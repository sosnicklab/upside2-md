#include "main.h"
#include "martini.h"
#include "monte_carlo_sampler.h"
#include "h5_support.h"
#include <tclap/CmdLine.h>
#include "deriv_engine.h"
#include "timing.h"
#include "thermostat.h"
#include <chrono>
#include <algorithm>
#include <set>
#include "random.h"
#include "state_logger.h"
#include <csignal>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "box.h"

using namespace std;
using namespace h5;

// If any stop signal is received (currently we trap sigterm and sigint)
// we increment any_stop_signal_received.
constexpr sig_atomic_t NO_SIGNAL = -1;  // FIXME is this a valid sentinel value?
volatile sig_atomic_t received_signal = NO_SIGNAL;

// if any stop signal is received, attempt to dump buffered state immediately.
// Perhaps I should also dump the exactly last position for each replica as
// well in case of an error, to aid debugging of say segfaults.  I need to be
// careful to stop each thread before draining the state to avoid race
// conditions.  This is all somewhat complicated, but worth it to have good
// stops no matter water.  I should also drain state in response to sigusr1 so
// that I can dump state before 100 frames if I need to while continuing the
// simulation.

static void abort_like_handler(int signal) {
    // NOTE TO THE INEXPERIENCED:
    //     This is a signal handler called in response to things like
    //     Ctrl-C.  These functions are very special, and there are very
    //     few things you are allowed to do in such a function without
    //     causing problems.  Basically, just set a global variable of type
    //     volatile sig_atomic_t and then notice the flag is set later in
    //     the code.  Even most forms of exiting the program are disallowed
    //     in a signal handler.  Please do not edit this function unless
    //     you have read extensively about signal handlers, especially for
    //     multithreaded code.

    // This signal handler is intended for signals that indicate the user
    // or OS wishes the program to terminate but after possibly writing any
    // remaining state.  There is some danger that a later SIGKILL could
    // corrupt the file being written to dump the last small amount of
    // state.  I think this is a risk worth taking.

    // Note that repeated Ctrl-C will not terminate the program if it is
    // still writing to the (possibly hanging) filesystem.  This is also
    // annoying, but the user can terminate with SIGKILL.

    received_signal = signal;
}

struct SignalHandlerHandler {
    // class to handle replacing signal handlers with orderly termination
    // If upside_main is used as a Python function, we must use RAII to
    // ensure than any termination will restore the Python signal handlers.
    typedef void (*signal_handler_t)(int);

    int signum;
    signal_handler_t old_handler;

    SignalHandlerHandler(int signum_, signal_handler_t handler):
        signum(signum_), old_handler(SIG_ERR)
    {
        // I will assume that I can use the standard signal system.  I am not 
        // sure how this will interact with Python, which can use sigaction.

        old_handler = signal(signum, handler);
        if(old_handler == SIG_ERR)
            fprintf(stderr, "Warning: problem installing signal handler."
                    " Does not affect correctness of simulation.\n");
    }

    virtual ~SignalHandlerHandler() {
        if(old_handler != SIG_ERR)
            if(signal(signum, old_handler) == SIG_ERR)
                fprintf(stderr, "Warning: problem restoring signal handler."
                    " Does not affect correctness of simulation.\n");
    }
};

    

struct System {
    int n_atom;
    uint32_t random_seed;
    float initial_temperature;
    float temperature;
    H5Obj config;
    H5Obj input;
    H5Obj output;
    shared_ptr<H5Logger> logger;
    DerivEngine engine;
    MultipleMonteCarloSampler mc_samplers;
    VecArrayStorage mom; // momentum
    OrnsteinUhlenbeckThermostat thermostat;
    uint64_t round_num;
    bool martini_hybrid_progress;
    vector<int> protein_rg_atom_indices;
    vector<unsigned char> protein_atom_mask;
    vector<unsigned char> lipid_atom_mask;
    int replica_index;  // the replica this system holds, which replica exchange moves between systems
    System(): round_num(0), martini_hybrid_progress(false), replica_index(0) {}

    void set_temperature(float new_temp) {
        temperature = new_temp;
        thermostat.set_temp(temperature);
    }
};

static vector<unsigned char> read_particle_class_mask(
        hid_t config_root,
        int n_atom,
        const string& target_class) {
    vector<unsigned char> out(static_cast<size_t>(n_atom), 0u);
    if(!h5_exists(config_root, "/input/particle_class")) return out;
    auto shape = get_dset_size(1, config_root, "/input/particle_class");
    if(static_cast<int>(shape[0]) != n_atom)
        throw string("/input/particle_class length must match n_atom");
    traverse_string_dset<1>(config_root, "/input/particle_class",
            [&](size_t i, const string& value) {
                if(value.find(target_class) != string::npos) out[i] = 1u;
            });
    return out;
}

static vector<int> collect_hybrid_ca_indices(hid_t config_root, int n_atom) {
    vector<int> out;
    if(n_atom <= 0) return out;
    if(!h5_exists(config_root, "/input/hybrid_bb_map/atom_indices")) return out;
    if(!h5_exists(config_root, "/input/hybrid_bb_map/atom_mask")) return out;

    try {
        auto idx_shape = get_dset_size(2, config_root, "/input/hybrid_bb_map/atom_indices");
        auto mask_shape = get_dset_size(2, config_root, "/input/hybrid_bb_map/atom_mask");
        if(idx_shape.size() != 2u || mask_shape.size() != 2u) return out;
        if(idx_shape[0] == 0u || idx_shape[1] < 2u) return out;
        if(mask_shape[0] != idx_shape[0] || mask_shape[1] != idx_shape[1]) return out;

        const size_t n_row = idx_shape[0];
        const size_t n_col = idx_shape[1];
        vector<int> atom_indices(n_row*n_col, -1);
        vector<int> atom_mask(n_row*n_col, 0);
        traverse_dset<2,int>(config_root, "/input/hybrid_bb_map/atom_indices",
                [&](size_t i, size_t j, int v) { atom_indices[i*n_col + j] = v; });
        traverse_dset<2,int>(config_root, "/input/hybrid_bb_map/atom_mask",
                [&](size_t i, size_t j, int v) { atom_mask[i*n_col + j] = v; });

        vector<unsigned char> seen(static_cast<size_t>(n_atom), 0u);
        out.reserve(n_row);
        for(size_t i = 0; i < n_row; ++i) {
            if(atom_mask[i*n_col + 1] == 0) continue;
            int ai = atom_indices[i*n_col + 1];
            if(ai < 0 || ai >= n_atom) continue;
            if(seen[static_cast<size_t>(ai)]) continue;
            seen[static_cast<size_t>(ai)] = 1u;
            out.push_back(ai);
        }
    } catch(...) {
        out.clear();
    }
    return out;
}

static double compute_protein_rg(const System& sys) {
    if(sys.protein_rg_atom_indices.empty()) return -1.0;
    float3 com = make_vec3(0.f, 0.f, 0.f);
    int n_backbone = 0;
    for(int ai: sys.protein_rg_atom_indices) {
        if(ai < 0 || ai >= sys.n_atom) continue;
        com += load_vec<3>(sys.engine.pos->output, ai);
        ++n_backbone;
    }
    if(n_backbone <= 0) return -1.0;
    com *= 1.f/n_backbone;
    double rg2 = 0.0;
    for(int ai: sys.protein_rg_atom_indices) {
        if(ai < 0 || ai >= sys.n_atom) continue;
        rg2 += mag2(load_vec<3>(sys.engine.pos->output, ai)-com);
    }
    return sqrt(rg2/n_backbone);
}

static double compute_all_atom_rg(const System& sys) {
    float3 com = make_vec3(0.f, 0.f, 0.f);
    for(int na=0; na<sys.n_atom; ++na)
        com += load_vec<3>(sys.engine.pos->output, na);
    com *= 1.f/sys.n_atom;

    double rg2 = 0.0;
    for(int na=0; na<sys.n_atom; ++na)
        rg2 += mag2(load_vec<3>(sys.engine.pos->output, na)-com);
    return sqrt(rg2/sys.n_atom);
}

static inline bool is_dry_martini_term_name(const std::string& name) {
    return is_prefix("martini_potential", name) ||
           is_prefix("dist_spring", name) ||
           is_prefix("angle_spring", name) ||
           is_prefix("dihedral_spring", name) ||
           is_prefix("restraint_position", name);
}

static inline bool is_sc_env_interface_term_name(const std::string& name) {
    return is_prefix("martini_sc_table_1body", name);
}

struct HybridProteinPotentialComponents {
    double upside_protein_potential = 0.0;
    double sc_env_interface_potential = 0.0;
    double bb_env_interface_potential = 0.0;
    double protein_total() const {
        return upside_protein_potential + sc_env_interface_potential + bb_env_interface_potential;
    }
};

static HybridProteinPotentialComponents compute_hybrid_protein_potential_components(const DerivEngine& engine) {
    HybridProteinPotentialComponents out;
    for(const auto& n : engine.nodes) {
        if(!n.computation->potential_term) continue;
        const auto& pot_node = dynamic_cast<const PotentialNode&>(*n.computation.get());
        double p = static_cast<double>(pot_node.potential);
        if(is_sc_env_interface_term_name(n.name)) {
            out.sc_env_interface_potential += p;
            continue;
        }
        if(is_dry_martini_term_name(n.name)) continue;
        out.upside_protein_potential += p;
    }
    out.bb_env_interface_potential = martini_hybrid::get_last_bb_env_interface_potential(engine);
    return out;
}

// Kinetic energy of the atoms in atom_mask per three degrees of freedom, so that it reads 1.5 kT in equilibrium:
// fixed atoms carry no momentum and z-held atoms none along z, so those degrees of freedom are not counted, and a
// rigid group, whose momenta are its rigid-body motion, has six however many atoms it holds.
static double compute_logged_subset_kinetic_energy(
        System* sys,
        const vector<unsigned char>& atom_mask) {
    if(atom_mask.size() != static_cast<size_t>(sys->n_atom)) return 0.0;
    vector<unsigned char> n_dof(static_cast<size_t>(sys->n_atom), 3u);
    for(int atom : martini_fix_rigid::get_z_fixed_atoms(sys->engine)) {
        if(atom >= 0 && atom < sys->n_atom) n_dof[static_cast<size_t>(atom)] = 2u;
    }
    for(int atom : martini_fix_rigid::get_fixed_atoms(sys->engine)) {
        if(atom >= 0 && atom < sys->n_atom) n_dof[static_cast<size_t>(atom)] = 0u;
    }
    double sum_kin = 0.0;
    int n_dof_total = 0;
    for(const auto& group : martini_fix_rigid::get_rigid_groups(sys->engine)) {
        if(!atom_mask[static_cast<size_t>(group[0])]) continue;
        for(int atom : group) {
            sum_kin += mag2(load_vec<3>(sys->mom, atom)) / martini_masses::get_mass(&sys->engine, atom);
            n_dof[static_cast<size_t>(atom)] = 0u;
        }
        n_dof_total += 6;
    }
    for(int na=0; na<sys->n_atom; ++na) {
        if(!atom_mask[static_cast<size_t>(na)] || !n_dof[static_cast<size_t>(na)]) continue;
        sum_kin += mag2(load_vec<3>(sys->mom, na)) / martini_masses::get_mass(&sys->engine, na);
        n_dof_total += n_dof[static_cast<size_t>(na)];
    }
    return n_dof_total > 0 ? (1.5/n_dof_total) * sum_kin : 0.0;
}

static double compute_logged_kinetic_energy(System* sys) {
    return compute_logged_subset_kinetic_energy(sys, vector<unsigned char>(static_cast<size_t>(sys->n_atom), 1u));
}

static vector<float> read_atom_thermostat_timescales(hid_t config_root, int n_atom) {
    vector<float> out;
    if(!h5_exists(config_root, "/input/thermostat_timescale")) return out;

    auto shape = get_dset_size(1, config_root, "/input/thermostat_timescale");
    if(shape.size() != 1u || shape[0] != static_cast<hsize_t>(n_atom)) {
        throw string("/input/thermostat_timescale must have one positive value per atom");
    }

    out.assign(static_cast<size_t>(n_atom), 0.f);
    traverse_dset<1,float>(config_root, "/input/thermostat_timescale",
            [&](size_t i, float tau) {
                if(!(tau > 0.f)) throw string("/input/thermostat_timescale values must be positive");
                out[i] = tau;
            });
    return out;
}

// Seed of system ns: a splitmix64 hash of (seed, ns).  Adding ns to the seed made system 1 of seed 5 run the
// same random streams as system 0 of seed 6.
static uint32_t system_seed(uint32_t seed, int ns) {
    uint64_t z = ((uint64_t(seed) << 32) | uint32_t(ns)) + 0x9e3779b97f4a7c15ull;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return uint32_t((z ^ (z >> 31)) >> 32);
}

double stod_strict(const std::string& s) {
    size_t nchar = -1u;
    double x = stod(s, &nchar);
    if(nchar != s.size()) throw("invalid float '" + s + "'");
    return x;
}
int stoi_strict(const std::string& s) {
    size_t nchar = -1u;
    int i = stoi(s, &nchar);
    if(nchar != s.size()) throw("invalid integer '" + s + "'");
    return i;
}

vector<string> split_string(const string& src, const string& sep) {
    vector<string> ret;

    for(auto curr_pos = begin(src); curr_pos-begin(src) < ptrdiff_t(src.size()); ) {
        auto next_match_end = search(curr_pos, end(src), begin(sep), end(sep));
        ret.emplace_back(curr_pos, next_match_end);
        curr_pos = next_match_end + sep.size();
    }

    return ret;
}


struct ReplicaExchange {
    struct SwapPair {int sys1; int sys2; uint64_t n_attempt; uint64_t n_success;};
    vector<vector<SwapPair>> swap_sets;
    vector<int> replica_indices;
    vector<vector<SwapPair*>> participating_swaps;
    int n_swap_sets;

    ReplicaExchange(vector<System>& systems, vector<string> swap_sets_strings) {
        n_swap_sets = (int)swap_sets_strings.size();
        int n_system = systems.size();
        for(int ns: range(n_system)) {
            replica_indices.push_back(systems[ns].replica_index);
            participating_swaps.emplace_back();
        }
        auto sorted_indices = replica_indices;
        sort(sorted_indices.begin(), sorted_indices.end());
        for(int ns: range(n_system))
            if(sorted_indices[ns] != ns)
                throw string("the replica indices of the systems (/input/mom@replica_index) are not a permutation");

        for(string& set_string: swap_sets_strings) {
            swap_sets.emplace_back();
            auto& set = swap_sets.back();
            
            for(auto pair_string: split_string(set_string,",")) {
                set.emplace_back();
                auto& s = set.back();
                auto p = split_string(pair_string, "-");
                if(p.size() != 2u) throw string("invalid swap pair, because it contains " +
                        to_string(p.size()) + "elements but should contain 2");

                s.sys1 = stoi_strict(p[0]);
                s.sys2 = stoi_strict(p[1]);

                s.n_attempt = 0u;
                s.n_success = 0u;

                int n_system = systems.size();
                if(s.sys1 >= n_system || s.sys2 >= n_system) throw string("invalid system");
            }
        }

        for(auto& ss: swap_sets) {
            set<int> systems_in_set;
            for(auto& sw: ss) {
                if(systems_in_set.count(sw.sys1) ||
                   systems_in_set.count(sw.sys2) ||
                   sw.sys1==sw.sys2) {
                    throw string("Overlapping indices in swap set.  "
                            "No replica index can appear more than once in a swap set.  "
                            "You probably (but maybe not; I didn't look that closely) need more "
                            "swap sets to get non-overlapping pairs.");
                }
                systems_in_set.insert(sw.sys1);
                systems_in_set.insert(sw.sys2);

                participating_swaps[sw.sys1].push_back(&sw);
                participating_swaps[sw.sys2].push_back(&sw);
            }
        }

        // enable logging of replica events
        for(int ns: range(n_system)) {
            auto logger = systems[ns].logger;
            if(!logger) continue;
            if(static_cast<int>(logger->level) < static_cast<int>(LOG_BASIC)) continue;

            auto* replica_ptr = &replica_indices[ns];
            logger->add_logger<int>("replica_index", {1}, [replica_ptr](int* buffer){buffer[0] = *replica_ptr;});

            auto* swap_vectors = &participating_swaps[ns];
            logger->log_once<int>("replica_swap_partner", {int(swap_vectors->size())},[swap_vectors,ns](int* buffer){
                    for(int i: range(swap_vectors->size())) {
                        auto sw = (*swap_vectors)[i];
                        // write down the system that is *not* this system
                        buffer[i] = (sw->sys1!=ns ? sw->sys1 : sw->sys2);
                    }});

            logger->add_logger<int>("replica_cumulative_swaps", {int(swap_vectors->size()), 2}, 
                  [swap_vectors](int* buffer) {
                    for(int i: range(swap_vectors->size())) {
                        auto sw = (*swap_vectors)[i];
                        buffer[2*i+0] = sw->n_success;
                        buffer[2*i+1] = sw->n_attempt;
                    }});
        }
    }

    void reset_stats() {
        for(auto& ss: swap_sets)
            for(auto& sw: ss)
                sw.n_success = sw.n_attempt = 0u;
    }

    // The swap sets alternate by attempt number, which the caller counts from the round so that a restart
    // continues the alternation.
    void attempt_swaps(uint32_t seed, uint64_t round, uint64_t n_attempt, vector<System>& systems) {
        int n_system = systems.size();

        vector<float> beta;
        for(auto &sys: systems) beta.push_back(1.f/sys.temperature);

        // compute the boltzmann factors for everyone
        auto compute_log_boltzmann = [&]() {
            vector<float> result(n_system);
            for(int i=0; i<n_system; ++i) {
                systems[i].engine.compute(PotentialAndDerivMode);
                result[i] = -beta[i]*systems[i].engine.potential;
            }
            return result;
        };

        auto coord_swap = [&](int ns1, int ns2) {
            swap(systems[ns1].engine.pos->output, systems[ns2].engine.pos->output);
        };

        // An accepted exchange carries the momenta with the coordinates, rescaled by sqrt(T_dest/T_src) so each
        // replica keeps its own kinetic temperature.  The integrator leaves the momentum half a step behind the
        // coordinates and correlated with the force there, so a coordinate-only swap breaks that correlation and
        // heats the replicas (KE/1.5kT 1.03 with a swap every 10 rounds; findings 3.14).
        auto mom_swap = [&](int ns1, int ns2) {
            swap(systems[ns1].mom, systems[ns2].mom);
            const float scale1 = sqrtf(systems[ns1].temperature / systems[ns2].temperature);
            const float scale2 = sqrtf(systems[ns2].temperature / systems[ns1].temperature);
            for(int ns: {ns1, ns2}) {
                const float scale = (ns == ns1) ? scale1 : scale2;
                VecArray m = systems[ns].mom;
                for(int na = 0; na < systems[ns].n_atom; ++na)
                    for(int d = 0; d < 3; ++d) m(d, na) *= scale;
            }
            swap(replica_indices[ns1], replica_indices[ns2]);
        };

        RandomGenerator random(seed, REPLICA_EXCHANGE_RANDOM_STREAM, 0u, round);

        // The energy is computed both before and after the swap, so that Hamiltonian replica exchange
        // works as well as temperature replica exchange.
        auto& set = swap_sets[n_attempt%n_swap_sets];
        auto old_lboltz = compute_log_boltzmann();
        for(auto& swap_pair: set) coord_swap(swap_pair.sys1, swap_pair.sys2);
        auto new_lboltz = compute_log_boltzmann();

        // reverse all swaps that should not occur by metropolis criterion
        for(auto& swap_pair: set) {
            auto s1 = swap_pair.sys1;
            auto s2 = swap_pair.sys2;
            swap_pair.n_attempt++;

            float lboltz_diff = (new_lboltz[s1] + new_lboltz[s2]) - (old_lboltz[s1]+old_lboltz[s2]);
            if(lboltz_diff < 0.f && expf(lboltz_diff) < random.uniform_open_closed().x()) {
                coord_swap(s1,s2);
            } else {
                swap_pair.n_success++;
                mom_swap(s1,s2);
            }
        }
    }
};

// The membrane term written with --use-curvature, empty if there is none: it reads Const3D_curvature_center and the
// absolute x and y of its atoms
static string curved_membrane_node(System& sys) {
    for(auto &n: sys.engine.nodes)
        for(auto type: {"cb_membrane_potential", "hb_membrane_potential",
                        "cb_surf_membrane_potential", "hb_surf_membrane_potential"})
            if(is_prefix(type, n.name) && read_attribute<int>(sys.config.get(),
                        ("/input/potential/" + n.name).c_str(), "use_curvature", 0))
                return n.name;
    return string();
}

struct CurvatureChange {
    vector<CoordNode*> centers;   // the centre of curvature of each system
    vector<float> min_radius;     // the membrane's half-thickness, below which it has no inner leaflet

    CurvatureChange(vector<System>& systems) {
        for(auto& sys: systems) {
            CoordNode* center = nullptr;
            for(auto &n: sys.engine.nodes)
                if(is_prefix("Const3D_curvature_center", n.name))
                    center = dynamic_cast<CoordNode*>(n.computation.get());
            string membrane = curved_membrane_node(sys);
            if(!center || membrane.empty())
                throw string("--curvature-changer-interval requires a membrane potential written with --use-curvature");
            centers.push_back(center);
            min_radius.push_back(read_attribute<float>(sys.config.get(), ("/input/potential/" + membrane).c_str(),
                                                       "half_thickness"));
            if(!(fabsf(center->output(2,0)) > min_radius.back()))
                throw string("the curvature radius must start above the membrane's half-thickness");
            // the changer moves the centre, so it is logged with every frame and a continuation starts from it
            sys.logger->add_logger<float>("curvature_center", {1, 3}, [center](float* buffer) {
                    for(int d=0; d<3; ++d) buffer[d] = center->output(d,0);});
        }
    }

    void attempt_change(uint32_t seed, uint64_t round, int ns, System& sys, float relative_change) {
        RandomGenerator random(seed, CURVATURE_CHANGE_RANDOM_STREAM, 0u, round);
        float4 u = random.uniform_open_closed();
        VecArray center = centers[ns]->output;

        // The step is symmetric in log|z|; the term -log_scale makes the stationary measure flat in curvature 1/R.
        // The measure is that of radii above the half-thickness, so a proposal at or below it is rejected: the
        // membrane energy is bounded, and a measure reaching R = 0 could not be normalised.
        float old_z     = center(2,0);
        float log_scale = relative_change*(2.f*u.x()-1.f);
        float new_z     = old_z*expf(log_scale);
        if(fabsf(new_z) <= min_radius[ns]) return;

        sys.engine.compute(PotentialAndDerivMode);
        float old_potential = sys.engine.potential;
        center(2,0) = new_z;

        sys.engine.compute(PotentialAndDerivMode);
        float new_potential = sys.engine.potential;

        float lboltz_diff = -log_scale - (new_potential-old_potential)/sys.temperature;
        if(lboltz_diff < 0.f && expf(lboltz_diff) < u.y())
            center(2,0) = old_z;
    }
};


vector<float> potential_deriv_agreement(DerivEngine& engine) {
    vector<float> relative_error;
    int n_atom = engine.pos->n_elem;
    VecArray pos_array = engine.pos->output;

    vector<float> input(n_atom*3);
    for(int na=0; na<n_atom; ++na)
        for(int d=0; d<3; ++d)
            input[na*3+d] = pos_array(d,na);
    std::vector<float> output(1);

    auto do_compute = [&]() {
        for(int na=0; na<n_atom; ++na)
            for(int d=0; d<3; ++d)
                pos_array(d,na) = input[na*3+d];
        engine.compute(PotentialAndDerivMode);
        output[0] = engine.potential;
    };

    for(auto &n: engine.nodes) {
        if(n.computation->potential_term) {
            auto &v = dynamic_cast<PotentialNode&>(*n.computation.get()).potential;
            printf("%s: % 4.3f\n", n.name.c_str(), v);
        }
    }
    printf("\n\n");

    auto central_diff_jac = central_difference_deriviative(do_compute, input, output, 1e-3);
    vector<float> deriv_array;
    for(int na=0; na<n_atom; ++na)
        for(int d=0; d<3; ++d)
            deriv_array.push_back(engine.pos->sens(d,na));

    relative_error.push_back(
            relative_rms_deviation(central_diff_jac, deriv_array));
    return relative_error;
}

int upside_main(int argc, const char* const * argv, int verbose=1)
try {
    received_signal = NO_SIGNAL;  // an earlier in-process run may have been interrupted
    using namespace TCLAP;  // Templatized C++ Command Line Parser (tclap.sourceforge.net)
    CmdLine cmd("Using Protein Statistical Information for Dynamics Estimation (Upside)\n Author: John Jumper, Xiangda Peng, Nabil Faruk", 
            ' ', "2.0alpha");

    ValueArg<double> time_step_arg("", "time-step", "time step for integration (default 0.009)", 
            false, 0.009, "float", cmd);
    ValueArg<double> duration_arg("", "duration", "duration of simulation", 
            false, -1., "float", cmd);
    ValueArg<int> duration_steps_arg("", "duration-steps",
            "number of integration rounds to run (step count; overrides --duration when >=0)",
            false, -1, "int", cmd);
    ValueArg<double> time_lim_arg("", "time-limit", "Run time limit of simulation in seconds (default unlimited)", 
            false, -1., "float", cmd);
    ValueArg<unsigned long> seed_arg("", "seed", "random seed (default 42)", 
            false, 42l, "int", cmd);
    ValueArg<string> temperature_arg("", "temperature", "thermostat temperature, as a comma-separated list with "
            "one temperature per system, or a single temperature for all systems (required)", 
            false, "", "temperature_list", cmd);
    MultiArg<string> swap_set_args("","swap-set", "list like 0-1,2-3,6-7,4-8 of non-overlapping swaps for a replica "
            "exchange.  May be specified multiple times for multiple swap sets (non-overlapping is only required "
            "within a single swap set).", false, "h5_files");
    cmd.add(swap_set_args);
    ValueArg<double> anneal_factor_arg("", "anneal-factor", "annealing factor (0.1 means the final temperature "
            "will be 10% of the initial temperature)", 
            false, 1., "float", cmd);
    ValueArg<double> anneal_duration_arg("", "anneal-duration", "duration of annealing phase "
            "(default is duration of simulation)", 
            false, -1., "float", cmd);

    ValueArg<double> anneal_start_arg("", "anneal-start", "start annealing phase ",
            false, -1., "float", cmd);

    ValueArg<double> anneal_end_arg("", "anneal-end", "end annealing phase ", 
            false, -1., "float", cmd);

    ValueArg<double> frame_interval_arg("", "frame-interval", "simulation time between frames", 
            true, -1., "float", cmd);
    ValueArg<double> replica_interval_arg("", "replica-interval", 
            "simulation time between applications of replica exchange (0 means no replica exchange, default 0.)", 
            false, 0., "float", cmd);
    ValueArg<double> mc_interval_arg("", "monte-carlo-interval", 
            "simulation time between attempts to do Monte Carlo moves (0. means no MC moves, default 0.)", 
            false, 0., "float", cmd);
    ValueArg<double> thermostat_interval_arg("", "thermostat-interval", 
            "simulation time between applications of the thermostat", 
            false, -1., "float", cmd);
    ValueArg<double> thermostat_timescale_arg("", "thermostat-timescale", "timescale for the thermostat", 
            false, 5., "float", cmd);

    ValueArg<double> curvature_changer_interval_arg("", "curvature-changer-interval", 
            "simulation time between applications of curvature change (0 means no curvature change, default 0.)", 
            false, 0., "float", cmd);

    ValueArg<double> relative_curvature_radius_change_arg("", "relative-curvature-radius-change", 
            "Maximum amount of curvature radius adjustment (default 0.05)", 
            false, 0.05, "float", cmd);

    ValueArg<string> integrator_arg("", "integrator", 
            "Use this option to control which Integrator are used.  Available levels are v(verlet) or mv(multi-step verlet). "
            "Default is verlet.",
            false, "", "v, mv ", cmd);
    ValueArg<int> inner_step_arg("", "inner-step", "inner step for the integrator", false, 3, "int", cmd);
    SwitchArg enable_min_arg("", "minimize", "Run an energy minimization before MD", cmd, false);
    ValueArg<int> min_max_iter_arg("", "min-max-iter", "Minimization max iterations", false, 1000, "int", cmd);
    ValueArg<double> min_energy_tol_arg("", "min-energy-tol", "Minimization energy tolerance", false, 1e-6, "float", cmd);
    ValueArg<double> min_force_tol_arg("", "min-force-tol", "Minimization force tolerance", false, 1e-6, "float", cmd);
    ValueArg<double> min_step_arg("", "min-step", "Minimization initial step: the largest-force atom's displacement in A", false, 0.1, "float", cmd);

    ValueArg<string> input_arg("i", "input", "h5df input file for position", false, "not_Defined_By_user", "string", cmd);
    ValueArg<string> input_base_arg("", "input-base", "h5df input files base for positions", false, "not_Defined_By_user", "string_list", cmd);
    ValueArg<string> output_arg("o", "output", "h5df output log file", false, "not_Defined_By_user", "string", cmd);
    ValueArg<string> output_base_arg("", "output-base", "the base name of h5df output log files.", false, "not_Defined_By_user", "string", cmd);
    SwitchArg disable_recenter_arg("", "disable-recentering", 
            "Disable all recentering of protein in the universe", 
            cmd, false);
    SwitchArg disable_z_recenter_arg("", "disable-z-recentering", 
            "Disable z-recentering of protein in the universe", 
            cmd, false);
    SwitchArg raise_signal_on_exit_if_received_arg("", "re-raise-signal", 
            "(Developer use only) Used for obscure details of signal handling.  No effect on simulation.", 
            cmd, false);
    ValueArg<string> log_level_arg("", "log-level", 
            "Use this option to control which arrays are stored in /output.  Available levels are basic, detailed, "
            "or extensive.  Default is detailed.",
            false, "", "basic, detailed, extensive", cmd);
    SwitchArg potential_deriv_agreement_arg("", "potential-deriv-agreement",
            "(developer use only) check the agreement of the derivative with finite differences "
            "of the potential for the initial structure.  This may give strange answers for native structures "
            "(no steric clashes may given an agreement of NaN) or random structures (where bonds and angles are "
            "exactly at their equilibrium values).  Interpret these results at your own risk.", cmd, false);
    SwitchArg record_momentum_arg("", "record-momentum",
            "record the momentum (so that the trajectory can be exactly restarted)"
            " the momentum will be recorded in output.mom ",
            cmd, false);
    SwitchArg restart_using_momentum_arg("", "restart-using-momentum",
            "restart the trajectory by initializing using the momentum stored in input.mom,"
            " this is usually helpful when the simulation ends due to the wall time."
            " To restart the trajectory, make sure to copy the end momentum of last run to input.mom (similar to the pos treatment when restarting):"
            " check examples for continue the simulation",
            cmd, false);
    ValueArg<string> set_param_arg("", "set-param", "Developer use only", false, "", "param_arg", cmd);
    UnlabeledMultiArg<string> config_args("config_files","configuration .h5 files", true, "h5_files");
    cmd.add(config_args);
    cmd.parse(argc, argv);

    try {
        if(verbose) printf("invocation: ");
        std::string invocation(argv[0]);
        for(auto arg=argv+1; arg!=argv+argc; ++arg) invocation += string(" ") + *arg;
        if(verbose) printf("%s\n", invocation.c_str());

        map<string,vector<float>> set_param_map;
        if(set_param_arg.getValue().size()) {
            auto param_file = h5_obj(H5Fclose, H5Fopen(set_param_arg.getValue().c_str(),
                        H5F_ACC_RDONLY, H5P_DEFAULT));

            for(const string& node_name: node_names_in_group(param_file.get(), ".")) {
                set_param_map[node_name] = vector<float>();
                auto& values = set_param_map[node_name];
                traverse_dset<1,float>(param_file.get(), node_name.c_str(), [&](size_t i, float x) {
                        values.push_back(x);});
            }
        }


        float dt = time_step_arg.getValue();
        if(!(dt > 0.f)) throw string("--time-step must be positive");
        if(!(frame_interval_arg.getValue() > 0.)) throw string("--frame-interval must be positive");
        if(!(thermostat_timescale_arg.getValue() > 0.)) throw string("--thermostat-timescale must be positive");
        if(!(anneal_factor_arg.getValue() > 0.)) throw string("--anneal-factor must be positive");
        if(replica_interval_arg.getValue() < 0. || mc_interval_arg.getValue() < 0. ||
                curvature_changer_interval_arg.getValue() < 0.)
            throw string("--replica-interval, --monte-carlo-interval and --curvature-changer-interval must not be negative");
        double duration = duration_arg.getValue();
        int duration_steps = duration_steps_arg.getValue();
        double time_lim = time_lim_arg.getValue();
        bool passed_time_lim = false;
        bool use_duration_steps = duration_steps >= 0;

        // A round is 3 time steps of the 3-stage Verlet integrator, inner_step steps of multi-step Verlet,
        // and one step of the g-JF integrator that configurations with /input/brownian use.
        bool brownian = h5_exists(h5_obj(H5Fclose, H5Fopen(config_args.getValue()[0].c_str(),
                        H5F_ACC_RDONLY, H5P_DEFAULT)).get(), "/input/brownian");
        if(integrator_arg.getValue() != "" && integrator_arg.getValue() != "v" && integrator_arg.getValue() != "mv")
            throw string("--integrator must be v or mv");
        if(integrator_arg.getValue() == "mv" && inner_step_arg.getValue() < 1)
            throw string("--inner-step must be at least 1");
        int inner_step = brownian ? 1 : integrator_arg.getValue() == "mv" ? inner_step_arg.getValue() : 3;

        uint64_t n_round = 0;
        if(use_duration_steps) {
            n_round = static_cast<uint64_t>(duration_steps);
            duration = static_cast<double>(n_round) * inner_step * dt;
        } else {
            if(duration < 0.) {
                throw string("Must provide either --duration or --duration-steps");
            }
            n_round = round(duration / (inner_step*dt));
        }
        int thermostat_interval = max(1.,round(thermostat_interval_arg.getValue() / (inner_step*dt))); //  FIXME inner_step
        int frame_interval = max(1.,round(frame_interval_arg.getValue() / (inner_step*dt)));

        unsigned long big_prime = 4294967291ul;  // largest prime smaller than 2^32
        uint32_t base_random_seed = uint32_t(seed_arg.getValue() % big_prime);

        // initialize thermostat and thermalize momentum
        if(verbose) printf("random seed: %lu\n", (unsigned long)(base_random_seed));

        int mc_interval = mc_interval_arg.getValue() > 0. 
            ? max(1.,round(mc_interval_arg.getValue()/(inner_step*dt))) 
            : 0;

        double display_duration_total = duration;
        int duration_print_width = ceil(log(1+display_duration_total)/log(10));

        bool do_recenter = !disable_recenter_arg.getValue();
        bool xy_recenter_only = do_recenter && disable_z_recenter_arg.getValue();

        h5_noerr(H5Eset_auto(H5E_DEFAULT, nullptr, nullptr));
        vector<string> config_paths = config_args.getValue();
        vector<System> systems(config_paths.size());
        // drop the systems' MARTINI registry entries when this run ends, however it ends
        struct UnregisterSystems {
            vector<System>& systems;
            ~UnregisterSystems() {for(auto& sys: systems) martini_unregister_engine(sys.engine);}
        } unregister_systems{systems};

        auto temperature_strings = split_string(temperature_arg.getValue(), ",");
        if(temperature_strings.size() != 1u && temperature_strings.size() != systems.size()) 
            throw string("Received "+to_string(temperature_strings.size())+" temperatures but have "
                    +to_string(systems.size())+" systems");

        for(int ns: range(systems.size())) {
            float T = stod_strict(temperature_strings.size()>1u ? temperature_strings[ns] : temperature_strings[0]);
            if(!(T > 0.f)) throw string("--temperature values must be positive");
            systems[ns].initial_temperature = T;
        }

        double anneal_factor = anneal_factor_arg.getValue();
        double anneal_duration = anneal_duration_arg.getValue();
        if(anneal_duration == -1.) anneal_duration = duration;

        double anneal_start = anneal_start_arg.getValue();
        if (anneal_start < 0.0) anneal_start = 0.0;
        double anneal_end = anneal_end_arg.getValue();
        if (anneal_end > anneal_start) anneal_duration = anneal_end-anneal_start;
        else anneal_end = anneal_start+anneal_duration;

        // tighter spacing at the low end of temperatures because the variance is decreasing
        auto anneal_temp = [=](double initial_temperature, double time) {
            auto fraction = max(0., (time-anneal_start) / anneal_duration);
            double T0 = initial_temperature;
            double T1 = initial_temperature*anneal_factor;
            if (time <= anneal_start)
                return T0;
            else if (time > anneal_end)
                return T1;
            else 
                return sqr(sqrt(T0)*(1.-fraction) + sqrt(T1)*fraction);
        };

        vector<string> outputs(systems.size());
        bool user_defined_output = false;

        auto out_name = output_arg.getValue();
        if (out_name != "not_Defined_By_user") {
            auto output_strings = split_string(out_name, ",");

            if(output_strings.size() != systems.size())
                throw string("Received "+to_string(output_strings.size())+" outputs but have "
                        +to_string(systems.size())+" systems; each system needs its own output file");

            for(int ns: range(systems.size())) 
                outputs[ns] = output_strings[ns];

            user_defined_output = true;
        }

        auto out_base = output_base_arg.getValue();
        if (out_base != "not_Defined_By_user")  {
            for(int ns: range(systems.size())) 
                outputs[ns] = out_base + "_" + to_string(ns) + ".h5";
            user_defined_output = true;
        }
        
        vector<string> inputs(systems.size());
        bool user_defined_input = false;

        auto in_name = input_arg.getValue();
        if (in_name != "not_Defined_By_user") {
            auto input_strings = split_string(in_name, ",");

            if(input_strings.size() != 1u && input_strings.size() != systems.size()) 
                throw string("Received "+to_string(input_strings.size())+" inputs but have "
                        +to_string(systems.size())+" systems");

            for(int ns: range(systems.size())) 
                inputs[ns] = input_strings.size()>1u ? input_strings[ns] : input_strings[0];

            user_defined_input = true;
        }
        auto in_base = input_base_arg.getValue();
        if (in_base != "not_Defined_By_user")  {
            for(int ns: range(systems.size())) 
                inputs[ns] = in_base + "_" + to_string(ns) + ".h5";
            user_defined_input = true;
        }
        int replica_interval = 0;
        if(replica_interval_arg.getValue())
            replica_interval = max(1.,round(replica_interval_arg.getValue()/(inner_step*dt)));

        int curvature_changer_interval = 0;
        if(curvature_changer_interval_arg.getValue())
            curvature_changer_interval = max(1.,round(curvature_changer_interval_arg.getValue()/(inner_step*dt)));
        float relative_curvature_radius_change = relative_curvature_radius_change_arg.getValue();

        // system 0 is the minimum temperature
        int n_system = systems.size();

        // We are not allowed to exit an OpenMP critical section early.  For this reason, we must trap
        // all exceptions.  To avoid crashing callers, we simply record the presence of an exception
        // then exit immediately after the block.
        bool error_exit_omp = false;
        #pragma omp critical
        for(int ns=0; ns<n_system; ++ns) try {
            System* sys = &systems[ns];  // a pointer here makes later lambda's more natural
            sys->random_seed = system_seed(base_random_seed, ns);

            try {
                if (user_defined_output) {
                    sys->config = h5_obj(H5Fclose,
                        H5Fopen(config_paths[ns].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
                }
                else {
                    sys->config = h5_obj(H5Fclose,
                        H5Fopen(config_paths[ns].c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
                }
            } catch(string &s) {
                throw string("Unable to open configuration file at ") + config_paths[ns];
            }

            if (user_defined_input) {
                try {
                    sys->input = h5_obj(H5Fclose,
                            H5Fopen(inputs[ns].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
                } catch(string &s) {
                    throw string("Unable to open input position file ") + inputs[ns];
                }
            }

            if (user_defined_output) {
                sys->output = h5_obj(H5Fclose, 
                   H5Fcreate(outputs[ns].c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
            }
            else {
                if(h5_exists(sys->config.get(), "output")) {
                    // Note that it is not possible in HDF5 1.8.x to reclaim space by deleting
                    // datasets or groups.  Subsequent h5repack will reclaim space, however.
                    h5_noerr(H5Ldelete(sys->config.get(), "/output", H5P_DEFAULT));
                }
            }


            LogLevel log_level;
            if     (log_level_arg.getValue() == "")          log_level = LOG_DETAILED;
            else if(log_level_arg.getValue() == "basic")     log_level = LOG_BASIC;
            else if(log_level_arg.getValue() == "detailed")  log_level = LOG_DETAILED;
            else if(log_level_arg.getValue() == "extensive") log_level = LOG_EXTENSIVE;
            else throw string("Illegal value for --log-level");

            if (user_defined_output) 
                sys->logger = make_shared<H5Logger>(sys->output, "output", log_level);
            else
                sys->logger = make_shared<H5Logger>(sys->config, "output", log_level);

            default_logger = sys->logger;  // FIXME kind of a hack for the ugly global variable

            if (user_defined_output) 
                write_string_attribute(sys->output.get(), "output", "invocation", invocation);
            else
                write_string_attribute(sys->config.get(), "output", "invocation", invocation);

            auto pos_shape = get_dset_size(3, sys->config.get(), "/input/pos");
            if (user_defined_input) 
                pos_shape = get_dset_size(3, sys->input.get(), "/input/pos");

            sys->n_atom = pos_shape[0];

            if(pos_shape[1]!=3) throw string("invalid dimensions for initial position");
            if(pos_shape[2]!=1) throw string("must have n_system 1 from config");

            auto potential_group = open_group(sys->config.get(), "/input/potential");
            sys->engine = initialize_engine_from_hdf5(sys->n_atom, potential_group.get());
            // Register NPT barostat for this engine (reads settings from H5)
            simulation_box::npt::register_barostat_for_engine(sys->config.get(), sys->engine, sys->random_seed);
            // Load masses for MARTINI integrators (read from H5)
            martini_masses::load_masses_for_engine(&sys->engine, sys->config.get());
            // Register fix rigid settings for this engine (read from H5)
            martini_fix_rigid::register_fix_rigid_for_engine(sys->config.get(), sys->engine);
            // Register stage-specific parameters for this engine (read from H5)
            martini_stage_params::register_stage_params_for_engine(&sys->engine, sys->config.get());
            // Register hybrid MARTINI/Upside metadata for this engine (read from H5)
            martini_hybrid::register_hybrid_for_engine(sys->config.get(), sys->engine);
            // Single-step g-JF integrator for every physical MARTINI particle
            martini_brownian::register_brownian_for_engine(
                    &sys->engine, sys->config.get(), sys->random_seed);
            if(martini_brownian::has_brownian(&sys->engine) != brownian)
                throw string("all configurations must use the same integrator, with or without /input/brownian");
            sys->martini_hybrid_progress = martini_hybrid::is_hybrid_enabled(sys->engine);
            sys->protein_atom_mask = read_particle_class_mask(
                    sys->config.get(), sys->n_atom, "PROTEIN");
            sys->lipid_atom_mask = read_particle_class_mask(
                    sys->config.get(), sys->n_atom, "LIPID");
            if(sys->martini_hybrid_progress) {
                sys->protein_rg_atom_indices = collect_hybrid_ca_indices(sys->config.get(), sys->n_atom);
            } else {
                sys->protein_rg_atom_indices.clear();
            }
            if  (integrator_arg.getValue() == "mv" )
                sys->engine.build_integrator_levels(true, dt, inner_step );
            // A moving constant advances by its time_step at each force evaluation, so it needs time_step = dt and one
            // evaluation per time step: each stage of the Verlet integrator, or the fast level of the multi-step one.
            // The g-JF step evaluates the forces 1 + inner_steps times per time step, the slow level once per round.
            for(auto& n: sys->engine.nodes) {
                if(!is_prefix("MovingConst3D", n.name)) continue;
                if(read_attribute<float>(potential_group.get(), n.name.c_str(), "time_step") != dt)
                    throw string("node ") + n.name + " needs a time_step equal to the integrator time step (--time-step)";
                if(brownian)
                    throw string("node ") + n.name + " cannot follow the g-JF integrator's clock (/input/brownian)";
                if(integrator_arg.getValue() == "mv" and n.integrator_level != 0)
                    throw string("node ") + n.name + " must feed only fast-level terms of the multi-step integrator";
            }

            // Override parameters as instructed by users
            for(const auto& p: set_param_map) {
                auto& computation = sys->engine.get(p.first).computation;
                if(computation->get_param().size() != p.second.size())
                    throw string("--set-param ") + p.first + ": expected " +
                        to_string(computation->get_param().size()) + " parameters but got " + to_string(p.second.size());
                computation->set_param(p.second);
            }


            if (user_defined_input) 
                traverse_dset<3,float>(sys->input.get(), "/input/pos", [&](size_t na, size_t d, size_t ns, float x) { 
                        sys->engine.pos->output(d,na) = x;});
            else
                traverse_dset<3,float>(sys->config.get(), "/input/pos", [&](size_t na, size_t d, size_t ns, float x) { 
                        sys->engine.pos->output(d,na) = x;});

            if(verbose) printf("%s\nn_atom %i\n\n", config_paths[ns].c_str(), sys->n_atom);

            if(potential_deriv_agreement_arg.getValue()){
                sys->engine.compute(PotentialAndDerivMode);
                if(verbose) printf("Initial potential:\n");
                auto relative_error = potential_deriv_agreement(sys->engine);
                if(verbose) printf("overall potential relative error: ");
                for(auto r: relative_error) printf(" %.5f", r);
                if(verbose) printf("\n");
            }

            sys->thermostat = OrnsteinUhlenbeckThermostat(
                    sys->random_seed,
                    thermostat_timescale_arg.getValue(),
                    1.,
                    1e8);
            sys->set_temperature(sys->initial_temperature);
            martini_brownian::set_brownian_temperature(&sys->engine, sys->temperature);
            sys->thermostat.set_delta_t(thermostat_interval*inner_step*dt);  // set true thermostat interval  //  FIXME inner_step
            sys->thermostat.set_atom_timescale(read_atom_thermostat_timescales(sys->config.get(), sys->n_atom));

            sys->mom.reset(3, sys->n_atom);
            // The replica this system holds.  A continuation carries it whether or not it resumes the momentum
            // (a position-only restart, as the hybrid REMD driver makes, keeps its walkers too).
            {
                hid_t restart_source = user_defined_input ? sys->input.get() : sys->config.get();
                sys->replica_index = h5_exists(restart_source, "input/mom")
                    ? read_attribute<int>(restart_source, "/input/mom", "replica_index", ns) : ns;
            }
            if (restart_using_momentum_arg.getValue()) { // initialize momentum using input.mom if requested
                if (user_defined_input) {
                    if (h5_exists(sys->input.get(), "input/mom")) {
                        check_size(sys->input.get(), "/input/mom", sys->n_atom, 3, 1);
                        traverse_dset<3,float>(sys->input.get(), "/input/mom", [&](size_t na, size_t d, size_t ns, float x) { 
                            sys->mom(d,na) = x;});
                    }
                    else throw  string("input h5 file doesn't have input.mom group, can't restart using the momentum!");
                }
                else {
                    if (h5_exists(sys->config.get(), "input/mom")) {
                        check_size(sys->config.get(), "/input/mom", sys->n_atom, 3, 1);
                        traverse_dset<3,float>(sys->config.get(), "/input/mom", [&](size_t na, size_t d, size_t ns, float x) { 
                                sys->mom(d,na) = x;});
                    }
                    else throw  string("input h5 file doesn't have input.mom group, can't restart using the momentum!");
                }

                // Resume the counters of the run that recorded the momentum, so that the continuation draws
                // the random numbers of the continuous run.  The initial thermalization is invocation 0.
                // Momentum recorded before 2026-10-19 carries no round and restarts at round 0, as it did then.
                int recorded_round = read_attribute<int>(user_defined_input ? sys->input.get() : sys->config.get(),
                        "/input/mom", "round", 0);
                if(recorded_round < 0) throw string("/input/mom@round must not be negative");
                sys->round_num = recorded_round;
                sys->thermostat.set_n_invocations(1 + (sys->round_num + thermostat_interval - 1)/thermostat_interval);
                martini_brownian::set_brownian_invocations(&sys->engine, sys->round_num);
                // and the annealing temperature that the last thermostat round before this one set
                if(anneal_factor != 1. && sys->round_num) {
                    uint64_t last_thermostat_round = (sys->round_num-1)/thermostat_interval*thermostat_interval;
                    sys->set_temperature(anneal_temp(sys->initial_temperature, inner_step*dt*(last_thermostat_round+1)));
                    martini_brownian::set_brownian_temperature(&sys->engine, sys->temperature);
                }
            }
            else {
                for(int d: range(3)) for(int na: range(sys->n_atom)) sys->mom(d,na) = 0.f;
                sys->thermostat.apply(sys->mom, sys->n_atom, &sys->engine, true, true); // initial thermalization if it's a fresh start
            }

            // the PO4 height hold of production runs on the round clock a continuation from the momentum carries
            martini_hybrid::refresh_transition_holds_for_engine(sys->engine, sys->round_num,
                                                                restart_using_momentum_arg.getValue());
            // Hybrid O/BB virtual sites are derived from the carrier set and
            // must not retain independently thermalized momentum.
            martini_fix_rigid::apply_fix_rigid_md(sys->engine, sys->engine.pos->output, sys->mom);


            // we must capture the sys pointer by value here so that it is available later
            sys->logger->add_logger<float>("pos", {1, sys->n_atom, 3}, [sys](float* pos_buffer) {
                    martini_hybrid::copy_positions_for_output(sys->engine, pos_buffer);
                    });
            if (record_momentum_arg.getValue()) { // record the momentum if requested, with the same frequency as the position recording
                sys->logger->add_logger<float>("mom", {1, sys->n_atom, 3}, [sys](float* mom_buffer) {
                        VecArray mom_array = sys->mom;
                        for(int na=0; na<sys->n_atom; ++na) 
                        for(int d=0; d<3; ++d) 
                        mom_buffer[na*3 + d] = mom_array(d,na);
                        });
                sys->logger->add_logger<int>("round", {}, [sys](int* round_buffer) {
                        *round_buffer = sys->round_num;});
            }
            sys->logger->add_logger<double>("kinetic", {1}, [sys](double* kin_buffer) {
                    kin_buffer[0] = compute_logged_kinetic_energy(sys);
                    });
            if(sys->martini_hybrid_progress) {
                sys->logger->add_logger<double>("protein_kinetic", {1}, [sys](double* kin_buffer) {
                        kin_buffer[0] = compute_logged_subset_kinetic_energy(sys, sys->protein_atom_mask);
                        });
                sys->logger->add_logger<double>("lipid_kinetic", {1}, [sys](double* kin_buffer) {
                        kin_buffer[0] = compute_logged_subset_kinetic_energy(sys, sys->lipid_atom_mask);
                        });
            }
            sys->logger->add_logger<double>("potential", {1}, [sys](double* pot_buffer) {
                    sys->engine.compute(PotentialAndDerivMode);
                    pot_buffer[0] = sys->engine.potential;});
            sys->logger->add_logger<double>("time", {}, [sys,dt,inner_step](double* time_buffer) {
                    *time_buffer=inner_step*double(dt)*sys->round_num;});

            if(mc_interval) {
                sys->mc_samplers = MultipleMonteCarloSampler{open_group(sys->config.get(), "/input").get(), *sys->logger,
                                                              sys->n_atom};
            }

            // Re-centering translates the system every frame, so it is refused for terms that depend on absolute
            // positions.  The names are those the config writers give these terms.
            auto has_node = [&](initializer_list<const char*> types) {
                for(auto &n: sys->engine.nodes)
                    for(auto type: types)
                        if(is_prefix(type, n.name)) return true;
                return false;
            };
            // a curved membrane, tethers, walls, pulling and restraints depend on x and y too
            bool curved_membrane = !curved_membrane_node(*sys).empty();
            if(do_recenter && (curved_membrane ||
                    has_node({"WallSpring_wall", "WallSpring_cavity_radial", "cavity_radial", "Spring_pulling",
                              "Spring_spring", "Spring_nail", "Spring_group_const", "Spring_group_nail",
                              "restraint_position"})))
                throw string("You have re-centering and a potential on absolute positions turned on.  "
                        "This is not what you want.  Consider --disable-recentering.");

            if(do_recenter && !xy_recenter_only &&
                    has_node({"cb_membrane_potential", "hb_membrane_potential",
                              "cb_surf_membrane_potential", "hb_surf_membrane_potential", "surface", "tension"}))
                throw string("You have z-centering and a z-dependent potential turned on.  "
                        "This is not what you want.  Consider --disable-z-recentering "
                        "or --disable-recentering.");

            if(curved_membrane && has_node({"surface"}))
                throw string("the surface node finds the lipid-facing residues in a flat slab, so it cannot serve a "
                             "curved membrane");

            if(do_recenter && (martini_fix_rigid::get_fixed_atoms(sys->engine).size() ||
                               martini_fix_rigid::get_z_fixed_atoms(sys->engine).size()))
                throw string("You have re-centering and fixed atoms turned on.  "
                        "This is not what you want.  Consider --disable-recentering.");


        } catch(const string &e) {
            fprintf(stderr, "\n\nERROR: %s\n", e.c_str());
            error_exit_omp = true;
        } catch(...) {
            fprintf(stderr, "\n\nERROR: unknown error\n");
            error_exit_omp = true;
        }
        // We have just left the critical section.  The global logger is released before any return: held, it
        // keeps the refused system's file open in a caller of the library, whose later engines would build their
        // node loggers into that file.
        default_logger = shared_ptr<H5Logger>();  // FIXME kind of a hack for the ugly global variable
        if(error_exit_omp) return 2;

        // A continuation runs the requested rounds on from the round its momentum was recorded at
        uint64_t start_round = systems[0].round_num;
        for(System& sys: systems)
            if(sys.round_num != start_round) throw string("all systems must restart from the same round");
        n_round += start_round;
        display_duration_total = n_round*double(dt*inner_step);
        duration_print_width = ceil(log(1+display_duration_total)/log(10));

        if(mc_interval) {
            for(System& sys: systems)
                if(martini_fix_rigid::get_fixed_atoms(sys.engine).size() || martini_fix_rigid::get_z_fixed_atoms(sys.engine).size())
                    throw string("Monte Carlo moves would move fixed atoms; do not combine them with --monte-carlo-interval");
        }
        // Half kick by which the momentum trails the coordinates, carried to the new coordinates of an accepted
        // Monte Carlo move (the g-JF momentum is synchronous with the coordinates)
        float mc_half_kick = brownian ? 0.f : 0.5f*dt;
        float mc_slow_half_kick = !brownian && integrator_arg.getValue() == "mv" ? 0.5f*inner_step*dt : mc_half_kick;

        unique_ptr<ReplicaExchange> replex;
        if(replica_interval) {
            if(verbose) printf("initializing replica exchange\n");
            replex.reset(new ReplicaExchange(systems, swap_set_args.getValue()));
            if(!replex->swap_sets.size()) throw string("replica exchange requested but no swap sets proposed");
        }


        if(replica_interval) {
            int n_atom = systems[0].n_atom;
            for(System& sys: systems) {
                if(simulation_box::npt::is_enabled(sys.engine))
                    throw string("Replica exchange does not exchange the barostat box; run it at constant volume");
                if(sys.n_atom != n_atom)
                    throw string("Replica exchange requires all systems have the same number of atoms");
            }
        }

        unique_ptr<CurvatureChange> curvature_changer;
        if(curvature_changer_interval) {
            if(verbose) printf("initializing curvature changer\n");
            curvature_changer.reset(new CurvatureChange(systems));
        }

        if(verbose) printf("\n");
        for(int ns: range(systems.size())) {
            if(verbose) printf("%i %.2f\n", ns, systems[ns].temperature);
            float* temperature_pointer = &(systems[ns].temperature);
            systems[ns].logger->add_logger<double>("temperature", {1}, [temperature_pointer](double* temperature_buffer) {
                    temperature_buffer[0] = *temperature_pointer;});

            // Add NPT-specific loggers if barostat is enabled
            if(simulation_box::npt::is_enabled(systems[ns].engine)) {
                // Log box dimensions
                systems[ns].logger->add_logger<float>("box", {3}, [ns, &systems](float* buffer) {
                    float bx, by, bz;
                    simulation_box::npt::get_current_box(systems[ns].engine, bx, by, bz);
                    buffer[0] = bx;
                    buffer[1] = by;
                    buffer[2] = bz;
                });

                // Log volume
                systems[ns].logger->add_logger<float>("volume", {1}, [ns, &systems](float* buffer) {
                    float vol = simulation_box::npt::get_volume(systems[ns].engine);
                    buffer[0] = vol;
                });
            }

        }
        if(verbose) printf("\n");

        if(verbose) printf("Initial potential energy:");
        for(System& sys: systems) {
            sys.engine.compute(PotentialAndDerivMode);
            if(verbose) printf(" %.2f", sys.engine.potential);
        }
        if(verbose) printf("\n");

        // Optional pre-run minimization
        if(enable_min_arg.getValue()) {
            int it = min_max_iter_arg.getValue();
            double etol = min_energy_tol_arg.getValue();
            double ftol = min_force_tol_arg.getValue();
            double mstep = min_step_arg.getValue();
            if(verbose) printf("\nMINIMIZATION: starting...\n");
            for(System& sys: systems) {
                martini_run_minimization(sys.engine, it, etol, ftol, mstep, verbose);
                // Save a frame of the minimised state when no MD round follows to save it, so that /output/pos
                // exists even if duration is 0
                if(sys.round_num == n_round) {
                    sys.engine.compute(PotentialAndDerivMode);
                    sys.logger->collect_samples();
                }
            }
            if(verbose) printf("MINIMIZATION: done\n\n");
        }


        // Install signal handlers to dump state only when the simulation has really started.  This is intended to prevent
        // loss of buffered data and to present final statistics.  It is especially useful when being killed due to running 
        // out of time on a cluster.
        SignalHandlerHandler sigint_handler (SIGINT,  abort_like_handler);
        SignalHandlerHandler sigterm_handler(SIGTERM, abort_like_handler);

        // we need to run everyone until the next synchronization event
        // a little care is needed if we are multiplexing the events
        auto tstart = chrono::high_resolution_clock::now();
        while(systems[0].round_num < n_round && received_signal==NO_SIGNAL) {
            // An exception must not escape an OpenMP structured block: doing so calls std::terminate
            // and the message is destroyed, which is exactly how a mid-run failure reported nothing but
            // "Abort trap: 6". Trap per system, keep the first message, and rethrow once serial so the
            // top-level handler can print why the run died.
            std::string omp_error;
            bool omp_failed = false;
            #pragma omp parallel for schedule(static,1)
            for(int ns=0; ns<int(systems.size()); ++ns) try {
                System& sys = systems[ns];
                for(bool do_break=false; (!do_break) && (sys.round_num<n_round); ++sys.round_num) {
                    int nr = sys.round_num;

                    // The frame of a round is its state before any of the round's operations, so that a
                    // continuation from it applies each of them once.  Frames lie on the frame grid, and every run
                    // also logs the round it starts from, so that a continuation's first frame is always the frame
                    // it continues from, on the grid or off it (a step-count run's final capture)
                    bool on_frame_grid = !(nr%frame_interval);
                    bool do_print = on_frame_grid || sys.round_num == start_round;
                    if(do_print) {
                        sys.engine.compute(PotentialAndDerivMode);
                        sys.logger->collect_samples();

                        if(verbose) {
                            double display_elapsed = nr*double(dt*inner_step);
                            if(sys.martini_hybrid_progress) {
                                double rg = compute_protein_rg(sys);
                                auto prot_components = compute_hybrid_protein_potential_components(sys.engine);
                                double protein_potential = prot_components.protein_total();
                                printf("step %d / %llu time %.3g system %2i temperature %.2f hbonds %5.1f",
                                       nr,
                                       static_cast<unsigned long long>(n_round),
                                       display_elapsed,
                                       ns,
                                       sys.temperature,
                                       get_n_hbond(sys.engine));
                                if(rg >= 0.0) printf(" %5.1f Rg,", rg);
                                else          printf("  N/A Rg,");
                                printf(" protein_potential % .2f total_potential % .2f\n",
                                       protein_potential,
                                       sys.engine.potential);
                            } else {
                                double rg = compute_all_atom_rg(sys);
                                printf("%*.0f / %*.0f elapsed %2i system %.2f temp %5.1f hbonds, Rg %5.1f A, potential % 8.2f\n",
                                       duration_print_width, display_elapsed,
                                       duration_print_width, display_duration_total,
                                       ns, sys.temperature,
                                       get_n_hbond(sys.engine),
                                       rg, sys.engine.potential);
                            }
                        }
                        fflush(stdout);
                    }

                    if(on_frame_grid && do_recenter) recenter(sys.engine.pos->output, xy_recenter_only, sys.n_atom);

                    // Don't pivot at t=0 so that a partially strained system may relax before the
                    // first pivot
                    if(nr && mc_interval && !(nr%mc_interval)) 
                        sys.mc_samplers.execute(sys.random_seed, nr, sys.temperature, sys.engine,
                                sys.mom, mc_half_kick, mc_slow_half_kick);

                    bool thermostat_round = !(nr%thermostat_interval);
                    if(thermostat_round) {
                        // Handle simulated annealing if applicable
                        if(anneal_factor != 1.)
                            sys.set_temperature(anneal_temp(sys.initial_temperature, inner_step*dt*(sys.round_num+1)));
                        martini_brownian::set_brownian_temperature(&sys.engine, sys.temperature);
                        // g-JF momenta are at the coordinates; the leapfrog's trail them by half a kick, so its
                        // integration cycle takes the thermostat step where it has the forces
                        if(brownian) sys.thermostat.apply(sys.mom, sys.n_atom, &sys.engine);
                    }

                    // Enforce fixed-in-space constraints before integration so fixed atoms
                    // do not receive a first-step kick from initialized/thermostatted momenta.
                    martini_fix_rigid::apply_fix_rigid_md(sys.engine, sys.engine.pos->output, sys.mom);

                    // MARTINI configurations list every dynamic particle in the single-step g-JF integrator
                    auto thermostat = thermostat_round ? &sys.thermostat : nullptr;
                    if(brownian) {
                        sys.engine.compute(DerivMode);
                        martini_brownian::apply_langevin_step(&sys.engine, sys.mom, dt);
                    } else if(integrator_arg.getValue() == "mv") {
                        sys.engine.integration_cycle(sys.mom, dt, inner_step, thermostat);
                    } else {
                        sys.engine.integration_cycle(sys.mom, dt, thermostat);
                    }

                    // Apply fix rigid constraints after integration
                    martini_fix_rigid::apply_fix_rigid_md(sys.engine, sys.engine.pos->output, sys.mom);

                    // Apply NPT barostat if enabled (skip step 0 to allow initial relaxation)
                    if(nr > 0)
                        simulation_box::npt::maybe_apply_barostat(sys.engine, sys.n_atom, nr, sys.temperature);

                    // Re-apply fixed-in-space constraints after barostat/box updates.
                    martini_fix_rigid::apply_fix_rigid_md(sys.engine, sys.engine.pos->output, sys.mom);
                    martini_hybrid::refresh_transition_holds_for_engine(sys.engine, nr + 1,
                                                                        restart_using_momentum_arg.getValue());

                    if(curvature_changer_interval && !(sys.round_num % curvature_changer_interval))
                        curvature_changer->attempt_change(sys.random_seed, sys.round_num, ns, sys, relative_curvature_radius_change);

                    // A block of rounds ends at the next replica exchange, and after at most 8 rounds so that a stop
                    // signal or the time limit, tested between blocks, stops every system at the same round
                    do_break = !((nr+1)%8) || (replica_interval && !((nr+1)%replica_interval));
                }
            } catch(const string &e) {
                #pragma omp critical
                if(!omp_failed) { omp_failed = true; omp_error = e; }
            } catch(const std::exception &e) {
                #pragma omp critical
                if(!omp_failed) { omp_failed = true; omp_error = string("std::exception: ") + e.what(); }
            } catch(...) {
                #pragma omp critical
                if(!omp_failed) { omp_failed = true; omp_error = "unknown exception"; }
            }
            // Here we are running in serial again
            if(omp_failed)
                throw string("integration failed at round ") + to_string(systems[0].round_num) +
                      ": " + omp_error;
            if(replica_interval && !(systems[0].round_num % replica_interval))
                replex->attempt_swaps(base_random_seed, systems[0].round_num,
                                      systems[0].round_num/replica_interval - 1, systems);

            if(time_lim > 0. && chrono::duration<double>(chrono::high_resolution_clock::now() - tstart).count() > time_lim) {
                passed_time_lim = true;
                break;
            }

        }

        if(use_duration_steps && n_round > 0 && !(enable_min_arg.getValue() && systems[0].round_num == start_round)) {
            // MARTINI step-count runs may end off a frame boundary; capture the final frame.  A run of no rounds
            // after a minimisation has logged this state already.
            for(System& sys: systems) {
                sys.engine.compute(PotentialAndDerivMode);
                sys.logger->collect_samples();
            }
        }


        if(received_signal!=NO_SIGNAL) {fprintf(stderr, "Received early termination signal\n");}
        if(passed_time_lim) {fprintf(stderr, "Passed time limit\n");}
        for(auto& sys: systems) sys.logger = shared_ptr<H5Logger>(); // release shared_ptr, which also flushes data during destructor

        auto elapsed = chrono::duration<double>(std::chrono::high_resolution_clock::now() - tstart).count();
        if(verbose)
            printf("\n\nfinished in %.1f seconds (%.2f us/systems/step, %.1e simulation_time_unit/hour)\n",
                elapsed,
                elapsed*1e6/systems.size()/(systems[0].round_num-start_round)/inner_step, 
                (systems[0].round_num-start_round)*inner_step*dt/elapsed * 3600.); 

        if(verbose) printf("\navg_kinetic_energy/1.5kT");
        for(auto& sys: systems) {
            double sum_kinetic = 0.;
            long n_kinetic = 0l;

            // each frame against its own temperature, which annealing changes during the run
            hid_t out_file = user_defined_output ? sys.output.get() : sys.config.get();
            vector<float> frame_temperature;
            traverse_dset<2,float>(out_file, "/output/temperature", [&](size_t nf, size_t ns, float x){
                frame_temperature.push_back(x);});
            size_t tot_frames = frame_temperature.size();

            traverse_dset<2,float>(out_file, "/output/kinetic", [&](size_t nf, size_t ns, float x){
                if(nf>=tot_frames/2){ sum_kinetic+=x/(1.5*frame_temperature[nf]); n_kinetic++; }
            });

            if(verbose) printf(" % .3f", sum_kinetic/n_kinetic);
        }
        if(verbose) printf("\n");

        // FIXME this code should be moved into MC sampler code
        try {
            if(mc_interval) {
                if(verbose)printf("pivot_success:\n");
                for(auto& sys: systems) {
                    std::vector<int64_t> ps(2,0);

                    if (user_defined_output) 
                        traverse_dset<2,int>(sys.output.get(), "/output/pivot_stats", [&](size_t nf, int d, int x) {ps[d] += x;});
                    else
                        traverse_dset<2,int>(sys.config.get(), "/output/pivot_stats", [&](size_t nf, int d, int x) {ps[d] += x;});

                    if(verbose)printf(" % .4f", double(ps[0])/double(ps[1]));
                }
                if(verbose)printf("\n");
            }
        } catch(...) {}  // stats reporting is optional

        try {
            if(mc_interval) {
                if(verbose)printf("jump_success:\n");
                for(auto& sys: systems) {
                    std::vector<int64_t> ps(2,0);

                    if (user_defined_output) 
                        traverse_dset<2,int>(sys.output.get(), "/output/jump_stats", [&](size_t nf, int d, int x) {ps[d] += x;});
                    else 
                        traverse_dset<2,int>(sys.config.get(), "/output/jump_stats", [&](size_t nf, int d, int x) {ps[d] += x;});

                    if(verbose)printf(" % .4f", double(ps[0])/double(ps[1]));
                }
                if(verbose)printf("\n");
            }
        } catch(...) {}  // stats reporting is optional

#ifdef COLLECT_PROFILE
        if(verbose) {
            printf("\n");
            global_time_keeper.print_report(inner_step*systems[0].round_num+1); //  FIXME inner_step
            printf("\n");
        }
#endif
    } catch(const string &e) {
        fprintf(stderr, "\n\nERROR: %s\n", e.c_str());
        return 1;
    } catch(...) {
        fprintf(stderr, "\n\nERROR: unknown error\n");
        return 1;
    }

    // By this point in the program, the signal handlers have been returned to the handlers installed
    // by the caller if Upside is running as a shared library function.  When we re-raise the signal,
    // the caller's handler will be able to take the signal.
    if(raise_signal_on_exit_if_received_arg.getValue() && received_signal!=NO_SIGNAL)
        raise(received_signal);

    return 0;
} catch(const TCLAP::ArgException &e) { 
    fprintf(stderr, "\n\nERROR: %s for argument %s\n", e.error().c_str(), e.argId().c_str());
    return 1;
} catch(const string &e) {
    fprintf(stderr, "\n\nERROR: %s\n", e.c_str());
    return 1;
}

int main(int argc, const char* const * argv) {
    return upside_main(argc, argv);
}
