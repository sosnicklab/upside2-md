#include "martini_brownian.h"
#include "martini.h"
#include "random.h"
#include "timing.h"
#include <map>
#include <vector>
#include <memory>
#include <string>
#include <cmath>

using namespace std;
using namespace h5;

namespace martini_brownian {

static vector<int> read_int_dataset(hid_t grp, const char* name) {
    vector<hsize_t> sz = get_dset_size(1, grp, name);
    vector<int> out(sz[0]);
    traverse_dset<1, int32_t>(grp, name, [&](size_t i, int32_t v) { out[i] = v; });
    return out;
}

struct BrownianRuntime {
    vector<int>   atom_index;      // beads integrated by Brownian dynamics
    vector<unsigned char> mask;    // per-atom (n_atom): 1 if Brownian
    float    temperature   = 1.f;  // current kT in E_up
    float    gamma_ref     = 1.f;  // friction at t_ref (E_up * T_up / Ang^2)
    float    t_ref         = 1.f;  // reference kT in E_up
    float    activation_energy = 0.f; // Ea in E_up (DOPC lateral-diffusion activation energy)
    int      n_substep     = 1;    // RESPA inner sub-steps of the lipids per outer step
    int      langevin      = 0;    // 0 = overdamped BD + RESPA (default); 1 = g-JF inertial Langevin single-step
    vector<float> friction_scale;  // per-bead multiplier on gamma(T), aligned with atom_index
                                   // (=1 for lipids; >1 for the protein, whose slower diffusion
                                   //  keeps it embedded and whose LJ collisions need sub-stepping)
    uint32_t random_seed   = 0;
    uint64_t n_invocations = 0;

    // Arrhenius friction gamma(T) = gamma_ref*(kT/t_ref)*exp(Ea*(1/kT - 1/t_ref)),
    // calibrated once at t_ref so lipid diffusion is temperature-transferable.
    float gamma() const {
        return gamma_ref * (temperature/t_ref)
             * expf(activation_energy * (1.f/temperature - 1.f/t_ref));
    }
};

static map<DerivEngine*, unique_ptr<BrownianRuntime>>& runtime_map() {
    static map<DerivEngine*, unique_ptr<BrownianRuntime>> m;
    return m;
}

void register_brownian_for_engine(DerivEngine* engine, hid_t config_root, uint32_t random_seed) {
    if(!engine) return;
    runtime_map().erase(engine);
    if(!h5_exists(config_root, "/input/brownian")) return;

    H5Obj grp_obj = open_group(config_root, "/input/brownian");
    hid_t grp = grp_obj.get();

    vector<int> atom_index = read_int_dataset(grp, "atom_index");
    if(atom_index.empty())
        throw string("/input/brownian/atom_index must contain at least one atom");

    float gamma_ref = read_attribute<float>(config_root, "/input/brownian", "gamma_ref");
    float t_ref     = read_attribute<float>(config_root, "/input/brownian", "t_ref");
    float ea        = read_attribute<float>(config_root, "/input/brownian", "activation_energy", 0.f);
    if(gamma_ref <= 0.f) throw string("/input/brownian gamma_ref must be positive");
    if(t_ref     <= 0.f) throw string("/input/brownian t_ref must be positive");

    auto rt = unique_ptr<BrownianRuntime>(new BrownianRuntime());
    rt->atom_index        = atom_index;
    rt->gamma_ref         = gamma_ref;
    rt->t_ref             = t_ref;
    rt->activation_energy = ea;
    rt->temperature       = t_ref;   // until set_brownian_temperature is called
    rt->n_substep   = read_attribute<int>(config_root, "/input/brownian", "n_substep", 1);
    if(rt->n_substep < 1) rt->n_substep = 1;
    rt->langevin    = read_attribute<int>(config_root, "/input/brownian", "langevin", 0);

    // Optional per-bead friction multiplier on gamma (default 1 for every bead). Lipids use 1
    // (full g-JF friction); the protein uses 0, which reduces the g-JF step to plain
    // velocity-Verlet (see apply_brownian_step) so it stays on the same integrator as the lipids
    // and is thermostatted separately by the Ornstein-Uhlenbeck thermostat.
    rt->friction_scale.assign(atom_index.size(), 1.f);
    if(h5_exists(grp, "friction_scale")) {
        vector<hsize_t> sz = get_dset_size(1, grp, "friction_scale");
        if(sz[0] != atom_index.size())
            throw string("/input/brownian friction_scale must match atom_index length");
        traverse_dset<1, float>(grp, "friction_scale", [&](size_t i, float v) {
            rt->friction_scale[i] = (v >= 0.f) ? v : 1.f;
        });
    }
    rt->random_seed = random_seed;
    int n_atom = engine->pos->n_atom;
    rt->mask.assign(static_cast<size_t>(n_atom), 0);
    for(int idx : atom_index) {
        if(idx < 0 || idx >= n_atom)
            throw string("/input/brownian atom_index out of range");
        rt->mask[static_cast<size_t>(idx)] = 1;
    }
    runtime_map()[engine] = std::move(rt);
}

bool has_brownian(DerivEngine* engine) {
    return runtime_map().count(engine) > 0;
}

const vector<unsigned char>& brownian_mask(DerivEngine* engine) {
    static const vector<unsigned char> empty;
    auto it = runtime_map().find(engine);
    return it == runtime_map().end() ? empty : it->second->mask;
}

void set_brownian_temperature(DerivEngine* engine, float temperature) {
    auto it = runtime_map().find(engine);
    if(it != runtime_map().end()) it->second->temperature = temperature;
}

void apply_brownian_step(DerivEngine* engine, VecArray mom, float dt) {
    auto it = runtime_map().find(engine);
    if(it == runtime_map().end()) return;
    BrownianRuntime& rt = *it->second;
    Timer timer(string("brownian"));

    const float kT = rt.temperature;
    const vector<int>& idx = rt.atom_index;
    size_t nb = idx.size();
    VecArray pos_out = engine->pos->output;
    VecArray sens    = engine->pos->sens;
    float gamma = rt.gamma();              // Arrhenius friction at the current temperature

    if(rt.langevin) {
        // Inertial Langevin, Gronbech-Jensen & Farago (Mol. Phys. 2013), single step (no RESPA).
        // alpha = gamma sets D = kT/alpha, so the overdamped calibration carries over.
        const bool have_mass = martini_masses::has_masses(engine);
        vector<float> fn(3*nb);        // f_n stored across the recompute; beta is redrawn (RNG is deterministic)
        for(size_t i=0;i<nb;++i) {
            const float alpha = gamma * rt.friction_scale[i];
            const float m     = have_mass ? martini_masses::get_mass(engine, idx[i]) : 1.f;
            const float b     = 1.f/(1.f + alpha*dt/(2.f*m));
            const float bnoise= sqrtf(2.f*alpha*kT*dt);
            RandomGenerator random(rt.random_seed, BROWNIAN_RANDOM_STREAM,
                                   static_cast<uint32_t>(idx[i]), rt.n_invocations);
            auto beta = bnoise*random.normal3();
            auto f    = (-1.f)*load_vec<3>(sens, idx[i]);   // f_n = -dU/dx
            auto p    = load_vec<3>(mom, idx[i]);           // p_n
            auto x    = load_vec<3>(pos_out, idx[i]);
            store_vec(pos_out, idx[i], x + (b*dt/m)*p + (b*dt*dt/(2.f*m))*f + (b*dt/(2.f*m))*beta);
            fn[3*i]=f.x(); fn[3*i+1]=f.y(); fn[3*i+2]=f.z();
        }
        engine->compute(DerivMode);                          // f_{n+1} at the updated positions
        for(size_t i=0;i<nb;++i) {
            const float alpha = gamma * rt.friction_scale[i];
            const float m     = have_mass ? martini_masses::get_mass(engine, idx[i]) : 1.f;
            const float half  = alpha*dt/(2.f*m);
            const float b     = 1.f/(1.f + half);
            const float a     = (1.f - half)/(1.f + half);
            const float bnoise= sqrtf(2.f*alpha*kT*dt);
            RandomGenerator random(rt.random_seed, BROWNIAN_RANDOM_STREAM,
                                   static_cast<uint32_t>(idx[i]), rt.n_invocations);
            auto beta = bnoise*random.normal3();             // identical draw to pass 1
            auto f_n  = make_vec3(fn[3*i], fn[3*i+1], fn[3*i+2]);
            auto f_n1 = (-1.f)*load_vec<3>(sens, idx[i]);
            auto p    = load_vec<3>(mom, idx[i]);            // still p_n (untouched in pass 1)
            store_vec(mom, idx[i], a*p + (dt*0.5f)*(a*f_n + f_n1) + b*beta);
        }
        rt.n_invocations++;
        return;
    }

    // RESPA sub-stepping: M overdamped inner steps of dt/M per outer step, refreshing forces
    // between them, so the friction sets the correct lateral diffusion (both lipids and protein).
    const int   M         = rt.n_substep;
    const float dt_inner  = dt / static_cast<float>(M);
    const float mob_base   = dt_inner / gamma;                 // per-bead: divided by friction_scale
    const float noise_base = sqrtf(2.f * kT * dt_inner / gamma); // per-bead: divided by sqrt(scale)

    // Per-inner-step displacement cap: a numerical safeguard against rare deep-core overshoots;
    // with enough sub-steps it never fires on thermal motion, so it does not bias diffusion.
    const float dmax = 2.0f;
    for(int m=0; m<M; ++m) {
        // Refresh forces between inner steps (m=0 reuses the outer compute).
        if(m > 0) engine->compute(DerivMode);
        for(size_t i=0;i<nb;++i) {
            float sc  = rt.friction_scale[i];               // gamma_i = gamma * sc
            float mob = mob_base / sc;
            float noise_scale = noise_base / sqrtf(sc);
            auto xi = load_vec<3>(pos_out, idx[i]);
            auto g  = load_vec<3>(sens,    idx[i]);
            RandomGenerator random(rt.random_seed, BROWNIAN_RANDOM_STREAM,
                                   static_cast<uint32_t>(idx[i]), rt.n_invocations);
            auto step = (-mob)*g + noise_scale*random.normal3();
            float sn = sqrtf(mag2(step));
            if(sn > dmax) step = step*(dmax/sn);
            store_vec(pos_out, idx[i], xi + step);
        }
        rt.n_invocations++;                         // independent noise per inner step
    }

    // Overdamped beads carry no momentum.
    for(size_t i=0;i<nb;++i) store_vec(mom, idx[i], make_zero<3>());
}

}
