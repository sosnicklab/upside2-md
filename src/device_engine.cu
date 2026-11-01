#include "device_engine.h"
#include "thermostat.h"
#include "random.h"
#include "martini.h"
#include "martini_brownian.h"
#include <typeinfo>

using namespace std;

DeviceCoord::DeviceCoord(int elem_width_, int row_width_, const vector<int>& n_elem_of_system):
    elem_width(elem_width_), row_width(row_width_), n_elem(0)
{
    for(int n: n_elem_of_system) {offset.push_back(n_elem); n_elem += n;}
    offset.push_back(n_elem);
    output = DeviceArray<float>(size_t(n_elem)*row_width);
    output.zero();
    sens = DeviceAccum(size_t(n_elem)*row_width);
}

vector<int> DeviceCoord::system_of_element() const {
    vector<int> ret;
    for(int ns=0; ns+1<int(offset.size()); ++ns)
        for(int i=offset[ns]; i<offset[ns+1]; ++i) ret.push_back(ns);
    return ret;
}


DeviceGraph::DeviceGraph(const vector<DerivEngine*>& engines_):
    n_system(engines_.size()), engines(engines_), potential(engines_.size()*engines_[0]->nodes.size()), error(1)
{
    error.zero();
    if(!n_system) throw string("the device engine needs at least one system");
    auto& e0 = *engines[0];
    int n_node = e0.nodes.size();
    for(auto e: engines) {
        if(int(e->nodes.size()) != n_node)
            throw string("the device engine runs systems of one potential graph, but the systems' node counts differ");
        for(int i=0; i<n_node; ++i) {
            const DerivComputation &c = *e->nodes[i].computation, &c0 = *e0.nodes[i].computation;
            if(e->nodes[i].name != e0.nodes[i].name || typeid(c) != typeid(c0))
                throw string("the device engine runs systems of one potential graph, but node ") + e0.nodes[i].name +
                    " differs between systems";
        }
        for(int i=0; i<n_node; ++i) node_index[e->nodes[i].computation.get()] = i;
    }

    coords.resize(n_node);
    nodes .resize(n_node);
    for(int i=0; i<n_node; ++i) {
        if(e0.nodes[i].computation->potential_term) continue;
        vector<int> n_elem;
        for(auto e: engines) {
            auto& c = static_cast<CoordNode&>(*e->nodes[i].computation);
            if(c.elem_width != static_cast<CoordNode&>(*e0.nodes[i].computation).elem_width)
                throw string("node ") + e0.nodes[i].name + " has a different width in different systems";
            n_elem.push_back(c.n_elem);
        }
        auto& c0 = static_cast<CoordNode&>(*e0.nodes[i].computation);
        coords[i].reset(new DeviceCoord(c0.elem_width, c0.output.row_width, n_elem));
    }

    string missing;
    for(int i=1; i<n_node; ++i) {  // node 0 is pos, which the engine itself writes
        vector<DerivComputation*> systems;
        for(auto e: engines) systems.push_back(e->nodes[i].computation.get());
        nodes[i].reset(systems[0]->make_device(*this, systems));
        if(!nodes[i]) missing += " " + e0.nodes[i].name;
    }
    if(missing.size()) throw string("these nodes have no device implementation:") + missing;
}

DeviceCoord& DeviceGraph::coord(const CoordNode& cpu_node) {
    auto it = node_index.find(&cpu_node);
    if(it == node_index.end() || !coords[it->second])
        throw string("the device engine has no coordinate node for this input");
    return *coords[it->second];
}

AccumView DeviceGraph::potential_of(const DerivComputation& cpu_node) {
    auto it = node_index.find(&cpu_node);
    if(it == node_index.end() || !cpu_node.potential_term)
        throw string("the device engine has no potential node for this input");
    auto v = potential.view();
    return AccumView{v.sum + it->second*n_system, v.bad + it->second*n_system};
}

void DeviceGraph::check_errors() const {
    int code = error.download()[0];
    if(code) throw checks.at(code-1);
}

vector<float> DeviceGraph::get_potential() const {
    auto p = potential.download();
    vector<float> ret(n_system, 0.f);
    for(int i: engines[0]->germ_exec_levels)
        if(engines[0]->nodes[i].computation->potential_term)
            for(int ns=0; ns<n_system; ++ns) ret[ns] += p[i*n_system+ns];
    return ret;
}

vector<float> DeviceGraph::get_potential(int node) const {
    auto p = potential.download();
    return vector<float>(p.begin()+node*n_system, p.begin()+(node+1)*n_system);
}

void DeviceGraph::run(ComputeMode mode, int integrator_level, bool all_levels) {
    if(mode == PotentialAndDerivMode) potential.zero();
    auto& e0 = *engines[0];
    auto at_level = [&](int i) {
        int lvl = e0.nodes[i].integrator_level;
        return all_levels || lvl == integrator_level || lvl == -1;};

    for(int i: e0.germ_exec_levels) {
        if(coords[i]) coords[i]->sens.zero();
        if(nodes[i] && at_level(i)) nodes[i]->compute_value(mode);
    }
    for(int i: e0.deriv_exec_levels)
        if(nodes[i] && at_level(i)) nodes[i]->propagate_deriv();
}

void DeviceGraph::compute(ComputeMode mode)                       {run(mode, 0, true);}
void DeviceGraph::compute(ComputeMode mode, int integrator_level) {run(mode, integrator_level, false);}


DeviceDynamics::DeviceDynamics(const vector<DerivEngine*>& engines):
    graph(engines), brownian(martini_brownian::has_brownian(engines[0])), n_inner_steps(1), brownian_time_step(0.f)
{
    auto& pos = *graph.coords[0];
    int n = pos.n_elem, nsys = graph.n_system;
    mom = DeviceArray<float>(pos.output.n);
    lag = DeviceArray<float>(pos.output.n);
    mom.zero();
    lag.zero();
    vector<int> local;
    for(int ns=0; ns<nsys; ++ns)
        for(int na=0; na<pos.offset[ns+1]-pos.offset[ns]; ++na) local.push_back(na);
    atom_system = DeviceArray<int>(pos.system_of_element());
    atom_local  = DeviceArray<int>(local);
    thermostat_seed        = DeviceArray<uint32_t>(nsys);
    thermostat_invocation  = DeviceArray<unsigned long long>(nsys);
    thermostat_mom_scale   = DeviceArray<float>(nsys);
    thermostat_noise_scale = DeviceArray<float>(nsys);
    thermostat_temp        = DeviceArray<float>(nsys);
    thermostat_delta_t     = DeviceArray<float>(nsys);

    vector<float> m;
    for(auto e: engines) {
        auto me = martini_masses::get_masses(e);
        if(martini_masses::has_masses(e) != martini_masses::has_masses(engines[0]))
            throw string("the device engine needs masses in every system or in none");
        m.insert(m.end(), me.begin(), me.end());
    }
    if(m.size()) mass = DeviceArray<float>(m);
    fixed   = DeviceArray<unsigned char>(n);
    z_fixed = DeviceArray<unsigned char>(n);
    vector<vector<int>> fa, za;
    for(auto e: engines) {
        if(martini_fix_rigid::get_rigid_groups(*e).size())
            throw string("the device engine has no rigid groups (the hybrid's preparation stages)");
        fa.push_back(martini_fix_rigid::get_fixed_atoms(*e));
        za.push_back(martini_fix_rigid::get_z_fixed_atoms(*e));
    }
    set_fixed_atoms(fa, za);

    vector<unsigned char> stoch(n, 0u), lst(n, 0u);
    vector<float> fr(n, 0.f), ref;
    vector<uint32_t> seed;
    if(brownian) {
        if(!m.size()) throw string("the g-JF integrator needs masses");
        for(int ns=0; ns<nsys; ++ns) {
            if(martini_brownian::has_brownian(engines[ns]) != brownian)
                throw string("the device engine needs the g-JF integrator in every system or in none");
            auto st = martini_brownian::get_settings(engines[ns]);
            if(ns == 0) {n_inner_steps = st.n_inner_steps; brownian_time_step = st.time_step;}
            if(st.n_inner_steps != n_inner_steps || st.time_step != brownian_time_step)
                throw string("the device engine needs one g-JF step and substep count in every system");
            for(size_t i=0; i<st.atom_index.size(); ++i) {
                int a = pos.offset[ns] + st.atom_index[i];
                lst[a] = 1u;
                fr[a] = st.friction[i];
                if(st.friction[i] > 0.f) stoch[a] = 1u;
            }
            ref.push_back(st.reference_temperature);
            seed.push_back(st.random_seed);
        }
        brownian_seed = DeviceArray<uint32_t>(seed);
        reference_temperature = DeviceArray<float>(ref);
        brownian_temperature = DeviceArray<float>(nsys);
        brownian_invocation = DeviceArray<unsigned long long>(nsys);
        force_n = DeviceArray<float>(pos.output.n);
    }
    stochastic = DeviceArray<unsigned char>(stoch);
    listed     = DeviceArray<unsigned char>(lst);
    friction   = DeviceArray<float>(fr);
}

void DeviceDynamics::set_fixed_atoms(const vector<vector<int>>& fixed_, const vector<vector<int>>& z_fixed_) {
    if(fixed_ == fixed_atoms && z_fixed_ == z_fixed_atoms) return;
    fixed_atoms = fixed_;
    z_fixed_atoms = z_fixed_;
    auto& pos = *graph.coords[0];
    vector<unsigned char> f(pos.n_elem, 0u), z(pos.n_elem, 0u);
    for(int ns=0; ns<graph.n_system; ++ns) {
        for(int a: fixed_atoms  [ns]) f[pos.offset[ns]+a] = 1u;
        for(int a: z_fixed_atoms[ns]) z[pos.offset[ns]+a] = 1u;
    }
    fixed.upload(f);
    z_fixed.upload(z);
}

void DeviceDynamics::set_state(const vector<VecArray>& pos, const vector<VecArray>& mom_) {
    auto& c = *graph.coords[0];
    vector<float> x(c.output.n, 0.f), p(c.output.n, 0.f);
    for(int ns=0; ns<graph.n_system; ++ns)
        for(int na=0; na<c.offset[ns+1]-c.offset[ns]; ++na)
            for(int d=0; d<3; ++d) {
                x[(c.offset[ns]+na)*c.row_width+d] = pos [ns](d,na);
                p[(c.offset[ns]+na)*c.row_width+d] = mom_[ns](d,na);
            }
    c.output.upload(x);
    mom.upload(p);
}

void DeviceDynamics::get_state(const vector<VecArray>& pos, const vector<VecArray>& mom_) const {
    auto& c = *graph.coords[0];
    auto x = c.output.download();
    auto p = mom.download();
    for(int ns=0; ns<graph.n_system; ++ns)
        for(int na=0; na<c.offset[ns+1]-c.offset[ns]; ++na)
            for(int d=0; d<3; ++d) {
                VecArray xs = pos[ns], ps = mom_[ns];
                xs(d,na) = x[(c.offset[ns]+na)*c.row_width+d];
                ps(d,na) = p[(c.offset[ns]+na)*c.row_width+d];
            }
}

// A kick of the momentum by -kick dV/dx and, unless drift is 0, a drift of the coordinates by drift p/m, as
// deriv_engine.cpp's leapfrog_step: fixed atoms keep zero momentum and z-fixed atoms none along z
void DeviceDynamics::leapfrog_step(float kick, float drift) {
    auto& c = *graph.coords[0];
    float* x = c.output.ptr;
    float* p = mom.ptr;
    auto sens = c.sens.view();
    int w = c.row_width;
    const float* m = mass.n ? mass.ptr : nullptr;
    const unsigned char *fx = fixed.ptr, *zf = z_fixed.ptr;
    parallel_for(c.n_elem, DEVICE_LAMBDA (int na) {
        if(fx[na]) {for(int d=0; d<3; ++d) p[na*w+d] = 0.f; return;}
        bool z_fixed = zf[na];
        float scale = m ? drift/m[na] : drift;
        for(int d=0; d<3; ++d) {
            bool held = z_fixed && d == 2;
            float pd = held ? 0.f : p[na*w+d] - kick*sens.get(na*w+d);
            p[na*w+d] = pd;
            if(drift != 0.f && !held) x[na*w+d] += scale*pd;
        }
    });
}

// OrnsteinUhlenbeckThermostat::apply, which skips the fixed atoms.  at_coordinates: the leapfrog's step on the
// momentum at the coordinates, mom - lag (deriv_engine.cpp's thermostat_at_coordinates); otherwise the g-JF momenta,
// skipping the g-JF bath's atoms too.  An atom's own timescale and its mass enter as in apply.
void DeviceDynamics::thermostat_step(const vector<OrnsteinUhlenbeckThermostat*>& thermostats, bool at_coordinates) {
    if(int(thermostats.size()) != graph.n_system) throw string("the device thermostat needs one thermostat per system");
    auto& c = *graph.coords[0];
    vector<uint32_t> seed;
    vector<unsigned long long> n_invocations;
    vector<float> mom_scale, noise_scale, temp, delta_t, tau;
    // apply takes an atom's own timescale where the thermostat holds one per atom
    auto own = [&](int ns) {return int(thermostats[ns]->atom_timescale.size()) == c.offset[ns+1]-c.offset[ns];};
    bool own_timescales = own(0);
    for(int ns=0; ns<graph.n_system; ++ns) {
        auto t = thermostats[ns];
        if(own(ns) != own_timescales)
            throw string("the device thermostat needs per-atom timescales in every system or in none");
        seed.push_back(t->random_seed);
        n_invocations.push_back(t->get_n_invocations());
        mom_scale.push_back(t->mom_scale);
        noise_scale.push_back(t->noise_scale);
        temp.push_back(t->temp);
        delta_t.push_back(t->delta_t);
        if(own_timescales) tau.insert(tau.end(), t->atom_timescale.begin(), t->atom_timescale.end());
    }
    thermostat_seed.upload(seed);
    thermostat_invocation.upload(n_invocations);
    thermostat_mom_scale.upload(mom_scale);
    thermostat_noise_scale.upload(noise_scale);
    thermostat_temp.upload(temp);
    thermostat_delta_t.upload(delta_t);
    // read from the configs, so they hold for the run
    if(own_timescales && atom_timescale.n != tau.size()) atom_timescale = DeviceArray<float>(tau);

    float* p = mom.ptr;
    const float* l = lag.ptr;
    const int* sys = atom_system.ptr;
    const int* local = atom_local.ptr;
    const uint32_t* sd = thermostat_seed.ptr;
    const unsigned long long* inv = thermostat_invocation.ptr;
    const float *ms = thermostat_mom_scale.ptr, *nsc = thermostat_noise_scale.ptr;
    const float *tp = thermostat_temp.ptr, *dt = thermostat_delta_t.ptr;
    const float* tau_d = own_timescales ? atom_timescale.ptr : nullptr;
    const float* mass_d = mass.n ? mass.ptr : nullptr;
    const unsigned char *fx = fixed.ptr, *st = stochastic.ptr;
    bool with_lag = at_coordinates;
    int w = c.row_width;
    parallel_for(c.n_elem, DEVICE_LAMBDA (int na) {
        if(fx[na] || (!with_lag && st[na])) return;
        int ns = sys[na];
        RandomGenerator random(sd[ns], THERMOSTAT_RANDOM_STREAM, local[na], inv[ns]);
        float atom_mom_scale = ms[ns], atom_noise_scale = nsc[ns];
        if(tau_d) {
            atom_mom_scale = expf(-dt[ns]/tau_d[na]);
            atom_noise_scale = sqrtf(tp[ns] * (1.f-atom_mom_scale*atom_mom_scale));
        }
        if(mass_d) atom_noise_scale *= sqrtf(mass_d[na]);
        auto noise = random.normal3();
        for(int d=0; d<3; ++d) {
            float pc = with_lag ? p[na*w+d] - l[na*w+d] : p[na*w+d];
            float pn = atom_mom_scale*pc + atom_noise_scale*noise[d];
            p[na*w+d] = with_lag ? pn + l[na*w+d] : pn;
        }
    });
    for(auto t: thermostats) t->set_n_invocations(t->get_n_invocations()+1);
}

void DeviceDynamics::zero_fixed_momenta() {
    auto& c = *graph.coords[0];
    float* p = mom.ptr;
    const unsigned char *fx = fixed.ptr, *zf = z_fixed.ptr;
    int w = c.row_width;
    parallel_for(c.n_elem, DEVICE_LAMBDA (int na) {
        if(fx[na]) for(int d=0; d<3; ++d) p[na*w+d] = 0.f;
        if(zf[na]) p[na*w+2] = 0.f;
    });
}

void DeviceDynamics::langevin_cycle(float dt, const vector<float>& temperature, uint64_t n_invocations) {
    if(!brownian) throw string("the systems have no g-JF integrator");
    if(fabsf(dt - brownian_time_step) > 1e-6f*fmaxf(dt, brownian_time_step))
        throw string("runtime timestep does not match /input/brownian numerical_time_step");
    brownian_temperature.upload(temperature);
    brownian_invocation.upload(vector<unsigned long long>(graph.n_system, n_invocations));

    auto& c = *graph.coords[0];
    float *x = c.output.ptr, *p = mom.ptr, *fn = force_n.ptr;
    auto sens = c.sens.view();
    const int *sys = atom_system.ptr, *local = atom_local.ptr;
    const unsigned char *lst = listed.ptr, *fx = fixed.ptr, *zf = z_fixed.ptr;
    const float *fr = friction.ptr, *ms = mass.ptr, *kt = brownian_temperature.ptr, *ref = reference_temperature.ptr;
    const uint32_t* sd = brownian_seed.ptr;
    const unsigned long long* inv = brownian_invocation.ptr;
    int w = c.row_width, N = n_inner_steps;
    const float dt_i = dt / static_cast<float>(N);

    graph.compute(DerivMode);
    for(int inner = 0; inner < N; inner++) {
        parallel_for(c.n_elem, DEVICE_LAMBDA (int na) {
            if(!lst[na]) return;
            if(fx[na]) {for(int d=0; d<3; ++d) p[na*w+d] = 0.f; return;}
            int ns = sys[na];
            const float kT = kt[ns];
            // Hold the calibrated diffusion, not the calibrated drag: D = kT/alpha, so alpha tracks the running
            // temperature
            const float alpha = fr[na]*(kT/ref[ns]);
            const float m = ms[na];
            const float b = 1.f/(1.f + alpha*dt_i/(2.f*m));
            const float bnoise = sqrtf(2.f*alpha*kT*dt_i);
            RandomGenerator random(sd[ns], BROWNIAN_RANDOM_STREAM, uint32_t(local[na]),
                                   inv[ns]*uint64_t(N) + uint64_t(inner));
            auto beta = bnoise*random.normal3();
            auto f = make_vec3(-sens.get(na*w+0), -sens.get(na*w+1), -sens.get(na*w+2));
            auto pv = load_vec<3>(p+na*w);
            if(zf[na]) {beta.z() = 0.f; f.z() = 0.f; pv.z() = 0.f;}
            // x + b dt/m (p + dt/2 f + beta/2), one rounding at the coordinate, as apply_langevin_step (findings 11g)
            auto xv = load_vec<3>(x+na*w);
            store_vec(x+na*w, xv + (b*dt_i/m)*(pv + (0.5f*dt_i)*f + 0.5f*beta));
            store_vec(fn+na*w, f);
        });
        graph.compute(DerivMode);
        parallel_for(c.n_elem, DEVICE_LAMBDA (int na) {
            if(!lst[na] || fx[na]) return;
            int ns = sys[na];
            const float kT = kt[ns];
            const float alpha = fr[na]*(kT/ref[ns]);
            const float m = ms[na];
            const float half = alpha*dt_i/(2.f*m);
            const float b = 1.f/(1.f + half);
            const float a = (1.f - half)/(1.f + half);
            const float bnoise = sqrtf(2.f*alpha*kT*dt_i);
            RandomGenerator random(sd[ns], BROWNIAN_RANDOM_STREAM, uint32_t(local[na]),
                                   inv[ns]*uint64_t(N) + uint64_t(inner));
            auto beta = bnoise*random.normal3();
            auto f_n  = load_vec<3>(fn+na*w);
            auto f_n1 = make_vec3(-sens.get(na*w+0), -sens.get(na*w+1), -sens.get(na*w+2));
            auto pv = load_vec<3>(p+na*w);
            if(zf[na]) {beta.z() = 0.f; f_n1.z() = 0.f; pv.z() = 0.f;}
            store_vec(p+na*w, a*pv + (dt_i*0.5f)*(a*f_n + f_n1) + b*beta);
        });
    }
    graph.check_errors();
    for(auto e: graph.engines) martini_brownian::set_brownian_invocations(e, n_invocations+1);
}

void DeviceDynamics::integration_cycle(float dt, const vector<OrnsteinUhlenbeckThermostat*>& thermostats) {
    auto& c = *graph.coords[0];
    float* l = lag.ptr;
    auto sens = c.sens.view();
    int w = c.row_width;
    for(int stage=0; stage<3; ++stage) {
        graph.compute(DerivMode);
        if(!stage && thermostats.size()) {
            float scale = 0.5f*dt;
            parallel_for(c.n_elem, DEVICE_LAMBDA (int na) {
                for(int d=0; d<3; ++d) l[na*w+d] = scale*sens.get(na*w+d);});
            thermostat_step(thermostats, true);
        }
        leapfrog_step(dt, dt);
    }
    graph.check_errors();
}

void DeviceDynamics::integration_cycle(float dt, int inner_step,
        const vector<OrnsteinUhlenbeckThermostat*>& thermostats) {
    auto& c = *graph.coords[0];
    float* l = lag.ptr;
    auto sens = c.sens.view();
    int w = c.row_width;

    graph.compute(DerivMode, 1);
    // after the slow kick the stored momentum leads the one at the coordinates by half a slow kick and trails it
    // by half a fast one
    if(thermostats.size()) {
        float scale = -0.5f*inner_step*dt;
        parallel_for(c.n_elem, DEVICE_LAMBDA (int na) {
            for(int d=0; d<3; ++d) l[na*w+d] = scale*sens.get(na*w+d);});
    }
    leapfrog_step(inner_step*dt, 0.f);

    for(int i=0; i<inner_step; ++i) {
        graph.compute(DerivMode, 0);
        if(!i && thermostats.size()) {
            float scale = 0.5f*dt;
            parallel_for(c.n_elem, DEVICE_LAMBDA (int na) {
                for(int d=0; d<3; ++d) l[na*w+d] += scale*sens.get(na*w+d);});
            thermostat_step(thermostats, true);
        }
        leapfrog_step(dt, dt);
    }
    graph.check_errors();
}
