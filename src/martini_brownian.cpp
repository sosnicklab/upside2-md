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
    vector<int>   atom_index;
    vector<unsigned char> stochastic_mask;
    float    temperature = 1.f;
    float    time_step = 0.f;
    vector<float> friction;
    uint32_t random_seed   = 0;
    uint64_t n_invocations = 0;
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

    auto rt = unique_ptr<BrownianRuntime>(new BrownianRuntime());
    rt->atom_index = atom_index;
    rt->time_step = read_attribute<float>(config_root, "/input/brownian", "numerical_time_step");
    if(!(rt->time_step > 0.f))
        throw string("/input/brownian numerical_time_step must be positive");

    if(!h5_exists(grp, "friction"))
        throw string("/input/brownian/friction is required");
    vector<hsize_t> friction_size = get_dset_size(1, grp, "friction");
    if(friction_size[0] != atom_index.size())
        throw string("/input/brownian friction must match atom_index length");
    rt->friction.resize(atom_index.size());
    traverse_dset<1, float>(grp, "friction", [&](size_t i, float value) {
        if(value < 0.f)
            throw string("/input/brownian friction must be nonnegative");
        rt->friction[i] = value;
    });
    rt->random_seed = random_seed;
    int n_atom = engine->pos->n_atom;
    rt->stochastic_mask.assign(static_cast<size_t>(n_atom), 0);
    for(size_t i=0; i<atom_index.size(); ++i) {
        int idx = atom_index[i];
        if(idx < 0 || idx >= n_atom)
            throw string("/input/brownian atom_index out of range");
        if(rt->friction[i] > 0.f) rt->stochastic_mask[static_cast<size_t>(idx)] = 1;
    }
    runtime_map()[engine] = std::move(rt);
}

bool has_brownian(DerivEngine* engine) {
    return runtime_map().count(engine) > 0;
}

const vector<unsigned char>& stochastic_mask(DerivEngine* engine) {
    static const vector<unsigned char> empty;
    auto it = runtime_map().find(engine);
    return it == runtime_map().end() ? empty : it->second->stochastic_mask;
}

void set_brownian_temperature(DerivEngine* engine, float temperature) {
    auto it = runtime_map().find(engine);
    if(it != runtime_map().end()) it->second->temperature = temperature;
}

void apply_langevin_step(DerivEngine* engine, VecArray mom, float dt) {
    auto it = runtime_map().find(engine);
    if(it == runtime_map().end()) return;
    BrownianRuntime& rt = *it->second;
    Timer timer(string("brownian"));
    if(fabsf(dt - rt.time_step) > 1e-6f*fmaxf(dt, rt.time_step))
        throw string("runtime timestep does not match /input/brownian numerical_time_step");

    const float kT = rt.temperature;
    const vector<int>& idx = rt.atom_index;
    size_t nb = idx.size();
    VecArray pos_out = engine->pos->output;
    VecArray sens    = engine->pos->sens;
    const bool have_mass = martini_masses::has_masses(engine);
    vector<unsigned char> fixed_mask(static_cast<size_t>(engine->pos->n_atom), 0u);
    vector<unsigned char> z_fixed_mask(static_cast<size_t>(engine->pos->n_atom), 0u);
    for(int atom : martini_fix_rigid::get_fixed_atoms(*engine)) {
        if(atom >= 0 && atom < engine->pos->n_atom) fixed_mask[static_cast<size_t>(atom)] = 1u;
    }
    for(int atom : martini_fix_rigid::get_z_fixed_atoms(*engine)) {
        if(atom >= 0 && atom < engine->pos->n_atom) z_fixed_mask[static_cast<size_t>(atom)] = 1u;
    }
    vector<float> fn(3*nb);
    for(size_t i=0;i<nb;++i) {
        int atom = idx[i];
        if(fixed_mask[static_cast<size_t>(atom)]) {
            store_vec(mom, atom, make_zero<3>());
            continue;
        }
        const float alpha = rt.friction[i];
        const float m = have_mass ? martini_masses::get_mass(engine, atom) : 1.f;
        const float b = 1.f/(1.f + alpha*dt/(2.f*m));
        const float bnoise = sqrtf(2.f*alpha*kT*dt);
        RandomGenerator random(rt.random_seed, BROWNIAN_RANDOM_STREAM,
                               static_cast<uint32_t>(atom), rt.n_invocations);
        auto beta = bnoise*random.normal3();
        auto f = (-1.f)*load_vec<3>(sens, atom);
        auto p = load_vec<3>(mom, atom);
        if(z_fixed_mask[static_cast<size_t>(atom)]) {
            beta.z() = 0.f;
            f.z() = 0.f;
            p.z() = 0.f;
        }
        auto x = load_vec<3>(pos_out, atom);
        store_vec(pos_out, atom, x + (b*dt/m)*p + (b*dt*dt/(2.f*m))*f + (b*dt/(2.f*m))*beta);
        fn[3*i] = f.x();
        fn[3*i+1] = f.y();
        fn[3*i+2] = f.z();
    }
    engine->compute(DerivMode);
    for(size_t i=0;i<nb;++i) {
        int atom = idx[i];
        if(fixed_mask[static_cast<size_t>(atom)]) continue;
        const float alpha = rt.friction[i];
        const float m = have_mass ? martini_masses::get_mass(engine, atom) : 1.f;
        const float half = alpha*dt/(2.f*m);
        const float b = 1.f/(1.f + half);
        const float a = (1.f - half)/(1.f + half);
        const float bnoise = sqrtf(2.f*alpha*kT*dt);
        RandomGenerator random(rt.random_seed, BROWNIAN_RANDOM_STREAM,
                               static_cast<uint32_t>(atom), rt.n_invocations);
        auto beta = bnoise*random.normal3();
        auto f_n = make_vec3(fn[3*i], fn[3*i+1], fn[3*i+2]);
        auto f_n1 = (-1.f)*load_vec<3>(sens, atom);
        auto p = load_vec<3>(mom, atom);
        if(z_fixed_mask[static_cast<size_t>(atom)]) {
            beta.z() = 0.f;
            f_n1.z() = 0.f;
            p.z() = 0.f;
        }
        store_vec(mom, atom, a*p + (dt*0.5f)*(a*f_n + f_n1) + b*beta);
    }
    rt.n_invocations++;
}

}
