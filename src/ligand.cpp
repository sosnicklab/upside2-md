#include "ligand.h"
#include "affine.h"
#include "deriv_engine.h"
#include "timing.h"
#include <array>
#include <cmath>

// Reuse residue-local sidechain heavy-atom geometry construction for atomistic
// ligand/rotamer interactions without widening changes outside src/ligand.*.
#include "generate_from_rotamer.cpp"

using namespace h5;
using namespace std;
using namespace ligand;

namespace ligand {

vector<int> read_ints(hid_t grp, const char* name) {
    auto sz = get_dset_size(1, grp, name);
    vector<int> values(sz[0], 0);
    traverse_dset<1, int>(grp, name, [&](size_t i, int x) { values[i] = x; });
    return values;
}

vector<float> read_floats(hid_t grp, const char* name) {
    auto sz = get_dset_size(1, grp, name);
    vector<float> values(sz[0], 0.f);
    traverse_dset<1, float>(grp, name, [&](size_t i, float x) { values[i] = x; });
    return values;
}

vector<string> read_strings(hid_t grp, const char* name) {
    auto sz = get_dset_size(1, grp, name);
    vector<string> values(sz[0]);
    traverse_string_dset<1>(grp, name, [&](size_t i, const string& x) { values[i] = x; });
    return values;
}

void require_size(const char* name, size_t expected, size_t actual) {
    if(expected != actual) {
        throw string("dataset '") + name + "' has size " + to_string(actual) +
            " but expected " + to_string(expected);
    }
}

float clamp_unit(float x) {
    return max(-1.f, min(1.f, x));
}

float pair_energy_and_deriv_over_r(
        float r2,
        float sigma,
        float epsilon,
        float qq,
        float cutoff2,
        float* deriv_over_r)
{
    if(r2 >= cutoff2) {
        if(deriv_over_r) *deriv_over_r = 0.f;
        return 0.f;
    }

    const float safe_r2 = max(r2, 1e-6f);
    const float inv_r2 = 1.f / safe_r2;
    const float sigma2_over_r2 = (sigma * sigma) * inv_r2;
    const float sr6 = sigma2_over_r2 * sigma2_over_r2 * sigma2_over_r2;
    const float sr12 = sr6 * sr6;
    const float lj = 4.f * epsilon * (sr12 - sr6);

    const float inv_r = sqrtf(inv_r2);
    const float coulomb = qq * inv_r;

    if(deriv_over_r) {
        const float lj_deriv_over_r = 24.f * epsilon * (sr6 - 2.f * sr12) * inv_r2;
        const float coulomb_deriv_over_r = -qq * inv_r * inv_r2;
        *deriv_over_r = lj_deriv_over_r + coulomb_deriv_over_r;
    }

    return lj + coulomb;
}

float charmm_pair_energy_and_deriv_over_r(
        float r2,
        float rmin,
        float epsilon,
        float qq,
        float cutoff2,
        float* deriv_over_r)
{
    if(r2 >= cutoff2) {
        if(deriv_over_r) *deriv_over_r = 0.f;
        return 0.f;
    }

    const float safe_r2 = max(r2, 1e-6f);
    const float inv_r2 = 1.f / safe_r2;
    const float rmin2_over_r2 = (rmin * rmin) * inv_r2;
    const float sr6 = rmin2_over_r2 * rmin2_over_r2 * rmin2_over_r2;
    const float sr12 = sr6 * sr6;
    const float lj = epsilon * (sr12 - 2.f * sr6);

    const float inv_r = sqrtf(inv_r2);
    const float coulomb = qq * inv_r;

    if(deriv_over_r) {
        const float lj_deriv_over_r = 12.f * epsilon * (sr6 - sr12) * inv_r2;
        const float coulomb_deriv_over_r = -qq * inv_r * inv_r2;
        *deriv_over_r = lj_deriv_over_r + coulomb_deriv_over_r;
    }

    return lj + coulomb;
}

}  // namespace ligand

namespace {

struct CharmmAtomParam {
    float rmin_half;
    float epsilon;
    float charge;
};

struct CharmmAtom {
    float3 local_pos;
    float rmin_half;
    float epsilon;
    float charge;
};

void load_prepared_rotamer_atoms(
        hid_t grp,
        size_t n_state,
        vector<int>& atom_start,
        vector<int>& atom_count,
        vector<CharmmAtom>& atoms)
{
    atom_start = read_ints(grp, "atom_start");
    atom_count = read_ints(grp, "atom_count");
    require_size("atom_start", n_state, atom_start.size());
    require_size("atom_count", n_state, atom_count.size());

    auto pos_size = get_dset_size(2, grp, "atom_local_pos");
    if(pos_size[1] != 3) {
        throw string("dataset 'atom_local_pos' must have width 3");
    }

    const size_t n_atom = pos_size[0];
    auto atom_rmin_half = read_floats(grp, "atom_rmin_half");
    auto atom_epsilon = read_floats(grp, "atom_epsilon");
    auto atom_charge = read_floats(grp, "atom_charge");
    require_size("atom_rmin_half", n_atom, atom_rmin_half.size());
    require_size("atom_epsilon", n_atom, atom_epsilon.size());
    require_size("atom_charge", n_atom, atom_charge.size());

    vector<float> atom_xyz(n_atom * 3, 0.f);
    traverse_dset<2, float>(grp, "atom_local_pos", [&](size_t i, size_t j, float x) {
        atom_xyz[i * 3 + j] = x;
    });

    atoms.resize(n_atom);
    for(size_t i = 0; i < n_atom; ++i) {
        atoms[i] = CharmmAtom{
            make_vec3(atom_xyz[i * 3 + 0], atom_xyz[i * 3 + 1], atom_xyz[i * 3 + 2]),
            atom_rmin_half[i],
            atom_epsilon[i],
            atom_charge[i],
        };
    }
}

const array<string, 20>& restype_order() {
    static const array<string, 20> names = {{
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    }};
    return names;
}

CharmmAtomParam aliphatic_carbon(float charge = 0.f) {
    return CharmmAtomParam{2.00f, 0.07f, charge};
}

CharmmAtomParam aromatic_carbon(float charge = 0.f) {
    return CharmmAtomParam{1.99f, 0.07f, charge};
}

CharmmAtomParam carbonyl_carbon(float charge = 0.f) {
    return CharmmAtomParam{2.00f, 0.07f, charge};
}

CharmmAtomParam polar_nitrogen(float charge = 0.f) {
    return CharmmAtomParam{1.85f, 0.20f, charge};
}

CharmmAtomParam oxygen_atom(float charge = 0.f) {
    return CharmmAtomParam{1.70f, 0.12f, charge};
}

CharmmAtomParam sulfur_atom(float charge = 0.f) {
    return CharmmAtomParam{2.00f, 0.25f, charge};
}

const vector<CharmmAtomParam>& residue_atom_params(const string& restype) {
    static const map<string, vector<CharmmAtomParam>> params = {
        {"ALA", {aliphatic_carbon()}},
        {"ARG", {
            aliphatic_carbon(), aliphatic_carbon(), aliphatic_carbon(),
            polar_nitrogen(0.00f), carbonyl_carbon(0.50f),
            polar_nitrogen(0.25f), polar_nitrogen(0.25f)}},
        {"ASN", {
            aliphatic_carbon(), carbonyl_carbon(0.40f),
            oxygen_atom(-0.40f), polar_nitrogen(0.00f)}},
        {"ASP", {
            aliphatic_carbon(), carbonyl_carbon(0.00f),
            oxygen_atom(-0.50f), oxygen_atom(-0.50f)}},
        {"CYS", {aliphatic_carbon(), sulfur_atom()}},
        {"GLN", {
            aliphatic_carbon(), aliphatic_carbon(), carbonyl_carbon(0.40f),
            oxygen_atom(-0.40f), polar_nitrogen(0.00f)}},
        {"GLU", {
            aliphatic_carbon(), aliphatic_carbon(), carbonyl_carbon(0.00f),
            oxygen_atom(-0.50f), oxygen_atom(-0.50f)}},
        {"GLY", {}},
        {"HIS", {
            aliphatic_carbon(), aromatic_carbon(), polar_nitrogen(0.00f),
            aromatic_carbon(), aromatic_carbon(), polar_nitrogen(0.00f)}},
        {"ILE", {
            aliphatic_carbon(), aliphatic_carbon(), aliphatic_carbon(), aliphatic_carbon()}},
        {"LEU", {
            aliphatic_carbon(), aliphatic_carbon(), aliphatic_carbon(), aliphatic_carbon()}},
        {"LYS", {
            aliphatic_carbon(), aliphatic_carbon(), aliphatic_carbon(),
            aliphatic_carbon(), polar_nitrogen(1.00f)}},
        {"MET", {aliphatic_carbon(), aliphatic_carbon(), sulfur_atom(), aliphatic_carbon()}},
        {"PHE", {
            aliphatic_carbon(), aromatic_carbon(), aromatic_carbon(), aromatic_carbon(),
            aromatic_carbon(), aromatic_carbon(), aromatic_carbon()}},
        {"PRO", {aliphatic_carbon(), aliphatic_carbon(), aliphatic_carbon()}},
        {"SER", {aliphatic_carbon(), oxygen_atom(0.00f)}},
        {"THR", {aliphatic_carbon(), oxygen_atom(0.00f), aliphatic_carbon()}},
        {"TRP", {
            aliphatic_carbon(), aromatic_carbon(), aromatic_carbon(), aromatic_carbon(),
            polar_nitrogen(0.00f), aromatic_carbon(), aromatic_carbon(),
            aromatic_carbon(), aromatic_carbon(), aromatic_carbon()}},
        {"TYR", {
            aliphatic_carbon(), aromatic_carbon(), aromatic_carbon(), aromatic_carbon(),
            aromatic_carbon(), aromatic_carbon(), aromatic_carbon(), oxygen_atom(0.00f)}},
        {"VAL", {aliphatic_carbon(), aliphatic_carbon(), aliphatic_carbon()}},
    };
    auto it = params.find(restype);
    if(it == params.end()) throw string("No CHARMM atom params for residue type ") + restype;
    return it->second;
}

array<float, 4> fill_missing_chi(const string& restype, const array<float, 4>& chi_in) {
    array<float, 4> chi = chi_in;
    const float trans = M_PI_F;
    if(restype == "GLN" || restype == "GLU") {
        chi[2] = 0.f;
    } else if(restype == "LYS") {
        chi[2] = trans;
        chi[3] = trans;
    } else if(restype == "ARG") {
        chi[2] = trans;
        chi[3] = 0.f;
    } else if(restype == "MET") {
        chi[2] = trans;
    }
    return chi;
}

vector<CharmmAtom> build_sidechain_atoms(int restype_index, const array<float, 4>& chi_in) {
    const auto& names = restype_order();
    if(restype_index < 0 || restype_index >= int(names.size())) {
        throw string("invalid restype index for CHARMM sidechain reconstruction");
    }

    const string& restype = names[restype_index];
    const auto& params = residue_atom_params(restype);
    if(params.empty()) return {};

    MatrixX3f pos;
    auto chi = fill_missing_chi(restype, chi_in);
    res_func_map().at(restype)(pos, 0.f, chi);

    if(int(params.size()) != pos.rows() - 4) {
        throw string("sidechain atom geometry/parameter mismatch for residue type ") + restype;
    }

    vector<CharmmAtom> atoms;
    atoms.reserve(params.size());
    for(int i = 0; i < int(params.size()); ++i) {
        atoms.push_back(CharmmAtom{
            make_vec3(pos(i + 4, 0), pos(i + 4, 1), pos(i + 4, 2)),
            params[i].rmin_half,
            params[i].epsilon,
            params[i].charge,
        });
    }
    return atoms;
}

struct LigandBondPotential : public PotentialNode {
    struct Param {
        int a;
        int b;
        float length;
        float spring;
    };

    CoordNode& pos;
    vector<Param> params;

    LigandBondPotential(hid_t grp, CoordNode& pos_):
        pos(pos_)
    {
        check_elem_width(pos, 3);
        auto atom1 = read_ints(grp, "atom1");
        auto atom2 = read_ints(grp, "atom2");
        auto length = read_floats(grp, "length");
        auto spring = read_floats(grp, "spring");

        require_size("atom2", atom1.size(), atom2.size());
        require_size("length", atom1.size(), length.size());
        require_size("spring", atom1.size(), spring.size());

        params.resize(atom1.size());
        for(size_t i = 0; i < params.size(); ++i) {
            params[i] = Param{atom1[i], atom2[i], length[i], spring[i]};
        }
    }

    void compute_value(ComputeMode mode) override {
        Timer timer("ligand_bond");
        VecArray out = pos.output;
        VecArray sens = pos.sens;
        if(mode == PotentialAndDerivMode) potential = 0.f;

        for(const auto& p : params) {
            auto r1 = load_vec<3>(out, p.a);
            auto r2 = load_vec<3>(out, p.b);
            auto disp = r1 - r2;
            float dist = mag(disp);
            float safe_dist = max(dist, 1e-6f);
            float delta = safe_dist - p.length;
            float deriv_mag = p.spring * delta / safe_dist;
            auto deriv = disp * deriv_mag;

            update_vec(sens, p.a, deriv);
            update_vec(sens, p.b, -deriv);
            if(mode == PotentialAndDerivMode) potential += 0.5f * p.spring * delta * delta;
        }
    }
};

struct LigandAnglePotential : public PotentialNode {
    struct Param {
        int a;
        int b;
        int c;
        float theta0;
        float spring;
    };

    CoordNode& pos;
    vector<Param> params;

    LigandAnglePotential(hid_t grp, CoordNode& pos_):
        pos(pos_)
    {
        check_elem_width(pos, 3);
        auto atom1 = read_ints(grp, "atom1");
        auto atom2 = read_ints(grp, "atom2");
        auto atom3 = read_ints(grp, "atom3");
        auto theta0 = read_floats(grp, "theta0");
        auto spring = read_floats(grp, "spring");

        require_size("atom2", atom1.size(), atom2.size());
        require_size("atom3", atom1.size(), atom3.size());
        require_size("theta0", atom1.size(), theta0.size());
        require_size("spring", atom1.size(), spring.size());

        params.resize(atom1.size());
        for(size_t i = 0; i < params.size(); ++i) {
            params[i] = Param{atom1[i], atom2[i], atom3[i], theta0[i], spring[i]};
        }
    }

    void compute_value(ComputeMode mode) override {
        Timer timer("ligand_angle");
        VecArray out = pos.output;
        VecArray sens = pos.sens;
        if(mode == PotentialAndDerivMode) potential = 0.f;

        for(const auto& p : params) {
            auto r1 = load_vec<3>(out, p.a);
            auto r2 = load_vec<3>(out, p.b);
            auto r3 = load_vec<3>(out, p.c);

            auto x1 = r1 - r2;
            auto x2 = r3 - r2;
            float inv_d1 = inv_mag(x1);
            float inv_d2 = inv_mag(x2);
            auto x1h = x1 * inv_d1;
            auto x2h = x2 * inv_d2;
            float dp = clamp_unit(dot(x1h, x2h));
            float theta = acosf(dp);
            float delta = theta - p.theta0;
            float denom = max(sqrtf(max(1.f - dp * dp, 1e-6f)), 1e-6f);
            float dE_ddp = -p.spring * delta / denom;

            auto deriv1 = (x2h - x1h * dp) * inv_d1;
            auto deriv3 = (x1h - x2h * dp) * inv_d2;
            auto deriv2 = -deriv1 - deriv3;

            update_vec(sens, p.a, dE_ddp * deriv1);
            update_vec(sens, p.b, dE_ddp * deriv2);
            update_vec(sens, p.c, dE_ddp * deriv3);
            if(mode == PotentialAndDerivMode) potential += 0.5f * p.spring * delta * delta;
        }
    }
};

struct LigandTorsionPotential : public PotentialNode {
    struct Param {
        int a;
        int b;
        int c;
        int d;
        float periodicity;
        float phase;
        float amplitude;
    };

    CoordNode& pos;
    vector<Param> params;

    LigandTorsionPotential(hid_t grp, CoordNode& pos_):
        pos(pos_)
    {
        check_elem_width(pos, 3);
        auto atom1 = read_ints(grp, "atom1");
        auto atom2 = read_ints(grp, "atom2");
        auto atom3 = read_ints(grp, "atom3");
        auto atom4 = read_ints(grp, "atom4");
        auto periodicity = read_floats(grp, "periodicity");
        auto phase = read_floats(grp, "phase");
        auto amplitude = read_floats(grp, "amplitude");

        require_size("atom2", atom1.size(), atom2.size());
        require_size("atom3", atom1.size(), atom3.size());
        require_size("atom4", atom1.size(), atom4.size());
        require_size("periodicity", atom1.size(), periodicity.size());
        require_size("phase", atom1.size(), phase.size());
        require_size("amplitude", atom1.size(), amplitude.size());

        params.resize(atom1.size());
        for(size_t i = 0; i < params.size(); ++i) {
            params[i] = Param{
                atom1[i], atom2[i], atom3[i], atom4[i],
                periodicity[i], phase[i], amplitude[i]
            };
        }
    }

    void compute_value(ComputeMode mode) override {
        Timer timer("ligand_torsion");
        VecArray out = pos.output;
        VecArray sens = pos.sens;
        if(mode == PotentialAndDerivMode) potential = 0.f;

        for(const auto& p : params) {
            auto r1 = load_vec<3>(out, p.a);
            auto r2 = load_vec<3>(out, p.b);
            auto r3 = load_vec<3>(out, p.c);
            auto r4 = load_vec<3>(out, p.d);

            float3 d1, d2, d3, d4;
            float phi = dihedral_germ(r1, r2, r3, r4, d1, d2, d3, d4);
            float arg = p.periodicity * phi - p.phase;
            float dE_dphi = -p.amplitude * p.periodicity * sinf(arg);

            update_vec(sens, p.a, dE_dphi * d1);
            update_vec(sens, p.b, dE_dphi * d2);
            update_vec(sens, p.c, dE_dphi * d3);
            update_vec(sens, p.d, dE_dphi * d4);
            if(mode == PotentialAndDerivMode) potential += p.amplitude * (1.f + cosf(arg));
        }
    }
};

struct LigandImproperPotential : public PotentialNode {
    struct Param {
        int a;
        int b;
        int c;
        int d;
        float theta0;
        float spring;
    };

    CoordNode& pos;
    vector<Param> params;

    LigandImproperPotential(hid_t grp, CoordNode& pos_):
        pos(pos_)
    {
        check_elem_width(pos, 3);
        auto atom1 = read_ints(grp, "atom1");
        auto atom2 = read_ints(grp, "atom2");
        auto atom3 = read_ints(grp, "atom3");
        auto atom4 = read_ints(grp, "atom4");
        auto theta0 = read_floats(grp, "theta0");
        auto spring = read_floats(grp, "spring");

        require_size("atom2", atom1.size(), atom2.size());
        require_size("atom3", atom1.size(), atom3.size());
        require_size("atom4", atom1.size(), atom4.size());
        require_size("theta0", atom1.size(), theta0.size());
        require_size("spring", atom1.size(), spring.size());

        params.resize(atom1.size());
        for(size_t i = 0; i < params.size(); ++i) {
            params[i] = Param{atom1[i], atom2[i], atom3[i], atom4[i], theta0[i], spring[i]};
        }
    }

    void compute_value(ComputeMode mode) override {
        Timer timer("ligand_improper");
        VecArray out = pos.output;
        VecArray sens = pos.sens;
        if(mode == PotentialAndDerivMode) potential = 0.f;

        for(const auto& p : params) {
            auto r1 = load_vec<3>(out, p.a);
            auto r2 = load_vec<3>(out, p.b);
            auto r3 = load_vec<3>(out, p.c);
            auto r4 = load_vec<3>(out, p.d);

            float3 d1, d2, d3, d4;
            float phi = dihedral_germ(r1, r2, r3, r4, d1, d2, d3, d4);
            float delta = phi - p.theta0;
            while(delta > M_PI_F) delta -= 2.f * M_PI_F;
            while(delta < -M_PI_F) delta += 2.f * M_PI_F;
            float dE_dphi = p.spring * delta;

            update_vec(sens, p.a, dE_dphi * d1);
            update_vec(sens, p.b, dE_dphi * d2);
            update_vec(sens, p.c, dE_dphi * d3);
            update_vec(sens, p.d, dE_dphi * d4);
            if(mode == PotentialAndDerivMode) potential += 0.5f * p.spring * delta * delta;
        }
    }
};

struct LigandNonbondedPotential : public PotentialNode {
    CoordNode& pos;
    vector<int> atom1;
    vector<int> atom2;
    vector<float> sigma;
    vector<float> epsilon;
    vector<float> charge;
    float cutoff2;

    LigandNonbondedPotential(hid_t grp, CoordNode& pos_):
        pos(pos_),
        atom1(read_ints(grp, "atom1")),
        atom2(read_ints(grp, "atom2")),
        sigma(read_floats(grp, "sigma")),
        epsilon(read_floats(grp, "epsilon")),
        charge(read_floats(grp, "charge")),
        cutoff2(powf(read_attribute<float>(grp, ".", "cutoff", 12.f), 2.f))
    {
        check_elem_width(pos, 3);
        require_size("atom2", atom1.size(), atom2.size());
        require_size("sigma", atom1.size(), sigma.size());
        require_size("epsilon", atom1.size(), epsilon.size());
        require_size("charge", atom1.size(), charge.size());
    }

    void compute_value(ComputeMode mode) override {
        Timer timer("ligand_nonbonded");
        VecArray out = pos.output;
        VecArray sens = pos.sens;
        if(mode == PotentialAndDerivMode) potential = 0.f;

        for(size_t i = 0; i < atom1.size(); ++i) {
            auto r1 = load_vec<3>(out, atom1[i]);
            auto r2 = load_vec<3>(out, atom2[i]);
            auto disp = r1 - r2;
            float deriv_over_r = 0.f;
            float energy = pair_energy_and_deriv_over_r(
                mag2(disp), sigma[i], epsilon[i], charge[i], cutoff2, &deriv_over_r);
            auto deriv = disp * deriv_over_r;
            update_vec(sens, atom1[i], deriv);
            update_vec(sens, atom2[i], -deriv);
            if(mode == PotentialAndDerivMode) potential += energy;
        }
    }
};

struct LigandNonbondedCharmmPotential : public PotentialNode {
    CoordNode& pos;
    vector<int> atom1;
    vector<int> atom2;
    vector<float> rmin;
    vector<float> epsilon;
    vector<float> charge;
    float cutoff2;

    LigandNonbondedCharmmPotential(hid_t grp, CoordNode& pos_):
        pos(pos_),
        atom1(read_ints(grp, "atom1")),
        atom2(read_ints(grp, "atom2")),
        rmin(read_floats(grp, "rmin")),
        epsilon(read_floats(grp, "epsilon")),
        charge(read_floats(grp, "charge")),
        cutoff2(powf(read_attribute<float>(grp, ".", "cutoff", 12.f), 2.f))
    {
        check_elem_width(pos, 3);
        require_size("atom2", atom1.size(), atom2.size());
        require_size("rmin", atom1.size(), rmin.size());
        require_size("epsilon", atom1.size(), epsilon.size());
        require_size("charge", atom1.size(), charge.size());
    }

    void compute_value(ComputeMode mode) override {
        Timer timer("charmm_ligand_nonbonded");
        VecArray out = pos.output;
        VecArray sens = pos.sens;
        if(mode == PotentialAndDerivMode) potential = 0.f;

        for(size_t i = 0; i < atom1.size(); ++i) {
            auto r1 = load_vec<3>(out, atom1[i]);
            auto r2 = load_vec<3>(out, atom2[i]);
            auto disp = r1 - r2;
            float deriv_over_r = 0.f;
            float energy = charmm_pair_energy_and_deriv_over_r(
                mag2(disp), rmin[i], epsilon[i], charge[i], cutoff2, &deriv_over_r);
            auto deriv = disp * deriv_over_r;
            update_vec(sens, atom1[i], deriv);
            update_vec(sens, atom2[i], -deriv);
            if(mode == PotentialAndDerivMode) potential += energy;
        }
    }
};

struct LigandBackboneCharmmPotential : public PotentialNode {
    CoordNode& ligand_pos;
    CoordNode& protein_pos;
    vector<int> protein_index;
    vector<float> ligand_rmin_half;
    vector<float> ligand_epsilon;
    vector<float> ligand_charge;
    vector<float> protein_rmin_half;
    vector<float> protein_epsilon;
    vector<float> protein_charge;
    float cutoff2;

    LigandBackboneCharmmPotential(hid_t grp, CoordNode& ligand_pos_, CoordNode& protein_pos_):
        ligand_pos(ligand_pos_),
        protein_pos(protein_pos_),
        protein_index(read_ints(grp, "protein_index")),
        ligand_rmin_half(read_floats(grp, "ligand_rmin_half")),
        ligand_epsilon(read_floats(grp, "ligand_epsilon")),
        ligand_charge(read_floats(grp, "ligand_charge")),
        protein_rmin_half(read_floats(grp, "protein_rmin_half")),
        protein_epsilon(read_floats(grp, "protein_epsilon")),
        protein_charge(read_floats(grp, "protein_charge")),
        cutoff2(powf(read_attribute<float>(grp, ".", "cutoff", 12.f), 2.f))
    {
        check_elem_width(ligand_pos, 3);
        check_elem_width_lower_bound(protein_pos, 3);
        require_size("ligand_rmin_half", ligand_pos.n_elem, ligand_rmin_half.size());
        require_size("ligand_epsilon", ligand_pos.n_elem, ligand_epsilon.size());
        require_size("ligand_charge", ligand_pos.n_elem, ligand_charge.size());
        require_size("protein_rmin_half", protein_index.size(), protein_rmin_half.size());
        require_size("protein_epsilon", protein_index.size(), protein_epsilon.size());
        require_size("protein_charge", protein_index.size(), protein_charge.size());
    }

    void compute_value(ComputeMode mode) override {
        Timer timer("ligand_backbone");
        VecArray lout = ligand_pos.output;
        VecArray pout = protein_pos.output;
        VecArray lsens = ligand_pos.sens;
        VecArray psens = protein_pos.sens;
        if(mode == PotentialAndDerivMode) potential = 0.f;

        for(int nl = 0; nl < ligand_pos.n_elem; ++nl) {
            auto r_lig = load_vec<3>(lout, nl);
            for(size_t np = 0; np < protein_index.size(); ++np) {
                auto r_pro = make_vec3(
                    pout(0, protein_index[np]),
                    pout(1, protein_index[np]),
                    pout(2, protein_index[np]));
                auto disp = r_lig - r_pro;
                const float rmin = ligand_rmin_half[nl] + protein_rmin_half[np];
                const float epsilon_mix = sqrtf(max(ligand_epsilon[nl] * protein_epsilon[np], 0.f));
                const float qq = ligand_charge[nl] * protein_charge[np];
                float deriv_over_r = 0.f;
                float energy = charmm_pair_energy_and_deriv_over_r(
                    mag2(disp), rmin, epsilon_mix, qq, cutoff2, &deriv_over_r);
                auto deriv = disp * deriv_over_r;
                update_vec(lsens, nl, deriv);
                update_vec(psens, protein_index[np], -deriv);
                if(mode == PotentialAndDerivMode) potential += energy;
            }
        }
    }
};

struct LigandRotamerCharmmEnergy : public CoordNode {
    CoordNode& alignment;
    CoordNode& ligand_pos;
    vector<int> affine_residue;
    vector<int> atom_start;
    vector<int> atom_count;
    vector<CharmmAtom> atoms;
    vector<float> ligand_rmin_half;
    vector<float> ligand_epsilon;
    vector<float> ligand_charge;
    float cutoff2;

    LigandRotamerCharmmEnergy(hid_t grp, CoordNode& alignment_, CoordNode& ligand_pos_):
        CoordNode(get_dset_size(1, grp, "affine_residue")[0], 1),
        alignment(alignment_),
        ligand_pos(ligand_pos_),
        affine_residue(read_ints(grp, "affine_residue")),
        atom_start(affine_residue.size(), 0),
        atom_count(affine_residue.size(), 0),
        ligand_rmin_half(read_floats(grp, "ligand_rmin_half")),
        ligand_epsilon(read_floats(grp, "ligand_epsilon")),
        ligand_charge(read_floats(grp, "ligand_charge")),
        cutoff2(powf(read_attribute<float>(grp, ".", "cutoff", 12.f), 2.f))
    {
        check_elem_width(alignment, 7);
        check_elem_width(ligand_pos, 3);
        require_size("ligand_rmin_half", ligand_pos.n_elem, ligand_rmin_half.size());
        require_size("ligand_epsilon", ligand_pos.n_elem, ligand_epsilon.size());
        require_size("ligand_charge", ligand_pos.n_elem, ligand_charge.size());

        if(h5_exists(grp, "atom_local_pos")) {
            load_prepared_rotamer_atoms(grp, affine_residue.size(), atom_start, atom_count, atoms);
        } else {
            auto restype_index = read_ints(grp, "restype_index");
            require_size("restype_index", affine_residue.size(), restype_index.size());

            auto chi_size = get_dset_size(2, grp, "chi");
            if(chi_size[0] != affine_residue.size() || chi_size[1] != 4) {
                throw string("dataset 'chi' has wrong shape for ligand_rotamer_charmm_1body");
            }
            vector<array<float, 4>> chi(affine_residue.size());
            traverse_dset<2, float>(grp, "chi", [&](size_t i, size_t j, float x) { chi[i][j] = x; });

            for(size_t i = 0; i < affine_residue.size(); ++i) {
                atom_start[i] = atoms.size();
                auto local_atoms = build_sidechain_atoms(restype_index[i], chi[i]);
                atom_count[i] = local_atoms.size();
                atoms.insert(atoms.end(), local_atoms.begin(), local_atoms.end());
            }
        }
    }

    void compute_value(ComputeMode mode) override {
        Timer timer("ligand_rotamer_charmm");
        VecArray aout = alignment.output;
        VecArray lout = ligand_pos.output;

        for(int ns = 0; ns < n_elem; ++ns) {
            auto aff = load_vec<7>(aout, affine_residue[ns]);
            auto t = extract<0, 3>(aff);
            float U[9];
            quat_to_rot(U, aff.v + 3);

            float energy = 0.f;
            const int start = atom_start[ns];
            const int stop = start + atom_count[ns];
            for(int na = start; na < stop; ++na) {
                auto r_sc = apply_affine(U, t, atoms[na].local_pos);
                for(int nl = 0; nl < ligand_pos.n_elem; ++nl) {
                    auto r_lig = load_vec<3>(lout, nl);
                    auto disp = r_sc - r_lig;
                    const float rmin = atoms[na].rmin_half + ligand_rmin_half[nl];
                    const float epsilon_mix = sqrtf(max(atoms[na].epsilon * ligand_epsilon[nl], 0.f));
                    const float qq = atoms[na].charge * ligand_charge[nl];
                    energy += charmm_pair_energy_and_deriv_over_r(
                        mag2(disp), rmin, epsilon_mix, qq, cutoff2, nullptr);
                }
            }
            output(0, ns) = energy;
        }
    }

    void propagate_deriv() override {
        Timer timer("d_ligand_rotamer_charmm");
        VecArray aout = alignment.output;
        VecArray asens = alignment.sens;
        VecArray lout = ligand_pos.output;
        VecArray lsens = ligand_pos.sens;

        for(int ns = 0; ns < n_elem; ++ns) {
            const float upstream = sens(0, ns);
            if(upstream == 0.f) continue;

            auto aff = load_vec<7>(aout, affine_residue[ns]);
            auto t = extract<0, 3>(aff);
            float U[9];
            quat_to_rot(U, aff.v + 3);

            float3 com_deriv = make_zero<3>();
            float3 torque = make_zero<3>();

            const int start = atom_start[ns];
            const int stop = start + atom_count[ns];
            for(int na = start; na < stop; ++na) {
                auto r_sc = apply_affine(U, t, atoms[na].local_pos);
                float3 atom_sens = make_zero<3>();
                for(int nl = 0; nl < ligand_pos.n_elem; ++nl) {
                    auto r_lig = load_vec<3>(lout, nl);
                    auto disp = r_sc - r_lig;
                    const float rmin = atoms[na].rmin_half + ligand_rmin_half[nl];
                    const float epsilon_mix = sqrtf(max(atoms[na].epsilon * ligand_epsilon[nl], 0.f));
                    const float qq = atoms[na].charge * ligand_charge[nl];
                    float deriv_over_r = 0.f;
                    charmm_pair_energy_and_deriv_over_r(
                        mag2(disp), rmin, epsilon_mix, qq, cutoff2, &deriv_over_r);
                    auto deriv = disp * (upstream * deriv_over_r);
                    atom_sens += deriv;
                    update_vec(lsens, nl, -deriv);
                }
                com_deriv += atom_sens;
                torque += cross(r_sc - t, atom_sens);
            }

            update_vec(&asens(0, affine_residue[ns]), com_deriv);
            update_vec(&asens(3, affine_residue[ns]), torque);
        }
    }
};

struct LigandRotamerEnergy : public CoordNode {
    CoordNode& sc_pos;
    CoordNode& ligand_pos;
    vector<float> sc_sigma;
    vector<float> sc_epsilon;
    vector<float> sc_charge;
    vector<float> ligand_sigma;
    vector<float> ligand_epsilon;
    vector<float> ligand_charge;
    float cutoff2;

    LigandRotamerEnergy(hid_t grp, CoordNode& sc_pos_, CoordNode& ligand_pos_):
        CoordNode(sc_pos_.n_elem, 1),
        sc_pos(sc_pos_),
        ligand_pos(ligand_pos_),
        sc_sigma(read_floats(grp, "sc_sigma")),
        sc_epsilon(read_floats(grp, "sc_epsilon")),
        sc_charge(read_floats(grp, "sc_charge")),
        ligand_sigma(read_floats(grp, "ligand_sigma")),
        ligand_epsilon(read_floats(grp, "ligand_epsilon")),
        ligand_charge(read_floats(grp, "ligand_charge")),
        cutoff2(powf(read_attribute<float>(grp, ".", "cutoff", 12.f), 2.f))
    {
        check_elem_width_lower_bound(sc_pos, 3);
        check_elem_width(ligand_pos, 3);
        require_size("sc_sigma", sc_pos.n_elem, sc_sigma.size());
        require_size("sc_epsilon", sc_pos.n_elem, sc_epsilon.size());
        require_size("sc_charge", sc_pos.n_elem, sc_charge.size());
        require_size("ligand_sigma", ligand_pos.n_elem, ligand_sigma.size());
        require_size("ligand_epsilon", ligand_pos.n_elem, ligand_epsilon.size());
        require_size("ligand_charge", ligand_pos.n_elem, ligand_charge.size());
    }

    void compute_value(ComputeMode mode) override {
        Timer timer("ligand_rotamer");
        VecArray sout = sc_pos.output;
        VecArray lout = ligand_pos.output;

        for(int ns = 0; ns < n_elem; ++ns) {
            auto r_sc = make_vec3(sout(0, ns), sout(1, ns), sout(2, ns));
            float energy = 0.f;
            for(int nl = 0; nl < ligand_pos.n_elem; ++nl) {
                auto r_lig = load_vec<3>(lout, nl);
                auto disp = r_sc - r_lig;
                float sigma_mix = 0.5f * (sc_sigma[ns] + ligand_sigma[nl]);
                float epsilon_mix = sqrtf(max(sc_epsilon[ns] * ligand_epsilon[nl], 0.f));
                float qq = sc_charge[ns] * ligand_charge[nl];
                energy += pair_energy_and_deriv_over_r(
                    mag2(disp), sigma_mix, epsilon_mix, qq, cutoff2, nullptr);
            }
            output(0, ns) = energy;
        }
    }

    void propagate_deriv() override {
        Timer timer("d_ligand_rotamer");
        VecArray sout = sc_pos.output;
        VecArray lout = ligand_pos.output;
        VecArray ssens = sc_pos.sens;
        VecArray lsens = ligand_pos.sens;

        for(int ns = 0; ns < n_elem; ++ns) {
            float upstream = sens(0, ns);
            if(upstream == 0.f) continue;
            auto r_sc = make_vec3(sout(0, ns), sout(1, ns), sout(2, ns));
            for(int nl = 0; nl < ligand_pos.n_elem; ++nl) {
                auto r_lig = load_vec<3>(lout, nl);
                auto disp = r_sc - r_lig;
                float sigma_mix = 0.5f * (sc_sigma[ns] + ligand_sigma[nl]);
                float epsilon_mix = sqrtf(max(sc_epsilon[ns] * ligand_epsilon[nl], 0.f));
                float qq = sc_charge[ns] * ligand_charge[nl];
                float deriv_over_r = 0.f;
                pair_energy_and_deriv_over_r(
                    mag2(disp), sigma_mix, epsilon_mix, qq, cutoff2, &deriv_over_r);
                auto deriv = disp * (upstream * deriv_over_r);
                ssens(0, ns) += deriv[0];
                ssens(1, ns) += deriv[1];
                ssens(2, ns) += deriv[2];
                update_vec(lsens, nl, -deriv);
            }
        }
    }
};

}  // namespace

static RegisterNodeType<LigandBondPotential, 1> ligand_bond_node("ligand_harmonic_bond");
static RegisterNodeType<LigandAnglePotential, 1> ligand_angle_node("ligand_angle");
static RegisterNodeType<LigandTorsionPotential, 1> ligand_torsion_node("ligand_torsion");
static RegisterNodeType<LigandImproperPotential, 1> ligand_improper_node("ligand_improper");
static RegisterNodeType<LigandNonbondedPotential, 1> ligand_nonbonded_node("ligand_nonbonded");
static RegisterNodeType<LigandNonbondedCharmmPotential, 1> ligand_nonbonded_charmm_node("charmm_ligand_nonbonded");
static RegisterNodeType<LigandBackboneCharmmPotential, 2> ligand_backbone_node("ligand_backbone_nonbonded");
static RegisterNodeType<LigandBackboneCharmmPotential, 2> ligand_virtual_backbone_node("ligand_virtual_backbone_nonbonded");
static RegisterNodeType<LigandRotamerCharmmEnergy, 2> ligand_rotamer_charmm_node("ligand_rotamer_charmm_1body");
static RegisterNodeType<LigandRotamerEnergy, 2> ligand_rotamer_node("ligand_rotamer_1body");
