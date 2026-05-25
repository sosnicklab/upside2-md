#include "martini_internal.h"
#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include <mutex>
#include "spline.h"
#include <iostream>
#include <H5Apublic.h> // for H5Aexists
#include <cmath> // For pow, cosf, sinf, acosf
#include <cctype>
#include <cstdint>
#include <set> // For std::set
#include <array>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <memory>
#include <limits>
#include "box.h" // For PBC minimum_image function

using namespace h5;
using namespace std;

inline bool attribute_exists(hid_t loc_id, const char* obj_name, const char* attr_name) {
    // H5Oopen(loc, ".") may alias the same ID in some HDF5 builds.
    // Closing that alias can invalidate the caller's handle.
    if(!obj_name || obj_name[0] == '\0' || (obj_name[0] == '.' && obj_name[1] == '\0')) {
        htri_t exists = H5Aexists(loc_id, attr_name);
        return exists > 0;
    }

    hid_t obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT);
    if (obj_id < 0) return false;
    htri_t exists = H5Aexists(obj_id, attr_name);
    H5Oclose(obj_id);
    return exists > 0;
}

struct DihedralSpring : public PotentialNode
{
    struct Params {
        index_t atom[4];
        float equil_dihedral;
        float spring_constant;
        int dihedral_type;  // 1=periodic, 2=harmonic
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    
    // Box dimensions for minimum image
    float box_x, box_y, box_z;
    
    // Spline interpolation for dihedral potential
    LayeredClampedSpline1D<1> dihedral_potential_spline;
    
    // Spline parameters
    float dihedral_min, dihedral_max;
    float max_spring;  // Store max spring constant for scaling

    DihedralSpring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(2, grp, "id")[0]), pos(pos_), params(n_elem),
        dihedral_potential_spline(1, 1000)
    {
        int n_dep = 4;  // number of atoms that each term depends on 
        check_size(grp, "id",           n_elem, n_dep);
        check_size(grp, "equil_dist",   n_elem);
        check_size(grp, "spring_const", n_elem);

        auto& p = params;
        traverse_dset<2,int>  (grp, "id",           [&](size_t i, size_t j, int   x) {p[i].atom[j]  =x;});
        traverse_dset<1,float>(grp, "equil_dist",   [&](size_t i,           float x) {p[i].equil_dihedral =x;});
        traverse_dset<1,float>(grp, "spring_const", [&](size_t i,           float x) {p[i].spring_constant=x;});
        // Read dihedral type (1=periodic, 2=harmonic), default to periodic if not present
        if(attribute_exists(grp, ".", "dihedral_type")) {
            traverse_dset<1,int>(grp, "dihedral_type", [&](size_t i, int x) {p[i].dihedral_type = x;});
        } else {
            // Default to periodic dihedrals for backward compatibility
            for(auto& param : params) param.dihedral_type = 1;
        }

        // Read box dimensions (same as MartiniPotential)
        // Support both individual dimensions and legacy wall_box_size
        if(attribute_exists(grp, ".", "x_len") && attribute_exists(grp, ".", "y_len") && attribute_exists(grp, ".", "z_len")) {
            // New format: separate x, y, z dimensions
            box_x = read_attribute<float>(grp, ".", "x_len");
            box_y = read_attribute<float>(grp, ".", "y_len");
            box_z = read_attribute<float>(grp, ".", "z_len");
        } else {
            // Legacy format: wall boundaries
            float wall_xlo = read_attribute<float>(grp, ".", "wall_xlo");
            float wall_xhi = read_attribute<float>(grp, ".", "wall_xhi");
            float wall_ylo = read_attribute<float>(grp, ".", "wall_ylo");
            float wall_yhi = read_attribute<float>(grp, ".", "wall_yhi");
            float wall_zlo = read_attribute<float>(grp, ".", "wall_zlo");
            float wall_zhi = read_attribute<float>(grp, ".", "wall_zhi");
            box_x = wall_xhi - wall_xlo;
            box_y = wall_yhi - wall_ylo;
            box_z = wall_zhi - wall_zlo;
        }
        
        // Initialize spline parameters for dihedral potential
        // Find the range of equilibrium dihedrals and spring constants
        float min_dihedral = std::numeric_limits<float>::max();
        float max_dihedral = std::numeric_limits<float>::lowest();
        float max_spring = std::numeric_limits<float>::lowest();
        
        for(const auto& p : params) {
            min_dihedral = std::min(min_dihedral, p.equil_dihedral);
            max_dihedral = std::max(max_dihedral, p.equil_dihedral);
            max_spring = std::max(max_spring, p.spring_constant);
        }
        
        // Store max_spring as member variable
        this->max_spring = max_spring;
        
        // Set spline range for dihedral angles (typically -pi to pi radians)
        dihedral_min = -M_PI_F;
        dihedral_max = M_PI_F;
        
        // Use periodic dihedral potential: V = k * (1 + cos(n*phi - phi0))
        // For MARTINI, n=1 (multiplicity=1) and phi0 is the equilibrium angle
        // We'll use a canonical spline for V = 1 + cos(phi - phi0)
        std::vector<double> dihedral_pot_data(1000);
        for(int i = 0; i < 1000; ++i) {
            float phi = dihedral_min + i * (dihedral_max - dihedral_min) / 999.0f;
            dihedral_pot_data[i] = 1.0 + cos(phi);  // Canonical form: 1 + cos(phi)
        }
        
        // Fit spline
        dihedral_potential_spline.fit_spline(dihedral_pot_data.data());

    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("dihedral_spring"));

        float* posc = pos.output.x.get();
        float* pos_sens = pos.sens.x.get();
        auto hybrid_state = martini_hybrid::get_state_for_coord(pos);
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;

        for(int nt=0; nt<n_elem; ++nt) {
            const auto& p = params[nt];
            if(hybrid_state && !martini_hybrid::allow_multibody_term_if_active(*hybrid_state, p.atom, 4)) {
                continue;
            }
            Float4 x_orig[4];
            for(int na: range(4)) x_orig[na] = Float4(posc + 4*params[nt].atom[na]);

            // Apply minimum image convention for periodic boundaries
            // Reconstruct chain: use atom 1 as reference, apply PBC to get other atoms
            Float4 x[4];
            x[1] = x_orig[1]; // Use atom 1 as reference
            
            // Apply PBC minimum image to reconstruct unbroken chain
            // Atom 0 relative to atom 1
            float dx0 = x_orig[0].x() - x_orig[1].x();
            float dy0 = x_orig[0].y() - x_orig[1].y();
            float dz0 = x_orig[0].z() - x_orig[1].z();
            simulation_box::minimum_image_scalar(dx0, dy0, dz0, box_x, box_y, box_z);
            {
                alignas(16) float disp0[4] = {dx0, dy0, dz0, 0.f};
                x[0] = x[1] + Float4(disp0);
            }
            
            // Atom 2 relative to atom 1
            float dx2 = x_orig[2].x() - x_orig[1].x();
            float dy2 = x_orig[2].y() - x_orig[1].y();
            float dz2 = x_orig[2].z() - x_orig[1].z();
            simulation_box::minimum_image_scalar(dx2, dy2, dz2, box_x, box_y, box_z);
            {
                alignas(16) float disp2[4] = {dx2, dy2, dz2, 0.f};
                x[2] = x[1] + Float4(disp2);
            }
            
            // Atom 3 relative to atom 2 (chain: 0-1-2-3)
            float dx3 = x_orig[3].x() - x_orig[2].x();
            float dy3 = x_orig[3].y() - x_orig[2].y();
            float dz3 = x_orig[3].z() - x_orig[2].z();
            simulation_box::minimum_image_scalar(dx3, dy3, dz3, box_x, box_y, box_z);
            {
                alignas(16) float disp3[4] = {dx3, dy3, dz3, 0.f};
                x[3] = x[2] + Float4(disp3);
            }

            Float4 d[4];
            float dihedral = dihedral_germ(x[0],x[1],x[2],x[3], d[0],d[1],d[2],d[3]).x();

            // Use appropriate dihedral potential based on type
            float delta_phi = dihedral - p.equil_dihedral;
            if(delta_phi > M_PI_F) delta_phi -= 2.0f * M_PI_F;
            if(delta_phi < -M_PI_F) delta_phi += 2.0f * M_PI_F;
            
            if(p.dihedral_type == 1) {
                // Periodic dihedral: V = k * (1 + cos(phi - phi0))
                if(delta_phi >= dihedral_min && delta_phi <= dihedral_max) {
                    // Transform delta_phi to spline coordinate [0, 999]
                    float phi_coord = (delta_phi - dihedral_min) / (dihedral_max - dihedral_min) * 999.0f;
                    
                    // Get potential and force from spline
                    float result[2];
                    dihedral_potential_spline.evaluate_value_and_deriv(result, 0, phi_coord);
                    float dihedral_pot = result[1]; // Index 1 is the value (1 + cos(delta_phi))
                    float dE_ddelta_spline = result[0]; // Index 0 is the derivative
                    
                    // Convert derivative from spline coordinate to physical coordinate
                    float coord_scale = 999.0f / (dihedral_max - dihedral_min);
                    float dE_ddelta = dE_ddelta_spline * coord_scale;
                    float dihedral_force_mag = p.spring_constant * dE_ddelta;
                    
                    // Use periodic dihedral potential: V = k * (1 + cos(phi - phi0))
                    if(pot) *pot += p.spring_constant * dihedral_pot;

                    // Apply force with mass scaling (dihedral_force_mag corresponds to dE/dphi)
                    auto s = Float4(dihedral_force_mag);
                    for(int na: range(4)) d[na].scale_update(s, pos_sens + 4*params[nt].atom[na]);
                } else {
                    // Fallback to direct calculation for out-of-range periodic dihedrals
                    if(pot) *pot += p.spring_constant * (1.0f + cos(delta_phi));
                    
                    // Force is the derivative of V = k * (1 + cos(phi - phi0))
                    // dV/dphi = k * (-sin(phi - phi0)) = -k * sin(delta_phi)
                    auto s = Float4(-p.spring_constant * sin(delta_phi));
                    for(int na: range(4)) d[na].scale_update(s, pos_sens + 4*params[nt].atom[na]);
                }
            } else if(p.dihedral_type == 2) {
                // Harmonic dihedral: V = 0.5 * k * (phi - phi0)^2
                if(pot) *pot += 0.5f * p.spring_constant * delta_phi * delta_phi;
                
                // Force is the derivative of V = 0.5 * k * (phi - phi0)^2
                // dV/dphi = k * (phi - phi0) = k * delta_phi
                auto s = Float4(p.spring_constant * delta_phi);
                for(int na: range(4)) d[na].scale_update(s, pos_sens + 4*params[nt].atom[na]);
            } else {
                // Unknown dihedral type - use harmonic as fallback
                if(pot) *pot += 0.5f * p.spring_constant * delta_phi * delta_phi;
                
                auto s = Float4(p.spring_constant * delta_phi);
                for(int na: range(4)) d[na].scale_update(s, pos_sens + 4*params[nt].atom[na]);
            }
        }
    }
    
};
static RegisterNodeType<DihedralSpring,1> dihedral_spring_node("dihedral_spring");

// MARTINI potential using spline interpolation for LJ and Coulomb calculations
struct MartiniPotential : public PotentialNode
{
    int n_atom;
    CoordNode& pos;

    struct PairParam {
        float eps;
        float sig;
        float qi;
        float qj;
        float qq;
        const LayeredClampedSpline1D<1>* combined_spline;

        PairParam():
            eps(0.f), sig(0.f), qi(0.f), qj(0.f), qq(0.f),
            combined_spline(nullptr) {}
    };

    vector<PairParam> param_table;
    vector<uint32_t> pair_param_index;
    vector<pair<int,int>> pairs;
    
    float epsilon, sigma, lj_cutoff, coul_cutoff;
    bool force_cap;
    float energy_conversion_kj_per_eup;
    float length_conversion_angstrom_per_nm;
    float coulomb_constant_native_kj_mol_nm_e2;
    float coulomb_k;

    // Box dimensions used for minimum-image pair displacements under PBC/NPT.
    float box_x, box_y, box_z;
    
    // Combined LJ+Coulomb spline, one spline per (eps, sig, qq) triple.
    std::map<std::tuple<float, float, float>, LayeredClampedSpline1D<1>> combined_splines;

    float r_min, r_max;
    float r_shift, r_scale;

    float cache_buffer;
    float pairlist_cutoff;
    bool pairlist_valid;
    float cached_box_x;
    float cached_box_y;
    float cached_box_z;
    vector<float> cached_pos;
    vector<int32_t> active_pair_indices;

    inline bool pairlist_needs_rebuild(const VecArray& pos1) const {
        if(!pairlist_valid) return true;
        if(cache_buffer <= 0.f) return true;
        if(fabsf(box_x - cached_box_x) > 1.0e-6f ||
           fabsf(box_y - cached_box_y) > 1.0e-6f ||
           fabsf(box_z - cached_box_z) > 1.0e-6f) return true;

        float max_cache_dist2 = sqr(0.5f * cache_buffer);
        for(int na = 0; na < n_atom; ++na) {
            float dx = pos1(0, na) - cached_pos[3*na + 0];
            float dy = pos1(1, na) - cached_pos[3*na + 1];
            float dz = pos1(2, na) - cached_pos[3*na + 2];
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f) {
                simulation_box::minimum_image_scalar(dx, dy, dz, box_x, box_y, box_z);
            }
            if(dx*dx + dy*dy + dz*dz > max_cache_dist2) return true;
        }
        return false;
    }

    void rebuild_pairlist(const VecArray& pos1) {
        Timer timer(string("martini_pairlist_rebuild"));
        float active_cutoff = pairlist_cutoff + std::max(0.f, cache_buffer);
        float active_cutoff2 = sqr(active_cutoff);

        active_pair_indices.clear();
        active_pair_indices.reserve(std::min<size_t>(pairs.size(), size_t(1u << 20)));

        for(size_t np = 0; np < pairs.size(); ++np) {
            const auto& param = param_table[pair_param_index[np]];
            if(param.eps == 0.f && param.sig == 0.f && param.qi == 0.f && param.qj == 0.f) continue;

            int i = pairs[np].first;
            int j = pairs[np].second;
            auto dr = load_vec<3>(pos1, i) - load_vec<3>(pos1, j);
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f) {
                dr = simulation_box::minimum_image(dr, box_x, box_y, box_z);
            }
            if(mag2(dr) < active_cutoff2) {
                active_pair_indices.push_back(int32_t(np));
            }
        }

        for(int na = 0; na < n_atom; ++na) {
            cached_pos[3*na + 0] = pos1(0, na);
            cached_pos[3*na + 1] = pos1(1, na);
            cached_pos[3*na + 2] = pos1(2, na);
        }
        cached_box_x = box_x;
        cached_box_y = box_y;
        cached_box_z = box_z;
        pairlist_valid = true;
    }
    
    MartiniPotential(hid_t grp, CoordNode& pos_):
        PotentialNode(), n_atom(pos_.n_elem), pos(pos_),
        cache_buffer(1.f),
        pairlist_cutoff(0.f),
        pairlist_valid(false),
        cached_box_x(0.f),
        cached_box_y(0.f),
        cached_box_z(0.f),
        cached_pos(3*n_atom, 0.f)
    {
        check_size(grp, "atom_indices", n_atom);
        check_size(grp, "charges", n_atom);
        check_size(grp, "/input/type", n_atom);
        
        epsilon     = read_attribute<float>(grp, ".", "epsilon");
        sigma       = read_attribute<float>(grp, ".", "sigma");  
        lj_cutoff   = read_attribute<float>(grp, ".", "lj_cutoff");
        coul_cutoff = read_attribute<float>(grp, ".", "coul_cutoff");
        cache_buffer = read_attribute<float>(grp, ".", "cache_buffer", 1.f);
        if(cache_buffer < 0.f) cache_buffer = 0.f;
        pairlist_cutoff = max(lj_cutoff, coul_cutoff);
        if(!attribute_exists(grp, ".", "energy_conversion_kj_per_eup") ||
           !attribute_exists(grp, ".", "length_conversion_angstrom_per_nm") ||
           !attribute_exists(grp, ".", "coulomb_constant_native_kj_mol_nm_e2")) {
            throw string("martini_potential requires explicit unit-conversion attrs: "
                         "energy_conversion_kj_per_eup, length_conversion_angstrom_per_nm, "
                         "coulomb_constant_native_kj_mol_nm_e2");
        }
        energy_conversion_kj_per_eup =
            read_attribute<float>(grp, ".", "energy_conversion_kj_per_eup");
        length_conversion_angstrom_per_nm =
            read_attribute<float>(grp, ".", "length_conversion_angstrom_per_nm");
        coulomb_constant_native_kj_mol_nm_e2 =
            read_attribute<float>(grp, ".", "coulomb_constant_native_kj_mol_nm_e2");
        if(!(energy_conversion_kj_per_eup > 0.f) || !(length_conversion_angstrom_per_nm > 0.f)) {
            throw string("martini_potential unit-conversion attrs must be positive");
        }
        coulomb_k =
            coulomb_constant_native_kj_mol_nm_e2 *
            (length_conversion_angstrom_per_nm / energy_conversion_kj_per_eup);

        force_cap = true;
        if(attribute_exists(grp, ".", "force_cap")) {
            force_cap = read_attribute<int>(grp, ".", "force_cap") != 0;
        }
        
        // Read box dimensions for minimum-image displacements.
        if(attribute_exists(grp, ".", "x_len") && attribute_exists(grp, ".", "y_len") && attribute_exists(grp, ".", "z_len")) {
            box_x = read_attribute<float>(grp, ".", "x_len");
            box_y = read_attribute<float>(grp, ".", "y_len");
            box_z = read_attribute<float>(grp, ".", "z_len");
        } else if(attribute_exists(grp, ".", "wall_xlo") && attribute_exists(grp, ".", "wall_xhi") &&
                  attribute_exists(grp, ".", "wall_ylo") && attribute_exists(grp, ".", "wall_yhi") &&
                  attribute_exists(grp, ".", "wall_zlo") && attribute_exists(grp, ".", "wall_zhi")) {
            float wall_xlo = read_attribute<float>(grp, ".", "wall_xlo");
            float wall_xhi = read_attribute<float>(grp, ".", "wall_xhi");
            float wall_ylo = read_attribute<float>(grp, ".", "wall_ylo");
            float wall_yhi = read_attribute<float>(grp, ".", "wall_yhi");
            float wall_zlo = read_attribute<float>(grp, ".", "wall_zlo");
            float wall_zhi = read_attribute<float>(grp, ".", "wall_zhi");
            box_x = wall_xhi - wall_xlo;
            box_y = wall_yhi - wall_ylo;
            box_z = wall_zhi - wall_zlo;
        } else {
            box_x = box_y = box_z = 0.f;
        }

        auto n_pair = get_dset_size(2, grp, "pairs")[0];
        
        // Check if we have the optimized format (coefficient_indices) or the original format
        bool optimized_format = false;
        if(attribute_exists(grp, ".", "optimized_format")) {
            optimized_format = read_attribute<int>(grp, ".", "optimized_format") != 0;
        }

        // Declare unique_coeff at function scope so it's available later
        std::vector<array<float,4>> unique_coeff;

        if(optimized_format) {
            // Load unique coefficients
            auto n_unique_coeff = get_dset_size(2, grp, "coefficients")[0];
            check_size(grp, "coefficient_indices", n_pair);

            unique_coeff.resize(n_unique_coeff);
            traverse_dset<2,float>(grp, "coefficients", [&](size_t nc, size_t d, float x) {
                unique_coeff[nc][d] = x;
            });
            
            // Load coefficient indices
            std::vector<long> coeff_indices;
            coeff_indices.resize(n_pair);
            traverse_dset<1,long>(grp, "coefficient_indices", [&](size_t np, long x) {
                coeff_indices[np] = x;
            });
            
            // Load pairs
            pairs.resize(n_pair);
            traverse_dset<2,int>(grp, "pairs", [&](size_t np, size_t d, int x) {
                if(d == 0) pairs[np].first = x;
                else pairs[np].second = x;
            });

            pair_param_index.resize(n_pair);
            for(size_t np = 0; np < n_pair; ++np) {
                long idx = coeff_indices[np];
                if(idx >= 0 && idx < (long)n_unique_coeff) {
                    pair_param_index[np] = uint32_t(idx);
                } else {
                    throw string("Invalid coefficient index in martini_potential");
                }
            }
        } else {
            check_size(grp, "coefficients", n_pair, 4);
            
            pairs.resize(n_pair);
            traverse_dset<2,int>(grp, "pairs", [&](size_t np, size_t d, int x) {
                if(d == 0) pairs[np].first = x;
                else pairs[np].second = x;
            });

            vector<array<float,4>> coeff_rows(n_pair);
            traverse_dset<2,float>(grp, "coefficients", [&](size_t np, size_t d, float x) {
                coeff_rows[np][d] = x;
            });

            std::map<array<float,4>, uint32_t> coeff_to_index;
            pair_param_index.resize(n_pair);
            for(size_t np = 0; np < n_pair; ++np) {
                auto it = coeff_to_index.find(coeff_rows[np]);
                if(it == coeff_to_index.end()) {
                    uint32_t idx = uint32_t(unique_coeff.size());
                    coeff_to_index.insert(std::make_pair(coeff_rows[np], idx));
                    unique_coeff.push_back(coeff_rows[np]);
                    pair_param_index[np] = idx;
                } else {
                    pair_param_index[np] = it->second;
                }
            }
        }
        
        // Find all epsilon/sigma pairs for separate LJ splines
        std::set<std::pair<float, float>> unique_lj_params;
        for(const auto& c : unique_coeff) {
            float eps = c[0];
            float sig = c[1];
            if(eps != 0.f && sig != 0.f) {
                unique_lj_params.insert({eps, sig});
            }
        }

        // Combined spline domain [0, 12] Angstroms, 1000 points
        r_min = 0.0f;
        r_max = 12.0f;
        r_shift = -r_min;
        r_scale = 999.0f / (r_max + r_shift);

        // Load combined LJ+Coulomb energy grids from injected per-run martini table
        {
            auto dims = get_dset_size(2, grp, "combined_energy_grids");
            size_t n_triples = dims[0];
            size_t n_pts = dims[1];
            if (n_pts != 1000)
                throw string("combined_energy_grids has " + to_string(n_pts) + " points, expected 1000");

            check_size(grp, "unique_eps_eup", n_triples);
            check_size(grp, "unique_sig_ang", n_triples);
            check_size(grp, "unique_charge_product", n_triples);

            std::vector<float> h5_eps(n_triples), h5_sig(n_triples), h5_qq(n_triples);
            traverse_dset<1,double>(grp, "unique_eps_eup", [&](size_t i, double x) { h5_eps[i] = (float)x; });
            traverse_dset<1,double>(grp, "unique_sig_ang", [&](size_t i, double x) { h5_sig[i] = (float)x; });
            traverse_dset<1,double>(grp, "unique_charge_product", [&](size_t i, double x) { h5_qq[i] = (float)x; });

            std::vector<double> all_grids(n_triples * 1000);
            traverse_dset<2,double>(grp, "combined_energy_grids", [&](size_t row, size_t col, double x) {
                all_grids[row * 1000 + col] = x;
            });

            std::map<std::tuple<float,float,float>, size_t> h5_index;
            for (size_t i = 0; i < n_triples; ++i)
                h5_index[{h5_eps[i], h5_sig[i], h5_qq[i]}] = i;

            for (const auto& c : unique_coeff) {
                float eps = c[0];
                float sig = c[1];
                float qq  = c[2] * c[3];
                auto key = std::make_tuple(eps, sig, qq);
                auto it = h5_index.find(key);
                if (it == h5_index.end()) {
                    // Fallback: nearest match within tolerance
                    for (size_t i = 0; i < n_triples; ++i) {
                        if (std::fabs(h5_eps[i] - eps) < 1e-5f &&
                            std::fabs(h5_sig[i] - sig) < 1e-3f &&
                            std::fabs(h5_qq[i] - qq) < 1e-5f) {
                            it = h5_index.find({h5_eps[i], h5_sig[i], h5_qq[i]});
                            break;
                        }
                    }
                }
                if (it == h5_index.end())
                    throw string("Combined params (" + to_string(eps) + ", " + to_string(sig) + ", " + to_string(qq) + ") not found in martini table");

                size_t row = it->second;
                auto inserted = combined_splines.emplace(std::piecewise_construct,
                                                        std::forward_as_tuple(eps, sig, qq),
                                                        std::forward_as_tuple(LayeredClampedSpline1D<1>(1, 1000)));
                if (inserted.second)
                    inserted.first->second.fit_spline(all_grids.data() + row * 1000);
            }
        }

        param_table.resize(unique_coeff.size());
        for(size_t ip = 0; ip < unique_coeff.size(); ++ip) {
            auto& param = param_table[ip];
            param.eps = unique_coeff[ip][0];
            param.sig = unique_coeff[ip][1];
            param.qi  = unique_coeff[ip][2];
            param.qj  = unique_coeff[ip][3];
            param.qq  = param.qi * param.qj;
            auto it = combined_splines.find(std::make_tuple(param.eps, param.sig, param.qq));
            if(it != combined_splines.end()) param.combined_spline = &it->second;
        }

    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("martini_potential"));
        
        VecArray pos1      = pos.output;
        VecArray pos1_sens = pos.sens;
        auto hybrid_state = martini_hybrid::get_state_for_coord(pos);
        auto* mutable_hybrid = static_cast<martini_hybrid::HybridRuntimeState*>(nullptr);
        if(hybrid_state) {
            martini_hybrid::refresh_bb_positions_if_active(*hybrid_state, pos1, n_atom);
            mutable_hybrid = const_cast<martini_hybrid::HybridRuntimeState*>(hybrid_state.get());
            if(mutable_hybrid) {
                mutable_hybrid->bb_env_interface_potential = 0.f;
                martini_hybrid::initialize_sc_env_po4_z_reference(*mutable_hybrid, pos1, n_atom);
            }
        }

        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;

        // Compute particle-particle interactions
        const float kMinDistance = 1.0e-6f;
        const bool active_hybrid_startup = (
            hybrid_state &&
            hybrid_state->enabled &&
            hybrid_state->active);
        float sc_force_uncap_mix = 1.f;
        float sc_backbone_feedback_mix = 1.f;
        if(active_hybrid_startup && mutable_hybrid) {
            sc_force_uncap_mix = martini_hybrid::compute_sc_force_uncap_mix(*mutable_hybrid);
            sc_backbone_feedback_mix = martini_hybrid::compute_sc_backbone_feedback_mix(*mutable_hybrid);
        }

        auto cap_force_vector = [&](Vec<3>& f, float cap_mag) {
            if(!(cap_mag > 0.f)) return;
            float f2 = mag2(f);
            float c2 = cap_mag * cap_mag;
            if(f2 > c2 && f2 > 0.f) {
                f *= (cap_mag / sqrtf(f2));
            }
        };

        auto eval_pair_force = [&](const Vec<3>& pa,
                                   const Vec<3>& pb,
                                   const PairParam& param,
                                   float interaction_scale,
                                   float& pair_potential,
                                   Vec<3>& pair_force,
                                   float force_cap_mag,
                                   float capped_to_regular_mix) -> bool {
            auto dr = pa - pb;
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f) {
                dr = simulation_box::minimum_image(dr, box_x, box_y, box_z);
            }
            auto dist2 = mag2(dr);
            auto dist = sqrtf(max(dist2, kMinDistance));
            pair_potential = 0.f;
            pair_force = make_zero<3>();

            float cutoff = max(lj_cutoff, coul_cutoff);
            if(dist > cutoff) return false;

            float cap_mix = std::max(0.f, std::min(capped_to_regular_mix, 1.f));

            if(param.combined_spline) {
                float r_coord = (dist - r_min) / (r_max - r_min) * 999.0f;
                float result[2];
                param.combined_spline->evaluate_value_and_deriv(result, 0, r_coord);
                float pot = result[1];
                float deriv_spline = result[0];
                float coord_scale = 999.0f / (r_max - r_min);
                float dE_dr = deriv_spline * coord_scale;
                float force_mag = -dE_dr;
                if(std::isfinite(pot) && std::isfinite(force_mag)) {
                    Vec<3> force_uncapped = (force_mag/dist) * dr;
                    Vec<3> force = force_uncapped;
                    if(force_cap_mag > 0.f && cap_mix < 1.f) {
                        cap_force_vector(force, force_cap_mag);
                        if(cap_mix > 0.f) {
                            force = ((1.f - cap_mix) * force) + (cap_mix * force_uncapped);
                        }
                    }
                    pot *= interaction_scale;
                    force *= interaction_scale;
                    pair_potential += pot;
                    pair_force += force;
                }
            }
            return !(pair_potential == 0.f && mag2(pair_force) == 0.f);
        };

        if(pairlist_needs_rebuild(pos1)) {
            rebuild_pairlist(pos1);
        }

        for(int32_t active_idx : active_pair_indices) {
            size_t np = size_t(active_idx);
            int i = pairs[np].first;
            int j = pairs[np].second;
            if(hybrid_state && martini_hybrid::skip_pair_if_intra_protein(*hybrid_state, i, j)) {
                continue;
            }

            const auto& param = param_table[pair_param_index[np]];

            auto p1 = load_vec<3>(pos1, i);
            auto p2 = load_vec<3>(pos1, j);
            bool i_is_protein = false;
            bool j_is_protein = false;
            if(hybrid_state && i >= 0 && j >= 0 &&
               i < (int)hybrid_state->protein_membership.size() &&
               j < (int)hybrid_state->protein_membership.size()) {
                i_is_protein = hybrid_state->protein_membership[i] >= 0;
                j_is_protein = hybrid_state->protein_membership[j] >= 0;
            }
            auto i_role = hybrid_state
                              ? martini_hybrid::atom_role_class_at(*hybrid_state, i)
                              : martini_hybrid::ROLE_OTHER;
            auto j_role = hybrid_state
                              ? martini_hybrid::atom_role_class_at(*hybrid_state, j)
                              : martini_hybrid::ROLE_OTHER;
            if(hybrid_state && hybrid_state->enabled && hybrid_state->active) {
                bool i_carrier_env_pair =
                    i_is_protein &&
                    !j_is_protein &&
                    martini_hybrid::atom_is_backbone_carrier_at(*hybrid_state, i) &&
                    martini_hybrid::bb_map_index_for_proxy(*hybrid_state, i) >= 0;
                bool j_carrier_env_pair =
                    j_is_protein &&
                    !i_is_protein &&
                    martini_hybrid::atom_is_backbone_carrier_at(*hybrid_state, j) &&
                    martini_hybrid::bb_map_index_for_proxy(*hybrid_state, j) >= 0;
                if(i_carrier_env_pair || j_carrier_env_pair) {
                    continue;
                }
            }
            if(active_hybrid_startup &&
               ((i_is_protein && i_role == martini_hybrid::ROLE_SC) ||
                (j_is_protein && j_role == martini_hybrid::ROLE_SC))) {
                continue;
            }

            Vec<3> force = make_zero<3>();
            float pair_pot = 0.f;
            float startup_force_cap = 0.f;
            float startup_cap_mix = 1.f;
            if(hybrid_state &&
               martini_hybrid::deterministic_startup_pair_cap_enabled(
                   *hybrid_state, i_is_protein, j_is_protein, i_role, j_role)) {
                startup_force_cap = max(hybrid_state->sc_env_lj_force_cap,
                                        hybrid_state->sc_env_coul_force_cap);
                startup_cap_mix = sc_force_uncap_mix;
            }
            float interface_scale = hybrid_state
                                        ? martini_hybrid::active_interface_interaction_scale(
                                              *hybrid_state, i_is_protein, j_is_protein)
                                        : 1.f;
            if(!eval_pair_force(p1, p2, param, interface_scale, pair_pot, force,
                                startup_force_cap,
                                startup_cap_mix)) {
                continue;
            }
            if(pot) *pot += pair_pot;
            if(mutable_hybrid && (i_is_protein != j_is_protein)) {
                mutable_hybrid->bb_env_interface_potential += pair_pot;
            }

            auto gi = -force;
            auto gj = force;
            if(active_hybrid_startup && i_is_protein && i_role == martini_hybrid::ROLE_BB) {
                gi *= sc_backbone_feedback_mix;
            }
            if(active_hybrid_startup && j_is_protein && j_role == martini_hybrid::ROLE_BB) {
                gj *= sc_backbone_feedback_mix;
            }
            bool i_projected = false;
            bool j_projected = false;
            if(hybrid_state && i_is_protein && i_role == martini_hybrid::ROLE_BB &&
               !martini_hybrid::atom_is_backbone_carrier_at(*hybrid_state, i) &&
               martini_hybrid::bb_map_index_for_proxy(*hybrid_state, i) >= 0) {
                martini_hybrid::project_bb_proxy_gradient_if_active(*hybrid_state, pos1_sens, n_atom, i, gi);
                i_projected = true;
            }
            if(hybrid_state && j_is_protein && j_role == martini_hybrid::ROLE_BB &&
               !martini_hybrid::atom_is_backbone_carrier_at(*hybrid_state, j) &&
               martini_hybrid::bb_map_index_for_proxy(*hybrid_state, j) >= 0) {
                martini_hybrid::project_bb_proxy_gradient_if_active(*hybrid_state, pos1_sens, n_atom, j, gj);
                j_projected = true;
            }
            if(!i_projected) {
                update_vec<3>(pos1_sens, i, gi);
            }
            if(!j_projected) {
                update_vec<3>(pos1_sens, j, gj);
            }
        }

    }
};
static RegisterNodeType<MartiniPotential, 1> martini_potential_node("martini_potential");

struct MartiniScTablePotential : public PotentialNode
{
    int n_cb;
    int n_env;
    int n_restype;
    int n_target;
    int n_grid;
    int n_layer;

    CoordNode& pos;
    CoordNode& cb_pos;

    vector<int> cb_index;
    vector<int> residue_table_index;
    vector<int> env_atom_index;
    vector<int> env_target_index;

    float energy_conversion_kj_per_eup;
    float length_conversion_angstrom_per_nm;
    float box_x;
    float box_y;
    float box_z;

    float grid_start_ang;
    float grid_step_ang;
    float cutoff_ang;
    int n_angle;
    float cos_start;
    float cos_step;
    float cos_end;

    vector<float> radial_table;
    vector<float> angular_table;
    vector<float> radial_left_value;
    vector<float> radial_left_slope;
    vector<float> angular_left_value;
    vector<float> angular_left_slope;
    vector<float> angular_profile_table;

    inline int radial_index(int layer, int grid_idx) const {
        return layer * n_grid + grid_idx;
    }

    inline int profile_index(int layer, int angle_idx) const {
        return layer * n_angle + angle_idx;
    }

    inline void evaluate_component_value_and_deriv(
            float& value,
            float& dVdr,
            const vector<float>& table,
            const vector<float>& left_value,
            const vector<float>& left_slope,
            int layer,
            float dist) const {
        if(dist <= grid_start_ang) {
            value = left_value[layer] + left_slope[layer] * (dist - grid_start_ang);
            dVdr = left_slope[layer];
            return;
        }

        float radial_coord = (dist - grid_start_ang) / grid_step_ang;
        int grid_idx = int(floorf(radial_coord));
        if(grid_idx < 0) grid_idx = 0;
        if(grid_idx > n_grid - 2) grid_idx = n_grid - 2;
        float frac = radial_coord - float(grid_idx);
        float e0 = table[radial_index(layer, grid_idx)];
        float e1 = table[radial_index(layer, grid_idx + 1)];
        value = (1.f - frac) * e0 + frac * e1;
        dVdr = (e1 - e0) / grid_step_ang;
    }

    inline void evaluate_angular_profile_and_deriv(
            float& value,
            float& dVdcoord,
            int layer,
            float angular_coord) const {
        if(n_angle <= 1) {
            value = angular_profile_table[profile_index(layer, 0)];
            dVdcoord = 0.f;
            return;
        }

        if(angular_coord <= cos_start) {
            value = angular_profile_table[profile_index(layer, 0)];
            dVdcoord = 0.f;
            return;
        }
        if(angular_coord >= cos_end) {
            value = angular_profile_table[profile_index(layer, n_angle - 1)];
            dVdcoord = 0.f;
            return;
        }

        float angle_coord = (angular_coord - cos_start) / cos_step;
        int angle_idx = int(floorf(angle_coord));
        if(angle_idx < 0) angle_idx = 0;
        if(angle_idx > n_angle - 2) angle_idx = n_angle - 2;
        float frac = angle_coord - float(angle_idx);
        float value_lo = angular_profile_table[profile_index(layer, angle_idx)];
        float value_hi = angular_profile_table[profile_index(layer, angle_idx + 1)];
        value = (1.f - frac) * value_lo + frac * value_hi;
        dVdcoord = (value_hi - value_lo) / cos_step;
    }

    inline void evaluate_angular_profile_and_deriv(
            float& value,
            float& dVdcoord,
            int layer,
            float angular_coord,
            const vector<float>& profile_table) const {
        if(n_angle <= 1) {
            value = profile_table[profile_index(layer, 0)];
            dVdcoord = 0.f;
            return;
        }
        if(angular_coord <= cos_start) {
            value = profile_table[profile_index(layer, 0)];
            dVdcoord = 0.f;
            return;
        }
        if(angular_coord >= cos_end) {
            value = profile_table[profile_index(layer, n_angle - 1)];
            dVdcoord = 0.f;
            return;
        }
        float angle_coord = (angular_coord - cos_start) / cos_step;
        int angle_idx = int(floorf(angle_coord));
        if(angle_idx < 0) angle_idx = 0;
        if(angle_idx > n_angle - 2) angle_idx = n_angle - 2;
        float frac = angle_coord - float(angle_idx);
        float value_lo = profile_table[profile_index(layer, angle_idx)];
        float value_hi = profile_table[profile_index(layer, angle_idx + 1)];
        value = (1.f - frac) * value_lo + frac * value_hi;
        dVdcoord = (value_hi - value_lo) / cos_step;
    }

    MartiniScTablePotential(hid_t grp, CoordNode& pos_, CoordNode& cb_pos_):
        PotentialNode(),
        n_cb(get_dset_size(1, grp, "cb_index")[0]),
        n_env(get_dset_size(1, grp, "env_atom_index")[0]),
        n_restype(get_dset_size(3, grp, "radial_energy_kj_mol")[0]),
        n_target(get_dset_size(3, grp, "radial_energy_kj_mol")[1]),
        n_grid(get_dset_size(3, grp, "radial_energy_kj_mol")[2]),
        n_layer(n_restype * n_target),
        pos(pos_),
        cb_pos(cb_pos_),
        cb_index(n_cb),
        residue_table_index(n_cb),
        env_atom_index(n_env),
        env_target_index(n_env),
        energy_conversion_kj_per_eup(read_attribute<float>(grp, ".", "energy_conversion_kj_per_eup")),
        length_conversion_angstrom_per_nm(read_attribute<float>(grp, ".", "length_conversion_angstrom_per_nm")),
        box_x(read_attribute<float>(grp, ".", "x_len")),
        box_y(read_attribute<float>(grp, ".", "y_len")),
        box_z(read_attribute<float>(grp, ".", "z_len")),
        grid_start_ang(0.f),
        grid_step_ang(0.f),
        cutoff_ang(0.f),
        n_angle(get_dset_size(3, grp, "angular_profile")[2]),
        cos_start(0.f),
        cos_step(0.f),
        cos_end(0.f),
        radial_table(n_layer * n_grid, 0.f),
        angular_table(n_layer * n_grid, 0.f),
        radial_left_value(n_layer, 0.f),
        radial_left_slope(n_layer, 0.f),
        angular_left_value(n_layer, 0.f),
        angular_left_slope(n_layer, 0.f),
        angular_profile_table(n_layer * n_angle, 0.f)
    {
        check_elem_width_lower_bound(pos, 3);
        check_elem_width_lower_bound(cb_pos, 6);

        check_size(grp, "residue_table_index", n_cb);
        check_size(grp, "env_target_index", n_env);
        check_size(grp, "grid_nm", n_grid);
        check_size(grp, "cos_theta_grid", n_angle);
        check_size(grp, "angular_energy_kj_mol", n_restype, n_target, n_grid);
        check_size(grp, "angular_profile", n_restype, n_target, n_angle);

        if(!(energy_conversion_kj_per_eup > 0.f) || !(length_conversion_angstrom_per_nm > 0.f)) {
            throw string("martini_sc_table_potential unit-conversion attrs must be positive");
        }
        if(n_restype <= 0 || n_target <= 0 || n_grid < 2 || n_angle < 1) {
            throw string("martini_sc_table_potential requires non-empty residue/target/grid dimensions");
        }

        traverse_dset<1,int>(grp, "cb_index", [&](size_t i, int x) { cb_index[i] = x; });
        traverse_dset<1,int>(grp, "residue_table_index", [&](size_t i, int x) { residue_table_index[i] = x; });
        traverse_dset<1,int>(grp, "env_atom_index", [&](size_t i, int x) { env_atom_index[i] = x; });
        traverse_dset<1,int>(grp, "env_target_index", [&](size_t i, int x) { env_target_index[i] = x; });

        vector<float> grid_nm(n_grid, 0.f);
        traverse_dset<1,float>(grp, "grid_nm", [&](size_t i, float x) { grid_nm[i] = x; });
        vector<float> cos_theta_grid(n_angle, 0.f);
        traverse_dset<1,float>(grp, "cos_theta_grid", [&](size_t i, float x) { cos_theta_grid[i] = x; });

        float grid_step_nm = grid_nm[1] - grid_nm[0];
        if(!(grid_step_nm > 0.f)) {
            throw string("martini_sc_table_potential grid_nm must be strictly increasing");
        }
        for(int i = 2; i < n_grid; ++i) {
            float step = grid_nm[i] - grid_nm[i-1];
            if(fabsf(step - grid_step_nm) > 1e-4f * std::max(1.f, fabsf(grid_step_nm))) {
                throw string("martini_sc_table_potential requires a uniform radial grid");
            }
        }
        grid_start_ang = grid_nm[0] * length_conversion_angstrom_per_nm;
        grid_step_ang = grid_step_nm * length_conversion_angstrom_per_nm;
        cutoff_ang = grid_nm[n_grid-1] * length_conversion_angstrom_per_nm;
        if(!(grid_step_ang > 0.f) || !(cutoff_ang > grid_start_ang)) {
            throw string("martini_sc_table_potential converted radial grid is invalid");
        }

        if(n_angle > 1) {
            cos_step = cos_theta_grid[1] - cos_theta_grid[0];
            if(!(cos_step > 0.f)) {
                throw string("martini_sc_table_potential cos_theta_grid must be strictly increasing");
            }
            for(int i = 2; i < n_angle; ++i) {
                float step = cos_theta_grid[i] - cos_theta_grid[i-1];
                if(fabsf(step - cos_step) > 1e-4f * std::max(1.f, fabsf(cos_step))) {
                    throw string("martini_sc_table_potential requires a uniform cos_theta_grid");
                }
            }
            cos_start = cos_theta_grid[0];
            cos_end = cos_theta_grid[n_angle - 1];
        } else {
            cos_start = cos_end = cos_theta_grid[0];
            cos_step = 1.f;
        }

        for(int i = 0; i < n_cb; ++i) {
            if(cb_index[i] < 0 || cb_index[i] >= cb_pos.n_elem) {
                throw string("martini_sc_table_potential cb_index out of bounds");
            }
            if(residue_table_index[i] < 0 || residue_table_index[i] >= n_restype) {
                throw string("martini_sc_table_potential residue_table_index out of bounds");
            }
        }
        for(int i = 0; i < n_env; ++i) {
            if(env_atom_index[i] < 0 || env_atom_index[i] >= pos.n_elem) {
                throw string("martini_sc_table_potential env_atom_index out of bounds");
            }
            if(env_target_index[i] < 0 || env_target_index[i] >= n_target) {
                throw string("martini_sc_table_potential env_target_index out of bounds");
            }
        }

        vector<float> radial_native(n_layer * n_grid, 0.f);
        vector<float> angular_native(n_layer * n_grid, 0.f);
        traverse_dset<3,float>(grp, "radial_energy_kj_mol", [&](size_t ir, size_t it, size_t ig, float x) {
            radial_native[radial_index(int(ir * n_target + it), int(ig))] = x;
        });
        traverse_dset<3,float>(grp, "angular_energy_kj_mol", [&](size_t ir, size_t it, size_t ig, float x) {
            angular_native[radial_index(int(ir * n_target + it), int(ig))] = x;
        });
        traverse_dset<3,float>(grp, "angular_profile", [&](size_t ir, size_t it, size_t ia, float x) {
            angular_profile_table[profile_index(int(ir * n_target + it), int(ia))] = x;
        });

        for(int layer = 0; layer < n_layer; ++layer) {
            float radial_tail = radial_native[radial_index(layer, n_grid - 1)];
            float angular_tail = angular_native[radial_index(layer, n_grid - 1)];
            for(int ig = 0; ig < n_grid; ++ig) {
                radial_table[radial_index(layer, ig)] =
                    (radial_native[radial_index(layer, ig)] - radial_tail) / energy_conversion_kj_per_eup;
                angular_table[radial_index(layer, ig)] =
                    (angular_native[radial_index(layer, ig)] - angular_tail) / energy_conversion_kj_per_eup;
            }
            radial_left_value[layer] = radial_table[radial_index(layer, 0)];
            radial_left_slope[layer] =
                (radial_table[radial_index(layer, 1)] - radial_table[radial_index(layer, 0)]) / grid_step_ang;
            angular_left_value[layer] = angular_table[radial_index(layer, 0)];
            angular_left_slope[layer] =
                (angular_table[radial_index(layer, 1)] - angular_table[radial_index(layer, 0)]) / grid_step_ang;
        }

    }

    virtual void update_box_dimensions_anisotropic(float scale_xy, float scale_z) override {
        box_x *= scale_xy;
        box_y *= scale_xy;
        box_z *= scale_z;
    }

    virtual void compute_value(ComputeMode mode) override {
        (void)mode;
        Timer timer(string("martini_sc_table_potential"));
        potential = 0.f;

        auto hybrid_state = martini_hybrid::get_state_for_coord(pos);
        if(!hybrid_state || !hybrid_state->enabled || !hybrid_state->active) return;
        float interface_scale = hybrid_state->protein_env_interface_scale;
        float protein_feedback_mix = martini_hybrid::compute_sc_backbone_feedback_mix(*hybrid_state);

        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;
        VecArray cbc = cb_pos.output;
        VecArray cb_sens = cb_pos.sens;
        constexpr float kMinDistance = 1.0e-8f;

        for(int icb = 0; icb < n_cb; ++icb) {
            int cb_idx = cb_index[icb];
            int residue_idx = residue_table_index[icb];
            Vec<6> cb_site = load_vec<6>(cbc, cb_idx);
            Vec<3> cbp = extract<0,3>(cb_site);
            Vec<3> cbv = extract<3,6>(cb_site);

            for(int ienv = 0; ienv < n_env; ++ienv) {
                int atom_idx = env_atom_index[ienv];
                int target_idx = env_target_index[ienv];
                int layer = residue_idx * n_target + target_idx;
                if(layer < 0 || layer >= n_layer) continue;

                Vec<3> envp = load_vec<3>(posc, atom_idx);
                Vec<3> dr = cbp - envp;
                if(box_x > 0.f && box_y > 0.f && box_z > 0.f) {
                    dr = simulation_box::minimum_image(dr, box_x, box_y, box_z);
                }

                float dist2 = mag2(dr);
                float dist = sqrtf(std::max(dist2, kMinDistance));
                if(dist >= cutoff_ang) continue;

                Vec<3> displace_unitvec = (1.f / dist) * dr;
                float cos_theta = dot(displace_unitvec, cbv);
                float angular_coord = -cos_theta;
                float value, dVdr, dVdcoord;

                float radial_value = 0.f, radial_dVdr = 0.f;
                float angular_value = 0.f, angular_dVdr = 0.f;
                float ang1_value = 0.f, dAng1dcoord = 0.f;
                evaluate_component_value_and_deriv(
                    radial_value, radial_dVdr, radial_table, radial_left_value, radial_left_slope, layer, dist);
                evaluate_component_value_and_deriv(
                    angular_value, angular_dVdr, angular_table, angular_left_value, angular_left_slope, layer, dist);
                evaluate_angular_profile_and_deriv(ang1_value, dAng1dcoord, layer, angular_coord);
                value = radial_value + ang1_value * angular_value;
                dVdr = radial_dVdr + ang1_value * angular_dVdr;
                dVdcoord = dAng1dcoord * angular_value;

                if(!std::isfinite(value) || !std::isfinite(dVdr) || !std::isfinite(dVdcoord)) continue;
                value *= interface_scale;
                dVdr *= interface_scale;
                dVdcoord *= interface_scale;
                potential += value;

                if(dist <= 1.0e-6f) continue;
                Vec<3> point_grad = dVdr * displace_unitvec +
                                    (-dVdcoord / dist) * (cbv - cos_theta * displace_unitvec);
                Vec<3> vector_grad = -dVdcoord * displace_unitvec;
                martini_hybrid::cap_force_vector(point_grad, hybrid_state->sc_env_lj_force_cap);
                martini_hybrid::cap_force_vector(vector_grad, hybrid_state->sc_env_lj_force_cap);
                Vec<6> grad_cb_full;
                store<0,3>(grad_cb_full, point_grad);
                store<3,6>(grad_cb_full, vector_grad);
                Vec<6> grad_cb = protein_feedback_mix * grad_cb_full;
                Vec<3> grad_env = -point_grad;

                update_vec<6>(cb_sens, cb_idx, grad_cb);
                update_vec<3>(pos_sens, atom_idx, grad_env);
            }
        }
    }
};
static RegisterNodeType<MartiniScTablePotential, 2> martini_sc_table_potential_node("martini_sc_table_potential");

struct MartiniScTableOneBody : public CoordNode
{
    int n_row;
    int n_env;
    int n_restype;
    int n_rotamer_max;
    int n_target;
    int n_grid;
    int n_layer;

    CoordNode& pos;
    CoordNode& cb_pos;

    vector<int> row_residue_index;
    vector<int> row_rotamer_index;
    vector<int> row_residue_table_index;
    vector<int> row_group_count;
    vector<int> env_atom_index;
    vector<int> env_target_index;
    vector<int> rotamer_count;

    float energy_conversion_kj_per_eup;
    float length_conversion_angstrom_per_nm;
    float box_x;
    float box_y;
    float box_z;

    float grid_start_ang;
    float grid_step_ang;
    float cutoff_ang;
    int n_angle;
    float cos_start;
    float cos_step;
    float cos_end;

    vector<float> radial_table;
    vector<float> angular_table;
    vector<float> radial_left_value;
    vector<float> radial_left_slope;
    vector<float> angular_left_value;
    vector<float> angular_left_slope;
    vector<float> angular_profile_table;

    struct ActiveContact {
        int row;
        int env;
    };

    float cache_buffer;
    bool active_contacts_valid;
    float cached_box_x;
    float cached_box_y;
    float cached_box_z;
    vector<float> cached_cb_point;
    vector<float> cached_env_pos;
    vector<ActiveContact> active_contacts;

    inline int radial_index(int layer, int grid_idx) const {
        return layer * n_grid + grid_idx;
    }

    inline int profile_index(int layer, int angle_idx) const {
        return layer * n_angle + angle_idx;
    }

    inline int layer_index(int residue_idx, int rotamer_idx, int target_idx) const {
        return ((residue_idx * n_rotamer_max) + rotamer_idx) * n_target + target_idx;
    }

    inline bool active_contacts_need_rebuild(const VecArray& posc, const VecArray& cbc) const {
        if(!active_contacts_valid) return true;
        if(cache_buffer <= 0.f) return true;
        if(fabsf(box_x - cached_box_x) > 1.0e-6f ||
           fabsf(box_y - cached_box_y) > 1.0e-6f ||
           fabsf(box_z - cached_box_z) > 1.0e-6f) return true;

        float max_cache_dist2 = sqr(0.5f * cache_buffer);

        for(int icb = 0; icb < cb_pos.n_elem; ++icb) {
            float dx = cbc(0, icb) - cached_cb_point[3*icb + 0];
            float dy = cbc(1, icb) - cached_cb_point[3*icb + 1];
            float dz = cbc(2, icb) - cached_cb_point[3*icb + 2];
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f) {
                simulation_box::minimum_image_scalar(dx, dy, dz, box_x, box_y, box_z);
            }
            if(dx*dx + dy*dy + dz*dz > max_cache_dist2) return true;
        }

        for(int ienv = 0; ienv < n_env; ++ienv) {
            int atom_idx = env_atom_index[ienv];
            float dx = posc(0, atom_idx) - cached_env_pos[3*ienv + 0];
            float dy = posc(1, atom_idx) - cached_env_pos[3*ienv + 1];
            float dz = posc(2, atom_idx) - cached_env_pos[3*ienv + 2];
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f) {
                simulation_box::minimum_image_scalar(dx, dy, dz, box_x, box_y, box_z);
            }
            if(dx*dx + dy*dy + dz*dz > max_cache_dist2) return true;
        }

        return false;
    }

    void rebuild_active_contacts(const VecArray& posc, const VecArray& cbc) {
        Timer timer(string("martini_sc_env_pairlist_rebuild"));
        float active_cutoff = cutoff_ang + std::max(0.f, cache_buffer);
        float active_cutoff2 = sqr(active_cutoff);

        active_contacts.clear();
        active_contacts.reserve(std::min<size_t>(size_t(n_row) * size_t(n_env), size_t(1u << 18)));

        for(int irow = 0; irow < n_row; ++irow) {
            int cb_idx = row_residue_index[irow];
            Vec<3> cbp = extract<0,3>(load_vec<6>(cbc, cb_idx));

            for(int ienv = 0; ienv < n_env; ++ienv) {
                int atom_idx = env_atom_index[ienv];
                Vec<3> envp = load_vec<3>(posc, atom_idx);
                Vec<3> dr = cbp - envp;
                if(box_x > 0.f && box_y > 0.f && box_z > 0.f) {
                    dr = simulation_box::minimum_image(dr, box_x, box_y, box_z);
                }
                if(mag2(dr) < active_cutoff2) {
                    active_contacts.push_back({irow, ienv});
                }
            }
        }

        for(int icb = 0; icb < cb_pos.n_elem; ++icb) {
            cached_cb_point[3*icb + 0] = cbc(0, icb);
            cached_cb_point[3*icb + 1] = cbc(1, icb);
            cached_cb_point[3*icb + 2] = cbc(2, icb);
        }
        for(int ienv = 0; ienv < n_env; ++ienv) {
            int atom_idx = env_atom_index[ienv];
            cached_env_pos[3*ienv + 0] = posc(0, atom_idx);
            cached_env_pos[3*ienv + 1] = posc(1, atom_idx);
            cached_env_pos[3*ienv + 2] = posc(2, atom_idx);
        }
        cached_box_x = box_x;
        cached_box_y = box_y;
        cached_box_z = box_z;
        active_contacts_valid = true;
    }

    inline void evaluate_component_value_and_deriv(
            float& value,
            float& dVdr,
            const vector<float>& table,
            const vector<float>& left_value,
            const vector<float>& left_slope,
            int layer,
            float dist) const {
        if(dist <= grid_start_ang) {
            value = left_value[layer] + left_slope[layer] * (dist - grid_start_ang);
            dVdr = left_slope[layer];
            return;
        }

        float radial_coord = (dist - grid_start_ang) / grid_step_ang;
        int grid_idx = int(floorf(radial_coord));
        if(grid_idx < 0) grid_idx = 0;
        if(grid_idx > n_grid - 2) grid_idx = n_grid - 2;
        float frac = radial_coord - float(grid_idx);
        float e0 = table[radial_index(layer, grid_idx)];
        float e1 = table[radial_index(layer, grid_idx + 1)];
        value = (1.f - frac) * e0 + frac * e1;
        dVdr = (e1 - e0) / grid_step_ang;
    }

    inline void evaluate_angular_profile_and_deriv(
            float& value,
            float& dVdcoord,
            int layer,
            float angular_coord) const {
        if(n_angle <= 1) {
            value = angular_profile_table[profile_index(layer, 0)];
            dVdcoord = 0.f;
            return;
        }

        if(angular_coord <= cos_start) {
            value = angular_profile_table[profile_index(layer, 0)];
            dVdcoord = 0.f;
            return;
        }
        if(angular_coord >= cos_end) {
            value = angular_profile_table[profile_index(layer, n_angle - 1)];
            dVdcoord = 0.f;
            return;
        }

        float angle_coord = (angular_coord - cos_start) / cos_step;
        int angle_idx = int(floorf(angle_coord));
        if(angle_idx < 0) angle_idx = 0;
        if(angle_idx > n_angle - 2) angle_idx = n_angle - 2;
        float frac = angle_coord - float(angle_idx);
        float value_lo = angular_profile_table[profile_index(layer, angle_idx)];
        float value_hi = angular_profile_table[profile_index(layer, angle_idx + 1)];
        value = (1.f - frac) * value_lo + frac * value_hi;
        dVdcoord = (value_hi - value_lo) / cos_step;
    }

    inline void evaluate_angular_profile_and_deriv(
            float& value,
            float& dVdcoord,
            int layer,
            float angular_coord,
            const vector<float>& profile_table) const {
        if(n_angle <= 1) {
            value = profile_table[profile_index(layer, 0)];
            dVdcoord = 0.f;
            return;
        }
        if(angular_coord <= cos_start) {
            value = profile_table[profile_index(layer, 0)];
            dVdcoord = 0.f;
            return;
        }
        if(angular_coord >= cos_end) {
            value = profile_table[profile_index(layer, n_angle - 1)];
            dVdcoord = 0.f;
            return;
        }
        float angle_coord = (angular_coord - cos_start) / cos_step;
        int angle_idx = int(floorf(angle_coord));
        if(angle_idx < 0) angle_idx = 0;
        if(angle_idx > n_angle - 2) angle_idx = n_angle - 2;
        float frac = angle_coord - float(angle_idx);
        float value_lo = profile_table[profile_index(layer, angle_idx)];
        float value_hi = profile_table[profile_index(layer, angle_idx + 1)];
        value = (1.f - frac) * value_lo + frac * value_hi;
        dVdcoord = (value_hi - value_lo) / cos_step;
    }

    MartiniScTableOneBody(hid_t grp, CoordNode& pos_, CoordNode& cb_pos_):
        CoordNode(get_dset_size(1, grp, "row_residue_index")[0], 1),
        n_row(get_dset_size(1, grp, "row_residue_index")[0]),
        n_env(get_dset_size(1, grp, "env_atom_index")[0]),
        n_restype(get_dset_size(4, grp, "rotamer_radial_energy_kj_mol")[0]),
        n_rotamer_max(get_dset_size(4, grp, "rotamer_radial_energy_kj_mol")[1]),
        n_target(get_dset_size(4, grp, "rotamer_radial_energy_kj_mol")[2]),
        n_grid(get_dset_size(4, grp, "rotamer_radial_energy_kj_mol")[3]),
        n_layer(n_restype * n_rotamer_max * n_target),
        pos(pos_),
        cb_pos(cb_pos_),
        row_residue_index(n_row),
        row_rotamer_index(n_row),
        row_residue_table_index(n_row),
        row_group_count(n_row, 1),
        env_atom_index(n_env),
        env_target_index(n_env),
        rotamer_count(n_restype, 0),
        energy_conversion_kj_per_eup(read_attribute<float>(grp, ".", "energy_conversion_kj_per_eup")),
        length_conversion_angstrom_per_nm(read_attribute<float>(grp, ".", "length_conversion_angstrom_per_nm")),
        box_x(read_attribute<float>(grp, ".", "x_len")),
        box_y(read_attribute<float>(grp, ".", "y_len")),
        box_z(read_attribute<float>(grp, ".", "z_len")),
        grid_start_ang(0.f),
        grid_step_ang(0.f),
        cutoff_ang(0.f),
        n_angle(get_dset_size(4, grp, "rotamer_angular_profile")[3]),
        cos_start(0.f),
        cos_step(0.f),
        cos_end(0.f),
        radial_table(n_layer * n_grid, 0.f),
        angular_table(n_layer * n_grid, 0.f),
        radial_left_value(n_layer, 0.f),
        radial_left_slope(n_layer, 0.f),
        angular_left_value(n_layer, 0.f),
        angular_left_slope(n_layer, 0.f),
        angular_profile_table(n_layer * n_angle, 0.f),
        cache_buffer(read_attribute<float>(grp, ".", "cache_buffer", 1.f)),
        active_contacts_valid(false),
        cached_box_x(0.f),
        cached_box_y(0.f),
        cached_box_z(0.f),
        cached_cb_point(3*cb_pos_.n_elem, 0.f),
        cached_env_pos(3*n_env, 0.f)
    {
        check_elem_width_lower_bound(pos, 3);
        check_elem_width_lower_bound(cb_pos, 6);
        if(cache_buffer < 0.f) cache_buffer = 0.f;

        check_size(grp, "row_rotamer_index", n_row);
        check_size(grp, "row_residue_table_index", n_row);
        check_size(grp, "env_target_index", n_env);
        check_size(grp, "grid_nm", n_grid);
        check_size(grp, "cos_theta_grid", n_angle);
        check_size(grp, "rotamer_count", n_restype);
        check_size(grp, "rotamer_angular_energy_kj_mol", n_restype, n_rotamer_max, n_target, n_grid);
        check_size(grp, "rotamer_angular_profile", n_restype, n_rotamer_max, n_target, n_angle);

        if(!(energy_conversion_kj_per_eup > 0.f) || !(length_conversion_angstrom_per_nm > 0.f)) {
            throw string("martini_sc_table_1body unit-conversion attrs must be positive");
        }
        if(n_restype <= 0 || n_rotamer_max <= 0 || n_target <= 0 || n_grid < 2 || n_angle < 1) {
            throw string("martini_sc_table_1body requires non-empty residue/rotamer/target/grid dimensions");
        }

        traverse_dset<1,int>(grp, "row_residue_index", [&](size_t i, int x) { row_residue_index[i] = x; });
        traverse_dset<1,int>(grp, "row_rotamer_index", [&](size_t i, int x) { row_rotamer_index[i] = x; });
        traverse_dset<1,int>(grp, "row_residue_table_index", [&](size_t i, int x) { row_residue_table_index[i] = x; });
        traverse_dset<1,int>(grp, "env_atom_index", [&](size_t i, int x) { env_atom_index[i] = x; });
        traverse_dset<1,int>(grp, "env_target_index", [&](size_t i, int x) { env_target_index[i] = x; });
        traverse_dset<1,int>(grp, "rotamer_count", [&](size_t i, int x) { rotamer_count[i] = x; });

        std::unordered_map<uint64_t, int> group_counts;
        group_counts.reserve(static_cast<size_t>(n_row) * 2u);
        for(int i = 0; i < n_row; ++i) {
            uint64_t key = (uint64_t(uint32_t(row_residue_index[i])) << 32) | uint32_t(row_rotamer_index[i]);
            group_counts[key] += 1;
        }
        for(int i = 0; i < n_row; ++i) {
            uint64_t key = (uint64_t(uint32_t(row_residue_index[i])) << 32) | uint32_t(row_rotamer_index[i]);
            auto it = group_counts.find(key);
            row_group_count[i] = (it == group_counts.end() || it->second < 1) ? 1 : it->second;
        }

        vector<float> grid_nm(n_grid, 0.f);
        traverse_dset<1,float>(grp, "grid_nm", [&](size_t i, float x) { grid_nm[i] = x; });
        vector<float> cos_theta_grid(n_angle, 0.f);
        traverse_dset<1,float>(grp, "cos_theta_grid", [&](size_t i, float x) { cos_theta_grid[i] = x; });

        float grid_step_nm = grid_nm[1] - grid_nm[0];
        if(!(grid_step_nm > 0.f)) {
            throw string("martini_sc_table_1body grid_nm must be strictly increasing");
        }
        for(int i = 2; i < n_grid; ++i) {
            float step = grid_nm[i] - grid_nm[i-1];
            if(fabsf(step - grid_step_nm) > 1e-4f * std::max(1.f, fabsf(grid_step_nm))) {
                throw string("martini_sc_table_1body requires a uniform radial grid");
            }
        }
        grid_start_ang = grid_nm[0] * length_conversion_angstrom_per_nm;
        grid_step_ang = grid_step_nm * length_conversion_angstrom_per_nm;
        cutoff_ang = grid_nm[n_grid - 1] * length_conversion_angstrom_per_nm;
        if(!(grid_step_ang > 0.f) || !(cutoff_ang > grid_start_ang)) {
            throw string("martini_sc_table_1body converted radial grid is invalid");
        }

        if(n_angle > 1) {
            cos_step = cos_theta_grid[1] - cos_theta_grid[0];
            if(!(cos_step > 0.f)) {
                throw string("martini_sc_table_1body cos_theta_grid must be strictly increasing");
            }
            for(int i = 2; i < n_angle; ++i) {
                float step = cos_theta_grid[i] - cos_theta_grid[i - 1];
                if(fabsf(step - cos_step) > 1e-4f * std::max(1.f, fabsf(cos_step))) {
                    throw string("martini_sc_table_1body requires a uniform cos_theta_grid");
                }
            }
            cos_start = cos_theta_grid[0];
            cos_end = cos_theta_grid[n_angle - 1];
        } else {
            cos_start = cos_end = cos_theta_grid[0];
            cos_step = 1.f;
        }

        for(int i = 0; i < n_row; ++i) {
            if(row_residue_index[i] < 0 || row_residue_index[i] >= cb_pos.n_elem) {
                throw string("martini_sc_table_1body row_residue_index out of bounds");
            }
            if(row_residue_table_index[i] < 0 || row_residue_table_index[i] >= n_restype) {
                throw string("martini_sc_table_1body row_residue_table_index out of bounds");
            }
            int res_idx = row_residue_table_index[i];
            if(row_rotamer_index[i] < 0 || row_rotamer_index[i] >= rotamer_count[res_idx]) {
                throw string("martini_sc_table_1body row_rotamer_index out of bounds");
            }
        }
        for(int i = 0; i < n_env; ++i) {
            if(env_atom_index[i] < 0 || env_atom_index[i] >= pos.n_elem) {
                throw string("martini_sc_table_1body env_atom_index out of bounds");
            }
            if(env_target_index[i] < 0 || env_target_index[i] >= n_target) {
                throw string("martini_sc_table_1body env_target_index out of bounds");
            }
        }

        vector<float> radial_native(n_layer * n_grid, 0.f);
        vector<float> angular_native(n_layer * n_grid, 0.f);
        traverse_dset<4,float>(grp, "rotamer_radial_energy_kj_mol",
                [&](size_t ir, size_t iro, size_t it, size_t ig, float x) {
                    radial_native[radial_index(layer_index(int(ir), int(iro), int(it)), int(ig))] = x;
                });
        traverse_dset<4,float>(grp, "rotamer_angular_energy_kj_mol",
                [&](size_t ir, size_t iro, size_t it, size_t ig, float x) {
                    angular_native[radial_index(layer_index(int(ir), int(iro), int(it)), int(ig))] = x;
                });
        traverse_dset<4,float>(grp, "rotamer_angular_profile",
                [&](size_t ir, size_t iro, size_t it, size_t ia, float x) {
                    angular_profile_table[profile_index(layer_index(int(ir), int(iro), int(it)), int(ia))] = x;
                });

        for(int layer = 0; layer < n_layer; ++layer) {
            float radial_tail = radial_native[radial_index(layer, n_grid - 1)];
            float angular_tail = angular_native[radial_index(layer, n_grid - 1)];
            for(int ig = 0; ig < n_grid; ++ig) {
                radial_table[radial_index(layer, ig)] =
                    (radial_native[radial_index(layer, ig)] - radial_tail) / energy_conversion_kj_per_eup;
                angular_table[radial_index(layer, ig)] =
                    (angular_native[radial_index(layer, ig)] - angular_tail) / energy_conversion_kj_per_eup;
            }
            radial_left_value[layer] = radial_table[radial_index(layer, 0)];
            radial_left_slope[layer] =
                (radial_table[radial_index(layer, 1)] - radial_table[radial_index(layer, 0)]) / grid_step_ang;
            angular_left_value[layer] = angular_table[radial_index(layer, 0)];
            angular_left_slope[layer] =
                (angular_table[radial_index(layer, 1)] - angular_table[radial_index(layer, 0)]) / grid_step_ang;
        }

    }

    virtual void compute_value(ComputeMode mode) override {
        (void)mode;
        Timer timer(string("martini_sc_table_1body"));
        fill(output, 0.f);

        auto hybrid_state = martini_hybrid::get_state_for_coord(pos);
        if(!hybrid_state || !hybrid_state->enabled || !hybrid_state->active) return;
        float interface_scale = hybrid_state->protein_env_interface_scale;

        VecArray posc = pos.output;
        VecArray cbc = cb_pos.output;
        constexpr float kMinDistance = 1.0e-8f;

        if(active_contacts_need_rebuild(posc, cbc)) {
            rebuild_active_contacts(posc, cbc);
        }

        for(const auto& active : active_contacts) {
            int irow = active.row;
            int ienv = active.env;
            int cb_idx = row_residue_index[irow];
            int residue_idx = row_residue_table_index[irow];
            int rotamer_idx = row_rotamer_index[irow];
            int atom_idx = env_atom_index[ienv];
            int target_idx = env_target_index[ienv];
            int layer = layer_index(residue_idx, rotamer_idx, target_idx);

            Vec<6> cb_site = load_vec<6>(cbc, cb_idx);
            Vec<3> cbp = extract<0,3>(cb_site);
            Vec<3> cbv = extract<3,6>(cb_site);
            Vec<3> envp = load_vec<3>(posc, atom_idx);
            Vec<3> dr = cbp - envp;
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f) {
                dr = simulation_box::minimum_image(dr, box_x, box_y, box_z);
            }

            float dist2 = mag2(dr);
            float dist = sqrtf(std::max(dist2, kMinDistance));
            if(dist >= cutoff_ang) continue;

            Vec<3> displace_unitvec = (1.f / dist) * dr;
            float angular_coord = -dot(displace_unitvec, cbv);

            float radial_value = 0.f, radial_dVdr = 0.f;
            float angular_value = 0.f, angular_dVdr = 0.f;
            float ang1_value = 0.f, dAng1dcoord = 0.f;
            evaluate_component_value_and_deriv(
                radial_value, radial_dVdr, radial_table, radial_left_value, radial_left_slope, layer, dist);
            evaluate_component_value_and_deriv(
                angular_value, angular_dVdr, angular_table, angular_left_value, angular_left_slope, layer, dist);
            evaluate_angular_profile_and_deriv(ang1_value, dAng1dcoord, layer, angular_coord);
            output(0, irow) += radial_value + ang1_value * angular_value;
        }

        for(int irow = 0; irow < n_row; ++irow) {
            int group_count = std::max(1, row_group_count[irow]);
            output(0, irow) = interface_scale * (output(0, irow) / float(group_count));
        }
    }

    virtual void propagate_deriv() override {
        Timer timer(string("d_martini_sc_table_1body"));

        auto hybrid_state = martini_hybrid::get_state_for_coord(pos);
        if(!hybrid_state || !hybrid_state->enabled || !hybrid_state->active) return;
        float interface_scale = hybrid_state->protein_env_interface_scale;
        float protein_feedback_mix = martini_hybrid::compute_sc_backbone_feedback_mix(*hybrid_state);

        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;
        VecArray cbc = cb_pos.output;
        VecArray cb_sens = cb_pos.sens;
        constexpr float kMinDistance = 1.0e-8f;

        if(active_contacts_need_rebuild(posc, cbc)) {
            rebuild_active_contacts(posc, cbc);
        }

        for(const auto& active : active_contacts) {
            int irow = active.row;
            int ienv = active.env;
            float row_scale = sens(0, irow);
            if(row_scale == 0.f) continue;

            int cb_idx = row_residue_index[irow];
            int residue_idx = row_residue_table_index[irow];
            int rotamer_idx = row_rotamer_index[irow];
            int group_count = std::max(1, row_group_count[irow]);
            row_scale = interface_scale * (row_scale / float(group_count));

            Vec<6> cb_site = load_vec<6>(cbc, cb_idx);
            Vec<3> cbp = extract<0,3>(cb_site);
            Vec<3> cbv = extract<3,6>(cb_site);

            int atom_idx = env_atom_index[ienv];
            int target_idx = env_target_index[ienv];
            int layer = layer_index(residue_idx, rotamer_idx, target_idx);

            Vec<3> envp = load_vec<3>(posc, atom_idx);
            Vec<3> dr = cbp - envp;
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f) {
                dr = simulation_box::minimum_image(dr, box_x, box_y, box_z);
            }

            float dist2 = mag2(dr);
            float dist = sqrtf(std::max(dist2, kMinDistance));
            if(dist >= cutoff_ang || dist <= 1.0e-6f) continue;

            Vec<3> displace_unitvec = (1.f / dist) * dr;
            float cos_theta = dot(displace_unitvec, cbv);
            float angular_coord = -cos_theta;
            float dVdr, dVdcoord;

            float radial_value = 0.f, radial_dVdr = 0.f;
            float angular_value = 0.f, angular_dVdr = 0.f;
            float ang1_value = 0.f, dAng1dcoord = 0.f;
            evaluate_component_value_and_deriv(
                radial_value, radial_dVdr, radial_table, radial_left_value, radial_left_slope, layer, dist);
            evaluate_component_value_and_deriv(
                angular_value, angular_dVdr, angular_table, angular_left_value, angular_left_slope, layer, dist);
            evaluate_angular_profile_and_deriv(ang1_value, dAng1dcoord, layer, angular_coord);
            dVdr = radial_dVdr + ang1_value * angular_dVdr;
            dVdcoord = dAng1dcoord * angular_value;
            if(!std::isfinite(dVdr) || !std::isfinite(dVdcoord)) continue;

            Vec<3> point_grad = row_scale * (
                dVdr * displace_unitvec +
                (-dVdcoord / dist) * (cbv - cos_theta * displace_unitvec));
            Vec<3> vector_grad = row_scale * (-dVdcoord * displace_unitvec);
            martini_hybrid::cap_force_vector(point_grad, hybrid_state->sc_env_lj_force_cap);
            martini_hybrid::cap_force_vector(vector_grad, hybrid_state->sc_env_lj_force_cap);
            Vec<6> grad_cb_full;
            store<0,3>(grad_cb_full, point_grad);
            store<3,6>(grad_cb_full, vector_grad);
            Vec<6> grad_cb = protein_feedback_mix * grad_cb_full;
            Vec<3> grad_env = -point_grad;

            update_vec<6>(cb_sens, cb_idx, grad_cb);
            update_vec<3>(pos_sens, atom_idx, grad_env);
        }
    }
};
static RegisterNodeType<MartiniScTableOneBody, 2> martini_sc_table_1body_node("martini_sc_table_1body");

// Bond potential using spline interpolation
struct DistSpring : public PotentialNode
{
    struct Params {
        index_t atom[2];
        float equil_dist;
        float spring_constant;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    vector<int> bonded_atoms;
    
    // Box dimensions for minimum image
    float box_x, box_y, box_z;
    
    // Spline interpolation for bond potential
    LayeredClampedSpline1D<1> bond_potential_spline;
    
    // Spline parameters
    float bond_r_min, bond_r_max;
    // Canonical delta-r spline domain (shared across bonds)
    float bond_delta_min, bond_delta_max;
    float max_spring;  // Store max spring constant for scaling

    DistSpring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(2, grp, "id")[0]), pos(pos_), params(n_elem),
        bond_potential_spline(1, 1000)
    {
        int n_dep = 2;  // number of atoms that each term depends on 
        check_size(grp, "id",           n_elem, n_dep);
        check_size(grp, "equil_dist",   n_elem);
        check_size(grp, "spring_const", n_elem);
        check_size(grp, "bonded_atoms", n_elem);

        auto& p = params;
        traverse_dset<2,int>  (grp, "id",           [&](size_t i, size_t j, int   x) {p[i].atom[j] = x;});
        traverse_dset<1,float>(grp, "equil_dist",   [&](size_t i,           float x) {p[i].equil_dist = x;});
        traverse_dset<1,float>(grp, "spring_const", [&](size_t i,           float x) {p[i].spring_constant = x;});
        traverse_dset<1,int>  (grp, "bonded_atoms", [&](size_t i,           int   x) {bonded_atoms.push_back(x);});

        // Read box dimensions
        if(attribute_exists(grp, ".", "x_len") && attribute_exists(grp, ".", "y_len") && attribute_exists(grp, ".", "z_len")) {
            box_x = read_attribute<float>(grp, ".", "x_len");
            box_y = read_attribute<float>(grp, ".", "y_len");
            box_z = read_attribute<float>(grp, ".", "z_len");
        } else {
            float wall_xlo = read_attribute<float>(grp, ".", "wall_xlo");
            float wall_xhi = read_attribute<float>(grp, ".", "wall_xhi");
            float wall_ylo = read_attribute<float>(grp, ".", "wall_ylo");
            float wall_yhi = read_attribute<float>(grp, ".", "wall_yhi");
            float wall_zlo = read_attribute<float>(grp, ".", "wall_zlo");
            float wall_zhi = read_attribute<float>(grp, ".", "wall_zhi");
            box_x = wall_xhi - wall_xlo;
            box_y = wall_yhi - wall_ylo;
            box_z = wall_zhi - wall_zlo;
        }

        // Initialize spline parameters for bond potential
        // Find the range of equilibrium distances and spring constants
        float min_equil = std::numeric_limits<float>::max();
        float max_equil = std::numeric_limits<float>::lowest();
        float max_spring = std::numeric_limits<float>::lowest();
        
        for(const auto& p : params) {
            min_equil = std::min(min_equil, p.equil_dist);
            max_equil = std::max(max_equil, p.equil_dist);
            max_spring = std::max(max_spring, p.spring_constant);
        }
        
        // Store max_spring as member variable
        this->max_spring = max_spring;
        
        // Define the bond r-range to start at 0 and extend to r_max
        bond_r_min = 0.0f;
        // Prefer an attribute if provided; otherwise choose a conservative default
        if(attribute_exists(grp, ".", "bond_r_max")) {
            bond_r_max = read_attribute<float>(grp, ".", "bond_r_max");
        } else {
            // heuristic default: at least 2x the largest equilibrium distance or 5 A
            bond_r_max = std::max(2.0f * max_equil, 5.0f);
        }
        // Define and store a GLOBAL delta-r domain shared by all bonds corresponding to r in [0, bond_r_max]
        bond_delta_min = bond_r_min - max_equil; // smallest possible (r - r0) with r_min=0
        bond_delta_max = bond_r_max - min_equil; // largest possible (r - r0)

        // Generate canonical spline in delta-r: pot=0.5*delta^2
        std::vector<double> bond_pot_data(1000);
        for(int i = 0; i < 1000; ++i) {
            float delta_r = bond_delta_min + i * (bond_delta_max - bond_delta_min) / 999.0f;
            bond_pot_data[i] = 0.5 * delta_r * delta_r;
        }
        
        // Fit spline
        bond_potential_spline.fit_spline(bond_pot_data.data());

    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("dist_spring"));

        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;
        auto hybrid_state = martini_hybrid::get_state_for_coord(pos);
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;

        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            if(hybrid_state && !martini_hybrid::allow_intra_protein_pair_if_active(*hybrid_state, p.atom[0], p.atom[1])) {
                continue;
            }

            auto x1 = load_vec<3>(posc, p.atom[0]);
            auto x2 = load_vec<3>(posc, p.atom[1]);

            // Apply PBC minimum image convention for distance calculation
            auto disp = simulation_box::minimum_image(x1 - x2, box_x, box_y, box_z);
            float dist = mag(disp);
            
            // Use spline interpolation for bond potential and force (in delta-r space)
            float delta_r = dist - p.equil_dist;
            if(delta_r >= bond_delta_min && delta_r <= bond_delta_max) {
                // Transform delta_r to spline coordinate [0, 999] (same as LJ and Coulomb splines)
                float r_coord = (delta_r - bond_delta_min) / (bond_delta_max - bond_delta_min) * 999.0f;
                
                // Get potential and force from single spline
                float result[2];
                bond_potential_spline.evaluate_value_and_deriv(result, 0, r_coord);
                float bond_pot = result[1]; // Index 1 is the value
                float bond_deriv_spline = result[0]; // Index 0 is the derivative w.r.t. spline coordinate
                
                
                // Convert derivative from spline coordinate to physical coordinate (dE/d(delta_r))
                // dE/d(delta_r) = dE/d(coord) * d(coord)/d(delta_r)
                float coord_scale = 999.0f / (bond_delta_max - bond_delta_min);
                float bond_deriv = bond_deriv_spline * coord_scale;
                
                // Scale by spring constant to match parameterized harmonic
                if(pot) *pot += p.spring_constant * bond_pot;
                
                // Apply force: F = -dE/dr, here dE/dr = k * dE/d(delta_r) * d(delta_r)/dr with d(delta_r)/dr = 1
                // The derivative is now correctly dE/d(delta_r), so we use it directly
                auto deriv = (p.spring_constant * bond_deriv / dist) * disp;
                update_vec(pos_sens, p.atom[0],  deriv);
                update_vec(pos_sens, p.atom[1], -deriv);
            } else {
                // Fallback to direct calculation for out-of-range distances
                auto deriv = p.spring_constant * (1.f - p.equil_dist*inv_mag(disp)) * disp;
                if(pot) *pot += 0.5f * p.spring_constant * sqr(mag(disp) - p.equil_dist);
                update_vec(pos_sens, p.atom[0],  deriv);
                update_vec(pos_sens, p.atom[1], -deriv);
            }
        }
    }
    
};
static RegisterNodeType<DistSpring, 1> dist_spring_node("dist_spring");

// Angle potential using spline interpolation
struct AngleSpring : public PotentialNode
{
    struct Params {
        index_t atom[3];
        float equil_angle_deg; // store equilibrium angle in degrees
        float spring_constant;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    
    // Box dimensions for minimum image
    float box_x, box_y, box_z;
    
    // Spline interpolation for angle potential
    LayeredClampedSpline1D<1> angle_potential_spline;
    
    // Spline parameters
    float angle_cos_min, angle_cos_max;
    float max_spring;  // Store max spring constant for scaling

    AngleSpring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(2, grp, "id")[0]), pos(pos_), params(n_elem),
        angle_potential_spline(1, 1000)
    {
        int n_dep = 3;  // number of atoms that each term depends on 
        check_size(grp, "id",              n_elem, n_dep);
        check_size(grp, "equil_angle_deg", n_elem);
        check_size(grp, "spring_const",    n_elem);

        auto& p = params;
        traverse_dset<2,int>  (grp, "id",              [&](size_t i, size_t j, int   x) { p[i].atom[j] = x;});
        traverse_dset<1,float>(grp, "equil_angle_deg", [&](size_t i,           float x) { p[i].equil_angle_deg = x;});
        traverse_dset<1,float>(grp, "spring_const",    [&](size_t i,           float x) { p[i].spring_constant = x;});

        // Read box dimensions
        if(attribute_exists(grp, ".", "x_len") && attribute_exists(grp, ".", "y_len") && attribute_exists(grp, ".", "z_len")) {
            box_x = read_attribute<float>(grp, ".", "x_len");
            box_y = read_attribute<float>(grp, ".", "y_len");
            box_z = read_attribute<float>(grp, ".", "z_len");
        } else {
            float wall_xlo = read_attribute<float>(grp, ".", "wall_xlo");
            float wall_xhi = read_attribute<float>(grp, ".", "wall_xhi");
            float wall_ylo = read_attribute<float>(grp, ".", "wall_ylo");
            float wall_yhi = read_attribute<float>(grp, ".", "wall_yhi");
            float wall_zlo = read_attribute<float>(grp, ".", "wall_zlo");
            float wall_zhi = read_attribute<float>(grp, ".", "wall_zhi");
            box_x = wall_xhi - wall_xlo;
            box_y = wall_yhi - wall_ylo;
            box_z = wall_zhi - wall_zlo;
        }

        // Initialize spline parameters for angle potential
        // Find the range of equilibrium angles and spring constants
        float min_angle = std::numeric_limits<float>::max();
        float max_angle = std::numeric_limits<float>::lowest();
        float max_spring = std::numeric_limits<float>::lowest();
        
        for(const auto& p : params) {
            min_angle = std::min(min_angle, p.equil_angle_deg);
            max_angle = std::max(max_angle, p.equil_angle_deg);
            max_spring = std::max(max_spring, p.spring_constant);
        }
        
        // Store max_spring as member variable
        this->max_spring = max_spring;
        
        // Set spline range for delta_cos = cos(theta) - cos(theta0)
        // Find the range of cos(equilibrium angles) to set proper delta_cos bounds
        float min_cos_equil = 1.0f;
        float max_cos_equil = -1.0f;
        for(const auto& p : params) {
            float cos_equil = cosf(p.equil_angle_deg * M_PI / 180.0f);
            min_cos_equil = std::min(min_cos_equil, cos_equil);
            max_cos_equil = std::max(max_cos_equil, cos_equil);
        }

        // Set conservative bounds for delta_cos
        angle_cos_min = -2.0f;  // Allow for angles deviating significantly from equilibrium
        angle_cos_max = 2.0f;

        // Canonical spline for angles: pot = 0.5*(delta_cos)^2 (cosine-based harmonic potential)
        std::vector<double> angle_pot_data(1000);
        for(int i = 0; i < 1000; ++i) {
            float delta_cos = angle_cos_min + i * (angle_cos_max - angle_cos_min) / 999.0f;
            angle_pot_data[i] = 0.5 * delta_cos * delta_cos;
        }
        
        // Fit spline
        angle_potential_spline.fit_spline(angle_pot_data.data());

    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("angle_spring"));

        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;
        auto hybrid_state = martini_hybrid::get_state_for_coord(pos);
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;

        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            if(hybrid_state && !martini_hybrid::allow_multibody_term_if_active(*hybrid_state, p.atom, 3)) {
                continue;
            }
            
            // Step 1: Load atomic positions
            auto x_orig1 = load_vec<3>(posc, p.atom[0]);
            auto x_orig2 = load_vec<3>(posc, p.atom[1]);
            auto x_orig3 = load_vec<3>(posc, p.atom[2]);

            auto x2 = x_orig2;
            auto x1 = x_orig1;
            auto x3 = x_orig3;

            // Step 2: Calculate vectors with PBC minimum image convention
            auto disp1 = simulation_box::minimum_image(x1 - x2, box_x, box_y, box_z);
            auto disp2 = simulation_box::minimum_image(x3 - x2, box_x, box_y, box_z);

            float norm1_sq = mag2(disp1);
            float norm2_sq = mag2(disp2);
            if(norm1_sq == 0.f || norm2_sq == 0.f) continue;

            float norm1 = sqrtf(norm1_sq);
            float norm2 = sqrtf(norm2_sq);
            float dot_p = dot(disp1, disp2);
            float dp = dot_p / (norm1 * norm2); // This is cos(theta)
            dp = std::max(-1.0f, std::min(1.0f, dp)); // Clamp for safety

            // Step 3: Calculate potential and derivative based on V = 1/2 * k * (cos(theta) - cos(theta_0))^2
            // Convert equilibrium angle from degrees to radians, then to cosine
            float equil_angle_rad = p.equil_angle_deg * M_PI / 180.0f;
            float cos_theta0 = cosf(equil_angle_rad);
            float delta_cos = dp - cos_theta0; // dp is cos(theta), cos_theta0 is cos(theta_0)

            if (pot) *pot += 0.5f * p.spring_constant * delta_cos * delta_cos;

            // Step 4: Calculate forces using the chain rule: F = -dV/dx = -(dV/d(cos(theta))) * (d(cos(theta))/dx)
            float dV_dcos = p.spring_constant * delta_cos;
            
            // Use standard formulas to get forces on atoms from dV/d(cos(theta))
            float3 r1 = disp1;
            float3 r2 = disp2;
            float inv_norm1 = 1.0f / norm1;
            float inv_norm2 = 1.0f / norm2;

            float3 dcos_dr1 = (r2 * inv_norm1 - r1 * (dot_p * inv_norm1 / norm1_sq)) * inv_norm2;
            float3 dcos_dr2 = (r1 * inv_norm2 - r2 * (dot_p * inv_norm2 / norm2_sq)) * inv_norm1;
            
            // Gradient on atoms 1 and 3 (dE/dx)
            float3 grad_a = dV_dcos * dcos_dr1;
            float3 grad_c = dV_dcos * dcos_dr2;
            // Gradient on atom 2
            float3 grad_b = -grad_a - grad_c;

            update_vec(pos_sens, p.atom[0], grad_a);
            update_vec(pos_sens, p.atom[1], grad_b);
            update_vec(pos_sens, p.atom[2], grad_c);
        }
    }
    
};
static RegisterNodeType<AngleSpring, 1> angle_spring_node("angle_spring");

// Position restraint potential for restraining atoms to reference positions
// Applies harmonic penalty: V = 0.5 * k * (r - r_ref)^2
struct PositionRestraint : public PotentialNode
{
    struct Params {
        index_t atom_index;
        Vec<3,float> ref_pos;
        Vec<3,float> spring_const_xyz;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;

    PositionRestraint(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(1, grp, "restraint_indices")[0]), pos(pos_), params(n_elem)
    {
        check_size(grp, "restraint_indices", n_elem);
        check_size(grp, "ref_pos", n_elem, 3);
        bool has_spring_const_xyz = h5_exists(grp, "spring_const_xyz");
        if(has_spring_const_xyz) {
            check_size(grp, "spring_const_xyz", n_elem, 3);
        } else {
            check_size(grp, "spring_const", n_elem);
        }

        auto& p = params;
        traverse_dset<1,int>(grp, "restraint_indices", [&](size_t i, int x) {p[i].atom_index = x;});
        traverse_dset<2,float>(grp, "ref_pos", [&](size_t i, size_t j, float x) {
            if(j == 0) p[i].ref_pos.x() = x;
            else if(j == 1) p[i].ref_pos.y() = x;
            else if(j == 2) p[i].ref_pos.z() = x;
        });
        if(has_spring_const_xyz) {
            traverse_dset<2,float>(grp, "spring_const_xyz", [&](size_t i, size_t j, float x) {
                if(j == 0) p[i].spring_const_xyz.x() = x;
                else if(j == 1) p[i].spring_const_xyz.y() = x;
                else if(j == 2) p[i].spring_const_xyz.z() = x;
            });
        } else {
            traverse_dset<1,float>(grp, "spring_const", [&](size_t i, float x) {
                p[i].spring_const_xyz = make_vec3(x, x, x);
            });
        }

    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("position_restraint"));

        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;

        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];

            auto x = load_vec<3>(posc, p.atom_index);
            auto disp = x - p.ref_pos;

            // Harmonic potential with anisotropic force constants:
            // V = 0.5 * (kx*dx^2 + ky*dy^2 + kz*dz^2)
            if(pot) {
                *pot += 0.5f * (
                    p.spring_const_xyz.x() * disp.x() * disp.x() +
                    p.spring_const_xyz.y() * disp.y() * disp.y() +
                    p.spring_const_xyz.z() * disp.z() * disp.z()
                );
            }

            // Gradient (sens): dV/dx = kx*dx, dV/dy = ky*dy, dV/dz = kz*dz
            auto deriv = make_vec3(
                p.spring_const_xyz.x() * disp.x(),
                p.spring_const_xyz.y() * disp.y(),
                p.spring_const_xyz.z() * disp.z()
            );
            update_vec(pos_sens, p.atom_index, deriv);
        }
    }
};
static RegisterNodeType<PositionRestraint, 1> position_restraint_node("restraint_position");

// Explicit registrar to ensure node types are available at runtime
// Even if some linkers strip unused static objects, this guarantees registration
namespace {
void update_martini_node_boxes(DerivEngine& engine, float scale_xy, float scale_z) {
    if(!std::isfinite(scale_xy) || !std::isfinite(scale_z)) return;
    if(!(scale_xy > 0.f) || !(scale_z > 0.f)) return;

    for(auto& n : engine.nodes) {
        if(!n.computation) continue;
        if(is_prefix("martini_potential", n.name)) {
            if(auto* node = dynamic_cast<MartiniPotential*>(n.computation.get())) {
                node->box_x *= scale_xy;
                node->box_y *= scale_xy;
                node->box_z *= scale_z;
            }
            continue;
        }
        if(is_prefix("martini_sc_table_potential", n.name)) {
            if(auto* node = dynamic_cast<MartiniScTablePotential*>(n.computation.get())) {
                node->box_x *= scale_xy;
                node->box_y *= scale_xy;
                node->box_z *= scale_z;
            }
            continue;
        }
        if(is_prefix("martini_sc_table_1body", n.name)) {
            if(auto* node = dynamic_cast<MartiniScTableOneBody*>(n.computation.get())) {
                node->box_x *= scale_xy;
                node->box_y *= scale_xy;
                node->box_z *= scale_z;
            }
            continue;
        }
        if(is_prefix("dist_spring", n.name)) {
            if(auto* node = dynamic_cast<DistSpring*>(n.computation.get())) {
                node->box_x *= scale_xy;
                node->box_y *= scale_xy;
                node->box_z *= scale_z;
            }
            continue;
        }
        if(is_prefix("angle_spring", n.name)) {
            if(auto* node = dynamic_cast<AngleSpring*>(n.computation.get())) {
                node->box_x *= scale_xy;
                node->box_y *= scale_xy;
                node->box_z *= scale_z;
            }
            continue;
        }
        if(is_prefix("dihedral_spring", n.name)) {
            if(auto* node = dynamic_cast<DihedralSpring*>(n.computation.get())) {
                node->box_x *= scale_xy;
                node->box_y *= scale_xy;
                node->box_z *= scale_z;
            }
            continue;
        }
        if(auto* node = dynamic_cast<PotentialNode*>(n.computation.get())) {
            node->update_box_dimensions_anisotropic(scale_xy, scale_z);
        }
    }
}

struct MartiniNodeRegistrar {
    MartiniNodeRegistrar() {
        auto& m = node_creation_map();
        simulation_box::npt::register_node_box_updater(update_martini_node_boxes);
        if(m.find("martini_potential") == m.end()) {
            add_node_creation_function("martini_potential", [](hid_t grp, const ArgList& args) {
                check_arguments_length(args,1);
                return new MartiniPotential(grp, *args[0]);
            });
        }
        if(m.find("martini_sc_table_potential") == m.end()) {
            add_node_creation_function("martini_sc_table_potential", [](hid_t grp, const ArgList& args) {
                check_arguments_length(args,2);
                return new MartiniScTablePotential(grp, *args[0], *args[1]);
            });
        }
        if(m.find("martini_sc_table_1body") == m.end()) {
            add_node_creation_function("martini_sc_table_1body", [](hid_t grp, const ArgList& args) {
                check_arguments_length(args,2);
                return new MartiniScTableOneBody(grp, *args[0], *args[1]);
            });
        }
        if(m.find("dist_spring") == m.end()) {
            add_node_creation_function("dist_spring", [](hid_t grp, const ArgList& args) {
                check_arguments_length(args,1);
                return new DistSpring(grp, *args[0]);
            });
        }
        if(m.find("angle_spring") == m.end()) {
            add_node_creation_function("angle_spring", [](hid_t grp, const ArgList& args) {
                check_arguments_length(args,1);
                return new AngleSpring(grp, *args[0]);
            });
        }
        if(m.find("dihedral_spring") == m.end()) {
            add_node_creation_function("dihedral_spring", [](hid_t grp, const ArgList& args) {
                check_arguments_length(args,1);
                return new DihedralSpring(grp, *args[0]);
            });
        }
    }
};
static MartiniNodeRegistrar s_martini_node_registrar;
}

// Lightweight pre-run minimization using gradient descent with backtracking
void martini_run_minimization(DerivEngine& engine,
                                          int max_iterations,
                                          double energy_tolerance,
                                          double force_tolerance,
                                          double initial_step_size,
                                          int verbose)
{
    martini_fix_rigid::apply_fix_rigid_projection(engine, engine.pos->output);
    engine.compute(PotentialAndDerivMode);
    martini_fix_rigid::apply_fix_rigid_minimization(engine, engine.pos->output, engine.pos->sens);
    double prev_potential = engine.potential;
    if(verbose) printf("MIN: Initial potential %.6f\n", prev_potential);

    const int n_atom = engine.pos->n_elem;
    VecArray position = engine.pos->output;

    auto grad_norm = [&]() -> double {
        double sum = 0.0;
        for(int i=0;i<n_atom;++i){
            float3 g = load_vec<3>(engine.pos->sens, i);
            sum += double(g.x())*g.x() + double(g.y())*g.y() + double(g.z())*g.z();
        }
        return sqrt(sum);
    };

    std::vector<float> saved_pos(3*n_atom);
    std::vector<float> descent(3*n_atom);

    int iter = 0;
    double step = initial_step_size;
    while(iter < max_iterations){
        // Apply fix rigid constraints to forces before building descent
        martini_fix_rigid::apply_fix_rigid_minimization(engine, position, engine.pos->sens);
        
        // Save positions and build descent = -grad
        for(int i=0;i<n_atom;++i){
            float3 p = load_vec<3>(position, i);
            saved_pos[3*i+0] = p.x();
            saved_pos[3*i+1] = p.y();
            saved_pos[3*i+2] = p.z();
            float3 g = load_vec<3>(engine.pos->sens, i);
            descent[3*i+0] = -g.x();
            descent[3*i+1] = -g.y();
            descent[3*i+2] = -g.z();
        }

        double gnorm = grad_norm();
        if(gnorm < force_tolerance){
            if(verbose) printf("MIN: Converged by force tol |grad|=%.3g at iter %d\n", gnorm, iter);
            break;
        }

        // Backtracking line search
        double alpha = step;
        const double min_alpha = 1e-6;
        bool accepted = false;
        double best_pot = prev_potential;
        while(alpha >= min_alpha){
            for(int i=0;i<n_atom;++i){
                auto p = make_vec3(
                    saved_pos[3*i+0] + float(alpha)*descent[3*i+0],
                    saved_pos[3*i+1] + float(alpha)*descent[3*i+1],
                    saved_pos[3*i+2] + float(alpha)*descent[3*i+2]
                );
                store_vec<3>(position, i, p);
            }
            martini_fix_rigid::apply_fix_rigid_projection(engine, position);
            engine.compute(PotentialAndDerivMode);
            martini_fix_rigid::apply_fix_rigid_minimization(engine, position, engine.pos->sens);
            if(engine.potential <= best_pot){
                best_pot = engine.potential;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if(!accepted){
            // restore and stop
            for(int i=0;i<n_atom;++i){
                auto p = make_vec3(saved_pos[3*i+0], saved_pos[3*i+1], saved_pos[3*i+2]);
                store_vec<3>(position, i, p);
            }
            if(verbose) printf("MIN: Line search failed; stopping at iter %d\n", iter);
            break;
        }

        double dE = prev_potential - best_pot;
        if(verbose && !(iter%10)){
            printf("MIN: iter %4d step %.3g E %.6f dE %.3g |grad| %.3g\n", iter, alpha, best_pot, dE, gnorm);
        }
        if(dE >= 0.0 && fabs(dE) < energy_tolerance){
            if(verbose) printf("MIN: Converged by energy tol dE=%.3g at iter %d\n", dE, iter);
            break;
        }

        prev_potential = best_pot;
        step = std::min(2.0*alpha, step*1.5);
        ++iter;
    }
    martini_fix_rigid::apply_fix_rigid_projection(engine, position);
    engine.compute(PotentialAndDerivMode);
    martini_fix_rigid::apply_fix_rigid_minimization(engine, position, engine.pos->sens);
    prev_potential = engine.potential;
    if(verbose) printf("MIN: Final potential %.6f after %d iterations\n", prev_potential, iter);
}
