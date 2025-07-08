#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include <iostream>
#include <H5Apublic.h> // for H5Aexists

using namespace h5;
using namespace std;
//Bond, Angle and Dihedral the same format in MARTINI 10.1021/jp071097f
//Missing: Proper Dihedral from 10.1021/ct700324x (might not need if only exist in protein model)
//Coulomb interactions can be softened using Slater potential: V(r) = q1*q2/r * (1 - (1 + αr/2) * exp(-αr))
//Reference: LAMMPS pair_coul_slater_cut implementation

// Helper to check if an HDF5 attribute exists
inline bool attribute_exists(hid_t loc_id, const char* obj_name, const char* attr_name) {
    hid_t obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT);
    if (obj_id < 0) return false;
    htri_t exists = H5Aexists(obj_id, attr_name);
    H5Oclose(obj_id);
    return exists > 0;
}

// Helper for minimum image convention in a rectangular box
inline Vec<3,float> minimum_image_rect(const Vec<3,float>& dr, float box_x, float box_y, float box_z) {
    Vec<3,float> out = dr;
    if (out.x() >  0.5f * box_x) out.x() -= box_x;
    if (out.x() < -0.5f * box_x) out.x() += box_x;
    if (out.y() >  0.5f * box_y) out.y() -= box_y;
    if (out.y() < -0.5f * box_y) out.y() += box_y;
    if (out.z() >  0.5f * box_z) out.z() -= box_z;
    if (out.z() < -0.5f * box_z) out.z() += box_z;
    return out;
}

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

    DistSpring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(2, grp, "id")[0]), pos(pos_), params(n_elem)
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

        if(logging(LOG_DETAILED))
            default_logger->add_logger<float>("nonbonded_spring_energy", {1}, [&](float* buffer) {
                    float pot = 0.f;
                    VecArray x = pos.output;

                    for(int nt=0; nt<n_elem; ++nt) {
                    if(bonded_atoms[nt]) continue;  // don't count bonded spring energy

                    auto p = params[nt];
                    float dmag = mag(load_vec<3>(x, p.atom[0]) - load_vec<3>(x, p.atom[1]));
                    pot += 0.5f * p.spring_constant * sqr(dmag - p.equil_dist);
                    }
                    buffer[0] = pot;
                    });
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("dist_spring"));

        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;

        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];

            auto x1 = load_vec<3>(posc, p.atom[0]);
            auto x2 = load_vec<3>(posc, p.atom[1]);

            auto disp = x1 - x2;
            auto deriv = p.spring_constant * (1.f - p.equil_dist*inv_mag(disp)) * disp;
            if(pot) *pot += 0.5f * p.spring_constant * sqr(mag(disp) - p.equil_dist);

            update_vec(pos_sens, p.atom[0],  deriv);
            update_vec(pos_sens, p.atom[1], -deriv);
        }
    }
};
static RegisterNodeType<DistSpring,1> dist_spring_node("dist_spring");


struct AngleSpring : public PotentialNode
{
    struct Params {
        index_t atom[3];
        float equil_dp;
        float spring_constant;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;

    AngleSpring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(2, grp, "id")[0]), pos(pos_), params(n_elem)
    {
        int n_dep = 3;  // number of atoms that each term depends on 
        check_size(grp, "id",              n_elem, n_dep);
        check_size(grp, "equil_dist",        n_elem);
        check_size(grp, "spring_const", n_elem);

        auto& p = params;
        traverse_dset<2,int>  (grp, "id",           [&](size_t i, size_t j, int   x) { p[i].atom[j] = x;});
        traverse_dset<1,float>(grp, "equil_dist",   [&](size_t i,           float x) { p[i].equil_dp = x;});
        traverse_dset<1,float>(grp, "spring_const", [&](size_t i,           float x) { p[i].spring_constant = x;});
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("angle_spring"));

        float* posc = pos.output.x.get();
        float* pos_sens = pos.sens.x.get();
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;

        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            auto atom1 = Float4(posc + 4*p.atom[0]);
            auto atom2 = Float4(posc + 4*p.atom[1]);
            auto atom3 = Float4(posc + 4*p.atom[2]);

            auto x1 = atom1 - atom3; auto inv_d1 = inv_mag(x1); auto x1h = x1*inv_d1;
            auto x2 = atom2 - atom3; auto inv_d2 = inv_mag(x2); auto x2h = x2*inv_d2;

            auto dp = dot(x1h, x2h);
            auto force_prefactor = Float4(p.spring_constant) * (dp - Float4(p.equil_dp));

            auto d1 = force_prefactor * (x2h - x1h*dp) * inv_d1;
            auto d2 = force_prefactor * (x1h - x2h*dp) * inv_d2;
            auto d3 = -d1-d2;

            d1.update(pos_sens + 4*p.atom[0]);
            d2.update(pos_sens + 4*p.atom[1]);
            d3.update(pos_sens + 4*p.atom[2]);

            if(pot) *pot += 0.5f * p.spring_constant * sqr(dp.x()-p.equil_dp);
        }
    }
};
static RegisterNodeType<AngleSpring,1> angle_spring_node("angle_spring");


struct DihedralSpring : public PotentialNode
{
    struct Params {
        index_t atom[4];
        float equil_dihedral;
        float spring_constant;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;

    DihedralSpring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(2, grp, "id")[0]), pos(pos_), params(n_elem)
    {
        int n_dep = 4;  // number of atoms that each term depends on 
        check_size(grp, "id",           n_elem, n_dep);
        check_size(grp, "equil_dist",   n_elem);
        check_size(grp, "spring_const", n_elem);

        auto& p = params;
        traverse_dset<2,int>  (grp, "id",           [&](size_t i, size_t j, int   x) {p[i].atom[j]  =x;});
        traverse_dset<1,float>(grp, "equil_dist",   [&](size_t i,           float x) {p[i].equil_dihedral =x;});
        traverse_dset<1,float>(grp, "spring_const", [&](size_t i,           float x) {p[i].spring_constant=x;});
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("dihedral_spring"));

        float* posc = pos.output.x.get();
        float* pos_sens = pos.sens.x.get();
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;

        for(int nt=0; nt<n_elem; ++nt) {
            const auto& p = params[nt];
            Float4 x[4];
            for(int na: range(4)) x[na] = Float4(posc + 4*params[nt].atom[na]);

            Float4 d[4];
            float dihedral = dihedral_germ(x[0],x[1],x[2],x[3], d[0],d[1],d[2],d[3]).x();

            // determine minimum periodic image (can be off by at most 2pi)
            float displacement = dihedral - p.equil_dihedral;
            displacement = (displacement> M_PI_F) ? displacement-2.f*M_PI_F : displacement;
            displacement = (displacement<-M_PI_F) ? displacement+2.f*M_PI_F : displacement;

            auto s = Float4(p.spring_constant * displacement);
            for(int na: range(4)) d[na].scale_update(s, pos_sens + 4*params[nt].atom[na]);

            if(pot) *pot += 0.5f * p.spring_constant * sqr(displacement);
        }
    }
};
static RegisterNodeType<DihedralSpring,1> dihedral_spring_node("dihedral_spring");

struct PeriodicBoundaryPotential : public PotentialNode
{
    int n_atom;
    CoordNode& pos;
    
    // Box parameters
    float box_x, box_y, box_z;
    
    PeriodicBoundaryPotential(hid_t grp, CoordNode& pos_):
        PotentialNode(), n_atom(pos_.n_elem), pos(pos_)
    {
        // Read box dimensions from attributes
        // Support both individual dimensions and legacy wall_box_size
        if(attribute_exists(grp, ".", "x_len") && attribute_exists(grp, ".", "y_len") && attribute_exists(grp, ".", "z_len")) {
            // New format: separate x, y, z dimensions
            box_x = read_attribute<float>(grp, ".", "x_len");
            box_y = read_attribute<float>(grp, ".", "y_len");
            box_z = read_attribute<float>(grp, ".", "z_len");
            std::cout << "PERIODIC: Using new box format - X=" << box_x << ", Y=" << box_y << ", Z=" << box_z << std::endl;
        } else {
            // Legacy format: wall boundaries
            float wall_xlo = read_attribute<float>(grp, ".", "wall_xlo");
            float wall_xhi = read_attribute<float>(grp, ".", "wall_xhi");
            float wall_ylo = read_attribute<float>(grp, ".", "wall_ylo");
            float wall_yhi = read_attribute<float>(grp, ".", "wall_yhi");
            float wall_zlo = read_attribute<float>(grp, ".", "wall_zlo");
            float wall_zhi = read_attribute<float>(grp, ".", "wall_zhi");
            
            // Calculate box dimensions
            box_x = wall_xhi - wall_xlo;
            box_y = wall_yhi - wall_ylo;
            box_z = wall_zhi - wall_zlo;
            std::cout << "PERIODIC: Using legacy box format - X=" << box_x << ", Y=" << box_y << ", Z=" << box_z << std::endl;
        }
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("periodic_boundary_potential"));
        
        VecArray pos1 = pos.output;
        
        // Apply periodic boundary conditions
        for(int na=0; na<n_atom; ++na) {
            auto p = load_vec<3>(pos1, na);
            bool wrapped = false;
            
            // X dimension
            if (p.x() < 0) {
                p.x() += box_x;
                wrapped = true;
            }
            if (p.x() >= box_x) {
                p.x() -= box_x;
                wrapped = true;
            }
            
            // Y dimension
            if (p.y() < 0) {
                p.y() += box_y;
                wrapped = true;
            }
            if (p.y() >= box_y) {
                p.y() -= box_y;
                wrapped = true;
            }
            
            // Z dimension
            if (p.z() < 0) {
                p.z() += box_z;
                wrapped = true;
            }
            if (p.z() >= box_z) {
                p.z() -= box_z;
                wrapped = true;
            }
            
            // Store wrapped position
            if (wrapped) {
                store_vec<3>(pos1, na, p);
            }
        }
        
        // No potential energy contribution from periodic boundaries
        potential = 0.f;
    }
};
static RegisterNodeType<PeriodicBoundaryPotential,1> periodic_boundary_potential_node("periodic_boundary_potential");

struct MartiniPotential : public PotentialNode
{
    int n_atom;
    CoordNode& pos;
    
    vector<array<float,4>> coeff;
    vector<pair<int,int>> pairs;
    
    float epsilon, sigma, lj_cutoff, coul_cutoff, dielectric;
    bool force_cap;
    bool coulomb_soften;
    float slater_alpha;

    // Box dimensions for minimum image
    float box_x, box_y, box_z;
    
    // Debug variables
    bool debug_mode;
    int debug_step_count;
    int max_debug_interactions;
    
    MartiniPotential(hid_t grp, CoordNode& pos_):
        PotentialNode(), n_atom(pos_.n_elem), pos(pos_), debug_step_count(0), max_debug_interactions(10)
    {
        check_size(grp, "atom_indices", n_atom);
        check_size(grp, "charges", n_atom);
        
        epsilon     = read_attribute<float>(grp, ".", "epsilon");
        sigma       = read_attribute<float>(grp, ".", "sigma");  
        lj_cutoff   = read_attribute<float>(grp, ".", "lj_cutoff");
        coul_cutoff = read_attribute<float>(grp, ".", "coul_cutoff");
        dielectric  = read_attribute<float>(grp, ".", "dielectric");
        force_cap = true;
        if(attribute_exists(grp, ".", "force_cap")) {
            force_cap = read_attribute<int>(grp, ".", "force_cap") != 0;
        }
        
        // Coulomb softening parameters
        coulomb_soften = false;
        if(attribute_exists(grp, ".", "coulomb_soften")) {
            coulomb_soften = read_attribute<int>(grp, ".", "coulomb_soften") != 0;
        }
        
        slater_alpha = 0.0f;
        if(coulomb_soften) {
            if(attribute_exists(grp, ".", "slater_alpha")) {
                slater_alpha = read_attribute<float>(grp, ".", "slater_alpha");
            } else {
                // Default value if not specified
                slater_alpha = 1.0f;
            }
            std::cout << "[DEBUG] SOFTENED COULOMB POTENTIAL IS ENABLED! slater_alpha=" << slater_alpha << std::endl;
        }
        
        // Debug mode - enable for first few steps
        debug_mode = false;
        if(attribute_exists(grp, ".", "debug_mode")) {
            debug_mode = read_attribute<int>(grp, ".", "debug_mode") != 0;
        }
        

        
        // Read box dimensions (same as PeriodicBoundaryPotential)
        // Support both individual dimensions and legacy wall_box_size
        if(attribute_exists(grp, ".", "x_len") && attribute_exists(grp, ".", "y_len") && attribute_exists(grp, ".", "z_len")) {
            // New format: separate x, y, z dimensions
            box_x = read_attribute<float>(grp, ".", "x_len");
            box_y = read_attribute<float>(grp, ".", "y_len");
            box_z = read_attribute<float>(grp, ".", "z_len");
            std::cout << "MARTINI: Using new box format - X=" << box_x << ", Y=" << box_y << ", Z=" << box_z << std::endl;
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
            std::cout << "MARTINI: Using legacy box format - X=" << box_x << ", Y=" << box_y << ", Z=" << box_z << std::endl;
        }
        
        auto n_pair = get_dset_size(2, grp, "pairs")[0];
        check_size(grp, "coefficients", n_pair, 4);
        
        pairs.resize(n_pair);
        traverse_dset<2,int>(grp, "pairs", [&](size_t np, size_t d, int x) {
            if(d == 0) pairs[np].first = x;
            else pairs[np].second = x;
        });
        
        coeff.resize(n_pair);
        traverse_dset<2,float>(grp, "coefficients", [&](size_t np, size_t d, float x) {
            coeff[np][d] = x;
        });
        
        // Print debug info if enabled
        if(debug_mode) {
            std::cout << "MARTINI DEBUG: Initialized with " << n_atom << " atoms and " << n_pair << " pairs" << std::endl;
            std::cout << "MARTINI DEBUG: Parameters - epsilon=" << epsilon << ", sigma=" << sigma 
                      << ", lj_cutoff=" << lj_cutoff << ", coul_cutoff=" << coul_cutoff 
                      << ", dielectric=" << dielectric << std::endl;
        }
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("martini_potential"));
        
        VecArray pos1      = pos.output;
        VecArray pos1_sens = pos.sens;
        
        fill(pos1_sens, 3, n_atom, 0.f);
        
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;
        
        // Debug: Print particle positions for first few steps
        if(debug_mode && debug_step_count < 3) {
            std::cout << "\n=== STEP " << debug_step_count << " ===" << std::endl;
            std::cout << "Coordinates:" << std::endl;
            for(int i = 0; i < n_atom; ++i) {
                auto p = load_vec<3>(pos1, i);
                std::cout << "  Particle " << i << ": (" << p.x() << ", " << p.y() << ", " << p.z() << ")" << std::endl;
            }
        }
        
        // Compute particle-particle interactions
        int debug_interaction_count = 0;
        for(size_t np=0; np<pairs.size(); ++np) {
            int i = pairs[np].first;
            int j = pairs[np].second;
            
            auto eps   = coeff[np][0];
            auto sig   = coeff[np][1];
            auto qi    = coeff[np][2];
            auto qj    = coeff[np][3];
            
            if(eps==0.f && sig==0.f && qi==0.f && qj==0.f) continue;
            
            auto p1 = load_vec<3>(pos1, i);
            auto p2 = load_vec<3>(pos1, j);
            // Apply minimum image convention for rectangular box
            auto dr = minimum_image_rect(p1 - p2, box_x, box_y, box_z);
            auto dist2 = mag2(dr);
            auto dist = sqrtf(dist2);
            
            // Apply distance cutoff for computational efficiency
            if(dist > max(lj_cutoff, coul_cutoff)) continue;
            
            // Debug: Print detailed interaction info for first few interactions
            if(debug_mode && debug_step_count < 3 && debug_interaction_count < max_debug_interactions) {
                std::cout << "\nInteraction " << debug_interaction_count << " (particles " << i << " and " << j << "):" << std::endl;
                std::cout << "  Distance: " << dist << " Angstroms" << std::endl;
                std::cout << "  Parameters: eps=" << eps << ", sig=" << sig << ", qi=" << qi << ", qj=" << qj << std::endl;
            }
            
            Vec<3> force = make_zero<3>();
            float lj_pot_contrib = 0.f;
            float coul_pot_contrib = 0.f;
            
            // Lennard-Jones potential
            if(eps != 0.f && sig != 0.f && dist < lj_cutoff) {
                // Use bare LJ potential without any protection
                auto sig_r = sig / dist;
                auto sig_r6 = sig_r * sig_r * sig_r * sig_r * sig_r * sig_r;
                auto sig_r12 = sig_r6 * sig_r6;
                
                // Check for overflow in LJ terms
                if(std::isfinite(sig_r6) && std::isfinite(sig_r12)) {
                    auto lj_pot = 4.f * eps * (sig_r12 - sig_r6);
                    // Corrected LJ force derivative: 24*epsilon/sigma*(2*(sigma/r)^13-(sigma/r)^7)
                    auto deriv = 24.f * eps / sig * (2.f * sig_r12 * sig_r - sig_r6 * sig_r);
                    
                    // Only apply if potential and derivative are finite
                    if(std::isfinite(lj_pot) && std::isfinite(deriv)) {
                        if(pot) *pot += lj_pot;
                        lj_pot_contrib = lj_pot;
                        force += (deriv/dist) * (-dr);  // FIXED: Use unit vector, not displacement vector
                        
                        // Debug: Print LJ details
                        if(debug_mode && debug_step_count < 3 && debug_interaction_count < max_debug_interactions) {
                            std::cout << "  LJ Energy: " << lj_pot << " UPSIDE units" << std::endl;
                            std::cout << "  LJ Force: (" << (deriv * (-dr)).x() << ", " << (deriv * (-dr)).y() << ", " << (deriv * (-dr)).z() << ") UPSIDE units" << std::endl;
                        }
                    }
                }
            }
            
            // Coulomb potential
            if(qi != 0.f && qj != 0.f && dist < coul_cutoff) {
                const float ke = 332.0636f;  // Coulomb constant in UPSIDE units
                float coul_pot, deriv;
                
                if(coulomb_soften && slater_alpha > 0.0f) {
                    // Slater softened Coulomb potential
                    // V(r) = q1*q2/r * (1 - (1 + αr/2) * exp(-αr))
                    float alpha_r = slater_alpha * dist;
                    float exp_alpha_r = expf(-alpha_r);
                    float soften_factor = 1.0f - (1.0f + alpha_r * 0.5f) * exp_alpha_r;
                    
                    coul_pot = ke * qi * qj * soften_factor / (dielectric * dist);
                    
                    // Derivative: dV/dr = q1*q2 * [-(1 - (1 + αr/2) * exp(-αr))/r^2 + α/2 * exp(-αr)/r]
                    float d_soften_dr = slater_alpha * 0.5f * exp_alpha_r / dist;
                    deriv = -ke * qi * qj * (soften_factor / dist2 - d_soften_dr) / dielectric;
                } else {
                    // Standard Coulomb potential
                    coul_pot = ke * qi * qj / (dielectric * dist);
                    deriv = -coul_pot / dist2;
                }
                
                // Only apply if potential and derivative are finite
                if(std::isfinite(coul_pot) && std::isfinite(deriv)) {
                    if(pot) *pot += coul_pot;
                    coul_pot_contrib = coul_pot;
                    force += (deriv/dist) * (-dr);  // FIXED: Use unit vector, not displacement vector
                    
                    // Debug: Print Coulomb details
                    if(debug_mode && debug_step_count < 3 && debug_interaction_count < max_debug_interactions) {
                        std::cout << "  Coulomb Energy: " << coul_pot << " UPSIDE units" << std::endl;
                        std::cout << "  Coulomb Force: (" << (deriv * dr).x() << ", " << (deriv * dr).y() << ", " << (deriv * dr).z() << ") UPSIDE units" << std::endl;
                    }
                }
            }
            
            // Debug: Print total interaction summary
            if(debug_mode && debug_step_count < 3 && debug_interaction_count < max_debug_interactions) {
                std::cout << "  Total Energy: " << (lj_pot_contrib + coul_pot_contrib) << " UPSIDE units" << std::endl;
                std::cout << "  Total Force: (" << force.x() << ", " << force.y() << ", " << force.z() << ") UPSIDE units" << std::endl;
            }
            
            update_vec<3>(pos1_sens, i,  force);
            update_vec<3>(pos1_sens, j, -force);
            
            debug_interaction_count++;
        }
        
        // Debug: Print total potential for this step
        if(debug_mode && debug_step_count < 3) {
            std::cout << "Total Energy: " << (pot ? *pot : 0.f) << " UPSIDE units" << std::endl;
            std::cout << "=== END STEP " << debug_step_count << " ===\n" << std::endl;
        }
        
        debug_step_count++;
    }
};
static RegisterNodeType<MartiniPotential,1> martini_potential_node("martini_potential");

// Conjugate Gradient Minimizer based on LAMMPS implementation
struct ConjugateGradientMinimizer : public PotentialNode
{
    CoordNode& pos;
    int n_atom;
    int max_iterations;
    float energy_tolerance;
    float force_tolerance;
    float step_size;
    bool verbose;
    
    // Minimization state variables
    vector<Vec<3,float>> forces;
    vector<Vec<3,float>> search_direction;
    vector<Vec<3,float>> old_forces;
    float old_energy;
    int iteration_count;
    bool converged;
    
    ConjugateGradientMinimizer(hid_t grp, CoordNode& pos_):
        PotentialNode(), pos(pos_), n_atom(pos_.n_elem), iteration_count(0), converged(false)
    {
        // Read minimization parameters
        max_iterations = 1000;
        if(attribute_exists(grp, ".", "max_iterations")) {
            max_iterations = read_attribute<int>(grp, ".", "max_iterations");
        }
        
        energy_tolerance = 1e-6f;
        if(attribute_exists(grp, ".", "energy_tolerance")) {
            energy_tolerance = read_attribute<float>(grp, ".", "energy_tolerance");
        }
        
        force_tolerance = 1e-6f;
        if(attribute_exists(grp, ".", "force_tolerance")) {
            force_tolerance = read_attribute<float>(grp, ".", "force_tolerance");
        }
        
        step_size = 0.1f;
        if(attribute_exists(grp, ".", "step_size")) {
            step_size = read_attribute<float>(grp, ".", "step_size");
        }
        
        verbose = true;
        if(attribute_exists(grp, ".", "verbose")) {
            verbose = read_attribute<int>(grp, ".", "verbose") != 0;
        }
        
        // Initialize vectors
        forces.resize(n_atom);
        search_direction.resize(n_atom);
        old_forces.resize(n_atom);
        
        std::cout << "[MINIMIZER] Conjugate Gradient Minimizer initialized:" << std::endl;
        std::cout << "  Max iterations: " << max_iterations << std::endl;
        std::cout << "  Energy tolerance: " << energy_tolerance << std::endl;
        std::cout << "  Force tolerance: " << force_tolerance << std::endl;
        std::cout << "  Step size: " << step_size << std::endl;
    }
    
    // Compute force magnitude for convergence check
    float compute_force_magnitude() {
        float total_force = 0.0f;
        for(int i = 0; i < n_atom; ++i) {
            total_force += mag2(forces[i]);
        }
        return sqrtf(total_force);
    }
    
    // Line search to find optimal step size (simplified version)
    float line_search() {
        float alpha = step_size;
        float energy_0 = potential;
        
        // Store current positions
        vector<Vec<3,float>> original_pos(n_atom);
        VecArray pos_array = pos.output;
        for(int i = 0; i < n_atom; ++i) {
            original_pos[i] = load_vec<3>(pos_array, i);
        }
        
        // Try different step sizes
        for(int attempt = 0; attempt < 10; ++attempt) {
            // Move along search direction
            for(int i = 0; i < n_atom; ++i) {
                auto new_pos = original_pos[i] + alpha * search_direction[i];
                store_vec<3>(pos_array, i, new_pos);
            }
            
            // Recompute energy (this would need to trigger a full energy calculation)
            // For now, we'll use a simplified approach
            float new_energy = energy_0; // Placeholder - would need actual recomputation
            
            if(new_energy < energy_0) {
                return alpha;
            }
            
            alpha *= 0.5f; // Reduce step size
        }
        
        // Restore original positions
        for(int i = 0; i < n_atom; ++i) {
            store_vec<3>(pos_array, i, original_pos[i]);
        }
        
        return step_size * 0.1f; // Return small step size if no improvement found
    }
    
    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("conjugate_gradient_minimizer"));
        
        if(converged || iteration_count >= max_iterations) {
            if(verbose) {
                std::cout << "[MINIMIZER] Minimization " << (converged ? "converged" : "reached max iterations") << std::endl;
            }
            return;
        }
        
        VecArray pos_array = pos.output;
        VecArray pos_sens = pos.sens;
        
        // Extract forces from sensitivity (gradient)
        for(int i = 0; i < n_atom; ++i) {
            forces[i] = -load_vec<3>(pos_sens, i); // Force = -gradient
        }
        
        float current_energy = potential;
        float force_magnitude = compute_force_magnitude();
        
        if(verbose && iteration_count % 10 == 0) {
            std::cout << "[MINIMIZER] Iteration " << iteration_count 
                      << " Energy: " << current_energy 
                      << " Force: " << force_magnitude << std::endl;
        }
        
        // Check convergence
        if(iteration_count > 0) {
            float energy_change = fabsf(current_energy - old_energy);
            if(energy_change < energy_tolerance && force_magnitude < force_tolerance) {
                converged = true;
                if(verbose) {
                    std::cout << "[MINIMIZER] CONVERGED! Energy change: " << energy_change 
                              << " Force magnitude: " << force_magnitude << std::endl;
                }
                return;
            }
        }
        
        // Compute search direction using conjugate gradient method
        if(iteration_count == 0) {
            // First iteration: use steepest descent
            for(int i = 0; i < n_atom; ++i) {
                search_direction[i] = forces[i];
            }
        } else {
            // Conjugate gradient: Polak-Ribiere formula
            float numerator = 0.0f;
            float denominator = 0.0f;
            
            for(int i = 0; i < n_atom; ++i) {
                Vec<3,float> force_diff = forces[i] - old_forces[i];
                numerator += dot(forces[i], force_diff);
                denominator += dot(old_forces[i], old_forces[i]);
            }
            
            float beta = (denominator > 0.0f) ? numerator / denominator : 0.0f;
            
            // Ensure beta is positive (Fletcher-Reeves fallback)
            if(beta < 0.0f) {
                beta = 0.0f;
            }
            
            for(int i = 0; i < n_atom; ++i) {
                search_direction[i] = forces[i] + beta * search_direction[i];
            }
        }
        
        // Normalize search direction
        float search_magnitude = 0.0f;
        for(int i = 0; i < n_atom; ++i) {
            search_magnitude += mag2(search_direction[i]);
        }
        search_magnitude = sqrtf(search_magnitude);
        
        if(search_magnitude > 0.0f) {
            for(int i = 0; i < n_atom; ++i) {
                search_direction[i] = search_direction[i] / search_magnitude;
            }
        }
        
        // Perform line search to find optimal step size
        float alpha = line_search();
        
        // Update positions
        for(int i = 0; i < n_atom; ++i) {
            auto current_pos = load_vec<3>(pos_array, i);
            auto new_pos = current_pos + alpha * search_direction[i];
            store_vec<3>(pos_array, i, new_pos);
        }
        
        // Store current state for next iteration
        old_energy = current_energy;
        for(int i = 0; i < n_atom; ++i) {
            old_forces[i] = forces[i];
        }
        
        iteration_count++;
        
        // Set potential to current energy for logging
        potential = current_energy;
    }
};

static RegisterNodeType<ConjugateGradientMinimizer,1> cg_minimizer_node("conjugate_gradient_minimizer");


