#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include <iostream>

using namespace h5;
using namespace std;
//Bond, Angle and Dihedral the same format in MARTINI 10.1021/jp071097f
//Missing: Proper Dihedral from 10.1021/ct700324x (might not need if only exist in protein model)
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
    
    MartiniPotential(hid_t grp, CoordNode& pos_):
        PotentialNode(), n_atom(pos_.n_elem), pos(pos_)
    {
        check_size(grp, "atom_indices", n_atom);
        check_size(grp, "charges", n_atom);
        
        epsilon     = read_attribute<float>(grp, ".", "epsilon");
        sigma       = read_attribute<float>(grp, ".", "sigma");  
        lj_cutoff   = read_attribute<float>(grp, ".", "lj_cutoff");
        coul_cutoff = read_attribute<float>(grp, ".", "coul_cutoff");
        dielectric  = read_attribute<float>(grp, ".", "dielectric");
        
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
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("martini_potential"));
        
        VecArray pos1      = pos.output;
        VecArray pos1_sens = pos.sens;
        
        fill(pos1_sens, 3, n_atom, 0.f);
        
        float pot = 0.f;
        
        // Compute particle-particle interactions
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
            auto dr = p1 - p2;
            auto dist2 = mag2(dr);
            auto dist = sqrtf(dist2);
            
            if(dist > max(lj_cutoff, coul_cutoff)) continue;
            
            Vec<3> force = make_zero<3>();
            
            // Lennard-Jones potential
            if(eps != 0.f && sig != 0.f && dist < lj_cutoff) {
                // Enhanced overflow protection: prevent core overlap
                // Set minimum distance to 1.12*sigma (LJ minimum) to prevent unphysical overlap
                float min_dist = 1.12f * sig;  // 2^(1/6) ≈ 1.12 is the LJ minimum
                float effective_dist = max(dist, min_dist);
                
                auto sig_r = sig / effective_dist;
                auto sig_r6 = sig_r * sig_r * sig_r * sig_r * sig_r * sig_r;
                auto sig_r12 = sig_r6 * sig_r6;
                
                // Check for overflow in LJ terms
                if(std::isfinite(sig_r6) && std::isfinite(sig_r12)) {
                    auto lj_pot = 4.f * eps * (sig_r12 - sig_r6);
                    auto deriv = 24.f * eps / (effective_dist * effective_dist) * (2.f * sig_r12 - sig_r6);
                    
                    // Cap maximum force magnitude to prevent integration instability
                    float max_force_magnitude = 1000.0f * eps;  // Reasonable force limit
                    float force_mag = mag(deriv * dr);
                    if(force_mag > max_force_magnitude) {
                        deriv = deriv * (max_force_magnitude / force_mag);
                    }
                    
                    // Only apply if potential and derivative are finite
                    if(std::isfinite(lj_pot) && std::isfinite(deriv)) {
                        pot += lj_pot;
                        force += deriv * dr;
                    }
                }
            }
            
            // Coulomb potential
            if(qi != 0.f && qj != 0.f && dist < coul_cutoff) {
                const float ke = 332.0636f;  // Coulomb constant in UPSIDE units
                auto coul_pot = ke * qi * qj / (dielectric * dist);
                auto deriv = -coul_pot / dist2;
                
                // Only apply if potential and derivative are finite
                if(std::isfinite(coul_pot) && std::isfinite(deriv)) {
                    pot += coul_pot;
                    force += deriv * dr;
                }
            }
            
            update_vec<3>(pos1_sens, i,  force);
            update_vec<3>(pos1_sens, j, -force);
        }
        
        potential = pot;
    }
};
static RegisterNodeType<MartiniPotential,1> martini_potential_node("martini_potential");


