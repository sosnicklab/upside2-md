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

struct MartiniPotential : public PotentialNode {
    CoordNode& pos;
    int n_elem;
    vector<int> atom_indices;
    
    // LJ parameters
    float epsilon;
    float sigma;
    float lj_cutoff;
    
    // Coulombic parameters
    vector<float> charges;
    float coul_cutoff;
    float dielectric;
    
    // Extra particles
    vector<int> extra_indices;
    
    MartiniPotential(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        pos(pos_),
        n_elem(get_dset_size(1, grp, "atom_indices")[0]),
        atom_indices(n_elem),
        epsilon(read_attribute<float>(grp, ".", "epsilon")),
        sigma(read_attribute<float>(grp, ".", "sigma")),
        lj_cutoff(read_attribute<float>(grp, ".", "lj_cutoff")),
        charges(n_elem),
        coul_cutoff(read_attribute<float>(grp, ".", "coul_cutoff")),
        dielectric(read_attribute<float>(grp, ".", "dielectric")) {
        
        // Read atom indices
        traverse_dset<1,int>(grp, "atom_indices", 
            [&](size_t i, int idx) {atom_indices[i] = idx;});
            
        // Read charges
        traverse_dset<1,float>(grp, "charges",
            [&](size_t i, float q) {charges[i] = q;});
            
        // Read extra particle indices if they exist
        if(H5Lexists(grp, "extra_indices", H5P_DEFAULT)) {
            int n_extra = get_dset_size(1, grp, "extra_indices")[0];
            extra_indices.resize(n_extra);
            traverse_dset<1,int>(grp, "extra_indices",
                [&](size_t i, int idx) {extra_indices[i] = idx;});
        }
    }
    
    virtual void compute_value(ComputeMode mode) {
        VecArray pos_array = pos.output;
        VecArray pos_sens = pos.sens;
        
        potential = 0.f;
        
        // Combine regular and extra particles
        vector<int> all_indices = atom_indices;
        all_indices.insert(all_indices.end(), extra_indices.begin(), extra_indices.end());
        
        for(int i=0; i<all_indices.size(); ++i) {
            for(int j=i+1; j<all_indices.size(); ++j) {
                int idx1 = all_indices[i];
                int idx2 = all_indices[j];
                
                // Calculate distance
                float dx = pos_array(0,idx2) - pos_array(0,idx1);
                float dy = pos_array(1,idx2) - pos_array(1,idx1);
                float dz = pos_array(2,idx2) - pos_array(2,idx1);
                float r2 = dx*dx + dy*dy + dz*dz;
                float r = sqrt(r2);
                
                // LJ potential
                if(r2 <= lj_cutoff*lj_cutoff) {
                    float r6 = pow(sigma/r, 6);
                    float r12 = r6*r6;
                    float lj_energy = 4*epsilon*(r12 - r6);
                    potential += lj_energy;
                    
                    float lj_force = 24*epsilon*(2*r12 - r6)/(r*r);
                    float fx = lj_force * dx/r;
                    float fy = lj_force * dy/r;
                    float fz = lj_force * dz/r;
                    
                    pos_sens(0,idx1) -= fx;
                    pos_sens(1,idx1) -= fy;
                    pos_sens(2,idx1) -= fz;
                    pos_sens(0,idx2) += fx;
                    pos_sens(1,idx2) += fy;
                    pos_sens(2,idx2) += fz;
                }
                
                // Coulombic potential
                if(r2 <= coul_cutoff*coul_cutoff) {
                    float coul_energy = charges[i]*charges[j]/(4*M_PI*dielectric*r);
                    potential += coul_energy;
                    
                    float coul_force = charges[i]*charges[j]/(4*M_PI*dielectric*r2);
                    float fx = coul_force * dx/r;
                    float fy = coul_force * dy/r;
                    float fz = coul_force * dz/r;
                    
                    pos_sens(0,idx1) -= fx;
                    pos_sens(1,idx1) -= fy;
                    pos_sens(2,idx1) -= fz;
                    pos_sens(0,idx2) += fx;
                    pos_sens(1,idx2) += fy;
                    pos_sens(2,idx2) += fz;
                }
            }
        }
    }
};
static RegisterNodeType<MartiniPotential,1> martini_potential_node("martini_potential");


