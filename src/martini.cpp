#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include "spline.h"
#include <iostream>
#include <H5Apublic.h> // for H5Aexists
#include <fstream> // For file writing
#include <cmath> // For pow, cosf, sinf, acosf
#include <set> // For std::set

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
    
    // Box dimensions for minimum image
    float box_x, box_y, box_z;
    float mass_scale;  // Mass scaling factor for forces (1/mass)

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
        
        // Mass scaling for MARTINI force field (mass = 72.0)
        mass_scale = 1.0f / 72.0f;  // Default MARTINI mass
        if(attribute_exists(grp, ".", "mass_scale")) {
            mass_scale = read_attribute<float>(grp, ".", "mass_scale");
        }
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("dihedral_spring"));

        float* posc = pos.output.x.get();
        float* pos_sens = pos.sens.x.get();
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;

        for(int nt=0; nt<n_elem; ++nt) {
            const auto& p = params[nt];
            Float4 x_orig[4];
            for(int na: range(4)) x_orig[na] = Float4(posc + 4*params[nt].atom[na]);

            // Apply minimum image convention for periodic boundaries
            // For dihedral, we need to ensure all bonds use minimum image
            Float4 x[4];
            x[1] = x_orig[1]; // Use atom 1 as reference
            
            // Apply minimum image for bonds 0-1, 2-1, 3-2
            auto disp01 = minimum_image_rect(make_vec3(x_orig[0].x(), x_orig[0].y(), x_orig[0].z()) - make_vec3(x[1].x(), x[1].y(), x[1].z()), box_x, box_y, box_z);
            float temp01[4] = {x[1].x() + disp01.x(), x[1].y() + disp01.y(), x[1].z() + disp01.z(), 0.0f};
            x[0] = Float4(temp01);
            
            auto disp21 = minimum_image_rect(make_vec3(x_orig[2].x(), x_orig[2].y(), x_orig[2].z()) - make_vec3(x[1].x(), x[1].y(), x[1].z()), box_x, box_y, box_z);
            float temp21[4] = {x[1].x() + disp21.x(), x[1].y() + disp21.y(), x[1].z() + disp21.z(), 0.0f};
            x[2] = Float4(temp21);
            
            auto disp32 = minimum_image_rect(make_vec3(x_orig[3].x(), x_orig[3].y(), x_orig[3].z()) - make_vec3(x[2].x(), x[2].y(), x[2].z()), box_x, box_y, box_z);
            float temp32[4] = {x[2].x() + disp32.x(), x[2].y() + disp32.y(), x[2].z() + disp32.z(), 0.0f};
            x[3] = Float4(temp32);

            Float4 d[4];
            float dihedral = dihedral_germ(x[0],x[1],x[2],x[3], d[0],d[1],d[2],d[3]).x();

            // determine minimum periodic image (can be off by at most 2pi)
            float displacement = dihedral - p.equil_dihedral;
            displacement = (displacement> M_PI_F) ? displacement-2.f*M_PI_F : displacement;
            displacement = (displacement<-M_PI_F) ? displacement+2.f*M_PI_F : displacement;

            auto s = Float4(p.spring_constant * displacement * mass_scale);
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
    // Static box dimensions for PBC wrapping
    static float static_box_x, static_box_y, static_box_z;
    
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
        float half_x = 0.5f * box_x;
        float half_y = 0.5f * box_y;
        float half_z = 0.5f * box_z;
        // Apply periodic boundary conditions (centered box convention)
        for(int na=0; na<n_atom; ++na) {
            auto p = load_vec<3>(pos1, na);
            bool wrapped = false;
            // X dimension: wrap into [-L/2, +L/2]
            while (p.x() < -half_x) { p.x() += box_x; wrapped = true; }
            while (p.x() >= half_x) { p.x() -= box_x; wrapped = true; }
            // Y dimension: wrap into [-L/2, +L/2]
            while (p.y() < -half_y) { p.y() += box_y; wrapped = true; }
            while (p.y() >= half_y) { p.y() -= box_y; wrapped = true; }
            // Z dimension: wrap into [-L/2, +L/2]
            while (p.z() < -half_z) { p.z() += box_z; wrapped = true; }
            while (p.z() >= half_z) { p.z() -= box_z; wrapped = true; }
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



// Static method to apply PBC wrapping to positions
static bool apply_pbc_wrapping(VecArray pos, int n_atoms) {
    bool any_wrapped = false;
    float half_x = 0.5f * PeriodicBoundaryPotential::static_box_x;
    float half_y = 0.5f * PeriodicBoundaryPotential::static_box_y;
    float half_z = 0.5f * PeriodicBoundaryPotential::static_box_z;
    for(int na=0; na<n_atoms; ++na) {
        auto p = load_vec<3>(pos, na);
        bool wrapped = false;
        // X dimension: wrap into [-L/2, +L/2]
        while (p.x() < -half_x) { p.x() += PeriodicBoundaryPotential::static_box_x; wrapped = true; }
        while (p.x() >= half_x) { p.x() -= PeriodicBoundaryPotential::static_box_x; wrapped = true; }
        // Y dimension: wrap into [-L/2, +L/2]
        while (p.y() < -half_y) { p.y() += PeriodicBoundaryPotential::static_box_y; wrapped = true; }
        while (p.y() >= half_y) { p.y() -= PeriodicBoundaryPotential::static_box_y; wrapped = true; }
        // Z dimension: wrap into [-L/2, +L/2]
        while (p.z() < -half_z) { p.z() += PeriodicBoundaryPotential::static_box_z; wrapped = true; }
        while (p.z() >= half_z) { p.z() -= PeriodicBoundaryPotential::static_box_z; wrapped = true; }
        if (wrapped) {
            store_vec<3>(pos, na, p);
            any_wrapped = true;
        }
    }
    return any_wrapped;
}

// MARTINI potential using spline interpolation for LJ and Coulomb calculations
struct MartiniPotential : public PotentialNode
{
    int n_atom;
    CoordNode& pos;
    
    vector<array<float,4>> coeff;
    vector<pair<int,int>> pairs;
    
    float epsilon, sigma, lj_cutoff, coul_cutoff, dielectric;
    float coulomb_constant;
    bool force_cap;
    bool coulomb_soften;
    float slater_alpha;
    float mass_scale;  // Mass scaling factor for forces (1/mass)

    // Box dimensions for minimum image
    float box_x, box_y, box_z;
    
    // Spline interpolation for LJ potential
    LayeredClampedSpline1D<1> lj_potential_spline;
    LayeredClampedSpline1D<1> lj_force_spline;
    
    // Spline interpolation for Coulomb potential
    LayeredClampedSpline1D<1> coul_potential_spline;
    LayeredClampedSpline1D<1> coul_force_spline;
    
    // Spline parameters
    float lj_r_min, lj_r_max, lj_r_scale;
    int lj_n_knots;
    float coul_r_min, coul_r_max, coul_r_scale;
    int coul_n_knots;
    
    // Debug variables
    bool debug_mode;
    int debug_step_count;
    int max_debug_interactions;
    
    MartiniPotential(hid_t grp, CoordNode& pos_):
        PotentialNode(), n_atom(pos_.n_elem), pos(pos_), debug_step_count(0), max_debug_interactions(10),
        lj_potential_spline(1, 1000), lj_force_spline(1, 1000),
        coul_potential_spline(1, 1000), coul_force_spline(1, 1000)
    {
        check_size(grp, "atom_indices", n_atom);
        check_size(grp, "charges", n_atom);
        
        epsilon     = read_attribute<float>(grp, ".", "epsilon");
        sigma       = read_attribute<float>(grp, ".", "sigma");  
        lj_cutoff   = read_attribute<float>(grp, ".", "lj_cutoff");
        coul_cutoff = read_attribute<float>(grp, ".", "coul_cutoff");
        dielectric  = read_attribute<float>(grp, ".", "dielectric");
        
        // Read Coulomb constant - use converted value from Python if available, otherwise use unit conversion value
        coulomb_constant = 476.6f; // Default value (21.831807297541² for unit conversion with 4πϵ₀ = 1)
        if(attribute_exists(grp, ".", "coulomb_constant")) {
            coulomb_constant = read_attribute<float>(grp, ".", "coulomb_constant");
        }
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
        }
        
        // Debug mode - enable for first few steps
        debug_mode = false;
        if(attribute_exists(grp, ".", "debug_mode")) {
            debug_mode = read_attribute<int>(grp, ".", "debug_mode") != 0;
        }
        
        // Mass scaling for MARTINI force field (mass = 72.0)
        mass_scale = 1.0f / 72.0f;  // Default MARTINI mass
        if(attribute_exists(grp, ".", "mass_scale")) {
            mass_scale = read_attribute<float>(grp, ".", "mass_scale");
        }
        
        // Read box dimensions (same as MartiniPotential)
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
        
        // Print box size for debug
        std::cout << "[DEBUG] MartiniPotential box_x=" << box_x << " box_y=" << box_y << " box_z=" << box_z << std::endl;

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
        
        // Initialize spline parameters
        lj_r_min = 0.5f * sigma;  // Minimum distance for spline
        lj_r_max = lj_cutoff;     // Maximum distance for spline
        lj_r_scale = 999.0f / (lj_r_max - lj_r_min);
        
        coul_r_min = 0.5f;        // Minimum distance for Coulomb spline (Angstroms)
        coul_r_max = coul_cutoff; // Maximum distance for Coulomb spline
        coul_r_scale = 999.0f / (coul_r_max - coul_r_min);
        
        // Generate spline data for LJ potential
        std::vector<double> lj_pot_data(1000);
        std::vector<double> lj_force_data(1000);
        
        for(int i = 0; i < 1000; ++i) {
            float r = lj_r_min + i * (lj_r_max - lj_r_min) / 999.0f;
            float sig_r = sigma / r;
            float sig_r6 = sig_r * sig_r * sig_r * sig_r * sig_r * sig_r;
            float sig_r12 = sig_r6 * sig_r6;
            
            // LJ potential: 4*epsilon*((sigma/r)^12 - (sigma/r)^6)
            lj_pot_data[i] = 4.0 * epsilon * (sig_r12 - sig_r6);
            
            // LJ force derivative: 24*epsilon/sigma*(2*(sigma/r)^13-(sigma/r)^7)
            lj_force_data[i] = 24.0 * epsilon / sigma * (2.0 * sig_r12 * sig_r - sig_r6 * sig_r);
        }
        
        // Generate spline data for Coulomb potential (for unit charges)
        std::vector<double> coul_pot_data(1000);
        std::vector<double> coul_force_data(1000);
        
        for(int i = 0; i < 1000; ++i) {
            float r = coul_r_min + i * (coul_r_max - coul_r_min) / 999.0f;
            
            // Coulomb potential: k*q1*q2/(dielectric*r) for unit charges (q1=q2=1)
            coul_pot_data[i] = coulomb_constant / (dielectric * r);
            
            // Coulomb force derivative: k*q1*q2/(dielectric*r^2) for unit charges (positive for repulsion)
            coul_force_data[i] = coulomb_constant / (dielectric * r * r);
        }
        
        // Fit splines
        lj_potential_spline.fit_spline(lj_pot_data.data());
        lj_force_spline.fit_spline(lj_force_data.data());
        coul_potential_spline.fit_spline(coul_pot_data.data());
        coul_force_spline.fit_spline(coul_force_data.data());

        // Debug: Write all unique spline tables to a single file if debug_mode is enabled
        if (debug_mode) {
            std::ofstream out("all_splines.txt", std::ios::app);
            out << "# All spline values below are divided by 72.0 (mass correction)\n";
            // --- LJ splines for each unique (epsilon, sigma) ---
            std::set<std::pair<float, float>> lj_params;
            for (const auto& c : coeff) {
                if (c[0] != 0.f && c[1] != 0.f) lj_params.insert({c[0], c[1]});
            }
            for (const auto& p : lj_params) {
                float epsilon = p.first, sigma = p.second;
                out << "# LJ Spline\n# epsilon=" << epsilon << ", sigma=" << sigma << ", r_min=" << 0.5f * sigma << ", r_max=" << lj_cutoff << "\n";
                out << "# r potential force\n";
                int n_pts = 10;
                for (int i = 0; i < n_pts; ++i) {
                    float r = 0.5f * sigma + i * (lj_cutoff - 0.5f * sigma) / (n_pts - 1);
                    float sig_r = sigma / r;
                    float sig_r6 = pow(sig_r, 6);
                    float sig_r12 = sig_r6 * sig_r6;
                    float pot = 4.0 * epsilon * (sig_r12 - sig_r6) / 72.0;
                    float force = 24.0 * epsilon / sigma * (2.0 * sig_r12 * sig_r - sig_r6 * sig_r) / 72.0;
                    out << r << " " << pot << " " << force << "\n";
                }
                out << "\n";
            }
            // --- Coulomb splines for each unique q1*q2 ---
            std::set<float> qq_params;
            for (const auto& c : coeff) {
                float qq = c[2] * c[3];
                if (qq != 0.f) qq_params.insert(qq);
            }
            for (float qq : qq_params) {
                out << "# Coulomb Spline\n# q1q2=" << qq << ", dielectric=" << dielectric << ", coulomb_constant=" << coulomb_constant << ", r_min=0.5, r_max=" << coul_cutoff << "\n";
                out << "# r potential force\n";
                int n_pts = 10;
                for (int i = 0; i < n_pts; ++i) {
                    float r = 0.5f + i * (coul_cutoff - 0.5f) / (n_pts - 1);
                    float pot = coulomb_constant * qq / (dielectric * r) / 72.0;
                    float force = coulomb_constant * qq / (dielectric * r * r) / 72.0;
                    out << r << " " << pot << " " << force << "\n";
                }
                out << "\n";
            }
            out.close();
        }
        
        std::cout << "MARTINI: Initialized splines with 1000 knots" << std::endl;
        std::cout << "  LJ range: " << lj_r_min << " to " << lj_r_max << " Angstroms" << std::endl;
        std::cout << "  Coulomb range: " << coul_r_min << " to " << coul_r_max << " Angstroms" << std::endl;
        std::cout << "  Coulomb constant: " << coulomb_constant << ", Dielectric: " << dielectric << std::endl;
        std::cout << "  Mass scale: " << mass_scale << " (mass = " << (1.0f/mass_scale) << ")" << std::endl;
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("martini_potential"));
        
        VecArray pos1      = pos.output;
        VecArray pos1_sens = pos.sens;
        
        fill(pos1_sens, 3, n_atom, 0.f);
        
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;
        
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
            // Apply minimum image convention for rectangular box
            auto dr = minimum_image_rect(p1 - p2, box_x, box_y, box_z);
            auto dist2 = mag2(dr);
            auto dist = sqrtf(dist2);
            // DEBUG: Print pairwise info for first few steps
            if (debug_mode && debug_step_count < 5) {
                std::cout << "[DEBUG] Pair " << i << "-" << j << " dr=(" << dr.x() << "," << dr.y() << "," << dr.z() << ") dist=" << dist << std::endl;
            }
            
            // Apply distance cutoff for computational efficiency
            if(dist > max(lj_cutoff, coul_cutoff)) continue;
            
            Vec<3> force = make_zero<3>();
            
            // Lennard-Jones potential using spline interpolation
            if(eps != 0.f && sig != 0.f && dist < lj_cutoff && dist >= lj_r_min) {
                // Use spline interpolation for LJ potential and force
                float r_coord = (dist - lj_r_min) * lj_r_scale;
                
                float lj_result[2];
                lj_potential_spline.evaluate_value_and_deriv(lj_result, 0, r_coord);
                float lj_pot = lj_result[1];
                
                float force_result[2];
                lj_force_spline.evaluate_value_and_deriv(force_result, 0, r_coord);
                float lj_force_mag = force_result[1];
                
                // Only apply if potential and force are finite
                if(std::isfinite(lj_pot) && std::isfinite(lj_force_mag)) {
                    if(pot) *pot += lj_pot;
                    force += (lj_force_mag/dist) * (-dr);  // Use unit vector
                }
            }
            
            // Coulomb potential using spline interpolation
            if(qi != 0.f && qj != 0.f && dist < coul_cutoff && dist >= coul_r_min) {
                float coul_pot, coul_force_mag;
                
                if(coulomb_soften && slater_alpha > 0.0f) {
                    // Slater softened Coulomb potential (use original calculation for now)
                    // V(r) = q1*q2/r * (1 - (1 + αr/2) * exp(-αr))
                    float alpha_r = slater_alpha * dist;
                    float exp_alpha_r = expf(-alpha_r);
                    float soften_factor = 1.0f - (1.0f + alpha_r * 0.5f) * exp_alpha_r;
                    
                    coul_pot = coulomb_constant * qi * qj * soften_factor / (dielectric * dist);
                    
                    // Derivative: dV/dr = q1*q2 * [-(1 - (1 + αr/2) * exp(-αr))/r^2 + α/2 * exp(-αr)/r]
                    float d_soften_dr = slater_alpha * 0.5f * exp_alpha_r / dist;
                    coul_force_mag = -coulomb_constant * qi * qj * (soften_factor / dist2 - d_soften_dr) / dielectric;
                } else {
                    // Use spline interpolation for standard Coulomb potential
                    float r_coord = (dist - coul_r_min) * coul_r_scale;
                    
                    float coul_result[2];
                    coul_potential_spline.evaluate_value_and_deriv(coul_result, 0, r_coord);
                    coul_pot = coul_result[1] * qi * qj;  // Scale by actual charges
                    
                    float force_result[2];
                    coul_force_spline.evaluate_value_and_deriv(force_result, 0, r_coord);
                    coul_force_mag = force_result[1] * qi * qj;  // Scale by actual charges
                }
                
                // Only apply if potential and force are finite
                if(std::isfinite(coul_pot) && std::isfinite(coul_force_mag)) {
                    if(pot) *pot += coul_pot;
                    force += (coul_force_mag/dist) * (-dr);  // Fixed: add negative sign for correct direction
                }
            }
            
            // Apply mass scaling to forces (divide by mass)
            update_vec<3>(pos1_sens, i,  force * mass_scale);
            update_vec<3>(pos1_sens, j, -force * mass_scale);
        }
        
        debug_step_count++;
    }
};
static RegisterNodeType<MartiniPotential,1> martini_potential_node("martini_potential");

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
    LayeredClampedSpline1D<1> bond_force_spline;
    
    // Spline parameters
    float bond_r_min, bond_r_max, bond_r_scale;
    float max_spring;  // Store max spring constant for scaling
    float mass_scale;  // Mass scaling factor for forces (1/mass)
    bool debug_mode;   // Debug flag for writing splines

    DistSpring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(2, grp, "id")[0]), pos(pos_), params(n_elem),
        bond_potential_spline(1, 1000), bond_force_spline(1, 1000)
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
        
        // Mass scaling for MARTINI force field (mass = 72.0)
        mass_scale = 1.0f / 72.0f;  // Default MARTINI mass
        if(attribute_exists(grp, ".", "mass_scale")) {
            mass_scale = read_attribute<float>(grp, ".", "mass_scale");
        }
        // Debug mode for spline output
        debug_mode = false;
        if(attribute_exists(grp, ".", "debug_mode")) {
            debug_mode = read_attribute<int>(grp, ".", "debug_mode") != 0;
        }
        
        // Set spline range to cover typical bond length variations
        bond_r_min = std::max(0.1f, min_equil * 0.5f);  // Minimum distance
        bond_r_max = max_equil * 2.0f;                  // Maximum distance
        bond_r_scale = 999.0f / (bond_r_max - bond_r_min);
        
        // Generate spline data for bond potential (using average spring constant)
        std::vector<double> bond_pot_data(1000);
        std::vector<double> bond_force_data(1000);
        
        for(int i = 0; i < 1000; ++i) {
            float r = bond_r_min + i * (bond_r_max - bond_r_min) / 999.0f;
            
            // Use average equilibrium distance and spring constant for spline
            float avg_equil = (min_equil + max_equil) * 0.5f;
            float avg_spring = max_spring * 0.5f;  // Use half of max for safety
            
            // Bond potential: 0.5 * k * (r - r0)^2
            float delta_r = r - avg_equil;
            bond_pot_data[i] = 0.5 * avg_spring * delta_r * delta_r;
            
            // Bond force: k * (r - r0)
            bond_force_data[i] = avg_spring * delta_r;
        }
        
        // Fit splines
        bond_potential_spline.fit_spline(bond_pot_data.data());
        bond_force_spline.fit_spline(bond_force_data.data());

        // Debug: Write all unique bond splines to a single file if debug_mode is enabled
        if (debug_mode) {
            std::ofstream out("bond_splines.txt", std::ios::app);
            out << "# All spline values below are divided by 72.0 (mass correction)\n";
            // Collect unique (k, r0)
            std::set<std::pair<float, float>> bond_params;
            for (const auto& p : params) bond_params.insert({p.spring_constant, p.equil_dist});
            for (const auto& bp : bond_params) {
                float k = bp.first, r0 = bp.second;
                out << "# Bond Spline\n# k=" << k << ", r0=" << r0 << "\n";
                out << "# r potential force\n";
                int n_pts = 10;
                for (int i = 0; i < n_pts; ++i) {
                    float r = std::max(0.1f, r0 * 0.5f) + i * (r0 * 2.0f - std::max(0.1f, r0 * 0.5f)) / (n_pts - 1);
                    float pot = 0.5f * k * (r - r0) * (r - r0) / 72.0;
                    float force = k * (r - r0) / 72.0;
                    out << r << " " << pot << " " << force << "\n";
                }
                out << "\n";
            }
            out.close();
        }
        
        std::cout << "BONDS: Initialized splines with 1000 knots" << std::endl;
        std::cout << "  Bond range: " << bond_r_min << " to " << bond_r_max << " Angstroms" << std::endl;
        std::cout << "  Equilibrium range: " << min_equil << " to " << max_equil << " Angstroms" << std::endl;
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

            // Apply minimum image convention for periodic boundaries
            auto disp = minimum_image_rect(x1 - x2, box_x, box_y, box_z);
            float dist = mag(disp);
            
            // Use spline interpolation for bond potential and force
            if(dist >= bond_r_min && dist <= bond_r_max) {
                float r_coord = (dist - bond_r_min) * bond_r_scale;
                
                // Get potential and force from spline
                float pot_result[2];
                bond_potential_spline.evaluate_value_and_deriv(pot_result, 0, r_coord);
                float bond_pot = pot_result[1];
                
                float force_result[2];
                bond_force_spline.evaluate_value_and_deriv(force_result, 0, r_coord);
                float bond_force_mag = force_result[1];
                
                // Scale by actual spring constant and equilibrium distance
                float scale_factor = p.spring_constant / (max_spring * 0.5f);  // Normalize to spline
                float delta_r = dist - p.equil_dist;
                float actual_pot = 0.5f * p.spring_constant * delta_r * delta_r;
                float actual_force = p.spring_constant * delta_r;
                
                if(pot) *pot += actual_pot;
                
                // Apply force with mass scaling
                auto deriv = (actual_force/dist) * disp * mass_scale;
                update_vec(pos_sens, p.atom[0],  deriv);
                update_vec(pos_sens, p.atom[1], -deriv);
            } else {
                // Fallback to direct calculation for out-of-range distances
                auto deriv = p.spring_constant * (1.f - p.equil_dist*inv_mag(disp)) * disp * mass_scale;
                if(pot) *pot += 0.5f * p.spring_constant * sqr(mag(disp) - p.equil_dist);
                update_vec(pos_sens, p.atom[0],  deriv);
                update_vec(pos_sens, p.atom[1], -deriv);
            }
        }
    }
};
static RegisterNodeType<DistSpring,1> dist_spring_node("dist_spring");

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
    LayeredClampedSpline1D<1> angle_force_spline;
    
    // Spline parameters
    float angle_cos_min, angle_cos_max, angle_cos_scale;
    float max_spring;  // Store max spring constant for scaling
    float mass_scale;  // Mass scaling factor for forces (1/mass)
    bool debug_mode;   // Debug flag for writing splines

    AngleSpring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(2, grp, "id")[0]), pos(pos_), params(n_elem),
        angle_potential_spline(1, 1000), angle_force_spline(1, 1000)
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
        
        // Mass scaling for MARTINI force field (mass = 72.0)
        mass_scale = 1.0f / 72.0f;  // Default MARTINI mass
        if(attribute_exists(grp, ".", "mass_scale")) {
            mass_scale = read_attribute<float>(grp, ".", "mass_scale");
        }
        // Debug mode for spline output
        debug_mode = false;
        if(attribute_exists(grp, ".", "debug_mode")) {
            debug_mode = read_attribute<int>(grp, ".", "debug_mode") != 0;
        }
        
        // Set spline range for cosine of angles (cos ranges from -1 to 1)
        angle_cos_min = -1.0f;  // cos(180°) = -1
        angle_cos_max = 1.0f;   // cos(0°) = 1
        angle_cos_scale = 999.0f / (angle_cos_max - angle_cos_min);
        
        // Generate spline data for angle potential (using average spring constant)
        std::vector<double> angle_pot_data(1000);
        std::vector<double> angle_force_data(1000);
        
        for(int i = 0; i < 1000; ++i) {
            float cos_theta = angle_cos_min + i * (angle_cos_max - angle_cos_min) / 999.0f;
            
            // Use average equilibrium angle and spring constant for spline
            float avg_angle_deg = (min_angle + max_angle) * 0.5f;
            float avg_spring = max_spring * 0.5f;  // Use half of max for safety
            
            // Convert to radians and calculate cos
            float avg_angle_rad = avg_angle_deg * M_PI / 180.0f;
            float cos_avg = cosf(avg_angle_rad);
            
            // Angle potential: 0.5 * k * (cos(theta) - cos(theta0))^2
            float delta_cos = cos_theta - cos_avg;
            angle_pot_data[i] = 0.5 * avg_spring * delta_cos * delta_cos;
            
            // Angle force: k * (cos(theta) - cos(theta0))
            angle_force_data[i] = avg_spring * delta_cos;
        }
        
        // Fit splines
        angle_potential_spline.fit_spline(angle_pot_data.data());
        angle_force_spline.fit_spline(angle_force_data.data());

        // Debug: Write all unique angle splines to a single file if debug_mode is enabled
        if (debug_mode) {
            std::ofstream out("angle_splines.txt", std::ios::app);
            out << "# All spline values below are divided by 72.0 (mass correction)\n";
            // Collect unique (k, theta0)
            std::set<std::pair<float, float>> angle_params;
            for (const auto& p : params) angle_params.insert({p.spring_constant, p.equil_angle_deg});
            for (const auto& ap : angle_params) {
                float k = ap.first, theta0 = ap.second;
                out << "# Angle Spline\n# k=" << k << ", theta0_deg=" << theta0 << "\n";
                out << "# theta_deg potential force\n";
                int n_pts = 10;
                for (int i = 0; i < n_pts; ++i) {
                    float theta = 180.0f * i / (n_pts - 1);
                    float delta = cosf(theta * M_PI / 180.0f) - cosf(theta0 * M_PI / 180.0f);
                    float pot = 0.5f * k * delta * delta / 72.0;
                    float force = k * delta / 72.0;
                    out << theta << " " << pot << " " << force << "\n";
                }
                out << "\n";
            }
            out.close();
        }
        
        std::cout << "ANGLES: Initialized splines with 1000 knots" << std::endl;
        std::cout << "  Angle range: " << min_angle << " to " << max_angle << " degrees" << std::endl;
        std::cout << "  Cosine range: " << angle_cos_min << " to " << angle_cos_max << std::endl;
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

            // Apply minimum image convention for periodic boundaries
            auto disp1 = minimum_image_rect(make_vec3(atom1.x(), atom1.y(), atom1.z()) - make_vec3(atom2.x(), atom2.y(), atom2.z()), box_x, box_y, box_z);
            auto disp2 = minimum_image_rect(make_vec3(atom3.x(), atom3.y(), atom3.z()) - make_vec3(atom2.x(), atom2.y(), atom2.z()), box_x, box_y, box_z);

            // Normalize vectors
            float norm1 = sqrtf(disp1.x()*disp1.x() + disp1.y()*disp1.y() + disp1.z()*disp1.z());
            float norm2 = sqrtf(disp2.x()*disp2.x() + disp2.y()*disp2.y() + disp2.z()*disp2.z());
            if(norm1 == 0.f || norm2 == 0.f) continue; // avoid division by zero
            float dp = (disp1.x()*disp2.x() + disp1.y()*disp2.y() + disp1.z()*disp2.z()) / (norm1 * norm2);
            dp = std::max(-1.0f, std::min(1.0f, dp)); // clamp for safety
            
            // Use spline interpolation for angle potential and force
            if(dp >= angle_cos_min && dp <= angle_cos_max) {
                float cos_coord = (dp - angle_cos_min) * angle_cos_scale;
                
                // Get potential and force from spline
                float pot_result[2];
                angle_potential_spline.evaluate_value_and_deriv(pot_result, 0, cos_coord);
                float angle_pot = pot_result[1];
                
                float force_result[2];
                angle_force_spline.evaluate_value_and_deriv(force_result, 0, cos_coord);
                float angle_force_mag = force_result[1];
                
                // Scale by actual spring constant and equilibrium angle
                float scale_factor = p.spring_constant / (max_spring * 0.5f);  // Normalize to spline
                float theta0_rad = p.equil_angle_deg * M_PI / 180.0f;
                float cos_theta0 = cosf(theta0_rad);
                float delta_cos = dp - cos_theta0;
                float actual_pot = 0.5f * p.spring_constant * delta_cos * delta_cos;
                float actual_force = p.spring_constant * delta_cos;
                
                if(pot) *pot += actual_pot;
                
                // Force calculation (same as original)
                float theta_rad = acosf(dp);
                float dE_dtheta = actual_force * (-sinf(theta_rad));
                float dtheta_ddp = -1.0f / sqrtf(1.0f - dp*dp);
                float dE_ddp = dE_dtheta * dtheta_ddp;

                // Distribute forces to atoms
                float3 r1 = disp1;
                float3 r2 = disp2;
                float inv_norm1 = 1.0f / norm1;
                float inv_norm2 = 1.0f / norm2;
                float inv_norm1_sq = inv_norm1 * inv_norm1;
                float inv_norm2_sq = inv_norm2 * inv_norm2;
                float inv_norm1_norm2 = inv_norm1 * inv_norm2;
                float3 dcos_dr1 = (r2 - r1 * (dp * inv_norm1_sq)) * inv_norm1_norm2;
                float3 dcos_dr2 = (r1 - r2 * (dp * inv_norm2_sq)) * inv_norm1_norm2;
                float3 fa = dE_ddp * dcos_dr1;
                float3 fc = dE_ddp * dcos_dr2;
                float3 fb = -fa - fc;

                // Update forces with mass scaling
                pos_sens[4*p.atom[0]+0] += fa.x() * mass_scale;
                pos_sens[4*p.atom[0]+1] += fa.y() * mass_scale;
                pos_sens[4*p.atom[0]+2] += fa.z() * mass_scale;
                pos_sens[4*p.atom[1]+0] += fb.x() * mass_scale;
                pos_sens[4*p.atom[1]+1] += fb.y() * mass_scale;
                pos_sens[4*p.atom[1]+2] += fb.z() * mass_scale;
                pos_sens[4*p.atom[2]+0] += fc.x() * mass_scale;
                pos_sens[4*p.atom[2]+1] += fc.y() * mass_scale;
                pos_sens[4*p.atom[2]+2] += fc.z() * mass_scale;
            } else {
                // Fallback to direct calculation for out-of-range angles
                float theta_rad = acosf(dp);
                float theta0_rad = p.equil_angle_deg * M_PI / 180.0f;
                float cos_theta0_rad = cosf(theta0_rad);
                float delta_cos = dp - cos_theta0_rad;

                if(pot) *pot += 0.5f * p.spring_constant * delta_cos * delta_cos;

                float dE_dtheta = p.spring_constant * delta_cos * (-sinf(theta_rad));
                float dtheta_ddp = -1.0f / sqrtf(1.0f - dp*dp);
                float dE_ddp = dE_dtheta * dtheta_ddp;

                float3 r1 = disp1;
                float3 r2 = disp2;
                float inv_norm1 = 1.0f / norm1;
                float inv_norm2 = 1.0f / norm2;
                float inv_norm1_sq = inv_norm1 * inv_norm1;
                float inv_norm2_sq = inv_norm2 * inv_norm2;
                float inv_norm1_norm2 = inv_norm1 * inv_norm2;
                float3 dcos_dr1 = (r2 - r1 * (dp * inv_norm1_sq)) * inv_norm1_norm2;
                float3 dcos_dr2 = (r1 - r2 * (dp * inv_norm2_sq)) * inv_norm1_norm2;
                float3 fa = dE_ddp * dcos_dr1;
                float3 fc = dE_ddp * dcos_dr2;
                float3 fb = -fa - fc;

                pos_sens[4*p.atom[0]+0] += fa.x() * mass_scale;
                pos_sens[4*p.atom[0]+1] += fa.y() * mass_scale;
                pos_sens[4*p.atom[0]+2] += fa.z() * mass_scale;
                pos_sens[4*p.atom[1]+0] += fb.x() * mass_scale;
                pos_sens[4*p.atom[1]+1] += fb.y() * mass_scale;
                pos_sens[4*p.atom[1]+2] += fb.z() * mass_scale;
                pos_sens[4*p.atom[2]+0] += fc.x() * mass_scale;
                pos_sens[4*p.atom[2]+1] += fc.y() * mass_scale;
                pos_sens[4*p.atom[2]+2] += fc.z() * mass_scale;
            }
        }
    }
};
static RegisterNodeType<AngleSpring,1> angle_spring_node("angle_spring");

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


