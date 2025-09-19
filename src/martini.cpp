#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include "spline.h"
#include <iostream>
#include <H5Apublic.h> // for H5Aexists
#include <fstream> // For file writing
#include <cmath> // For pow, cosf, sinf, acosf
#include <set> // For std::set
#include <complex> // For complex numbers in PME
#include <vector> // For PME grid operations
#include <algorithm> // For PME algorithms

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
    
    // Spline interpolation for dihedral potential
    LayeredClampedSpline1D<1> dihedral_potential_spline;
    
    // Spline parameters
    float dihedral_min, dihedral_max;
    float max_spring;  // Store max spring constant for scaling
    bool debug_mode;   // Debug flag for writing splines

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
        

        
        // Debug mode for spline output
        debug_mode = false;
        if(attribute_exists(grp, ".", "debug_mode")) {
            debug_mode = read_attribute<int>(grp, ".", "debug_mode") != 0;
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
        
        // Set spline range for dihedral angles (typically -π to π radians)
        dihedral_min = -M_PI_F;  // -180°
        dihedral_max = M_PI_F;   // +180°
        
        // Use a parameterless canonical spline: pot = 0.5 * (delta_phi)^2
        // with delta_phi in [-pi, pi]
        std::vector<double> dihedral_pot_data(1000);
        for(int i = 0; i < 1000; ++i) {
            float delta_phi = dihedral_min + i * (dihedral_max - dihedral_min) / 999.0f;
            dihedral_pot_data[i] = 0.5 * delta_phi * delta_phi;
        }
        
        // Fit spline
        dihedral_potential_spline.fit_spline(dihedral_pot_data.data());

        // Debug: Write all unique dihedral splines to a single file if debug_mode is enabled
        if (debug_mode) {
            std::ofstream out("dihedral_splines.txt", std::ios::out | std::ios::trunc);
            out << "# Canonical dihedral spline (delta_phi = phi - phi0)\n";
            out << "# Forces are calculated as analytical derivatives of the potential\n";
            // Collect unique (k, phi0)
            std::set<std::pair<float, float>> dihedral_params;
            for (const auto& p : params) dihedral_params.insert({p.spring_constant, p.equil_dihedral});
            for (const auto& dp : dihedral_params) {
                float k = dp.first, phi0 = dp.second;
                out << "# Dihedral Spline\n# k=" << k << ", phi0_rad=" << phi0 << "\n";
                out << "# phi_rad potential\n";
                int n_pts = 10;
                for (int i = 0; i < n_pts; ++i) {
                    float phi = -M_PI_F + i * 2.0f * M_PI_F / (n_pts - 1);
                    float delta_phi = phi - phi0;
                    if(delta_phi > M_PI_F) delta_phi -= 2.0f * M_PI_F;
                    if(delta_phi < -M_PI_F) delta_phi += 2.0f * M_PI_F;
                    float pot = 0.5f * k * delta_phi * delta_phi;
                    out << phi << " " << pot << "\n";
                }
                out << "\n";
            }
            out.close();
        }
        
        std::cout << "DIHEDRALS: Initialized splines with 1000 knots" << std::endl;
        std::cout << "  Dihedral range: " << min_dihedral << " to " << max_dihedral << " radians" << std::endl;
        std::cout << "  Spline range: " << dihedral_min << " to " << dihedral_max << " radians" << std::endl;
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

            // Use spline interpolation for dihedral potential and force
            // Use delta phi relative to equilibrium
            float delta_phi = dihedral - p.equil_dihedral;
            if(delta_phi > M_PI_F) delta_phi -= 2.0f * M_PI_F;
            if(delta_phi < -M_PI_F) delta_phi += 2.0f * M_PI_F;
            if(delta_phi >= dihedral_min && delta_phi <= dihedral_max) {
                float phi_coord = delta_phi;
                
                // Get potential and force from single spline
                float result[2];
                dihedral_potential_spline.evaluate_value_and_deriv(result, 0, phi_coord);
                float dihedral_pot = result[1]; // Index 1 is the value
                float dE_ddelta = result[0]; // Index 0 is the derivative w.r.t. delta_phi
                float dihedral_force_mag = p.spring_constant * dE_ddelta;
                
                // Use spline-evaluated potential and force
                if(pot) *pot += p.spring_constant * dihedral_pot;

                // Apply force with mass scaling (dihedral_force_mag corresponds to dE/dphi)
                auto s = Float4(dihedral_force_mag);
                for(int na: range(4)) d[na].scale_update(s, pos_sens + 4*params[nt].atom[na]);
            } else {
                // Fallback to direct calculation for out-of-range dihedrals
                float displacement = dihedral - p.equil_dihedral;
                displacement = (displacement> M_PI_F) ? displacement-2.f*M_PI_F : displacement;
                displacement = (displacement<-M_PI_F) ? displacement+2.f*M_PI_F : displacement;

                if(pot) *pot += 0.5f * p.spring_constant * sqr(displacement);
                
                auto s = Float4(p.spring_constant * displacement);
                for(int na: range(4)) d[na].scale_update(s, pos_sens + 4*params[nt].atom[na]);
            }
        }
    }
    
    // Method to update box dimensions for NPT barostat
    void update_box_dimensions(float scale_factor) override {
        box_x *= scale_factor;
        box_y *= scale_factor;
        box_z *= scale_factor;
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
        
        // Initialize static variables for global access
        static_box_x = box_x;
        static_box_y = box_y;
        static_box_z = box_z;
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



// Static member variable definitions
float PeriodicBoundaryPotential::static_box_x = 0.0f;
float PeriodicBoundaryPotential::static_box_y = 0.0f;
float PeriodicBoundaryPotential::static_box_z = 0.0f;

// Function to update box dimensions for NPT barostat
void update_box_dimensions(float scale_factor) {
    // Update static box dimensions used by PeriodicBoundaryPotential
    PeriodicBoundaryPotential::static_box_x *= scale_factor;
    PeriodicBoundaryPotential::static_box_y *= scale_factor;
    PeriodicBoundaryPotential::static_box_z *= scale_factor;
    
    // Note: Individual potential classes (MartiniPotential, DistSpring, AngleSpring)
    // also store their own box dimensions that need to be updated.
    // This requires access to the potential instances, which is not available here.
    // The proper solution would be to add an update_box_dimensions method to each
    // potential class and call it from the main simulation loop.
}

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

// Simple 1D FFT implementation for PME
class SimpleFFT {
public:
    static void fft1d(std::vector<std::complex<float>>& data) {
        int n = data.size();
        if (n <= 1) return;
        
        // Cooley-Tukey FFT algorithm
        std::vector<std::complex<float>> even(n/2), odd(n/2);
        for (int i = 0; i < n/2; ++i) {
            even[i] = data[2*i];
            odd[i] = data[2*i + 1];
        }
        
        fft1d(even);
        fft1d(odd);
        
        for (int i = 0; i < n/2; ++i) {
            std::complex<float> t = std::polar(1.0f, -2.0f * M_PI_F * i / n) * odd[i];
            data[i] = even[i] + t;
            data[i + n/2] = even[i] - t;
        }
    }
    
    static void ifft1d(std::vector<std::complex<float>>& data) {
        // Conjugate the input
        for (auto& x : data) {
            x = std::conj(x);
        }
        
        // Apply forward FFT
        fft1d(data);
        
        // Conjugate and normalize
        int n = data.size();
        for (auto& x : data) {
            x = std::conj(x) / static_cast<float>(n);
        }
    }
    
    static void fft3d(std::vector<std::vector<std::vector<std::complex<float>>>>& data) {
        int nx = data.size();
        int ny = data[0].size();
        int nz = data[0][0].size();
        
        // FFT along x direction
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                std::vector<std::complex<float>> x_slice(nx);
                for (int i = 0; i < nx; ++i) {
                    x_slice[i] = data[i][j][k];
                }
                fft1d(x_slice);
                for (int i = 0; i < nx; ++i) {
                    data[i][j][k] = x_slice[i];
                }
            }
        }
        
        // FFT along y direction
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nz; ++k) {
                std::vector<std::complex<float>> y_slice(ny);
                for (int j = 0; j < ny; ++j) {
                    y_slice[j] = data[i][j][k];
                }
                fft1d(y_slice);
                for (int j = 0; j < ny; ++j) {
                    data[i][j][k] = y_slice[j];
                }
            }
        }
        
        // FFT along z direction
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                std::vector<std::complex<float>> z_slice(nz);
                for (int k = 0; k < nz; ++k) {
                    z_slice[k] = data[i][j][k];
                }
                fft1d(z_slice);
                for (int k = 0; k < nz; ++k) {
                    data[i][j][k] = z_slice[k];
                }
            }
        }
    }
    
    static void ifft3d(std::vector<std::vector<std::vector<std::complex<float>>>>& data) {
        int nx = data.size();
        int ny = data[0].size();
        int nz = data[0][0].size();
        
        // IFFT along z direction
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                std::vector<std::complex<float>> z_slice(nz);
                for (int k = 0; k < nz; ++k) {
                    z_slice[k] = data[i][j][k];
                }
                ifft1d(z_slice);
                for (int k = 0; k < nz; ++k) {
                    data[i][j][k] = z_slice[k];
                }
            }
        }
        
        // IFFT along y direction
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nz; ++k) {
                std::vector<std::complex<float>> y_slice(ny);
                for (int j = 0; j < ny; ++j) {
                    y_slice[j] = data[i][j][k];
                }
                ifft1d(y_slice);
                for (int j = 0; j < ny; ++j) {
                    data[i][j][k] = y_slice[j];
                }
            }
        }
        
        // IFFT along x direction
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                std::vector<std::complex<float>> x_slice(nx);
                for (int i = 0; i < nx; ++i) {
                    x_slice[i] = data[i][j][k];
                }
                ifft1d(x_slice);
                for (int i = 0; i < nx; ++i) {
                    data[i][j][k] = x_slice[i];
                }
            }
        }
    }
};

// B-spline interpolation for PME charge assignment
class BSplineInterpolation {
public:
    static float bspline3(float x) {
        // Third-order B-spline
        x = fabsf(x);
        if (x < 1.0f) {
            return (2.0f/3.0f) - x*x + 0.5f*x*x*x;
        } else if (x < 2.0f) {
            float t = 2.0f - x;
            return t*t*t / 6.0f;
        }
        return 0.0f;
    }
    
    static void assign_charges_to_grid(const std::vector<Vec<3,float>>& positions,
                                     const std::vector<float>& charges,
                                     std::vector<std::vector<std::vector<float>>>& grid,
                                     float box_x, float box_y, float box_z,
                                     int nx, int ny, int nz) {
        // Clear grid
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    grid[i][j][k] = 0.0f;
                }
            }
        }
        
        float dx = box_x / nx;
        float dy = box_y / ny;
        float dz = box_z / nz;
        
        // Assign charges to grid using B-spline interpolation
        for (size_t p = 0; p < positions.size(); ++p) {
            const auto& pos = positions[p];
            float q = charges[p];
            
            // Convert to grid coordinates
            float gx = (pos.x() + 0.5f * box_x) / dx;
            float gy = (pos.y() + 0.5f * box_y) / dy;
            float gz = (pos.z() + 0.5f * box_z) / dz;
            
            // Find grid points within B-spline support
            int ix_start = static_cast<int>(floorf(gx - 2.0f));
            int iy_start = static_cast<int>(floorf(gy - 2.0f));
            int iz_start = static_cast<int>(floorf(gz - 2.0f));
            
            for (int ix = ix_start; ix <= ix_start + 4; ++ix) {
                for (int iy = iy_start; iy <= iy_start + 4; ++iy) {
                    for (int iz = iz_start; iz <= iz_start + 4; ++iz) {
                        // Apply periodic boundary conditions
                        int i = ((ix % nx) + nx) % nx;
                        int j = ((iy % ny) + ny) % ny;
                        int k = ((iz % nz) + nz) % nz;
                        
                        // Calculate B-spline weights
                        float wx = bspline3(gx - ix);
                        float wy = bspline3(gy - iy);
                        float wz = bspline3(gz - iz);
                        
                        // Add charge contribution
                        grid[i][j][k] += q * wx * wy * wz;
                    }
                }
            }
        }
    }
    
    static void interpolate_forces_from_grid(const std::vector<Vec<3,float>>& positions,
                                           const std::vector<std::vector<std::vector<Vec<3,float>>>>& force_grid,
                                           std::vector<Vec<3,float>>& forces,
                                           float box_x, float box_y, float box_z,
                                           int nx, int ny, int nz) {
        float dx = box_x / nx;
        float dy = box_y / ny;
        float dz = box_z / nz;
        
        for (size_t p = 0; p < positions.size(); ++p) {
            const auto& pos = positions[p];
            Vec<3,float> force = make_zero<3>();
            
            // Convert to grid coordinates
            float gx = (pos.x() + 0.5f * box_x) / dx;
            float gy = (pos.y() + 0.5f * box_y) / dy;
            float gz = (pos.z() + 0.5f * box_z) / dz;
            
            // Find grid points within B-spline support
            int ix_start = static_cast<int>(floorf(gx - 2.0f));
            int iy_start = static_cast<int>(floorf(gy - 2.0f));
            int iz_start = static_cast<int>(floorf(gz - 2.0f));
            
            for (int ix = ix_start; ix <= ix_start + 4; ++ix) {
                for (int iy = iy_start; iy <= iy_start + 4; ++iy) {
                    for (int iz = iz_start; iz <= iz_start + 4; ++iz) {
                        // Apply periodic boundary conditions
                        int i = ((ix % nx) + nx) % nx;
                        int j = ((iy % ny) + ny) % ny;
                        int k = ((iz % nz) + nz) % nz;
                        
                        // Calculate B-spline weights
                        float wx = bspline3(gx - ix);
                        float wy = bspline3(gy - iy);
                        float wz = bspline3(gz - iz);
                        
                        // Add force contribution
                        force += force_grid[i][j][k] * wx * wy * wz;
                    }
                }
            }
            
            forces[p] = force;
        }
    }
};

// Particle Mesh Ewald (PME) implementation for long-range Coulomb interactions
struct ParticleMeshEwald : public PotentialNode
{
    int n_atom;
    CoordNode& pos;
    
    // Box dimensions
    float box_x, box_y, box_z;
    float volume;
    
    // PME parameters
    float alpha;           // Ewald parameter (screening parameter)
    float rcut;            // Real space cutoff
    int nx, ny, nz;        // Grid dimensions
    int order;             // B-spline interpolation order
    
    // Charge data
    vector<float> charges;
    
    // PME grids
    vector<vector<vector<float>>> charge_grid;
    vector<vector<vector<complex<float>>>> potential_grid;
    vector<vector<vector<Vec<3,float>>>> force_grid;
    
    // Precomputed factors
    float coulomb_k;       // Coulomb constant
    float alpha_sqrt_pi;   // α/sqrt(π)
    
    ParticleMeshEwald(hid_t grp, CoordNode& pos_):
        PotentialNode(), n_atom(pos_.n_elem), pos(pos_), coulomb_k(31.775347952181f)
    {
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
        
        volume = box_x * box_y * box_z;
        
        // Read charges
        check_size(grp, "charges", n_atom);
        charges.resize(n_atom);
        traverse_dset<1,float>(grp, "charges", [&](size_t i, float q) {
            charges[i] = q;
        });
        
        // PME parameters
        alpha = 0.2f;  // Default Ewald parameter
        if(attribute_exists(grp, ".", "pme_alpha")) {
            alpha = read_attribute<float>(grp, ".", "pme_alpha");
        }
        
        rcut = 10.0f;  // Default real space cutoff
        if(attribute_exists(grp, ".", "pme_rcut")) {
            rcut = read_attribute<float>(grp, ".", "pme_rcut");
        }
        
        // Grid dimensions (should be powers of 2 for efficient FFT)
        nx = 32;  // Default grid size
        if(attribute_exists(grp, ".", "pme_nx")) {
            nx = read_attribute<int>(grp, ".", "pme_nx");
        }
        
        ny = 32;
        if(attribute_exists(grp, ".", "pme_ny")) {
            ny = read_attribute<int>(grp, ".", "pme_ny");
        }
        
        nz = 32;
        if(attribute_exists(grp, ".", "pme_nz")) {
            nz = read_attribute<int>(grp, ".", "pme_nz");
        }
        
        order = 4;  // B-spline interpolation order
        if(attribute_exists(grp, ".", "pme_order")) {
            order = read_attribute<int>(grp, ".", "pme_order");
        }
        
        // Precompute factors
        alpha_sqrt_pi = alpha / sqrtf(M_PI);
        
        // Initialize grids
        charge_grid.resize(nx, vector<vector<float>>(ny, vector<float>(nz, 0.0f)));
        potential_grid.resize(nx, vector<vector<complex<float>>>(ny, vector<complex<float>>(nz, 0.0f)));
        force_grid.resize(nx, vector<vector<Vec<3,float>>>(ny, vector<Vec<3,float>>(nz, make_zero<3>())));
        
        std::cout << "PME: Initialized with alpha=" << alpha << ", rcut=" << rcut << std::endl;
        std::cout << "PME: Box dimensions: " << box_x << " x " << box_y << " x " << box_z << std::endl;
        std::cout << "PME: Grid dimensions: " << nx << " x " << ny << " x " << nz << std::endl;
        std::cout << "PME: B-spline order: " << order << std::endl;
    }
    
    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("particle_mesh_ewald"));
        
        VecArray pos1 = pos.output;
        VecArray pos1_sens = pos.sens;
        
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;
        
        // Extract positions
        vector<Vec<3,float>> positions(n_atom);
        for(int i = 0; i < n_atom; ++i) {
            positions[i] = load_vec<3>(pos1, i);
        }
        
        // Step 1: Assign charges to grid using B-spline interpolation
        BSplineInterpolation::assign_charges_to_grid(positions, charges, charge_grid, 
                                                    box_x, box_y, box_z, nx, ny, nz);
        
        // Step 2: Convert charge grid to complex format for FFT
        for(int i = 0; i < nx; ++i) {
            for(int j = 0; j < ny; ++j) {
                for(int k = 0; k < nz; ++k) {
                    potential_grid[i][j][k] = complex<float>(charge_grid[i][j][k], 0.0f);
                }
            }
        }
        
        // Step 3: Forward 3D FFT
        SimpleFFT::fft3d(potential_grid);
        
        // Step 4: Solve Poisson's equation in reciprocal space
        float dx = box_x / nx;
        float dy = box_y / ny;
        float dz = box_z / nz;
        
        for(int i = 0; i < nx; ++i) {
            for(int j = 0; j < ny; ++j) {
                for(int k = 0; k < nz; ++k) {
                    // Calculate k-vector components
                    float kx = (i <= nx/2) ? 2.0f * M_PI * i / box_x : 2.0f * M_PI * (i - nx) / box_x;
                    float ky = (j <= ny/2) ? 2.0f * M_PI * j / box_y : 2.0f * M_PI * (j - ny) / box_y;
                    float kz = (k <= nz/2) ? 2.0f * M_PI * k / box_z : 2.0f * M_PI * (k - nz) / box_z;
                    
                    float k_sq = kx*kx + ky*ky + kz*kz;
                    
                    if(k_sq > 0.0f) {
                        // Green's function: G(k) = 1/k² * exp(-k²/(4α²))
                        float green = expf(-k_sq / (4.0f * alpha * alpha)) / k_sq;
                        potential_grid[i][j][k] *= green;
                    } else {
                        // Handle k=0 case (should be zero for neutral systems)
                        potential_grid[i][j][k] = 0.0f;
                    }
                }
            }
        }
        
        // Step 5: Inverse 3D FFT
        SimpleFFT::ifft3d(potential_grid);
        
        // Step 6: Calculate potential energy
        float recip_pot = 0.0f;
        for(int i = 0; i < nx; ++i) {
            for(int j = 0; j < ny; ++j) {
                for(int k = 0; k < nz; ++k) {
                    recip_pot += charge_grid[i][j][k] * potential_grid[i][j][k].real();
                }
            }
        }
        recip_pot *= 0.5f * coulomb_k / volume;
        
        // Step 7: Calculate forces on grid
        for(int i = 0; i < nx; ++i) {
            for(int j = 0; j < ny; ++j) {
                for(int k = 0; k < nz; ++k) {
                    // Calculate gradients using finite differences
                    int ip = (i + 1) % nx;
                    int im = (i - 1 + nx) % nx;
                    int jp = (j + 1) % ny;
                    int jm = (j - 1 + ny) % ny;
                    int kp = (k + 1) % nz;
                    int km = (k - 1 + nz) % nz;
                    
                    float fx = (potential_grid[ip][j][k].real() - potential_grid[im][j][k].real()) / (2.0f * dx);
                    float fy = (potential_grid[i][jp][k].real() - potential_grid[i][jm][k].real()) / (2.0f * dy);
                    float fz = (potential_grid[i][j][kp].real() - potential_grid[i][j][km].real()) / (2.0f * dz);
                    
                    force_grid[i][j][k] = make_vec3(-fx, -fy, -fz) * coulomb_k;
                }
            }
        }
        
        // Step 8: Interpolate forces back to particles
        vector<Vec<3,float>> recip_forces(n_atom, make_zero<3>());
        BSplineInterpolation::interpolate_forces_from_grid(positions, force_grid, recip_forces,
                                                          box_x, box_y, box_z, nx, ny, nz);
        
        // Step 9: Self-interaction correction
        float self_correction = 0.0f;
        for(int i = 0; i < n_atom; ++i) {
            self_correction += charges[i] * charges[i];
        }
        self_correction *= alpha_sqrt_pi;
        
        // Total reciprocal space potential
        recip_pot -= self_correction;
        
        if(pot) *pot += recip_pot;
        
        // Apply forces
        for(int i = 0; i < n_atom; ++i) {
            update_vec<3>(pos1_sens, i, -recip_forces[i]);
        }
    }
    
    // Method to update box dimensions for NPT barostat
    void update_box_dimensions(float scale_factor) override {
        box_x *= scale_factor;
        box_y *= scale_factor;
        box_z *= scale_factor;
        volume = box_x * box_y * box_z;
    }
};
static RegisterNodeType<ParticleMeshEwald, 1> pme_node("particle_mesh_ewald");

// MARTINI potential using spline interpolation for LJ and Coulomb calculations
struct MartiniPotential : public PotentialNode
{
    int n_atom;
    CoordNode& pos;
    
    vector<array<float,4>> coeff;
    vector<pair<int,int>> pairs;
    
    float epsilon, sigma, lj_cutoff, coul_cutoff;
    bool force_cap;
    bool coulomb_soften;
    float slater_alpha;
    bool lj_soften;
    float lj_soften_alpha;
    bool overwrite_spline_tables;
    
    // PME parameters
    bool use_pme;
    float pme_alpha;
    float pme_rcut;

    // Box dimensions for minimum image
    float box_x, box_y, box_z;
    
    // Spline interpolation for LJ potential - single spline for each epsilon/sigma pair
    std::map<std::pair<float, float>, LayeredClampedSpline1D<1>> lj_splines;
    
    // Spline interpolation for Coulomb potential - single spline for each charge product
    std::map<float, LayeredClampedSpline1D<1>> coulomb_splines;
    
    // Spline parameters
    float lj_r_min, lj_r_max;
    float lj_r_shift, lj_r_scale;  // Coordinate transformation parameters
    float coul_r_min, coul_r_max;
    float coul_r_shift, coul_r_scale;  // Coordinate transformation parameters
    int coul_n_knots;
    
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
        // dielectric constant is now included in the Coulomb k constant (31.775347952181)
        
        // Coulomb constant is now hardcoded as 31.775347952181 in the potential calculation
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
        
        // LJ softening parameters (soft-core LJ)
        lj_soften = false;
        if(attribute_exists(grp, ".", "lj_soften")) {
            lj_soften = read_attribute<int>(grp, ".", "lj_soften") != 0;
        }
        lj_soften_alpha = 0.0f;
        if(lj_soften) {
            if(attribute_exists(grp, ".", "lj_soften_alpha")) {
                lj_soften_alpha = read_attribute<float>(grp, ".", "lj_soften_alpha");
            } else {
                lj_soften_alpha = 1.0f; // sensible default, dimensionless added to (r/sigma)^6
            }
        }

        // Debug mode - enable for first few steps
        debug_mode = false;
        if(attribute_exists(grp, ".", "debug_mode")) {
            debug_mode = read_attribute<int>(grp, ".", "debug_mode") != 0;
        }
        
        // Optionally overwrite existing spline table output files (debug text files)
        overwrite_spline_tables = false;
        if(attribute_exists(grp, ".", "overwrite_spline_tables")) {
            overwrite_spline_tables = read_attribute<int>(grp, ".", "overwrite_spline_tables") != 0;
        }
        
        // PME parameters
        use_pme = false;
        if(attribute_exists(grp, ".", "use_pme")) {
            use_pme = read_attribute<int>(grp, ".", "use_pme") != 0;
        }
        
        pme_alpha = 0.2f;  // Default PME parameter
        if(attribute_exists(grp, ".", "pme_alpha")) {
            pme_alpha = read_attribute<float>(grp, ".", "pme_alpha");
        }
        
        pme_rcut = 10.0f;  // Default real space cutoff
        if(attribute_exists(grp, ".", "pme_rcut")) {
            pme_rcut = read_attribute<float>(grp, ".", "pme_rcut");
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
        
        // Check if we have the optimized format (coefficient_indices) or the original format
        bool optimized_format = false;
        if(attribute_exists(grp, ".", "optimized_format")) {
            optimized_format = read_attribute<int>(grp, ".", "optimized_format") != 0;
        }

        // Declare unique_coeff at function scope so it's available later
        std::vector<array<float,4>> unique_coeff;

        if(optimized_format) {
            // Optimized format: unique coefficients + indices
            std::cout << "MARTINI: Using optimized interaction table format" << std::endl;
            
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
            
            // Reconstruct full coefficient array from indices
            coeff.resize(n_pair);
            for(size_t np = 0; np < n_pair; ++np) {
                long idx = coeff_indices[np];
                if(idx >= 0 && idx < (long)n_unique_coeff) {
                    coeff[np] = unique_coeff[idx];
                } else {
                    std::cerr << "ERROR: Invalid coefficient index " << idx << " for pair " << np << std::endl;
                    coeff[np] = {0.0f, 0.0f, 0.0f, 0.0f};
                }
            }
            
            // Debug: Print first few coefficients to verify reconstruction
            std::cout << "MARTINI: Debug - First 5 reconstructed coefficients:" << std::endl;
            for(size_t i = 0; i < std::min(size_t(5), size_t(n_pair)); ++i) {
                std::cout << "  Pair " << i << " (atoms " << pairs[i].first << "-" << pairs[i].second 
                          << "): idx=" << (int)coeff_indices[i] << " coeff=[" 
                          << coeff[i][0] << ", " << coeff[i][1] << ", " << coeff[i][2] << ", " << coeff[i][3] << "]" << std::endl;
            }
            
            std::cout << "MARTINI: Loaded " << n_unique_coeff << " unique coefficients for " << n_pair << " pairs" << std::endl;
            std::cout << "MARTINI: Compression ratio: " << (float)n_pair / n_unique_coeff << "x" << std::endl;
            
        } else {
            // Original format: full coefficient array
            std::cout << "MARTINI: Using original interaction table format" << std::endl;
            
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
        
        // Find all epsilon/sigma pairs for separate LJ splines
        std::set<std::pair<float, float>> unique_lj_params;

        if(optimized_format) {
            // In optimized format, iterate over unique coefficients
            for(const auto& c : unique_coeff) {
                float eps = c[0];
                float sig = c[1];
                if(eps != 0.f && sig != 0.f) {
                    unique_lj_params.insert({eps, sig});
                }
            }
        } else {
            // In original format, iterate over all coefficients
            for(const auto& c : coeff) {
                float eps = c[0];
                float sig = c[1];
                if(eps != 0.f && sig != 0.f) {
                    unique_lj_params.insert({eps, sig});
                }
            }
        }

        std::cout << "MARTINI: Generating separate LJ splines for " << unique_lj_params.size() << " unique epsilon/sigma pairs from coefficients array" << std::endl;
        for(const auto& params : unique_lj_params) {
            std::cout << "  epsilon=" << params.first << ", sigma=" << params.second << std::endl;
        }

        // Initialize spline parameters for LJ - fixed domain [0, 12]
        lj_r_min = 0.0f;
        lj_r_max = 12.0f;  // Fixed domain for all LJ interactions
        
        // Coordinate transformation parameters (following membrane spline pattern)
        lj_r_shift = -lj_r_min;  // shift to make coordinate system start at 0
        lj_r_scale = 999.0f / (lj_r_max + lj_r_shift);  // scale to [0, 999]

        coul_r_min = 0.0f;        // Minimum distance for Coulomb spline (Angstroms)
        coul_r_max = 12.0f;       // Fixed domain for all Coulomb interactions
        
        // Coordinate transformation parameters for Coulomb
        coul_r_shift = -coul_r_min;  // shift to make coordinate system start at 0
        coul_r_scale = 999.0f / (coul_r_max + coul_r_shift);  // scale to [0, 999]

        // Generate separate LJ splines for each unique epsilon/sigma pair
        for(const auto& params : unique_lj_params) {
            float eps = params.first;
            float sig = params.second;

            std::vector<double> lj_pot_data(1000);

            for(int i = 0; i < 1000; ++i) {
                float r = lj_r_min + i * (lj_r_max - lj_r_min) / 999.0f;
                if(r == 0.0f) r = 1.0e-6f;

                // Check for numerical stability - avoid very small r values for regular LJ
                if(!lj_soften || lj_soften_alpha <= 0.0f) {
                    if(r < 0.1f * sig) r = 0.1f * sig;  // Minimum r = 0.1 * sigma to avoid numerical issues
                }

                if(lj_soften && lj_soften_alpha > 0.0f) {
                    // Soft-core LJ: t = (r/sigma)^6 + alpha; V = 4*epsilon*(1/t^2 - 1/t)
                    float x = r / sig;
                    float x2 = x * x;
                    float x3 = x2 * x;
                    float x6 = x3 * x3; // (r/sigma)^6
                    float t = x6 + lj_soften_alpha;
                    float inv_t = 1.0f / t;
                    float inv_t2 = inv_t * inv_t;
                    // potential for softened Lennard-Jones
                    lj_pot_data[i] = 4.0 * eps * (inv_t2 - inv_t);
                } else {
                    // Regular LJ potential: V = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)
                    float r2 = r * r;
                    float r3 = r2 * r;
                    float r6 = r3 * r3;
                    float sig2 = sig * sig;
                    float sig6 = sig2 * sig2 * sig2;
                    float sig12 = sig6 * sig6;
                    float inv_r6 = sig6 / r6;
                    float inv_r12 = sig12 / (r6 * r6);
                    // LJ potential
                    lj_pot_data[i] = 4.0 * eps * (inv_r12 - inv_r6);
                }
            }

            // Create single spline for this epsilon/sigma combination
            auto [it, inserted] = lj_splines.emplace(std::piecewise_construct,
                                                     std::forward_as_tuple(eps, sig),
                                                     std::forward_as_tuple(LayeredClampedSpline1D<1>(1, 1000)));
            auto& spline = it->second;

            // Only compute spline data if this is a new parameter set
            if (inserted) {
                // Initialize the spline with potential data
                spline.fit_spline(lj_pot_data.data());
            }
            
        }
        

        // Generate separate Coulomb splines for each unique charge product
        std::set<float> unique_charge_products;
        if(optimized_format) {
            // In optimized format, look at unique coefficients
            for(const auto& c : unique_coeff) {
                float qq = c[2] * c[3];
                if(std::abs(qq) > 1e-10f) {  // Use small epsilon for floating point comparison
                    unique_charge_products.insert(qq);
                }
            }
        } else {
            // In original format, look at all coefficients
            for(const auto& c : coeff) {
                float qq = c[2] * c[3];
                if(std::abs(qq) > 1e-10f) {  // Use small epsilon for floating point comparison
                    unique_charge_products.insert(qq);
                }
            }
        }

        std::cout << "MARTINI: Generating separate Coulomb splines for " << unique_charge_products.size() << " unique charge products" << std::endl;

        // Generate separate Coulomb splines for each unique charge product
        for(float qq : unique_charge_products) {
            std::vector<double> coul_pot_data_for_spline(1000 * 1);  // 1 layer, 1000 points, 1 value per point

            for(int i = 0; i < 1000; ++i) {
                float r = coul_r_min + i * (coul_r_max - coul_r_min) / 999.0f;
                if(r == 0.0f) r = 1.0e-6f;

                // Coulomb potential: V = k * qq / r, where k = 31.775347952181 (includes epsilon_r=15)
                float coulomb_k = 31.775347952181f;
                float potential = coulomb_k * qq / r;

                // Apply softening if enabled
                if(coulomb_soften) {
                    // Slater softening: V(r) = k*qq/r * (1 - (1 + αr/2) * exp(-αr))
                    float alpha_r = slater_alpha * r;
                    float exp_term = expf(-alpha_r);
                    float soft_factor = 1.0f - (1.0f + alpha_r * 0.5f) * exp_term;

                    // Softened potential
                    coul_pot_data_for_spline[i] = potential * soft_factor;
                } else {
                    coul_pot_data_for_spline[i] = potential;
                }
            }

            // Create single spline for this charge product
            auto [coulomb_it, coulomb_inserted] = coulomb_splines.emplace(std::piecewise_construct,
                                                                          std::forward_as_tuple(qq),
                                                                          std::forward_as_tuple(LayeredClampedSpline1D<1>(1, 1000)));
            auto& coulomb_spline = coulomb_it->second;

            // Initialize the spline with potential data
            coulomb_spline.fit_spline(coul_pot_data_for_spline.data());
        }

        std::cout << "MARTINI: Generated " << coulomb_splines.size() << " Coulomb splines" << std::endl;

        // Debug: Write all unique spline tables to a single file
        {
            static bool s_truncated_all_splines = false;
            std::ios_base::openmode mode = std::ios::app;
            if (overwrite_spline_tables && !s_truncated_all_splines) {
                mode = std::ios::out | std::ios::trunc;
                s_truncated_all_splines = true;
            }
            std::ofstream out("all_splines.txt", mode);
            out << "# LJ splines: Separate tables for each unique epsilon/sigma pair\n";
            out << "# Each spline contains only the potential - forces calculated as analytical derivatives\n\n";
            out << "# Coulomb splines: Separate tables for each unique charge product\n";
            out << "# Each spline contains only the potential - forces calculated as analytical derivatives\n\n";
            out << "# Bond splines: Harmonic potential for bond distances\n";
            out << "# Angle splines: Cosine-based potential for bond angles\n";
            out << "# Dihedral splines: Harmonic potential for dihedral angles\n";
            out << "# All forces are calculated as analytical derivatives of the potential splines\n";

            // --- LJ splines for each unique (epsilon, sigma) ---
            for (const auto& spline_pair : lj_splines) {
                float epsilon = spline_pair.first.first;
                float sigma = spline_pair.first.second;
                const auto& spline = spline_pair.second;

                out << "# LJ Spline\n# epsilon=" << epsilon << ", sigma=" << sigma << ", r_min=" << lj_r_min << ", r_max=" << lj_r_max
                    << ", softened=" << (lj_soften?1:0) << ", lj_soften_alpha=" << lj_soften_alpha << "\n";
                out << "# r potential\n";

                int n_pts = 10;
                for (int i = 0; i < n_pts; ++i) {
                    float r = lj_r_min + i * (lj_r_max - lj_r_min) / (n_pts - 1);
                    if(r == 0.0f) r = 1.0e-6f;
                    float r_coord = (r - lj_r_min) / (lj_r_max - lj_r_min) * 999.0f;

                    float result[2];
                    spline.evaluate_value_and_deriv(result, 0, r_coord);
                    float pot = result[1];  // Index 1 is the value

                    out << r << " " << pot << "\n";
                }
                out << "\n";
            }
            // --- Coulomb splines for each unique charge product ---
            for (const auto& coulomb_pair : coulomb_splines) {
                float qq = coulomb_pair.first;
                const auto& spline = coulomb_pair.second;

                out << "# Coulomb Spline\n# q1q2=" << qq << ", k=31.775347952181, r_min=" << coul_r_min << ", r_max=" << coul_r_max << ", softened=" << (coulomb_soften?1:0) << ", slater_alpha=" << slater_alpha << "\n";
                out << "# r potential\n";

                int n_pts = 10;
                for (int i = 0; i < n_pts; ++i) {
                    float r = coul_r_min + i * (coul_r_max - coul_r_min) / (n_pts - 1);
                    if(r == 0.0f) r = 1.0e-6f;
                    float r_coord = (r - coul_r_min) / (coul_r_max - coul_r_min) * 999.0f;

                    float result[2];
                    spline.evaluate_value_and_deriv(result, 0, r_coord);
                    float pot = result[1];  // Index 1 is the value

                    out << r << " " << pot << "\n";
                }
                out << "\n";
            }
            
            // Add bond and angle splines to the main output file
            // Read and append bond splines if they exist
            std::ifstream bond_file("bond_splines.txt");
            if (bond_file.is_open()) {
                out << "\n# === BOND SPLINES ===\n";
                std::string line;
                while (std::getline(bond_file, line)) {
                    out << line << "\n";
                }
                bond_file.close();
            }
            
            // Read and append angle splines if they exist
            std::ifstream angle_file("angle_splines.txt");
            if (angle_file.is_open()) {
                out << "\n# === ANGLE SPLINES ===\n";
                std::string line;
                while (std::getline(angle_file, line)) {
                    out << line << "\n";
                }
                angle_file.close();
            }
            
            out.close();
        }

        std::cout << "MARTINI: Initialized splines with 1000 knots" << std::endl;
        std::cout << "  LJ range: " << lj_r_min << " to " << lj_r_max << " Angstroms" << std::endl;
        std::cout << "  Coulomb range: " << coul_r_min << " to " << coul_r_max << " Angstroms" << std::endl;
        std::cout << "  Coulomb k: 31.775347952181 (includes epsilon_r=15)" << std::endl;
        if(use_pme) {
            std::cout << "  PME enabled: alpha=" << pme_alpha << ", rcut=" << pme_rcut << std::endl;
        } else {
            std::cout << "  PME disabled: using standard Coulomb cutoff" << std::endl;
        }

    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("martini_potential"));
        
        VecArray pos1      = pos.output;
        VecArray pos1_sens = pos.sens;
        
        // --- REMOVED: fill(pos1_sens, 3, n_atom, 0.f); ---
        // This line was incorrectly zeroing the force array, erasing all bonded forces
        // Force array initialization is now handled centrally in DerivEngine::compute()
        
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
            // if (debug_mode && debug_step_count < 5) {
            //     std::cout << "[DEBUG] Pair " << i << "-" << j << " dr=(" << dr.x() << "," << dr.y() << "," << dr.z() << ") dist=" << dist << std::endl;
            // }
            
            // Apply distance cutoff for computational efficiency
            if(dist > max(lj_cutoff, coul_cutoff)) continue;
            
            Vec<3> force = make_zero<3>();
            
            // Lennard-Jones potential using single spline for each epsilon/sigma pair
            if(eps != 0.f && sig != 0.f && dist < lj_cutoff) {
                // Look up the appropriate spline for this epsilon/sigma pair
                auto spline_it = lj_splines.find({eps, sig});
                if(spline_it != lj_splines.end()) {
                    // Use spline interpolation for LJ potential and force
                    // Transform physical distance to spline coordinate [0, 999]
                    float r_coord = (dist - lj_r_min) / (lj_r_max - lj_r_min) * 999.0f;

                    float lj_result[2];
                    // CORRECT: Call evaluate_value_and_deriv only on the potential spline
                    spline_it->second.evaluate_value_and_deriv(lj_result, 0, r_coord);

                    float lj_pot = lj_result[1];           // Index 1 is the value
                    float lj_deriv_spline = lj_result[0];   // Index 0 is the derivative w.r.t. spline coordinate

                    // Convert derivative from spline coordinate to physical coordinate (dE/dr)
                    // dE/dr = dE/d(coord) * d(coord)/dr
                    float coord_scale = 999.0f / (lj_r_max - lj_r_min);
                    float dE_dr = lj_deriv_spline * coord_scale;
                    
                    // The force is the negative gradient: F = -dE/dr
                    float lj_force_mag = -dE_dr;

                    // Only apply if potential and force are finite and reasonable
                    if(std::isfinite(lj_pot) && std::isfinite(lj_force_mag)) {
                        if(pot) *pot += lj_pot;
                        // lj_force_mag is now correctly calculated as -dV/dr (physical force)
                        // Accumulate the physical force directly
                        force += (lj_force_mag/dist) * dr;
                    }
                }
            }
            
            // Coulomb potential using single spline for each charge product
            if(qi != 0.f && qj != 0.f && dist < coul_cutoff) {
                float qq = qi * qj;
                // DEBUG: Check charge product
                // if(debug_step_count < 5) std::cout << "COULOMB: Looking for q1*q2=" << qq << " (qi=" << qi << ", qj=" << qj << ")" << std::endl;

                float coul_pot = 0.0f;
                float coul_force_mag = 0.0f;
                
                if(use_pme && dist < pme_rcut) {
                    // Real space part of PME: V = k * qq * erfc(α*r) / r
                    float coulomb_k = 31.775347952181f;
                    float alpha_r = pme_alpha * dist;
                    float erfc_alpha_r = erfc(alpha_r);
                    
                    // Potential: V = k * qq * erfc(α*r) / r
                    coul_pot = coulomb_k * qq * erfc_alpha_r / dist;
                    
                    // Force: F = -k * qq * (erfc(α*r)/r² + 2*α*exp(-α²r²)/(√π*r))
                    float exp_term = 2.0f * pme_alpha * expf(-alpha_r * alpha_r) / sqrtf(M_PI);
                    coul_force_mag = -coulomb_k * qq * (erfc_alpha_r / (dist * dist) + exp_term / dist);
                } else {
                    // Standard Coulomb potential (short-range only or when Ewald is disabled)
                    // Look up the appropriate spline for this charge product
                    auto coulomb_it = coulomb_splines.find(qq);
                    if(coulomb_it != coulomb_splines.end()) {
                        // DEBUG: Uncomment to check Coulomb computation
                        // std::cout << "COULOMB: Using spline for q1*q2=" << qq << " at r=" << dist << std::endl;
                        // Use spline interpolation for Coulomb potential and force
                        // Transform physical distance to spline coordinate [0, 999]
                        float r_coord = (dist - coul_r_min) / (coul_r_max - coul_r_min) * 999.0f;

                        float coul_result[2];
                        // CORRECT: Call evaluate_value_and_deriv only on the potential spline
                        coulomb_it->second.evaluate_value_and_deriv(coul_result, 0, r_coord);

                        coul_pot = coul_result[1];           // Index 1 is the value
                        float coul_deriv_spline = coul_result[0];   // Index 0 is the derivative w.r.t. spline coordinate

                        // Convert derivative from spline coordinate to physical coordinate (dE/dr)
                        // dE/dr = dE/d(coord) * d(coord)/dr
                        float coord_scale = 999.0f / (coul_r_max - coul_r_min);
                        float dE_dr = coul_deriv_spline * coord_scale;
                        
                        // The force is the negative gradient: F = -dE/dr
                        coul_force_mag = -dE_dr;
                    }
                }

                // Only apply if potential and force are finite and reasonable
                if(std::isfinite(coul_pot) && std::isfinite(coul_force_mag)) {
                    if(pot) *pot += coul_pot;
                    // coul_force_mag is now correctly calculated as -dV/dr (physical force)
                    // Accumulate the physical force directly
                    force += (coul_force_mag/dist) * dr;
                }
            }
            
            // Apply mass scaling to forces (divide by mass)
            // Store gradient (∇E = -F) in pos_sens for UPSIDE integrator
            update_vec<3>(pos1_sens, i, -force);
            update_vec<3>(pos1_sens, j,  force);
        }
        
        debug_step_count++;
    }
    
    // Method to update box dimensions for NPT barostat
    void update_box_dimensions(float scale_factor) override {
        box_x *= scale_factor;
        box_y *= scale_factor;
        box_z *= scale_factor;
    }
};
static RegisterNodeType<MartiniPotential, 1> martini_potential_node("martini_potential");

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
    bool debug_mode;   // Debug flag for writing splines

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
        

        // Debug mode for spline output
        debug_mode = false;
        if(attribute_exists(grp, ".", "debug_mode")) {
            debug_mode = read_attribute<int>(grp, ".", "debug_mode") != 0;
        }
        
        // Define the bond r-range to start at 0 and extend to r_max
        bond_r_min = 0.0f;
        // Prefer an attribute if provided; otherwise choose a conservative default
        if(attribute_exists(grp, ".", "bond_r_max")) {
            bond_r_max = read_attribute<float>(grp, ".", "bond_r_max");
        } else {
            // heuristic default: at least 2x the largest equilibrium distance or 5 Å
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

        // Debug: Write all unique bond splines to a single file if debug_mode is enabled
        if (debug_mode) {
            std::ofstream out("bond_splines.txt", std::ios::out | std::ios::trunc);
            out << "# Bond splines: Two-particle interactions (atom1-atom2)\n";
            out << "# Forces are calculated as analytical derivatives of the potential\n";
            // Collect unique (k, r0)
            std::set<std::pair<float, float>> bond_params;
            for (const auto& p : params) bond_params.insert({p.spring_constant, p.equil_dist});
            for (const auto& bp : bond_params) {
                float k = bp.first, r0 = bp.second;
                out << "# Bond Spline\n# k=" << k << ", r0=" << r0 << "\n";
                out << "# r potential\n";
                int n_pts = 10;
                for (int i = 0; i < n_pts; ++i) {
                    float r = std::max(0.1f, r0 * 0.5f) + i * (r0 * 2.0f - std::max(0.1f, r0 * 0.5f)) / (n_pts - 1);
                    float pot = 0.5f * k * (r - r0) * (r - r0);
                    out << r << " " << pot << "\n";
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
            
            // Use spline interpolation for bond potential and force (in delta-r space)
            float delta_r = dist - p.equil_dist;
            if(delta_r >= bond_delta_min && delta_r <= bond_delta_max) {
                float r_coord = delta_r;
                
                // Get potential and force from single spline
                float result[2];
                bond_potential_spline.evaluate_value_and_deriv(result, 0, r_coord);
                float bond_pot = result[1]; // Index 1 is the value
                float bond_deriv = result[0]; // Index 0 is the derivative w.r.t. delta_r
                
                // Scale by spring constant to match parameterized harmonic
                if(pot) *pot += p.spring_constant * bond_pot;
                
                // Apply force: F = -dE/dr, here dE/dr = k * dE/ddelta * ddelta/dr with ddelta/dr = 1
                // The derivative is already dE/d(delta_r), so we use it directly
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
    
    // Method to update box dimensions for NPT barostat
    void update_box_dimensions(float scale_factor) override {
        box_x *= scale_factor;
        box_y *= scale_factor;
        box_z *= scale_factor;
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
    bool debug_mode;   // Debug flag for writing splines

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
        

        // Debug mode for spline output
        debug_mode = false;
        if(attribute_exists(grp, ".", "debug_mode")) {
            debug_mode = read_attribute<int>(grp, ".", "debug_mode") != 0;
        }
        
        // Set spline range for delta_cos = cos(θ) - cos(θ₀)
        // Find the range of cos(equilibrium angles) to set proper delta_cos bounds
        float min_cos_equil = 1.0f;  // cos(0°) = 1
        float max_cos_equil = -1.0f; // cos(180°) = -1
        for(const auto& p : params) {
            float cos_equil = cosf(p.equil_angle_deg * M_PI / 180.0f);
            min_cos_equil = std::min(min_cos_equil, cos_equil);
            max_cos_equil = std::max(max_cos_equil, cos_equil);
        }

        // Set conservative bounds for delta_cos
        angle_cos_min = -2.0f;  // Allow for angles deviating significantly from equilibrium
        angle_cos_max = 2.0f;

        // Canonical spline for angles: pot = 0.5*(delta_cos)^2
        std::vector<double> angle_pot_data(1000);
        for(int i = 0; i < 1000; ++i) {
            float delta_cos = angle_cos_min + i * (angle_cos_max - angle_cos_min) / 999.0f;
            angle_pot_data[i] = 0.5 * delta_cos * delta_cos;
        }
        
        // Fit spline
        angle_potential_spline.fit_spline(angle_pot_data.data());

        // Debug: Write all unique angle splines to a single file if debug_mode is enabled
        if (debug_mode) {
            std::ofstream out("angle_splines.txt", std::ios::out | std::ios::trunc);
            out << "# Angle splines: Three-particle interactions (atom1-atom2-atom3)\n";
            out << "# Forces are calculated as analytical derivatives of the potential\n";
            // Collect unique (k, theta0)
            std::set<std::pair<float, float>> angle_params;
            for (const auto& p : params) angle_params.insert({p.spring_constant, p.equil_angle_deg});
            for (const auto& ap : angle_params) {
                float k = ap.first, theta0 = ap.second;
                out << "# Angle Spline\n# k=" << k << ", theta0_deg=" << theta0 << "\n";
                out << "# theta_deg potential\n";
                int n_pts = 10;
                for (int i = 0; i < n_pts; ++i) {
                    float theta = 180.0f * i / (n_pts - 1);
                    float cos_theta = cosf(theta * M_PI / 180.0f);
                    float cos_theta0 = cosf(theta0 * M_PI / 180.0f);
                    float delta_cos = cos_theta - cos_theta0;
                    float pot = 0.5f * k * delta_cos * delta_cos;
                    
                    out << theta << " " << pot << "\n";
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
            // Calculate delta_cos = cos(θ) - cos(θ₀)
            float cos_equil = cosf(p.equil_angle_deg * M_PI / 180.0f);
            float delta_cos = dp - cos_equil;
            if(delta_cos >= angle_cos_min && delta_cos <= angle_cos_max) {
                float cos_coord = delta_cos;
                
                // Get potential and force from single spline
                float result[2];
                angle_potential_spline.evaluate_value_and_deriv(result, 0, cos_coord);
                float angle_pot = result[1]; // Index 1 is the value
                float angle_deriv = result[0]; // Index 0 is the derivative w.r.t. delta_cos
                
                // Use spline-evaluated potential
                if(pot) *pot += p.spring_constant * angle_pot;

                // Convert spline dE/d(cos) to forces on atoms
                // The derivative is dE/d(delta_cos), so we use it directly
                float dE_ddp = p.spring_constant * angle_deriv;

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

                // Update forces
                pos_sens[4*p.atom[0]+0] += fa.x();
                pos_sens[4*p.atom[0]+1] += fa.y();
                pos_sens[4*p.atom[0]+2] += fa.z();
                pos_sens[4*p.atom[1]+0] += fb.x();
                pos_sens[4*p.atom[1]+1] += fb.y();
                pos_sens[4*p.atom[1]+2] += fb.z();
                pos_sens[4*p.atom[2]+0] += fc.x();
                pos_sens[4*p.atom[2]+1] += fc.y();
                pos_sens[4*p.atom[2]+2] += fc.z();
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

                pos_sens[4*p.atom[0]+0] += fa.x();
                pos_sens[4*p.atom[0]+1] += fa.y();
                pos_sens[4*p.atom[0]+2] += fa.z();
                pos_sens[4*p.atom[1]+0] += fb.x();
                pos_sens[4*p.atom[1]+1] += fb.y();
                pos_sens[4*p.atom[1]+2] += fb.z();
                pos_sens[4*p.atom[2]+0] += fc.x();
                pos_sens[4*p.atom[2]+1] += fc.y();
                pos_sens[4*p.atom[2]+2] += fc.z();
            }
        }
    }
    
    // Method to update box dimensions for NPT barostat
    void update_box_dimensions(float scale_factor) override {
        box_x *= scale_factor;
        box_y *= scale_factor;
        box_z *= scale_factor;
    }
};
static RegisterNodeType<AngleSpring, 1> angle_spring_node("angle_spring");

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

static RegisterNodeType<ConjugateGradientMinimizer, 1> cg_minimizer_node("conjugate_gradient_minimizer");



