#include "protein_backbone.h"
#include <vector>
#include <cmath>
#include <iostream>

namespace ProteinBackbone {

// Calculate center of mass of backbone atoms (N, CA, C, O)
void calculate_backbone_center(const std::vector<double>& positions, 
                              const std::vector<int>& bb_atom_indices,
                              double* bb_center) {
    bb_center[0] = 0.0;
    bb_center[1] = 0.0;
    bb_center[2] = 0.0;
    
    for (int atom_idx : bb_atom_indices) {
        bb_center[0] += positions[atom_idx * 3];
        bb_center[1] += positions[atom_idx * 3 + 1];
        bb_center[2] += positions[atom_idx * 3 + 2];
    }
    
    // Divide by number of atoms (4 for backbone: N, CA, C, O)
    double inv_n = 1.0 / bb_atom_indices.size();
    bb_center[0] *= inv_n;
    bb_center[1] *= inv_n;
    bb_center[2] *= inv_n;
}

// Calculate forces on backbone atoms from external forces on virtual particle
void distribute_backbone_forces(const double* virtual_force,
                               const std::vector<int>& bb_atom_indices,
                               std::vector<double>& forces) {
    // Distribute force equally among all backbone atoms
    double force_per_atom[3] = {
        virtual_force[0] / bb_atom_indices.size(),
        virtual_force[1] / bb_atom_indices.size(),
        virtual_force[2] / bb_atom_indices.size()
    };
    
    for (int atom_idx : bb_atom_indices) {
        forces[atom_idx * 3] += force_per_atom[0];
        forces[atom_idx * 3 + 1] += force_per_atom[1];
        forces[atom_idx * 3 + 2] += force_per_atom[2];
    }
}

// Calculate LJ force between two particles
void calculate_lj_force(const double* pos1, const double* pos2,
                       double epsilon, double sigma, double cutoff,
                       double* force) {
    double dx = pos2[0] - pos1[0];
    double dy = pos2[1] - pos1[1];
    double dz = pos2[2] - pos1[2];
    
    double r2 = dx*dx + dy*dy + dz*dz;
    double r = std::sqrt(r2);
    
    if (r > cutoff || r < 1e-6) {
        force[0] = force[1] = force[2] = 0.0;
        return;
    }
    
    double inv_r = 1.0 / r;
    double sigma_r = sigma * inv_r;
    double sigma_r6 = sigma_r * sigma_r * sigma_r * sigma_r * sigma_r * sigma_r;
    double sigma_r12 = sigma_r6 * sigma_r6;
    
    double force_magnitude = 24.0 * epsilon * inv_r * (2.0 * sigma_r12 - sigma_r6);
    
    force[0] = force_magnitude * dx * inv_r;
    force[1] = force_magnitude * dy * inv_r;
    force[2] = force_magnitude * dz * inv_r;
}

// Calculate Coulomb force between two charged particles
void calculate_coulomb_force(const double* pos1, const double* pos2,
                           double q1, double q2, double cutoff, double coulomb_constant,
                           double* force) {
    double dx = pos2[0] - pos1[0];
    double dy = pos2[1] - pos1[1];
    double dz = pos2[2] - pos1[2];
    
    double r2 = dx*dx + dy*dy + dz*dz;
    double r = std::sqrt(r2);
    
    if (r > cutoff || r < 1e-6) {
        force[0] = force[1] = force[2] = 0.0;
        return;
    }
    
    double inv_r = 1.0 / r;
    double inv_r3 = inv_r * inv_r * inv_r;
    
    double force_magnitude = coulomb_constant * q1 * q2 * inv_r3;
    
    force[0] = force_magnitude * dx;
    force[1] = force_magnitude * dy;
    force[2] = force_magnitude * dz;
}

// Main function to calculate protein backbone forces
void calculate_protein_backbone_forces(const std::vector<double>& positions,
                                     const std::vector<double>& velocities,
                                     const std::vector<int>& backbone_groups,
                                     const std::vector<std::string>& backbone_bead_types,
                                     const std::vector<double>& all_charges,
                                     const std::vector<std::string>& all_atom_types,
                                     const std::map<std::pair<std::string, std::string>, std::pair<double, double>>& lj_params,
                                     double lj_cutoff, double coulomb_cutoff, double coulomb_constant,
                                     std::vector<double>& forces) {
    
    int n_backbones = backbone_groups.size() / 4; // Each backbone has 4 atoms
    
    for (int bb_idx = 0; bb_idx < n_backbones; bb_idx++) {
        // Get backbone atom indices
        std::vector<int> bb_atoms(4);
        for (int i = 0; i < 4; i++) {
            bb_atoms[i] = backbone_groups[bb_idx * 4 + i];
        }
        
        // Calculate virtual backbone center
        double bb_center[3];
        calculate_backbone_center(positions, bb_atoms, bb_center);
        
        // Get backbone bead type
        std::string bb_type = backbone_bead_types[bb_idx];
        
        // Calculate forces from virtual particle to all other particles
        for (size_t other_idx = 0; other_idx < positions.size() / 3; other_idx++) {
            // Skip if this is one of the backbone atoms
            bool is_backbone_atom = false;
            for (int bb_atom : bb_atoms) {
                if (other_idx == static_cast<size_t>(bb_atom)) {
                    is_backbone_atom = true;
                    break;
                }
            }
            if (is_backbone_atom) continue;
            
            // Get other particle position and type
            double other_pos[3] = {
                positions[other_idx * 3],
                positions[other_idx * 3 + 1],
                positions[other_idx * 3 + 2]
            };
            
            std::string other_type = all_atom_types[other_idx];
            double other_charge = all_charges[other_idx];
            
            // Calculate LJ force
            double lj_force[3];
            auto lj_key = std::make_pair(bb_type, other_type);
            if (lj_params.find(lj_key) != lj_params.end()) {
                auto params = lj_params.at(lj_key);
                double epsilon = params.first;
                double sigma = params.second;
                calculate_lj_force(bb_center, other_pos, epsilon, sigma, lj_cutoff, lj_force);
            } else {
                lj_force[0] = lj_force[1] = lj_force[2] = 0.0;
            }
            
            // Calculate Coulomb force
            double coulomb_force[3];
            double bb_charge = 0.0; // Backbone particles are typically neutral
            calculate_coulomb_force(bb_center, other_pos, bb_charge, other_charge, 
                                  coulomb_cutoff, coulomb_constant, coulomb_force);
            
            // Total force on virtual particle
            double total_force[3] = {
                lj_force[0] + coulomb_force[0],
                lj_force[1] + coulomb_force[1],
                lj_force[2] + coulomb_force[2]
            };
            
            // Distribute force to backbone atoms
            distribute_backbone_forces(total_force, bb_atoms, forces);
            
            // Apply equal and opposite force to other particle
            forces[other_idx * 3] -= total_force[0];
            forces[other_idx * 3 + 1] -= total_force[1];
            forces[other_idx * 3 + 2] -= total_force[2];
        }
    }
}

// Logging function for debugging
void log_backbone_calculation(int bb_idx, const double* bb_center, 
                            const std::string& bb_type, double total_force[3]) {
    std::cout << "Backbone " << bb_idx << " (type: " << bb_type 
              << ") center: [" << bb_center[0] << ", " << bb_center[1] << ", " << bb_center[2] << "]"
              << " total force: [" << total_force[0] << ", " << total_force[1] << ", " << total_force[2] << "]"
              << std::endl;
}

} // namespace ProteinBackbone

