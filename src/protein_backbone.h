#ifndef PROTEIN_BACKBONE_H
#define PROTEIN_BACKBONE_H

#include <vector>
#include <string>
#include <map>
#include <utility>

namespace ProteinBackbone {

// Calculate center of mass of backbone atoms (N, CA, C, O)
void calculate_backbone_center(const std::vector<double>& positions, 
                              const std::vector<int>& bb_atom_indices,
                              double* bb_center);

// Calculate forces on backbone atoms from external forces on virtual particle
void distribute_backbone_forces(const double* virtual_force,
                               const std::vector<int>& bb_atom_indices,
                               std::vector<double>& forces);

// Calculate LJ force between two particles
void calculate_lj_force(const double* pos1, const double* pos2,
                       double epsilon, double sigma, double cutoff,
                       double* force);

// Calculate Coulomb force between two charged particles
void calculate_coulomb_force(const double* pos1, const double* pos2,
                           double q1, double q2, double cutoff, double coulomb_constant,
                           double* force);

// Main function to calculate protein backbone forces
void calculate_protein_backbone_forces(const std::vector<double>& positions,
                                     const std::vector<double>& velocities,
                                     const std::vector<int>& backbone_groups,
                                     const std::vector<std::string>& backbone_bead_types,
                                     const std::vector<double>& all_charges,
                                     const std::vector<std::string>& all_atom_types,
                                     const std::map<std::pair<std::string, std::string>, std::pair<double, double>>& lj_params,
                                     double lj_cutoff, double coulomb_cutoff, double coulomb_constant,
                                     std::vector<double>& forces);

// Logging function for debugging
void log_backbone_calculation(int bb_idx, const double* bb_center, 
                            const std::string& bb_type, double total_force[3]);

} // namespace ProteinBackbone

#endif // PROTEIN_BACKBONE_H

