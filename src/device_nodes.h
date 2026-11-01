#ifndef DEVICE_NODES_H
#define DEVICE_NODES_H

// The device implementations of the nodes (device_<family>.cu), made by each CPU node's make_device from its
// constructed data.  Index arrays come already offset into their input's concatenated array (DeviceGraph); an
// element's system, where a node needs it, comes in `system`.  A potential node adds its energy to `potential`
// (DeviceGraph::potential_of) at its element's system in PotentialAndDerivMode.

#include "device_engine.h"
#include "interaction_graph.h"
#include "spline.h"
#include <vector>

//! An InteractionGraph's data for every system, as device_interaction_graph.h's DeviceInteractionGraph is built
//! from: positions offset into the inputs' concatenated arrays, the other arrays concatenated over the systems.  For
//! a symmetric graph the second set is the first, and its arrays are empty.
struct DeviceInteractionData {
    DeviceCoord* pos1;
    DeviceCoord* pos2;
    std::vector<int> n_elem1, n_elem2;  // per system
    std::vector<int> loc1, loc2, type1, type2, id1, id2, sequence_id1, sequence_id2;
    std::vector<float> interaction_param;
    int n_type1, n_type2, max_n_edge;
    float cutoff, cache_buffer;
};

//! The DeviceInteractionData of the same InteractionGraph of every system; the systems must share its parameters
template <typename IType>
DeviceInteractionData device_interaction_data(DeviceGraph& g, const std::vector<const InteractionGraph<IType>*>& igs) {
    auto& ig0 = *igs.at(0);
    DeviceInteractionData d;
    d.pos1 = &g.coord(*ig0.pos_node1);
    d.pos2 = IType::symmetric ? nullptr : &g.coord(*ig0.pos_node2);
    for(int ns=0; ns<g.n_system; ++ns) {
        auto& ig = *igs[ns];
        if(ig.get_param() != ig0.get_param() || ig.cutoff != ig0.cutoff || ig.max_n_edge != ig0.max_n_edge)
            throw std::string("the device engine shares an interaction's parameters between systems, "
                    "but the systems' parameters differ");
        d.n_elem1.push_back(ig.n_elem1);
        for(int i=0; i<ig.n_elem1; ++i) {
            d.loc1.push_back(ig.loc1[i] + d.pos1->offset[ns]);
            d.type1.push_back(ig.types1[i]);
            d.id1.push_back(ig.id1[i]);
            d.sequence_id1.push_back(ig.sequence_id1[i]);
        }
        if(IType::symmetric) continue;
        d.n_elem2.push_back(ig.n_elem2);
        for(int i=0; i<ig.n_elem2; ++i) {
            d.loc2.push_back(ig.loc2[i] + d.pos2->offset[ns]);
            d.type2.push_back(ig.types2[i]);
            d.id2.push_back(ig.id2[i]);
            d.sequence_id2.push_back(ig.sequence_id2[i]);
        }
    }
    d.interaction_param = ig0.get_param();
    d.n_type1 = ig0.n_type1;
    d.n_type2 = ig0.n_type2;
    d.max_n_edge = ig0.max_n_edge;
    d.cutoff = ig0.cutoff;
    d.cache_buffer = ig0.pairlist.get_cache_buffer();
    return d;
}

// device_coord.cu
DeviceNode* device_distance(DeviceCoord& out, DeviceCoord& pos1, DeviceCoord& pos2,
        const std::vector<int>& atom1, const std::vector<int>& atom2);
DeviceNode* device_angle(DeviceCoord& out, DeviceCoord& pos, const std::vector<int>& atoms);  // 3 per element
// 4 per element; an element whose first or last atom is -1 is a dummy angle, its other atoms unused
DeviceNode* device_dihedral(DeviceCoord& out, DeviceCoord& pos, const std::vector<int>& atoms);
DeviceNode* device_paste_rama(DeviceCoord& out, DeviceCoord& pos1, DeviceCoord& pos2,
        const std::vector<int>& index1, const std::vector<int>& index2);
// element i copies n_dim components of element index[i] of pos2 if from2[i], else of pos1
DeviceNode* device_cat_pos(DeviceCoord& out, DeviceCoord& pos1, DeviceCoord& pos2, int n_dim,
        const std::vector<int>& index, const std::vector<int>& from2);

// device_spring.cu
DeviceNode* device_spring(AccumView potential, DeviceCoord& pos, int dim1, bool pbc, float box_len,
        const std::vector<int>& id, const std::vector<float>& equil_dist, const std::vector<float>& spring_const,
        const std::vector<int>& system);

// device_hbond.cu
// sites are the donors then the acceptors of each system; atom holds 3 per site
DeviceNode* device_infer_h_o(DeviceCoord& out, DeviceCoord& pos, const std::vector<int>& atom,
        const std::vector<float>& bond_length);
// n_donor of each system; the interaction's first set is the donors, its second the acceptors, both of `infer`
DeviceNode* device_protein_hbond(DeviceGraph& graph, DeviceCoord& out, DeviceCoord& infer,
        const DeviceInteractionData& ig, const std::vector<int>& n_donor);
// per residue: its rama element, and the protein_hbond elements whose occupancy it counts (donors then acceptors,
// in index order) and whose sensitivity takes its energy (its first donor and first acceptor, -1 for none); per
// peptide (residue i to i+1 of a system): its residue i, and the elements it counts; params as hbond_energy's
DeviceNode* device_hbond_energy(AccumView potential, DeviceCoord& protein_hbond, DeviceCoord& rama,
        const std::vector<float>& params, const std::vector<int>& residue_class,
        const std::vector<int>& residue_rama, const std::vector<int>& residue_first_donor,
        const std::vector<int>& residue_first_acceptor, const std::vector<int>& residue_count_start,
        const std::vector<int>& residue_count, const std::vector<int>& peptide_residue,
        const std::vector<int>& peptide_count_start, const std::vector<int>& peptide_count,
        const std::vector<int>& residue_system);

// device_rama.cu
// the spline coefficients of LayeredPeriodicSpline2D<1>, shared by every system
DeviceNode* device_rama_map_pot(AccumView potential, DeviceCoord& rama, const std::vector<int>& residue,
        const std::vector<int>& map_id, const std::vector<float>& coefficients, int nx, int ny,
        const std::vector<int>& system);

// device_affine.cu
// atom holds 3 per residue, ref_geom 9 (atom, xyz)
DeviceNode* device_affine_alignment(DeviceCoord& out, DeviceCoord& pos, const std::vector<int>& atom,
        const std::vector<float>& ref_geom);
// ref_pos holds 12 per residue (4 atoms, xyz); n_residue per system
DeviceNode* device_backbone_pairs(DeviceGraph& graph, AccumView potential, DeviceCoord& alignment,
        const std::vector<int>& residue, const std::vector<int>& n_atom, const std::vector<float>& ref_pos,
        const std::vector<int>& sequence_id, const std::vector<int>& n_residue, float dist_cutoff,
        float cache_buffer, int max_n_edge);

// device_placement.cu
// kinds: each output component group's PlaceT (0 scalar, 1 vector, 2 point), n_pos_dim values per element.  The
// values are table[layer*n_pos_dim + d] (fixed placement, rama null), or LayeredPeriodicSpline2D<n_pos_dim>'s
// coefficients table at rama element rama_residue (nx, ny its grid)
DeviceNode* device_placement(DeviceCoord& out, DeviceCoord& alignment, DeviceCoord* rama, const std::vector<int>& kinds,
        int n_pos_dim, const std::vector<int>& affine_residue, const std::vector<int>& layer,
        const std::vector<float>& table, int nx, int ny, const std::vector<int>& rama_residue);

// device_environment.cu
// coverage nodes (environment.cpp's environment_coverage, hb_environment_coverage, hbbb_coverage; hbond.cpp's
// hbond_coverage): the value of edge (i1, i2) of the interaction adds to output element base1[i1] + add2[i2]
DeviceNode* device_environment_coverage(DeviceGraph& graph, DeviceCoord& out, const DeviceInteractionData& ig,
        const std::vector<int>& base1, const std::vector<int>& add2);
DeviceNode* device_hbond_environment_coverage(DeviceGraph& graph, DeviceCoord& out, const DeviceInteractionData& ig,
        const std::vector<int>& base1, const std::vector<int>& add2);
DeviceNode* device_hbond_backbone_coverage(DeviceGraph& graph, DeviceCoord& out, const DeviceInteractionData& ig,
        const std::vector<int>& base1, const std::vector<int>& add2);
DeviceNode* device_hbond_coverage(DeviceGraph& graph, DeviceCoord& out, const DeviceInteractionData& ig,
        const std::vector<int>& base1, const std::vector<int>& add2);
DeviceNode* device_weighted_pos(DeviceCoord& out, DeviceCoord& pos, DeviceCoord& energy,
        const std::vector<int>& index_pos, const std::vector<int>& index_weight);
// per residue: its coupling type and first input element (n_restype in a row); scale, center, sharpness per type,
// weights n_restype per type
DeviceNode* device_sigmoid_coupling(AccumView potential, DeviceCoord& input, int n_restype,
        const std::vector<float>& scale, const std::vector<float>& center, const std::vector<float>& sharpness,
        const std::vector<float>& weights, const std::vector<int>& coupling_type,
        const std::vector<int>& residue_input, const std::vector<int>& residue_system);
// per site: its first side-chain coverage element (n_restype in a row) and its H-bond coverage element
DeviceNode* device_bb_sigmoid_coupling(AccumView potential, DeviceCoord& env_sc, DeviceCoord& env_hb, int n_restype,
        float scale, float center, float sharpness, float hbond_weight, const std::vector<float>& weights,
        const std::vector<int>& site_sc, const std::vector<int>& site_hb, const std::vector<int>& site_system);

// device_rotamer.cu
// rotamer.cpp's RotamerSidechain.  Per bead (the interaction's elements): its residue's rotamer count (1, 3 or 6),
// the residue's index among the residues of that count in its system (the id's residue field), and its rotamer;
// n_residue holds, per system, the residues of each count (3 per system: 1, 3, 6).  onebody are the nodes whose
// output is each bead's 1-body energy.
DeviceNode* device_rotamer(DeviceGraph& graph, AccumView potential, const DeviceInteractionData& ig,
        const std::vector<DeviceCoord*>& onebody, const std::vector<int>& bead_n_rot,
        const std::vector<int>& bead_residue, const std::vector<int>& bead_rot, const std::vector<int>& n_residue,
        float damping, int max_iter, float tol, int iteration_chunk_size);

// device_martini.cu
// martini_hybrid_position.  Per BB site: its BB and O atoms, the infer_H_O element of its O (-1 for none), its N, CA,
// C atoms (3 per site), the weights of N, CA, C, O that build BB and those that carry its sensitivity back (4 per
// site, 0 for a term left out), and the factor that renormalises BB (1 where the weights sum to 1); derived marks
// every site's BB and O.
DeviceNode* device_hybrid_position(DeviceCoord& out, DeviceCoord& source, DeviceCoord& hbond,
        const std::vector<unsigned char>& derived, const std::vector<int>& bb, const std::vector<int>& o,
        const std::vector<int>& o_elem, const std::vector<int>& src, const std::vector<float>& weight,
        const std::vector<float>& back_weight, const std::vector<float>& rescale);
// martini_potential, with the pair table every system shares: atom i's partners of higher index are entries
// [row_start[i], row_start[i+1]), each with its partner and its spline (an index into splines)
DeviceNode* device_martini_potential(AccumView potential, DeviceCoord& pos, const std::vector<int>& row_start,
        const std::vector<int>& neighbor, const std::vector<int>& spline,
        const std::vector<LayeredClampedSpline1D<1>>& splines, float r_min, float r_max, float cutoff,
        float cache_buffer, float box_x, float box_y, float box_z);
// martini_sc_table_1body, with the tables every system shares.  Per row: its CB element, its first layer (that of
// target 0; -1 for a row without a table), its group count and its system; per environment entry its atom and target,
// each system's entries starting at env_start (then their total).
DeviceNode* device_sc_table_1body(DeviceCoord& out, DeviceCoord& pos, DeviceCoord& cb,
        const std::vector<float>& spline_coeff, int n_grid, int n_theta, float grid_step_ang, float theta_step,
        float cutoff_ang, float cache_buffer, float box_x, float box_y, float box_z, const std::vector<int>& row_cb,
        const std::vector<int>& row_layer, const std::vector<int>& row_group, const std::vector<int>& row_system,
        const std::vector<int>& env_atom, const std::vector<int>& env_target, const std::vector<int>& env_start);
// dist_spring (2 atoms per term) and angle_spring (3), the terms the hybrid keeps, with the canonical spline of
// delta^2/2 every system shares
DeviceNode* device_dist_spring(AccumView potential, DeviceCoord& pos, const std::vector<int>& atom,
        const std::vector<float>& equil_dist, const std::vector<float>& spring_const, const std::vector<int>& system,
        const LayeredClampedSpline1D<1>& spline, float delta_min, float coord_scale,
        float box_x, float box_y, float box_z);
DeviceNode* device_angle_spring(AccumView potential, DeviceCoord& pos, const std::vector<int>& atom,
        const std::vector<float>& cos_theta0, const std::vector<float>& spring_const, const std::vector<int>& system,
        const LayeredClampedSpline1D<1>& spline, float delta_min, float coord_scale,
        float box_x, float box_y, float box_z);

#endif
