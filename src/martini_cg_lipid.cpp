#include "deriv_engine.h"
#include "interaction_graph.h"
#include "bead_interaction.h"
#include "state_logger.h"

#include <string>
#include <vector>

using namespace std;
using namespace h5;

// ===================================================================
// ComposeVector6D — Combines 3D positions with static 3D directions
// ===================================================================
struct ComposeVector6D : public CoordNode {
    CoordNode& pos;
    vector<int> elem_index;
    vector<float> direction;

    ComposeVector6D(hid_t grp, CoordNode& pos_node);

    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
};

ComposeVector6D::ComposeVector6D(hid_t grp, CoordNode& pos_node)
    : CoordNode(h5::get_dset_size(1, grp, "elem_index")[0], 6)
    , pos(pos_node)
{
    check_elem_width(pos_node, 3);

    int n = n_elem;
    elem_index.resize(n);
    traverse_dset<1, int32_t>(grp, "elem_index", [&](size_t i, int32_t v) {
        elem_index[i] = v;
    });

    direction.resize(n * 3);
    traverse_dset<2, float>(grp, "direction", [&](size_t i, size_t d, float v) {
        direction[i * 3 + d] = v;
    });

    fill(output, 0.f);
}

void ComposeVector6D::compute_value(ComputeMode) {
    VecArray posc = pos.output;
    for (int i = 0; i < n_elem; i++) {
        int ai = elem_index[i];
        output(0, i) = posc(0, ai);
        output(1, i) = posc(1, ai);
        output(2, i) = posc(2, ai);
        output(3, i) = direction[i * 3];
        output(4, i) = direction[i * 3 + 1];
        output(5, i) = direction[i * 3 + 2];
    }
}

void ComposeVector6D::propagate_deriv() {
    VecArray pos_sens = pos.sens;
    for (int i = 0; i < n_elem; i++) {
        int ai = elem_index[i];
        pos_sens(0, ai) += sens(0, i);
        pos_sens(1, ai) += sens(1, i);
        pos_sens(2, ai) += sens(2, i);
    }
    // Direction gradients (components 3-5) are zeroed — direction is static
}

// ===================================================================
// CGLipidPairPotential — CG lipid ↔ CG lipid via quadspline
// ===================================================================
struct CGLipidPairPotential : public PotentialNode {
    InteractionGraph<PosQuadSplineInteraction> igraph;

    CGLipidPairPotential(hid_t grp, CoordNode& cg_pos);
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
    virtual void set_param(const vector<float>& new_param) override {
        igraph.set_param(new_param);
    }
};

CGLipidPairPotential::CGLipidPairPotential(hid_t grp, CoordNode& cg_pos)
    : PotentialNode()
    , igraph(open_group(grp, "pair_interaction").get(), &cg_pos)
{}

void CGLipidPairPotential::compute_value(ComputeMode mode) {
    igraph.compute_edges();
    float total = 0.f;
    for (int ne = 0; ne < igraph.n_edge; ne++)
        total += igraph.edge_value[ne];
    potential = total;
    (void)mode;
}

void CGLipidPairPotential::propagate_deriv() {
    for (int ne = 0; ne < igraph.n_edge; ne++)
        igraph.edge_sensitivity[ne] = 1.f;
    igraph.propagate_derivatives();
}

// ===================================================================
// CGLipidSCPotential — CG lipid ↔ Sidechain via quadspline
// ===================================================================
struct CGLipidSCPotential : public PotentialNode {
    InteractionGraph<PosQuadSplineInteractionAsymmetric> igraph;

    CGLipidSCPotential(hid_t grp, CoordNode& sc_pos, CoordNode& cg_pos);
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
    virtual void set_param(const vector<float>& new_param) override {
        igraph.set_param(new_param);
    }
};

CGLipidSCPotential::CGLipidSCPotential(hid_t grp, CoordNode& sc_pos, CoordNode& cg_pos)
    : PotentialNode()
    , igraph(open_group(grp, "pair_interaction").get(), &sc_pos, &cg_pos)
{}

void CGLipidSCPotential::compute_value(ComputeMode mode) {
    igraph.compute_edges();
    float total = 0.f;
    for (int ne = 0; ne < igraph.n_edge; ne++)
        total += igraph.edge_value[ne];
    potential = total;
    (void)mode;
}

void CGLipidSCPotential::propagate_deriv() {
    for (int ne = 0; ne < igraph.n_edge; ne++)
        igraph.edge_sensitivity[ne] = 1.f;
    igraph.propagate_derivatives();
}

// ===================================================================
// Node registration
// ===================================================================
static RegisterNodeType<ComposeVector6D, 1> _reg_cv("compose_vector6d");
static RegisterNodeType<CGLipidPairPotential, 1> _reg_cg_pair("cg_lipid_pair");
static RegisterNodeType<CGLipidSCPotential, 2> _reg_cg_sc("cg_lipid_sc");
