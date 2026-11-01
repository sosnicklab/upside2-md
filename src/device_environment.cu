// Device implementations of environment.cpp's coverage, weighted_pos and coupling nodes, and hbond.cpp's
// hbond_coverage
#include "device_interaction_graph.h"
#include "environment_interaction.h"
#include "hbond_interaction.h"

using namespace std;

namespace {

// A coverage node: each output element sums the values of its edges, base1[i1] + add2[i2] naming an edge's element
template <typename IType>
struct DeviceCoverage : public DeviceNode {
    DeviceCoord& out;
    DeviceInteractionGraph<IType> igraph;
    DeviceArray<int> base1, add2;
    DeviceAccum sum;

    DeviceCoverage(DeviceGraph& graph, const string& name, DeviceCoord& out_, const DeviceInteractionData& ig,
            const vector<int>& base1_, const vector<int>& add2_):
        out(out_), igraph(graph, name, ig), base1(base1_), add2(add2_), sum(out_.n_elem) {}

    void compute_value(ComputeMode mode) override {
        igraph.compute_edges();
        sum.zero();
        auto sv = sum.view();
        const int *e1 = igraph.pairlist.edge1.ptr, *e2 = igraph.pairlist.edge2.ptr, *ne = igraph.pairlist.n_edge.ptr;
        const int *b1 = base1.ptr, *a2 = add2.ptr;
        const float* val = igraph.edge_value.ptr;
        int cap = igraph.max_n_edge;
        slots_for(igraph.n_system, cap, ne, DEVICE_LAMBDA (int e) {
            sv.add(b1[e1[e]] + a2[e2[e]], val[e]);
        });
        float* o = out.output.ptr;
        int w = out.row_width;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int i) {o[i*w] = sv.get(i);});
    }

    void propagate_deriv() override {
        auto s = out.sens.view();
        int w = out.row_width;
        const int *e1 = igraph.pairlist.edge1.ptr, *e2 = igraph.pairlist.edge2.ptr, *ne = igraph.pairlist.n_edge.ptr;
        const int *b1 = base1.ptr, *a2 = add2.ptr;
        float* es = igraph.edge_sensitivity.ptr;
        int cap = igraph.max_n_edge;
        slots_for(igraph.n_system, cap, ne, DEVICE_LAMBDA (int e) {
            es[e] = s.get((b1[e1[e]] + a2[e2[e]])*w);
        });
        igraph.propagate_derivatives();
    }
};


struct DeviceWeightedPos : public DeviceNode {
    DeviceCoord &out, &pos, &energy;
    DeviceArray<int> index_pos, index_weight;

    DeviceWeightedPos(DeviceCoord& out_, DeviceCoord& pos_, DeviceCoord& energy_, const vector<int>& index_pos_,
            const vector<int>& index_weight_):
        out(out_), pos(pos_), energy(energy_), index_pos(index_pos_), index_weight(index_weight_) {}

    void compute_value(ComputeMode mode) override {
        VecArray x = pos.value(), en = energy.value();
        float* o = out.output.ptr;
        int w = out.row_width;
        const int *ip = index_pos.ptr, *iw = index_weight.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int ne) {
            o[ne*w+0] = x(0,ip[ne]);
            o[ne*w+1] = x(1,ip[ne]);
            o[ne*w+2] = x(2,ip[ne]);
            o[ne*w+3] = expf(-en(0,iw[ne]));
        });
    }

    void propagate_deriv() override {
        auto s = out.sens.view(), ps = pos.sens.view(), es = energy.sens.view();
        const float* o = out.output.ptr;
        int w = out.row_width, wp = pos.row_width, we = energy.row_width;
        const int *ip = index_pos.ptr, *iw = index_weight.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int ne) {
            for(int d=0; d<3; ++d) ps.add(ip[ne]*wp+d, s.get(ne*w+d));
            es.add(iw[ne]*we, -(o[ne*w+3]*s.get(ne*w+3)));  // exponential derivative
        });
    }
};


struct DeviceSigmoidCoupling : public DeviceNode {
    AccumView potential;
    DeviceCoord& input;
    int n_restype;
    DeviceArray<float> scale, center, sharpness, weights;
    DeviceArray<int> coupling_type, residue_input, residue_system;

    DeviceSigmoidCoupling(AccumView potential_, DeviceCoord& input_, int n_restype_, const vector<float>& scale_,
            const vector<float>& center_, const vector<float>& sharpness_, const vector<float>& weights_,
            const vector<int>& coupling_type_, const vector<int>& residue_input_, const vector<int>& residue_system_):
        potential(potential_), input(input_), n_restype(n_restype_), scale(scale_), center(center_),
        sharpness(sharpness_), weights(weights_), coupling_type(coupling_type_), residue_input(residue_input_),
        residue_system(residue_system_) {}

    void compute_value(ComputeMode mode) override {
        VecArray in = input.value();
        auto is = input.sens.view(), pot = potential;
        bool want_pot = mode == PotentialAndDerivMode;
        int nt = n_restype, wi = input.row_width;
        const float *sc = scale.ptr, *ce = center.ptr, *sh = sharpness.ptr, *wt = weights.ptr;
        const int *ct = coupling_type.ptr, *ri = residue_input.ptr, *sys = residue_system.ptr;
        parallel_for(int(coupling_type.n), DEVICE_LAMBDA (int nr) {
            float wnumber = 0.0;
            int ctype = ct[nr];
            for(int aa=0; aa<nt; ++aa)
                wnumber += wt[ctype*nt+aa] * in(0, ri[nr]+aa);
            auto out = compact_sigmoid(wnumber-ce[ctype], sh[ctype]);

            if(want_pot) pot.add(sys[nr], sc[ctype] * out.x());
            for(int aa=0; aa<nt; ++aa)
                is.add((ri[nr]+aa)*wi, wt[ctype*nt+aa] * sc[ctype] * out.y());
        });
    }
};


struct DeviceBackboneSigmoidCoupling : public DeviceNode {
    AccumView potential;
    DeviceCoord &env_sc, &env_hb;
    int n_restype;
    float scale, center, sharpness, hbond_weight;
    DeviceArray<float> weights;
    DeviceArray<int> site_sc, site_hb, site_system;

    DeviceBackboneSigmoidCoupling(AccumView potential_, DeviceCoord& env_sc_, DeviceCoord& env_hb_, int n_restype_,
            float scale_, float center_, float sharpness_, float hbond_weight_, const vector<float>& weights_,
            const vector<int>& site_sc_, const vector<int>& site_hb_, const vector<int>& site_system_):
        potential(potential_), env_sc(env_sc_), env_hb(env_hb_), n_restype(n_restype_), scale(scale_),
        center(center_), sharpness(sharpness_), hbond_weight(hbond_weight_), weights(weights_),
        site_sc(site_sc_), site_hb(site_hb_), site_system(site_system_) {}

    void compute_value(ComputeMode mode) override {
        VecArray sc = env_sc.value(), hb = env_hb.value();
        auto ss = env_sc.sens.view(), hs = env_hb.sens.view(), pot = potential;
        bool want_pot = mode == PotentialAndDerivMode;
        int nt = n_restype, wsc = env_sc.row_width, whb = env_hb.row_width;
        float sca = scale, cen = center, sha = sharpness, hbw = hbond_weight;
        const float* wt = weights.ptr;
        const int *isc = site_sc.ptr, *ihb = site_hb.ptr, *sys = site_system.ptr;
        parallel_for(int(site_sc.n), DEVICE_LAMBDA (int nr) {
            float bl = 0.0;
            for(int aa=0; aa<nt; ++aa)
                bl += wt[aa] * sc(0, isc[nr]+aa);
            bl += hbw*hb(0, ihb[nr]);
            auto out = compact_sigmoid(bl-cen, sha);

            if(want_pot) pot.add(sys[nr], sca * out.x());
            auto d_pot  = sca * out.y();
            for(int aa=0; aa<nt; ++aa)
                ss.add((isc[nr]+aa)*wsc, wt[aa] * d_pot);
            hs.add(ihb[nr]*whb, hbw* d_pot);
        });
    }
};

}

DeviceNode* device_environment_coverage(DeviceGraph& graph, DeviceCoord& out, const DeviceInteractionData& ig,
        const vector<int>& base1, const vector<int>& add2) {
    return new DeviceCoverage<EnvironmentCoverageInteraction>(graph, "environment_coverage", out, ig, base1, add2);
}
DeviceNode* device_hbond_environment_coverage(DeviceGraph& graph, DeviceCoord& out, const DeviceInteractionData& ig,
        const vector<int>& base1, const vector<int>& add2) {
    return new DeviceCoverage<HbondEnvironmentCoverageInteraction>(graph, "hb_environment_coverage", out, ig,
            base1, add2);
}
DeviceNode* device_hbond_backbone_coverage(DeviceGraph& graph, DeviceCoord& out, const DeviceInteractionData& ig,
        const vector<int>& base1, const vector<int>& add2) {
    return new DeviceCoverage<HbondEnvironmentCoverageInteraction2>(graph, "hbbb_coverage", out, ig, base1, add2);
}
DeviceNode* device_hbond_coverage(DeviceGraph& graph, DeviceCoord& out, const DeviceInteractionData& ig,
        const vector<int>& base1, const vector<int>& add2) {
    return new DeviceCoverage<HBondCoverageInteraction>(graph, "hbond_coverage", out, ig, base1, add2);
}
DeviceNode* device_weighted_pos(DeviceCoord& out, DeviceCoord& pos, DeviceCoord& energy,
        const vector<int>& index_pos, const vector<int>& index_weight) {
    return new DeviceWeightedPos(out, pos, energy, index_pos, index_weight);
}
DeviceNode* device_sigmoid_coupling(AccumView potential, DeviceCoord& input, int n_restype,
        const vector<float>& scale, const vector<float>& center, const vector<float>& sharpness,
        const vector<float>& weights, const vector<int>& coupling_type,
        const vector<int>& residue_input, const vector<int>& residue_system) {
    return new DeviceSigmoidCoupling(potential, input, n_restype, scale, center, sharpness, weights, coupling_type,
            residue_input, residue_system);
}
DeviceNode* device_bb_sigmoid_coupling(AccumView potential, DeviceCoord& env_sc, DeviceCoord& env_hb, int n_restype,
        float scale, float center, float sharpness, float hbond_weight, const vector<float>& weights,
        const vector<int>& site_sc, const vector<int>& site_hb, const vector<int>& site_system) {
    return new DeviceBackboneSigmoidCoupling(potential, env_sc, env_hb, n_restype, scale, center, sharpness,
            hbond_weight, weights, site_sc, site_hb, site_system);
}
