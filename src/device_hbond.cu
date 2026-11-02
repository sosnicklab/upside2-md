// Device implementations of hbond.cpp's infer_H_O, protein_hbond and hbond_energy
#include "device_interaction_graph.h"
#include "hbond_interaction.h"

using namespace std;

namespace {

// hbond_energy's energy per H-bond of a residue at rama r and its derivatives by phi and psi: the residue's helix,
// sheet and turn energies e_alpha, e_beta, e_other weighted by its Ramachandran region (hbond_energy's parameters 4
// to 11, its boundaries and sharpnesses).  The arithmetic of HBondEnergy::compute_value, written out again because
// moving the CPU's into a function changes the executable's code generation (findings 11g).
UPSIDE_HD inline void hbond_rama_energy(float& E, float& dE_dphi, float& dE_dpsi, const Vec<2>& r,
        float e_alpha, float e_beta, float e_other,
        float boundary_turn1,  float sharpness_turn1,  float boundary_turn2,  float sharpness_turn2,
        float boundary_helix1, float sharpness_helix1, float boundary_helix2, float sharpness_helix2)
{
    auto phiv1 = compact_sigmoid(r.v[0]-boundary_turn1, sharpness_turn1);
    auto phiv2 = compact_sigmoid(r.v[0]-boundary_turn2, sharpness_turn2);
    auto psiv1 = compact_sigmoid(r.v[1]-boundary_helix1, sharpness_helix1);
    auto psiv2 = compact_sigmoid(r.v[1]-boundary_helix2, sharpness_helix2);

    float turn_score   = (1.-phiv1.x())*phiv2.x();
    float psi_score    = (1.-psiv1.x())*psiv2.x();
    float n_turn_score = 1. - turn_score;
    float n_psi_score  = 1. - psi_score;
    float helix_score  = n_turn_score * psi_score;
    float sheet_score  = n_turn_score * n_psi_score;

    float dturn_dphi  = phiv2.y() - (phiv1.y()*phiv2.x() + phiv1.x()*phiv2.y());
    float dpsis_dpsi  = psiv2.y() - (psiv1.y()*psiv2.x() + psiv1.x()*psiv2.y());
    float dhelix_dphi = -dturn_dphi * psi_score;
    float dsheet_dphi = -dturn_dphi * n_psi_score;
    float dhelix_dpsi = n_turn_score *  dpsis_dpsi;
    float dsheet_dpsi = n_turn_score * -dpsis_dpsi;

    E       = e_alpha * helix_score + e_beta * sheet_score + e_other * turn_score;
    dE_dphi = e_alpha * dhelix_dphi + e_beta * dsheet_dphi + e_other * dturn_dphi;
    dE_dpsi = e_alpha * dhelix_dpsi + e_beta * dsheet_dpsi;
}

struct DeviceInferHO : public DeviceNode {
    DeviceCoord &out, &pos;
    DeviceArray<int> atom;          // 3 per site: previous, current and next backbone atom
    DeviceArray<float> bond_length;
    DeviceArray<float> data;        // unit vectors to prev and next, and their normalised sum, each with its 1/|v|

    DeviceInferHO(DeviceCoord& out_, DeviceCoord& pos_, const vector<int>& atom_, const vector<float>& bond_length_):
        out(out_), pos(pos_), atom(atom_), bond_length(bond_length_), data(12*bond_length_.size()) {}

    void compute_value(ComputeMode mode) override {
        VecArray x = pos.value();
        float* o = out.output.ptr;
        int w = out.row_width;
        const int* a = atom.ptr;
        const float* bl = bond_length.ptr;
        float* dd = data.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nv) {
            // For output, first three components are position of H/O and second three are HN/OC bond direction
            // Bond direction is a unit vector; curr is the N or C atom.  This algorithm assumes perfect 120 degree
            // bond angles.
            auto prev_c = load_vec<3>(x, a[3*nv+0]);
            auto curr_c = load_vec<3>(x, a[3*nv+1]);
            auto next_c = load_vec<3>(x, a[3*nv+2]);

            auto prev = prev_c - curr_c; auto prev_invmag = inv_mag(prev); prev *= prev_invmag;
            auto next = next_c - curr_c; auto next_invmag = inv_mag(next); next *= next_invmag;
            auto disp = prev   + next  ; auto disp_invmag = inv_mag(disp); disp *= disp_invmag;

            auto hbond_dir = -disp;
            auto hbond_pos = bl[nv]*hbond_dir + curr_c;

            store_vec(dd+12*nv+0, prev); dd[12*nv+ 3] = prev_invmag;
            store_vec(dd+12*nv+4, next); dd[12*nv+ 7] = next_invmag;
            store_vec(dd+12*nv+8, disp); dd[12*nv+11] = disp_invmag;

            store_vec(o+nv*w+0, hbond_pos);
            store_vec(o+nv*w+3, hbond_dir);
        });
    }

    void propagate_deriv() override {
        auto s = out.sens.view(), ps = pos.sens.view();
        int w = out.row_width, wp = pos.row_width;
        const int* a = atom.ptr;
        const float* bl = bond_length.ptr;
        const float* dd = data.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nv) {
            auto sens_pos = make_vec3(s.get(nv*w+0), s.get(nv*w+1), s.get(nv*w+2));
            auto sens_dir = make_vec3(s.get(nv*w+3), s.get(nv*w+4), s.get(nv*w+5));

            auto sens_neg_unitdisp = sens_dir + bl[nv]*sens_pos;

            auto prev = load_vec<3>(dd+12*nv+0); float prev_invmag = dd[12*nv+ 3];
            auto next = load_vec<3>(dd+12*nv+4); float next_invmag = dd[12*nv+ 7];
            auto disp = load_vec<3>(dd+12*nv+8); float disp_invmag = dd[12*nv+11];

            auto sens_nonunit_disp =   disp_invmag *(dot(disp,sens_neg_unitdisp)*disp - sens_neg_unitdisp);
            auto sens_nonunit_prev = (-prev_invmag)*(dot(prev,sens_nonunit_disp)*prev - sens_nonunit_disp);
            auto sens_nonunit_next = (-next_invmag)*(dot(next,sens_nonunit_disp)*next - sens_nonunit_disp);

            auto sens_curr = sens_pos - sens_nonunit_prev - sens_nonunit_next;
            for(int d=0; d<3; ++d) {
                ps.add(a[3*nv+0]*wp+d, sens_nonunit_prev[d]);
                ps.add(a[3*nv+1]*wp+d, sens_curr[d]);
                ps.add(a[3*nv+2]*wp+d, sens_nonunit_next[d]);
            }
        });
    }
};


struct DeviceProteinHBond : public DeviceNode {
    DeviceCoord &out, &infer;
    DeviceInteractionGraph<ProteinHBondInteraction> igraph;
    DeviceArray<int> elem1, elem2;   // the output element of each donor and acceptor of the interaction
    DeviceAccum hb_log;              // per output element: the sum of -log(1-hb) over its H-bonds
    DeviceArray<float> sens_scaled;

    DeviceProteinHBond(DeviceGraph& graph, DeviceCoord& out_, DeviceCoord& infer_, const DeviceInteractionData& ig,
            const vector<int>& n_donor):
        out(out_), infer(infer_), igraph(graph, "protein_hbond", ig), hb_log(out_.n_elem), sens_scaled(out_.n_elem)
    {
        if(out.offset != infer.offset)
            throw string("protein_hbond needs one row per infer_H_O row in every system");
        vector<int> e1, e2;
        for(int ns=0; ns<graph.n_system; ++ns) {
            if(ig.n_elem1[ns] != n_donor[ns] || ig.n_elem1[ns]+ig.n_elem2[ns] != out.offset[ns+1]-out.offset[ns])
                throw string("protein_hbond's interaction must hold its donors then its acceptors");
            for(int nd=0; nd<ig.n_elem1[ns]; ++nd) e1.push_back(out.offset[ns] + nd);
            for(int na=0; na<ig.n_elem2[ns]; ++na) e2.push_back(out.offset[ns] + n_donor[ns] + na);
        }
        elem1 = DeviceArray<int>(e1);
        elem2 = DeviceArray<int>(e2);
    }

    void compute_value(ComputeMode mode) override {
        VecArray ho = infer.value();
        float* vs = out.output.ptr;
        int w = out.row_width;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nv) {for(int d=0; d<6; ++d) vs[nv*w+d] = ho(d,nv);});

        // Compute protein hbonding score and its derivative
        igraph.compute_edges();
        hb_log.zero();
        auto hl = hb_log.view();
        const int *e1 = igraph.pairlist.edge1.ptr, *e2 = igraph.pairlist.edge2.ptr, *ne = igraph.pairlist.n_edge.ptr;
        const int *el1 = elem1.ptr, *el2 = elem2.ptr;
        const float* val = igraph.edge_value.ptr;
        int cap = igraph.max_n_edge;
        slots_for(igraph.n_system, cap, ne, DEVICE_LAMBDA (int e) {
            hl.add(el1[e1[e]], val[e]);
            hl.add(el2[e2[e]], val[e]);
        });
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nv) {vs[nv*w+6] = -expm1f(-hl.get(nv));});  // as the CPU's
    }

    void propagate_deriv() override {
        // we accumulated derivatives for z = 1-exp(-log(no_hb)), so we need to convert back with z_sens*(1.f-hb)
        auto s = out.sens.view();
        const float* vs = out.output.ptr;
        float* ss = sens_scaled.ptr;
        int w = out.row_width;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nv) {ss[nv] = s.get(nv*w+6) * (1.f-vs[nv*w+6]);});

        // Push protein HBond derivatives
        const int *e1 = igraph.pairlist.edge1.ptr, *e2 = igraph.pairlist.edge2.ptr, *ne = igraph.pairlist.n_edge.ptr;
        const int *el1 = elem1.ptr, *el2 = elem2.ptr;
        float* es = igraph.edge_sensitivity.ptr;
        int cap = igraph.max_n_edge;
        slots_for(igraph.n_system, cap, ne, DEVICE_LAMBDA (int e) {
            es[e] = ss[el1[e1[e]]] + ss[el2[e2[e]]];
        });
        igraph.propagate_derivatives();

        // pass through derivatives on all other components
        auto is = infer.sens.view();
        int wi = infer.row_width;
        const int *l1 = igraph.loc1.ptr, *l2 = igraph.loc2.ptr;
        parallel_for(igraph.n_elem1, DEVICE_LAMBDA (int i) {
            for(int d=0; d<6; ++d) is.add(l1[i]*wi+d, s.get(el1[i]*w+d));});
        parallel_for(igraph.n_elem2, DEVICE_LAMBDA (int i) {
            for(int d=0; d<6; ++d) is.add(l2[i]*wi+d, s.get(el2[i]*w+d));});
    }
};


struct DeviceHBondEnergy : public DeviceNode {
    AccumView potential;
    DeviceCoord &pp, &rama;
    int n_res, n_pep;
    DeviceArray<float> params;
    DeviceArray<int> residue_class, residue_rama, residue_first_donor, residue_first_acceptor;
    DeviceArray<int> residue_count_start, residue_count, peptide_residue, peptide_count_start, peptide_count;
    DeviceArray<int> residue_system;

    DeviceHBondEnergy(AccumView potential_, DeviceCoord& pp_, DeviceCoord& rama_, const vector<float>& params_,
            const vector<int>& residue_class_, const vector<int>& residue_rama_,
            const vector<int>& residue_first_donor_, const vector<int>& residue_first_acceptor_,
            const vector<int>& residue_count_start_, const vector<int>& residue_count_,
            const vector<int>& peptide_residue_, const vector<int>& peptide_count_start_,
            const vector<int>& peptide_count_, const vector<int>& residue_system_):
        potential(potential_), pp(pp_), rama(rama_), n_res(residue_rama_.size()), n_pep(peptide_residue_.size()),
        params(params_), residue_class(residue_class_), residue_rama(residue_rama_),
        residue_first_donor(residue_first_donor_), residue_first_acceptor(residue_first_acceptor_),
        residue_count_start(residue_count_start_), residue_count(residue_count_),
        peptide_residue(peptide_residue_), peptide_count_start(peptide_count_start_), peptide_count(peptide_count_),
        residue_system(residue_system_) {}

    void compute_value(ComputeMode mode) override {
        VecArray ramac = rama.value(), ppv = pp.value();
        auto rs = rama.sens.view(), ps = pp.sens.view(), pot = potential;
        bool want_pot = mode == PotentialAndDerivMode;
        int wr = rama.row_width, wp = pp.row_width;
        const float* p = params.ptr;
        const int *cls = residue_class.ptr, *rr = residue_rama.ptr, *fd = residue_first_donor.ptr;
        const int *fa = residue_first_acceptor.ptr, *rcs = residue_count_start.ptr, *rc = residue_count.ptr;
        const int *pr = peptide_residue.ptr, *pcs = peptide_count_start.ptr, *pc = peptide_count.ptr;
        const int *sys = residue_system.ptr;

        parallel_for(n_res, DEVICE_LAMBDA (int i) {
            // offset of the residue's class on basin energy k (0 alpha, 1 beta, 2 other)
            auto class_offset = [&](int k) {return cls[i] ? p[12 + 3*(cls[i]-1) + k] : 0.f;};
            float E, dPhi, dPsi;
            hbond_rama_energy(E, dPhi, dPsi, load_vec<2>(ramac, rr[i]),
                    p[0] + class_offset(0), p[1] + class_offset(1), p[2] + class_offset(2),
                    p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11]);

            // hbond number in the residue, its donors' occupancies then its acceptors'
            float hb_number1 = 0.f;
            for(int k=rcs[i]; k<rcs[i+1]; ++k) hb_number1 += ppv(6, rc[k]);
            if(want_pot) pot.add(sys[i], hb_number1*E);
            if(fd[i]>=0) ps.add(fd[i]*wp+6, E);
            if(fa[i]>=0) ps.add(fa[i]*wp+6, E);
            rs.add(rr[i]*wr+0, hb_number1*dPhi);
            rs.add(rr[i]*wr+1, hb_number1*dPsi);
        });

        parallel_for(n_pep, DEVICE_LAMBDA (int k) {
            float hb_number2 = 0.f;  // hbond number in the peptide
            for(int j=pcs[k]; j<pcs[k+1]; ++j) hb_number2 += ppv(6, pc[j]);
            if (hb_number2 <= 1.0) return;
            int i = pr[k];
            auto dE       = compact_sigmoid(1.5-hb_number2, 2.0);
            if(want_pot) pot.add(sys[i], dE.x() * p[3]);
            auto dE_dHB   = -dE.y() * p[3];
            if(fd[i+1]>=0) ps.add(fd[i+1]*wp+6, dE_dHB);
            if(fa[i  ]>=0) ps.add(fa[i  ]*wp+6, dE_dHB);
        });
    }
};

}

DeviceNode* device_infer_h_o(DeviceCoord& out, DeviceCoord& pos, const vector<int>& atom,
        const vector<float>& bond_length) {
    return new DeviceInferHO(out, pos, atom, bond_length);
}

DeviceNode* device_protein_hbond(DeviceGraph& graph, DeviceCoord& out, DeviceCoord& infer,
        const DeviceInteractionData& ig, const vector<int>& n_donor) {
    return new DeviceProteinHBond(graph, out, infer, ig, n_donor);
}

DeviceNode* device_hbond_energy(AccumView potential, DeviceCoord& protein_hbond, DeviceCoord& rama,
        const vector<float>& params, const vector<int>& residue_class,
        const vector<int>& residue_rama, const vector<int>& residue_first_donor,
        const vector<int>& residue_first_acceptor, const vector<int>& residue_count_start,
        const vector<int>& residue_count, const vector<int>& peptide_residue,
        const vector<int>& peptide_count_start, const vector<int>& peptide_count,
        const vector<int>& residue_system) {
    return new DeviceHBondEnergy(potential, protein_hbond, rama, params, residue_class, residue_rama,
            residue_first_donor, residue_first_acceptor, residue_count_start, residue_count, peptide_residue,
            peptide_count_start, peptide_count, residue_system);
}
