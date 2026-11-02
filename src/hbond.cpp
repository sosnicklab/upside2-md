#include "deriv_engine.h"
#include <string>
#include "timing.h"
#include <math.h>
#include <map>
#include <algorithm>
#include "state_logger.h"
#include "interaction_graph.h"
#include "spline.h"
#include "bead_interaction.h"
#include "hbond_interaction.h"
#include "device_nodes.h"
#include "vector_math.h"

#define N_KNOT_MID_NC 7
#define KNOT_SPACING_MID_NC 0.2
#define OFFSET_MID_NC 3.5

using namespace h5;
using namespace std;

struct Infer_H_O : public CoordNode
{
    struct Params {
        index_t atom[3];
        float bond_length;
    };

    CoordNode& pos;
    int n_donor, n_acceptor, n_virtual;
    vector<Params> params;
    unique_ptr<float[]> data_for_deriv;

    Infer_H_O(hid_t grp, CoordNode& pos_):
        CoordNode(
                get_dset_size(2, grp, "donors/id")[0]+get_dset_size(2, grp, "acceptors/id")[0], 6),
        pos(pos_),
        n_donor(get_dset_size(2, grp, "donors/id")[0]),
        n_acceptor(get_dset_size(2, grp, "acceptors/id")[0]),
        n_virtual(n_donor+n_acceptor), params(n_virtual),
        data_for_deriv(new_aligned<float>(n_virtual*3*4, 4))
    {
        auto don = open_group(grp, "donors");
        auto acc = open_group(grp, "acceptors");

        check_size(don.get(), "id",          n_donor,    3);
        check_size(don.get(), "bond_length", n_donor);
        check_size(acc.get(), "id",          n_acceptor, 3);
        check_size(acc.get(), "bond_length", n_acceptor);

        traverse_dset<2,int  >(don.get(),"id",          [&](size_t i,size_t j, int   x){params[        i].atom[j]=x;});
        traverse_dset<1,float>(don.get(),"bond_length", [&](size_t i,          float x){params[        i].bond_length  =x;});
        traverse_dset<2,int  >(acc.get(),"id",          [&](size_t i,size_t j, int   x){params[n_donor+i].atom[j]=x;});
        traverse_dset<1,float>(acc.get(),"bond_length", [&](size_t i,          float x){params[n_donor+i].bond_length  =x;});

        if(logging(LOG_EXTENSIVE)) {
            default_logger->add_logger<float>("virtual", {n_elem, 3}, [&](float* buffer) {
                    for(int nv=0; nv<n_virtual; ++nv) {
                        auto x = load_vec<6>(output, nv);
                        for(int d=0; d<3; ++d) buffer[nv*3 + d] = x[d];
                    }
                });
        }
    }


    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("infer_H_O"));

        VecArray posc  = pos.output;
        for(int nv=0; nv<n_virtual; ++nv) {
            // For output, first three components are position of H/O and second three are HN/OC bond direction
            // Bond direction is a unit vector
            // curr is the N or C atom
            // This algorithm assumes perfect 120 degree bond angles

            auto& p = params[nv];

            auto prev_c = Float4(&posc(0,p.atom[0]));
            auto curr_c = Float4(&posc(0,p.atom[1]));
            auto next_c = Float4(&posc(0,p.atom[2]));

            auto prev = prev_c - curr_c; auto prev_invmag = inv_mag(prev); prev *= prev_invmag;
            auto next = next_c - curr_c; auto next_invmag = inv_mag(next); next *= next_invmag;
            auto disp = prev   + next  ; auto disp_invmag = inv_mag(disp); disp *= disp_invmag;

            auto hbond_dir = -disp;
            auto hbond_pos = fmadd(Float4(p.bond_length),hbond_dir, curr_c);

            // store derived values for derivatives later
            prev.blend<0,0,0,1>(prev_invmag).store(data_for_deriv + nv*3*4 + 0);
            next.blend<0,0,0,1>(next_invmag).store(data_for_deriv + nv*3*4 + 4);
            disp.blend<0,0,0,1>(disp_invmag).store(data_for_deriv + nv*3*4 + 8);

            // write pos
            hbond_pos.store(&output(0,nv));
            hbond_dir.store(&output(3,nv), Alignment::unaligned);
        }
    }

    virtual void propagate_deriv() override {
        Timer timer(string("infer_H_O_deriv"));
        VecArray pos_sens = pos.sens;

        for(int nv=0; nv<n_virtual; ++nv) {
            const auto& p = params[nv];

            auto sens_pos = Float4(&sens(0,nv)).zero_entries<0,0,0,1>(); // last entry should be zero
            auto sens_dir = Float4(&sens(3,nv), Alignment::unaligned);

            auto sens_neg_unitdisp = sens_dir + Float4(p.bond_length)*sens_pos;

            // loading: first 3 entries are unitvec and last is inv_mag
            auto prev4 = Float4(data_for_deriv + nv*3*4 + 0);  auto prev_invmag = prev4.broadcast<3>();
            auto next4 = Float4(data_for_deriv + nv*3*4 + 4);  auto next_invmag = next4.broadcast<3>();
            auto disp4 = Float4(data_for_deriv + nv*3*4 + 8);  auto disp_invmag = disp4.broadcast<3>();

            // use dot3 here so we don't have to zero the last component of disp4, etc
            auto sens_nonunit_disp =   disp_invmag *fmsub(dot3(disp4,sens_neg_unitdisp),disp4, sens_neg_unitdisp);
            auto sens_nonunit_prev = (-prev_invmag)*fmsub(dot3(prev4,sens_nonunit_disp),prev4, sens_nonunit_disp);
            auto sens_nonunit_next = (-next_invmag)*fmsub(dot3(next4,sens_nonunit_disp),next4, sens_nonunit_disp);

            sens_nonunit_prev                                 .update(&pos_sens(0,p.atom[0]));
            (sens_pos - sens_nonunit_prev - sens_nonunit_next).update(&pos_sens(0,p.atom[1]));
            sens_nonunit_next                                 .update(&pos_sens(0,p.atom[2]));
        }
    }

    virtual DeviceNode* make_device(DeviceGraph& g, const vector<DerivComputation*>& systems) override {
        auto param = [&](int ns) -> const vector<Params>& {return static_cast<Infer_H_O&>(*systems[ns]).params;};
        return device_infer_h_o(g.coord(*this), g.coord(pos),
                g.offset_index(pos, [&](int ns) {
                    vector<int> r; for(auto& p: param(ns)) r.insert(r.end(), p.atom, p.atom+3); return r;}),
                g.concat<float>([&](int ns) {vector<float> r; for(auto& p: param(ns)) r.push_back(p.bond_length); return r;}));
    }
};
static RegisterNodeType<Infer_H_O,1> infer_node("infer_H_O");


struct ProteinHBond : public CoordNode
{
    CoordNode& infer;
    InteractionGraph<ProteinHBondInteraction> igraph;
    int n_donor, n_acceptor, n_virtual;
    unique_ptr<float[]> sens_scaled;

    ProteinHBond(hid_t grp, CoordNode& infer_):
        CoordNode(get_dset_size(1,grp,"index1")[0]+get_dset_size(1,grp,"index2")[0], 7),
        infer(infer_),
        igraph(grp, &infer, &infer) ,
        n_donor   (igraph.n_elem1),
        n_acceptor(igraph.n_elem2),
        n_virtual (n_donor+n_acceptor),
        sens_scaled(new_aligned<float>(n_virtual,4))
    {
        if(logging(LOG_DETAILED)) {
            default_logger->add_logger<float>("hbond", {n_donor+n_acceptor}, [&](float* buffer) {
                   for(int nv: range(n_donor+n_acceptor))
                       buffer[nv] = output(6,nv);});
        }
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("protein_hbond"));

        int n_virtual = n_donor + n_acceptor;
        VecArray vs = output;
        VecArray ho = infer.output;

        for(int nv: range(n_virtual)) {
            Float4(&ho(0,nv)).store(&vs(0,nv));
            Float4(&ho(4,nv)).store(&vs(4,nv)); // result is already zero padded
        }

        // Compute protein hbonding score and its derivative
        igraph.compute_edges();
        for(int ne=0; ne<igraph.n_edge; ++ne) {
            int nd = igraph.edge_indices1[ne];
            int na = igraph.edge_indices2[ne];
            float hb_log = igraph.edge_value[ne];
            vs(6,nd)         += hb_log;
            vs(6,na+n_donor) += hb_log;
        }

        // 1-exp(-x) by expm1f, since the subtraction keeps only whole steps of 2^-24 of a weak H-bond's score
        for(int nv: range(n_virtual)) vs(6,nv) = -expm1f(-vs(6,nv));
    }

    virtual void propagate_deriv() override {
        Timer timer(string("protein_hbond_deriv"));

        // we accumulated derivatives for z = 1-exp(-log(no_hb))
        // so we need to convert back with z_sens*(1.f-hb)
        for(int nv=0; nv<n_virtual; ++nv)
            sens_scaled[nv] = sens(6,nv) * (1.f-output(6,nv));

        // Push protein HBond derivatives
        for(int ne: range(igraph.n_edge)) {
            auto don_sens = sens_scaled[igraph.edge_indices1[ne]];
            auto acc_sens = sens_scaled[igraph.edge_indices2[ne]+n_donor];

            igraph.edge_sensitivity[ne] = don_sens + acc_sens;
        }
        igraph.propagate_derivatives();

        // pass through derivatives on all other components
        VecArray pd1 = igraph.pos_node1->sens;
        for(int nd=0; nd<n_donor; ++nd) {
            // the last component is taken care of by the edge loop
            update_vec(pd1, igraph.loc1[nd], load_vec<6>(sens, nd));
        }
        VecArray pd2 = igraph.pos_node2->sens;
        for(int na=0; na<n_acceptor; ++na) {  // acceptor loop
            // the last component is taken care of by the edge loop
            update_vec(pd2, igraph.loc2[na], load_vec<6>(sens, na+n_donor));
        }
    }

    virtual DeviceNode* make_device(DeviceGraph& g, const vector<DerivComputation*>& systems) override {
        vector<const InteractionGraph<ProteinHBondInteraction>*> igs;
        vector<int> n_donors;
        for(auto s: systems) {
            igs.push_back(&static_cast<ProteinHBond&>(*s).igraph);
            n_donors.push_back(static_cast<ProteinHBond&>(*s).n_donor);
        }
        return device_protein_hbond(g, g.coord(*this), g.coord(infer), device_interaction_data(g, igs), n_donors);
    }
};
static RegisterNodeType<ProteinHBond,1> hbond_node("protein_hbond");


struct HBondCoverage : public CoordNode {
    InteractionGraph<HBondCoverageInteraction> igraph;
    int n_sc;

    HBondCoverage(hid_t grp, CoordNode& infer_, CoordNode& sidechains_) try :
        CoordNode(get_dset_size(1,grp,"index2")[0], 1),
        igraph(grp, &infer_, &sidechains_),
        n_sc(igraph.n_elem2) {
    } catch (runtime_error& e) {
        cerr << e.what() << endl;
    } 

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("hbond_coverage"));

        // Compute coverage and its derivative
        igraph.compute_edges();

        fill(output, 0.f);
        for(int ne=0; ne<igraph.n_edge; ++ne) {
            output(0, igraph.edge_indices2[ne]) += igraph.edge_value[ne];
        }
    }

    virtual void propagate_deriv() override {
        Timer timer(string("hbond_coverage_deriv"));

        for(int ne: range(igraph.n_edge))
            igraph.edge_sensitivity[ne] = sens(0,igraph.edge_indices2[ne]);
        igraph.propagate_derivatives();
    }

    virtual std::vector<float> get_param() const override {return igraph.get_param();}
#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {return igraph.get_param_deriv();}
#endif
    virtual void set_param(const std::vector<float>& new_param) override {igraph.set_param(new_param);}

    virtual vector<float> get_value_by_name(const char* log_name) override {
        if(!strcmp(log_name, "count_edges_by_type")) {
            return igraph.count_edges_by_type();
        } else {
            throw string("Value ") + log_name + string(" not implemented");
        }
    }

    virtual DeviceNode* make_device(DeviceGraph& g, const vector<DerivComputation*>& systems) override {
        auto& out = g.coord(*this);
        vector<const InteractionGraph<HBondCoverageInteraction>*> igs;
        vector<int> base1, add2;  // an edge adds to its second element's output
        for(int ns: range(g.n_system)) {
            auto& n = static_cast<HBondCoverage&>(*systems[ns]);
            igs.push_back(&n.igraph);
            base1.insert(base1.end(), n.igraph.n_elem1, 0);
            for(int i: range(n.igraph.n_elem2)) add2.push_back(out.offset[ns] + i);
        }
        return device_hbond_coverage(g, out, device_interaction_data(g, igs), base1, add2);
    }
};
static RegisterNodeType<HBondCoverage,2> coverage_node("hbond_coverage");


// parameters holds the 12 shared values, then optionally 3 per extra residue class: offsets on
// E_alpha, E_beta and E_other for the residues that residue_class assigns to it (class 0 uses the
// shared energies alone, so a 12-value config has no residue_class and behaves as before).
struct HBondEnergy : public HBondCounter
{
    CoordNode& protein_hbond;
    CoordNode& rama;
    int n_res;
    int n_virtual;
    int n_donor;
    int n_acceptor;
    int n_param;
    std::vector<float> params;
    std::vector<int>   res_class;
    std::vector<float> hb_number1; // in one residue
    std::vector<float> hb_number2; // in one peptide
    std::vector<int>   donor_resid;
    std::vector<int>   acceptor_resid;
    std::map<int,int>  resid_donor;
    std::map<int,int>  resid_acceptor;

    std::vector<float> Ehbond;
    std::vector<float> dPhi;
    std::vector<float> dPsi;
    float E_beta;
    float E_alpha;
    float E_other;
    float boundary_turn1;
    float sharpness_turn1;
    float boundary_turn2;
    float sharpness_turn2;
    float boundary_helix1;
    float sharpness_helix1;
    float boundary_helix2;
    float sharpness_helix2;
    float E_bias;

    HBondEnergy(hid_t grp, CoordNode& protein_hbond_, CoordNode& rama_):
        HBondCounter(),
        protein_hbond(protein_hbond_),
        rama(rama_),
        n_res(get_dset_size(1, grp, "rama_resid")[0]),
        n_virtual(protein_hbond.n_elem),
        n_donor(get_dset_size(1, grp, "donor_resid")[0]),
        n_acceptor(get_dset_size(1, grp, "acceptor_resid")[0]),
        n_param(get_dset_size(1, grp, "parameters")[0]),
        params(n_param),
        res_class(n_res, 0),
        hb_number1(n_res),
        hb_number2(n_res-1),
        donor_resid(n_donor),
        acceptor_resid(n_acceptor),
        Ehbond(n_res),
        dPhi(n_res),
        dPsi(n_res)
    {
        if(n_donor+n_acceptor != n_virtual)
            throw string("hbond_energy donors and acceptors must be the rows of protein_hbond");
        if(n_res != rama.n_elem)
            throw string("hbond_energy rama_resid must have one row per rama_coord row");
        traverse_dset<1, int>(grp, "rama_resid", [&](size_t i, int x) {
                if(x != int(i)) throw string("hbond_energy rama_resid must be the rama_coord rows in order");});
        if(n_param < 12 || (n_param-12)%3)
            throw string("hbond_energy parameters must hold 12 values plus 3 per residue class");
        traverse_dset<1, float> (grp, "parameters",     [&](size_t i, float x) {params[i] = x;});
        if(n_param > 12) {
            int n_class = 1 + (n_param-12)/3;
            check_size(grp, "residue_class", n_res);
            traverse_dset<1, int>(grp, "residue_class", [&](size_t i, int x) {
                    if(x<0 || x>=n_class) throw string("hbond_energy residue_class out of range");
                    res_class[i] = x;});
        }
        // a donor's N-H closes the peptide before its residue, an acceptor's C=O the one after it
        traverse_dset<1, int>   (grp, "donor_resid",    [&](size_t i, int x)   {
                if(x<1 || x>=n_res) throw string("hbond_energy donor_resid outside 1..n_res-1");
                donor_resid[i] = x;    resid_donor.insert(pair<int, int>(x, (int)i));  });
        traverse_dset<1, int>   (grp, "acceptor_resid", [&](size_t i, int x)   {
                if(x<0 || x>=n_res-1) throw string("hbond_energy acceptor_resid outside 0..n_res-2");
                acceptor_resid[i] = x; resid_acceptor.insert(pair<int, int>(x, (int)i));});
        E_alpha         = params[0];
        E_beta          = params[1];
        E_other         = params[2];
        E_bias          = params[3];
        boundary_turn1  = params[4];
        sharpness_turn1 = params[5];
        boundary_turn2  = params[6];
        sharpness_turn2 = params[7];
        boundary_helix1  = params[8];
        sharpness_helix1 = params[9];
        boundary_helix2  = params[10];
        sharpness_helix2 = params[11];
    }

    // offset of residue i's class on basin energy k (0 alpha, 1 beta, 2 other)
    float class_offset(int i, int k) const {
        int c = res_class[i];
        return c ? params[12 + 3*(c-1) + k] : 0.f;
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("hbond_energy"));
        VecArray pp        = protein_hbond.output;
        VecArray pp_sens   = protein_hbond.sens;
        VecArray ramac     = rama.output;
        VecArray rama_sens = rama.sens;

        for( int i=0; i<n_res; i++) {
            auto r = load_vec<2>(ramac, i);

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

            float e_alpha = E_alpha + class_offset(i,0);
            float e_beta  = E_beta  + class_offset(i,1);
            float e_other = E_other + class_offset(i,2);

            Ehbond[i] = e_alpha * helix_score + e_beta * sheet_score + e_other * turn_score;
            dPhi[i]   = e_alpha * dhelix_dphi + e_beta * dsheet_dphi + e_other * dturn_dphi;
            dPsi[i]   = e_alpha * dhelix_dpsi + e_beta * dsheet_dpsi;

        }

        // count hbond number in one residue
        for( int i =0; i<n_res; i++ ) 
            hb_number1[i] = 0.0;
        for(int nd=0; nd<n_donor; ++nd) {
            int resid = donor_resid[nd];
            hb_number1[resid] += pp(6,nd);
        }
        for(int na=0; na<n_acceptor; ++na) {
            int resid = acceptor_resid[na];
            hb_number1[resid] += pp(6,na+n_donor);
        }

        // count hbond number in one peptide
        for( int i =0; i<n_res-1; i++ ) 
            hb_number2[i] = 0.0;
        for(int nd=0; nd<n_donor; ++nd) {
            int pepid = donor_resid[nd]-1;
            hb_number2[pepid] += pp(6,nd);
        }
        for(int na=0; na<n_acceptor; ++na) {
            int pepid = acceptor_resid[na];
            hb_number2[pepid] += pp(6,na+n_donor);
        }

        // the energy is linear in each residue's count, so its derivative by an occupancy is the residue's energy
        // per H-bond at any count, zero included
        n_hbond = 0.0;
        potential = 0.0;
        for( int i =0; i<n_res; i++ ) {
            n_hbond      += 0.5f*hb_number1[i];   // a bond counts at its donor's and its acceptor's residue
            potential    += hb_number1[i]*Ehbond[i];

            auto iter = resid_donor.find(i);
            if(iter!=resid_donor.end()) {
                pp_sens(6,iter->second) += Ehbond[i];
            }
            iter = resid_acceptor.find(i);
            if(iter!=resid_acceptor.end()) {
                pp_sens(6,iter->second+n_donor) += Ehbond[i];
            }

            rama_sens(0, i) += hb_number1[i]*dPhi[i];
            rama_sens(1, i) += hb_number1[i]*dPsi[i];
        }

        for( int i =0; i<n_res-1; i++ ) {
            if (hb_number2[i] <= 1.0)
                continue;
            auto dE       = compact_sigmoid(1.5-hb_number2[i], 2.0);
            potential    += dE.x() * E_bias;
            auto dE_dHB   = -dE.y() * E_bias ;
            auto iter = resid_donor.find(i+1);
            if(iter!=resid_donor.end()) {
                pp_sens(6,iter->second) += dE_dHB;
            }
            iter = resid_acceptor.find(i);
            if(iter!=resid_acceptor.end()) {
                pp_sens(6,iter->second+n_donor) += dE_dHB;
            }
        }
    }

    virtual std::vector<float> get_param() const override {return params;}
#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        std::vector<float> dparams(n_param, 0.f);
        VecArray ramac = rama.output;
        VecArray pp    = protein_hbond.output;

        for( int i =0; i<n_res; i++ ) 
            hb_number1[i] = 0.0;
        for(int nd=0; nd<n_donor; ++nd) {
            int resid = donor_resid[nd];
            hb_number1[resid] += pp(6,nd);
        }
        for(int na=0; na<n_acceptor; ++na) {
            int resid = acceptor_resid[na];
            hb_number1[resid] += pp(6,na+n_donor);
        }

        for( int i=0; i<n_res; i++) {
            float h = hb_number1[i];
            if (h == 0.0) continue;
            auto r = load_vec<2>(ramac, i);

            auto phiv1 = compact_sigmoid(r.v[0]-boundary_turn1, sharpness_turn1);
            auto phiv2 = compact_sigmoid(r.v[0]-boundary_turn2, sharpness_turn2);
            auto psiv1 = compact_sigmoid(r.v[1]-boundary_helix1, sharpness_helix1);
            auto psiv2 = compact_sigmoid(r.v[1]-boundary_helix2, sharpness_helix2);

            float turn_score   = (1.-phiv1.x())*phiv2.x();
            float psi_score    = (1.-psiv1.x())*psiv2.x();
            float helix_score  = (1.-turn_score) * psi_score;
            float sheet_score  = (1.-turn_score) * (1.-psi_score);

            dparams[0] += h*helix_score;
            dparams[1] += h*sheet_score;
            dparams[2] += h*turn_score;
            if(res_class[i]) {  // the residue class's offsets, after the shared values
                int o = 12 + 3*(res_class[i]-1);
                dparams[o+0] += h*helix_score;
                dparams[o+1] += h*sheet_score;
                dparams[o+2] += h*turn_score;
            }

            // the rama boundaries and sharpnesses act through turn_score and psi_score; a sigmoid's value
            // depends on its boundary b and sharpness s through (angle-b)*s
            float e_alpha  = E_alpha + class_offset(i,0);
            float e_beta   = E_beta  + class_offset(i,1);
            float e_other  = E_other + class_offset(i,2);
            float dE_dturn = e_other - e_alpha*psi_score - e_beta*(1.f-psi_score);
            float dE_dpsi  = (1.f-turn_score) * (e_alpha-e_beta);
            float x1 = r.v[0]-boundary_turn1,  x2 = r.v[0]-boundary_turn2;
            float y1 = r.v[1]-boundary_helix1, y2 = r.v[1]-boundary_helix2;
            dparams[4]  += h*dE_dturn *  phiv1.y()*phiv2.x();
            dparams[5]  += h*dE_dturn * -x1/sharpness_turn1*phiv1.y()*phiv2.x();
            dparams[6]  += h*dE_dturn * -(1.f-phiv1.x())*phiv2.y();
            dparams[7]  += h*dE_dturn *  (1.f-phiv1.x())*x2/sharpness_turn2*phiv2.y();
            dparams[8]  += h*dE_dpsi  *  psiv1.y()*psiv2.x();
            dparams[9]  += h*dE_dpsi  * -y1/sharpness_helix1*psiv1.y()*psiv2.x();
            dparams[10] += h*dE_dpsi  * -(1.f-psiv1.x())*psiv2.y();
            dparams[11] += h*dE_dpsi  *  (1.f-psiv1.x())*y2/sharpness_helix2*psiv2.y();
        }

        for( int i =0; i<n_res-1; i++ ) 
            hb_number2[i] = 0.0;
        for(int nd=0; nd<n_donor; ++nd) {
            int pepid = donor_resid[nd]-1;
            hb_number2[pepid] += pp(6,nd);
        }
        for(int na=0; na<n_acceptor; ++na) {
            int pepid = acceptor_resid[na];
            hb_number2[pepid] += pp(6,na+n_donor);
        }
        for( int i =0; i<n_res-1; i++ ) {
            auto dE     = compact_sigmoid(1.5-hb_number2[i], 2.0);
            dparams[3] += dE.x();
        }

        return dparams;
    }
#endif
    virtual DeviceNode* make_device(DeviceGraph& g, const vector<DerivComputation*>& systems) override {
        auto &cp = g.coord(protein_hbond), &cr = g.coord(rama);
        vector<int> cls, res_rama, first_donor, first_acceptor, count_start(1,0), count, pep_res, pep_start(1,0), pep_count, sys;
        for(int ns: range(g.n_system)) {
            auto& n = static_cast<HBondEnergy&>(*systems[ns]);
            if(n.params != params) throw string("the device engine shares hbond_energy's parameters between systems");
            int r0 = res_rama.size();  // this system's first residue
            // a residue's or peptide's elements, its donors then its acceptors, each in index order
            auto counted = [&](vector<int>& out, function<bool(int)> donor, function<bool(int)> acceptor) {
                for(int nd: range(n.n_donor))    if(donor(nd))    out.push_back(cp.offset[ns] + nd);
                for(int na: range(n.n_acceptor)) if(acceptor(na)) out.push_back(cp.offset[ns] + n.n_donor + na);
            };
            for(int i: range(n.n_res)) {
                cls.push_back(n.res_class[i]);
                res_rama.push_back(cr.offset[ns] + i);
                auto d = n.resid_donor.find(i), a = n.resid_acceptor.find(i);
                first_donor.push_back(d == n.resid_donor.end() ? -1 : cp.offset[ns] + d->second);
                first_acceptor.push_back(a == n.resid_acceptor.end() ? -1 : cp.offset[ns] + n.n_donor + a->second);
                counted(count, [&](int nd) {return n.donor_resid[nd]==i;}, [&](int na) {return n.acceptor_resid[na]==i;});
                count_start.push_back(count.size());
                sys.push_back(ns);
            }
            for(int i: range(n.n_res-1)) {
                pep_res.push_back(r0 + i);
                counted(pep_count, [&](int nd) {return n.donor_resid[nd]-1==i;}, [&](int na) {return n.acceptor_resid[na]==i;});
                pep_start.push_back(pep_count.size());
            }
        }
        return device_hbond_energy(g.potential_of(*this), cp, cr, params, cls, res_rama, first_donor, first_acceptor,
                count_start, count, pep_res, pep_start, pep_count, sys);
    }

    virtual void set_param(const std::vector<float>& new_param) override {
        if(int(new_param.size()) != n_param)
            throw string("hbond_energy set_param: wrong number of parameters");
        params = new_param;
        E_alpha          = params[0];
        E_beta           = params[1];
        E_other          = params[2];
        E_bias           = params[3];
        boundary_turn1   = params[4];
        sharpness_turn1  = params[5];
        boundary_turn2   = params[6];
        sharpness_turn2  = params[7];
        boundary_helix1  = params[8];
        sharpness_helix1 = params[9];
        boundary_helix2  = params[10];
        sharpness_helix2 = params[11];
    }
};
static RegisterNodeType<HBondEnergy,2> hbond_energy_node("hbond_energy");

