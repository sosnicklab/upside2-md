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
};
static RegisterNodeType<Infer_H_O,1> infer_node("infer_H_O");


#define radial_cutoff2 (3.5f*3.5f)
// angular cutoff at 90 degrees
#define angular_cutoff (0.f)

template <typename S>
Vec<2,S> hbond_radial_potential(const S& input,
        const S& inner_barrier, const S& inv_inner_width,
        const S& outer_barrier, const S& inv_outer_width
        )
{
    Vec<2,S> outer_sigmoid = sigmoid((outer_barrier-input)*inv_outer_width);
    Vec<2,S> inner_sigmoid = sigmoid((input-inner_barrier)*inv_inner_width);

    return make_vec2( outer_sigmoid.x() * inner_sigmoid.x(),
            - inv_outer_width * outer_sigmoid.y() * inner_sigmoid.x()
            + inv_inner_width * inner_sigmoid.y() * outer_sigmoid.x());
}


template <typename S>
Vec<2,S> hbond_angular_potential(const S& dotp, const S& wall_dp, const S& inv_dp_width)
{
    Vec<2,S> v = sigmoid((dotp-wall_dp)*inv_dp_width);
    return make_vec2(v.x(), inv_dp_width*v.y());
}


namespace {
    struct ProteinHBondInteraction {
        // inner_barrier, inner_scale, outer_barrier, outer_scale, wall_dp, inv_dp_width
        // first group is donors; second group is acceptors
        constexpr static const bool symmetric=false;
        constexpr static const int n_param=8, n_dim1=6, n_dim2=6, simd_width=1;

        static float cutoff(const float* p) {
            return sqrtf(radial_cutoff2); // FIXME make parameter dependent
        }

        static Int4 acceptable_id_pair(const Int4& id1, const Int4& id2) {
            //return Int4() == Int4();  // No exclusions (all true)
            auto sequence_exclude = Int4(1);
            return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
        }

        static Float4 compute_edge(Vec<n_dim1,Float4> &d1, Vec<n_dim2,Float4> &d2, const float* p[4],
                const Vec<n_dim1,Float4> &x1, const Vec<n_dim2,Float4> &x2) {
            auto one = Float4(1.f);

            auto  H = extract<0,3>(x1);
            auto  O = extract<0,3>(x2);
            auto  rHN = extract<3,6>(x1);
            auto  rOC = extract<3,6>(x2);

            auto HO = H-O;

            auto magHO2 = mag2(HO) + Float4(1e-6f); // a bit of paranoia to avoid division by zero later
            auto invHOmag = rsqrt(magHO2);
            auto magHO    = magHO2 * invHOmag;  // avoid a sqrtf later

            auto rHO = HO*invHOmag;

            auto dotHOC =  dot(rHO,rOC);
            auto dotOHN = -dot(rHO,rHN);

            Vec<3,Float4> dH,dO,drHN,drOC;
            Float4 hb;
            auto within_angular_cutoff = (Float4(angular_cutoff) < dotHOC) & (Float4(angular_cutoff) < dotOHN);
            if(none(within_angular_cutoff)) {
                dH=dO=drHN=drOC=make_zero<3,Float4>();
                hb = zero<Float4>();
            } else {
                // FIXME I have to load up 4 of these rather pointlessly since they will all be the same
                // FIXME I don't know if I guarantee alignment on parameters
                // I expand to 8 to be sure
                auto p0 = Float4(p[0]);
                auto p1 = Float4(p[1]);
                auto p2 = Float4(p[2]);
                auto p3 = Float4(p[3]);   transpose4(p0,p1,p2,p3);
                auto p4 = Float4(p[0]+4);
                auto p5 = Float4(p[1]+4);
                auto p6 = Float4(p[2]+4);
                auto p7 = Float4(p[3]+4); transpose4(p4,p5,p6,p7);

                auto radial   = hbond_radial_potential (magHO , p0, p1, p2, p3);  // x has val, y has deriv
                auto angular1 = hbond_angular_potential(dotHOC, p4, p5);
                auto angular2 = hbond_angular_potential(dotOHN, p4, p5);

                hb      =  radial.x() * angular1.x() * angular2.x();
                auto c0 =  radial.y() * angular1.x() * angular2.x();
                auto c1 =  radial.x() * angular1.y() * angular2.x();
                auto c2 = -radial.x() * angular1.x() * angular2.y();

                drOC = c1*rHO;
                drHN = c2*rHO;

                dH = c0*rHO + (c1*invHOmag)*(rOC-dotHOC*rHO) + (c2*invHOmag)*(rHN+dotOHN*rHO);
                dO = -dH;
            }

            auto hb_log = ternary(one<=hb, Float4(100.f), -logf(one-hb));  // work in multiplicative space

            auto deriv_prefactor = min(rcp(one-hb),Float4(1e5f)); // FIXME this is a mess
            store<0,3>(d1, dH   * deriv_prefactor);
            store<3,6>(d1, drHN * deriv_prefactor);
            store<0,3>(d2, dO   * deriv_prefactor);
            store<3,6>(d2, drOC * deriv_prefactor);

            return hb_log;
        }

        static void param_deriv(Vec<n_param> &d_param, const float* p,
                const Vec<n_dim1> &x1, const Vec<n_dim2> &x2) {
            for(int np: range(n_param)) d_param[np] = -1.f;
        }

        static bool is_compatible(const float* p1, const float* p2) {return true;};
    };


    struct HBondCoverageInteraction {
        // radius scale angular_width angular_scale
        // first group is hb; second group is sc

        constexpr static bool  symmetric = false;
        constexpr static int   n_knot = N_KNOT_SC_BB, n_knot_angular=N_KNOT_ANGULAR;
        constexpr static int   n_param=2*n_knot_angular+2*n_knot, n_dim1=7, n_dim2=6, simd_width=1;
        constexpr static float inv_dx = 1.f/KNOT_SPACING_BB, inv_dtheta = (n_knot_angular-3)/2.f;

        static float cutoff(const float* p) {
            return (n_knot-2-1e-6)/inv_dx;  // 1e-6 insulates from roundoff
        }

        static Int4 acceptable_id_pair(const Int4& id1, const Int4& id2) {
            // return Int4() == Int4();  // No exclusions (all true)
            // return id1 != id2; // exclude interactions on the same residue
            auto sequence_exclude = Int4(2);  // exclude i,i, i,i+1, and i,i+2
            return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
        }

        static Float4 compute_edge(Vec<n_dim1,Float4> &d1, Vec<n_dim2,Float4> &d2, const float* p[4],
                const Vec<n_dim1,Float4> &hb_pos, const Vec<n_dim2,Float4> &sc_pos) {
            Float4 one(1.f);

            // print_vector("hb_pos[0]", hb_pos[0]);
            // print_vector("sc_pos[0]", sc_pos[0]);
            // print_vector("hbond_dist", mag(extract<0,3>(hb_pos)-extract<0,3>(sc_pos)));
            auto coverage = quadspline<n_knot_angular, n_knot>(d1,d2, inv_dtheta,inv_dx,p, hb_pos,sc_pos);

            auto prefactor = sqr(one-hb_pos[6]);
            d1 *= prefactor;
            d2 *= prefactor;
            d1[6] = -coverage * (one-hb_pos[6])*Float4(2.f);

            return prefactor * coverage;
        }

        static void param_deriv(Vec<n_param> &d_param, const float* p,
                const Vec<n_dim1> &hb_pos, const Vec<n_dim2> &sc_pos) {
            quadspline_param_deriv<n_knot_angular, n_knot>(d_param, inv_dtheta,inv_dx,p, hb_pos,sc_pos);
            auto prefactor = sqr(1.f-hb_pos[6]);
            d_param *= prefactor;
        }

        static bool is_compatible(const float* p1, const float* p2) {return true;};
    };
}

struct CTerHBond : public CoordNode {
    CoordNode& infer;
    CoordNode& bb;
    int n_donor;
    float cutoff;
    std::vector<int> id1, id2;
    std::vector<index_t> index1, index2;
    float inner_barrier, inner_scale, outer_barrier, outer_scale;

    CTerHBond(hid_t grp, CoordNode& infer_, CoordNode& bb_):
        CoordNode(get_dset_size(1,grp,"index1")[0], 1),
        infer(infer_), bb(bb_),
        n_donor(get_dset_size(1, grp, "index1")[0]),
        cutoff(5.5),
        id1(n_donor), id2(1),
        index1(n_donor), index2(1),
        inner_barrier(read_attribute<float>(grp, ".", "inner_barrier")),
        inner_scale(read_attribute<float>(grp, ".", "inner_scale")),
        outer_barrier(read_attribute<float>(grp, ".", "outer_barrier")),
        outer_scale(read_attribute<float>(grp, ".", "outer_scale"))
    {

        check_size(grp, "id1",      n_donor);
        check_size(grp, "index1",   n_donor);
        check_size(grp, "id2",      1);
        check_size(grp, "index2",   1);

        traverse_dset<1,int>(grp, "id1",    [&](size_t i, int x){id1[i] = x;});
        traverse_dset<1,int>(grp, "id2",    [&](size_t i, int x){id2[i] = x;});
        traverse_dset<1,int>(grp, "index1", [&](size_t i, int x){index1[i] = x;});
        traverse_dset<1,int>(grp, "index2", [&](size_t i, int x){index2[i] = x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("CTer_hbond"));

        VecArray vs = output;

        auto C = load_vec<3>(bb.output, index2[0]);
        for(int n=0; n<n_donor; ++n) {

            auto NH = load_vec<6>(infer.output, index1[n]);
            auto H = extract<0,3>(NH);
            auto disp = H - C;
            auto magHC = mag(disp);
            if (magHC <  cutoff and abs(id1[n]-id2[0]) > 1) {
                //auto radial = hbond_radial_potential(magHC, inner_barrier, inner_scale, outer_barrier, outer_scale);
                //vs(0,n) = radial.x();
                vs(0,n) = magHC;
            }
            else
                vs(0,n) = cutoff;
       }
    }

    virtual void propagate_deriv() override { }
};
static RegisterNodeType<CTerHBond,2> cter_hbond_node("cter_hbond");


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

        for(int nv: range(n_virtual)) vs(6,nv) = 1.f-expf(-vs(6,nv));
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
};
static RegisterNodeType<HBondCoverage,2> coverage_node("hbond_coverage");


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
        hb_number1(n_res),
        hb_number2(n_res-1),
        donor_resid(n_donor),
        acceptor_resid(n_acceptor),
        Ehbond(n_res),
        dPhi(n_res),
        dPsi(n_res)
    {
        assert(n_donor+n_acceptor==n_virtual);
        traverse_dset<1, float> (grp, "parameters",     [&](size_t i, float x) {params[i] = x;});
        traverse_dset<1, int>   (grp, "donor_resid",    [&](size_t i, int x)   {donor_resid[i] = x;    resid_donor.insert(pair<int, int>(x, (int)i));  });
        traverse_dset<1, int>   (grp, "acceptor_resid", [&](size_t i, int x)   {acceptor_resid[i] = x; resid_acceptor.insert(pair<int, int>(x, (int)i));});
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

            Ehbond[i] = E_alpha * helix_score + E_beta * sheet_score + E_other * turn_score;
            dPhi[i]   = E_alpha * dhelix_dphi + E_beta * dsheet_dphi + E_other * dturn_dphi;
            dPsi[i]   = E_alpha * dhelix_dpsi + E_beta * dsheet_dpsi;

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

        n_hbond = 0.0;
        potential = 0.0;
        for( int i =0; i<n_res; i++ ) {
            if (hb_number1[i] == 0.0)
                continue;
            n_hbond      += hb_number1[i];
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
        std::vector<float> dparams(12);

        std::vector<float> dE_other(n_res);
        std::vector<float> dE_alpha(n_res);
        std::vector<float> dE_beta(n_res);

        for( int i =0; i<12; i++ ) 
            dparams[i] = 0.0;

        VecArray ramac = rama.output;
        VecArray pp    = protein_hbond.output;
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

            dE_alpha[i] = helix_score;
            dE_beta[i]  = sheet_score;
            dE_other[i] = turn_score;
        }

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

        for( int i =0; i<n_res; i++ ) {
            dparams[0] += hb_number1[i]*dE_alpha[i];
            dparams[1] += hb_number1[i]*dE_beta[i];
            dparams[2] += hb_number1[i]*dE_other[i];
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
    virtual void set_param(const std::vector<float>& new_param) override {
        // FIXME add a check to the size of parameters
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

