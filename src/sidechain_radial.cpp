#include "deriv_engine.h"
#include <string>
#include "timing.h"
#include "affine.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include "interaction_graph.h"
#include "spline.h"
#include "state_logger.h"

using namespace std;
using namespace h5;

namespace {
template <bool is_symmetric>
struct RadialHelper {
    // spline-based distance interaction
    // n_knot is the number of basis splines (including those required to get zero
    //   derivative in the clamped spline)
    // spline is constant over [0,dx] to avoid funniness at origin

    // Please obey these 4 conditions:
    // p[0] = 1./dx, that is the inverse of the knot spacing
    // should have p[1] == p[3] for origin clamp (p[0] is inv_dx)
    // should have p[-3] == p[-1] (negative indices from the end, Python-style) for terminal clamping
    // should have (1./6.)*p[-3] + (2./3.)*p[-2] + (1./6.)*p[-1] == 0. for continuity at cutoff

    constexpr static bool  symmetric = is_symmetric;
    constexpr static int   n_knot=16, n_param=1+n_knot, n_dim1=3, n_dim2=3, simd_width=1;

    static float cutoff(const float* p) {
        const float inv_dx = p[0];
        return (n_knot-2-1e-6)/inv_dx;  // 1e-6 just insulates us from round-off error
    }

    static bool is_compatible(const float* p1, const float* p2) {
        if(symmetric) for(int i: range(n_param)) if(p1[i]!=p2[i]) return false;
        return true;
    }

    static Int4 acceptable_id_pair(const Int4& id1, const Int4& id2) {
        auto sequence_exclude = Int4(2);
        return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
    }

    static Float4 compute_edge(Vec<n_dim1,Float4> &d1, Vec<n_dim2,Float4> &d2, const float* p[4], 
            const Vec<n_dim1,Float4> &x1, const Vec<n_dim2,Float4> &x2) {
        alignas(16) const float inv_dx_data[4] = {p[0][0], p[1][0], p[2][0], p[3][0]};

        auto inv_dx     = Float4(inv_dx_data, Alignment::aligned);
        auto disp       = x1-x2;
        auto dist2      = mag2(disp);
        auto inv_dist   = rsqrt(dist2+Float4(1e-7f));  // 1e-7 is divergence protection
        auto dist_coord = dist2*(inv_dist*inv_dx);

        const float* pp[4] = {p[0]+1, p[1]+1, p[2]+1, p[3]+1};
        auto en = clamped_deBoor_value_and_deriv(pp, dist_coord, n_knot);
        d1 = disp*(inv_dist*inv_dx*en.y());
        d2 = -d1;
        return en.x();
    }

    static void param_deriv(Vec<n_param> &d_param, const float* p,
            const Vec<n_dim1> &x1, const Vec<n_dim2> &x2) {
        d_param = make_zero<n_param>();
        float inv_dx = p[0];
        auto dist_coord = inv_dx*mag(x1-x2); // to convert to spline coords of interger grid of knots
        auto dV_dinv_dx = clamped_deBoor_value_and_deriv(p+1, dist_coord, n_knot).y()*mag(x1-x2);
        d_param[0] = dV_dinv_dx;
 
        int starting_bin;
        float result[4];
        clamped_deBoor_coeff_deriv(&starting_bin, result, dist_coord, n_knot);
        for(int i: range(4)) d_param[1+starting_bin+i] = result[i];
   }
};



struct SidechainRadialPairs : public PotentialNode
{

    InteractionGraph<RadialHelper<true>> igraph;

    SidechainRadialPairs(hid_t grp, CoordNode& bb_point_):
        PotentialNode(),
        igraph(grp, &bb_point_)
    {};

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("radial_pairs"));

        igraph.compute_edges();
        for(int ne=0; ne<igraph.n_edge; ++ne) igraph.edge_sensitivity[ne] = 1.f;
        igraph.propagate_derivatives();

        if(mode==PotentialAndDerivMode) {
            potential = 0.f;
            for(int ne=0; ne<igraph.n_edge; ++ne) 
                potential += igraph.edge_value[ne];
        }
    }
};

struct HBondSidechainRadialPairs : public PotentialNode
{

    InteractionGraph<RadialHelper<false>> igraph;

    HBondSidechainRadialPairs(hid_t grp, CoordNode& hb_point_, CoordNode& bb_point_):
        PotentialNode(),
        igraph(grp, &hb_point_, &bb_point_)
    {};

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("hbond_sc_radial_pairs"));

        igraph.compute_edges();
        for(int ne=0; ne<igraph.n_edge; ++ne) igraph.edge_sensitivity[ne] = 1.f;
        igraph.propagate_derivatives();

        if(mode==PotentialAndDerivMode) {
            potential = 0.f;
            for(int ne=0; ne<igraph.n_edge; ++ne) 
                potential += igraph.edge_value[ne];
        }
    }

    virtual std::vector<float> get_param() const override {return igraph.get_param();}
#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {return igraph.get_param_deriv();}
#endif
    virtual void set_param(const std::vector<float>& new_param) override {igraph.set_param(new_param);}
};

template <bool is_symmetric>
struct ContactInteraction {

    constexpr static bool  symmetric = is_symmetric;
    constexpr static int   n_param = 4, n_dim1=3, n_dim2=3, simd_width=1;

    static float cutoff(const float* p) {
        return p[3]+1e-6;
    }

    static Int4 acceptable_id_pair(const Int4& id1, const Int4& id2) {
        auto sequence_exclude = Int4(2);  // exclude i,i, i,i+1, and i,i+2
        return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
    }

    static bool is_compatible(const float* p1, const float* p2) { return true;}

    static Float4 compute_edge(Vec<n_dim1,Float4> &d1, Vec<n_dim2,Float4> &d2, const float* p[4], 
                const Vec<n_dim1,Float4> &x1, const Vec<n_dim2,Float4> &x2) {

        auto disp       = x1-x2;
        auto dist2      = mag2(disp);
        auto inv_dist   = rsqrt(dist2+Float4(1e-7f));  // 1e-7 is divergence protection
        auto dist       = dist2*inv_dist;
        auto scale      = Float4(p[0]);
        auto center     = Float4(p[1]);
        auto sharpness  = Float4(p[2]);
        auto empty      = Float4(p[3]);
        transpose4(scale,  center,  sharpness,  empty);
        auto sig  = compact_sigmoid(dist-center,  sharpness);
        d1 = disp*inv_dist*scale*sig.y();
        d2 = -d1;
        return scale*sig.x();
    }

    static void param_deriv(Vec<n_param> &d_param, const float* p,
            const Vec<n_dim1> &x1, const Vec<n_dim2> &x2) {
        d_param     = make_zero<n_param>();
        auto dist   = mag(x1-x2);
        auto out    = compact_sigmoid(dist-p[1], p[2]);
        d_param[0] +=  out.x();
        d_param[1] += -p[0] * out.y();
        d_param[2] +=  p[0] * (dist-p[1]) * out.y()/p[2];
   }

};

struct SymmContactNListEnergy : public PotentialNode
{
    InteractionGraph<ContactInteraction<true>> igraph;
    SymmContactNListEnergy(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        igraph(grp, &pos_)
    {};
    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("fast_contact"));
        igraph.compute_edges();
        for(int ne=0; ne<igraph.n_edge; ++ne) igraph.edge_sensitivity[ne] = 1.f;
        igraph.propagate_derivatives();
        if(mode==PotentialAndDerivMode) {
            potential = 0.f;
            for(int ne=0; ne<igraph.n_edge; ++ne) 
                potential += igraph.edge_value[ne];
        }
    }
    virtual std::vector<float> get_param() const override {return igraph.get_param();}
#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {return igraph.get_param_deriv();}
#endif
    virtual void set_param(const std::vector<float>& new_param) override {igraph.set_param(new_param);}
};

struct NoSymmContactNListEnergy : public PotentialNode
{
    InteractionGraph<ContactInteraction<false>> igraph;
    NoSymmContactNListEnergy(hid_t grp, CoordNode& posA_, CoordNode& posB_):
        PotentialNode(),
        igraph(grp, &posA_, &posB_)
    {};
    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("fast contact"));
        igraph.compute_edges();
        for(int ne=0; ne<igraph.n_edge; ++ne) igraph.edge_sensitivity[ne] = 1.f;
        igraph.propagate_derivatives();
        if(mode==PotentialAndDerivMode) {
            potential = 0.f;
            for(int ne=0; ne<igraph.n_edge; ++ne) 
                potential += igraph.edge_value[ne];
        }
    }
    virtual std::vector<float> get_param() const override {return igraph.get_param();}
#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {return igraph.get_param_deriv();}
#endif
    virtual void set_param(const std::vector<float>& new_param) override {igraph.set_param(new_param);}
};

struct WallSpring : public PotentialNode
{
    struct Params {
        index_t atom[2];
        float equil_dist;
        float spring_constant;
        float sign;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;

    WallSpring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(2, grp, "id")[0]), pos(pos_), params(n_elem)
    {
        int n_dep = 2;  // number of atoms that each term depends on 
        check_size(grp, "id",           n_elem, n_dep);
        check_size(grp, "signs",        n_elem);
        check_size(grp, "equil_dist",   n_elem);
        check_size(grp, "spring_const", n_elem);

        auto& p = params;
        traverse_dset<2,int>  (grp, "id",           [&](size_t i, size_t j, int   x) {p[i].atom[j] = x;});
        traverse_dset<1,float>(grp, "signs",        [&](size_t i,           float x) {p[i].sign = x;});
        traverse_dset<1,float>(grp, "equil_dist",   [&](size_t i,           float x) {p[i].equil_dist = x;});
        traverse_dset<1,float>(grp, "spring_const", [&](size_t i,           float x) {p[i].spring_constant = x;});

        if(logging(LOG_DETAILED))
            default_logger->add_logger<float>("wall_energy", {1}, [&](float* buffer) {
                    float pot = 0.f;
                    VecArray x = pos.output;

                    for(int nt=0; nt<n_elem; ++nt) {

                        auto p = params[nt];
                        float dmag = mag(load_vec<3>(x, p.atom[0]) - load_vec<3>(x, p.atom[1]));

                        auto dist = dmag - p.equil_dist;
                        if ( p.sign*dist >= 0.0) pot += 0.5f * p.spring_constant * sqr(dist);
                    }
                    buffer[0] = pot;
                    });
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("walls"));

        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;

        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];

            auto x1 = load_vec<3>(posc, p.atom[0]);
            auto x2 = load_vec<3>(posc, p.atom[1]);

            auto disp = x1 - x2;
            auto dist = mag(disp);

            if (dist*p.sign > 0.0) {

                auto deriv = p.spring_constant * (1.f - p.equil_dist*inv_mag(disp)) * disp;
                if(pot) *pot += 0.5f * p.spring_constant * sqr(mag(disp) - p.equil_dist);

                update_vec(pos_sens, p.atom[0],  deriv);
                update_vec(pos_sens, p.atom[1], -deriv);
            }
        }
    }
};

struct ContactEnergy : public PotentialNode
{
    struct Param {
        index_t   loc[2];
        float     energy;
        float     dist;
        float     scale;  // 1.f/width
        float     cutoff;
    };

    int n_contact;
    CoordNode& bead_pos;
    vector<Param> params;
    float cutoff;

    ContactEnergy(hid_t grp, CoordNode& bead_pos_):
        PotentialNode(),
        n_contact(get_dset_size(2, grp, "id")[0]),
        bead_pos(bead_pos_), 
        params(n_contact)
    {
        check_size(grp, "id",       n_contact, 2);
        check_size(grp, "energy",   n_contact);
        check_size(grp, "distance", n_contact);
        check_size(grp, "width",    n_contact);

        traverse_dset<2,int  >(grp, "id",       [&](size_t nc, size_t i, int x){params[nc].loc[i] = x;});
        traverse_dset<1,float>(grp, "distance", [&](size_t nc, float x){params[nc].dist = x;});
        traverse_dset<1,float>(grp, "energy",   [&](size_t nc, float x){params[nc].energy = x;});
        traverse_dset<1,float>(grp, "width",    [&](size_t nc, float x){params[nc].scale = 1.f/x;});
        for(auto &p: params) p.cutoff = p.dist + 1.f/p.scale;

        if(logging(LOG_DETAILED)) {
            default_logger->add_logger<float>("contact_energy", {bead_pos.n_elem}, 
                    [&](float* buffer) {
                       fill_n(buffer, bead_pos.n_elem, 0.f);
                       VecArray pos  = bead_pos.output;

                       for(const auto &p: params) {
                           auto dist = mag(load_vec<3>(pos, p.loc[0]) - load_vec<3>(pos, p.loc[1]));
                           float en = p.energy * compact_sigmoid(dist-p.dist, p.scale)[0];
                           buffer[p.loc[0]] += 0.5f*en;
                           buffer[p.loc[1]] += 0.5f*en;
                       }});
        }
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("contact_energy"));
        VecArray pos  = bead_pos.output;
        VecArray sens = bead_pos.sens;
        potential = 0.f;

        for(int nc=0; nc<n_contact; ++nc) {
            const auto& p = params[nc];
            auto disp = load_vec<3>(pos, p.loc[0]) - load_vec<3>(pos, p.loc[1]);
            auto dist = mag(disp);
            if(dist>=p.cutoff) continue;

            Vec<2> contact = compact_sigmoid(dist-p.dist, p.scale);
            potential += p.energy*contact.x();
            auto deriv = (p.energy*contact.y()*rcp(dist)) * disp;
            update_vec(sens, p.loc[0],  deriv);
            update_vec(sens, p.loc[1], -deriv);
        }
    }
};

struct CooperationContacts : public PotentialNode
{
    struct Param {
        index_t   loc[2];
        float     energy;
        float     dist;
        float     scale;  // 1.f/width
        float     cutoff;
    };

    float energy;
    int n_contact;
    CoordNode& bead_pos;
    vector<Param> params;
    vector<float> dists;
    vector<Vec<2> > contacts;
    float cutoff;

    CooperationContacts(hid_t grp, CoordNode& bead_pos_):
        PotentialNode(),
        n_contact(get_dset_size(2, grp, "id")[0]),
        bead_pos(bead_pos_), 
        params(n_contact),
        dists(n_contact),
        contacts(n_contact)
    {
        check_size(grp, "id",       n_contact, 2);
        check_size(grp, "distance", n_contact);
        check_size(grp, "width",    n_contact);
        check_size(grp, "energy",   n_contact);

        traverse_dset<2,int  >(grp, "id",       [&](size_t nc, size_t i, int x){params[nc].loc[i] = x;});
        traverse_dset<1,float>(grp, "distance", [&](size_t nc, float x){params[nc].dist = x;});
        traverse_dset<1,float>(grp, "width",    [&](size_t nc, float x){params[nc].scale = 1.f/x;});
        traverse_dset<1,float>(grp, "energy",   [&](size_t nc, float x){params[nc].energy = x;});

        for(auto &p: params) p.cutoff = p.dist + 1.f/p.scale;

        energy = 1.0;
        for(auto &p: params) energy *= p.energy;
        energy = pow(energy, 1.f/n_contact);

        if(logging(LOG_DETAILED)) {
            default_logger->add_logger<float>("cooperation_contacts", {bead_pos.n_elem}, 
                    [&](float* buffer) {
                       fill_n(buffer, bead_pos.n_elem, 0.f);
                       VecArray pos  = bead_pos.output;

                       float en = 1.f;
                       for(const auto &p: params) {
                           auto dist = mag(load_vec<3>(pos, p.loc[0]) - load_vec<3>(pos, p.loc[1]));
                           en *= compact_sigmoid(dist-p.dist, p.scale)[0];
                       }
                       en = energy * pow(en, 1.f/n_contact);

                       for(const auto &p: params) {
                           buffer[p.loc[0]] += 0.5f*en;
                           buffer[p.loc[1]] += 0.5f*en;
                       }
                  }
             );
         }
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("cooperation_contacts"));

        VecArray pos  = bead_pos.output;
        VecArray sens = bead_pos.sens;

        float product_contact = 1.f;

        // calc potential energy
        for(int nc=0; nc<n_contact; ++nc) {
            const auto& p = params[nc];
            auto disp = load_vec<3>(pos, p.loc[0]) - load_vec<3>(pos, p.loc[1]);
            auto dist = mag(disp);
            dists[nc] = dist;
            Vec<2> contact = compact_sigmoid(dist-p.dist, p.scale);
            contacts[nc] = contact;
            product_contact *= contact.x();
        }

        // S**(1/N)
        float exponent = 1.f/n_contact;
        product_contact = pow(product_contact, exponent);
        potential = product_contact * energy;

        // calc force
        if( product_contact > 0.f and product_contact<1.f) {
            for(int nc=0; nc<n_contact; ++nc) {
                const auto& p = params[nc];
                auto dist = dists[nc];

                auto disp = load_vec<3>(pos, p.loc[0]) - load_vec<3>(pos, p.loc[1]);
                auto contact = contacts[nc];
                auto deriv = exponent * (potential*contact.y()*rcp(contact.x())*rcp(dist)) * disp;

                update_vec(sens, p.loc[0],  deriv);
                update_vec(sens, p.loc[1], -deriv);
            }
        }
    }
};

}

static RegisterNodeType<ContactEnergy,1>             contact_node("contact");
static RegisterNodeType<SymmContactNListEnergy,1>    sfast_contact_node("symm_fast_contact");
static RegisterNodeType<NoSymmContactNListEnergy,2>  nsfast_contact_node("nosymm_fast_contact");
static RegisterNodeType<CooperationContacts,1>       cooperation_contacts_node("cooperation_contacts");
static RegisterNodeType<SidechainRadialPairs,1>      radial_node ("radial");
static RegisterNodeType<HBondSidechainRadialPairs,2> hbond_sc_radial_node ("hbond_sc_radial");
