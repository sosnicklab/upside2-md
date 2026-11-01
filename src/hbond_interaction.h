#ifndef HBOND_INTERACTION_H
#define HBOND_INTERACTION_H

// The pair potentials of hbond.cpp's protein_hbond and hbond_coverage, written in the lane type S (edge_params.h)
// for the CPU's InteractionGraph and the device's

#include "vector_math.h"
#include "edge_params.h"
#include "bead_interaction.h"
#include <string>

#define radial_cutoff2 (3.5f*3.5f)

// internal linkage, as in the .cpp files they came from, so the CPU code generation is unchanged
namespace {

template <typename S>
UPSIDE_HD Vec<2,S> hbond_radial_potential(const S& input,
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
UPSIDE_HD Vec<2,S> hbond_angular_potential(const S& dotp, const S& wall_dp, const S& inv_dp_width)
{
    Vec<2,S> v = sigmoid((dotp-wall_dp)*inv_dp_width);
    return make_vec2(v.x(), inv_dp_width*v.y());
}


struct ProteinHBondInteraction {
    // inner_barrier, inner_scale, outer_barrier, outer_scale, wall_dp, inv_dp_width
    // first group is donors; second group is acceptors
    constexpr static const bool symmetric=false;
    constexpr static const int n_param=8, n_dim1=6, n_dim2=6, simd_width=1;

    static float cutoff(const float* p) {
        return sqrtf(radial_cutoff2); // FIXME make parameter dependent
    }

    template <typename I>  // Int4 on the CPU, int on the device
    UPSIDE_HD static I acceptable_id_pair(const I& id1, const I& id2) {
        //return Int4() == Int4();  // No exclusions (all true)
        auto sequence_exclude = I(1);
        return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
    }

    template <typename S>
    UPSIDE_HD static S compute_edge(Vec<n_dim1,S> &d1, Vec<n_dim2,S> &d2, EdgeParams<S> p,
            const Vec<n_dim1,S> &x1, const Vec<n_dim2,S> &x2) {
        auto one = S(1.f);

        auto  H = extract<0,3>(x1);
        auto  O = extract<0,3>(x2);
        auto  rHN = extract<3,6>(x1);
        auto  rOC = extract<3,6>(x2);

        auto HO = H-O;

        auto magHO2 = mag2(HO) + S(1e-6f); // a bit of paranoia to avoid division by zero later
        auto invHOmag = rsqrt(magHO2);
        auto magHO    = magHO2 * invHOmag;  // avoid a sqrtf later

        auto rHO = HO*invHOmag;

        auto dotHOC =  dot(rHO,rOC);
        auto dotOHN = -dot(rHO,rHN);

        // FIXME I have to load up 4 of these rather pointlessly since they will all be the same
        S p0,p1,p2,p3; p.load4(0, p0,p1,p2,p3);
        S p4,p5,p6,p7; p.load4(4, p4,p5,p6,p7);

        auto radial   = hbond_radial_potential (magHO , p0, p1, p2, p3);  // x has val, y has deriv
        auto angular1 = hbond_angular_potential(dotHOC, p4, p5);
        auto angular2 = hbond_angular_potential(dotOHN, p4, p5);

        auto hb =  radial.x() * angular1.x() * angular2.x();
        auto c0 =  radial.y() * angular1.x() * angular2.x();
        auto c1 =  radial.x() * angular1.y() * angular2.x();
        auto c2 = -radial.x() * angular1.x() * angular2.y();

        auto drOC = c1*rHO;
        auto drHN = c2*rHO;

        auto dH = c0*rHO + (c1*invHOmag)*(rOC-dotHOC*rHO) + (c2*invHOmag)*(rHN+dotOHN*rHO);
        auto dO = -dH;

        auto hb_log = ternary(one<=hb, S(100.f), -logf(one-hb));  // work in multiplicative space

        auto deriv_prefactor = lane_min(rcp(one-hb),S(1e5f)); // FIXME this is a mess
        store<0,3>(d1, dH   * deriv_prefactor);
        store<3,6>(d1, drHN * deriv_prefactor);
        store<0,3>(d2, dO   * deriv_prefactor);
        store<3,6>(d2, drOC * deriv_prefactor);

        return hb_log;
    }

    static void param_deriv(Vec<n_param> &d_param, const float* p,
            const Vec<n_dim1> &x1, const Vec<n_dim2> &x2) {
        throw std::string("the protein_hbond parameter derivative is not implemented");
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

    template <typename I>
    UPSIDE_HD static I acceptable_id_pair(const I& id1, const I& id2) {
        // return Int4() == Int4();  // No exclusions (all true)
        // return id1 != id2; // exclude interactions on the same residue
        auto sequence_exclude = I(2);  // exclude i,i, i,i+1, and i,i+2
        return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
    }

    template <typename S>
    UPSIDE_HD static S compute_edge(Vec<n_dim1,S> &d1, Vec<n_dim2,S> &d2, EdgeParams<S> p,
            const Vec<n_dim1,S> &hb_pos, const Vec<n_dim2,S> &sc_pos) {
        S one(1.f);

        // print_vector("hb_pos[0]", hb_pos[0]);
        // print_vector("sc_pos[0]", sc_pos[0]);
        // print_vector("hbond_dist", mag(extract<0,3>(hb_pos)-extract<0,3>(sc_pos)));
        auto coverage = quadspline<n_knot_angular, n_knot>(d1,d2, inv_dtheta,inv_dx,p, hb_pos,sc_pos);

        auto prefactor = sqr(one-hb_pos[6]);
        d1 *= prefactor;
        d2 *= prefactor;
        d1[6] = -coverage * (one-hb_pos[6])*S(2.f);

        return prefactor * coverage;
    }

#ifndef __CUDACC__  // parameter derivatives are the CPU's (upside3.md 10.2)
    static void param_deriv(Vec<n_param> &d_param, const float* p,
            const Vec<n_dim1> &hb_pos, const Vec<n_dim2> &sc_pos) {
        quadspline_param_deriv<n_knot_angular, n_knot>(d_param, inv_dtheta,inv_dx,p, hb_pos,sc_pos);
        auto prefactor = sqr(1.f-hb_pos[6]);
        d_param *= prefactor;
    }
#endif

    static bool is_compatible(const float* p1, const float* p2) {return true;};
};
}

#endif
