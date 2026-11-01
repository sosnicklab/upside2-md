#ifndef ENVIRONMENT_INTERACTION_H
#define ENVIRONMENT_INTERACTION_H

// The pair potentials of environment.cpp's coverage nodes, written in the lane type S (edge_params.h) for the CPU's
// InteractionGraph and the device's

#include "vector_math.h"
#include "edge_params.h"
#include <string>

// internal linkage, as in the .cpp files they came from, so the CPU code generation is unchanged
namespace {

struct EnvironmentCoverageInteraction {
    // parameters are r0,r_sharpness, dot0,dot_sharpness

    constexpr static bool  symmetric = false;
    constexpr static int   n_param=4, n_dim1=6, n_dim2=4, simd_width=1;

    static float cutoff(const float* p) {
        return p[0] + compact_sigmoid_cutoff(p[1]);
    }

    template <typename I>  // Int4 on the CPU, int on the device
    UPSIDE_HD static I acceptable_id_pair(const I& id1, const I& id2) {
        auto sequence_exclude = I(2);  // exclude i,i, i,i+1, and i,i+2
        return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
    }

    template <typename S>
    UPSIDE_HD static S compute_edge(Vec<n_dim1,S> &d1, Vec<n_dim2,S> &d2, EdgeParams<S> p,
            const Vec<n_dim1,S> &cb_pos, const Vec<n_dim2,S> &sc_pos) {
        S one(1.f);
        auto displace = extract<0,3>(sc_pos)-extract<0,3>(cb_pos);
        auto rvec1    = extract<3,6>(cb_pos);
        auto prob     = sc_pos[3];

        auto dist2 = mag2(displace);
        auto inv_dist = rsqrt(dist2);
        auto dist = dist2*inv_dist;
        auto displace_unitvec = inv_dist*displace;

        // read parameters then transpose
        S r0, r_sharpness, dot0, dot_sharpness;
        p.load4(0, r0,r_sharpness,dot0,dot_sharpness);

        auto dp = dot(displace_unitvec,rvec1);
        auto radial_sig  = compact_sigmoid(dist-r0, r_sharpness);
        auto angular_sig = compact_sigmoid(dot0-dp, dot_sharpness);

        // now we compute derivatives (minus sign is from the derivative of angular_sig)
        auto d_displace = prob*(radial_sig.y()*angular_sig.x() * displace_unitvec -
                                radial_sig.x()*angular_sig.y()* inv_dist*(rvec1 - dp*displace_unitvec));

        store<3,6>(d1, -prob*radial_sig.x()*angular_sig.y()*displace_unitvec);
        store<0,3>(d1, -d_displace);
        store<0,3>(d2,  d_displace);
        auto score = radial_sig.x() * angular_sig.x();
        d2[3] = score;
        return prob * score;
    }

    static void param_deriv(Vec<n_param> &d_param, const float* p,
            const Vec<n_dim1> &hb_pos, const Vec<n_dim2> &sc_pos) {
        throw std::string("the coverage parameter derivative is not implemented");
    }

    static bool is_compatible(const float* p1, const float* p2) {return true;};
};

struct HbondEnvironmentCoverageInteraction {
    // parameters are r0,r_sharpness, dot0,dot_sharpness

    constexpr static bool  symmetric = false;
    constexpr static int   n_param=4, n_dim1=6, n_dim2=6, simd_width=1;

    static float cutoff(const float* p) {
	return p[0] + compact_sigmoid_cutoff(p[1]);
    }

    template <typename I>  // Int4 on the CPU, int on the device
    UPSIDE_HD static I acceptable_id_pair(const I& id1, const I& id2) {
        auto sequence_exclude = I(2);  // exclude i,i, i,i+1, and i,i+2
        return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
        //return Int4() == Int4();  // No exclusions (all true)
    }

    template <typename S>
    UPSIDE_HD static S compute_edge(Vec<n_dim1,S> &d1, Vec<n_dim2,S> &d2, EdgeParams<S> p,
	    const Vec<n_dim1,S> &x1, const Vec<n_dim2,S> &x2) {

	S one(1.f);
	auto displace = extract<0,3>(x2)-extract<0,3>(x1);
	auto rvec1 = extract<3,6>(x1);

	auto dist2 = mag2(displace);
	auto inv_dist = rsqrt(dist2);
	auto dist = dist2*inv_dist;
	auto displace_unitvec = inv_dist*displace;

	// read parameters then transpose
	S r0, r_sharpness, dot0, dot_sharpness;
	p.load4(0, r0,r_sharpness,dot0,dot_sharpness);

	auto dp = dot(displace_unitvec,rvec1);
	auto radial_sig  = compact_sigmoid(dist-r0, r_sharpness);
	auto angular_sig = compact_sigmoid(dot0-dp, dot_sharpness);

	// now we compute derivatives (minus sign is from the derivative of angular_sig)
	auto d_displace = radial_sig.y()*angular_sig.x() * displace_unitvec -
	       		radial_sig.x()*angular_sig.y()* inv_dist*(rvec1 - dp*displace_unitvec);

	store<3,6>(d1, -radial_sig.x()*angular_sig.y()*displace_unitvec);
	store<0,3>(d1, -d_displace);
	store<0,3>(d2,  d_displace);
	store<3,6>(d2,  make_zero<3,S>());  // the coverage does not depend on the second site's direction
	auto score = radial_sig.x() * angular_sig.x();
	return score;
    }

    static void param_deriv(Vec<n_param> &d_param, const float* p,
            const Vec<n_dim1> &hb_pos, const Vec<n_dim2> &sc_pos) {
        throw std::string("the coverage parameter derivative is not implemented");
    }

    static bool is_compatible(const float* p1, const float* p2) {return true;};
};

struct HbondEnvironmentCoverageInteraction2 {
    // parameters are r0,r_sharpness, dot0,dot_sharpness

    constexpr static bool  symmetric = false;
    constexpr static int   n_param=4, n_dim1=6, n_dim2=3, simd_width=1;

    static float cutoff(const float* p) {
	return p[0] + compact_sigmoid_cutoff(p[1]);
    }

    template <typename I>  // Int4 on the CPU, int on the device
    UPSIDE_HD static I acceptable_id_pair(const I& id1, const I& id2) {
        auto sequence_exclude = I(2);  // exclude i,i, i,i+1, and i,i+2
        return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
        //return Int4() == Int4();  // No exclusions (all true)
    }

    template <typename S>
    UPSIDE_HD static S compute_edge(Vec<n_dim1,S> &d1, Vec<n_dim2,S> &d2, EdgeParams<S> p,
	    const Vec<n_dim1,S> &x1, const Vec<n_dim2,S> &x2) {

	S one(1.f);
	auto displace = extract<0,3>(x2)-extract<0,3>(x1);
	auto rvec1 = extract<3,6>(x1);

	auto dist2 = mag2(displace);
	auto inv_dist = rsqrt(dist2);
	auto dist = dist2*inv_dist;
	auto displace_unitvec = inv_dist*displace;

	// read parameters then transpose
	S r0, r_sharpness, dot0, dot_sharpness;
	p.load4(0, r0,r_sharpness,dot0,dot_sharpness);

	auto dp = dot(displace_unitvec,rvec1);
	auto radial_sig  = compact_sigmoid(dist-r0, r_sharpness);
	auto angular_sig = compact_sigmoid(dot0-dp, dot_sharpness);

	// now we compute derivatives (minus sign is from the derivative of angular_sig)
	auto d_displace = radial_sig.y()*angular_sig.x() * displace_unitvec -
	       		radial_sig.x()*angular_sig.y()* inv_dist*(rvec1 - dp*displace_unitvec);

	store<3,6>(d1, -radial_sig.x()*angular_sig.y()*displace_unitvec);
	store<0,3>(d1, -d_displace);
	store<0,3>(d2,  d_displace);
	auto score = radial_sig.x() * angular_sig.x();
	return score;
    }

    static void param_deriv(Vec<n_param> &d_param, const float* p,
            const Vec<n_dim1> &hb_pos, const Vec<n_dim2> &sc_pos) {
        throw std::string("the coverage parameter derivative is not implemented");
    }

    static bool is_compatible(const float* p1, const float* p2) {return true;};
};
}

#endif
