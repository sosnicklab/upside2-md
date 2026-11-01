#ifndef EDGE_PARAMS_H
#define EDGE_PARAMS_H

// The parameters an InteractionGraph edge's interaction (an IType) reads, for the two lane types its compute_edge
// is written in: Float4, four edges at once with one parameter pointer each (the CPU), and float, one edge (the
// device).  The pair potentials are written once, in S, against these.

#include "vector_math.h"
#include "spline.h"

template <typename S> struct EdgeParams;

template <> struct EdgeParams<Float4> {  // host code; InteractionGraph names it in device files too
    const float* p[4];

    EdgeParams<Float4> operator+(int o) const {return EdgeParams<Float4>{{p[0]+o, p[1]+o, p[2]+o, p[3]+o}};}

    //! the parameters o to o+3 of each edge, parameter k in its own Float4 across the edges
    void load4(int o, Float4& a, Float4& b, Float4& c, Float4& d) const {
        a = Float4(p[0]+o); b = Float4(p[1]+o); c = Float4(p[2]+o); d = Float4(p[3]+o);
        transpose4(a,b,c,d);
    }
};

#ifndef __CUDACC__  // the four-lane spline evaluators (spline.h) are the CPU's
inline Vec<2,Float4> deBoor_value_and_deriv(EdgeParams<Float4> p, const Float4& x) {
    return deBoor_value_and_deriv(p.p, x);
}
inline Vec<2,Float4> clamped_deBoor_value_and_deriv(EdgeParams<Float4> p, const Float4& x, int n_knot) {
    return clamped_deBoor_value_and_deriv(p.p, x, n_knot);
}
#endif

template <> struct EdgeParams<float> {
    const float* p;

    UPSIDE_HD EdgeParams<float> operator+(int o) const {return EdgeParams<float>{p+o};}

    UPSIDE_HD void load4(int o, float& a, float& b, float& c, float& d) const {
        a = p[o]; b = p[o+1]; c = p[o+2]; d = p[o+3];
    }
};

//! a<b ? a : b, as SSE's min of two Float4
UPSIDE_HD inline float lane_min(float a, float b) {return a<b ? a : b;}
#ifndef __CUDACC__
inline Float4 lane_min(const Float4& a, const Float4& b) {return min(a,b);}
#endif

UPSIDE_HD inline Vec<2> deBoor_value_and_deriv(EdgeParams<float> p, float x) {
    return deBoor_lane(p.p, x);
}
UPSIDE_HD inline Vec<2> clamped_deBoor_value_and_deriv(EdgeParams<float> p, float x, int n_knot) {
    return clamped_deBoor_lane(p.p, x, n_knot);
}

#endif
