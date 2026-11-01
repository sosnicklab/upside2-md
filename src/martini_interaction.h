#ifndef MARTINI_INTERACTION_H
#define MARTINI_INTERACTION_H

// The terms of martini_potential.cpp's nodes that they share with their device implementations (device_martini.cu).
// martini_sc_table_1body's table and the two springs keep their own code on the CPU: shared from here, the same
// lines compiled to other roundings there (findings 11g), so the device has its own copy of those.

#include "vector_math.h"
#include "spline.h"
#include "box.h"

// internal linkage, as in the .cpp file they came from, so the CPU code generation is unchanged
namespace {

// martini_potential's pair of atoms at pa and pb: the combined LJ and Coulomb energy, a LayeredClampedSpline1D<1>'s
// arrays tabulated on n_grid points from r_min to r_max, at their minimum-image distance, and the force on the first
// atom (on the second its negative).  False, with no energy and force, past the cutoff or where both vanish.
UPSIDE_HD inline bool martini_pair_force(float& pair_potential, Vec<3>& pair_force, const Vec<3>& pa,
        const Vec<3>& pb, const float* coefficients, const float* left_clamped_value,
        const float* right_clamped_value, float cutoff, float r_min, float r_max, int n_grid,
        float box_x, float box_y, float box_z) {
    const float kMinDistance = 1.0e-6f;
    auto dr = pa - pb;
    dr = simulation_box::minimum_image(dr, box_x, box_y, box_z);
    auto dist2 = mag2(dr);
    auto dist = sqrtf(dist2 < kMinDistance ? kMinDistance : dist2);
    pair_potential = 0.f;
    pair_force = make_zero<3>();

    if(dist > cutoff) return false;

    float r_coord = (dist - r_min) / (r_max - r_min) * float(n_grid - 1);
    float result[2];
    clamped_spline_value_and_deriv<1>(result, coefficients, left_clamped_value, right_clamped_value, n_grid, 0,
            r_coord);
    float coord_scale = float(n_grid - 1) / (r_max - r_min);
    float dE_dr = result[0] * coord_scale;
    pair_potential = result[1];
    pair_force = (-dE_dr/dist) * dr;
    return !(pair_potential == 0.f && mag2(pair_force) == 0.f);
}

// Uniform cubic B-spline weights and their derivatives at fraction t of a knot interval.
UPSIDE_HD inline void cubic_bspline_weights(float w[4], float dw[4], float t) {
    float u = 1.f - t;
    w[0] = (1.f/6.f)*u*u*u;
    w[1] = (1.f/6.f)*(4.f + t*t*(3.f*t - 6.f));
    w[2] = (1.f/6.f)*(1.f + 3.f*t*(1.f + t*(1.f - t)));
    w[3] = (1.f/6.f)*t*t*t;
    dw[0] = -0.5f*u*u;
    dw[1] = t*(1.5f*t - 2.f);
    dw[2] = 0.5f + t*(1.f - 1.5f*t);
    dw[3] = 0.5f*t*t;
}

}

#endif
