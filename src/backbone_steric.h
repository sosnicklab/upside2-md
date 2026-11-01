#ifndef BACKBONE_STERIC_H
#define BACKBONE_STERIC_H

// The steric wall between backbone atoms of backbone_pairs (backbone_steric.cpp), for the CPU and the device

#include "vector_math.h"

// internal linkage, as in the .cpp file it came from, so the CPU code generation is unchanged
namespace {

constexpr float nonbonded_atom_cutoff2 = 3.f*3.f + 0.1f*3.f;

template <bool return_deriv>
UPSIDE_HD inline float
nonbonded_kernel_or_deriv_over_r(float r_mag2)
{
    const float energy_scale = 4.f;
    const float wall = 3.0f;  // corresponds to vdW *diameter*
    const float wall_squared = wall*wall;  
    const float width = 0.10f;
    const float sharpness = 1.f/(wall*width);  // ensure character

    const Vec<2> V = energy_scale*compact_sigmoid(r_mag2-wall_squared, sharpness);
    return return_deriv ? 2.f*V.y() : V.x();
}

template <typename I>  // Int4 on the CPU, int on the device
UPSIDE_HD I acceptable_backbone_pair(const I& id1, const I& id2) {
        auto sequence_exclude = I(1);
        return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
}
}

#endif
