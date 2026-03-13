#ifndef LIGAND_H
#define LIGAND_H

#include "h5_support.h"
#include "vector_math.h"
#include <string>
#include <vector>

namespace ligand {

std::vector<int> read_ints(hid_t grp, const char* name);
std::vector<float> read_floats(hid_t grp, const char* name);
std::vector<std::string> read_strings(hid_t grp, const char* name);

void require_size(const char* name, size_t expected, size_t actual);
float clamp_unit(float x);

float pair_energy_and_deriv_over_r(
        float r2,
        float sigma,
        float epsilon,
        float qq,
        float cutoff2,
        float* deriv_over_r);

}  // namespace ligand

#endif
