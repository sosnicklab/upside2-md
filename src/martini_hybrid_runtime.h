#pragma once

#include "deriv_engine.h"

namespace martini_hybrid {

bool project_bb_proxy_gradient_for_coord(
        const CoordNode& coord,
        VecArray pos_sens,
        int n_atom,
        int atom_idx,
        const Vec<3>& grad);

}  // namespace martini_hybrid
