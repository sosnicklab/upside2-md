#ifndef EXTRA_PARTICLES_H
#define EXTRA_PARTICLES_H

#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"

using namespace h5;
using namespace std;

struct ExtraParticles : public PotentialNode {
    CoordNode& pos;
    int n_extra;
    vector<int> extra_indices;
    
    ExtraParticles(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        pos(pos_),
        n_extra(get_dset_size(1, grp, "extra_indices")[0]),
        extra_indices(n_extra) {
        
        // Read extra particle indices
        traverse_dset<1,int>(grp, "extra_indices", 
            [&](size_t i, int idx) {extra_indices[i] = idx;});
    }
    
    virtual void compute_value(ComputeMode mode) {
        // This is a placeholder - the actual forces will be computed by MARTINI
        // We just need to make sure these particles exist in the system
    }
};

static RegisterNodeType<ExtraParticles,1> extra_particles_node("extra_particles");

#endif // EXTRA_PARTICLES_H 