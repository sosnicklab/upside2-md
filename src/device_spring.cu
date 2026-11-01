// Device implementation of spring.cpp's Spring
#include "device_nodes.h"

using namespace std;

namespace {

struct DeviceSpring : public DeviceNode {
    AccumView potential;
    DeviceCoord& pos;
    int dim1;
    bool pbc;
    float box_len;
    DeviceArray<int> id, system;
    DeviceArray<float> equil_dist, spring_const;

    DeviceSpring(AccumView potential_, DeviceCoord& pos_, int dim1_, bool pbc_, float box_len_,
            const vector<int>& id_, const vector<float>& equil_dist_, const vector<float>& spring_const_,
            const vector<int>& system_):
        potential(potential_), pos(pos_), dim1(dim1_), pbc(pbc_), box_len(box_len_),
        id(id_), system(system_), equil_dist(equil_dist_), spring_const(spring_const_) {}

    void compute_value(ComputeMode mode) override {
        VecArray x = pos.value();
        auto ps = pos.sens.view();
        auto pot = potential;
        bool want_pot = mode == PotentialAndDerivMode;
        int w = pos.row_width, d1 = dim1;
        bool periodic = pbc;
        float box = box_len;
        const int *ids = id.ptr, *sys = system.ptr;
        const float *x0 = equil_dist.ptr, *k = spring_const.ptr;
        parallel_for(int(id.n), DEVICE_LAMBDA (int nt) {
            float dist   = x(d1, ids[nt]);
            float excess = dist-x0[nt];
            float sqr_excess = sqr(excess);

            if (excess == 0.f) return;

            if (periodic) {
                float excess1 = dist-x0[nt]-box;
                float excess2 = dist-x0[nt]+box;
                float sqr_excess1 = sqr(excess1);
                float sqr_excess2 = sqr(excess2);
                if (sqr_excess1 < sqr_excess) {
                    sqr_excess = sqr_excess1;
                    excess     = excess1;
                }
                if (sqr_excess2 < sqr_excess) {
                    sqr_excess = sqr_excess2;
                    excess     = excess2;
                }
            }

            if(want_pot) pot.add(sys[nt], 0.5f * k[nt] * sqr_excess);
            ps.add(ids[nt]*w+d1, k[nt] * excess);
        });
    }
};

}

DeviceNode* device_spring(AccumView potential, DeviceCoord& pos, int dim1, bool pbc, float box_len,
        const vector<int>& id, const vector<float>& equil_dist, const vector<float>& spring_const,
        const vector<int>& system) {
    return new DeviceSpring(potential, pos, dim1, pbc, box_len, id, equil_dist, spring_const, system);
}
