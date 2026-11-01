// Device implementation of rama_map_pot.cpp's RamaMapPot
#include "device_nodes.h"

using namespace std;

namespace {

struct DeviceRamaMapPot : public DeviceNode {
    AccumView potential;
    DeviceCoord& rama;
    int nx, ny;
    DeviceArray<int> residue, map_id, system;
    DeviceArray<float> coefficients;  // LayeredPeriodicSpline2D<1>'s, 16 per (layer, x bin, y bin)

    DeviceRamaMapPot(AccumView potential_, DeviceCoord& rama_, const vector<int>& residue_, const vector<int>& map_id_,
            const vector<float>& coefficients_, int nx_, int ny_, const vector<int>& system_):
        potential(potential_), rama(rama_), nx(nx_), ny(ny_), residue(residue_), map_id(map_id_), system(system_),
        coefficients(coefficients_) {}

    void compute_value(ComputeMode mode) override {
        VecArray ramac = rama.value();
        auto rs = rama.sens.view(), pot = potential;
        bool want_pot = mode == PotentialAndDerivMode;
        int w = rama.row_width, gx = nx, gy = ny;
        const int *res = residue.ptr, *map = map_id.ptr, *sys = system.ptr;
        const float* coeff = coefficients.ptr;

        // add a litte paranoia to make sure there are no rounding problems
        const float scale = nx * (0.5f/M_PI_F - 1e-7f);
        const float shift = M_PI_F;

        parallel_for(int(residue.n), DEVICE_LAMBDA (int nr) {
            auto r = load_vec<2>(ramac, res[nr]);
            float x = (r.v[0]+shift)*scale, y = (r.v[1]+shift)*scale;

            // LayeredPeriodicSpline2D<1>::evaluate_value_and_deriv
            int x_bin = int(x);
            int y_bin = int(y);
            float fx = x - x_bin;
            float fy = y - y_bin;
            const float* c = coeff + (map[nr]*gx*gy + x_bin*gy + y_bin)*16;
            float value, dx, dy;
            spline_value_and_deriv(value, dx, dy, c, fx, fy);

            if(want_pot) pot.add(sys[nr], value);
            rs.add(res[nr]*w+0, dx * scale);
            rs.add(res[nr]*w+1, dy * scale);
        });
    }
};

}

DeviceNode* device_rama_map_pot(AccumView potential, DeviceCoord& rama, const vector<int>& residue,
        const vector<int>& map_id, const vector<float>& coefficients, int nx, int ny, const vector<int>& system) {
    return new DeviceRamaMapPot(potential, rama, residue, map_id, coefficients, nx, ny, system);
}
