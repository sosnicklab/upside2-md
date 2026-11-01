// Device implementation of placement.cpp's PlacementNode, for its fixed and Ramachandran-dependent placements
#include "device_nodes.h"
#include "affine.h"

using namespace std;

namespace {

constexpr int max_pos_dim = 8;
enum DevicePlaceT {SCALAR=0, VECTOR=1, POINT=2};  // PlaceT's order

struct DevicePlacement : public DeviceNode {
    DeviceCoord &out, &alignment;
    DeviceCoord* rama;
    int n_pos_dim, n_kind, nx, ny;
    int kinds[3];
    DeviceArray<int> affine_residue, layer, rama_residue;
    DeviceArray<float> table;
    DeviceArray<float> rama_deriv;  // per element: the values' derivatives by phi, then by psi

    DevicePlacement(DeviceCoord& out_, DeviceCoord& alignment_, DeviceCoord* rama_, const vector<int>& kinds_,
            int n_pos_dim_, const vector<int>& affine_residue_, const vector<int>& layer_, const vector<float>& table_,
            int nx_, int ny_, const vector<int>& rama_residue_):
        out(out_), alignment(alignment_), rama(rama_), n_pos_dim(n_pos_dim_), n_kind(kinds_.size()), nx(nx_), ny(ny_),
        affine_residue(affine_residue_), layer(layer_), rama_residue(rama_residue_), table(table_),
        rama_deriv(rama_ ? 2*n_pos_dim_*affine_residue_.size() : 0)
    {
        if(n_kind > 3 || n_pos_dim > max_pos_dim) throw string("placement signature too long for the device");
        for(int i=0; i<n_kind; ++i) kinds[i] = kinds_[i];
    }

    void compute_value(ComputeMode mode) override {
        VecArray al = alignment.value();
        VecArray rc = rama ? rama->value() : VecArray();
        float* o = out.output.ptr;
        int w = out.row_width, nd = n_pos_dim, nk = n_kind, gx = nx, gy = ny;
        int k0 = kinds[0], k1 = kinds[1], k2 = kinds[2];
        bool from_rama = rama;
        const int *res = affine_residue.ptr, *lay = layer.ptr, *rr = rama_residue.ptr;
        const float* tab = table.ptr;
        float* rd = rama_deriv.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int ne) {
            auto aff = load_vec<7>(al, res[ne]);
            auto t   = extract<0,3>(aff);
            float U[9]; quat_to_rot(U, aff.v+3);

            float val[max_pos_dim];
            if(from_rama) {
                // RamaPlacement::evaluate, through LayeredPeriodicSpline2D::evaluate_value_and_deriv
                const float scale_x = gx * (0.5f/M_PI_F - 1e-7f);
                const float scale_y = gy * (0.5f/M_PI_F - 1e-7f);
                const float shift = M_PI_F;
                auto r = load_vec<2>(rc, rr[ne]);
                float x = (r[0]+shift)*scale_x, y = (r[1]+shift)*scale_y;
                int x_bin = int(x);
                int y_bin = int(y);
                float fx = x - x_bin;
                float fy = y - y_bin;
                const float* c = tab + (lay[ne]*gx*gy + x_bin*gy + y_bin)*16*nd;
                for(int id=0; id<nd; ++id)
                    spline_value_and_deriv(val[id], rd[(2*ne+0)*nd+id], rd[(2*ne+1)*nd+id], c+id*16, fx, fy);
            } else {
                for(int d=0; d<nd; ++d) val[d] = tab[lay[ne]*nd+d];
            }

            // do_transformations
            int offset = 0;
            for(int k=0; k<nk; ++k) {
                int kind = k==0 ? k0 : k==1 ? k1 : k2;
                if(kind == SCALAR) {
                    o[ne*w+offset] = val[offset];
                    offset += 1;
                } else {
                    auto v = make_vec3(val[offset],val[offset+1],val[offset+2]);
                    store_vec(o+ne*w+offset, kind==VECTOR ? apply_rotation(U, v) : apply_affine(U,t, v));
                    offset += 3;
                }
            }
        });
    }

    void propagate_deriv() override {
        VecArray al = alignment.value();
        auto s = out.sens.view(), as = alignment.sens.view();
        AccumView rs = rama ? rama->sens.view() : AccumView{nullptr, nullptr};
        const float* o = out.output.ptr;
        int w = out.row_width, wa = alignment.row_width, wr = rama ? rama->row_width : 0;
        int nd = n_pos_dim, nk = n_kind, gx = nx, gy = ny;
        int k0 = kinds[0], k1 = kinds[1], k2 = kinds[2];
        bool from_rama = rama;
        const int *res = affine_residue.ptr, *rr = rama_residue.ptr;
        const float* rd = rama_deriv.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int ne) {
            auto aff = load_vec<7>(al, res[ne]);
            auto t   = extract<0,3>(aff);
            float U[9]; quat_to_rot(U, aff.v+3);

            float ref_sens[max_pos_dim];
            Vec<3> com_deriv = make_zero<3>();
            Vec<3> torque    = make_zero<3>();

            // do_sens_transformations
            int offset = 0;
            for(int k=0; k<nk; ++k) {
                int kind = k==0 ? k0 : k==1 ? k1 : k2;
                if(kind == SCALAR) {
                    ref_sens[offset] = s.get(ne*w+offset);
                    offset += 1;
                    continue;
                }
                auto sv = make_vec3(s.get(ne*w+offset), s.get(ne*w+offset+1), s.get(ne*w+offset+2));
                auto xv = make_vec3(o[ne*w+offset], o[ne*w+offset+1], o[ne*w+offset+2]);
                store_vec(ref_sens+offset, apply_inverse_rotation(U,sv));
                if(kind == POINT) {
                    com_deriv += sv;
                    torque    += cross(xv-t, sv);
                } else {
                    torque += cross(xv, sv);
                }
                offset += 3;
            }

            if(from_rama) {
                // RamaPlacement::propagate_deriv
                const float scale_x = gx * (0.5f/M_PI_F - 1e-7f);
                const float scale_y = gy * (0.5f/M_PI_F - 1e-7f);
                float dphi = 0.f, dpsi = 0.f;
                for(int id=0; id<nd; ++id) {
                    dphi += ref_sens[id]*rd[(2*ne+0)*nd+id];
                    dpsi += ref_sens[id]*rd[(2*ne+1)*nd+id];
                }
                rs.add(rr[ne]*wr+0, scale_x*dphi);
                rs.add(rr[ne]*wr+1, scale_y*dpsi);
            }

            for(int d=0; d<3; ++d) {
                as.add(res[ne]*wa+d,   com_deriv[d]);
                as.add(res[ne]*wa+3+d, torque[d]);
            }
        });
    }
};

}

DeviceNode* device_placement(DeviceCoord& out, DeviceCoord& alignment, DeviceCoord* rama, const vector<int>& kinds,
        int n_pos_dim, const vector<int>& affine_residue, const vector<int>& layer,
        const vector<float>& table, int nx, int ny, const vector<int>& rama_residue) {
    return new DevicePlacement(out, alignment, rama, kinds, n_pos_dim, affine_residue, layer, table, nx, ny,
            rama_residue);
}
