// Device implementations of coord_basic.cpp's and coord_map.cpp's nodes
#include "device_nodes.h"

using namespace std;

namespace {

struct DeviceDistance : public DeviceNode {
    DeviceCoord &out, &pos1, &pos2;
    DeviceArray<int> atom1, atom2;
    DeviceArray<float> deriv;  // unit vector from atom 2 to atom 1, 3 per element

    DeviceDistance(DeviceCoord& out_, DeviceCoord& pos1_, DeviceCoord& pos2_,
            const vector<int>& atom1_, const vector<int>& atom2_):
        out(out_), pos1(pos1_), pos2(pos2_), atom1(atom1_), atom2(atom2_), deriv(3*atom1_.size()) {}

    void compute_value(ComputeMode mode) override {
        VecArray x1 = pos1.value(), x2 = pos2.value();
        float* o = out.output.ptr;
        int w = out.row_width;
        const int *a1 = atom1.ptr, *a2 = atom2.ptr;
        float* dv = deriv.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nt) {
            auto disp = load_vec<3>(x1, a1[nt]) - load_vec<3>(x2, a2[nt]);
            auto dist = mag(disp);
            o[nt*w] = dist;
            store_vec(dv+3*nt, dist != 0. ? disp*inv_mag(disp) : disp*0.f);
        });
    }

    void propagate_deriv() override {
        auto s = out.sens.view(), s1 = pos1.sens.view(), s2 = pos2.sens.view();
        int w = out.row_width, w1 = pos1.row_width, w2 = pos2.row_width;
        const int *a1 = atom1.ptr, *a2 = atom2.ptr;
        const float* dv = deriv.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nt) {
            float d = s.get(nt*w);
            for(int k=0; k<3; ++k) {
                s1.add(a1[nt]*w1+k,  dv[3*nt+k]*d);
                s2.add(a2[nt]*w2+k, -dv[3*nt+k]*d);
            }
        });
    }
};

struct DeviceAngle : public DeviceNode {
    DeviceCoord &out, &pos;
    DeviceArray<int> atoms;
    DeviceArray<float> deriv;  // derivative of the cosine by each of the 3 atoms, 9 per element

    DeviceAngle(DeviceCoord& out_, DeviceCoord& pos_, const vector<int>& atoms_):
        out(out_), pos(pos_), atoms(atoms_), deriv(3*atoms_.size()) {}

    void compute_value(ComputeMode mode) override {
        VecArray x = pos.value();
        float* o = out.output.ptr;
        int w = out.row_width;
        const int* a = atoms.ptr;
        float* dv = deriv.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nt) {
            auto atom1 = load_vec<3>(x, a[3*nt+0]);
            auto atom2 = load_vec<3>(x, a[3*nt+1]);
            auto atom3 = load_vec<3>(x, a[3*nt+2]);

            auto x1 = atom1 - atom3; auto inv_d1 = inv_mag(x1); auto x1h = x1*inv_d1;
            auto x2 = atom2 - atom3; auto inv_d2 = inv_mag(x2); auto x2h = x2*inv_d2;
            auto dp = dot(x1h, x2h);
            o[nt*w] = dp;
            auto deriv1 = (x2h - x1h*dp) * inv_d1;
            auto deriv2 = (x1h - x2h*dp) * inv_d2;
            store_vec(dv+9*nt+0, deriv1);
            store_vec(dv+9*nt+3, deriv2);
            store_vec(dv+9*nt+6, -deriv1-deriv2);
        });
    }

    void propagate_deriv() override {
        auto s = out.sens.view(), ps = pos.sens.view();
        int w = out.row_width, wp = pos.row_width;
        const int* a = atoms.ptr;
        const float* dv = deriv.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nt) {
            float d = s.get(nt*w);
            for(int j=0; j<3; ++j)
                for(int k=0; k<3; ++k) ps.add(a[3*nt+j]*wp+k, d*dv[9*nt+3*j+k]);
        });
    }
};

struct DeviceDihedral : public DeviceNode {
    DeviceCoord &out, &pos;
    DeviceArray<int> atoms;
    DeviceArray<float> jac;  // derivative of the angle by each of the 4 atoms, 12 per element

    DeviceDihedral(DeviceCoord& out_, DeviceCoord& pos_, const vector<int>& atoms_):
        out(out_), pos(pos_), atoms(atoms_), jac(3*atoms_.size()) {}

    void compute_value(ComputeMode mode) override {
        VecArray x = pos.value();
        float* o = out.output.ptr;
        int w = out.row_width;
        const int* a = atoms.ptr;
        float* jv = jac.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nt) {
            if(a[4*nt+0] == -1 || a[4*nt+3] == -1) {
                o[nt*w] = -1.3963f;   // -80 degrees if dummy angle; its atom -1 must not be read
                return;
            }
            Vec<3> r[4], d[4];
            for(int na=0; na<4; ++na) r[na] = load_vec<3>(x, a[4*nt+na]);
            o[nt*w] = dihedral_germ(r[0],r[1],r[2],r[3], d[0],d[1],d[2],d[3]);
            for(int na=0; na<4; ++na) store_vec(jv+12*nt+3*na, d[na]);
        });
    }

    void propagate_deriv() override {
        auto s = out.sens.view(), ps = pos.sens.view();
        int w = out.row_width, wp = pos.row_width;
        const int* a = atoms.ptr;
        const float* jv = jac.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nt) {
            if(a[4*nt+0] == -1 || a[4*nt+3] == -1) return;
            float d = s.get(nt*w);
            for(int na=0; na<4; ++na)
                for(int k=0; k<3; ++k) ps.add(a[4*nt+na]*wp+k, d*jv[12*nt+3*na+k]);
        });
    }
};

struct DevicePasteRama : public DeviceNode {
    DeviceCoord &out, &pos1, &pos2;
    DeviceArray<int> index1, index2;

    DevicePasteRama(DeviceCoord& out_, DeviceCoord& pos1_, DeviceCoord& pos2_,
            const vector<int>& index1_, const vector<int>& index2_):
        out(out_), pos1(pos1_), pos2(pos2_), index1(index1_), index2(index2_) {}

    void compute_value(ComputeMode mode) override {
        VecArray x1 = pos1.value(), x2 = pos2.value();
        float* o = out.output.ptr;
        int w = out.row_width;
        const int *i1 = index1.ptr, *i2 = index2.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int ne) {
            o[ne*w+0] = x1(0, i1[ne]);
            o[ne*w+1] = x2(0, i2[ne]);
        });
    }

    void propagate_deriv() override {
        auto s = out.sens.view(), s1 = pos1.sens.view(), s2 = pos2.sens.view();
        int w = out.row_width, w1 = pos1.row_width, w2 = pos2.row_width;
        const int *i1 = index1.ptr, *i2 = index2.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int ne) {
            s1.add(i1[ne]*w1, s.get(ne*w+0));
            s2.add(i2[ne]*w2, s.get(ne*w+1));
        });
    }
};

struct DeviceCatPos : public DeviceNode {
    DeviceCoord &out, &pos1, &pos2;
    int n_dim;
    DeviceArray<int> index, from2;

    DeviceCatPos(DeviceCoord& out_, DeviceCoord& pos1_, DeviceCoord& pos2_, int n_dim_,
            const vector<int>& index_, const vector<int>& from2_):
        out(out_), pos1(pos1_), pos2(pos2_), n_dim(n_dim_), index(index_), from2(from2_) {}

    void compute_value(ComputeMode mode) override {
        VecArray x1 = pos1.value(), x2 = pos2.value();
        float* o = out.output.ptr;
        int w = out.row_width;
        const int *idx = index.ptr, *f2 = from2.ptr;
        int nd = n_dim;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int ne) {
            VecArray x = f2[ne] ? x2 : x1;
            for(int d=0; d<nd; ++d) o[ne*w+d] = x(d, idx[ne]);
        });
    }

    void propagate_deriv() override {
        auto s = out.sens.view(), s1 = pos1.sens.view(), s2 = pos2.sens.view();
        int w = out.row_width, w1 = pos1.row_width, w2 = pos2.row_width;
        const int *idx = index.ptr, *f2 = from2.ptr;
        int nd = n_dim;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int ne) {
            for(int d=0; d<nd; ++d) {
                if(f2[ne]) s2.add(idx[ne]*w2+d, s.get(ne*w+d));
                else       s1.add(idx[ne]*w1+d, s.get(ne*w+d));
            }
        });
    }
};

}

DeviceNode* device_distance(DeviceCoord& out, DeviceCoord& pos1, DeviceCoord& pos2,
        const vector<int>& atom1, const vector<int>& atom2) {
    return new DeviceDistance(out, pos1, pos2, atom1, atom2);
}
DeviceNode* device_angle(DeviceCoord& out, DeviceCoord& pos, const vector<int>& atoms) {
    return new DeviceAngle(out, pos, atoms);
}
DeviceNode* device_dihedral(DeviceCoord& out, DeviceCoord& pos, const vector<int>& atoms) {
    return new DeviceDihedral(out, pos, atoms);
}
DeviceNode* device_paste_rama(DeviceCoord& out, DeviceCoord& pos1, DeviceCoord& pos2,
        const vector<int>& index1, const vector<int>& index2) {
    return new DevicePasteRama(out, pos1, pos2, index1, index2);
}
DeviceNode* device_cat_pos(DeviceCoord& out, DeviceCoord& pos1, DeviceCoord& pos2, int n_dim,
        const vector<int>& index, const vector<int>& from2) {
    return new DeviceCatPos(out, pos1, pos2, n_dim, index, from2);
}
