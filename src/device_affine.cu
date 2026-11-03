// Device implementations of eig.cpp's affine_alignment and backbone_steric.cpp's backbone_pairs
#include "device_interaction_graph.h"
#include "eig4.h"
#include "affine.h"
#include "backbone_steric.h"

using namespace std;

namespace {

// The derivative of one residue's frame by its 3 atoms, as AffineAlignment::propagate_deriv computes it for each
// lane: sens3 is the sensitivity to the frame's origin, torque to its rotation, g the atom's ref_geom row; the
// derivative of the quaternion, the largest eigenvector, is the first-order perturbation by the other eigenvectors
// (the largest eigenvalue is taken to be nondegenerate, as for a reasonable alignment to a nearly rigid structure)
UPSIDE_HD inline Vec<3> affine_atom_deriv(const float* evals, const float* evecs, const float* sens3,
        const float* torque, const float* g)
{
    float inv_evals[4];
    for(int j=1; j<4; ++j) inv_evals[j] = rcp(evals[0]-evals[j]);

    // push back torque to affine derivatives (multiply by quaternion); the torque is in the tangent space of the
    // rotated frame, so it is pushed into the rotated space by a right multiply by the quaternion itself
    float quat_sens[4] = {
        2.f*(-torque[0]*evecs[1] - torque[1]*evecs[2] - torque[2]*evecs[3]),
        2.f*( torque[0]*evecs[0] + torque[1]*evecs[3] - torque[2]*evecs[2]),
        2.f*( torque[1]*evecs[0] + torque[2]*evecs[1] - torque[0]*evecs[3]),
        2.f*( torque[2]*evecs[0] + torque[0]*evecs[2] - torque[1]*evecs[1])};

    float quat_sens_diag_basis[4] = {0.f,0.f,0.f,0.f};
    for(int d=0; d<4; ++d) for(int i=0; i<4; ++i) quat_sens_diag_basis[d] += quat_sens[i]*evecs[d*4+i];

    // the perturbation term of the matrix element f (upper triangle, eig4.h's r)
    auto perturb = [&](float& ret, const float f[10]) {
        for(int k=1; k<4; ++k) {
            float c = 0.f;
            for(int i=0; i<4; ++i)
                for(int j=i; j<4; ++j)
                    c += f[r(i,j)] * ((i==j) ? evecs[k*4+i]*evecs[0*4+j]
                                             : evecs[k*4+i]*evecs[0*4+j]+evecs[k*4+j]*evecs[0*4+i]);
            ret += inv_evals[k]*c*quat_sens_diag_basis[k];
        }
    };

    Vec<3> deriv = (1.f/3.f)*make_vec3(sens3[0],sens3[1],sens3[2]);
    const float f0[10] = {g[0],  0.f, g[2],-g[1],
                                 g[0],g[1], g[2],
                                      -g[0],  0.f,
                                            -g[0]};
    const float f1[10] = {g[1],-g[2],  0.f, g[0],
                                -g[1], g[0],  0.f,
                                       g[1], g[2],
                                            -g[1]};
    const float f2[10] = {g[2], g[1],-g[0],  0.f,
                                -g[2],  0.f, g[0],
                                      -g[2], g[1],
                                             g[2]};
    perturb(deriv[0], f0);
    perturb(deriv[1], f1);
    perturb(deriv[2], f2);
    return deriv;
}

struct DeviceAffineAlignment : public DeviceNode {
    DeviceCoord &out, &pos;
    DeviceArray<int> atom;
    DeviceArray<float> ref_geom;
    DeviceArray<float> evals, evecs;  // 4 and 16 per residue, for the derivative

    DeviceAffineAlignment(DeviceCoord& out_, DeviceCoord& pos_, const vector<int>& atom_, const vector<float>& ref_geom_):
        out(out_), pos(pos_), atom(atom_), ref_geom(ref_geom_), evals(4*out_.n_elem), evecs(16*out_.n_elem) {}

    void compute_value(ComputeMode mode) override {
        VecArray x = pos.value();
        float* o = out.output.ptr;
        int w = out.row_width;
        const int* a = atom.ptr;
        const float* rg = ref_geom.ptr;
        float *ev = evals.ptr, *ec = evecs.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nr) {
            auto atom1 = load_vec<3>(x, a[3*nr+0]);
            auto atom2 = load_vec<3>(x, a[3*nr+1]);
            auto atom3 = load_vec<3>(x, a[3*nr+2]);

            auto center = (1.f/3.f)*(atom1+atom2+atom3);
            atom1 -= center;
            atom2 -= center;
            atom3 -= center;

            const float* g = rg + 9*nr;  // (atom, xyz)
            float R_[3][3];
            #define R(i,j) (R_[i][j])
            for(int i=0; i<3; ++i)
                for(int j=0; j<3; ++j)
                    R(i,j) = atom1[j] * g[0*3+i]
                           + atom2[j] * g[1*3+i]
                           + atom3[j] * g[2*3+i];

            float F[10] = {R(0,0)+R(1,1)+R(2,2), R(1,2)-R(2,1),         R(2,0)-R(0,2),         R(0,1)-R(1,0),
                                                 R(0,0)-R(1,1)-R(2,2),  R(0,1)+R(1,0),         R(0,2)+R(2,0),
                                                                       -R(0,0)+R(1,1)-R(2,2),  R(1,2)+R(2,1),
                                                                                              -R(0,0)-R(1,1)+R(2,2)};
            #undef R

            float* e = ev + 4*nr;
            float* v = ec + 16*nr;
            symm_QR_4x4(e, v, F, 2e-7f, 100);  // 1e-5 left frame errors that switched on and off: energy jumps of 1e-3

            // swap largest eigenvalue into location 0
            for(int i=1; i<4; ++i) {
                if(e[0] < e[i]) {
                    float t = e[0]; e[0] = e[i]; e[i] = t;
                    for(int d=0; d<4; ++d) {float u = v[d]; v[d] = v[i*4+d]; v[i*4+d] = u;}
                }
            }

            for(int j=0; j<3; ++j) o[nr*w+j]   = center[j];
            for(int j=0; j<4; ++j) o[nr*w+3+j] = v[j];
        });
    }

    void propagate_deriv() override {
        auto s = out.sens.view(), ps = pos.sens.view();
        int w = out.row_width, wp = pos.row_width;
        const int* a = atom.ptr;
        const float* rg = ref_geom.ptr;
        const float *ev = evals.ptr, *ec = evecs.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int nr) {
            float sens3 [3] = {s.get(nr*w+0), s.get(nr*w+1), s.get(nr*w+2)};
            float torque[3] = {s.get(nr*w+3), s.get(nr*w+4), s.get(nr*w+5)};
            for(int na=0; na<3; ++na) {
                auto d = affine_atom_deriv(ev+4*nr, ec+16*nr, sens3, torque, rg+9*nr+3*na);
                for(int k=0; k<3; ++k) ps.add(a[3*nr+na]*wp+k, d[k]);
            }
        });
    }
};


struct BackboneAccept {UPSIDE_HD static bool ok(int a, int b) {return acceptable_backbone_pair<int>(a, b);}};

struct DeviceBackbonePairs : public DeviceNode {
    AccumView potential;
    DeviceCoord& alignment;
    int n_residue;
    DeviceArray<int> residue, n_atom, system;
    DeviceArray<float> ref_pos;          // 12 per residue
    DeviceArray<float> coords, ref_pos_coords;  // each residue's frame origin (3), and its atoms placed (12)
    DevicePairlist<true> pairlist;

    DeviceBackbonePairs(DeviceGraph& graph, AccumView potential_, DeviceCoord& alignment_, const vector<int>& residue_,
            const vector<int>& n_atom_, const vector<float>& ref_pos_, const vector<int>& sequence_id,
            const vector<int>& n_residue_of_system, float dist_cutoff, float cache_buffer,
            const vector<int>& max_n_edge):
        potential(potential_), alignment(alignment_), n_residue(residue_.size()), residue(residue_), n_atom(n_atom_),
        system(DevicePairlist<true>::systems(n_residue_of_system)), ref_pos(ref_pos_),
        coords(3*residue_.size()), ref_pos_coords(12*residue_.size()),
        pairlist(graph, "backbone_pairs", n_residue_of_system, n_residue_of_system, sequence_id, sequence_id,
                dist_cutoff, cache_buffer, max_n_edge) {}

    void compute_value(ComputeMode mode) override {
        VecArray al = alignment.value();
        const int* res = residue.ptr;
        const float* rp = ref_pos.ptr;
        float *cd = coords.ptr, *rc = ref_pos_coords.ptr;
        parallel_for(n_residue, DEVICE_LAMBDA (int nr) {
            auto aff = load_vec<7>(al, res[nr]);
            float U[9]; quat_to_rot(U, aff.v+3);
            auto t = extract<0,3>(aff);
            store_vec(cd+3*nr, t);
            for(int na=0; na<4; ++na)
                store_vec(rc+12*nr+3*na, apply_affine(U,t, load_vec<3>(rp+12*nr+3*na)));
        });

        // acceptable_backbone_pair checks that the sequence ids differ by at least 2
        pairlist.find_edges<BackboneAccept>(cd, 3, cd, 3);

        auto as = alignment.sens.view(), pot = potential;
        bool want_pot = mode == PotentialAndDerivMode;
        int wa = alignment.row_width, cap = pairlist.max_n_edge;
        const int *e1 = pairlist.edge1.ptr, *e2 = pairlist.edge2.ptr, *ne = pairlist.n_edge.ptr;
        const int *na_ = n_atom.ptr, *sys = system.ptr;
        slots_for(pairlist.n_system, cap, ne, DEVICE_LAMBDA (int e) {
            int nr1 = e1[e], nr2 = e2[e];

            auto d1 = make_zero<3>(); auto torque1 = make_zero<3>(); auto t1 = load_vec<3>(cd+3*nr1);
            auto d2 = make_zero<3>(); auto torque2 = make_zero<3>(); auto t2 = load_vec<3>(cd+3*nr2);

            bool hit = false;
            float pot_sum = 0.f;
            for(int i1=0; i1<na_[nr1]; ++i1) {
                const Vec<3> x1 = load_vec<3>(rc+12*nr1+3*i1);

                for(int i2=0; i2<na_[nr2]; ++i2) {
                    const Vec<3> x2 = load_vec<3>(rc+12*nr2+3*i2);

                    const Vec<3> r = x1-x2;
                    const float r_mag2 = mag2(r);
                    if(r_mag2>nonbonded_atom_cutoff2) continue;
                    hit = true;
                    const float deriv_over_r  = nonbonded_kernel_or_deriv_over_r<true> (r_mag2);
                    pot_sum += nonbonded_kernel_or_deriv_over_r<false>(r_mag2);
                    const Vec<3> g = deriv_over_r*r;

                    d1 +=  g;  torque1 += cross(x1-t1,  g);
                    d2 += -g;  torque2 += cross(x2-t2, -g);
                }
            }

            if(hit) {
                if(want_pot) pot.add(sys[nr1], pot_sum);
                for(int d=0; d<3; ++d) {
                    as.add(res[nr1]*wa+d,   d1[d]);  as.add(res[nr1]*wa+3+d, torque1[d]);
                    as.add(res[nr2]*wa+d,   d2[d]);  as.add(res[nr2]*wa+3+d, torque2[d]);
                }
            }
        });
    }
};

}

DeviceNode* device_affine_alignment(DeviceCoord& out, DeviceCoord& pos, const vector<int>& atom,
        const vector<float>& ref_geom) {
    return new DeviceAffineAlignment(out, pos, atom, ref_geom);
}

DeviceNode* device_backbone_pairs(DeviceGraph& graph, AccumView potential, DeviceCoord& alignment,
        const vector<int>& residue, const vector<int>& n_atom, const vector<float>& ref_pos,
        const vector<int>& sequence_id, const vector<int>& n_residue, float dist_cutoff,
        float cache_buffer, const vector<int>& max_n_edge) {
    return new DeviceBackbonePairs(graph, potential, alignment, residue, n_atom, ref_pos, sequence_id, n_residue,
            dist_cutoff, cache_buffer, max_n_edge);
}
