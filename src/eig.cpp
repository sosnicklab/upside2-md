// Author: John Jumper
//
// RMSD alignment algorithm based on Coutsias, Seok, and Dill, Journal of
// Computational Chemistry, Vol. 25, No. 15
//
// Linear algebra based on algorithms in Golub and van Loan, Matrix
// Computations, chapters 5 and 8

#include "deriv_engine.h"
#include "timing.h"
#include "affine.h"
#include "eig4.h"
#include "device_nodes.h"

using namespace h5;
using namespace std;

typedef Float4 S;

struct AffineAlignment : public CoordNode
{
    struct Params {
        alignas(16) int32_t atom_offsets[3][4];
        alignas(16) float   ref_geom[3][3][4];
    };

    int n_group;
    CoordNode& pos;
    vector<Params> params;
    unique_ptr<Float4[]> evals_storage;
    unique_ptr<Float4[]> evecs_storage;

    AffineAlignment(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(2, grp, "atoms")[0], 7),
        n_group(round_up(n_elem,4)/4),
        
        pos(pos_), params(n_group),
        evals_storage(new_aligned<Float4>(n_group* 4)),
        evecs_storage(new_aligned<Float4>(n_group*16))
    {
        check_size(grp, "atoms",    n_elem, 3);
        check_size(grp, "ref_geom", n_elem, 3,3);  // (residue, atom, xyz)

        traverse_dset<2,int  >(grp,"atoms",   [&](size_t i,size_t j,          int   x){
                if(x<0 || x>=pos.n_elem) throw string("affine_alignment atoms outside its input");
                params[i/4].atom_offsets[j][i%4] = x*pos.output.row_width;});
        traverse_dset<3,float>(grp,"ref_geom",[&](size_t i,size_t na,size_t d,float x){
                params[i/4].ref_geom[na][d][i%4]=x;});

        // handle non-multiple of 4's by padding (needs testing)
        for(int i=n_elem; i<round_up(n_elem,4); ++i) {
            for(int j: range(3))
                params[i/4].atom_offsets[j][i%4] = params[i/4].atom_offsets[j][0];

            for(int na: range(3))
                for(int d: range(3))
                    params[i/4].ref_geom[na][d][i%4] = params[i/4].ref_geom[na][d][0];
        }

        // the frame's origin is the centroid of its atoms, and the derivative below takes the covariance's
        // derivative with respect to an atom as its ref_geom row, which holds only for a centred ref_geom
        for(int i: range(n_elem))
            for(int d: range(3))
                if(fabsf(params[i/4].ref_geom[0][d][i%4] + params[i/4].ref_geom[1][d][i%4] +
                         params[i/4].ref_geom[2][d][i%4]) > 3e-4f)
                    throw string("ref_geom of residue ") + to_string(i) + " is not centred";
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("affine_alignment"));

        VecArray rigid_body = output;
        float* posc = pos.output.x.get();

        for(int ng=0; ng<n_group; ++ng) {
            const auto& p = params[ng];

            auto atom1 = aligned_gather_vec<3>(posc, Int4(p.atom_offsets[0]));
            auto atom2 = aligned_gather_vec<3>(posc, Int4(p.atom_offsets[1]));
            auto atom3 = aligned_gather_vec<3>(posc, Int4(p.atom_offsets[2]));

            auto center = S(1.f/3.f)*(atom1+atom2+atom3);
            atom1 -= center;
            atom2 -= center;
            atom3 -= center;

            auto ref_geom1 = make_vec3(Float4(p.ref_geom[0][0]), Float4(p.ref_geom[0][1]), Float4(p.ref_geom[0][2]));
            auto ref_geom2 = make_vec3(Float4(p.ref_geom[1][0]), Float4(p.ref_geom[1][1]), Float4(p.ref_geom[1][2]));
            auto ref_geom3 = make_vec3(Float4(p.ref_geom[2][0]), Float4(p.ref_geom[2][1]), Float4(p.ref_geom[2][2]));

            S R_[3][3];
            #define R(i,j) (R_[i][j])
            for(int i=0; i<3; ++i)
                for(int j=0; j<3; ++j)
                    R(i,j) = atom1[j] * ref_geom1[i]
                           + atom2[j] * ref_geom2[i]
                           + atom3[j] * ref_geom3[i];

            S F[10] = {R(0,0)+R(1,1)+R(2,2), R(1,2)-R(2,1),         R(2,0)-R(0,2),         R(0,1)-R(1,0),
                                             R(0,0)-R(1,1)-R(2,2),  R(0,1)+R(1,0),         R(0,2)+R(2,0),
                                                                   -R(0,0)+R(1,1)-R(2,2),  R(1,2)+R(2,1),
                                                                                          -R(0,0)-R(1,1)+R(2,2)};
            #undef R

            // S evals[4], evecs[16];
            Float4* restrict evals = evals_storage.get() + ng* 4;
            Float4* restrict evecs = evecs_storage.get() + ng*16;

            symm_QR_4x4(evals, evecs, F, 2e-7f, 100);  // 1e-5 left frame errors that switched on and off: energy jumps of 1e-3

            // swap largest eigenvalue into location 0
            for(int i=1; i<4; ++i) {
                auto do_flip = evals[0] < evals[i];

                auto eval0 = ternary(do_flip, evals[i], evals[0]);
                auto evali = ternary(do_flip, evals[0], evals[i]);
                evals[0] = eval0;
                evals[i] = evali;

                for(int d=0; d<4; ++d) {
                    auto evec0 = ternary(do_flip, evecs[i*4+d], evecs[0*4+d]);
                    auto eveci = ternary(do_flip, evecs[0*4+d], evecs[i*4+d]);
                    evecs[0*4+d] = evec0;
                    evecs[i*4+d] = eveci;
                }
            }

            Vec<8,S> body; // really 7 components but I need the eight for the transpose
            for(int j=0; j<3; ++j) body[j] = center[j];
            for(int j=0; j<4; ++j) body[3+j] = evecs[0*4+j];

            transpose4(body[0],body[1],body[2],body[3]);
            for(int i=0; i<4; ++i) body[i].store(&rigid_body(0,4*ng+i));

            transpose4(body[4],body[5],body[6],body[7]);
            for(int i=0; i<4; ++i) body[4+i].store(&rigid_body(4,4*ng+i));
        }
    }

    virtual void propagate_deriv() {
        Timer timer(string("affine_alignment_deriv"));
        float* pos_sens = pos.sens.x.get();

        for(int ng=0; ng<n_group; ++ng) {
            const auto& p = params[ng];

            const Float4* restrict evals = evals_storage.get() + ng* 4;
            const Float4* restrict evecs = evecs_storage.get() + ng*16;

            // compute inverse eigenvalue differences
            S inv_evals[4];
            for(int j=1; j<4; ++j) inv_evals[j] = rcp(evals[0]-evals[j]);

            // push back torque to affine derivatives (multiply by quaternion)
            // the torque is in the tangent space of the rotated frame
            // to act on a tangent in the affine space, I need to push that tangent into the rotated space
            // this means a right multiply by the quaternion itself

            // evecs[0] is the rotation quaternion
            S sens3  [3] = {Float4(&sens(0,4*ng+0)), Float4(&sens(0,4*ng+1)),Float4(&sens(0,4*ng+2))};
            S torque [3] = {Float4(&sens(0,4*ng+3)), Float4(&sens(4,4*ng+0)),Float4(&sens(4,4*ng+1))};
            S padding[2] = {Float4(&sens(4,4*ng+2)), Float4(&sens(4,4*ng+3))};

            transpose4(sens3 [0],sens3 [1],sens3  [2],torque [0]);
            transpose4(torque[1],torque[2],padding[0],padding[1]);

            S quat_sens[4] = {
                two<S>()*(-torque[0]*evecs[1] - torque[1]*evecs[2] - torque[2]*evecs[3]),
                two<S>()*( torque[0]*evecs[0] + torque[1]*evecs[3] - torque[2]*evecs[2]),
                two<S>()*( torque[1]*evecs[0] + torque[2]*evecs[1] - torque[0]*evecs[3]),
                two<S>()*( torque[2]*evecs[0] + torque[0]*evecs[2] - torque[1]*evecs[1])};

            S quat_sens_diag_basis[4] = {zero<S>(),zero<S>(),zero<S>(),zero<S>()};
            for(int d=0; d<4; ++d) for(int i=0; i<4; ++i) quat_sens_diag_basis[d] += quat_sens[i]*evecs[d*4+i];

            // The derivative of an eigenvector is given by perturbation theory as 
            // a linear combination of the other eigenvectors weighted by interaction
            // matrix elements and the difference of eigenvalues

            // I will assume that the largest eigenvalue is nondegenerate, as would be expected
            // for a reasonable alignment to a nearly rigid structure

            // unsafe macro for perturbation term
            #define t(i,j) (f[r(i,j)] * ((i==j) \
                                            ? evecs[k*4+i]*evecs[0*4+j] \
                                            : evecs[k*4+i]*evecs[0*4+j]+evecs[k*4+j]*evecs[0*4+i]))
            #define perturb(ret, f00,f01,f02,f03,f11,f12,f13,f22,f23,f33) do { \
                const S f[10] = {f00,f01,f02,f03,f11,f12,f13,f22,f23,f33}; \
                for(int k=1; k<4; ++k) { \
                    S c = inv_evals[k]*(t(0,0)+t(0,1)+t(0,2)+t(0,3) \
                                              +t(1,1)+t(1,2)+t(1,3) \
                                                     +t(2,2)+t(2,3) \
                                                            +t(3,3)); \
                    (ret) += c*quat_sens_diag_basis[k]; \
                } \
            } while(0)


            for(int na=0; na<3; ++na) {
                Vec<3,S> deriv = S(1.f/3.f)*make_vec3(sens3[0],sens3[1],sens3[2]);
                auto g = make_vec3(Float4(p.ref_geom[na][0]), Float4(p.ref_geom[na][1]), Float4(p.ref_geom[na][2]));

                // quaternion rotation derivative
                perturb(deriv[0], g[0],  zero<S>(), g[2],-g[1], 
                                         g[0],      g[1], g[2], 
                                                   -g[0],  zero<S>(), 
                                                          -g[0]);

                perturb(deriv[1], g[1],-g[2],  zero<S>(), g[0], 
                                       -g[1],  g[0],      zero<S>(), 
                                               g[1],      g[2], 
                                                         -g[1]);

                perturb(deriv[2], g[2], g[1],  -g[0],      zero<S>(), 
                                       -g[2],   zero<S>(), g[0], 
                                               -g[2],      g[1], 
                                                           g[2]);

                aligned_scatter_update_vec_destructive(pos_sens, Int4(p.atom_offsets[na]), deriv);
            }
        }
    }

    virtual DeviceNode* make_device(DeviceGraph& g, const vector<DerivComputation*>& systems) {
        vector<int> atoms;
        vector<float> ref;
        auto& c = g.coord(pos);
        for(int ns: range(g.n_system)) {
            auto& n = static_cast<AffineAlignment&>(*systems[ns]);
            for(int i: range(n.n_elem)) {
                auto& p = n.params[i/4];
                for(int j: range(3)) atoms.push_back(p.atom_offsets[j][i%4]/n.pos.output.row_width + c.offset[ns]);
                for(int na: range(3)) for(int d: range(3)) ref.push_back(p.ref_geom[na][d][i%4]);
            }
        }
        return device_affine_alignment(g.coord(*this), c, atoms, ref);
    }
};

static RegisterNodeType<AffineAlignment,1>_8("affine_alignment");
