#ifndef EIG4_H
#define EIG4_H

// The symmetric 4x4 eigensolver of eig.cpp's affine_alignment (Householder tridiagonalisation, then implicit QR with
// Givens rotations; Golub and van Loan, Matrix Computations, chapters 5 and 8), written in the lane type S: Float4
// on the CPU, four residues at once, and float on the device.

#include "vector_math.h"

// internal linkage, as in the .cpp file it came from, so the CPU code generation is unchanged
namespace {


template <typename S>
UPSIDE_HD S house(
        int n,
        S* x)  // size (n,), replaced by v, v[0] = |x| (x[0] itself when the column is already reduced)
{
    S beta,mu,s;
    S sigma2 = zero<S>();

    for(int i=1; i<n; ++i) sigma2 += x[i]*x[i];

    // a column with no tail is already reduced, so its reflection is the identity (Golub and van Loan 5.1.1)
    auto reduced = sigma2 == zero<S>();

    mu = sqrtf(x[0]*x[0] + sigma2);   // |x|
    s = ternary(zero<S>()<x[0], -sigma2*rcp(x[0]+mu), x[0]-mu);
    beta = ternary(reduced, zero<S>(), two<S>()*s*s * rcp(sigma2+s*s));

    x[0] = ternary(reduced, x[0], mu);
    auto inv_s = ternary(reduced, zero<S>(), rcp(s));
    for(int i=1; i<n; ++i) x[i] *= inv_s;

    return beta;
}

template <typename S>
UPSIDE_HD void givens(S *c, S *s, S a, S b)
{
    S inv_r = rsqrt(a*a+b*b);

    auto trivial = b==zero<S>();
    *c = ternary(trivial, one <S>(),  a*inv_r);
    *s = ternary(trivial, zero<S>(), -b*inv_r);
}



UPSIDE_HD inline int r(int i, int j) {
    int ii = (i<j) ? i : j;
    int jj = (i<j) ? j : i;

    switch(ii) {
        case 0:
            return (0-0) + jj;
        case 1:
            return (4-1) + jj;
        case 2:
            return (7-2) + jj;
        case 3:
            return (9-3) + jj;
        default:
            return 1000; // error
    }
}


template <typename S>
UPSIDE_HD void
symmetric_tridiagonalize_4x4(
        S* restrict beta, // length 2 beta sequence of Householder sequence
        S* restrict A)   // upper triangle of 4x4 matrix, 10 elements
{
    S p[3];
    S w[3];

    for(int k=0; k<4-2; ++k) {
        int m = 4-(k+1);
        beta[k] = house(m, A+r(k,k+1));   // overwrite with householder

        #define v(j) ((j)==0 ? one<S>() : A[r(k,k+1+(j))])
        for(int i=0; i<m; ++i) {
            p[i] = zero<S>();
            for(int j=0; j<m; ++j)
                p[i] += A[r(i+k+1,j+k+1)] * v(j);
            p[i] *= beta[k];
        }

        S p_dot_v = zero<S>();
        for(int i=0; i<m; ++i) p_dot_v += p[i]*v(i);
        for(int i=0; i<m; ++i) w[i] = p[i] - (half<S>()*beta[k]*p_dot_v) * v(i);

        for(int i=0; i<m; ++i)
            for(int j=i; j<m; ++j)
                A[r(i+k+1,j+k+1)] -= v(i)*w[j] + v(j)*w[i];
        #undef v
    }
}


template <typename S>
UPSIDE_HD void
unpack_tridiagonalize_4x4(
        S * restrict d,  // length 4
        S * restrict u,  // length 3
        S * restrict rot_, // size 4x4 = 16, nonsymmetric
        S * restrict beta, // length 2
        S * restrict A)    // symmetric_4x4, length 10
{
    #define rot(i,j) (rot_[(i)*4+(j)])
    for(int i=0; i<4; ++i) d[i] = A[r(i,i)];
    for(int i=0; i<3; ++i) u[i] = A[r(i,i+1)];

    // first update with the symmetric part
    for(int i=0; i<4; ++i) for(int j=i; j<4; ++j) rot(i,j) = (i==j);   // identity matrix to start

    // unpack Householder reflection vectors
    #define v0(j) ((j)==0 ? one<S>() : A[r(0,1+(j))])
    #define v1(j) ((j)==0 ? one<S>() : A[r(1,2+(j))])

    // first Householder reflection operates on lower left 3x3 block
    for(int i=1; i<4; ++i)
        for(int j=i; j<4; ++j)
            rot(i,j) -= beta[0] * (v0(i-1)*v0(j-1));

    // second Householder reflection operates on the lowest 2x2 block
    for(int i=2; i<4; ++i)
        for(int j=i; j<4; ++j)
            rot(i,j) -= beta[1] * (v1(i-2)*v1(j-2));

    // fill in using symmetry
    for(int i=0; i<4; ++i) for(int j=i+1; j<4; ++j) rot(j,i) = rot(i,j);

    // now fill in non-symmetric Householder interaction (note leading zero for v1)
    // fills in lower left 2x3 block
    S coeff = beta[0]*beta[1] * (v0(1)*v1(0) + v0(2)*v1(1));
    for(int i=2; i<4; ++i)
        for(int j=1; j<4; ++j)
            rot(i,j) += coeff * v1(i-2)*v0(j-1);
    #undef rot
    #undef v0
    #undef v1
}


#define rot(i,j) (rot_[(i)*4 + (j)])

template <typename S>
UPSIDE_HD void
implicit_symm_QR_step_4x4(
        int n,   // length of diagonal
        S * restrict d,    // diagonal
        S * restrict u,   // upper superdiagonal
        S * restrict rot_) // rotation to accumulate Givens rotations
                     // must be at least nx4-sized
{
    S dval = half<S>()*(d[n-2] - d[n-1]) + S(1e-20f);  // paranoia against NaN
    S mu = d[n-1] - u[n-2]*u[n-2]*rcp(dval+copysignf(sqrtf(dval*dval + u[n-2]*u[n-2]),dval));

    S x = d[0] - mu;
    S z = u[0];
    
    for(int k=0; k<n-1; ++k){
        S c,s;
        givens(&c,&s,x,z);

        /* accumulate Givens rotation */
        for(int j=0; j<4; ++j) {
            S t1 = rot(k  ,j);
            S t2 = rot(k+1,j);
            rot(k  ,j) = c*t1 - s*t2;
            rot(k+1,j) = s*t1 + c*t2;
        }

        /* compute the T(k-1,:) row */
        if(k > 0) u[k-1] = c*x - s*z;

        /* mix the T(k,k),T(k,k+1),T(k+1,k+1) block */
        S T00 = d[k  ];
        S T11 = d[k+1];
        S T01 = u[k  ];
        
        d[k  ]  = T00*c*c - T01*two<S>()*c*s + T11*s*s;
        d[k+1]  = T00*s*s + T01*two<S>()*c*s + T11*c*c;
        u[k  ]  = (T00-T11)*c*s + T01*(c*c-s*s);
        x = u[k];

        /* compute the T(k+2,:) row */
        if(k<n-2) {
            z = -u[k+1]*s;
            u[k+1] *= c;
        }
    } 
}
#undef rot


template <typename S>
UPSIDE_HD int  __attribute__ ((noinline))
symm_QR_4x4(
        S* restrict d, // will contain eigenvalues, length 4
        S* restrict rot, // will contain eigenvectors in rows, size 4x4
        S* restrict A, // symmetric 4x4 matrix, 10 elements, overwritten
        float tol_,  // relative tolerance for off-diagonal entries (suggest something like 1e-5)
        int max_iter) 
{
    const int n = 4;
    S tol = S(tol_);

    S beta[2];
    symmetric_tridiagonalize_4x4(beta, A);

    S u[3];
    unpack_tridiagonalize_4x4(d, u, rot, beta, A);

    for(int k=0; k<max_iter; ++k) {
        // exactly zero off-diagonal elements that are nearly zero
        for(int i=0; i<n-1; ++i) {
            auto small_u = fabsf(u[i]) <= tol * (fabsf(d[i]) + fabsf(d[i+1]));
            u[i] = ternary(small_u, zero<S>(), u[i]);
        }

        // find largest diagonal lower left subblock
        int q;
        if      (lane_any(u[2] != zero<S>())) q = 0;
        else if (lane_any(u[1] != zero<S>())) q = 1;
        else if (lane_any(u[0] != zero<S>())) q = 2;
        else return k;   // q=3 and the matrix is diagonal and we are done

        // p is start of largest unreduced triangular submatrix adjacent to the q region
        // unreduced means no entries in the upper triangle are zero
        int p;
        for(p=n-q-1; p>0; --p) {
            if(lane_none(u[p-1] != zero<S>())) break;
        }

        // perform QR iteration and accumulate into rot
        implicit_symm_QR_step_4x4(n-q-p, d+p, u+p, rot+4*p);
    }
    return -1;  // non-convergence
}
}

#endif
