#include <math.h>
#include <string>
#include <algorithm>
#include <iostream>
#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include "vector_math.h"

using namespace h5;
using namespace std;

struct DistCoord : public CoordNode
{
    struct Params { index_t atom[2]; };

    int n_elem;
    CoordNode& pos1;
    CoordNode& pos2;
    vector<Params> params;
    vector<float3> deriv;

    DistCoord(hid_t grp, CoordNode& pos1_,  CoordNode& pos2_):
        CoordNode(get_dset_size(2, grp, "id")[0], 1),
        n_elem(get_dset_size(2, grp, "id")[0]), 
	pos1(pos1_), 
	pos2(pos2_), 
	params(n_elem),
	deriv(n_elem) 
    {
        int n_dep = 2;  // number of atoms that each term depends on 
        check_size(grp, "id", n_elem, n_dep);
        traverse_dset<2,int> (grp, "id", [&](size_t i, size_t j, int x) {params[i].atom[j] = x;});
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("distance"));
        VecArray posc1 = pos1.output;
        VecArray posc2 = pos2.output;
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            auto x1 = load_vec<3>(posc1, p.atom[0]);
            auto x2 = load_vec<3>(posc2, p.atom[1]);
            auto disp = x1 - x2;
	    auto dist = mag(disp);
	    output(0, nt) = dist;
	    if (dist != 0.)
                deriv[nt] = disp * inv_mag(disp);
	    else
                deriv[nt] = disp * 0.f;
        }
    }

    virtual void propagate_deriv() {
        Timer timer(string("distance_deriv"));
        VecArray pos_sens1 = pos1.sens;
        VecArray pos_sens2 = pos2.sens;
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            auto  d = sens(0, nt);
            update_vec(pos_sens1, p.atom[0],  deriv[nt]*d);
            update_vec(pos_sens2, p.atom[1], -deriv[nt]*d);
        }
    }
};
static RegisterNodeType<DistCoord,2> dist_coord_node("Distance3D");

struct Dist2DCoord : public CoordNode
{
    struct Params { index_t atom[2]; };

    int n_elem;
    int dim1;
    int dim2;
    CoordNode& pos1;
    CoordNode& pos2;
    vector<Params> params;
    vector<float2> deriv;

    Dist2DCoord(hid_t grp, CoordNode& pos1_,  CoordNode& pos2_):
        CoordNode(get_dset_size(2, grp, "id")[0], 1),
        n_elem(get_dset_size(2, grp, "id")[0]), 
        dim1(read_attribute<int>(grp, ".", "dim1")),
        dim2(read_attribute<int>(grp, ".", "dim2")),
	pos1(pos1_), 
	pos2(pos2_), 
	params(n_elem),
	deriv(n_elem) 
    {
        int n_dep = 2;  // number of atoms that each term depends on 
        check_size(grp, "id", n_elem, n_dep);
        traverse_dset<2,int> (grp, "id", [&](size_t i, size_t j, int x) {params[i].atom[j] = x;});
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("distance2d"));
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p       = params[nt];
            auto x1       = pos1.output(dim1, p.atom[0]);
            auto x2       = pos2.output(dim1, p.atom[1]);
            auto y1       = pos1.output(dim2, p.atom[0]);
            auto y2       = pos2.output(dim2, p.atom[1]);
	    auto disp     = make_vec2(x1-x2, y1-y2);
	    auto dist     = mag(disp);
	    output(0, nt) = dist;
	    if (dist != 0.)
                deriv[nt] = disp * inv_mag(disp);
	    else
                deriv[nt] = disp * 0.f;
        }
    }

    virtual void propagate_deriv() {
        Timer timer(string("distance2d_deriv"));
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            auto  d = sens(0, nt);
            pos1.sens(dim1, p.atom[0]) += deriv[nt][0]*d;
            pos2.sens(dim1, p.atom[1]) -= deriv[nt][0]*d;
            pos1.sens(dim2, p.atom[0]) += deriv[nt][1]*d;
            pos2.sens(dim2, p.atom[1]) -= deriv[nt][1]*d;

        }
    }
};
static RegisterNodeType<Dist2DCoord,2> dist2D_coord_node("Distance2D");

struct Dist1DCoord : public CoordNode
{
    struct Params { index_t atom[2]; };

    int n_elem;
    int dim1;
    CoordNode& pos1;
    CoordNode& pos2;
    vector<float> deriv;
    vector<Params> params;

    Dist1DCoord(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
        CoordNode(get_dset_size(2, grp, "id")[0], 1),
        n_elem(get_dset_size(2, grp, "id")[0]), 
        dim1(read_attribute<int>(grp, ".", "dim1")),
	pos1(pos1_), 
	pos2(pos2_), 
	deriv(n_elem),
	params(n_elem)
    {
        int n_dep = 2;  // number of atoms that each term depends on 
        check_size(grp, "id", n_elem, n_dep);
        traverse_dset<2,int> (grp, "id", [&](size_t i, size_t j, int x) {params[i].atom[j] = x;});
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("distance1d"));
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p       = params[nt];
            auto x1       = pos1.output(dim1, p.atom[0]);
            auto x2       = pos2.output(dim1, p.atom[1]);
	    auto disp     = x1-x2;
	    auto dist     = sqrtf(disp*disp);
	    output(0, nt) = dist;
	    if (dist == 0)
	        deriv[nt] = 0.0;
	    else
	        deriv[nt] = disp/dist;
        }
    }

    virtual void propagate_deriv() {
        Timer timer(string("distance1d_deriv"));
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            auto  d = sens(0, nt);
            pos1.sens(dim1, p.atom[0]) += d*deriv[nt];
            pos2.sens(dim1, p.atom[1]) -= d*deriv[nt];

        }
    }
};
static RegisterNodeType<Dist1DCoord,2> dist1D_coord_node("Distance1D");


struct AngleCoord : public CoordNode
{
    //struct alignas(16) Jac {float d[3][4];}; // padding for vector load/store
    struct Params { index_t atom[3];};

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    //unique_ptr<Jac[]> jac;
    vector<float3> deriv1;
    vector<float3> deriv2;
    vector<float3> deriv3;

    AngleCoord(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(2, grp, "id")[0], 1),
        n_elem(get_dset_size(2, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
	deriv1(n_elem),
	deriv2(n_elem),
	deriv3(n_elem)
        //jac(new_aligned<Jac>(n_elem,1))
    {
        int n_dep = 3;  // number of atoms that each term depends on 
        check_size(grp, "id",              n_elem, n_dep);
        traverse_dset<2,int> (grp, "id", [&](size_t i, size_t j, int x) { params[i].atom[j] = x;});
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("angle"));
        //float* posc = pos.output.x.get();
        VecArray posc = pos.output;
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];

            auto atom1 = load_vec<3>(posc, p.atom[0]);
            auto atom2 = load_vec<3>(posc, p.atom[1]);
            auto atom3 = load_vec<3>(posc, p.atom[2]);

            auto x1 = atom1 - atom3; auto inv_d1 = inv_mag(x1); auto x1h = x1*inv_d1;
            auto x2 = atom2 - atom3; auto inv_d2 = inv_mag(x2); auto x2h = x2*inv_d2;
            auto dp = dot(x1h, x2h);
	    output(0, nt) = dp;
            deriv1[nt] = (x2h - x1h*dp) * inv_d1;
            deriv2[nt] = (x1h - x2h*dp) * inv_d2;
            deriv3[nt] = -deriv1[nt]-deriv2[nt];
        }
    }

    virtual void propagate_deriv() {
        Timer timer(string("angle_deriv"));
        VecArray pos_sens = pos.sens;
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            update_vec(pos_sens, p.atom[0], sens(0,nt)*deriv1[nt]);
            update_vec(pos_sens, p.atom[1], sens(0,nt)*deriv2[nt]);
            update_vec(pos_sens, p.atom[2], sens(0,nt)*deriv3[nt]);
        }
    }
};
static RegisterNodeType<AngleCoord,1> angle_coord_node("Angle");

struct AngleCoord2 : public CoordNode
{
    //struct alignas(16) Jac {float d[3][4];}; // padding for vector load/store
    struct Params { index_t atom[4];};

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    //unique_ptr<Jac[]> jac;
    vector<float3> deriv1;
    vector<float3> deriv2;
    vector<float3> deriv3;
    vector<float3> deriv4;

    AngleCoord2(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(2, grp, "id")[0], 1),
        n_elem(get_dset_size(2, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
	deriv1(n_elem),
	deriv2(n_elem),
	deriv3(n_elem),
	deriv4(n_elem)
        //jac(new_aligned<Jac>(n_elem,1))
    {
        check_size(grp, "id", n_elem, 4);
        traverse_dset<2,int> (grp, "id", [&](size_t i, size_t j, int x) { params[i].atom[j] = x;});
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("vector_angle"));
        //float* posc = pos.output.x.get();
        VecArray posc = pos.output;
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];

            auto atom1 = load_vec<3>(posc, p.atom[0]);
            auto atom2 = load_vec<3>(posc, p.atom[1]);
            auto atom3 = load_vec<3>(posc, p.atom[2]);
            auto atom4 = load_vec<3>(posc, p.atom[3]);

            auto x1 = atom1 - atom2; auto inv_d1 = inv_mag(x1); auto x1h = x1*inv_d1;
            auto x2 = atom3 - atom4; auto inv_d2 = inv_mag(x2); auto x2h = x2*inv_d2;
            auto dp = dot(x1h, x2h);
	    output(0, nt) = dp;
            deriv1[nt] = (x2h - x1h*dp) * inv_d1;
            deriv3[nt] = (x1h - x2h*dp) * inv_d2;
            deriv2[nt] = -deriv1[nt];
            deriv4[nt] = -deriv3[nt];
        }
    }

    virtual void propagate_deriv() {
        Timer timer(string("vector_angle_deriv"));
        VecArray pos_sens = pos.sens;
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            update_vec(pos_sens, p.atom[0], sens(0,nt)*deriv1[nt]);
            update_vec(pos_sens, p.atom[1], sens(0,nt)*deriv2[nt]);
            update_vec(pos_sens, p.atom[2], sens(0,nt)*deriv3[nt]);
            update_vec(pos_sens, p.atom[3], sens(0,nt)*deriv4[nt]);
        }
    }
};
static RegisterNodeType<AngleCoord2,1> vector_angle_coord_node("VectorAngle");

struct DihedralCoord : public CoordNode
{
    struct alignas(16) Jac {float j[4][4];}; // padding for vector load/store
    struct Params {bool dummy_angle; index_t atom[4];};

    CoordNode& pos;
    vector<Params> params;
    unique_ptr<Jac[]> jac;

    DihedralCoord(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(2, grp, "id")[0], 1),
        pos(pos_),
        params(n_elem),
        jac(new_aligned<Jac>(n_elem,1))
    {
        check_size(grp, "id", n_elem, 4);
        traverse_dset<2,int>(grp, "id", [&](size_t nr, size_t na, int x) { params[nr].atom[na] = x;});
        for(auto& p: params) {
            // handle dummy angles uniformly (N-terminal phi and C-terminal psi)
            p.dummy_angle = p.atom[0] == index_t(-1) or p.atom[3] == index_t(-1);
        }
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("dihedral"));

        VecArray dihe_pos = output;
        float*   posv     = pos.output.x.get();

        for(int nt=0; nt<n_elem; ++nt) {
            const auto& p = params[nt];
            Float4 x[4];
            for(int na: range(4)) x[na] = Float4(posv + 4*p.atom[na]);

            Float4 d[4];
	    if (p.dummy_angle)
	        dihe_pos(0,nt) = -1.3963f;   // -80 degrees if dummy angle
	    else {
                dihe_pos(0,nt) = dihedral_germ(x[0],x[1],x[2],x[3], d[0],d[1],d[2],d[3]).x();
                for(int na: range(4)) d[na].store(jac[nt].j[na]);
            }
        }
    }

    virtual void propagate_deriv() {
        Timer timer(string("dihedral_deriv"));
        float* pos_sens = pos.sens.x.get();

        for(int nt=0; nt<n_elem; ++nt) {
            const auto& p = params[nt];

	    if (p.dummy_angle)
                continue;

            Float4 ps[4];
            for(int na: range(4))
                ps[na] = fmadd( Float4(sens(0,nt)), Float4(jac[nt].j[na]), Float4(pos_sens + 4*p.atom[na]) );

            for(int na: range(4))
                ps[na].store(pos_sens + 4*p.atom[na]);  // value was added above
        }
    }
};
static RegisterNodeType<DihedralCoord,1> dihedral_coord_node("Dihedral");

struct PosCenter : public CoordNode {
    CoordNode &pos;
    int n_dim;
    int n_group;
    int n_atom;

    vector<int> index_pos;
    vector<int> border;
    vector<float> scale;
    vector<int> index_dim;

    PosCenter(hid_t grp, CoordNode& pos_):
        CoordNode(read_attribute<int>(grp, ".", "n_group"), read_attribute<int>(grp, ".", "n_dim")),
        pos(pos_),
        n_dim(read_attribute<int>(grp, ".", "n_dim")),
        n_group(read_attribute<int>(grp, ".", "n_group")),
        n_atom(get_dset_size(1, grp, "index_pos")[0]),
        index_pos(n_atom),
	border(n_group+1),
	scale(n_group),
        index_dim(n_dim)
    {
        traverse_dset<1, int>(grp, "index_pos", [&](size_t np, int x) {index_pos[np] = x;});
        traverse_dset<1, int>(grp, "border",    [&](size_t np, int x) {border[np]    = x;});
        traverse_dset<1, int>(grp, "index_dim" ,[&](size_t ne, int x) {index_dim[ne] = x;});

	for (int i=0; i<n_group; ++i) {
            scale[i] = 1.f/(border[i+1]-border[i]);	
	}

        if(logging(LOG_DETAILED)) {
            default_logger->add_logger<float>("centers", {n_elem, n_dim}, [&](float* buffer) {
                for(int ne=0; ne<n_elem; ++ne) {
                    for(int nd=0; nd<n_dim; ++nd) 
                        buffer[ne*n_dim + nd] = output(nd, ne);
                }
            });
	}
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("pos_center");

        for(int nd=0; nd<n_dim; ++nd) {
            for(int g=0; g<n_group; ++g) {
                output(nd, g) = 0.0;
                for(int na=border[g]; na<border[g+1]; ++na) {
                    output(nd, g) += pos.output(index_dim[nd], index_pos[na]) * scale[g];
                }
            }
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_pos_center");

        for(int nd=0; nd<n_dim; ++nd) {
            for(int g=0; g<n_group; ++g) {
                auto scale_sens = sens(nd, g)*scale[g];
                for(int na=border[g]; na<border[g+1]; ++na) {
                    pos.sens(index_dim[nd], index_pos[na]) += scale_sens;
                }
            }
        }
    }
};

static RegisterNodeType<PosCenter,1> pos_center_node("GroupCenter");

