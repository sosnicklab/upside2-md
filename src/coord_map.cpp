#include <math.h>
#include <string>
#include <algorithm>
#include <iostream>
#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include "vector_math.h"

using namespace std;
using namespace h5;

struct PasteRama : public CoordNode {

    struct Param {
        index_t index_pos1;
        index_t index_pos2;
    };

    CoordNode &pos1;
    CoordNode &pos2;
    vector<Param> params;

    PasteRama(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
        CoordNode(get_dset_size(1,grp,"index_pos1")[0], 2),
        pos1(pos1_),
        pos2(pos2_),
        params(n_elem)
    {
        check_size(grp, "index_pos1", n_elem);
        check_size(grp, "index_pos2", n_elem);
        traverse_dset<1,int>(grp,"index_pos1" ,[&](size_t ne, int x){
                if(x<0 || x>=pos1.n_elem) throw string("rama_coord index_pos1 outside its input");
                params[ne].index_pos1 = x;});
        traverse_dset<1,int>(grp,"index_pos2" ,[&](size_t ne, int x){
                if(x<0 || x>=pos2.n_elem) throw string("rama_coord index_pos2 outside its input");
                params[ne].index_pos2 = x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("paste_rama");
        for(int ne=0; ne<n_elem; ++ne) {
            auto p = params[ne];
            output(0, ne) = pos1.output(0, p.index_pos1);
            output(1, ne) = pos2.output(0, p.index_pos2);
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_paste_rama");
        for(int ne=0; ne<n_elem; ++ne) {
            auto p = params[ne];
            pos1.sens(0, p.index_pos1) += sens(0,ne);
            pos2.sens(0, p.index_pos2) += sens(1,ne);
        }
    }
};
static RegisterNodeType<PasteRama,2> paste_rama_node("rama_coord");

struct CatPos : public CoordNode {
    CoordNode &pos1;
    CoordNode &pos2;
    int pos1_size;
    int pos2_size;
    int n_dim;

    vector<int> index;

    CatPos(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
        CoordNode(get_dset_size(1,grp,"index_pos1")[0]+get_dset_size(1,grp,"index_pos2")[0], read_attribute<int>(grp, ".", "n_dim")),
        pos1(pos1_),
        pos2(pos2_),
        pos1_size(get_dset_size(1,grp,"index_pos1")[0]),
        pos2_size(get_dset_size(1,grp,"index_pos2")[0]),
        n_dim(read_attribute<int>(grp, ".", "n_dim")),
        index(n_elem)
    {
        check_elem_width_lower_bound(pos1, n_dim);
        check_elem_width_lower_bound(pos2, n_dim);
        traverse_dset<1,int>(grp,"index_pos1" ,[&](size_t ne, int x){
                if(x<0 || x>=pos1.n_elem) throw string("cat_pos index_pos1 outside its input");
                index[ne] = x;});
        traverse_dset<1,int>(grp,"index_pos2" ,[&](size_t ne, int x){
                if(x<0 || x>=pos2.n_elem) throw string("cat_pos index_pos2 outside its input");
                index[ne+pos1_size] = x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("cat_pos");

        for(int ne=0; ne<pos1_size; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                output(nd, ne) = pos1.output(nd, index[ne]);
        }
        for(int ne=pos1_size; ne<n_elem; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                output(nd, ne) = pos2.output(nd, index[ne]);
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_cat_pos");

        for(int ne=0; ne<pos1_size; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                pos1.sens(nd,index[ne]) += sens(nd,ne);
        }
        for(int ne=pos1_size; ne<n_elem; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                pos2.sens(nd,index[ne]) += sens(nd,ne);
        }
    }
};

static RegisterNodeType<CatPos,2> cat_pos_node("cat_pos");

struct SigmoidCoord : public CoordNode {

    struct Params {
        index_t id;
        float c;
        float s;
	float sign;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    vector<float> deriv;
    int dim1;

    SigmoidCoord(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0],1),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
	deriv(n_elem),
        dim1(read_attribute<int>(grp, ".", "dim1"))
    {
        check_size(grp, "c",    n_elem);
        check_size(grp, "s",    n_elem);
        check_size(grp, "sign", n_elem);

        auto& p = params;
        traverse_dset<1,int>  (grp, "id",   [&](size_t i, int   x) {
                if(x<0 || x>=pos.n_elem) throw string("SigmoidCoord id outside its input");
                p[i].id   = x;});
        traverse_dset<1,float>(grp, "c",    [&](size_t i, float x) {p[i].c    = x;});
        traverse_dset<1,float>(grp, "s",    [&](size_t i, float x) {p[i].s    = x;});
        traverse_dset<1,float>(grp, "sign", [&](size_t i, float x) {p[i].sign = x;});
        if(dim1<0 || dim1>=pos.elem_width) throw string("SigmoidCoord dim1 outside its input");

        for(int ne=0; ne<n_elem; ++ne) {
            deriv[ne] = 0.f;
	}
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("sigmoid_coord");
        for(int ne=0; ne<n_elem; ++ne) {
            auto& p  = params[ne];
	    auto x   = pos.output(dim1, p.id);
	    auto sig = compact_sigmoid(x-p.c, p.s);
	    if (p.sign > 0) {
                output(0, ne) = sig.x();
		deriv[ne] = sig.y();
	    }
	    else {
                output(0, ne) = 1.f - sig.x();
		deriv[ne] = -sig.y();
            }
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_sigmoid_coord");
        VecArray pos_sens = pos.sens;
        for(int ne=0; ne<n_elem; ++ne) 
            pos_sens(dim1, params[ne].id) += sens(0, ne)*deriv[ne];
    }
};

static RegisterNodeType<SigmoidCoord,1> sigmoid_coord_node("SigmoidCoord");

