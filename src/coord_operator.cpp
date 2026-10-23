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

// every dimension inside its input, every input index inside its input and every output index inside the output
static void check_operator_indices(const CoordNode& pos1, const CoordNode& pos2, int n_size,
        const vector<int>& dim1, const vector<int>& dim2, const vector<int>& id_pos1, const vector<int>& id_pos2,
        const vector<int>& id1_out, const vector<int>& id2_out) {
    auto inside = [](const vector<int>& v, int n) {for(int i: v) if(i<0 || i>=n) return false; return true;};
    if(!inside(dim1, pos1.elem_width) || !inside(dim2, pos2.elem_width)) throw string("dim outside its input");
    if(!inside(id_pos1, pos1.n_elem) || !inside(id_pos2, pos2.n_elem)) throw string("id_pos outside its input");
    if(!inside(id1_out, n_size) || !inside(id2_out, n_size)) throw string("id_out outside the output");
}

struct Add : public CoordNode {
    CoordNode &pos1;
    CoordNode &pos2;
    int n_size;
    int n_size1;
    int n_size2;
    vector<int> id_pos1;
    vector<int> id_pos2;
    vector<int> id1_out;
    vector<int> id2_out;
    int n_dim;
    vector<int> dim1;
    vector<int> dim2;

    Add(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
        CoordNode(read_attribute<int>(grp, ".", "n_size"), read_attribute<int>(grp, ".", "n_dim")),
        pos1(pos1_),
        pos2(pos2_),
	n_size(read_attribute<int>(grp, ".", "n_size")),
	n_size1(get_dset_size(1,grp,"id_pos1")[0]),
	n_size2(get_dset_size(1,grp,"id_pos2")[0]),
	id_pos1(n_size1),
	id_pos2(n_size2),
	id1_out(n_size1),
	id2_out(n_size2),
	n_dim(read_attribute<int>(grp, ".", "n_dim")),
	dim1(n_dim),
	dim2(n_dim)
    {
        check_size(grp, "dim1", n_dim);
        check_size(grp, "dim2", n_dim);
        check_size(grp, "id1_out", n_size1);
        check_size(grp, "id2_out", n_size2);

        traverse_dset<1,int>(grp,"id_pos1" ,[&](size_t ne, int x){id_pos1[ne] = x;});
        traverse_dset<1,int>(grp,"id_pos2" ,[&](size_t ne, int x){id_pos2[ne] = x;});
        traverse_dset<1,int>(grp,"id1_out" ,[&](size_t ne, int x){id1_out[ne] = x;});
        traverse_dset<1,int>(grp,"id2_out" ,[&](size_t ne, int x){id2_out[ne] = x;});
        traverse_dset<1,int>(grp,"dim1"    ,[&](size_t ne, int x){dim1[ne]    = x;});
        traverse_dset<1,int>(grp,"dim2"    ,[&](size_t ne, int x){dim2[ne]    = x;});

	check_operator_indices(pos1, pos2, n_size, dim1, dim2, id_pos1, id_pos2, id1_out, id2_out);
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("Add");

	for(int i: range(n_size)) {
            for(int d: range(n_dim))
                output(d, i) = 0.f;
        }

	for(int i: range(n_size1)) {
            for(int d: range(n_dim))
                output(d, id1_out[i]) += pos1.output(dim1[d], id_pos1[i]);
        }

	for(int i: range(n_size2)) {
            for(int d: range(n_dim))
                output(d, id2_out[i]) += pos2.output(dim2[d], id_pos2[i]);
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_Add");
	
	for(int i: range(n_size1)) {
            for(int d: range(n_dim))
                pos1.sens(dim1[d], id_pos1[i]) += sens(d, id1_out[i]);
        }

	for(int i: range(n_size2)) {
            for(int d: range(n_dim))
                pos2.sens(dim2[d], id_pos2[i]) += sens(d, id2_out[i]);
        }
    }
};

static RegisterNodeType<Add, 2> add_node("Add");

struct Multiply: public CoordNode {
    CoordNode &pos1;
    CoordNode &pos2;
    int n_size;
    int n_size1;
    int n_size2;
    vector<int> id_pos1;
    vector<int> id_pos2;
    vector<int> id1_out;
    vector<int> id2_out;
    int n_dim;
    vector<int> dim1;
    vector<int> dim2;
    vector<int> count;
    vector<float> s1;
    vector<float> s2;

    Multiply(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
        CoordNode(read_attribute<int>(grp, ".", "n_size"), read_attribute<int>(grp, ".", "n_dim")),
        pos1(pos1_),
        pos2(pos2_),
	n_size(read_attribute<int>(grp, ".", "n_size")),
	n_size1(get_dset_size(1,grp,"id_pos1")[0]),
	n_size2(get_dset_size(1,grp,"id_pos2")[0]),
	id_pos1(n_size1),
	id_pos2(n_size2),
	id1_out(n_size1),
	id2_out(n_size2),
	n_dim(read_attribute<int>(grp, ".", "n_dim")),
	dim1(n_dim),
	dim2(n_dim),
	count(n_size),
	s1(n_size*n_dim),
	s2(n_size*n_dim)
    {
        check_size(grp, "dim1", n_dim);
        check_size(grp, "dim2", n_dim);
        check_size(grp, "id1_out", n_size1);
        check_size(grp, "id2_out", n_size2);

        traverse_dset<1,int>(grp,"id_pos1" ,[&](size_t ne, int x){id_pos1[ne] = x;});
        traverse_dset<1,int>(grp,"id_pos2" ,[&](size_t ne, int x){id_pos2[ne] = x;});
        traverse_dset<1,int>(grp,"id1_out" ,[&](size_t ne, int x){id1_out[ne] = x;});
        traverse_dset<1,int>(grp,"id2_out" ,[&](size_t ne, int x){id2_out[ne] = x;});
        traverse_dset<1,int>(grp,"dim1"    ,[&](size_t ne, int x){dim1[ne]    = x;});
        traverse_dset<1,int>(grp,"dim2"    ,[&](size_t ne, int x){dim2[ne]    = x;});

	check_operator_indices(pos1, pos2, n_size, dim1, dim2, id_pos1, id_pos2, id1_out, id2_out);

	// each output element is the product of at most one factor from each input
	vector<int> count1(n_size, 0), count2(n_size, 0);
	for(int i: id1_out) count1[i] += 1;
	for(int i: id2_out) count2[i] += 1;
	for(int i: range(n_size)) {
            if(count1[i]>1 or count2[i]>1) throw string("Multiply takes at most one factor from each input for an output element");
            count[i] = count1[i] + count2[i];
        }

	// s1 and s2 scale the sensitivity of each input's factor: they hold the other input's factor, or 1 without one
	for(int i: range(n_size*n_dim)) {
            s1[i] = 1.f;
            s2[i] = 1.f;
        }
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("Multiply");

	for(int i: range(n_size)) {
            for(int d: range(n_dim))
                output(d, i) = 1.f;
        }

	for(int i: range(n_size1)) {
            for(int d: range(n_dim)) {
                output(d, id1_out[i]) *= pos1.output(dim1[d], id_pos1[i]);
		s2[id1_out[i]*n_dim+d] = pos1.output(dim1[d], id_pos1[i]);
            }
        }

	for(int i: range(n_size2)) {
            for(int d: range(n_dim)) {
                output(d, id2_out[i]) *= pos2.output(dim2[d], id_pos2[i]);
		s1[id2_out[i]*n_dim+d] = pos2.output(dim2[d], id_pos2[i]);
            }
        }

	for(int i: range(n_size)) {
            if (count[i] > 0)
                continue;
            for(int d: range(n_dim)) 
                output(d, i) = 0.f;
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_Multiply");
	
	for(int i: range(n_size1)) {
            for(int d: range(n_dim))
                pos1.sens(dim1[d], id_pos1[i]) += sens(d, id1_out[i])* s1[id1_out[i]*n_dim+d];
        }

	for(int i: range(n_size2)) {
            for(int d: range(n_dim))
                pos2.sens(dim2[d], id_pos2[i]) += sens(d, id2_out[i])* s2[id2_out[i]*n_dim+d];
        }
    }
};

static RegisterNodeType<Multiply, 2> multiply_node("Multiply");
