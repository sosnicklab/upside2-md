#include <math.h>
#include <string>
#include <algorithm>
#include <iostream>
#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include "vector_math.h"
#include <Eigen/Dense>

using namespace std;
using namespace h5;
using namespace Eigen;

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

	bool check_s = true;
	for(int i: dim1) if(i>=pos1.output.row_width) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: dim2) if(i>=pos2.output.row_width) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id_pos1) if(i>=pos1.output.n_elem) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id_pos2) if(i>=pos2.output.n_elem) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id1_out) if(i>=output.n_elem) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id2_out) if(i>=output.n_elem) {check_s = false; break;};

	if (not check_s) cout << "The ID you gave is out of range!" << endl;
	assert(check_s);
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

struct Mean : public CoordNode {
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
    vector<float> s;

    Mean(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
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
	s(n_size)
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

	bool check_s = true;
	for(int i: dim1) if(i>=pos1.output.row_width) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: dim2) if(i>=pos2.output.row_width) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id_pos1) if(i>=pos1.output.n_elem) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id_pos2) if(i>=pos2.output.n_elem) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id1_out) if(i>=output.n_elem) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id2_out) if(i>=output.n_elem) {check_s = false; break;};

	if (not check_s) cout << "The ID you gave is out of range!" << endl;
	assert(check_s);

	for(int i: range(n_size)) count[i]  = 0;
	for(int i: id1_out)       count[i] += 1;
	for(int i: id2_out)       count[i] += 1;

	for(int i: range(n_size)) {
            if (count[i] == 0)
                s[i] = 1.f;
            else
                s[i] = 1.f/count[i];
        }
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("Mean");

	for(int i: range(n_size)) {
            for(int d: range(n_dim))
                output(d, i) = 0.f;
        }

	for(int i: range(n_size1)) {
            for(int d: range(n_dim))
                output(d, id1_out[i]) += pos1.output(dim1[d], id_pos1[i]) * s[id1_out[i]];
        }

	for(int i: range(n_size2)) {
            for(int d: range(n_dim))
                output(d, id2_out[i]) += pos2.output(dim2[d], id_pos2[i]) * s[id2_out[i]];
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_Mean");
	
	for(int i: range(n_size1)) {
            for(int d: range(n_dim))
                pos1.sens(dim1[d], id_pos1[i]) += sens(d, id1_out[i])* s[id1_out[i]];
        }

	for(int i: range(n_size2)) {
            for(int d: range(n_dim))
                pos2.sens(dim2[d], id_pos2[i]) += sens(d, id2_out[i])* s[id2_out[i]];
        }
    }
};

static RegisterNodeType<Mean, 2> mean_node("Mean");

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
	s1(n_size),
	s2(n_size)
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

	bool check_s = true;
	for(int i: dim1) if(i>=pos1.output.row_width) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: dim2) if(i>=pos2.output.row_width) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id_pos1) if(i>=pos1.output.n_elem) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id_pos2) if(i>=pos2.output.n_elem) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id1_out) if(i>=output.n_elem) {check_s = false; break;};
	assert(check_s);
	check_s = true;
	for(int i: id2_out) if(i>=output.n_elem) {check_s = false; break;};

	if (not check_s) cout << "The ID you gave is out of range!" << endl;
	assert(check_s);

	for(int i: range(n_size)) count[i]  = 0;
	for(int i: id1_out)       count[i] += 1;
	for(int i: id2_out)       count[i] += 1;

	for(int i: range(n_size)) {
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
		s2[id1_out[i]] = pos1.output(dim1[d], id_pos1[i]);
            }
        }

	for(int i: range(n_size2)) {
            for(int d: range(n_dim)) {
                output(d, id2_out[i]) *= pos2.output(dim2[d], id_pos2[i]);
		s1[id2_out[i]] = pos2.output(dim2[d], id_pos2[i]);
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
                pos1.sens(dim1[d], id_pos1[i]) += sens(d, id1_out[i])* s2[id1_out[i]];
        }

	for(int i: range(n_size2)) {
            for(int d: range(n_dim))
                pos2.sens(dim2[d], id_pos2[i]) += sens(d, id2_out[i])* s1[id2_out[i]];
        }
    }
};

static RegisterNodeType<Multiply, 2> multiply_node("Multiply");

struct Sum: public CoordNode {
    CoordNode &pos1;
    int n_size;
    vector<int> id_pos;
    int n_dim;
    vector<int> dim;
    vector<float> weight;

    Sum(hid_t grp, CoordNode& pos1_):
        CoordNode(1, read_attribute<int>(grp, ".", "n_dim")),
        pos1(pos1_),
	n_size(get_dset_size(1,grp,"id_pos")[0]),
	id_pos(n_size),
	n_dim(read_attribute<int>(grp, ".", "n_dim")),
	dim(n_dim),
	weight(n_size)
    {
        check_size(grp, "dim", n_dim);

        traverse_dset<1,int>  (grp, "id_pos" ,[&](size_t ne, int x)   {id_pos[ne] = x;});
        traverse_dset<1,int>  (grp, "dim"    ,[&](size_t ne, int x)   {dim[ne]    = x;});
        traverse_dset<1,float>(grp, "weight" ,[&](size_t ne, float x) {weight[ne] = x;});

	bool check_s = true;
	for(int i: dim) if(i>=pos1.output.row_width) {check_s = false; break;};
	if (not check_s) cout << "The dim you gave is out of range!" << endl;
	assert(check_s);
	check_s = true;
	for(int i: id_pos) if(i>=pos1.output.n_elem) {check_s = false; break;};
	if (not check_s) cout << "The id you gave is out of range!" << endl;
	assert(check_s);
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("Sum");

        for(int d: range(n_dim)) 
            output(d, 0) = 0.0;
	for(int i: range(n_size)) {
            for(int d: range(n_dim)) output(d, 0) += weight[i]*pos1.output(dim[d], id_pos[i]);
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_Sum");
	
	for(int i: range(n_size)) {
            for(int d: range(n_dim))
                pos1.sens(dim[d], id_pos[i]) += sens(d, 0)* weight[id_pos[i]];
        }
    }
};

static RegisterNodeType<Sum, 1> Sum_node("Sum");

