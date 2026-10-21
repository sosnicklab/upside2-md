#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include <iostream>

using namespace h5;
using namespace std;

struct SigmoidEnergy : public PotentialNode
{
    struct Params {
        index_t id;
        float a;
        float b;
        float c;
        float s;
	float min;
	float max;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    int n_dim;
    int dim1;
    bool use_cutoff;

    SigmoidEnergy(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
	n_dim(pos.output.row_width),
        dim1( read_attribute<int>(grp, ".", "dim1") ),
        use_cutoff( read_attribute<int>(grp, ".", "use_cutoff") )
    {
        check_size(grp, "a",   n_elem);
        check_size(grp, "b",   n_elem);
        check_size(grp, "c",   n_elem);
        check_size(grp, "s",   n_elem);
        check_size(grp, "min", n_elem);
        check_size(grp, "max", n_elem);

        auto& p = params;
        traverse_dset<1,int>  (grp, "id",  [&](size_t i, int   x) {p[i].id  = x;});
        traverse_dset<1,float>(grp, "a",   [&](size_t i, float x) {p[i].a   = x;});
        traverse_dset<1,float>(grp, "b",   [&](size_t i, float x) {p[i].b   = x;});
        traverse_dset<1,float>(grp, "c",   [&](size_t i, float x) {p[i].c   = x;});
        traverse_dset<1,float>(grp, "s",   [&](size_t i, float x) {p[i].s   = x;});
        traverse_dset<1,float>(grp, "min", [&](size_t i, float x) {p[i].min = x;});
        traverse_dset<1,float>(grp, "max", [&](size_t i, float x) {p[i].max = x;});
	assert(dim1<n_dim);

	if (not use_cutoff) {
            for(int nt=0; nt<n_elem; ++nt) {
                auto& p = params[nt];
	        p.min =  p.c - 1.f/p.s;
	        p.max =  p.c + 1.f/p.s;
            }
        }
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("sigmoid_energy"));
        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;

        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
	    float x = posc(dim1, p.id);

            if (use_cutoff and (x>=p.max or x<p.min)) {  // [min,max): adjacent windows share a border
                continue;
	    }
	    else if (x>p.max) {
	        if(pot) *pot += p.a;
	    }
	    else if (x<p.min) {
	        if(pot) *pot += p.a + p.b;
	    }
	    else {
	        Vec<2,float> sig = compact_sigmoid(x-p.c, p.s);
	        if(pot) *pot += p.a + p.b*sig.x();
                pos_sens(dim1, p.id) += p.b*sig.y();
            }
        }
    }
};
static RegisterNodeType<SigmoidEnergy,1> sigmoid_node("SigmoidEnergy");

