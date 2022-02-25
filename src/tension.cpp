#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"

using namespace h5;
using namespace std;

struct TensionPotential : public PotentialNode
{
    int n_elem;
    CoordNode& pos;
    struct Param {
        index_t atom;
        Vec<3>  tension_coeff;
    };
    vector<Param> params;

    TensionPotential(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(1, grp, "atom")[0]), pos(pos_), params(n_elem)
    {
        check_size(grp, "atom", n_elem);
        check_size(grp, "tension_coeff", n_elem, 3);

        auto &p = params;
        traverse_dset<1,int>  (grp,"atom",         [&](size_t i,         int   x){p[i].atom = x;});
        traverse_dset<2,float>(grp,"tension_coeff",[&](size_t i,size_t d,float x){p[i].tension_coeff[d] = x;});
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("tension"));

        VecArray pos_c = pos.output;
        VecArray pos_sens = pos.sens;

        float pot = 0.f;
        for(auto &p: params) {
            auto x = load_vec<3>(pos_c, p.atom);
            pot -= dot(x, p.tension_coeff);
            update_vec(pos_sens, p.atom, -p.tension_coeff);
        }
        potential = pot;
    }
};
static RegisterNodeType<TensionPotential,1> tension_node("tension");
