#include "deriv_engine.h"
#include <string>
#include "timing.h"
#include "affine.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include "interaction_graph.h"
#include "spline.h"
#include "state_logger.h"

using namespace std;
using namespace h5;

namespace {
struct ContactEnergy : public PotentialNode
{
    struct Param {
        index_t   loc[2];
        float     energy;
        float     dist;
        float     scale;  // 1.f/width
        float     cutoff;
    };

    int n_contact;
    CoordNode& bead_pos;
    vector<Param> params;
    float cutoff;

    ContactEnergy(hid_t grp, CoordNode& bead_pos_):
        PotentialNode(),
        n_contact(get_dset_size(2, grp, "id")[0]),
        bead_pos(bead_pos_), 
        params(n_contact)
    {
        check_size(grp, "id",       n_contact, 2);
        check_size(grp, "energy",   n_contact);
        check_size(grp, "distance", n_contact);
        check_size(grp, "width",    n_contact);

        traverse_dset<2,int  >(grp, "id",       [&](size_t nc, size_t i, int x){params[nc].loc[i] = x;});
        traverse_dset<1,float>(grp, "distance", [&](size_t nc, float x){params[nc].dist = x;});
        traverse_dset<1,float>(grp, "energy",   [&](size_t nc, float x){params[nc].energy = x;});
        traverse_dset<1,float>(grp, "width",    [&](size_t nc, float x){params[nc].scale = 1.f/x;});
        for(auto &p: params) p.cutoff = p.dist + 1.f/p.scale;

        if(logging(LOG_DETAILED)) {
            default_logger->add_logger<float>("contact_energy", {bead_pos.n_elem}, 
                    [&](float* buffer) {
                       fill_n(buffer, bead_pos.n_elem, 0.f);
                       VecArray pos  = bead_pos.output;

                       for(const auto &p: params) {
                           auto dist = mag(load_vec<3>(pos, p.loc[0]) - load_vec<3>(pos, p.loc[1]));
                           float en = p.energy * compact_sigmoid(dist-p.dist, p.scale)[0];
                           buffer[p.loc[0]] += 0.5f*en;
                           buffer[p.loc[1]] += 0.5f*en;
                       }});
        }
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("contact_energy"));
        VecArray pos  = bead_pos.output;
        VecArray sens = bead_pos.sens;
        potential = 0.f;

        for(int nc=0; nc<n_contact; ++nc) {
            const auto& p = params[nc];
            auto disp = load_vec<3>(pos, p.loc[0]) - load_vec<3>(pos, p.loc[1]);
            auto dist = mag(disp);
            if(dist>=p.cutoff) continue;

            Vec<2> contact = compact_sigmoid(dist-p.dist, p.scale);
            potential += p.energy*contact.x();
            auto deriv = (p.energy*contact.y()*rcp(dist)) * disp;
            update_vec(sens, p.loc[0],  deriv);
            update_vec(sens, p.loc[1], -deriv);
        }
    }
};

}

static RegisterNodeType<ContactEnergy,1> contact_node("contact");
