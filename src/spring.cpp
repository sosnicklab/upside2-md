#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include <iostream>

using namespace h5;
using namespace std;

struct Spring : public PotentialNode
{
    struct Params {
        index_t id;
        float equil_dist;
        float spring_const;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    int n_dim;
    int dim1;
    bool pbc;
    float box_len;

    Spring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
	n_dim(pos.output.row_width),
        dim1( read_attribute<int>(grp, ".", "dim1") ),
        pbc(  read_attribute<int>(grp, ".", "pbc") ),
        box_len(read_attribute<float>(grp, ".", "box_len") )
    {
        check_size(grp, "equil_dist",   n_elem);
        check_size(grp, "spring_const", n_elem);

        auto& p = params;
        traverse_dset<1,int>  (grp, "id",           [&](size_t i, int   x) {p[i].id = x;});
        traverse_dset<1,float>(grp, "equil_dist",   [&](size_t i, float x) {p[i].equil_dist = x;});
        traverse_dset<1,float>(grp, "spring_const", [&](size_t i, float x) {p[i].spring_const = x;});
	assert(dim1<n_dim);

	if (pbc) assert(box_len>0);
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("pos_spring_1d"));
        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            float dist   = posc(dim1, p.id);
	    float excess = dist-p.equil_dist;
	    float sqr_excess = sqr(excess);

	    if (excess == 0.f) continue;

	    if (pbc) {
                float excess1 = dist-p.equil_dist-box_len;
                float excess2 = dist-p.equil_dist+box_len;
	        float sqr_excess1 = sqr(excess1);
	        float sqr_excess2 = sqr(excess2);
		if (sqr_excess1 < sqr_excess) {
                    sqr_excess = sqr_excess1;
                    excess     = excess1;
		}
		if (sqr_excess2 < sqr_excess) {
		    sqr_excess = sqr_excess2;
                    excess     = excess2;
		}
	    }

            if(pot) *pot         += 0.5f * p.spring_const * sqr_excess;
            pos_sens(dim1, p.id) +=        p.spring_const * excess;
        }
    }
};
static RegisterNodeType<Spring,1> spring_node("Spring");

struct WallSpring : public PotentialNode
{
    struct Params {
        index_t id;
        float equil_dist;
        float spring_const;
	int wall_type;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    int n_dim;
    int dim1;

    WallSpring(hid_t grp, CoordNode& pos_):
        PotentialNode(),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
	n_dim(pos.output.row_width),
        dim1( read_attribute<int>(grp, ".", "dim1") )
    {
        check_size(grp, "equil_dist",   n_elem);
        check_size(grp, "spring_const", n_elem);
        check_size(grp, "wall_type",    n_elem);

        auto& p = params;
        traverse_dset<1,int>  (grp, "id",           [&](size_t i, int   x) {p[i].id = x;});
        traverse_dset<1,float>(grp, "equil_dist",   [&](size_t i, float x) {p[i].equil_dist = x;});
        traverse_dset<1,float>(grp, "spring_const", [&](size_t i, float x) {p[i].spring_const = x;});
        traverse_dset<1,int>  (grp, "wall_type",    [&](size_t i, int   x) {p[i].wall_type = x;});
	assert(dim1<n_dim);
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("pos_spring_1d"));
        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;
        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        if(pot) *pot = 0.f;
        for(int nt=0; nt<n_elem; ++nt) {
            auto& p = params[nt];
            float dist   = posc(dim1, p.id);
	    float excess = dist-p.equil_dist;
            if (p.wall_type==0 and excess>=0.f) 
                excess = 0.f;
	    else if (p.wall_type==1 and excess <= 0.f) 
                excess = 0.f;

	    if (excess != 0.f) {
                if(pot) *pot         += 0.5f*p.spring_const*sqr(excess);
                pos_sens(dim1, p.id) += p.spring_const*excess;
	    }
        }
    }
};
static RegisterNodeType<WallSpring,1> wall1_node("WallSpring");

//struct MovingSpring : public PotentialNode
//{
//    struct Params {
//        index_t id;
//        float starting_equil_dist;
//        float spring_const;
//	float pulling_vel;
//    };
//
//    int n_elem;
//    vector<Param> params;
//    CoordNode& pos;
//
//    float time_initial;
//    float time_step;      // WARNING: this should be the same as the global time step !
//    float time_estimate;  // WARNING: changes frequently
//    int round_num;        // WARNING: changes frequently
//    
//    MovingSpring(hid_t grp, CoordNode& pos_):
//        PotentialNode(),
//        n_elem(get_dset_size(1, grp, "atom")[0]),
//        params(n_elem),
//        pos(pos_),
//        time_initial(read_attribute<float>(grp, "pulling_vel", "time_initial")),
//        time_step(read_attribute<float>(grp, "pulling_vel", "time_step")),
//        time_estimate(time_initial),
//        round_num(0),
//    {
//        check_size(grp, "atom",                n_elem);
//        check_size(grp, "spring_const",        n_elem);
//        check_size(grp, "starting_equil_dist", n_elem);
//        check_size(grp, "pulling_vel",         n_elem);
//        
//        auto &p = params;
//        traverse_dset<1,int>  (grp,"atom",                [&](size_t i,int   x){p[i].atom                = x;});
//        traverse_dset<1,float>(grp,"spring_const",        [&](size_t i,float x){p[i].spring_const        = x;});
//        traverse_dset<1,float>(grp,"starting_equil_dist", [&](size_t i,float x){p[i].starting_equil_dist = x;});
//        traverse_dset<1,float>(grp,"pulling_vel",         [&](size_t i,float x){p[i].pulling_vel         = x;});
//        
//        if(logging(LOG_BASIC)) {
//            default_logger->add_logger<float>("tip_pos", {n_elem, 3}, [&](float* buffer) {
//                
//                time_estimate = time_initial + float(time_step)*round_num;
//                
//                for(int nt=0; nt<n_elem; ++nt) {
//                    auto p = params[nt];
//                    auto tip_pos = p.starting_tip_pos + p.pulling_vel * time_estimate;
//                    for(int d=0; d<3; ++d) buffer[nt*3+d] = tip_pos[d];
//                }
//            });
//           
//            default_logger->add_logger<float>("time_estimate", {1}, [&](float* buffer) {
//                buffer[0] = time_estimate;
//            });
//        }
//    }
//
//    virtual void compute_value(ComputeMode mode) {
//        Timer timer(string("AFM"));
//        
//        if (mode == DerivMode) round_num += 1;
//        time_estimate = time_initial + float(time_step)*round_num;
//        
//        VecArray pos_c    = pos.output;
//        VecArray pos_sens = pos.sens;
//        
//        float pot = 0.f;
//        for(auto &p: params) {
//            auto x = load_vec<3>(pos_c, p.atom);
//            auto tip_pos = p.starting_tip_pos + p.pulling_vel * time_estimate;
//            auto x_diff = x - tip_pos;
//            pot += 0.5*p.spring_const*mag2(x_diff);
//            update_vec(pos_sens, p.atom, p.spring_const*x_diff);
//        }
//        potential = pot;
//    }
//};
//static RegisterNodeType<AFMPotential,1> AFM_node("AFM");
//
