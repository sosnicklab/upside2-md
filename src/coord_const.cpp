#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include "vector_math.h"
#include <iostream>

using namespace h5;
using namespace std;

struct ConstantCoord3D : public CoordNode
{
    ConstantCoord3D(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0], 3)
    {
        check_size(grp, "value", n_elem, 3);  // one row of value for each row of id
        traverse_dset<2,float>(grp, "value", [&](size_t ne, size_t nt, float x) { output(nt, ne) = x;});

        if(logging(LOG_DETAILED)) {
            if (default_logger) {
                bool suff_exists = true;
                int log_suffix = 0;
                while(suff_exists) {
                    if(h5_exists(default_logger->logging_group.get(), ("Const3D_" + to_string(log_suffix)).c_str())) {
                        log_suffix++;
                    } else {
                        suff_exists = false;
                    }
                }
                default_logger->add_logger<float>( ("Const3D_" + to_string(log_suffix)).c_str(), {n_elem, 3}, [&](float* buffer) {
                    for(int nt=0; nt<n_elem; ++nt) {
                        auto x = load_vec<3>(output, nt);
                        for(int d=0; d<3; ++d) buffer[nt*3 + d] = x[d];
                    }
                });
            }
        }
    }

    virtual void compute_value(ComputeMode mode) {}
    virtual void propagate_deriv() {}
};
static RegisterNodeType<ConstantCoord3D,1> const3d_node("Const3D");

struct MovingConstantCoord3D : public CoordNode
{
    struct Params { 
        float3 velocities;
        float3 start_pos;
    };

    vector<Params> params;

    float time_initial;
    float time_step;      // WARNING: this should be the same as the global time step !
    float time_estimate;  // WARNING: changes frequently
    long round_num;       // time steps taken since time_initial: @step_initial before this run, then one per force evaluation

    MovingConstantCoord3D(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0], 3),
	params(n_elem),
        time_initial(read_attribute<float>(grp, ".", "time_initial")),
        time_step(read_attribute<float>(grp,    ".", "time_step")),
        round_num(read_attribute<long>(grp, ".", "step_initial"))
    {
        // A continuation resumes the step count (its frame's /output/tip_step), so that its clock is the continuous
        // run's time_initial + time_step*steps, bit for bit; restarting the count at a shifted time_initial would
        // round differently
        if(round_num < 0) throw string("MovingConst3D step_initial must not be negative");
        time_estimate = time_initial + float(time_step)*round_num;
        auto &p = params;
        check_size(grp, "velocities",  n_elem, 3);
        check_size(grp, "start_pos",   n_elem, 3);
        traverse_dset<2, float>(grp,"velocities",  [&](size_t i, size_t j, float x) {p[i].velocities[j]  = x;});
        traverse_dset<2, float>(grp,"start_pos",   [&](size_t i, size_t j, float x) {p[i].start_pos[j]   = x;});

        if(logging(LOG_BASIC)) {
            default_logger->add_logger<float>("tip_pos", {n_elem, 3}, [&](float* buffer) {
                for(int nt=0; nt<n_elem; ++nt) {
                    auto p = params[nt];
                    for(int d=0; d<3; ++d) 
		        buffer[nt*3+d] = p.start_pos[d] + p.velocities[d] * time_estimate;
                }
            });
           
            default_logger->add_logger<float>("time_estimate", {1}, [&](float* buffer) {
                buffer[0] = time_estimate;
            });

            default_logger->add_logger<long>("tip_step", {1}, [&](long* buffer) {
                buffer[0] = round_num;
            });
        }
    }

    virtual void compute_value(ComputeMode mode) {
        if (mode == DerivMode) round_num += 1;
        time_estimate = time_initial + float(time_step)*round_num;
        for(int nt=0; nt<n_elem; ++nt) {
            auto& pa = params[nt];
	    for(int i: range(3)) {
                output(i, nt) = pa.start_pos[i] + pa.velocities[i] * time_estimate;
	    }
        }
    }
    virtual void propagate_deriv() {}
};
static RegisterNodeType<MovingConstantCoord3D,1> moving_const3d_node("MovingConst3D");

