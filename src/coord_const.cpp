#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include "vector_math.h"
#include <iostream>

using namespace h5;
using namespace std;

struct ConstantCoord1D : public CoordNode
{
    struct Params { index_t id; };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    bool initialized_by_coord;
    bool initialized;
    int dim;

    ConstantCoord1D(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0], 1),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
        initialized_by_coord(read_attribute<int>(grp, ".", "initialized_by_coord")),
        initialized(false),
        dim(read_attribute<int>(grp, ".", "dim"))
    {
	if (initialized_by_coord) {
	    assert(dim<pos.output.row_width);
            traverse_dset<1,int> (grp, "id", [&](size_t i, int x) {params[i].id = x;});
	}
	else{
            n_elem = get_dset_size(1, grp, "value")[0];
            check_size(grp, "value", n_elem, 1);
            traverse_dset<1,float>(grp, "value", [&](size_t ne, float x) { output(0,ne) = x;});
	}


    }

    virtual void compute_value(ComputeMode mode) {
	if (initialized_by_coord and not initialized) {
            VecArray posc = pos.output;
            for(int nt=0; nt<n_elem; ++nt) {
                auto& p = params[nt];
	        output(0, nt) = posc(dim, p.id);
            }
	    initialized = true;
	}
    }
    virtual void propagate_deriv() {}
};
static RegisterNodeType<ConstantCoord1D,1> const1d_node("Const1D");

struct ConstantCoord2D : public CoordNode
{
    struct Params { index_t id; };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    bool initialized_by_coord;
    bool initialized;
    int dim1;
    int dim2;

    ConstantCoord2D(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0], 2),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
        initialized_by_coord(read_attribute<int>(grp, ".", "initialized_by_coord")),
        initialized(false),
        dim1(read_attribute<int>(grp, ".", "dim1")),
        dim2(read_attribute<int>(grp, ".", "dim2"))
    {
	if (initialized_by_coord) {
	    assert(dim1<pos.output.row_width);
	    assert(dim2<pos.output.row_width);
            traverse_dset<1,int> (grp, "id", [&](size_t i, int x) {params[i].id = x;});
	}
	else{
            n_elem = get_dset_size(2, grp, "value")[0];
            check_size(grp, "value", n_elem, 2);
            traverse_dset<2,float>(grp, "value", [&](size_t ne, size_t nt, float x) { output(nt, ne) = x;});
	}
    }

    virtual void compute_value(ComputeMode mode) {
    
	if (initialized_by_coord and not initialized) {
            VecArray posc = pos.output;
            for(int nt=0; nt<n_elem; ++nt) {
                auto& p = params[nt];
	        output(0, nt) = posc(dim1, p.id);
	        output(1, nt) = posc(dim2, p.id);
            }
	    initialized = true;
	}
    
    }
    virtual void propagate_deriv() {}
};
static RegisterNodeType<ConstantCoord2D,1> const2d_node("Const2D");

struct ConstantCoord3D : public CoordNode
{
    struct Params { index_t id; };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    bool initialized_by_coord;
    bool initialized;
    int dim1;
    int dim2;
    int dim3;

    ConstantCoord3D(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0], 3),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
        initialized_by_coord(read_attribute<int>(grp, ".", "initialized_by_coord")),
        initialized(false),
        dim1(read_attribute<int>(grp, ".", "dim1")),
        dim2(read_attribute<int>(grp, ".", "dim2")),
        dim3(read_attribute<int>(grp, ".", "dim3"))
    {
	if (initialized_by_coord) {
	    assert(dim1<pos.output.row_width);
	    assert(dim2<pos.output.row_width);
	    assert(dim3<pos.output.row_width);
            traverse_dset<1,int> (grp, "id", [&](size_t i, int x) {params[i].id = x;});
	}
	else{
            n_elem = get_dset_size(2, grp, "value")[0];
            check_size(grp, "value", n_elem, 3);
            traverse_dset<2,float>(grp, "value", [&](size_t ne, size_t nt, float x) { output(nt, ne) = x;});
	}

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

    virtual void compute_value(ComputeMode mode) {
	if (initialized_by_coord and not initialized) {
            VecArray posc = pos.output;
            for(int nt=0; nt<n_elem; ++nt) {
                auto& p = params[nt];
	        output(0, nt) = posc(dim1, p.id);
	        output(1, nt) = posc(dim2, p.id);
	        output(2, nt) = posc(dim3, p.id);
            }
	    initialized = true;
	}
    }
    virtual void propagate_deriv() {}
};
static RegisterNodeType<ConstantCoord3D,1> const3d_node("Const3D");

struct MovingConstantCoord1D : public CoordNode
{
    struct Params { 
	index_t id; 
        float velocities;
        float start_pos;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    bool initialized_by_coord;
    bool initialized;
    int dim;

    float time_initial;
    float time_step;      // WARNING: this should be the same as the global time step !
    float time_estimate;  // WARNING: changes frequently
    int round_num;        // WARNING: changes frequently

    MovingConstantCoord1D(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0], 1),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
        initialized_by_coord(read_attribute<int>(grp, ".", "initialized_by_coord")),
        initialized(false),
        time_initial(read_attribute<float>(grp, ".", "time_initial")),
        time_step(read_attribute<float>(grp,    ".", "time_step")),
        time_estimate(time_initial),
        round_num(0)
       
    {
        auto &p = params;
        traverse_dset<1,int> (grp, "id", [&](size_t i, int x) {p[i].id = x;});

	if (initialized_by_coord) {
            dim = read_attribute<int>(grp, ".", "dim");
	    assert(dim<pos.output.row_width);
	}
	else{
            n_elem = get_dset_size(1, grp, "start_pos")[0];
            check_size(grp, "start_pos", n_elem);
            traverse_dset<1,float>(grp, "start_pos", [&](size_t i, float x) { output(0,i) = x; p[i].start_pos = x;});
	}

        check_size(grp, "velocities",  n_elem);
        traverse_dset<1,float>(grp, "velocities", [&](size_t i, float x) {p[i].velocities = x;});

    }

    virtual void compute_value(ComputeMode mode) {

	if (initialized_by_coord and not initialized) {
            VecArray posc = pos.output;
            for(int nt=0; nt<n_elem; ++nt) {
                auto& pa = params[nt];
	        output(0, nt) = posc(dim, pa.id);
		pa.start_pos  = posc(dim, pa.id);
            }
	    initialized = true;
	}

        if (mode == DerivMode) round_num += 1;
        time_estimate = time_initial + float(time_step)*round_num;
        for(auto &p: params) {
            output(0, p.id) = p.start_pos + p.velocities * time_estimate;
        }
    }
    virtual void propagate_deriv() {}
};
static RegisterNodeType<MovingConstantCoord1D,1> moving_const1d_node("MovingConst1D");

struct WhirlingConstantCoord1D : public CoordNode
{
    struct Params { 
	index_t id; 
        float whirling_vel;
        float start_angle;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    bool initialized_by_coord;
    bool initialized;
    int dim;

    float time_initial;
    float time_step;      // WARNING: this should be the same as the global time step !
    float time_estimate;  // WARNING: changes frequently
    int round_num;        // WARNING: changes frequently

    WhirlingConstantCoord1D(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0], 1),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
        initialized_by_coord(read_attribute<int>(grp, ".", "initialized_by_coord")),
        initialized(false),
        time_initial(read_attribute<float>(grp, ".", "time_initial")),
        time_step(read_attribute<float>(grp,    ".", "time_step")),
        time_estimate(time_initial),
        round_num(0)
       
    {
        auto &p = params;
        traverse_dset<1,int> (grp, "id", [&](size_t i, int x) {p[i].id = x;});

	if (initialized_by_coord) {
            dim = read_attribute<int>(grp, ".", "dim");
	    assert(dim<pos.output.row_width);
	}
	else{
            n_elem = get_dset_size(1, grp, "start_angle")[0];
            check_size(grp, "start_angle", n_elem);
            traverse_dset<1,float>(grp, "start_angle", [&](size_t i, float x) { output(0,i) = x; p[i].start_angle = x;});
	}

        check_size(grp, "whirling_vel",  n_elem);
        traverse_dset<1,float>(grp,"whirling_vel",  [&](size_t i, float x) {p[i].whirling_vel = x;});
    }

    virtual void compute_value(ComputeMode mode) {

	if (initialized_by_coord and not initialized) {
            VecArray posc = pos.output;
            for(int nt=0; nt<n_elem; ++nt) {
                auto& pa = params[nt];
	        output(0, nt)  = posc(dim, pa.id);
		pa.start_angle = posc(dim, pa.id);
            }
	    initialized = true;
	}

        if (mode == DerivMode) round_num += 1;
        time_estimate = time_initial + float(time_step)*round_num;
        for(auto &p: params) {
	    float angle = p.start_angle + p.whirling_vel * time_estimate;
            if (angle > M_PI_F)
                angle -= 2*M_PI_F;
            if (angle < -M_PI_F)
                angle += 2*M_PI_F;
            output(0, p.id) = angle;
        }
    }
    virtual void propagate_deriv() {}
};
static RegisterNodeType<WhirlingConstantCoord1D,1> whirling_const1d_node("WhirlingConst1D");

struct MovingConstantCoord2D : public CoordNode
{
    struct Params { 
	index_t id; 
        float2 velocities;
        float2 start_pos;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    bool initialized_by_coord;
    bool initialized;
    int dim1;
    int dim2;

    float time_initial;
    float time_step;      // WARNING: this should be the same as the global time step !
    float time_estimate;  // WARNING: changes frequently
    int round_num;        // WARNING: changes frequently

    MovingConstantCoord2D(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0], 2),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
        initialized_by_coord(read_attribute<int>(grp, ".", "initialized_by_coord")),
        initialized(false),
        time_initial(read_attribute<float>(grp, ".", "time_initial")),
        time_step(read_attribute<float>(grp,    ".", "time_step")),
        time_estimate(time_initial),
        round_num(0)
    {
        auto &p = params;
        traverse_dset<1,int> (grp, "id", [&](size_t i, int x) {p[i].id = x;});

	if (initialized_by_coord) {
            dim1 = read_attribute<int>(grp, ".", "dim1");
            dim2 = read_attribute<int>(grp, ".", "dim2");
	    assert(dim1<pos.output.row_width);
	    assert(dim2<pos.output.row_width);
	}
	else{
            n_elem = get_dset_size(1, grp, "start_pos")[0];
            check_size(grp, "start_pos", n_elem, 2);
            traverse_dset<2,float>(grp, "start_pos", [&](size_t i, size_t j, float x) { output(j,i) = x; p[i].start_pos[j] = x;});
	}

        check_size(grp, "velocities",  n_elem, 2);
        traverse_dset<2,float>(grp,"velocities",  [&](size_t i, size_t j, float x) {p[i].velocities[j]  = x;});

    }

    virtual void compute_value(ComputeMode mode) {

	if (initialized_by_coord and not initialized) {
            VecArray posc = pos.output;
            for(int nt=0; nt<n_elem; ++nt) {
                auto& pa = params[nt];
	        output(0, nt)  = posc(dim1, pa.id);
	        output(1, nt)  = posc(dim2, pa.id);
		pa.start_pos[0] = posc(dim1, pa.id);
		pa.start_pos[1] = posc(dim2, pa.id);
            }
	    initialized = true;
	}

        if (mode == DerivMode) round_num += 1;
        time_estimate = time_initial + float(time_step)*round_num;
        for(auto &p: params) {
	    for(int i: range(2)) {
                output(i, p.id) = p.start_pos[i] + p.velocities[i] * time_estimate;
	    }
        }
    }
    virtual void propagate_deriv() {}
};
static RegisterNodeType<MovingConstantCoord2D,1> moving_const2d_node("MovingConst2D");

struct MovingConstantCoord3D : public CoordNode
{
    struct Params { 
	index_t id; 
        float3 velocities;
        float3 start_pos;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    bool initialized_by_coord;
    bool initialized;
    int dim1;
    int dim2;
    int dim3;

    float time_initial;
    float time_step;      // WARNING: this should be the same as the global time step !
    float time_estimate;  // WARNING: changes frequently
    long round_num;        // WARNING: changes frequently

    MovingConstantCoord3D(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0], 3),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
        initialized_by_coord(read_attribute<int>(grp, ".", "initialized_by_coord")),
        initialized(false),
	dim1(0),
	dim2(1),
	dim3(2),
        time_initial(read_attribute<float>(grp, ".", "time_initial")),
        time_step(read_attribute<float>(grp,    ".", "time_step")),
        time_estimate(time_initial),
        round_num(0)
       
    {
        auto &p = params;
        traverse_dset<1,int> (grp, "id", [&](size_t i, int x) {p[i].id = x;});

        dim1 = read_attribute<int>(grp, ".", "dim1");
        dim2 = read_attribute<int>(grp, ".", "dim2");
        dim3 = read_attribute<int>(grp, ".", "dim3");

        check_size(grp, "velocities",  n_elem, 3);
        traverse_dset<2, float>(grp,"velocities",  [&](size_t i, size_t j, float x) {p[i].velocities[j]  = x;});

	if (initialized_by_coord) {
	    assert(dim1<pos.output.row_width);
	    assert(dim2<pos.output.row_width);
	    assert(dim3<pos.output.row_width);
	}
	else {
            check_size(grp, "start_pos", n_elem, 3);
            traverse_dset<2,float>(grp, "start_pos", [&](size_t i, size_t j, float x) {p[i].start_pos[j] = x;});
	}

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
        }
    }

    virtual void compute_value(ComputeMode mode) {

	if (initialized_by_coord and not initialized) {
            VecArray posc = pos.output;
            for(int nt=0; nt<n_elem; ++nt) {
                auto& pa = params[nt];
	        output(0, nt)   = posc(dim1, pa.id);
	        output(1, nt)   = posc(dim2, pa.id);
	        output(2, nt)   = posc(dim3, pa.id);
		pa.start_pos[0] = posc(dim1, pa.id);
		pa.start_pos[1] = posc(dim2, pa.id);
		pa.start_pos[2] = posc(dim3, pa.id);
            }

	    initialized = true;
	}

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

