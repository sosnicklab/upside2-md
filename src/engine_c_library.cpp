#include "engine_c_library.h"
#include "deriv_engine.h"
#include "martini.h"
#include <algorithm>

using namespace h5;
using namespace std;

DerivEngine* construct_deriv_engine(int n_atom, const char* potential_file) try {
    H5Obj config = h5_obj(H5Fclose, H5Fopen(potential_file, H5F_ACC_RDONLY, H5P_DEFAULT));
    auto potential_group = open_group(config.get(), "/input/potential");
    
    auto engine = new DerivEngine(initialize_engine_from_hdf5(n_atom, potential_group.get() ));

    // Register the same per-engine metadata main.cpp does. Without this the library silently evaluated a
    // different model than the simulator for a hybrid system: the BB proxy state was never registered, so
    // martini_hybrid_position fell through to a plain passthrough and the derived BB and O slots behaved
    // as free particles carrying environment force directly.
    martini_masses::load_masses_for_engine(engine, config.get());
    martini_fix_rigid::register_fix_rigid_for_engine(config.get(), *engine);
    martini_stage_params::register_stage_params_for_engine(engine, config.get());
    martini_hybrid::register_hybrid_for_engine(config.get(), *engine);

    engine->build_exec_levels();
    return engine;
} catch(const string& e) {
    fprintf(stderr, "\n\nERROR: %s\n", e.c_str());
    return 0;
} catch(...) {
    return 0;
}


void free_deriv_engine(DerivEngine* engine) {
    // like delete, a no-op on the NULL that construct_deriv_engine returns for a refused config
    if(engine) martini_unregister_engine(*engine);
    delete engine;
}


// 0 indicates success, anything else is failure
int evaluate_energy(float* energy, DerivEngine* engine, const float* pos) try {
    VecArray a = engine->pos->output;
    for(int na: range(engine->pos->n_atom))
        for(int d: range(3))
            a(d,na) = pos[na*3+d];
    engine->compute(PotentialAndDerivMode);
    *energy = engine->potential;
    return 0;
} catch(const char* e) {
    fprintf(stderr, "\n\nERROR: %s\n", e);
    return 1;
} catch(const string& e) {
    fprintf(stderr, "\n\nERROR: %s\n", e.c_str());
    return 1;
} catch(...) {
    return 1;
}

// 0 indicates success, anything else is failure
int evaluate_deriv(float* deriv, DerivEngine* engine, const float* pos) try {
    // result is size (n_atom,3)
    VecArray a = engine->pos->output;
    for(int na: range(engine->pos->n_atom))
        for(int d: range(3))
            a(d,na) = pos[na*3+d];

    engine->compute(PotentialAndDerivMode);

    VecArray b = engine->pos->sens;
    for(int na: range(engine->pos->n_atom))
        for(int d: range(3))
            deriv[na*3+d] = b(d,na);
    return 0;
} catch(const char* e) {
    fprintf(stderr, "\n\nERROR: %s\n", e);
    return 1;
} catch(const string& e) {
    fprintf(stderr, "\n\nERROR: %s\n", e.c_str());
    return 1;
} catch(...) {
    return 1;
}


int set_param(int n_param, const float* param, DerivEngine* engine, const char* node_name) try {
    auto& computation = engine->get(string(node_name)).computation;
    size_t n_expected = computation->get_param().size();
    if(n_expected != size_t(n_param))
        throw string("Wrong number of parameters, expected ") + to_string(n_expected) + " but got " + to_string(n_param);
    computation->set_param(vector<float>(param, param+n_param));
    return 0;
} catch(const string& s) {
    fprintf(stderr, "ERROR: %s\n", s.c_str());
    return 1;
} catch(...) {
    fprintf(stderr, "ERROR: %s\n", "unknown error");
    return 1;
}

int get_param(int n_param, float* param, DerivEngine* engine, const char* node_name) try {
    auto param_v = engine->get(string(node_name)).computation->get_param();
    if(param_v.size() != size_t(n_param)) 
        throw string("Wrong number of parameters, expected ") + to_string(param_v.size()) + " but got " + 
            to_string(n_param);
    copy(begin(param_v), end(param_v), param);
    return 0;
} catch(const string& s) {
    fprintf(stderr, "ERROR: %s\n", s.c_str());
    return 1;
} catch(...) {
    return 1;
}

int get_param_size(int* n_param, DerivEngine* engine, const char* node_name) try {
    *n_param = engine->get(string(node_name)).computation->get_param().size();
    return 0;
} catch(const string& s) {
    fprintf(stderr, "ERROR: %s\n", s.c_str());
    return 1;
} catch(...) {
    return 1;
}

int get_param_deriv(int n_param, float* deriv, DerivEngine* engine, const char* node_name) try {
#ifdef PARAM_DERIV
    auto deriv_v = engine->get(string(node_name)).computation->get_param_deriv();
    if(deriv_v.size() != size_t(n_param)) 
        throw string("Wrong number of parameters, expected ") + to_string(deriv_v.size()) + " but got " + 
            to_string(n_param);
    copy(begin(deriv_v), end(deriv_v), deriv);
    return 0;
#else
    return -1;
#endif
} catch(const string& s) {
    fprintf(stderr, "ERROR: %s\n", s.c_str());
    return 1;
} catch(...) {
    return 1;
}


int get_sens(int n_output, float* sens, DerivEngine* engine, const char* node_name) try {
    auto& dc = engine->get_computation<DerivComputation&>(string(node_name));

    if(dc.potential_term) {
        if(n_output!=1) throw string("wrong size for potential node");
        *sens = 1.f;  // the total potential's sensitivity to one of its terms
    } else {
        auto& c = dynamic_cast<CoordNode&>(dc);
        if(n_output != c.n_elem*c.elem_width) throw string("wrong size for CoordNode");
        VecArray a = c.sens;
        for(int ne: range(c.n_elem))
            for(int d: range(c.elem_width))
                sens[ne*c.elem_width + d] = a(d,ne);
    }
    return 0;
} catch(const string& s) {
    fprintf(stderr, "ERROR: %s\n", s.c_str());
    return 1;
} catch(...) {
    return 1;
}


int get_output(int n_output, float* output, DerivEngine* engine, const char* node_name) try {
    auto& dc = engine->get_computation<DerivComputation&>(string(node_name));

    if(dc.potential_term) {
        auto& p = dynamic_cast<PotentialNode&>(dc);
        if(n_output!=1) throw string("wrong size for potential node");
        *output = p.potential;
    } else {
        auto& c = dynamic_cast<CoordNode&>(dc);
        if(n_output != c.n_elem*c.elem_width) throw string("wrong size for CoordNode");
        VecArray a = c.output;
        for(int ne: range(c.n_elem))
            for(int d: range(c.elem_width))
                output[ne*c.elem_width + d] = a(d,ne);
    }
    return 0;
} catch(const string& s) {
    fprintf(stderr, "ERROR: %s\n", s.c_str());
    return 1;
} catch(...) {
    return 1;
}


int get_output_dims(int* n_elem, int* elem_width, DerivEngine* engine, const char* node_name) try {
    auto& dc = engine->get_computation<DerivComputation&>(string(node_name));

    if(dc.potential_term) {
        *n_elem = 1;
        *elem_width = 1;
    } else {
        auto& c = dynamic_cast<CoordNode&>(dc);
        *n_elem = c.n_elem;
        *elem_width = c.elem_width;
    }
    return 0;
} catch(const string& s) {
    fprintf(stderr, "ERROR: %s\n", s.c_str());
    return 1;
} catch(...) {
    return 1;
}


int get_value_by_name(int n_output, float* output, DerivEngine* engine, const char* node_name, const char* log_name)
try {
    auto& dc = engine->get_computation<DerivComputation&>(string(node_name));
    auto value = dc.get_value_by_name(log_name);
    if(n_output != int(value.size()))
        throw string("expected size (") + to_string(n_output) +
            " elements) inconsistent with actual size (" + to_string(value.size()) + ")";
    copy(begin(value),end(value), output);
    return 0;
} catch(const string& s) {
    fprintf(stderr, "ERROR: %s\n", s.c_str());
    return 1;
} catch(...) {
    return 1;
}
