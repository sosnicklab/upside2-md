#ifndef ENGINE_C_LIBRARY_H
#define ENGINE_C_LIBRARY_H

#include "main.h"

extern "C" {
    struct DerivEngine;
    struct DeviceGraph;

    // defined in main.h
    // int upside_main(int argc, const char* const char*)

    DerivEngine* construct_deriv_engine(int n_atom, const char* potential_file);
    void free_deriv_engine(DerivEngine* engine);

    int evaluate_energy(float* energy, DerivEngine* engine, const float* pos);
    int evaluate_deriv (float* deriv,  DerivEngine* engine, const float* pos);

    int set_param      (int n_param, const  float* param,  DerivEngine* engine, const char* node_name);

    int get_param_deriv(int n_param,  float* deriv,  DerivEngine* engine, const char* node_name);
    int get_param      (int n_param,  float* param,  DerivEngine* engine, const char* node_name);
    int get_param_size (int* n_param,                DerivEngine* engine, const char* node_name);
    int get_output_dims(int* n_elem, int* elem_width, DerivEngine* engine, const char* node_name);
    int get_output     (int n_output, float* output, DerivEngine* engine, const char* node_name);
    int get_sens       (int n_output, float* output, DerivEngine* engine, const char* node_name);

    int get_value_by_name(int n_output, float* output, DerivEngine* engine,
            const char* node_name, const char* log_name);

    // The device engine (device_engine.h) of the systems of n_engine engines.  Positions and derivatives are
    // (n_atom,3) per system, system after system; the potential is one per system.
    DeviceGraph* construct_device_graph(int n_engine, DerivEngine** engines);
    void free_device_graph(DeviceGraph* graph);
    int device_describe(int n_char, char* description);
    int device_evaluate(float* potential, float* deriv, DeviceGraph* graph, const float* pos);
    int device_get_output(int n_output, float* output, DeviceGraph* graph, int system, const char* node_name);
    int device_get_sens  (int n_output, float* sens,   DeviceGraph* graph, int system, const char* node_name);
}

#endif
