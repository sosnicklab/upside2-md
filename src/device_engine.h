#ifndef DEVICE_ENGINE_H
#define DEVICE_ENGINE_H

// The device engine: the potential graph of every system of a run on one device (device.h), built from the CPU
// engines' constructed nodes, so it reads the same parameters and spline coefficients.  The systems' elements are
// concatenated node by node and every index is offset into its parent's concatenated array, so one launch per node
// evaluates every system.  The nodes run in the CPU engine's order (DerivEngine::compute).

#include "deriv_engine.h"
#include "device.h"
#include <map>

struct OrnsteinUhlenbeckThermostat;

//! One CoordNode's output and sensitivity for every system, the systems' elements one after another
struct DeviceCoord {
    int elem_width;
    int row_width;                 // VecArrayStorage's: element i, component d at i*row_width+d
    int n_elem;                    // over all systems
    std::vector<int> offset;       // the first element of each system, then n_elem
    DeviceArray<float> output;
    DeviceAccum sens;

    DeviceCoord(int elem_width_, int row_width_, const std::vector<int>& n_elem_of_system);
    VecArray value() const {return VecArray(output.ptr, row_width);}  // for reading: a device lambda's copy is const
    std::vector<int> system_of_element() const;
};

//! The device implementation of one node, made by its CPU node's make_device
struct DeviceNode {
    virtual ~DeviceNode() {}
    //! A potential node adds its energy into the graph's per-system potential in PotentialAndDerivMode, and, as on
    //! the CPU, its parents' sensitivities in every mode
    virtual void compute_value(ComputeMode mode) = 0;
    virtual void propagate_deriv() {}
};

//! The graph of every system on the device; every system must have the same nodes
struct DeviceGraph {
    int n_system;
    std::vector<DerivEngine*> engines;
    std::vector<std::unique_ptr<DeviceCoord>> coords;  // per node; null for a potential node
    std::vector<std::unique_ptr<DeviceNode>>  nodes;   // per node; null for pos
    DeviceAccum potential;                             // per node and system, node i system s at i*n_system+s
    DeviceArray<int> error;                            // a failed precondition's code (add_check), 0 if none

    explicit DeviceGraph(const std::vector<DerivEngine*>& engines);

    //! The DeviceCoord of a CoordNode of any of the systems
    DeviceCoord& coord(const CoordNode& cpu_node);
    //! The per-system potential of a potential node of any of the systems, indexed by system
    AccumView potential_of(const DerivComputation& cpu_node);
    //! The concatenation over the systems of f(ns), each system's vector
    template <typename T, typename F>
    std::vector<T> concat(F f) {
        std::vector<T> ret;
        for(int ns=0; ns<n_system; ++ns)
            for(const T& x: f(ns)) ret.push_back(x);
        return ret;
    }
    //! The concatenation of index(ns), each system's indices into the node indexed, offset by the system's first
    //! element there; a negative index, a dummy, stays as it is
    template <typename F>
    std::vector<int> offset_index(const CoordNode& indexed, F index) {
        auto& c = coord(indexed);
        std::vector<int> ret;
        for(int ns=0; ns<n_system; ++ns)
            for(int i: index(ns)) ret.push_back(i<0 ? i : i + c.offset[ns]);
        return ret;
    }

    void compute(ComputeMode mode);
    void compute(ComputeMode mode, int integrator_level);
    //! The code a kernel writes into error when the condition message describes fails
    int add_check(const std::string& message) {checks.push_back(message); return int(checks.size());}
    //! Throws the message of a failed check; called after each integration cycle and evaluation
    void check_errors() const;

    //! Each system's potential, its nodes' potentials summed in the order of DerivEngine::compute
    std::vector<float> get_potential() const;
    //! A potential node's potential in each system
    std::vector<float> get_potential(int node) const;

    private:
        std::map<const DerivComputation*, int> node_index;  // the nodes of every system
        std::vector<std::string> checks;
        void run(ComputeMode mode, int integrator_level, bool all_levels);
};

//! Molecular dynamics of every system on the device: the integrators and thermostat of DerivEngine,
//! OrnsteinUhlenbeckThermostat and martini_brownian, with masses, fixed and z-fixed atoms.  The rigid groups of the
//! hybrid's preparation stages and the barostat stay on the CPU.
struct DeviceDynamics {
    DeviceGraph graph;
    DeviceArray<float> mom;          // the systems' momenta, the layout of pos
    DeviceArray<float> lag;          // the half kick by which the stored momentum trails the one at the coordinates
    DeviceArray<int> atom_system;    // each atom's system
    DeviceArray<int> atom_local;     // each atom's index in its system, which keys its random streams
    DeviceArray<float> mass;         // per atom, empty without masses
    DeviceArray<unsigned char> fixed, z_fixed;  // per atom, as martini_fix_rigid holds them (set_fixed_atoms)
    std::vector<std::vector<int>> fixed_atoms, z_fixed_atoms;  // of each system, the sets the masks hold
    // the thermostat: each system's random seed, invocation count, momentum and noise scales, temperature and step,
    // and per atom whether it is the g-JF bath's (skipped) and its own timescale (empty, or one per atom)
    DeviceArray<uint32_t> thermostat_seed;
    DeviceArray<unsigned long long> thermostat_invocation;
    DeviceArray<float> thermostat_mom_scale, thermostat_noise_scale, thermostat_temp, thermostat_delta_t;
    DeviceArray<unsigned char> stochastic;
    DeviceArray<float> atom_timescale;
    // the g-JF integrator (martini_brownian): each atom's listing and friction, each system's seed, reference
    // temperature and, per step, temperature and invocation count
    bool brownian;
    int n_inner_steps;
    float brownian_time_step;
    DeviceArray<unsigned char> listed;
    DeviceArray<float> friction, force_n;
    DeviceArray<uint32_t> brownian_seed;
    DeviceArray<float> reference_temperature, brownian_temperature;
    DeviceArray<unsigned long long> brownian_invocation;

    explicit DeviceDynamics(const std::vector<DerivEngine*>& engines);

    void set_state(const std::vector<VecArray>& pos, const std::vector<VecArray>& mom);
    void get_state(const std::vector<VecArray>& pos, const std::vector<VecArray>& mom) const;
    //! the fixed and z-fixed atoms of each system (martini_fix_rigid's, which the PO4 hold changes), uploaded when
    //! they differ from the sets held
    void set_fixed_atoms(const std::vector<std::vector<int>>& fixed_, const std::vector<std::vector<int>>& z_fixed_);

    //! As DerivEngine::integration_cycle; a non-empty thermostats (one per system) takes its step, and each
    //! thermostat advances its invocation count as OrnsteinUhlenbeckThermostat::apply does
    void integration_cycle(float dt, const std::vector<OrnsteinUhlenbeckThermostat*>& thermostats);
    void integration_cycle(float dt, int inner_step, const std::vector<OrnsteinUhlenbeckThermostat*>& thermostats);
    //! main.cpp's g-JF round: the evaluation at the coordinates, then martini_brownian::apply_langevin_step with each
    //! system's temperature, at the systems' invocation count, which it advances in each engine as that does
    void langevin_cycle(float dt, const std::vector<float>& temperature, uint64_t n_invocations);
    //! OrnsteinUhlenbeckThermostat::apply, one thermostat per system, which advances each one's invocation count.
    //! at_coordinates: within a leapfrog cycle (integration_cycle); otherwise on the g-JF momenta, which are at the
    //! coordinates, for the atoms outside the g-JF bath.
    void thermostat_step(const std::vector<OrnsteinUhlenbeckThermostat*>& thermostats, bool at_coordinates);
    //! martini_fix_rigid::apply_fix_rigid_md without rigid groups: fixed atoms at rest, z-fixed ones along z
    void zero_fixed_momenta();

    // a step of the cycles; public because a CUDA device lambda cannot sit in a private member function
    void leapfrog_step(float kick, float drift);
};

#endif
