#include "deriv_engine.h"
#include "timing.h"
#include "Plumed.h"
#include <iostream>
#include <cstring>
#include <array>
#include "state_logger.h"
#include "h5_support.h"

using namespace h5;
using namespace std;

struct PlumedForce : public PotentialNode
{
    plumed plumedmain;
    CoordNode& pos;
    float dt;
    float kbT;
    int n_atoms;
    bool just_print;
    string plumedFile;
    int loopsize;
    int frame_interval;

    double virial[9] = {0.0f}; // initialize to all zeros array as we don't care about pressure

    std::vector<float> masses;
    bool hasInitialized = false;
    int step;


    PlumedForce(hid_t grp, CoordNode& pos_): PotentialNode(), pos(pos_),
                  dt(read_attribute<float>(grp, ".", "dt")),
                  kbT(read_attribute<float>(grp, ".", "kbT")),
                  n_atoms(read_attribute<int>(grp, ".", "n_atoms")),
                  just_print((bool)read_attribute<int>(grp, ".", "just_print"))
    {
        loopsize = round(1.0/(3.0*dt))*3;

        int n_line = get_dset_size(1, grp, "plumedFile")[0];
        vector<string> plumedInput(n_line);
        traverse_string_dset<1>(grp, "plumedFile", [&](size_t nr, string s) { plumedInput[nr] = s;});
        plumedFile = plumedInput[0];

        plumedmain=plumed_create();                 // Create the plumed object
        
        // Calls to pass data to plumed
        int real_precision=4;
        plumed_cmd(plumedmain, "setRealPrecision",&real_precision);   // Pass a pointer to an integer containing the size of a real number (4 or 8)
        float conversionUnits = 1.0;
        plumed_cmd(plumedmain, "setMDEnergyUnits",&conversionUnits);  // Pass a pointer to the conversion factor between the energy unit used in your code and kJ mol-1
        plumed_cmd(plumedmain, "setMDTimeUnits", &conversionUnits);   // Pass a pointer to the conversion factor between the time unit used in your code and ps
        float lengthUnits=0.1;
        plumed_cmd(plumedmain, "setMDLengthUnits",&lengthUnits);      // Pass a pointer to the conversion factor between the length unit used in your code and nm 
        plumed_cmd(plumedmain, "setMDEngine","upside");               // Pass the name of your md engine to plumed (now it is just a label) 
        plumed_cmd(plumedmain, "setPlumedDat", plumedFile.c_str());   // Pass the name of the plumed input file from the md code to plumed
        plumed_cmd(plumedmain, "setNatoms",&n_atoms);                 // Pass a pointer to the number of atoms in the system to plumed
        plumed_cmd(plumedmain, "setTimestep", &conversionUnits);       // Pass a pointer to the molecular dynamics timestep to plumed
        plumed_cmd(plumedmain, "setKbT", &kbT);                        // Pointer to a real containing the value of kbT

        plumed_cmd(plumedmain, "init", NULL);

        masses.resize(n_atoms);
        for (int i = 0; i < n_atoms; i++) {
            masses[i] = 1.0;
        }
        step = 0;
        hasInitialized = true;
    }

    ~PlumedForce() {
        if (hasInitialized)
            plumed_finalize(plumedmain);
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("plumedforce")); 
        potential = 0.0;
        VecArray posc = pos.output;
        VecArray pos_sens = pos.sens;

        if (step%loopsize == 0) {
            frame_interval = step/loopsize;
        }
        else {
            frame_interval = -1;
        }
        plumed_cmd(plumedmain, "setStep", &frame_interval);
        plumed_cmd(plumedmain, "setMasses", &masses[0]);

        const int pfsize = n_atoms*3;  

        float *p = new float[pfsize];
        float *f = new float[pfsize];

        for (int i=0; i<n_atoms; i++) {
            auto x = load_vec<3>(posc, i);
            for (int d=0;d<3;d++) {
                p[d+i*3] = x[d];
                f[d+i*3] = 0.0;
            }
        }

        plumed_cmd(plumedmain, "setPositions", p);
        plumed_cmd(plumedmain, "setForces", f);
        plumed_cmd(plumedmain, "setVirial", &virial);

        // Calculate the forces and energy.
        plumed_cmd(plumedmain, "prepareCalc", NULL);
        plumed_cmd(plumedmain, "performCalcNoUpdate", NULL);
        if (mode==DerivMode) {
            step++;
            if (frame_interval >= 0) 
                plumed_cmd(plumedmain, "update", NULL);
        }

        if (!just_print) {
            if (mode==PotentialAndDerivMode)
                plumed_cmd(plumedmain, "getBias", &potential);

            for (int i=0; i<n_atoms; i++) {
                for (int d=0;d<3;d++) {
                     pos_sens.x[d+i*4] -= f[d+i*3];
                }
            }
        }

        delete p;
        delete f;
    }
};
static RegisterNodeType<PlumedForce,1> pos_spring_node("plumedforce");

