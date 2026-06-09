**Context & Role:**
You are an expert computational biophysicist and simulation software developer. Your task is to debug, develop, and iterate on a hybrid simulation interface combining the UPSIDE (vector-based side chain) and dry-MARTINI force fields.

**Current State & Objective:**
We currently use a dry-MARTINI `.itp` file to calculate interactions between an UPSIDE particle (typically a side chain, SC) and a dry-MARTINI particle. The current implementation yields imperfect results that seem to rely on "twisting" parameters, but it is partially functional. The current output data and simulation results are available for your review at `example/16.MARTINI/outputs/*`.

Your core objective is to take this a step further: train a force field (FF) that uses a single vector particle to represent a DOPC lipid (Coarse-Grained Lipid, CGL). Because the current model somewhat works, you must prioritize small, incremental modifications to the existing implementation based on the data in the outputs directory. However, every modification must remain strictly physical. Do not compromise thermodynamic or structural accuracy for convenience.

**Known Defects to Address:**
Based on the current outputs, the model exhibits two specific physical defects that you must resolve:

1. **Bilayer-Protein Gap:** There is an unphysical spatial gap between the CGL bilayer and the embedded protein. Your modifications to the CGL-SC interaction parameters must facilitate proper solvation and packing of the CGL particles against the protein side chains.
2. **Orientational Artifacts:** The CGL orientations are problematic both in the immediate vicinity of the protein and at the periodic box boundaries. You must ensure the orientational dependencies (especially when interacting with SCs or evaluating CGL-CGL interactions across boundaries) maintain a stable, physical bilayer structure without artificial flipping or distortion.

**Mathematical Framework & Required Modifications:**
For SC-SC interactions, UPSIDE uses the following potential:


$$V=\kappa(V_{radial}(r_{12})+Ang_1(-\mathbf{n}_1 \cdot \mathbf{n}_{12}) \cdot Ang_2(\mathbf{n}_2 \cdot \mathbf{n}_{12}) \cdot V_{angular}(r_{12}))$$

You must modify and implement this potential for our specific topological scenarios:

1. **SC-Particle Interaction:** Standard dry-MARTINI particles lack directionality. Modify the potential to require fewer parameters (removing orientation dependence on the particle side).
2. **CGL-CGL and CGL-SC Interactions:** The single-vector CGL representation is not purely radial. Formulate the extra parameters needed to capture the anisotropic nature of the lipid vector, specifically targeting the packing and orientational defects mentioned above.

**Strict Implementation Constraints:**
To achieve the mathematically proper result, you must strictly adhere to the following constraints throughout development and testing:

* **No Parameter Twisting:** The physics must dictate the stability, not arbitrary scaling factors.
* **No Force or Energy Capping:** The potential definitions must naturally avoid singularities or exploding forces.
* **No Thermodynamic Cheating:** It is strictly forbidden to lower the temperature, drastically reduce the timestep ($dt$), or alter thermostat/barostat controls to make the system artificially slow just to bypass stability tests. Stability must be a consequence of the correct potential geometry and parameters.
* **No Additional Orientation Potentials:** Do not add ad-hoc orientation potentials on CGL interactions; the core modified UPSIDE potential must handle it natively.
* **Decomposition Strategy:** Utilize Singular Value Decomposition (SVD) or similar algorithms to decompose the potentials into sums and products of single-parameter potentials. UPSIDE assumes that single-parameter potentials are smooth, which will drastically reduce the data points required for training.

**Execution Plan & Testing Workflow:**
Please execute the development and testing in the following specific phases:

* **Phase 1: Pure CGL Bilayer Debugging.** Debug and iterate on a CGL-only bilayer model. Do not move forward until this is physically correct. Criteria for success: The bilayer is stable without parameter twisting or thermodynamic cheating, boundary orientational artifacts are resolved, and the CGL orientations are physically realistic for a DOPC membrane.
* **Phase 2: Full System Integration.** Once the pure CGL bilayer works perfectly, apply the method to SC-CGL, CGL-particle, and SC-particle interactions. You must specifically verify that the bilayer-protein gap and near-protein orientational defects are resolved. Run and verify the results through the following specific workflow scripts:
* `example/16.MARTINI/run_sim_1afo.sh`
* `example/16.MARTINI/run_sim_1rkl.sh`
* `example/16.MARTINI/run_sim_1rkl_full.sh`
* `example/16.MARTINI/run_sim_1afo_full.sh`



**Suggested Training Methodologies:**

* **Directional Restraints:** Apply a directional restraint to two interacting molecules, allowing movement only along the line connecting their centers of mass (COM), and then relax the system. Use this to filter out unwanted orthogonal movements during training.
* **Time Integration:** Maintain the same timestep ($dt$) and core integration settings between the UPSIDE core and dry-MARTINI particles if possible. Provide a rationale if they must differ, keeping in mind the anti-cheating constraint above.