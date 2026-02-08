# Task: Implement Ewald Summation for MARTINI 3.0 Potentials

## 1. Objective
Implement a standard Ewald summation to handle periodic long-range electrostatics accurately within the current NPT/NVT framework. This involves updating `box.cpp` to handle reciprocal space, modifying `martini.cpp` to use screened real-space potentials, and updating the shell workflow in `run_sim_bilayer.sh`.

## 2. Technical Context & Unit Constraints
The implementation must adhere to the specific unit system defined in the existing MARTINI 3.0 port:

* **Length Units**: Angstroms ($\text{\AA}$).
* **Energy Units**: $E_{up}$ ($1 \, E_{up} \approx 2.915 \, \text{kJ/mol}$).
* **Coulomb Constant ($k_C$)**: $31.775347952181$.
    * This constant already incorporates the MARTINI relative permittivity ($\epsilon_r = 15$).
* **Force Units**: $E_{up}/\text{\AA}$.

---

## 3. Implementation Requirements

### A. Modifications to box.h & box.cpp
* **Data Structure**: Add an `EwaldState` struct to the `simulation_box::npt` namespace to store the splitting parameter $\alpha$, $k$-vectors, and precomputed structural factors.
* **Reciprocal Space Function**: Implement `compute_ewald_reciprocal` to calculate:
    1.  **Reciprocal Potential**: $V_{recip} = \frac{k_C}{2\pi V} \sum_{\mathbf{k} \neq 0} \frac{4\pi^2}{k^2} e^{-k^2/4\alpha^2} |\rho(\mathbf{k})|^2$.
    2.  **Forces**: Apply the reciprocal gradient to `engine.pos->sens`.
    3.  **Virial contribution**: Ensure the reciprocal sum contributes to the pressure tensor for NPT box scaling.
* **Self-Interaction Correction**: Add the constant term $V_{self} = -k_C \frac{\alpha}{\sqrt{\pi}} \sum_i q_i^2$.

### B. Modifications to martini.cpp
* **Real-Space Screening**: Update `MartiniPotential::compute_value` to utilize the screened Coulomb potential: $V_{real}(r) = \frac{k_C q_i q_j \text{erfc}(\alpha r)}{r}$.
* **Configuration Reading**: Update the engine to read `ewald_alpha` and `ewald_enabled` from the HDF5 path `/input/potential/martini_potential`.

### C. Modifications to run_sim_bilayer.sh
* **Environment Variables**: Add `UPSIDE_EWALD_ENABLE` and `UPSIDE_EWALD_ALPHA` to the user configuration section.
* **Workflow Integration**: Ensure the production stage (`npt_prod`) is updated to pass these parameters into the simulation.

---

## 4. Constraints
* **Direct Summation**: Use a direct $k$-space summation for this initial implementation (O(N*N_k)) to maintain simplicity and avoid external FFT dependencies.
* **Barostat Awareness**: The reciprocal space summation must be re-evaluated whenever the box dimensions are updated via `apply_semi_isotropic_scaling` in `box.cpp`.
* **Logging**: Add the Ewald potential contribution to the "potential" logger in `main.cpp`.