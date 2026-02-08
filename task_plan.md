# Task: Implement Ewald Summation for MARTINI 3.0 Potentials

## 1. Objective
Implement a standard Ewald summation for periodic long-range electrostatics. This involves implementing reciprocal space logic in `src/box.cpp`, updating the real-space potential in `src/martini.cpp`, ensuring the `src/main.cpp` simulation loop handles the updated force calls, and enabling the feature via `example/16.MARTINI/run_sim_bilayer.sh`.

## 2. Technical Context & Unit Constraints
The implementation must adhere to the specific unit system already established in the MARTINI 3.0 port:

* **Length Units**: Angstroms ($\text{\AA}$).
* **Energy Units**: $E_{up}$ ($1 \, E_{up} \approx 2.915 \, \text{kJ/mol}$).
* **Coulomb Constant ($k_C$)**: $31.775347952181$.
    * This constant incorporates the MARTINI relative permittivity ($\epsilon_r = 15$) and unit conversions.
* **Charge Units**: The current implementation in `src/martini.cpp` expects charges in elementary charge units ($e$).
    * **Note**: The workflow should set input charges to $1.0$ (or $-1.0$) as the conversion to $E_{up}$ is handled internally by $k_C$ during potential/force calculation.
* **Force Units**: $E_{up}/\text{\AA}$.

---

## 3. Implementation Requirements

### A. Modifications to `src/box.h` & `src/box.cpp`
* **Data Structure**: Add an `EwaldState` struct to the `simulation_box::npt` namespace to store the splitting parameter $\alpha$, $k$-vectors, and structural factors.
* **Reciprocal Space Function**: Implement `compute_ewald_reciprocal` to calculate:
    1.  **Reciprocal Potential**: $V_{recip} = \frac{k_C}{2\pi V} \sum_{\mathbf{k} \neq 0} \frac{4\pi^2}{k^2} e^{-k^2/4\alpha^2} |\rho(\mathbf{k})|^2$.
    2.  **Forces**: Apply the reciprocal gradient directly to the `engine.pos->sens` array.
    3.  **Virial contribution**: Ensure the reciprocal sum contributes to the pressure tensor for NPT box scaling.
* **Self-Interaction Correction**: Implement the constant term $V_{self} = -k_C \frac{\alpha}{\sqrt{\pi}} \sum_i q_i^2$.

### B. Modifications to `src/martini.cpp`
* **Real-Space Screening**: Update `MartiniPotential::compute_value` to utilize the screened Coulomb potential: $V_{real}(r) = \frac{k_C q_i q_j \text{erfc}(\alpha r)}{r}$.
* **Configuration Reading**: Update the initialization to read `ewald_alpha` and `ewald_enabled` from the HDF5 path `/input/potential/martini_potential`.

### C. Modifications to `src/main.cpp`
* **Namespace Hooks**: Update the `simulation_box::npt` namespace declarations to include `initialize_ewald` and `apply_ewald_reciprocal`.
* **Simulation Loop**: Trigger the Ewald reciprocal calculation during the force evaluation step, ensuring it correctly accounts for box dimensions if they were modified by the barostat.
* **Logging**: Update potential energy logging to include reciprocal and self-correction terms for accurate energy monitoring.

### D. Modifications to `example/16.MARTINI/run_sim_bilayer.sh`
* **Charge Configuration**: Ensure the preparation scripts maintain elementary charges (e.g., $+1.0, -1.0$) as the code handles the $k_C$ conversion.
* **Environment Variables**: Add `UPSIDE_EWALD_ENABLE` and `UPSIDE_EWALD_ALPHA` to the user configuration.
* **Workflow Integration**: Ensure the `npt_prod` stage writes the Ewald configuration to the HDF5 input file.

---

## 4. Constraints
* **Direct Summation**: Use a direct $k$-space summation ($O(N \cdot N_k)$) for this initial implementation to maintain simplicity and avoid adding external FFT dependencies.
* **Barostat Awareness**: Reciprocal sums must be re-evaluated whenever `apply_semi_isotropic_scaling` is called in `src/box.cpp`.
* **No Emoji**: Do not use emojis in code comments or log outputs.