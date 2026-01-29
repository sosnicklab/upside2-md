# Claude Code Plan: NPT Simulation Updates & Analysis Filtering

This plan outlines the steps for modifying the NPT simulation output and refining the analysis logic to exclude coarse-grained MARTINI particles.

---

## 1. Context Exploration
Before making changes, identify the core logic within the codebase:
- **Search for Output Logic:** Locate where the simulation status is printed (e.g., `log`, `print`, or `sys.stdout`).
- **Identify Particle Metadata:** Determine how the system distinguishes between all-atom and MARTINI particles (e.g., via atom names, masses, or a `is_martini` flag).
- **Locate Analysis Modules:** Find the specific functions responsible for calculating Hydrogen Bonds (hbond) and the Radius of Gyration ($R_g$).

---

## 2. Implementation Tasks

### Task A: Enhance NPT Print Output
Add high-verbosity debug information to the NPT simulation logging to monitor box stability and pressure fluctuations.
- **Variables to Add:** Current box dimensions (length, width, height) and instantaneous pressure.
- **Location:** The main integration loop or the dedicated `Logger` class.
- **Format:** Ensure the debug info is clearly labeled so it can be parsed or ignored by downstream analysis scripts.

### Task B: Filter MARTINI Particles from Analysis
Exclude MARTINI particles from $R_g$ and hbond calculations to ensure the metrics reflect only the all-atom components of the system.

1. **Hydrogen Bond (hbond) Filtering:**
   - Update the neighbor search or donor/acceptor identification logic.
   - Skip any particle identified as "MARTINI".
2. **Radius of Gyration ($R_g$) Filtering:**
   - Modify the mass-weighted calculation:
     $$R_g = \sqrt{\frac{1}{M} \sum_{i \in \{non-MARTINI\}} m_i (\mathbf{r}_i - \mathbf{R}_{cm})^2}$$
   - Ensure the total mass $M$ and Center of Mass $\mathbf{R}_{cm}$ are also calculated using only the filtered particle set.

---

## 3. Verification & Quality Control
- **Dry Run:** Execute a 100-step simulation to verify that the new debug strings appear in the logs without breaking the formatting.
- **Logic Check:** Compare the $R_g$ output against a previous run. The new $R_g$ should specifically represent the all-atom subset, typically resulting in a different value if MARTINI solvent or lipids were previously included.
- **Unit Test:** If a test suite exists, add a case that confirms the `calculate_hbond` function returns zero when provided with a system containing only MARTINI particles.

---

## 4. Execution Command for Claude Code
*Paste the following into your terminal to start the task:*

> "Act as an expert MD simulation developer. First, find the logging function and add debug info for box dimensions and pressure. Second, modify the hbond and Rg calculation functions to filter out any particles belonging to the MARTINI model. Ensure the Rg calculation uses a filtered center of mass and total mass. Show me the diffs before applying."
