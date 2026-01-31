# Claude Code Plan: NPT/NVT Workflow Correction & Analysis Filtering

This plan outlines the steps to fix the multi-stage simulation workflow (ensuring NPT actually runs), adding necessary debug logging for pressure coupling, and refining analysis to exclude coarse-grained MARTINI particles.

---

## 1. Workflow Correction (run_sim_bilayer.sh)
The current script may be running NVT twice instead of NPT then NVT.
- **Verify NPT Activation:** Check if `UPSIDE_NPT_ENABLE=1` is sufficient for the `upside` binary or if a flag (e.g., `--npt`) must be added to `CMD_NPT`.
- **Sequential Continuity:** Ensure `set_initial_position.py` correctly passes the *final* coordinates of the previous stage to the next stage's input file.
- **VTF Separation:** Ensure Stage 3 (NPT) and Stage 4 (NVT) produce distinct `.vtf` files representing their respective trajectories.

---

## 2. Source Code Implementation

### Task A: Enhance NPT Print Output
Add high-verbosity debug information to the NPT simulation logging to monitor box stability.
- **Variables:** Current box dimensions ($L_x, L_y, L_z$) and instantaneous pressure.
- **Location:** Locate the barostat application logic or the main integration loop in the C++/source files.
- **Output:** Ensure these values are printed to `stdout` at every `frame-interval`.

### Task B: Filter MARTINI Particles from Analysis
Exclude MARTINI particles from $R_g$ and hbond calculations so metrics reflect only all-atom components.
- **Identification:** Use particle metadata (names/masses) to identify "MARTINI" or "CG" beads.
- **Hbond Filtering:** Update donor/acceptor logic to skip MARTINI particles.
- **Radius of Gyration ($R_g$):** Modify the mass-weighted calculation:
  $$R_g = \sqrt{\frac{1}{M} \sum_{i \in \{non-MARTINI\}} m_i (\mathbf{r}_i - \mathbf{R}_{cm})^2}$$
  *Note: $M$ and $\mathbf{R}_{cm}$ must also be calculated using only the filtered subset.*

---

## 3. Verification & Quality Control
- **Stage Validation:** Run the shell script and check `logs/npt_equilibration.log`. If the box dimensions do not change over time, the NPT stage is not functioning.
- **Dry Run:** Execute a 100-step test to verify new debug strings appear in the logs.
- **Analysis Logic Check:** Compare $R_g$ output against a previous run. The new value should differ as it now excludes the MARTINI environment.

---

## 4. Execution Command for Claude Code
*Paste this into your terminal:*

> "Act as an expert MD developer. 
> 1. Fix `example/16.MARTINI/run_sim_bilayer.sh` to ensure Stage 3 actually runs NPT (check if environment variables or flags are needed) and that Stage 3 and 4 generate distinct VTF trajectories. 
> 2. Find the logging/barostat logic in the source code and add debug prints for box dimensions and pressure. 
> 3. Modify the hbond and Rg calculation functions to filter out MARTINI/CG particles, ensuring Rg uses a filtered center of mass and total mass. 
> Show me the diffs before applying."