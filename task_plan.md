# Claude Code Plan: MARTINI Workflow Correction & NPT Logging

This plan outlines the steps to fix the multi-stage MARTINI simulation workflow, ensuring the NPT stage correctly utilizes the barostat and provides necessary debug logging for box stability.

---

## 1. Workflow Sequence Correction (run_sim_bilayer.sh)
The script `example/16.MARTINI/run_sim_bilayer.sh` must be updated to follow a strict sequential dependency where each stage uses the result of the previous stage as its starting structure.

**Required Sequence:**
1. **Prepare Input:** Generate the initial topology and coordinates.
2. **Minimization:** Relax the system to remove high-energy clashes.
3. **NPT Equilibration:** Allow the box dimensions to fluctuate to reach target pressure.
4. **NVT Production:** Run production at a fixed volume determined by the NPT stage.
5. **VTF Extraction:** Generate distinct trajectory files for both NPT and NVT stages.

---

## 2. Implementation Tasks

### Task A: Fix Stage Logic & Continuity
- **Verify NPT Activation:** Confirm if `UPSIDE_NPT_ENABLE=1` is correctly recognized by the `upside` binary or if explicit command-line flags (e.g., `--npt`) are required in `CMD_NPT`.
- **Sequential Continuity:** Audit the `cp` and `set_initial_position.py` calls to ensure the final frame of the NPT equilibration is correctly passed as the starting state for NVT production.
- **VTF Separation:** Modify the extraction stage to ensure that Stage 3 (NPT) and Stage 4 (NVT) generate independent `.vtf` files for distinct analysis.

### Task B: Source Code Implementation (NPT Debugging)
To verify the stability of the NPT stage, the simulation must report periodic updates on the system volume.
- **Update Logging:** Enhance the print output during NPT simulations to include box dimensions.
- **Variables to Add:** Current box dimensions ($L_x, L_y, L_z$) and instantaneous pressure.
- **Location:** The main integration loop or the dedicated `Logger` class within the source code.

---

## 3. Verification & Quality Control
- **Logic Check:** Execute the workflow and inspect `npt_equilibration.log`. If box dimensions remain static despite NPT being "enabled," the barostat is not being triggered.
- **Dry Run:** Execute a 100-step simulation to verify that the new debug strings appear in the logs without breaking the formatting.
- **File Integrity:** Confirm that the output VTF files for the NPT and NVT stages contain the unique trajectories of their respective runs.

---

## 4. Execution Command for Claude Code
*Paste the following into your terminal to start the task:*

> "Act as an expert MD simulation developer. 
> 1. Modify `example/16.MARTINI/run_sim_bilayer.sh` to follow the strict sequence: Prepare -> Minimization -> NPT -> NVT -> VTF. 
> 2. Ensure Stage 3 (NPT) actually activates the barostat (check if flags like --npt are needed) and ensure Stage 4 (NVT) uses the NPT final state. 
> 3. Fix the VTF extraction to produce unique files for stages 3 and 4. 
> 4. In the source code, add debug prints for box dimensions and pressure during NPT. 
> Show me the diffs before applying."