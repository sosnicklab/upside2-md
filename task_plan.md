# Plan: Align UPSIDE-MARTINI with CHARMM-GUI Bilayer Workflow

This plan outlines the implementation of a multi-stage simulation workflow to replicate the GROMACS MARTINI protocol for lipid bilayers (Steps 6.0 through 7.0).

---

## 1. Implement Advanced Pressure Coupling (C++)
The current barostat is limited to the Berendsen algorithm. Production runs for bilayers require the Parrinello-Rahman barostat to ensure correct volume fluctuations and ensemble properties.
- [ ] **Modify `box.cpp`**:
    - Add `ParrinelloRahman` to the `BarostatType` enumeration.
    - Implement the **Parrinello-Rahman** barostat by replicating the standard algorithm used in GROMACS or LAMMPS.
    - Ensure `semi_isotropic` coupling is supported to scale $X/Y$ (lateral) together and $Z$ (normal) independently.
    - Use a default `compressibility` of $3 \times 10^{-4}$ bar⁻¹.

## 2. Implement Position Restraint Node (C++)
This is required to support lipid headgroup restraints (`BILAYER_LIPIDHEAD_FC`) during equilibration stages.
- [ ] **Modify `martini.cpp`**:
    - Create a `PositionRestraint` potential node.
    - Apply a harmonic penalty $V = \frac{1}{2} k (r - r_{ref})^2$ to specified atom indices.
    - Use initial coordinates from `/input/pos` as the reference positions $r_{ref}$.

## 3. Update Preparation Script (Python)
- [ ] **Modify `prepare_martini.py`**:
    - **Lipid Head Selection**: Automatically identify `PO4` and `GL1` beads in DOPC for restraints.
    - **Step 6.0 Soft-core**: Enable soft-core potential handling via `UPSIDE_SOFTEN_LJ` and `UPSIDE_SOFTEN_COULOMB`.
    - **HDF5 Stage Map**: Create an 8-stage sequence in the HDF5 input:
        - **Steps 6.0-6.1**: Minimization; Soft-core enabled for 6.0.
        - **Step 6.2**: MD; $dt=0.002$; $k_{lipid}=200$; Berendsen barostat.
        - **Steps 6.3-6.6**: MD; Gradually increase $dt$ (0.005 to 0.02) and decrease $k_{lipid}$ (100 to 10).
        - **Step 7.0**: MD Production; $dt=0.02$; $k=0$; **Parrinello-Rahman** barostat.



## 4. Workflow Orchestration (Bash)
- [ ] **Refactor `run_sim_bilayer.sh`**:
    - Iterate through all 8 stages sequentially.
    - Use `set_initial_position.py` to chain the output of one stage to the input of the next.
    - Switch the `UPSIDE_BAROSTAT_TYPE` to Parrinello-Rahman only for Step 7.0.

## 5. Global Physics Settings
- **Electrostatics**: Reaction-field with $\epsilon_r=15$ and 1.1 nm cutoff.
- **Thermostat**: Ornstein-Uhlenbeck (`v-rescale`) targeting 303.15 K.
- **Pressure**: 1.0 bar semi-isotropic coupling.