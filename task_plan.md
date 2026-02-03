# Plan: Align UPSIDE-MARTINI with CHARMM-GUI Bilayer Workflow

This plan outlines the implementation of a multi-stage simulation workflow to mirror the GROMACS MARTINI protocol for lipid bilayers (Steps 6.0 through 7.0), focusing on advanced pressure coupling and stage-wise restraint management.

## 1. Implement Advanced Pressure Coupling (C++)
The current barostat is limited to the Berendsen algorithm. Production runs for bilayers require the Parrinello-Rahman barostat to ensure correct volume fluctuations and ensemble properties.
- [ ] **Modify `box.cpp`**:
    - Add `ParrinelloRahman` to the `BarostatType` enumeration.
    - Implement box dynamical variables (velocities $v_{box\_xy}$, $v_{box\_z}$) to track box expansion/contraction rates.
    - Implement the Parrinello-Rahman integration logic: $v_{box}(t + dt) = v_{box}(t) + \frac{dt}{W} (P_{inst} - P_{target})$.
    - Set the default `compressibility` to $3 \times 10^{-4}$ bar⁻¹ to match bilayer standards.
    - Ensure `semi_isotropic` coupling scales $X/Y$ (lateral) together and $Z$ (normal) independently to maintain bilayer integrity.

## 2. Implement Position Restraint Node (C++)
To support the lipid headgroup restraints (`BILAYER_LIPIDHEAD_FC`) used in equilibration steps.
- [ ] **Modify `martini.cpp`**:
    - Create a `PositionRestraint` potential node.
    - The node should read a list of atom indices and apply a harmonic penalty: $V = \frac{1}{2} k (r - r_{ref})^2$.
    - Use the initial coordinates from `/input/pos` as the reference positions ($r_{ref}$).

## 3. Update Preparation Script (Python)
- [ ] **Modify `prepare_martini.py`**:
    - **Identify Lipid Heads**: Automatically select `PO4` and `GL1` beads in DOPC for the `BILAYER_LIPIDHEAD_FC` restraints.
    - **Step 6.0 Logic**: Configure a "Soft-core" potential mode using `UPSIDE_SOFTEN_LJ` and `UPSIDE_SOFTEN_COULOMB` to handle initial overlaps.
    - **HDF5 Stage Map**: Write a configuration group `/input/stages` containing the exact parameters for the 8-stage sequence:
        - **Steps 6.0-6.1**: Minimization; Soft-core enabled for 6.0; $dt$ N/A.
        - **Step 6.2**: MD; $dt=0.002$; $k_{lipid}=200$; Berendsen barostat.
        - **Steps 6.3-6.6**: MD; Gradually increase $dt$ (0.005 to 0.02) and decrease $k_{lipid}$ (100 to 10).
        - **Step 7.0**: MD Production; $dt=0.02$; $k=0$; **Parrinello-Rahman** barostat.



## 4. Workflow Orchestration (Bash)
- [ ] **Refactor `run_sim_bilayer.sh`**:
    - Iterate through the sequence of 8 checkpoint files.
    - Update `set_initial_position.py` to ensure box dimensions and velocities are passed between stages.
    - Switch the `UPSIDE_BAROSTAT_TYPE` from Berendsen (Equilibration) to Parrinello-Rahman (Production) at Step 7.0.

## 5. Global Physics Settings
- **Electrostatics**: Enforce reaction-field with $\epsilon_r=15$ and a 1.1 nm cutoff for all stages.
- **Thermostat**: Use the Ornstein-Uhlenbeck (`v-rescale`) thermostat targeting 303.15 K.
- **Pressure**: Target 1.0 bar semi-isotropically.
