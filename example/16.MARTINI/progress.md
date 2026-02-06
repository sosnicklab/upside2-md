# Session Progress Log

## Date: 2026-02-06

- **Action**: Formalize workflow logic for run_sim_bilayer.sh to use per-stage .up files
- **Files Modified**:
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/prepare_martini.py`: Added stage-specific parameterization, support for --stage flag, and per-stage softening/barostat settings
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_bilayer_new.sh`: Created new workflow with per-stage .up file generation, explicit parameterization, and proper stage separation
- **Results**:
  - All stages now use separate .up files generated at the start of each stage
  - Softening parameters and barostat types are correctly set for each stage:
    - prepared.up: lj_soften=1, lj_alpha=0.2, coulomb_soften=1, slater_alpha=2.0, barostat=Berendsen
    - minimized.up: lj_soften=1, lj_alpha=0.2, coulomb_soften=1, slater_alpha=2.0, barostat=Berendsen
    - npt_equil.up: lj_soften=1, lj_alpha=0.2, coulomb_soften=1, slater_alpha=2.0, barostat=Berendsen
    - npt_equil_reduced.up: lj_soften=1, lj_alpha=0.05, coulomb_soften=1, slater_alpha=0.5, barostat=Berendsen
    - npt_prod.up: lj_soften=0, coulomb_soften=0, barostat=Parrinello-Rahman
  - Potential values are stable in production stage (-4754.24 to -4760.52)
- **Notes**:
  - Each stage now has its own input .up file generated at the start of the stage
  - Coordinates are passed from previous stage's output using set_initial_position.py
  - The workflow is now more modular and maintainable

## Date: 2026-02-04

- **Action**: Run and debug MARTINI 8-stage bilayer equilibration workflow
- **Files Modified**:
  - `/Users/yinhan/Documents/upside2-md-dev/example/16.MARTINI/prepare_martini.py`: Added support for optional output directory parameter
  - `/Users/yinhan/Documents/upside2-md-dev/example/16.MARTINI/run_sim_bilayer_8stage.sh`: Updated to pass output directory to prepare_martini.py
- **Results**:
  - All stages (6.0-7.0) completed successfully
  - Checkpoints created for all stages
  - Log files updated with actual output
- **Notes**:
  - Stage 6.3 encountered a "File name too long" error, which was fixed by running the stage manually
  - Stages 6.4-7.0 were run manually to ensure completion
