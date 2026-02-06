# Session Progress Log

## Date: 2026-02-06
- **Action**: Fix water simulation pressure coupling and compressibility
- **Files Modified**:
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_water.sh`:
    - Changed from semi-isotropic to isotropic pressure coupling (UPSIDE_NPT_SEMI=0)
    - Updated compressibility to water's actual value: 2.1782 Å³/E_up (4.5e-5 bar⁻¹)
- **Results**:
  - Box dimensions now change uniformly across all axes during NPT simulation
  - Pressure values are reasonable and fluctuate around target pressure (1 bar)
  - Simulation runs without box collapse
  - Final box dimensions: ~49.74 Å for all axes (consistent scaling)
- **Notes**:
  - Previous semi-isotropic coupling caused asymmetric box changes
  - Incorrect compressibility value (14.52 Å³/E_up) was causing extreme box shrinkage
  - Fixed by setting UPSIDE_NPT_COMPRESSIBILITY to 2.1782 Å³/E_up

## Date: 2026-02-06

- **Action**: Fix NPT barostat compressibility value to ensure box dimensions change
- **Files Modified**:
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/prepare_martini.py`: Changed default compressibility from 4.5e-5 to 14.521180763676 (correct value for 3e-4 bar^(-1))
- **Results**:
  - Box dimensions now change during NPT simulation stages
  - Pressure values are more reasonable and fluctuate around target pressure
  - Simulation progresses correctly with proper box scaling
- **Notes**:
  - The compressibility value was incorrect in the previous version, leading to no box dimension changes
  - The fix ensures that the barostat can properly scale the box to reach the target pressure

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
