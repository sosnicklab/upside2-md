# Session Progress Log

## Date: 2026-02-04

- **Action**: Fix softened to hard particle switching issue in run_sim_bilayer.sh
- **Files Modified**:
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_bilayer.sh`: Added HDF5 attribute modification code, updated stage descriptions, added progressive softening stages
- **Results**:
  - Simulation now properly switches from softened to hard particles using progressive softening reduction
  - Potential values are stable in production stage (-4721.86 to -4762.31)
- **Notes**:
  - Implemented two NPT equilibration stages with reduced softening before production
  - Added code to directly modify HDF5 file attributes to ensure softening parameters are correctly updated

## Date: 2026-02-03

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
