# UPSIDE-MARTINI Bilayer Equilibration Implementation

This document describes the implementation of the CHARMM-GUI MARTINI protocol (Steps 6.0-7.0) for lipid bilayer equilibration in UPSIDE.

## Overview

The implementation adds support for:
1. **Parrinello-Rahman barostat** for proper NPT ensemble dynamics
2. **Position restraints** for lipid headgroup restraints during equilibration
3. **8-stage equilibration protocol** matching CHARMM-GUI workflow
4. **Soft-core potentials** for initial minimization

## Implementation Details

### 1. Parrinello-Rahman Barostat (C++)

**Files Modified:**
- `src/box.h`: Added `BarostatType` enum and state variables
- `src/box.cpp`: Implemented Parrinello-Rahman algorithm

**Features:**
- Semi-isotropic coupling (X/Y lateral, Z normal independent)
- Box velocity tracking for extended Lagrangian dynamics
- Damping to prevent oscillations
- Backward compatible with Berendsen barostat

**Configuration (HDF5):**
```
/input/barostat/
  - enable: 1
  - type: 0 (Berendsen) or 1 (Parrinello-Rahman)
  - target_p_xy: target lateral pressure
  - target_p_z: target normal pressure
  - tau_p: barostat time constant
  - compressibility: isothermal compressibility
  - semi_isotropic: 1 for semi-isotropic coupling
  - interval: steps between barostat applications
```

**Environment Variables:**
```bash
export UPSIDE_BAROSTAT_TYPE=1  # 0=Berendsen, 1=Parrinello-Rahman
```

### 2. Position Restraint Node (C++)

**Files Modified:**
- `src/martini.cpp`: Added `PositionRestraint` potential node

**Features:**
- Harmonic penalty: V = 0.5 * k * (r - r_ref)^2
- Per-atom force constants
- Reference positions from initial configuration

**Configuration (HDF5):**
```
/input/potential/restraint_position/
  - restraint_indices: [n_restraints] atom indices
  - ref_pos: [n_restraints x 3] reference positions
  - spring_const: [n_restraints] force constants
```

### 3. Preparation Script Updates (Python)

**Files Modified:**
- `example/16.MARTINI/prepare_martini.py`

**Changes:**
- Added `UPSIDE_BAROSTAT_TYPE` environment variable support
- Barostat type is written to HDF5 configuration
- Displays barostat type in output

**Usage:**
```bash
export UPSIDE_BAROSTAT_TYPE=0  # For equilibration stages
python3 prepare_martini.py bilayer
```

### 4. 8-Stage Workflow Script (Bash)

**Files Created:**
- `example/16.MARTINI/run_sim_bilayer_8stage.sh`

**Workflow Stages:**

| Stage | Type | Timestep | Steps | Barostat | Soft-Core | Notes |
|-------|------|----------|-------|----------|-----------|-------|
| 6.0 | Minimization | - | 500 | - | Yes | Initial relaxation |
| 6.1 | Minimization | - | 500 | - | No | Final minimization |
| 6.2 | MD | 0.002 | 2000 | Berendsen | No | Strong restraints |
| 6.3 | MD | 0.005 | 2000 | Berendsen | No | Gradual equilibration |
| 6.4 | MD | 0.010 | 2000 | Berendsen | No | Gradual equilibration |
| 6.5 | MD | 0.015 | 2000 | Berendsen | No | Gradual equilibration |
| 6.6 | MD | 0.020 | 2000 | Berendsen | No | Final equilibration |
| 7.0 | MD | 0.020 | 5000 | Parrinello-Rahman | No | Production |

**Features:**
- Automatic stage chaining using `set_initial_position.py`
- Soft-core potentials for Stage 6.0
- Gradual timestep increase (0.002 → 0.020 ps)
- Automatic barostat type switching (Berendsen → Parrinello-Rahman)
- Checkpoint files for each stage
- Detailed logging

**Usage:**
```bash
cd example/16.MARTINI
export PDB_ID=bilayer
./run_sim_bilayer_8stage.sh
```

## Usage Examples

### Basic 8-Stage Bilayer Equilibration

```bash
cd example/16.MARTINI
source ../../source.sh
source ../../.venv/bin/activate

# Set system name
export PDB_ID=bilayer

# Run 8-stage workflow
./run_sim_bilayer_8stage.sh
```

### Custom Stage Durations

```bash
# Adjust stage durations (in MD steps)
export MIN_STEPS_60=1000
export MIN_STEPS_61=1000
export MD_STEPS_62=5000
export MD_STEPS_63=5000
export MD_STEPS_64=5000
export MD_STEPS_65=5000
export MD_STEPS_66=5000
export MD_STEPS_70=10000

./run_sim_bilayer_8stage.sh
```

### Custom NPT Settings

```bash
# Adjust NPT parameters (values in UPSIDE units)
# Note: 1 bar = 0.000020659 E_up/Angstrom^3
#       1 atm = 0.000020933215 E_up/Angstrom^3
export UPSIDE_NPT_TARGET_PXY=0.000020659  # 1 bar lateral pressure
export UPSIDE_NPT_TARGET_PZ=0.000020659   # 1 bar normal pressure
export UPSIDE_NPT_TAU=1.0
export UPSIDE_NPT_COMPRESSIBILITY=3e-4
export UPSIDE_NPT_INTERVAL=10

./run_sim_bilayer_8stage.sh
```

## Output Files

The workflow creates the following directory structure:

```
outputs/martini_8stage/
├── checkpoints/
│   ├── bilayer.stage_6.0.up  # After minimization with soft-core
│   ├── bilayer.stage_6.1.up  # After minimization without soft-core
│   ├── bilayer.stage_6.2.up  # After MD dt=0.002
│   ├── bilayer.stage_6.3.up  # After MD dt=0.005
│   ├── bilayer.stage_6.4.up  # After MD dt=0.010
│   ├── bilayer.stage_6.5.up  # After MD dt=0.015
│   ├── bilayer.stage_6.6.up  # After MD dt=0.020
│   └── bilayer.stage_7.0.up  # Production (Parrinello-Rahman)
└── logs/
    ├── stage_6.0.log
    ├── stage_6.1.log
    ├── stage_6.2.log
    ├── stage_6.3.log
    ├── stage_6.4.log
    ├── stage_6.5.log
    ├── stage_6.6.log
    └── stage_7.0.log
```

## Technical Notes

### UPSIDE Unit Conversions

| Quantity | UPSIDE Unit | Standard Equivalent |
|----------|-------------|---------------------|
| Energy | 1 E_up | 2.914952774272 kJ/mol |
| Length | 1 Angstrom | 1 Angstrom |
| Mass | 1 m_up | 12 g/mol |
| Temperature | 1.0 T_up | 350.588235 Kelvin |
| Pressure | 0.000020659 E_up/Å³ | 1 bar |

### Barostat Comparison

**Berendsen:**
- Deterministic pressure coupling
- Fast equilibration
- Does not produce correct NPT ensemble
- Recommended for equilibration only

**Parrinello-Rahman:**
- Extended Lagrangian dynamics
- Correct NPT ensemble fluctuations
- Slower equilibration
- Recommended for production runs

### Soft-Core Potentials

Soft-core potentials prevent numerical instabilities during initial minimization:

**Lennard-Jones Softening:**
```
V_soft(r) = V_LJ(√(r² + α²))
```

**Coulomb Softening (Slater):**
```
V_soft(r) = q₁q₂/r * [1 - (1 + αr/2) * exp(-αr)]
```

**Environment Variables:**
```bash
export UPSIDE_SOFTEN_LJ=1
export UPSIDE_LJ_ALPHA=0.2
export UPSIDE_SOFTEN_COULOMB=1
export UPSIDE_SLATER_ALPHA=2.0
```

## Validation

To validate the implementation:

1. **Check barostat type switching:**
```bash
# Verify barostat type in checkpoint files
python3 -c "
import h5py
with h5py.File('outputs/martini_8stage/checkpoints/bilayer.stage_6.6.up', 'r') as f:
    print('Stage 6.6 barostat type:', f['/input/barostat'].attrs['type'])
with h5py.File('outputs/martini_8stage/checkpoints/bilayer.stage_7.0.up', 'r') as f:
    print('Stage 7.0 barostat type:', f['/input/barostat'].attrs['type'])
"
```

2. **Monitor pressure convergence:**
```bash
# Check log files for pressure values
grep "NPT" outputs/martini_8stage/logs/stage_7.0.log
```

3. **Verify box dimensions:**
```bash
# Extract box dimensions from trajectory
python3 -c "
import h5py
import numpy as np
with h5py.File('outputs/martini_8stage/checkpoints/bilayer.stage_7.0.up', 'r') as f:
    if '/output/volume' in f:
        volumes = f['/output/volume'][:]
        print(f'Volume: {np.mean(volumes):.2f} ± {np.std(volumes):.2f} Å³')
"
```

## Future Enhancements

Potential improvements for future versions:

1. **Lipid headgroup restraints:**
   - Automatic PO4/GL1 bead selection
   - Gradual restraint reduction (200 → 0 kJ/mol/nm²)

2. **Multi-stage HDF5 configuration:**
   - Store all 8 stages in single HDF5 file
   - Automatic stage switching

3. **Enhanced analysis:**
   - Automatic pressure/volume analysis
   - Area per lipid calculation
   - Membrane thickness monitoring

4. **Visualization:**
   - VTF trajectory generation for each stage
   - Automated VMD visualization scripts

## References

1. CHARMM-GUI MARTINI Maker: http://www.charmm-gui.org/
2. MARTINI Force Field: http://cgmartini.nl/
3. Parrinello-Rahman barostat: Parrinello & Rahman, J. Appl. Phys. 52, 7182 (1981)
4. GROMACS MARTINI tutorial: http://www.mdtutorials.com/gmx/martini/

## Support

For issues or questions:
- Check log files in `outputs/martini_8stage/logs/`
- Verify environment variables are set correctly
- Ensure UPSIDE is compiled with latest changes
- Review checkpoint files with h5py

## Changelog

### 2026-02-03
- Initial implementation of Parrinello-Rahman barostat
- Added position restraint potential node
- Created 8-stage bilayer equilibration workflow
- Updated preparation script with barostat type support
