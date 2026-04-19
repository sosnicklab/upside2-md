# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview
Upside is a molecular dynamics simulation package for protein folding and conformational dynamics. It combines a fast C++ core with Python scripts for configuration and analysis.

## Hybrid Interface RMSF Calibration

- For the hybrid `1rkl` membrane workflow, the current chosen production interface factor is:
  - `PROTEIN_ENV_INTERFACE_SCALE = 1.15`
- This value is used in:
  - `example/16.MARTINI/run_sim_1rkl.sh`
- Calibration target:
  - reproduce the protein fluctuation amplitude of the implicit membrane reference from `example/08.MembraneSimulation`
  - compare only the membrane-embedded protein region, not the solution-exposed residues
- Decision rule:
  - prefer direct reference-vs-hybrid RMSF overlays over weak fitted trend lines
  - current selected scale from the latest stable rerun is `1.15`

## Development Rules
## Development Rules
* **Backward Compatibility**: Modifications to C++ source files must not break existing function calls or the Python-to-C++ interface. Exception: This rule does not apply to `src/box.h`, `src/box.cpp`, `src/martini.cpp`, and any files under the `example/16.MARTINI/` directory.
* **Function Signatures**: When adding parameters to an existing function, the additional parameters must be optional (i.e., provide default values), except when modifying `src/box.h`, `src/box.cpp`, `src/martini.cpp`, or files within `example/16.MARTINI/`.
* **Master Branch Parity**: The `master` branch is the gold standard; all modifications must produce results identical to those of the `master` branch for existing simulation configurations.
* **Memory Layout**: Do not reorder existing member variables in classes accessed by Python to avoid memory corruption.
* **Deprecation**: Mark old functions as deprecated instead of removing them to support legacy scripts. Exception: Old functions may be removed directly in `src/box.h`, `src/box.cpp`, `src/martini.cpp`, and files under `example/16.MARTINI/`.

### Environment Setup
Crucial: You must run these commands from the project root before running anything in this project:
```bash
source .venv/bin/activate
source source.sh

```

### Slurm Environment Setup
For Slurm jobs on the cluster, do not rely on the Apple-Silicon `source.sh` bootstrap as the primary environment setup.

Use a self-contained Slurm setup from the repo root:
```bash
PROJECT_ROOT=/path/to/upside2-md

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load python/3.11.9 || true
  module load cmake || true
  module load openmpi || true
  module load hdf5/1.14.3 || true
fi

if [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
  source "$PROJECT_ROOT/.venv/bin/activate"
fi

export UPSIDE_HOME="$PROJECT_ROOT"
export PATH="$PROJECT_ROOT/obj:$PATH"
export PYTHONPATH="$PROJECT_ROOT/py${PYTHONPATH:+:$PYTHONPATH}"
```

Rules:
- For interactive local Mac work: `source .venv/bin/activate && source source.sh`.
- For Slurm jobs: prefer module load + repo `.venv` activation + explicit `UPSIDE_HOME/PATH/PYTHONPATH`.
- If a Slurm wrapper sets up the environment itself, it should set `UPSIDE_SKIP_SOURCE_SH=1` before invoking lower-level workflow scripts so they do not re-enter the local-only bootstrap path.
- A proper Slurm job for this project should do these steps in order:
  1. Resolve `PROJECT_ROOT` explicitly.
  2. Source `/etc/profile.d/modules.sh` when available.
  3. Load the required modules: Python, CMake, OpenMPI, and HDF5.
  4. Activate `PROJECT_ROOT/.venv` if it exists.
  5. Export `UPSIDE_HOME="$PROJECT_ROOT"`.
  6. Prepend `PROJECT_ROOT/obj` to `PATH`.
  7. Prepend `PROJECT_ROOT/py` to `PYTHONPATH`.
  8. Set `UPSIDE_SKIP_SOURCE_SH=1` if the wrapper is handing off to lower-level workflow scripts that would otherwise source the local Mac bootstrap.

Example Slurm wrapper skeleton:
```bash
#!/bin/bash
#SBATCH --time=36:00:00
set -euo pipefail

PROJECT_ROOT=/path/to/upside2-md

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load python/3.11.9 || true
  module load cmake || true
  module load openmpi || true
  module load hdf5/1.14.3 || true
fi

if [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
  source "$PROJECT_ROOT/.venv/bin/activate"
fi

export UPSIDE_HOME="$PROJECT_ROOT"
export PATH="$PROJECT_ROOT/obj:$PATH"
export PYTHONPATH="$PROJECT_ROOT/py${PYTHONPATH:+:$PYTHONPATH}"
export UPSIDE_SKIP_SOURCE_SH=1

bash "$PROJECT_ROOT/example/16.MARTINI/run_sim_1rkl_outlipid.sh"
```

### Installation & Build
```bash
# Install dependencies and compile C++ core
./install_M1.sh
./install_python_env.sh
```

### Upside Unit Conversions

| Quantity | Upside Unit | Standard Equivalent |
| :--- | :--- | :--- |
| **Energy** | 1 E_up | 2.914952774272 kJ/mol |
| **Length** | 1 Angstrom | 1 Angstrom |
| **Mass** | 1 m_up | 12 g/mol |
| **Temperature** | 1.0 T_up | 350.588235 Kelvin |
| **Pressure** | 0.000020933215 E_up / (Angstrom^3) | 1 atm |
| **Pressure** | 0.000020659477 E_up /(Angstrom^3) | 1 bar |
| **Compressibility** | 14.521180763676 Angstrom^3 / E_up | 3e-4 bar^(-1) |

### Dry-MARTINI Unit Contract

Training artifacts under `SC-training/` stay in native dry-MARTINI units.

- implementation rule:
  - training outputs and forcefield parameters are authored in native dry-MARTINI units (`nm`, `kJ/mol`, `e`);
  - simulation code must not bake dry-MARTINI to Upside conversion numbers into the training artifacts;
  - simulation code must receive the required unit-conversion factors as explicit parameters / attrs when converting native dry-MARTINI quantities into Upside runtime units.
