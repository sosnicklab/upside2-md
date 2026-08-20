## Overview
Upside is a molecular dynamics simulation package for protein folding and conformational dynamics. It combines a fast C++ core with Python scripts for configuration and analysis.

## Physical Model Integrity
**CRITICAL: Do not modify, scale to zero, or disable core physics interactions.**
* **Hybrid Interface Interactions**: The interaction potentials between the protein Side Chains and the dry-MARTINI environment (SC-env), as well as the protein Backbone and the dry-MARTINI environment (BB-env), must **NEVER** be turned off.
* **No "Debugging" Exclusions**: Do not disable or bypass these hybrid interface interactions to circumvent crashes, optimize performance, or troubleshoot workflow scripts. Disabling them completely breaks the physical model of the Upside simulation.
* **Strict Adherence**: Any generated script, parameter modification, or configuration must 100% respect this physical model.

## Purpose: Make The Result RIGHT, Not Make It "Run"

**The objective is always a correct scientific result, never a job that merely completes.** A run that
finishes, exits 0, or stops producing errors is not evidence of anything. Judge every change by whether
the physics is right, and report a job as working only after checking the observables that would reveal
it is not.

* **A green exit code means nothing.** Slurm has reported `COMPLETED` on runs whose every replica was
  destroyed. `isfinite` has passed on coordinates of 1e12 A. Check the physical observables — Rg,
  peptide C-N bond lengths, the sign and trend of the potential, `avg_kinetic_energy/1.5kT` — not the
  return code.
* **Never make something "work" by weakening the test.** Do not widen a threshold, disable a check,
  skip a frame, or relax a tolerance so a run proceeds. If a check fires, first establish whether it is
  right. Recalibrate a detection threshold only from measured data on the system it guards, and say so.
* **Never make something "work" by adjusting physics.** Timesteps, friction, masses, ion counts, box
  sizes and force-field parameters are determined by the model, not chosen to avoid a crash. Several
  are calibrated against each other (glpG's dt is locked to `/input/brownian` because the friction is
  tuned against it for a target lipid diffusion) — changing one silently invalidates the other.
* **Do not transfer settings, thresholds or analysis between different simulations.** They are separate
  systems with separate methods and integrators. A number derived from one is not evidence for another.
* **Diagnose before concluding.** Measure the mechanism; do not assert one by analogy. State plainly
  when a cause is unproven, and correct the record when a claim turns out to be wrong.

## Development Rules
* **Backward Compatibility**: Modifications to C++ source files must not break existing function calls or the Python-to-C++ interface.
* **Function Signatures**: When adding parameters to an existing function, the additional parameters must be optional (i.e., provide default values).
* **Master Branch Parity**: The `master` branch is the gold standard; all modifications must produce results identical to those of the `master` branch for existing simulation configurations.
* **Memory Layout**: Do not reorder existing member variables in classes accessed by Python to avoid memory corruption.
* **Deprecation**: Mark old functions as deprecated instead of removing them to support legacy scripts.
* **Physical Interactions**: Except for the Upside core, which uses a trained force field for protein dynamics, everything else should be physical. "Twisting parameters to make it work" is not allowed. No arbitrary capping, and no additional orientational potential for CGL. The CGL force field should not contain information of the bilayer. A stable bilayer should be a result of correct force field evolving, not something you twist parameters into.
* **Spline Table Only**: During simulations, all interactions computed by Upside should use a spline table. Even if two particles have a simple Lennard-Jones potential, their interaction potentials need to be written to an .h5 file before the simulation, and Upside needs to read them from the .h5 file.
* **Spline Tables Must Be The Original Potential**: A spline table is a representation, never a variant. Evaluated in native force-field units it must equal the published functional form exactly, including how that form reaches its cutoff — for dry-MARTINI that means reaction-field electrostatics (`epsilon_r = 15`, `epsilon_rf = 0`) and a potential-shifted Lennard-Jones, both going to zero at 1.2 nm. Do not tabulate a bare truncation, a re-fit, or a smoothed approximation. Assert the equivalence against the analytic form after unit conversion rather than assuming it.
* **NO GUARDS**: Never add a guard that hides, masks, bypasses, or works around a numerical problem. This is prohibited without exception, and includes: skipping a pair, term, or gradient because its value is non-finite; clamping, capping, or flooring a force, energy, or displacement to keep a run alive; catching a blow-up and continuing; and aborting on NaN in place of explaining it. A guard destroys the evidence needed to find the cause and silently corrupts the physics that survives it. When a run produces NaN, diverges, or blows up, that is the signal to find the real defect — instrument the code, localize the event, and fix the underlying force field, table, exclusion, or integration error. A validated precondition on a genuinely invalid domain (a box length must be positive; a required input must exist) is not a guard and is fine. If you believe an exception is warranted, stop and ask; do not add one.
* **H5 Force Field Files**: Do not make version numbers for h5 force field files. Backup the old ones and overwrite them.
* **No Hard-Coded Protein / System Identity**: Scripts under `py/` are shared infrastructure and must never name a specific protein, PDB id, or system as a default or fallback. A hard-coded id does not fail loudly — it silently attaches the wrong protein's metadata to another system's trajectory, or writes another system's outputs into a foreign run directory. Take the identity from an explicit argument or environment variable; derive every dependent path from it (`run_dir`, `runtime_pdb_id`, `protein_aa_pdb`, metadata PDB). If it is genuinely absent, either raise a clear error naming the missing option, or skip the optional feature that needed it — never substitute a guess. Sample-specific defaults belong in the per-example shell scripts under `example/`, not in `py/`. Fixed 2026-08-09: `martini_extract_vtf.infer_pdb_id` returned `"1rkl"` for any non-bilayer input (this mislabelled the 1AO6 nanoparticle VTFs), and `martini_prepare_system` defaulted `--pdb-id` to `1rkl` and `--run-dir` to `outputs/martini_test_1rkl_hybrid`.

### dryMARTINI Interface Refactoring Rules
* **Master Repository Path**: Use `/Users/yinhan/Documents/upside2-md-master` as the master repository reference for all file diffs and code comparisons.
* **Architectural Integrity**: Clean up and thoroughly rewrite the dryMARTINI interface code (including Python, C++, and MD scripts) to ensure a straight, logical, and cohesive architecture.
* **Code Quality**: Eliminate the current fragmented, patch-on-patch structure introduced by previous AI iterations.
* **Inactive Flags**: Completely remove any inactive, disabled, or unused configuration flags within the diffed scope. Do not leave them as commented-out code or dead toggles. 
* **Stylistic Matching**: Exactly match the formatting, naming conventions, and style of the human-written code found in the master repository.
* **Exclusions**: Completely ignore `/Users/yinhan/Documents/upside2-md-master/example/00.AnalysisScripts` for style reference, as it is AI-written and not a valid baseline.

### The "Clean Slate" Exception
The rules for **Backward Compatibility**, **Function Signatures**, and **Deprecation** DO NOT apply to the code actively being developed. You must determine this scope by diffing the current repository against the master repository. 

**CRITICAL RULE FOR DIFFED FILES:** Keep the actively modified or newly added interface files impeccably clean. Do not build layers of disabled code, do not leave commented-out legacy blocks, and do not write wrapper functions for old implementations. You must completely remove old or unused code and write the new implementations directly.

### Coherence of Edits
**All edits must be coherent.** Every change to a C++, Python, or Markdown file must leave that file reading as a single, unified, logical whole — not layers of patches on patches, disabled toggles, commented-out legacy blocks, or wrappers around old implementations. When you fix or extend something, overwrite or delete the old implementation and write the new one directly, so the file reads as though authored in one pass. This applies equally to source (`.cpp`/`.h`), Python (`.py`), and documentation (`.md`).

Two exceptions only:
1. **Master parity takes precedence.** Do not refactor for coherence when doing so would change code that must remain identical to the `master` branch. Where a coherent rewrite would break master parity, keep master parity and leave that code alone.
2. **Temporary debugging scaffolding.** Incoherent, throwaway edits (scaffolding, probes, temporary toggles) are allowed *while actively debugging*, but must be removed and the code returned to a coherent state before the task is considered complete.

### Environment Setup
Crucial: You must run these commands from the project root before running anything in this project:
```bash
source .venv/bin/activate
source source.sh

```

### Slurm Environment Setup

For Slurm jobs on the cluster, do not rely on the Apple Silicon `source.sh` bootstrap as the primary environment setup.

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

* For interactive local Mac work: `source .venv/bin/activate && source source.sh`.
* For Slurm jobs: prefer module load + repo `.venv` activation + explicit `UPSIDE_HOME/PATH/PYTHONPATH`.
* If a Slurm wrapper sets up the environment itself, it should set `UPSIDE_SKIP_SOURCE_SH=1` before invoking lower-level workflow scripts so they do not re-enter the local-only bootstrap path.
* A proper Slurm job for this project should complete these steps in order:

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

### Remote Job Records

**All remote job recording goes into `remote_jobs.md`, and nowhere else.** It is the single source of
truth for what is running on midway3: job ids, what each job is, its submit script, its log path, its
data directory, and the next action each one is waiting on. Do not record jobs in `plan.md`,
`progress.md` or `findings.md` — those track technical direction, execution history and knowledge
respectively, and a job table duplicated across them goes stale silently and then misleads.

Maintenance rules:
* Update the snapshot date and the job table whenever jobs are submitted, finish, or are cancelled.
* **Delete stale rows.** A finished or cancelled job is removed from the current-jobs table; it is kept
  only as a one-line entry under history if it carries a lesson worth reusing. A table listing jobs that
  no longer exist is worse than no table.
* Record the log path and data directory for every job, so a cold session can check it without guessing.
* Record what must not be forgotten alongside the jobs (STOP files, binaries carrying temporary
  instrumentation, undecided policy questions), since those are what break a resumed session.

### Installation & Build

```bash
# Install dependencies and compile C++ core
./install_M1.sh
./install_python_env.sh

```

### Compiling Upside on midway3 (login node)

To (re)compile the C++ core on the RCC midway3 login node, from the repo root:

```bash
# 1. Load build modules (names verified on midway3)
module load cmake
module load openmpi
module load hdf5/1.14.3

# 2. Provide EIGEN_HOME (and PATH/PYTHONPATH) by sourcing the repo's source_sh
source source_sh

# 3. Build: install.sh copies source_sh -> source.sh, then runs cmake + make
./install.sh
```

Notes:
* `install.sh` runs `cmake ../src/ -DEIGEN3_INCLUDE_DIR=$EIGEN_HOME && make` in `obj/`, so
  `EIGEN_HOME` (set by `source_sh`, e.g. `/software/eigen-3.4-el8-x86_64/include/eigen3`) must
  resolve to a valid Eigen 3.4 include directory before building.
* The build is Python-independent; only C++ recompilation requires these steps. `UPSIDE_HOME`
  in `source_sh` is unrelated to compilation (cmake uses the repo's own `src/`).
* This compiles on the login node (no Slurm job needed); it is CPU-only and takes a few minutes.

### Upside Unit Conversions

| Quantity | Upside Unit | Standard Equivalent |
| --- | --- | --- |
| **Energy** | 1 E_up | 2.914952774272 kJ/mol |
| **Length** | 1 Angstrom | 1 Angstrom |
| **Mass** | 1 m_up | 12 g/mol |
| **Temperature** | 1.0 T_up | 350.588235 Kelvin |
| **Pressure** | 0.000020933215 E_up / (Angstrom^3) | 1 atm |
| **Pressure** | 0.000020659477 E_up / (Angstrom^3) | 1 bar |
| **Compressibility** | 14.521180763676 Angstrom^3 / E_up | 3e-4 bar^(-1) |

### Dry-MARTINI Unit Contract

Training artifacts under `SC-training/` stay in native dry-MARTINI units.

* Implementation rules:
* Training outputs and forcefield parameters are authored in native dry-MARTINI units (`nm`, `kJ/mol`, `e`).
* The simulation code must not bake dry-MARTINI to Upside conversion numbers into the training artifacts.
* The native dry-MARTINI to Upside unit conversion happens ONCE, at h5-build time in Python, which receives the required conversion factors as explicit parameters. The runtime h5 force-field files and configs therefore store Upside-unit values (energies in `E_up`, lengths in Angstrom), and the simulation code (C++ engine) performs NO unit conversion — it reads the pre-converted spline tables directly.

### midway3 SSH — Self-Connection

**Do not ask the user to connect SSH.** Establish the ControlMaster socket yourself:

```bash
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw3_master.exp
```

This script reads the password from `~/.bin/ssh_mdw3`, sends it, selects Duo option 1 (push), and opens `~/.ssh/cm-mdw3.sock` with `ControlPersist=8h`. After it exits, reuse with:

```bash
ssh -S ~/.ssh/cm-mdw3.sock yinhanw@midway3.rcc.uchicago.edu '<command>'
```

Check socket liveness before running: `ssh -S ~/.ssh/cm-mdw3.sock -O check yinhanw@midway3.rcc.uchicago.edu`. If the socket is missing or dead, re-run the expect script. The user still approves the Duo push on their phone — but you initiate it; you never ask them to run the script themselves.

