## Overview
Upside is a molecular dynamics simulation package for protein folding and conformational dynamics. It combines a fast C++ core with Python scripts for configuration and analysis.

### File Format Reference

Before editing any `.up` simulation input file, force-field `.h5` parameter file, or preparation script in `py/`, read **`up.md`** in the repo root. It documents:
- The HDF5 structure of `.up` files (all `/input` datasets, all `/input/potential` node types and their dataset schemas, the `/output` trajectory layout)
- Every force-field parameter file (`sidechain.h5`, `hbond.h5`, `environment.h5`, `bb_env.dat`, `martini.h5`, `dopc.h5`, `rama.dat`, `sheet`, `membrane.h5`)
- Hybrid-specific groups (`hybrid_bb_map`, `hybrid_remap`, `cgl_gle`, `barostat`, `stage_parameters`)
- Critical conventions: angle storage as cosines, CGL target angle sign, GLY Ramachandran handedness (it lives in the library file, never in the config writer), spline table exactness, thermostat timescales

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
* **No Hard-Coded Protein / System Identity**: Scripts under `py/` are shared infrastructure and must never name a specific protein, PDB id, or system as a default or fallback. A hard-coded id does not fail loudly: it silently attaches the wrong protein's metadata to another system's trajectory, or writes another system's outputs into a foreign run directory. Take the identity from an explicit argument or environment variable; derive every dependent path from it (`run_dir`, `runtime_pdb_id`, `protein_aa_pdb`, metadata PDB). If it is genuinely absent, either raise a clear error naming the missing option, or skip the optional feature that needed it; never substitute a guess. Sample-specific defaults belong in the per-example shell scripts under `example/`, not in `py/`. The instances found and fixed are in `findings.md` 1.7.

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

### Shared Upside Deployment (midway2 + midway3)

**`/beagle3/trsosnic/yinhan/upside2-md` is the one deployment both clusters use.**

```bash
source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh
```

That single file gives an identical binary, force fields and Python on midway2 and midway3.
Do not build a second per-cluster tree; the whole point is that there is one. The ConDiv training
tree, `/project/trsosnic/yinhan/upside2-md-mdw2` (`$P` in `remote_jobs.md`), is separate: it runs
only on midway2, uses its own `training/env.sh`, and is never overwritten from the repo during a
campaign (`remote_jobs.md` §1b).

Why it works, each point measured rather than assumed:

* **`/beagle3` is visible and writable from the COMPUTE NODES of both clusters**, as are `/project`
  and `/project2`. **`/cds3` is login-node only** and cannot be used by jobs.
* **One binary serves both.** It is compiled on midway2 (Broadwell) where `-march=native`
  (`src/CMakeLists_Other.txt:8`) emits no AVX-512, so it also runs on midway3's Cascade Lake.
  **Always compile on midway2.** Building on midway3 produces AVX-512 and the binary dies with
  SIGILL on broadwl.
* **`hdf5/1.14.3+oneapi-2023.1` exists on both** and supplies the `libhdf5.so.310` the binary needs.
  midway2 *additionally* needs `gcc/10.1.0`, because its system libstdc++ lacks `GLIBCXX_3.4.20`.
* **The venv is portable because its interpreter lives on shared storage.** `pyrt/` is a copy of
  midway2's el7 python 3.9.18, and `.venv` is built from it, so nothing points at `/software`,
  which is **per-cluster** and is why earlier venvs worked on only one machine. An el7 interpreter
  runs on el8 (midway3) but not the reverse, the same compatibility direction as the binary.
  `env_shared.sh` puts `pyrt/lib` on `LD_LIBRARY_PATH` for `libpython3.9.so.1.0`.
  Carries numpy, pytables, prody and h5py; verified importing on both clusters.

**Filesystem headroom, and a trap.** Check the *fileset* or *group* quota, never `df` on a mount
point, which reports the whole device (`/project` shows 6.3 PB). From midway2, `df` on the
subdirectory (`df -h /project/trsosnic`) reports the fileset, and `rcchelp quota` reports only home,
scratch and the `/project2` group. Current headroom is tracked in `remote_jobs.md` (disk section);
`/project2` sits near its group quota and `/cds3` is unusable from compute nodes.

### Slurm Environment Setup

Cluster jobs do not use the Mac `source.sh` bootstrap. Each sources one environment file that loads
the modules, activates its tree's venv and sets `UPSIDE_HOME`, `PATH` and `PYTHONPATH`:

* **Shared deployment** (benchmarks, glpG, anything that may run on either cluster):
  `source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh`.
* **ConDiv training** (`$P`, midway2 only): the run directory's own copy of `$P/training/env.sh`,
  which `train_chain.sbatch` sources. The repo's `training/env.sh` is site-neutral (it only
  activates the repo `.venv`); the cluster copies carry the midway2 module loads.
* **HDX analysis** is the exception: it runs on midway3 from `.venv_el8_py311_bak` (pymbar,
  matplotlib), set by hand rather than through that venv's stale `activate` (memory
  `renamed-venv-activate-trap`).

Rules:

* For interactive local Mac work: `source .venv/bin/activate && source source.sh`.
* Do not add module loads on top of these files. Both cluster venvs are built from python 3.9.18
  (`env_shared.sh` from `pyrt/`, the training venv from midway2's `python/3.9.18` module), and each
  file's module loads and library paths are what let `libupside.so` and `libpython3.9.so.1.0` load.
* A wrapper that sets up the environment itself sets `UPSIDE_SKIP_SOURCE_SH=1` before invoking
  lower-level workflow scripts, so they do not re-enter the local-only bootstrap.
* Never pipe `source <env file>`: the pipe runs it in a subshell and the environment is lost.
* Submit through midway2 (Default Cluster below), and verify a new environment on a compute node
  (`command -v python3`, `python3 -V`, the imports the job needs) before submitting a fleet.

Example wrapper skeleton:

```bash
#!/bin/bash
#SBATCH --account=pi-trsosnic
#SBATCH --partition=broadwl
#SBATCH --time=36:00:00
set -euo pipefail

source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh
export UPSIDE_SKIP_SOURCE_SH=1

bash "$UPSIDE_HOME/example/16.MARTINI/<workflow>.sh"
```

### Remote Job Records

**"continue jobs":** a session asked to continue the jobs starts at `remote_jobs.md` "Resume here,
from any computer" and follows its checklist in order.

**All remote job recording goes into `remote_jobs.md`, and nowhere else.** It is the single source of
truth for what is running on midway2 and midway3: job ids, what each job is, its submit script, its
log path, its data directory, and the next action each one is waiting on. Do not record jobs in `plan.md`,
`progress.md` or `findings.md`: those track technical direction, execution history and knowledge
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

### Compiling Upside on the cluster (midway2 login node only)

**Always compile on midway2.** `-march=native` (`src/CMakeLists_Other.txt:8`) on midway3's Cascade
Lake emits AVX-512 and that binary dies with SIGILL on broadwl; a midway2 (Broadwell) build runs on
both clusters. Build into a fresh directory and install by rename, never with `install.sh` in a live
tree: `install.sh` empties `obj/` first, under jobs that have the binary loaded.

```bash
source /software/modules/init/bash
module load gcc/10.1.0 hdf5/1.14.3+oneapi-2023.1
TREE=/project/trsosnic/yinhan/upside2-md-mdw2      # or /beagle3/trsosnic/yinhan/upside2-md
mkdir -p $TREE/obj_<tag> && cd $TREE/obj_<tag>
/software/cmake-3.19-el7-x86_64/bin/cmake ../src -DEIGEN3_INCLUDE_DIR=/software/eigen-3.3-el7-x86_64/include/eigen3
make -j8
```

Then back up `obj/upside` and `obj/libupside.so` as `.bak_<tag>`, copy the new files in with `cp -p`
and `mv` (a deploy must carry file modes; `findings.md` 10.11), and check bitwise parity against the
installed build before any job uses it. The 2026-10-02 deploy is the worked example
(`/project/trsosnic/yinhan/checks/hbg_deploy_20261002/deploy_hbg.sh`, `finish_hbg.sh`). The build is
CPU-only and takes a few minutes on the login node; analysis heavier than a quick read goes to a
compute node.

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

### Default Cluster for Remote Jobs

**Unless the user explicitly specifies otherwise, all remote jobs go to midway2 (broadwl partition), not midway3.**

Midway3 (caslake) has very long Priority queue times and should only be used when midway2 is unavailable or when a job requires midway3-specific resources. Submitting from the midway3 login node routes jobs to caslake even if the script requests broadwl, so always submit training and simulation jobs through the **midway2 SSH socket** (`~/.ssh/cm-mdw2.sock`). Never send CPU work to a GPU partition (`amd`, `beagle3`) because it accepts the job: that bills the group's GPU allocation. midway3's login node may be used to move or read files on the shared `/project` and `/beagle3`.

```bash
# Submit to midway2:
ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu 'cd <project_dir> && sbatch $(cat /project/trsosnic/yinhan/slurm/midway2.args) <script>'
```

**midway2 node exclusions are one law** (`remote_jobs.md` §0d). `/project/trsosnic/yinhan/slurm/midway2.args`
holds `--partition=broadwl` and the excluded nodes; it is the only place a midway2 node list is written.
* Submit with `sbatch $(cat /project/trsosnic/yinhan/slurm/midway2.args) ...`, or through a
  `slurm.args` made with `ln -s /project/trsosnic/yinhan/slurm/midway2.args slurm.args`. Never write
  `--exclude=midway2-...` into a script, and never copy the list into a `slurm.args`.
* `~/bin/sbatch` on midway2 adds the law's `--exclude` to every other submission. A caller's own
  `--exclude` replaces it, so pass one only on purpose.
* Add a node only for its own record of failures. Edit the one line, record the evidence in
  `README.md` beside it, then run `update_pending.sh` there so queued jobs follow.

### Cluster SSH: Self-Connection

**Do not ask the user to connect SSH.** Open the ControlMaster socket yourself; the user only approves the Duo push on their phone.

| cluster | open the socket | socket |
| --- | --- | --- |
| midway2 (all job submission) | `expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw2_master.exp` | `~/.ssh/cm-mdw2.sock` |
| midway3 (file moves and reads) | `expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw3_master.exp` | `~/.ssh/cm-mdw3.sock` |

Both scripts read the password from `~/.bin/ssh_mdw3`, select Duo option 1 (push), and keep the socket with `ControlPersist=8h`. If midway2 refuses this IP, `scratchpad/mdw2_via_mdw3.exp` opens the midway2 socket through a live midway3 socket (`remote_jobs.md` §0).

Check the socket before every call and always pass `-o BatchMode=yes`:

```bash
if ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu >/dev/null 2>&1; then
    ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu '<command>'
fi
```

On a dead master a plain `ssh -S` falls back to password attempts, and a few of those get this IP throttled. If the check fails, run the expect script once; every launch sends a Duo push, so never retry it in a loop, and after a throttle wait at least 30 minutes. The full connection handbook is `remote_jobs.md` §0.
