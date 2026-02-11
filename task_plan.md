# Task: Debug 1UBQ MARTINI workflow failure

## 1. Objective
Run `example/16.MARTINI/run_sim_1ubq.sh` and fix the NPT equilibration instability (NaN potential after ~20 steps).

## 2. Architecture & Key Decisions
* **Scope**: Fix workflow stability for 1UBQ without changing core simulation architecture.
* **Source of truth**: Use existing MARTINI protein topology/mapping files already in the repo.
* **Minimal change**: Prefer physically justified parameter updates over ad hoc tweaks.
* **Gromacs alignment**: Mirror CHARMM-GUI Gromacs minimization controls and NPT barostat settings (tau_p, compressibility, isotropic) using Upside unit conversions.

## 3. Execution Phases
- [x] Phase 1: Inspect 1UBQ topology/mapping inputs and parser expectations.
- [x] Phase 2: Implement minimal mapping/topology fix for missing atoms.
- [ ] Phase 3: Re-run the workflow to confirm NPT equilibration is stable.

## 4. Known Errors / Blockers
- NPT equilibration produces NaN potential around step 20; barostat scales box upward with NaN pressure.
- Current run with Gromacs-aligned minimization/barostat is in progress; no NaNs observed through 100 steps of NPT equilibration.
