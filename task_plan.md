# Task: Migrate bilayer workflow to dry MARTINI protocol

## 1. Project Goal
Update `example/16.MARTINI/run_sim_bilayer.sh` so the bilayer workflow uses dry MARTINI setup and relaxation stages aligned with CHARMM-GUI, and ensure pressure coupling supports dry MARTINI membrane equilibration semantics.

## 2. Architecture & Key Decisions
- Scope: Keep changes focused on bilayer workflow, preparation wiring, and barostat behavior.
- Relaxation parity: Use explicit stages mirroring CHARMM-GUI sequence (`6.0`, `6.1`, `6.2`-`6.6`, `7.0`).
- Revised Decisions:
  - Implement axis-specific compressibility in NPT settings (`compressibility_xy`, `compressibility_z`) for membrane semi-isotropic coupling.
  - Keep backward compatibility by still honoring legacy single `compressibility` when new attributes are absent.
  - Fix pressure-convergence checks for zero target pressure to avoid divide-by-zero behavior.
  - Implement dry MARTINI lipid-head restraints from `BILAYER_LIPIDHEAD_FC` for stages `6.2`-`6.6` with stage-specific ramp (`200/100/50/20/10`) and disabled in production (`7.0`), matching CHARMM-GUI stage semantics.
  - Use NVT as the default for stage `7.0` production in this dry bilayer workflow (matching CHARMM-GUI `step7_production.mdp` `pcoupl = no` and the paper note for implicit-solvent production), while keeping optional NPT override via environment variables.
  - Keep a softened transition for early equilibration in UPSIDE (`6.2: npt_equil`, `6.3: npt_equil_reduced`) to avoid immediate post-minimization blow-up, then switch to hard interactions for `6.4`-`7.0`.
- Compatibility: Preserve existing interfaces; new attributes/env vars are optional extensions.

## 3. Execution Phases
- [x] Phase 1: Inspect current bilayer script, dry MARTINI inputs, and `.mdp` control parameters.
- [x] Phase 2: Determine gaps in pressure-coupling capabilities (`src/box.cpp` and related interface).
- [x] Phase 3: Implement code/script changes for dry MARTINI workflow and relaxation stages.
- [x] Phase 4: Run targeted validation checks and document outcomes.

## 4. Known Errors / Blockers
- Full runtime/physics validation requires executing the long bilayer workflow; only static checks were run in this change.
