# Task: Enforce Weighted Hybrid SC Coupling Without Pair-Kernel Fallback

## 1. Project Goal
Ensure hybrid SC probabilistic coupling always uses weighted SC-row interactions (including SC-BB and SC-env cases), while keeping interaction evaluation strictly on dry-MARTINI-side particles (no AA-carrier interaction fallback in the pair kernel).

## 2. Architecture & Key Decisions
- Scope:
  - Inspect and patch hybrid SC coupling logic in `src/martini.cpp` only.
  - Preserve existing hybrid schema and runtime node dependencies.
  - Avoid changes to stage-file generation unless required.
- Key decisions:
  - Keep per-row probabilistic model; make weight normalization robust by enforcing fallback uniform weights when probability sums are degenerate.
  - Keep probabilistic edge path active whenever SC proxy rows have valid weighted rows.
  - Remove deterministic MARTINI pair fallback for SC-configured proxy rows when no active row weights are available in that step.
  - Disable protein SC-SC nonbonded interactions (deterministic and probabilistic paths) to prevent instability.
  - Enforce strict rotamer-position locking for probabilistic SC rows during interaction evaluation (no displacement by env/SC interaction-driven relaxation) for SC-env paths.
- Revised Decisions:
  - Use normalized row weights (`sc_row_prob_norm`) consistently for edge eligibility checks.

## 3. Execution Phases
- [x] Phase 1: Keep strict placement-based SC-row mapping (no pair-kernel fallback to carrier/proxy coordinates).
- [x] Phase 2: Patch SC-row probability normalization to guarantee nonzero per-proxy row weights when valid rows exist.
- [x] Phase 3: Patch SC-edge selection to rely on normalized weights and preserve probabilistic routing.
- [x] Phase 4: Build and run focused validation.
- [x] Phase 5: Record outcomes in `progress.md` and task review.

## 4. Known Errors / Blockers
- None currently.

## 5. Review
- Implemented in `src/martini.cpp`:
  - Hardened per-proxy SC-row probability normalization so that if valid rows exist but summed probability is degenerate/non-positive, weights fall back to a uniform distribution over valid rows.
  - Updated SC probabilistic edge eligibility to use normalized per-row weights (`sc_row_prob_norm`) and valid-row masks.
  - Removed deterministic SC pair fallback inside probabilistic SC selection: when a proxy is SC-configured but has no active weighted rows for the current step, that pair is skipped instead of evaluated through unweighted deterministic MARTINI fallback.
  - Disabled SC-SC interactions by restoring role rule to `ROLE_SC/ROLE_SC -> false`.
  - Added an explicit pair-loop guard to skip protein SC-SC nonbonded interactions before deterministic/probabilistic SC routing.
  - Locked probabilistic SC row positions to their mapped rotamer targets during SC-env evaluation, disabling interaction-driven row displacement.
- Validation:
  - `source .venv/bin/activate && source source.sh && cmake --build obj -j4` succeeded.
  - Build emitted only existing warnings; no new errors.
