# Task: Fix Hybrid Mapping Atom-Count Mismatch During Stage Preparation

## 1. Project Goal
Resolve `Hybrid mapping n_atom mismatch` in `example/16.MARTINI/run_sim_1rkl.sh` by ensuring stage conversion reads the same runtime packed PDB/ITP that mapping export used.

## 2. Architecture & Key Decisions
- Scope:
  - Add runtime input-path overrides in `example/16.MARTINI/prepare_martini.py` via environment variables.
  - Wire `example/16.MARTINI/prepare_system.py` to set those env vars when invoking `prepare_martini.py`.
- Revised Decisions:
  - Keep existing positional CLI behavior of `prepare_martini.py`; only add optional env-based path overrides.
  - Preserve fallback behavior to legacy `pdb/<pdb_id>.MARTINI.pdb` and `pdb/<pdb_id>_proa.itp` when overrides are absent.

## 3. Execution Phases
- [x] Phase 1: Confirm mismatch root cause across `run_sim_1rkl.sh`, `prepare_system.py`, and `prepare_martini.py`.
- [x] Phase 2: Implement runtime path override support in `prepare_martini.py`.
- [x] Phase 3: Pass runtime override env vars from `prepare_system.py` subprocess call.
- [x] Phase 4: Validate syntax and run focused reproduction/smoke check.
- [x] Phase 5: Log outcomes in `progress.md`.

## 4. Known Errors / Blockers
- None currently.
