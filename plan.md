# Project Goal

Debug the recent regression that makes both `example/16.MARTINI/outputs/martini_1rkl_hybrid_full/1rkl.stage_7.0.vtf` and `example/16.MARTINI/outputs/martini_1rkl_hybrid/1rkl.stage_7.0.vtf` show unstable protein secondary structure compared with the last committed version.

# Architecture & Key Decisions

- Treat this as a regression against `HEAD`; use read-only git commands only.
- Do not disable or scale hybrid interface interactions to hide the issue.
- Compare behavior from HDF5/log metrics, not only VTF visualization.
- Focus first on changes shared by full-lipid and CGL workflows, because both outputs regress.
- Keep fixes minimal and local to the regression cause.

# Execution Phases

- [x] Phase 1: Inventory changed files and identify workflow/runtime changes common to both full and CGL 1RKL outputs.
- [x] Phase 2: Quantify protein secondary-structure/protein-stability metrics from current full and CGL stage-7 outputs.
- [x] Phase 3: Run controlled copied-checkpoint tests to isolate the temperature regression.
- [x] Phase 4: Implement the minimal fix.
- [x] Phase 5: Update findings/progress with the root cause, verification, and any remaining risk.
- [x] Phase 6: Fix hybrid progress logging to print actual MD step counts instead of rounded elapsed time.

# Known Errors / Blockers

- Existing `martini_1rkl_hybrid` and `martini_1rkl_hybrid_full` stage-7 outputs were generated before the fix with `--temperature 1.2`; they should be regenerated to replace the bad VTF/HDF5 artifacts.

# Review

- Root cause: `example/16.MARTINI/run_sim_hybrid.sh` accidentally changed the default `TEMPERATURE` from the committed `0.8647` to `1.2`. Both reported bad outputs were generated with `--temperature 1.2`.
- Current bad output metrics confirm secondary-structure loss: CGL hbond sum `23.92 -> 14.28`, CA Rg `12.73 -> 14.02 A`; full-lipid hbond sum `26.38 -> 15.45`, CA Rg `12.69 -> 13.73 A`.
- Controlled copied-checkpoint test from the same CGL stage-7 start isolated temperature as the cause. At `T=0.8647`, hbond final `33.48`, hbond last-20 mean `32.15`, CA RMSD max `1.03 A`; at `T=1.2`, hbond final `22.30`, hbond last-20 mean `21.55`, CA RMSD max `1.34 A`.
- Fix applied: restored `TEMPERATURE="${TEMPERATURE:-0.8647}"` in `run_sim_hybrid.sh`.
- Progress logging fix applied: hybrid verbose logs now print actual MD steps as `step nr / n_round time ...` instead of rounded elapsed time in the leading field. A copied 120-step smoke run printed `step 0 / 120`, `step 50 / 120`, and `step 100 / 120`.
