Project Goal
- Hard-code `icb_energy` and `ihb_energy` to the non-membrane baseline in `ConDiv_symlay` so the training constraint no longer learns a hidden inner membrane branch.

Architecture & Key Decisions
- Treat the non-membrane baseline for `icb/ihb` as zero membrane increment. The non-membrane protein physics already comes from the ordinary force-field terms outside `membrane.h5`.
- Freeze `icb/ihb` everywhere that matters:
  - constrained seed generation
  - per-minibatch projection
  - checkpoint/restart state sanitization
  - emitted `membrane.h5` files
  - optimizer gradients and updates
- Keep the outer trained channels unchanged:
  - `cb_energy`
  - `hb_energy`
- Make validation report deviation from the required zero inner baseline instead of silently treating arbitrary inner curves as valid.

Execution Phases
- [x] Update tracking files with the new inner-baseline task and the relevant finding/lesson.
- [x] Patch `ConDiv_symlay` to pin `icb/ihb` to zero baseline across seed, update, resume, and writeout.
- [x] Run focused validation to confirm zero inner curves in constrained outputs and a short local smoke run.

Known Errors / Blockers
- None.

Review
- Patched [ConDiv_mem.py](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/ConDiv_mem.py) so `icb/ihb` are pinned to the non-membrane baseline during init-load, restart sanitization, optimizer updates, and emitted force-field writeout.
- Patched [symlay_utils.py](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/symlay_utils.py) so constrained projection forces `icb/ihb` to zero and validator stats measure deviation from that required baseline.
- Updated [TRAINING_POTENTIAL_DESIGN.md](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/TRAINING_POTENTIAL_DESIGN.md), [findings.md](/Users/yinhan/Documents/upside2-md/findings.md), [progress.md](/Users/yinhan/Documents/upside2-md/progress.md), and [lessons.md](/Users/yinhan/Documents/upside2-md/lessons.md) to reflect the new inner-branch rule.
- Verification:
  - `python3 -m py_compile ConDiv_symlay/ConDiv_mem.py ConDiv_symlay/symlay_utils.py ConDiv_symlay/validate_symlay_constraints.py`
  - `bash -n ConDiv_symlay/run_init.sh`
  - controlled `run_init.sh` smokes in `/tmp/condiv_symlay_inner_baseline_smoke` and `/tmp/condiv_symlay_inner_baseline_smoke_fast`
  - direct HDF5 checks showing zero `icb/ihb` in seeded and restart-temporary `membrane.h5`
  - direct `_project_update_to_symlay(...)` probe showing injected nonzero inner values are projected back to zero
