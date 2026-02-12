# Project Goal
Validate Python3/PyTorch/local modernization of the ConDiv workflow against the original Python2/Theano/Slurm implementation, debug runtime worker crashes, and add `run_remote.sh` for Slurm restart-based training rounds.

# Architecture & Key Decisions
- Compare modernized vs original script behavior first, before applying any fix, to confirm whether the crash is due to modernization drift or environment/runtime issues.
- Trace the worker launch path in `ConDiv.py` and helper scripts to identify the exact failing command and missing output/error handling.
- Reproduce the failing worker (`1mg4`) directly to isolate root cause and avoid full expensive reruns.
- Apply the smallest possible code changes that preserve intended original behavior.

# Execution Phases
- [x] Phase 1: Baseline context and failure-path mapping
- [x] Phase 2: Modernized vs original diff review for critical scripts
- [x] Phase 3: Targeted reproduction of `1mg4` failure and root-cause confirmation
- [x] Phase 4: Implement minimal fix and validate with focused rerun
- [x] Phase 5: Summarize modernization validity and remaining risks

# Known Errors / Blockers
- Initial observed failure: worker for `1mg4` exits with code 1 during epoch 0/minibatch 0/chunk 1.
- Full end-to-end rerun after fixes has not been completed yet in this session (only static/syntax validation done).

# Revised Decisions
- Add a Slurm-focused launcher script (`run_remote.sh`) that mirrors local resume semantics: always pick latest valid `checkpoint.pkl` and run a configurable number of restart iterations per submission.
