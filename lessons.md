# Lessons

- 2026-02-23: For stage-transition workflows, do not introduce coordinate transforms (recentering, per-atom wrapping, carrier remapping) unless explicitly requested. Default to exact previous-stage last-frame handoff.
- 2026-02-23: When user states intended stage semantics ("only update tables/settings"), preserve that exactly and avoid convenience logic that changes physical state.
- 2026-02-23: If the user asks for stage parity (e.g., `6.6 -> 7.0`), validate by checking exact coordinate equality first, then evaluate energy differences as a separate effect of stage parameter changes.
- 2026-02-23: If user defines a canonical workflow entrypoint (here `prepare_system.py`), do not route execution through legacy prep scripts; keep helper code internal and preserve the requested entrypoint contract.
- 2026-02-23: When users ask to reduce script chaos, enforce one operational script entrypoint in workflows and convert legacy scripts into import-only helpers instead of subprocess targets.
- 2026-02-23: If user further requires one script + one lib, merge helper modules promptly and remove transitional helper files after dependency checks.
- 2026-02-25: When user distinguishes "convert to MARTINI interaction representation" from "fallback behavior", do not add fallback paths that mask mapping gaps; keep pair-kernel semantics strict and validate interactions remain within the intended particle model.
- 2026-02-25: When user asks for selective interaction re-enable (e.g., SC-SC but not same residue), implement the minimal rule-gate change in pair-allowance logic instead of broad interaction-policy changes.
- 2026-02-25: For probabilistic coupling requests with directional weighting examples, implement explicit per-side force scaling (counterpart-weighted) rather than assuming symmetric pair averaging.
- 2026-02-25: When user requires state-locked probabilistic coordinates (rotamer-locked SC rows), enforce it in the active interaction-evaluation path rather than relying on large restraint constants.
- 2026-02-25: In probabilistic SC mode, do not leave SC-SC interactions on direct proxy force paths, because proxy forces are later cleared by row projection; route SC-SC through row-space gradient accumulation so interactions survive projection and stay rotamer-locked.
- 2026-02-25: If SC-SC stabilization fails in practice, disable SC-SC at both policy and pair-loop gates immediately, instead of leaving partial SC-SC paths active.
