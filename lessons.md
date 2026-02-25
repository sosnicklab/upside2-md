# Lessons

- 2026-02-23: For stage-transition workflows, do not introduce coordinate transforms (recentering, per-atom wrapping, carrier remapping) unless explicitly requested. Default to exact previous-stage last-frame handoff.
- 2026-02-23: When user states intended stage semantics ("only update tables/settings"), preserve that exactly and avoid convenience logic that changes physical state.
- 2026-02-23: If the user asks for stage parity (e.g., `6.6 -> 7.0`), validate by checking exact coordinate equality first, then evaluate energy differences as a separate effect of stage parameter changes.
- 2026-02-23: If user defines a canonical workflow entrypoint (here `prepare_system.py`), do not route execution through legacy prep scripts; keep helper code internal and preserve the requested entrypoint contract.
- 2026-02-23: When users ask to reduce script chaos, enforce one operational script entrypoint in workflows and convert legacy scripts into import-only helpers instead of subprocess targets.
- 2026-02-23: If user further requires one script + one lib, merge helper modules promptly and remove transitional helper files after dependency checks.
- 2026-02-25: When user distinguishes "convert to MARTINI interaction representation" from "fallback behavior", do not add fallback paths that mask mapping gaps; keep pair-kernel semantics strict and validate interactions remain within the intended particle model.
