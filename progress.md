# Progress log

## Current phase (2026-08-10/11) — findings-88 fix deployed; NP dt corrected

### findings-88 fix — COMPLETE, production running
`HybridPositionNode::propagate_deriv` now routes the O-site sensitivity through the placement Jacobian
and BB's fourth weight is no longer dropped. `_write_nonbonded_pairs` restricts the protein side to BB
only (N/CA/C/O env pairs = 0). Seeds rebuilt for all four glpG variants and all six NP faces on the
corrected pair table and binary. Jobs:
- **glpG REMD**: jobs 53233848/49/51/52 (79HIS, 79HIS_S115T, 79ALA, 79ALA_S115T), PENDING as of 08-11.
- **NP**: job 53251911, RUNNING as of 08-11, dt=0.001, block 2.

**Still open:** avg_kinetic_energy/1.5kT is +2.1% above 1.000 after the fix. Second cause unidentified;
not blocking production. `--potential-deriv-agreement` is unusable as a gate on these systems (round-off
dominated at this potential scale).

### NP dt corrected — COMPLETE
Job 53234804 (dt=0.005) had runs 1 and 4 destroyed at t≈250 during protein unfolding. Root cause:
backbone spring bonds (k=48, ω=9.8 rad/t_u) are driven to 3–7× thermal amplitude by conformational
changes; at dt=0.005 the integration accuracy is insufficient to track force reversals at large amplitude.
Reverted to dt=0.001; /output from destroyed runs deleted; resubmitted as job 53251911. See findings.md
"NP-1AO6 blow-up at dt=0.005 after findings-88 fix".

## Environment notes (still valid)
- midway3 default shell is zsh (1-indexed arrays) — use explicit index mapping in upload loops.
- A micelle has no fixed normal (asphericity fluctuates 0.19–0.33); depth and tilt must be measured
  against the aggregate's instantaneous short principal axis, never box z.
- Cancelling a job is not always a no-op — a 9-second run was enough for reseed() to consume /output
  on all six NP configs and strand the next job.
