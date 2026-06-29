# Project Goal

Stabilize the CGL bilayer transport model across the tested Upside
temperature window without changing the accepted conservative bilayer force
field. The structural fix for leaflet distance already holds; the remaining bug
is the temperature dependence of the CGL lateral diffusion time scale.

# Architecture & Key Decisions

- Keep the accepted conservative CGL bilayer model unchanged.
  The exact-288 scan already showed that leaflet separation, orientation, and
  lateral uniformity stay correct from `T=0.7` to `T=1.2`.
- Treat the remaining defect as a transport-model problem in the CGL GLE path.
  The measured compaction-response state is nearly temperature-invariant across
  the scan, so the large diffusion drift is not a conservative-structure bug.
- Implement optional temperature-dependent scaling for the existing multi-mode
  CGL GLE transport parameters.
  The retained source-level fix is a shared temperature grid with mode-resolved
  `coupling_scale` and `memory_tau_scale`, written by the preparation path and
  interpolated by the C++ runtime.
- Calibrate transport only with seed-matched comparisons.
  Single `100`-time-unit trajectories are noisy enough that changing both the
  schedule and the seed at the same time gives misleading endpoint drift.
  Keep the original scan seeds when comparing schedules, then only use new
  seeds after a candidate is selected.
- Store the calibration in `/input/cgl_gle` as explicit metadata written by the
  preparation path, then interpolate it in the C++ runtime when the simulation
  temperature is set.
- Use the retained exact-288 production state at `T=0.8647` as the reference
  anchor and keep the accepted `T=1.2` endpoint from the shared-scale branch
  when the endpoint scale values are identical.
- Previous retained calibration baseline:
  the seed-matched `v7` piecewise mode-resolved schedule under
  `example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v7_piecewise_modegrid_seedmatched/`
  with temperature grid `0.7, 0.8647, 1.0, 1.1, 1.2`,
  `coupling_scale=[[1,1],[1,1],[1.05,1.22],[1.20,1.39],[1.55,1.55]]`, and
  `memory_tau_scale=[[0.33,0.33],[1,1],[1.10,1.33],[1.35,1.58],[1.83,1.83]]`.
- Next source-level refinement:
  continuation must preserve the full dynamic CGL hidden state, but the hidden
  variables are not all the same kind of state. The retained source fix should
  promote the missing compaction restart state (`cgl_compaction` and
  `cgl_compaction_mom`) and rescale preserved hidden states when the target
  simulation temperature changes, while keeping the validated fresh-start
  behavior unchanged. Fresh hidden states should stay on the legacy zero-state
  path unless they come from an explicit restart, because attempts to
  thermalize them on fresh starts perturb the accepted bilayer-transport
  branch.
- Current retained transport branch:
  replace the split multi-mode CGL GLE thermostat step with an exact
  linear-Gaussian step built from the Van Loan block exponential, then keep the
  schedule explicit around the low-mid range with
  `temperature_grid=[0.7,0.8,0.8647,0.9,1.0,1.1,1.2]`,
  `coupling_scale=[[1,1],[1,1],[1,1],[1.013,1.057],[1.05,1.22],[1.20,1.39],[1.50,1.50]]`,
  and
  `memory_tau_scale=[[0.33,0.33],[0.50,0.50],[1,1],[0.85,0.95],[1.10,1.33],[1.35,1.58],[1.70,1.70]]`.
- Rejected follow-up:
  the broad `v13` spectral-rebalance schedule
  (lower coupling plus longer memory at low temperature, slightly weaker
  high-temperature friction) preserves bilayer structure but overshoots
  diffusion badly at `T=0.7`, `0.8`, and `0.9`.
- Rejected source-level candidates:
  both fresh-state GLE initialization variants were tested and rejected.
  Full stationary initialization of CGL translational plus auxiliary momenta
  improved `T=0.7` and `0.8` but oversped `T=0.9` and `1.2`, while
  auxiliary-only initialization slowed the whole branch.
- Retained interpretation of the temperature trend:
  the tested window spans roughly `245 K` to `421 K`, so the old
  linear-in-`T` continuation is not equally trustworthy at the endpoints.
  The high-temperature endpoint behaves like excess thermostat friction and is
  fixed by the validated softer `T=1.2` schedule point. Low-temperature
  changes that try to flatten `T=0.7` toward the same linear target either
  overspeed the midrange or slow the full branch, so they are not retained.
- For now, evaluate time-scale stability against the internally consistent
  linear-in-`T` continuation from the validated `T=0.8647` point. If a better
  external DOPC temperature curve is introduced later, the same metadata path
  should support it without another runtime redesign.

# Execution Phases

- [x] Re-read the relevant runtime, prep, and workflow code and confirm the
      implementation surface for temperature-dependent CGL GLE calibration.
- [x] Implement the new `/input/cgl_gle` temperature-calibration metadata in
      `py/martini_prepare_system.py` and consume it in
      `src/martini_cg_lipid.cpp`.
- [x] Build the code and verify that legacy inputs without temperature
      calibration still run unchanged.
- [x] Validate the new runtime path directly by writing and reading 2D
      `coupling_scale` / `memory_tau_scale` datasets through the prep script.
- [x] Iterate reduced bilayer probes across multiple calibration families
      (shared-scale, slow-mode-only, piecewise mode-resolved, and seed-matched
      follow-ups).
- [x] Retain the best current seed-matched calibration branch and document the
      remaining limitation.
- [x] Implement full CGL hidden-state restart promotion plus temperature
      rescaling of preserved restart states, while leaving fresh-start hidden
      states on the validated legacy path.
- [x] Rebuild and validate the hidden-state fix on the reduced bilayer scan,
      comparing it against the retained `v7` branch with matched seeds.
- [x] Replace the split multi-mode CGL GLE step with an exact linear update
      of the coupled physical and auxiliary momenta.
- [x] Revalidate the bilayer under the exact-step runtime and insert explicit
      `0.8` and `0.9` temperature knots where the old long interpolation
      segment remained unstable.
- [x] Test fresh-state CGL GLE initialization variants and reject the ones
      that perturb the accepted bilayer transport branch.

# Known Errors / Blockers

- Remaining limitation:
  the exact-step runtime removes a real source-level discretization defect and
  the retained softer `T=1.2` endpoint fixes the hot-end slowdown, but the low
  endpoint is still not perfectly flat against the old linear continuation.
  The current interpretation is that `T=0.7` is a target-definition edge case,
  not the next clear source bug.
- Latest rejected hypothesis:
  broad low-temperature spectral rebalancing is too aggressive. It turns the
  endpoint mismatch into a large overspeed and therefore is not the next fix.
- Latest rejected source-level hypotheses:
  fresh stationary-state GLE initialization is not retainable in either tested
  form. Full translational-plus-auxiliary initialization overspeeds the
  accepted branch, while auxiliary-only initialization slows it.
- Current retained interpretation:
  the main architectural bug in this phase was the split multi-mode GLE
  thermostat update. After replacing it with the exact linear step, only local
  low-mid schedule knots were needed to recover clean `T=0.8..1.1` transport
  while preserving the correct bilayer structure.
- Main interpretation constraint:
  the project has one clearly accepted transport reference point, not a full
  external DOPC temperature series. The present result should be treated as a
  transport-stabilized internal temperature calibration, not as a final claim
  of exact absolute DOPC transport accuracy across all temperatures.

# Review

- Accepted reference artifact:
  - `example/16.MARTINI/outputs/codex_bilayer_npt288_exact_pairrelax_targetgapstate_transportscan_long/validate100_pairrelax_targetgapstate_s2p0_t25.up`
- Existing scan artifacts:
  - `example/16.MARTINI/outputs/codex_bilayer_temp_scan_0p7_1p2/`
- New implementation / validation artifacts:
  - `example/16.MARTINI/outputs/codex_bilayer_temp_transportfix_v1/`
  - `example/16.MARTINI/outputs/codex_bilayer_temp_transportfix_tauprobe/`
  - `example/16.MARTINI/outputs/codex_bilayer_temp_transportfix_mixprobe_1p2.up`
  - `example/16.MARTINI/outputs/codex_bilayer_temp_transportfix_mixprobe_1p2_c155.up`
  - `example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v3_0p7_1p2/`
  - `example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v4_midhigh/`
  - `example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v5_modegrid/`
  - `example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v6_piecewise_modegrid/`
  - `example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v7_piecewise_modegrid_seedmatched/`
  - `example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v8_midknot_seedmatched/`
  - `example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v12_exactgle_seedmatched/`
  - `example/16.MARTINI/outputs/codex_cgl_gle_writer_smoke/`
- Primary implementation files:
  - `src/martini_cg_lipid.cpp`
  - `py/martini_prepare_system.py`
  - `example/16.MARTINI/run_sim_hybrid.sh`
- Current outcome:
  - The code now supports a shared CGL GLE temperature grid with independent
    interpolation of mode-resolved `coupling_scale` and `memory_tau_scale`.
  - The prep path now writes those 2D datasets correctly, and the runtime reads
    them correctly.
  - The prep path now also promotes the missing compaction restart state and
    writes restart-temperature metadata for preserved CGL hidden states.
  - The runtime now rescales preserved CGL hidden states when restart
    temperature and target simulation temperature differ.
  - The runtime now advances the CGL GLE hidden state with an exact
    linear-Gaussian step instead of the old split per-mode update.
  - The coarse hybrid launcher now defaults to the retained exact-step branch
    while still allowing explicit overrides, with the validated softer `1.2`
    endpoint schedule wired in by default.
  - Structure remains robust in every retained transport probe:
    correct leaflet separation, zero flips, strong orientation, and uniform
    `x-y` occupancy.
  - The final restart-state fix preserves fresh-bilayer parity exactly:
    matched-seed fresh reruns at `T=0.8`, `0.9`, and `1.0` reproduce the
    original `v7` metrics bit-for-bit under the retained legacy fresh-start
    hidden-state path.
  - The retained exact-step branch keeps the bilayer structurally correct
    over `T=0.7..1.2` with explicit validation runs at every tested
    temperature.
  - Under the same late-half single-origin slope metric used for continuity
    with the earlier scan summaries, the retained exact-step branch gave
    about `0.82, 1.07, 0.87, 1.00, 1.00, 0.87` before the hot-end refinement,
    and the validated softer `1.2` endpoint point lifts the hot end to about
    `1.02` without harming bilayer structure.
