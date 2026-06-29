# Progress Log

## Current Task: Exact-Step CGL GLE Runtime With Low-Mid Temperature Knots

- Actions taken:
  - Diagnosed the residual temperature trend after the exact-step runtime fix
    and tested three follow-up families against the same reduced-bilayer
    metric:
    broad spectral rebalance (`v13`), full fresh-state GLE stationary
    initialization (`v14`), and auxiliary-only fresh-state initialization
    (`v15`).
  - Rejected both fresh-state source-level initialization variants:
    `v14` oversped the accepted mid/high branch, while `v15` slowed the whole
    branch.
  - Rejected the broad low-temperature spectral rebalance because it oversped
    `T=0.7..0.9`.
  - Kept one targeted improvement from that diagnosis:
    updated the default coarse-hybrid schedule so the validated `T=1.2`
    endpoint now uses the softer `v13` point
    `coupling_scale=1.50,1.50` and `memory_tau_scale=1.70,1.70`,
    while the rest of the retained exact-step schedule stays unchanged.
  - Replaced the split per-mode CGL GLE thermostat update in
    `/Users/yinhan/Documents/upside2-md/src/martini_cg_lipid.cpp`
    with an exact linear-Gaussian step for the coupled physical and auxiliary
    momenta using the Van Loan block exponential.
  - Rebuilt the code and validated that the new exact-step runtime runs cleanly
    on fresh bilayer-only probes.
  - Reused the seed-matched `v7` bilayer inputs as fresh-start exact-step
    probes at `T=0.7`, `0.8`, `0.9`, `1.0`, `1.1`, and `1.2`.
  - Added explicit low-mid schedule knots at `0.8` and `0.9` after the exact
    runtime fix showed that those points were still too dependent on the old
    long interpolation segment.
  - Updated
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_hybrid.sh`
    so the default coarse-CGL schedule now uses the retained exact-step branch:
    grid `0.7,0.8,0.8647,0.9,1.0,1.1,1.2`,
    `coupling_scale=1,1;1,1;1,1;1.013,1.057;1.05,1.22;1.20,1.39;1.55,1.55`,
    and
    `memory_tau_scale=0.33,0.33;0.50,0.50;1,1;0.85,0.95;1.10,1.33;1.35,1.58;1.83,1.83`.
  - Extended the CGL transport metadata path so `/input/cgl_gle` can now carry
    a shared `temperature_grid` plus independent
    `coupling_scale` and `memory_tau_scale`, with either shared 1D schedules
    or mode-resolved 2D schedules.
  - Updated the C++ runtime in
    `/Users/yinhan/Documents/upside2-md/src/martini_cg_lipid.cpp`
    so those scales are interpolated when the simulation temperature is set,
    mode by mode when a 2D table is present.
  - Updated the prep path in
    `/Users/yinhan/Documents/upside2-md/py/martini_prepare_system.py`
    so the new transport metadata can be written from env-backed workflow
    inputs.
  - Updated the coarse hybrid launcher in
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_hybrid.sh`
    so the retained `v7` schedule is now the default coarse-CGL temperature
    calibration.
  - Rebuilt the code and verified the updated Python prep module compiles.
  - Verified the writer path directly by generating a fresh bilayer-only
    `test.input.up` and confirming that the emitted `/input/cgl_gle` datasets
    have the expected 2D shapes and values.
  - Ran multiple reduced bilayer transport probes:
    slow-mode-only scaling, piecewise mode-resolved schedules, and
    seed-matched follow-up scans.
  - Kept the seed-matched `v7` piecewise mode-resolved schedule as the best
    current branch and rejected the `v8` mid-knot refinement.
  - Implemented restart-state preservation for CGL hidden variables:
    compaction coordinate and momentum promotion in
    `/Users/yinhan/Documents/upside2-md/py/martini_prepare_system.py`,
    plus restart-temperature-aware rescaling in
    `/Users/yinhan/Documents/upside2-md/src/martini_cg_lipid.cpp`.
  - Tested and rejected two fresh hidden-state thermalization variants:
    `v9` (thermalize all hidden states) and `v10` (keep zero fresh GLE memory
    but thermalize compaction/orientation momenta). Both disturbed the accepted
    fresh bilayer branch.
  - Retained the minimal `v11` policy instead:
    preserve/rescale only true restart states and leave fresh starts on the
    validated legacy zero-hidden-state path.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/src/martini_cg_lipid.cpp`
  - `/Users/yinhan/Documents/upside2-md/py/martini_prepare_system.py`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_hybrid.sh`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`

- New output artifacts:
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v13_exactgle_spectralbalance/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v14_exactgle_stationaryinit/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v15_exactgle_auxinit/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_transportfix_v1/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_transportfix_tauprobe/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_transportfix_mixprobe_1p2.up`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_transportfix_mixprobe_1p2_c155.up`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v3_0p7_1p2/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v4_midhigh/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v5_modegrid/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v6_piecewise_modegrid/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v7_piecewise_modegrid_seedmatched/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v8_midknot_seedmatched/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_cgl_gle_writer_smoke/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v9_hiddenstatefix_seedmatched/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v10_memzero_seedmatched/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v11_restartfix_only_seedmatched/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_bilayer_temp_scan_transportfix_v12_exactgle_seedmatched/`

- Verification:
  - The working tree is back on the retained exact-step runtime; the rejected
    fresh-state initialization edits were removed and the code rebuilds cleanly.
  - `v13` confirms that a softer `T=1.2` endpoint is safe and useful:
    with the same exact-step runtime, `T=1.2` improves from
    `D_ratio_linear≈0.87` to about `1.02` while keeping
    `sep_late≈16.51 A`, zero flips, and strong orientation.
  - `v14` and `v15` confirm that fresh-state GLE initialization is not the
    next retainable source fix.
  - The new metadata and runtime code build cleanly.
  - The exact-step runtime no longer throws covariance/build errors after
    fixing the Van Loan block assembly.
  - The prep path writes the new 2D transport metadata correctly.
  - The coarse hybrid launcher syntax-checks cleanly with the retained default
    schedule.
  - Structure stays acceptable in all retained transport probes:
    correct leaflet separation, zero flips, strong `n_z`, and uniform `6x6`
    occupancy.
  - Best retained transport scan:
    `v7` keeps `D_ratio_linear` at about
    `0.848, 0.921, 1.009, 0.726, 0.778, 1.082, 0.985`
    over `T=0.7, 0.8, 0.8647, 0.9, 1.0, 1.1, 1.2`.
  - Final retained restart fix does not perturb the accepted fresh bilayer
    branch: matched-seed `v11` reruns at `0.8`, `0.9`, and `1.0` reproduce the
    original `v7` structure and diffusion metrics exactly.
  - Restart promotion checks now pass on copied hybrid and bilayer smoke
    artifacts:
    `input/cgl_orientation_mom`, `input/cgl_gle` restart attrs, and
    `input/cgl_compaction_mom` are all written with restart-temperature
    metadata, and a short promoted dyncomp rerun at `T=0.9` completes.
  - The retained exact-step branch keeps the bilayer structurally correct at
    every tested temperature from `0.7` to `1.2`: late-half leaflet distance
    stays around `16.40..16.52 A`, `sep_max` stays below `18.1 A`, and flips
    stay at zero.
  - Under the same late-half single-origin continuity metric used in the
    earlier scan summaries, the retained exact-step branch now gives about
    `D_ratio_linear≈0.82, 1.07, 0.87, 1.00, 1.00, 0.87` at
    `T=0.7, 0.8, 0.9, 1.0, 1.1, 1.2`.
  - The remaining limitation is now confined to the endpoints:
    `0.7` and `1.2` are still modestly slow relative to the internal linear
    continuation target, while the midrange is materially improved.
