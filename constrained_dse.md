# Constrained DSE: the unfolded-state target as a constraint in ConDiv (designed 2026-11-11)

Repo plan.md Phase 38 item 20. User 11-11: "so let's say this is a design issue on training framework, what can we do
to improve it?", then "record this in md file, expect to run this on another computer tomorrow afternoon". This file
is the design and the run book; job ids go only into `remote_jobs.md`.

## Resume here (the other computer, 11-11 afternoon)

1. Read this file, then `remote_jobs.md` "Resume here" for what is running (ff3.3's last epochs and its validation
   waiter, the two ff2.1 checks `ff21_coilstart` and `ff21_jumperlr` and their TM reads). This Mac's watcher is a
   background loop of the session that wrote this file; it ends with that session, so start your own watch if the
   user asks for one.
2. Open the midway2 socket once (`expect scratchpad/mdw2_master.exp`, one Duo push, never in a loop).
3. Build the test run as in "The test run" below, in `/project/trsosnic/yinhan/checks/dsecon_20261111/` (exists,
   holding a copy of this file as `DESIGN.md`).
4. Run the checks before submission (no test job), then submit the first link and record it in `remote_jobs.md`.
5. At step 19, read the panel `dc_00` by the rule below, in TM-score and Ca-RMSD only (the folded fraction is not read,
   user 11-10).

## Why: the measured problem

- **The drop is real in Peng's and Jumper's metrics.** Released ff2.1 on Peng's benchmark (16 proteins, his scorer, our
  pipeline) reproduces his FF2: native 0.552, de novo 0.427, against 0.55 and 0.42. Continued by our FF2 trainer at
  lambda 0.3 (`ff21_ctrl_dt009`, driver `8977d7f7`), it scores after 37 steps 0.486 native (-0.066 [-0.135, -0.015])
  and 0.299 de novo (-0.128 [-0.200, -0.065]), below Peng's FF1 (0.45 / 0.37). After 19 steps on the 44 domains from
  the native, mean TM 0.644 against 0.682 (-0.038 [-0.047, -0.030], 42 of 44 lower) and Ca-RMSD +0.45 A (40 of 44
  higher) (`checks/fold_vs_tm_20261110/tm_refs.txt`; findings.md 1.42).
- **The two terms of the objective are not alike.** contrast = NSE + lambda DSE. The NSE (restrained native minus the
  free ensemble) vanishes when the model's native ensemble matches the restrained one: it limits itself. The DSE (our
  SARW minus the pooled unfolded frames near Tm) compares Upside's compact unfolded globules with a self-avoiding
  random walk the model never reaches, so it never vanishes. At ff2.1 its pull on the backbone burial scale is +20
  per protein against the NSE's -0.4; the balance lambda* is 0.018 at dt 0.009 (0.088 at 0.015), not 0.3. Adam then
  moves the scale at its full rate in the DSE's direction each step, and that scalar alone carries the folding loss
  (plan.md items 15-18).
- **Peng's 0.3 does not transfer**: it depends on how his code scales the term, and his trainer is not readable
  (findings 9v). His stated aim is a constrained problem: the DSE "serves as a regularizer that produces the minimum
  interaction energy that can still generate both a stable NSE and a disordered DSE" (JCTC 2022). A fixed weight
  enforces neither condition.
- **The DSE's job is real.** Trained without it, ff3.2 gains mdCATH TM (+0.018 [+0.012, +0.025]) but its unfolded state
  keeps more residual H-bonds than ff2.1's (h +0.112 [+0.038, +0.194], 18 of 24 higher; `unf_hb_ff_3.2_p57_09.txt`).
  Peng's own measure of the DSE's effect is the same quantity: the non-NSE to NSE H-bond fraction at Tm fell from
  38-50% to 20-31%.

## Design

### Objective and constraint

Native stability is the objective; the unfolded state's residual structure is a constraint, h <= h_target. The DSE
weight lambda becomes the constraint's Lagrange multiplier, updated by dual ascent: it rises only while the unfolded
state keeps more H-bonds than the target and returns to 0 once it does not, so it cannot go on trading native
stability for a target the model never reaches. This is the KKT setting for one inequality constraint, the multiplier
non-negative and zero when the constraint is inactive (`~/Documents/mdfiles/math/Fudan_NumericalOptimization.md`
[P19]). It does not depend on Peng's lambda, and it holds under either remaining explanation of the drop (our
objective scaled differently from his, or ff2.1 never trained near its balance point). The DSE's direction is kept
as it is (SARW minus unfolded, every field); only its weight changes.

### h, as the trainer computes it

Per protein with an unfolded ensemble, from the three ensembles the worker already builds (native = the restrained
replica, unf = the pooled bracketing replicas, SARW = ours), with S(x) the sum over hb[0:3] (E_alpha, E_beta,
E_other) of the ensemble's mean derivative of the H-bond energy:

    h = (S(unf) - S(SARW)) / (S(native) - S(SARW))      0 at the SARW, 1 at the restrained native

This is `checks/ff33_ramp_20261109/unf_hb_read.py`'s h (there summed over per-residue derivatives; the engine's
parameter derivative is that sum). The step's value is the mean over the minibatch's proteins with a DSE.

**The target, measured** (`checks/fold_vs_tm_20261110/h_target.py`, `h_target.txt`; released ff2.1 at the control's
step on its minibatches 0-1, 48 proteins, dt 0.009, every system from the native): mean 0.145 [0.116, 0.179],
median 0.099, range 0.039 to 0.649 (skewed, so the constraint reads a running mean, below); the two minibatches 0.150
and 0.139. On the 4 proteins shared with `unf_hb`'s visits: r 0.78, median |difference| 0.019. **h_target = 0.145.**

### The lambda update

After each step's workers return, the driver computes that step's mean h, appends it to the checkpointed history, and
sets the next step's weight

    lambda_next = max(0, lambda + eta * (mean of the last 5 steps' h - h_target))

- lambda_0 = 0 (constraint first; the NSE acts alone until h exceeds the target).
- eta = 0.3 per unit h per step, fixed before the run: an excess of 0.05 (about half ff3.2's measured excess) raises
  lambda by about 0.3 over one epoch of 19 steps, and the window's noise (per-protein spread about 0.12, so about 0.011
  over 120 protein visits) moves it by about 0.003 a step. No variant of eta or the window is tried on any outcome.
- No upper cap on lambda (no arbitrary caps, CLAUDE.md); the max(0, .) is the KKT sign condition, not a guard.
- The driver writes the step's lambda into the step directory (`dse_weight.txt`) before it launches the workers, and
  each worker reads it there, as `ConDiv_ramp57.py`'s `dse_weight_at` reads its schedule; each divergence file keeps
  its weight, its h and its three S values.

### What stays the same

Everything else of the control (`8977d7f7`): the SI's minibatch of 24, 19 steps an epoch, every system from the
native, the 3-replica NSE mix 0.6/0.3/0.1, dt 0.009, the solver and its rates, the initial checkpoint (ff2.1).

## The test run: ff2.1, one epoch (`ff21_dsecon`, panel `dc_00`)

The control with the lambda rule in place of `dse_weight = 0.3`, so the only difference from `c9_00` (the control's
epoch 0) is the rule. One epoch, steps 0-19, `$P/training/ff21_dsecon`.

### Build (template: `checks/coilstart_20261110/`, which built `ff21_coilstart` the same way)

- `make_dsecon.py`: the control's `run_output/ConDiv.py` (`8977d7f7`) -> `ConDiv_dsecon.py` by asserted replacements
  (as `make_coil.py`): `dse_weight = 0.3` replaced by `h_target = 0.145`, `dual_rate = 0.3`, `h_window = 5`; the
  worker reads its lambda from `dse_weight.txt` and records `h`, `S_native`, `S_sarw`, `S_unf` in the divergence
  dict; the driver's step computes the step's h, updates lambda, checkpoints `state['dse_weight']` and
  `state['h_hist']`, writes the next step's `dse_weight.txt`, and prints one line a step (lambda, step h, window h).
- `move_ckpt_ctrl.py` (copy from coilstart): the control's initial checkpoint into `ff21_dsecon`, run paths only, with
  lambda 0 and an empty h history added.
- `build_dsecon.sh` and `after_training_ff21_dsecon.sbatch` (copy and edit coilstart's): at step 19, panel `dc_00`
  (44 domains, dt 0.015, trajectories kept) and a TM read: `tm_panel.py domains runs dc_00 <ref> <C>/tm_runs "" 10 26`
  from `ff3_selection/` with C = `checks/fold_vs_tm_20261110` (both references cached there), and `panel_rmsd.py`.

### Checks before submission (no test job; findings 10.16)

- The driver's diff from `8977d7f7` is the lambda rule and the h record alone.
- The h function, fed the replica-resolved visit's ensemble means (`ff33_cause_20261110/replicas/ff_2.1`),
  reproduces `h_target.txt` (mean 0.145 over the 48).
- The update function on a made-up h series: lambda stays 0 while the window is under 0.145; rises by 0.3 x excess a
  step above it; returns to 0, never below.
- The moved checkpoint equals the control's but for paths and the two new fields; `--test-only` of the first link.

### Rule, fixed before the data (TM-score and Ca-RMSD only)

At step 19, `dc_00` paired over the 44 domains:
- **Passes**: mean TM against `ff21_released_tr` not wholly under 0, and the last 5 steps' mean h at most 0.179
  (ff2.1's upper 95% bound). The constrained objective keeps ff2.1's native stability and its unfolded state.
- **Partial**: mean TM wholly above `c9_00_tr` but wholly under `ff21_released_tr`: the rule removes part of the loss.
- **Fails**: mean TM not wholly above `c9_00_tr`: the weight is not where the loss comes from; the search turns to the
  NSE side or a code defect.
- Reported beside, deciding nothing: Ca-RMSD for each pair, lambda and h at every step, the backbone burial scale per
  step (the control's moved -0.301 to -0.404 by step 19).
- **Read the lambda trace with the verdict.** If lambda stays 0 for all 19 steps, the run is the control with the DSE
  off and h within ff2.1's level: a pass then says the DSE at 0.3 was pressure ff2.1 did not need over one epoch, and
  the constraint itself is first tested where it binds, from ff3.2 (next section).

### Cost and submission

24 tasks of 14 CPUs (12 broadwl nodes) at about 30-35 min a step: 19 steps in about 10-11 h, about 130 node-hours, plus
the panel (one node, 45 min) and TM (one node, about 10 min). Submit through midway2 with
`sbatch $(cat /project/trsosnic/yinhan/slurm/midway2.args) ...`, never to a GPU partition. The user's direction (11-11)
is to run it on the other computer in the afternoon.

## After the test

- **Pass or partial**: the same rule on ff3.2 (`ff32_pcd57` step 80, the driver `c0277e8d` with its persistent chains,
  57 a step), where ff3.2's h is above the target so lambda engages. Its h_target is recomputed on those proteins
  under that driver's protocol before the run; validated with the full suite (mdCATH TM, h, Peng's 32 arms, glpG),
  against ff3.2 and ff3.3.
- **Fail**: report to the user; the remaining leads are plan.md items 16-18 (`ff21_coilstart`, `ff21_jumperlr`, and
  Peng's answers on how many steps at 0.3 made ff2.1 and whether his DSE term is normalized differently from the NSE).

## Other options considered, in order

1. (This file) the constraint with a dual-ascent lambda.
2. A reachable, self-limiting DSE target: the model's own unfolded frames reweighted by maximum entropy to meet the
   target h, so the DSE contrast goes to zero once h is met (experiment-directed maximum entropy: Pitera and Chodera,
   JCTC 2012; White and Voth, JCTC 2014; the max-entropy distribution under a moment constraint,
   `~/Documents/mdfiles/math/CUHK_InformationTheory_Yeung.md` [P58]). Compatible with option 1, more work.
3. A native-stability check every epoch in TM (the 44-domain panel and `tm_panel.py`, under an hour), keeping the best
   checkpoint: protects the selection, does not fix the objective.
4. Removing the part of the DSE gradient that opposes the NSE's (PCGrad, Yu et al., NeurIPS 2020): cheap, but both
   gradients are noisy at 24 proteins a step.

Not done: freezing or capping the backbone burial scale, or tuning lambda by hand until TM looks right. Both hide the
conflict (CLAUDE.md: no twisting parameters, no guards).
