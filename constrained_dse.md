# Constrained DSE: the unfolded-state target as a constraint in ConDiv (designed 2026-11-11)

Repo plan.md Phase 38 item 20. User 11-11: "so let's say this is a design issue on training framework, what can we do
to improve it?", then "record this in md file, expect to run this on another computer tomorrow afternoon". This file
is the design and the run book; job ids go only into `remote_jobs.md`. Amended 11-11 before any data, after a
check that the design does not pick its result (user: "double check that it is not picking results", then "Yes"):
the rule's branches are exhaustive and ordered, a binding criterion decides whether the constraint was tested at all,
the h condition is the controller's check rather than evidence, `dc_00` is a screen, and the follow-up's rule is
written here before it runs.

## Resume here (state 11-11 10:20)

1. **Built, checked and submitted 11-11** (user: "go ahead and implement and test repaired training"), in
   `/project/trsosnic/yinhan/checks/dsecon_20261111/` (README; this file copied there as `DESIGN.md`, the version
   before the amendments as `DESIGN.md.bak_pre_amend_20261111`). Run `$P/training/ff21_dsecon`, driver md5
   `3150a435`, first link to step 19 in `remote_jobs.md`. The checks (`build_dsecon.log`, `static_test.log`) passed.
2. Per step, the link log's `constraint h <= 0.145:` line and the step directory's `constraint.txt` carry lambda and h.
3. At step 19 `after_training.sbatch` submits the panel `dc_00` and `dc_read.sbatch`; read `dc_read.txt`'s VERDICT
   by the rule below, in TM-score and Ca-RMSD only (the folded fraction is not read, user 11-10).
4. **The follow-up is built and gated** (user 11-11: "Do our planned debug assuming he cannot find it"), in
   `/project/trsosnic/yinhan/checks/dsecon_k32_20261111/` (README): run `$P/training/ff32_dsecon57` (driver `4fa75139`),
   its target visit (writes `run_output/h_target.txt`), the gate that reads `dc_read.txt` and starts the run with its
   validation waiter if the rule allows, and the validation at step 160 (`validate/`, dry run passed). Job ids in
   `remote_jobs.md`. If the gate stops, it records why in its `submitted.log`: report the screen's verdict to the user.

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

- lambda_0 = 0 (constraint first; the NSE acts alone until h exceeds the target). Until 5 steps exist the window
  is the steps there are. A step with no protein holding an unfolded ensemble has no h: lambda is left unchanged,
  and that step has no DSE term anyway.
- eta = 0.3 per unit h per step, sized by signal against noise: an excess of 0.05 (about half ff3.2's measured
  excess) moves lambda by 0.015 a step, while the window's noise (per-protein spread about 0.12, so about 0.011 over
  120 protein visits) moves it by about 0.003. A sustained excess of 0.05 therefore takes lambda to about 0.3 within
  an epoch, past the measured balance at ff2.1 (0.02-0.09), so a run can overshoot before h responds; the trace is
  read with the verdict for that reason. eta was kept at 0.3 when the rule was amended (11-11, before any data). No
  variant of eta or the window is tried on any outcome.
- No upper cap on lambda (no arbitrary caps, CLAUDE.md); the max(0, .) is the KKT sign condition, not a guard.
- The driver writes the step's lambda into the step directory (`dse_weight.txt`) before it launches the workers, and
  each worker reads it there, as `ConDiv_ramp57.py`'s `dse_weight_at` reads its schedule; each divergence file keeps
  its weight, its h and its three S values. After the update the driver writes `constraint.txt` in the step
  directory (the weight it used, the step's h and its protein count, the window's h, the next weight), the record
  the readout reads.

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
  (44 domains, dt 0.015, trajectories kept) and `dc_read.sbatch`: `tm_panel.py domains runs dc_00 <ref> <C>/tm_runs ""
  10 <n>` from `ff3_selection/` with C = `checks/fold_vs_tm_20261110` (both references cached there), `panel_rmsd.py`,
  then `dc_read.py`, which applies the rule below mechanically and writes `dc_read.txt` ending in its VERDICT line.

### Checks before submission (no test job; findings 10.16)

- The driver's diff from `8977d7f7` is the lambda rule and the h record alone.
- The h function, fed the replica-resolved visit's ensemble means (`ff33_cause_20261110/replicas/ff_2.1`),
  reproduces `h_target.txt` (mean 0.145 over the 48).
- The update function on a made-up h series: lambda stays 0 while the window is under 0.145; rises by 0.3 x excess a
  step above it; returns to 0, never below.
- The moved checkpoint equals the control's but for paths and the two new fields; `--test-only` of the first link.
- `dc_read.py`'s rule on every combination of the three tests (synthetic intervals and traces), and one dry run on
  real files with `c9_00_tr` standing in for `dc_00`.

### Rule, fixed before the data (TM-score and Ca-RMSD only)

At step 19, `dc_00` paired over the 44 domains, in mean TM-score: A = `dc_00` - `ff21_released_tr` and B = `dc_00` -
`c9_00_tr`, each with its 95% bootstrap interval. Steps are numbered 0-18 by minibatch; w_k is the weight step k used.
`dc_read.py` applies the three tests in order and prints one VERDICT.

1. **The controller's check** (not evidence about the force field: the controller enforces it). Held if the window's
   h after step 18 is at most 0.179, ff2.1's upper 95% bound. If not held, the TM class is reported and the verdict is
   "constraint not held": the controller is diagnosed before the TM is read as a property of the objective.
2. **Did the constraint bind?** Bound if w_k > 0 on at least 5 of the 14 steps 5-18, the steps whose weight comes
   from a full window. Reported beside: the summed weight over the epoch against the control's 19 x 0.3 = 5.7.
3. **The TM class**, exhaustive: *level* if A's interval is not wholly under 0; *partial* if A's is wholly under 0 and
   B's wholly above 0; *none* if A's is wholly under 0 and B's is not wholly above 0.

The verdicts, held and bound:
- level: **passes**, the constrained weight keeps ff2.1's TM over one epoch with the constraint active;
- partial: **partial**, the constraint removes part of the loss;
- none: **fails**, the constraint as built (target 0.145, eta 0.3, window 5) does not prevent the loss; no variant of
  it is tried on this outcome.

Held but not bound, the run is the control with the DSE effectively off, so it says nothing about the constraint:
- level: **constraint untested**, one epoch from ff2.1 with the DSE effectively off keeps ff2.1's TM; the
  constraint's test is the follow-up below;
- partial: **constraint untested**, with the DSE effectively off part of the loss remains: the DSE weight carries part
  of it and the NSE side the rest;
- none: **constraint untested**, with the DSE effectively off TM is not above the control: the loss is not the DSE
  weight's; the search turns to the NSE side or a code defect.

Reported beside, deciding nothing: Ca-RMSD for each pair, w_k and h at every step, the backbone burial scale and
E_alpha at step 19 (the control's scale moved -0.301 to -0.404).

`dc_00` is a screen: a pass or a partial decides only whether the follow-up runs. The claim that the fix works rests
on the follow-up's validation, on suites that decided nothing in this search (the 44 domains have decided about 40
panels since 10-31, the diagnosis of the burial scale among them).

### Cost and submission

24 tasks of 14 CPUs (12 broadwl nodes) at about 30-35 min a step: 19 steps in about 10-11 h, about 130 node-hours, plus
the panel (one node, 45 min) and TM (one node, about 10 min). Submit through midway2 with
`sbatch $(cat /project/trsosnic/yinhan/slurm/midway2.args) ...`, never to a GPU partition. The user's direction (11-11)
is to run it on the other computer in the afternoon.

## After the test

- **Pass, partial, or untested with level TM**: the follow-up, the user's go first. The same lambda rule on ff3.2:
  `ff32_pcd57`'s step-80 checkpoint (`ff_3.2_p57_09`) with its driver `c0277e8d` (persistent chains, 57 a step), the
  rule in place of its DSE weight, from step 80 to 160 as `ff33_ramp57` ran, so the two differ only in how the weight
  is set (the 0.1/0.2/0.3 schedule against the constraint). ff3.2's h is above ff2.1's (+0.112), so the constraint
  binds there. Its h_target is released ff2.1's files (`$P/parameters/ff_2.1`, never ff3.2's) visited once under
  `c0277e8d`'s protocol, by the follow-up driver's own worker, on the 114 proteins of its first two minibatches (as
  the screen's target took the control's first two), every system from the native, the mean over the proteins with
  an unfolded ensemble, fixed before the run in `run_output/h_target.txt`, which the driver reads; item 13's
  `dse_visits` values for released ff2.1 (24 proteins, same protocol) are reported beside it.
- **Started without the user in the loop** (user 11-11, Peng asked for his code: "Do our planned debug assuming he
  cannot find it"). A gate job reads `dc_read.txt` hourly: on a VERDICT of passes or partial, or untested with TM
  class level, it submits the follow-up's first link; on any other verdict it submits nothing and records why.
- **The follow-up's rule, written now (11-11, before any data).** At step 160, item 14's suite against released
  ff2.1: mdCATH mean TM (89 domains), Peng's native and de novo mean TM (paired over the 16 proteins, his scorer), h
  (item 13's 24 proteins, paired). It passes if each TM interval is not wholly under 0 and h's paired interval is not
  wholly above 0. glpG's TM4 (block 1 on Upside 3), ff3.2 (`p57_09`) and ff3.3 (`r33_19`) are reported beside and
  decide nothing; the user chooses among the force fields on these numbers. The 44-domain panels decide nothing.
- **Fail, or untested with TM below ff2.1**: report to the user; the remaining leads are plan.md items 16-18
  (`ff21_coilstart`, `ff21_jumperlr`, and Peng's answers on how many steps at 0.3 made ff2.1 and whether his DSE term
  is normalized differently from the NSE).

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
