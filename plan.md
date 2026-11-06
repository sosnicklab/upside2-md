# ff3.0: repair glpG TM4 by fixing glycine's Ramachandran map within ff2.1's own training workflow

## Project Goal

Train ff3.0 from ff2.1 with ff2.1's own training workflow (Peng et al. 2022 SI), modernised only.
The one design change is glycine's. The test is whether glpG's TM4 stays helical in the hybrid
dry-MARTINI glpG system, under the in-training force field (findings 10.19).

**Why glycine.** TM4 (residues 134-151) holds three helical glycines: GLY136, 143 and 149. Upside
applies NDRD's glycine Ramachandran map as a local energy at every glycine. That map is a
statistic over loop sites in folded proteins, so its alpha_L excess is fold selection (glycine is
placed wherever a fold needs alpha_L), not local physics: ln(aR/aL) is -1.19 in NDRD against -0.58
over all residues of the 456 training natives (findings 1.9). Data from folded proteins cannot
separate the energy from the selection (findings 1.15-1.16). So glycine's map comes from a source
where no fold selects anything, and its placement in folds has to come from the non-local terms.

**Fixed by the user.** ff2.1's training data and workflow are kept, so the old and the new force
field compare directly (findings 10.17). A drift is attributed against a control trained the same
way without the glycine change, never fixed by retuning the workflow. TM4 is judged only by the
hybrid glpG test; the soluble selection panel, the H-bond margin and the training set's composition
are not TM4 evidence (findings 10.19).

**TM1** (user, 10-09). glpG's TM1 loses structure under every force field, ff2.1 included
(findings 4.6). It is addressed in the hybrid itself, by correct physics and without retraining
(Phase 13).

**ff3.1** (user, 10-09, after talking with Jumper). Train the relative depth of every basin of every
Ramachandran map, not only glycine's pair, as offsets from NDRD under an explicit L2 prior, with no
new input simulations (route 1 of findings 1.28; Phase 14). ff3.2 (user, 10-15) was the same training
started from ff3.0 (Phase 17). Since 10-26 "ff3.1" means the merge of ff3.1.1 and ff3.2.1, the depth
training under a pooled prior with every other term held at ff2.1 (Phases 29-30); the earlier ff3.1 and
ff3.2 sets are retired as `.bak_pre_merge_20261026`.

**The basin-depth path, not the BioEmu map** (user, 10-11 23:50: "we decided to go with basin depth
training path (as it is expandable to all residues)"). **ff3.0 is ff_3.0_gdepth** (user, 10-12 00:25:
"since ff3.0 is working, it is the only branch we need (and its expansion to all residues)"). Dropped
10-11/12: the BioEmu-map line (ff30_bio, ff30_bio_dt009, ff30_bio_si; released ff_3.0_bio, ff_3.0_bio_si),
the SI-rate twin ff30_gdepth_si (stopped at step 72 before its gate and deleted), the H-bond-frozen arms,
polygly and the local TM4 glycine-flip test; the BioEmu line and the frozen arms were deleted with their inputs
and derived products (remote_jobs.md §1 "Finished runs"). The reason: on the fixed side-chain/lipid table
ff2.1 holds TM4 (findings 4.15). TM4 HDX is closed (user, 10-23; Phase 16).

## Architecture and Key Decisions

**Trainer.** `training/ConDiv.py` is O. Kleinmann's Python 3 port of Peng's FF2 dual-target trainer,
restored to the SI (findings 9v). Per protein and step it runs:
- one native-restrained replica, 12 free replicas at T 0.8 to 1.1 and one SARW replica;
- 8000 time units, with the second half analysed;
- contrast NSE + 0.3 * DSE;
- 456 proteins in 19 minibatches of 24.

Every difference from the port is justified in the file's docstring.
- **dt 0.009** (user, 10-06): Upside's standard step, as in master, the FF1 trainer and the Peng
  benchmark. The port's 0.015 destroyed free replicas, and dt 0.009 has destroyed none in
  2,136 protein-steps (findings 1.21).
- **Side-chain learning rate 10x smaller** (user, 10-02): the 31,420-coefficient pair table
  random-walks under Adam's normalised step (findings 1.17).
- **Learning rates otherwise the port's** (`alpha_scale` 0.5). This makes `hb`, `dhb`, `sheet` and
  the burial groups twice the rate the SI implies (fine-tuning x 0.25; findings 1.23). The SI-rate arm
  (every group but rot at the SI's rate, hb 0.005) lost as much folding: the rates set speed, not the
  loss (Phase 28).

**ff2.1's parameter set** (user, 09-24):
- rot: pair, coverage, hydrophobe;
- burial: 20 x 3 sigmoid parameters and 400 weights;
- the backbone-desolvation scale;
- the three H-bond energies and the second-H-bond term;
- the 20 sheet values.

Left untrained, as in ff2.1, because the engine returns no derivative for them: the backbone term's
center, sharpness and hbond weight, and `hbond.h5` entries 4-11.

**ff3.0's glycine design, gdepth** (user, 10-04/05; findings 1.19): ff2.1's own NDRD library with one
pooled glycine alpha_R / alpha_L depth pair on the 37 GLY|X maps. It starts where the maps hold BioEmu's
alpha_R / alpha_L weights, the only place BioEmu enters, and is trained per epoch by a damped Newton step
on free against restrained populations, with a Gaussian prior of 1 nat on NDRD. The claim tested:
training on the original data alone keeps TM4 stable.

**Control** (user, 10-06): `ff21_ctrl_dt009` has ff2.1's own `rama.dat`, the same trainer, and
target 38 (the half-trained point), with no gate. It separates the workflow's own drift from the
glycine change.

**End of training and release** (user, 10-01/05; end of training revised 11-02). Training runs FF2's fixed
budget of 76 steps and has no convergence test (user 11-02: "we do not need convergence, as stated by Peng
and Jumper"; Peng's SI: 76 iterations, no test; Jumper stopped early on purpose because contrastive
divergence peaks before it converges; findings 1.23). Step 76's checkpoint is the result, unless an earlier
epoch end folds significantly better on the selection panel (44 all-L CATH domains against all-atom; paired,
p < 0.05). mdCATH, Peng's 16 proteins and glpG judge and are never used to choose; the release is the
user's. Runs before ff3.2 (`ff32_dse0`) were stopped by `convergence_gate.py`, a fixed-point test that
trained on to step 114 or 13 epochs; their records stand.

**ff3.2 is the target; ff3.0 and ff3.1 are transitional** (user 11-05: "ff3.1 and ff3.0 no longer matters, if the
ff3.2 validation is positive. they are transitional"). ff3.2 (`ff_3.2_w32_03`) is judged against released ff2.1 alone,
by Phase 38 item 10 (c): mdCATH TM, Peng's TM and glpG's TM4. If that validation is positive, ff3.0 and ff3.1 drop out:
their unfinished work (ff3.1's Peng scoring and the scorer fix it needs) is not done and their columns are not
reported. Until the verdict they are the fallback, so their running arms and kept replicas stay.

**Retired round-4 arms and tests** (detail in findings):
* `bio`: glycine's maps fitted bottom-up to BioEmu's plain-MD octapeptides by iterative Boltzmann
  inversion in Upside (`parameters/common/rama31.dat`), frozen in ConDiv (findings 1.19).
* H-bond and sheet frozen (user, 10-07; `ff30_bio_fz`, `ff30_gdepth_fz`, `ff21_ctrl_fz`): stopped 10-09
  23:09, once the dt 0.009 training kept TM4 stable; they folded less than their twins (findings 1.21, 1.37).
* SI-rate arm (user, 10-07; `ff30_bio_si`, `ff30_gdepth_si`): above.
* The local TM4 test (findings 1.22, 1.25, 3.11): glpG 79HIS in the POPE/POPG hybrid seed patched by
  `patch_glpg.py`, T 0.80, 4000 tu, dt 0.009, 12 seeds; primary readout each seed's last-block DSSP alpha of
  135-151, two-sided Mann-Whitney at p < 0.05 with a bootstrap interval (user, 10-09 03:20). Recipe:
  remote_jobs.md "Round-4 watch" step 4. Retired 10-11 23:10 (Phase 11).

**Retired designs** (one line each; detail in findings):
* Full-map glycine row (Phase 2): failed the convergence gate at every checkpoint.
* Per-pair basin offsets on the NDRD maps, the ff_3.0 released 09-30: favoured alpha_L at every
  glycine; noise-limited per pair (1.12, 1.14-1.15).
* Physics glycine map (AWH dipeptides) plus glycine H-bond offsets, round 3 (`ff30_glyhb`): panels
  lost helix; the offsets went left-handed (1.17).
* Ting et al.'s product combining rule: breaks GLY|GLY symmetry (1.13).
* A design through H-bond terms for glycine (user, 10-04): a map defect is fixed on the map.
* Pre-proline correction: parked; no simulation failure behind it (1.16).

## Execution Phases

### Done or closed (detail in findings and progress.md)
- Phase 1: the trainer validated on ff2.1 (09-24; 9v).
- Phase 2: the full-map glycine row, cancelled 09-28 (9w).
- Phases 3-5: the basin-offset ff_3.0, released 09-30 and withdrawn the same day (1.14-1.15).
  - Its glpG validation ran on live seeds that kept the FF1-form ff_3.0's coverage tables (3.11),
    so its "TM4 glycines flipped" result is of a mixed force field.
- Phase 6: pre-proline, parked.
- Phase 7: the glycine handedness probe. From equal depth the data pull glycine back toward alpha_L
  (1.15).
- Phase 8: round 3, `ff30_glyhb`, stopped 10-04 (1.17).
- Phase 9: repo cleanup for redistribution, 10-01/02. `training/` has five files; the cluster keeps
  the old layout (remote_jobs.md "Where things are").
- Phase 10: the local glpG test of round 3. Its TM4 counts used pre-fix inputs and are invalid
  (3.11).

### Phase 11 - round 4 (BioEmu and frozen arms dropped 10-11, user)
Round 4 trained gdepth, bio, the frozen and SI-rate arms and the control on midway2 broadwl from 10-06 at
dt 0.009 (the dt 0.015 runs were stopped and restarted, otherwise identical). Record: findings 1.19-1.26,
1.31-1.32.
- **ff30_gdepth_dt009 converged and is ff_3.0** (gate 10-09 08:54, findings 1.26; user, 10-09 10:40: "pull
  ff30_gdepth_dt009 down locally as ff3.0"). `parameters/ff_3.0/` holds the six released files, its own
  `rama.dat` included, md5-identical to `$P/parameters/ff_3.0_gdepth` and to `d9_03`. `martini.h5` and
  `membrane.h5` are not trained and stay ff_2.1's.
- ff30_bio_si (gate 10-10 02:25) and ff30_bio (10-10 16:46) converged and were released as ff_3.0_bio_si and
  ff_3.0_bio, then dropped with the BioEmu line (findings 1.26, 1.31, 1.32). The frozen arms stopped 10-09
  23:09; their TM4 and panel results stand (findings 1.21-1.25).
- Built for round 4 and still used: the automatic validation (gate, `validate_ff.sh`, `after_training.sbatch`),
  `patch_glpg.py`'s coverage fix (findings 3.11), the per-step checks (`check_step.py`, the KE scan, the
  glycine readouts) and midway2's node exclusions as one law (`/project/trsosnic/yinhan/slurm/midway2.args`,
  remote_jobs.md §0d).
- The local TM4 test (findings 1.24-1.26): no in-training checkpoint is resolved from its start; d9_02 leaned
  highest (0.765 against gdepth_start's 0.657, p 0.58) and the released d9_03 is 0.726 (p 0.64); bz_02 against
  its twin b9_02 was the only resolved primary test, just (0.856 against 0.694, p 0.046, one of 32 tables). At
  12 seeds a true 0.10-0.15 difference is detected 18-31% of the time (findings 1.25), so a lean in one 12-seed
  set is a reason for 24 seeds (c9_00's two halves differed at p 0.009). On the fixed side-chain/lipid table
  ff2.1 itself holds TM4 (0.944 against 0.610 on the old table, p 0.002), every start reads alike and no
  released candidate improves on it (findings 4.15); the user retired the test 10-11 23:10 ("we have a working
  ff and we are trying to expand it").

### Phase 12 - poly-Gly: all-atom Ac-(Gly)20-NHMe as Upside's reference (STOPPED 2026-10-11, user; started 10-05, PI's suggestion)
With neutral caps the chain is achiral, so alpha_R = alpha_L exactly at equilibrium: a reference for the
glycine map's own handedness. amber99sb-ildn / TIP3P at 300 K on midway2 (files and README in
`/project/trsosnic/yinhan/polygly/`; Upside runs and analysis in the Mac's `scratchpad/polygly/`): collapse,
7.5 nm dodecahedron production boxes, then a planned 4 x 1 us. Upside G20 under ff21_released and bio_start is
chiral, ln(aR/aL) -0.72 and -0.45 (findings 1.20). Stopped by the user at ~97 ns per replica (10-11 23:05): the
question it served, the glycine map behind TM4's glycine flips, was settled by the side-chain/lipid table fix
(findings 4.15). The data are kept (`polygly/`, 213 MB), and the chain resumes from `prod.cpt` if the glycine
map's own handedness is ever taken up.

### Phase 13 - glpG TM1: hold it in the hybrid by correct physics (STARTED 2026-10-09, user)
Hold glpG's TM1 in the hybrid without retraining the force field (user, 10-09): no scaling or disabling of
BB-env or SC-env, no tuned parameters, no guards. Diagnosis (findings 4.6-4.7): TM1's middle (32-38) sits at the
bilayer centre and loses structure in every force field, ff2.1 included; lipid gain is not the driver. Opening it
costs only Upside's soluble-calibrated H-bond energy, since the hybrid has no cost for an unpaired backbone NH or
CO in the acyl core and its MARTINI backbone typing is fixed at the starting secondary structure. Local work on
the Mac Studio in `scratchpad/glpg_tm1/` (from midway2's `checks/tm1_hbmem_20261009/`).
1. [x] TM1 readout: each seed's last-block DSSP alpha of 29-48, two-sided Mann-Whitney against the 12 local d9_03
   seeds (`ff3_local_test/runs/`), resolved at p < 0.05, with a bootstrap 95% interval.
2. [x] The impulse localised (findings 4.9-4.10): the side-chain/lipid table (`martini.h5` `sc_table`) was
   sampled on 24 directions and one bead frame, and its one-node ridges give single-evaluation spikes up to
   33,260 E_up/A on CA34/CA35. Found with a read-only probe build (`scratchpad/glpg_tm1/src_probe`, `obj_probe/`).
3. [x] Fixed at the root: the table at converged azimuthal sampling (600 directions x 12 frames,
   `kick_ff30/sc_table_vectorized.py`); in production since Phase 15.
4. [x] Fixed ff3.0, 12 seeds (findings 4.13): TM1 0.931 against 0.685, p 0.002, resolved; TM4 0.924 against
   0.726, p 0.10; no evaluation above 1000 E_up/A.
5. [x] Option A (user, 10-09 13:45): Upside's membrane backbone H-bond term, its H-bond-state-dependent part only,
   through the engine's existing `hb_membrane_potential` node (half-thickness 15.9 A, the GL1/GL2 mid-plane of our
   own bilayer; `scratchpad/glpg_tm1_hbmem/add_hbmem.py`). On ff3.0 and the fixed table (findings 4.11, 4.14) TM1
   0.964 (p < 0.001) and TM4 0.959 (p 0.009); against the table alone TM1 +0.033 (p 0.049), TM4 p 0.21. ff2.1 +
   term's tears were the table's kicks on MET142 (findings 4.8). Not adopted (user, 10-10 11:20). Option B, each BB
   bead's MARTINI type following its H-bond state, would charge about +38 E_up for the opening, and the literature
   forbids combining it with A (findings 4.12).
6. [ ] Delete the probe build and its source.

### Phase 14 - ff3.1: every basin depth of every map, under an L2 prior (STARTED 2026-10-09, user)
Route 1 of findings 1.28 (user, 10-09: no new simulations for input data). Each coil map (central, direction,
neighbour) gets one offset per basin of `rama_basin.py`'s six (alpha_R, alpha_L, beta, pPII and the two phi > 0
extended regions; five free per map after renormalisation), trained on that map's own reads under its own Gaussian
prior centred on NDRD, with no pooling across maps (rule 1.8), by the damped Newton step with the prior, the MAP
(ridge) estimate, not Adam (findings 1.28: Adam moves noise-only coordinates at its full rate). Risks recorded
before the probe (findings 1.28): the depths can learn native placement and absorb H-bond and side-chain physics,
which a prior does not prevent; 46 held-out proteins sit at the sample-size floor (1.12), so the width needs
cross-validation over many folds. This set is retired: ff3.1 is now Phase 30's merge.
1. [x] Probe on ff21_ctrl_dt009's recorded basin populations, no job (findings 1.29): cross-validation finds a
   finite width (0.3 nat), but per-map depths transfer less well than pooled ones, and what transfers is a global
   alpha_R deepening.
2. [x] Trainer `training/ff31/` (findings 1.30).
3. [x] `ff31_depth_all` (5,040 depths, prior 0.3) and `ff31_depth_data` (1,442, prior 0.2) from ff2.1 with every
   offset at 0 (user, 10-09 23:05), target 76, no gate.
4. [x] Reviewed 10-13 (user; findings 1.34): `all` is worse than ff_3.0_gdepth (folded -0.058, paired), `data`
   ties it (-0.012, not resolved); the depth rounds closed 1.6% (`all`) and 1.3% (`data`) of the held-out gap
   (findings 1.30).
5. [x] **Steric watch** (user, 10-10 21:10, after Jumper: a moving map may open conformations where the side
   chain and backbone clash; findings 1.33). Upside has no excluded volume between a residue's CB and its
   neighbours' backbone, so the map is what keeps them apart. `training/ff31/steric_check.py` (deployed
   `checks/ff31_steric/`) scores each step's library against a hard-sphere Ramachandran map and reads the free
   ensembles' non-glycine phi > 0 populations. No flag but `data`'s step 39, shown to be scatter.
6. [x] fd_03 validated as `ff_3.1_depth_data` (user, 10-13 12:00: "Use Peng validation set and glpG TM4 to
   validate ff3.1"; `all` is not validated): TM4 helical in all four glpG chains (findings 1.35), Peng arms scored
   (findings 1.37). The chains were rerun after the placement fix (Phase 18); the set moved aside as
   `.bak_pre_merge_20261026` (Phase 30).

### Phase 15 - the side-chain/lipid table at converged sampling, in production (DONE 2026-10-10, user)
The user's word (10-10 08:30): apply Phase 13's table fix to the force-field files so future simulations run
with it. The fix is a sampling density, not a parameter: the same potential, Boltzmann average and grids, at 600
target directions x 12 side-chain bead frames instead of 24 x 1 (findings 4.10, 4.13, 4.15; up.md 2.6). The
membrane H-bond term is not part of it. Done:
- `py/martini_build_tables.py` computes the per-sample loop as arrays, 600 x 12 by default; at 24 x 1 it
  reproduces the old table.
- `parameters/ff_2.1/martini.h5` rebuilt (md5 `2e612261...`, old kept as `.bak_pre_sc600x12_20261010`); since
  rebuilt by the audit (Known Errors; Phase 24 item 5.6).
- `patch_glpg.py` (md5 `97ca3294...`) copies martini.h5's table into a seed; TM4 sets ran as `<tag>_sc`.
- Cluster: martini.h5 in `$P` and `/beagle3/.../upside2-md`, the builder in `/beagle3`, `patch_glpg.py` in
  `$P/training/` and `tm4_local/scripts/` (`checks/sc_table_20261010/`); `$P/py/martini_build_tables.py` stays
  old (campaign rule).
- TM4 references on the fixed table: ff2.1 holds TM4, every start reads alike, and none of the three released
  candidates is resolved from its start (findings 4.15).

### Phase 16 - glpG TM4 HDX on the fixed table (CLOSED by the user 2026-10-23)
The only pre-fix hybrid HDX set (ff2.1, 09-04) gave TM4 3.5-4.9 kcal/mol against the H145A peptide 140-144's
6.1-7.4, through protein protection failing (`pp_fail` 0.032 on 135-151 at the T 0.85 rung, findings 5.3); it
predates the 09-10 temperature fix and the 10-10 table (findings 5.5). Residue numbers here are HDX labels,
PDB - 1 (findings 10.38). Closing rule (user, 10-14): resolved if labels 137-147 (PDB 138-148) stay censored
above ~5.1 at 298 K and 148-151 fray no more than locally; HXMS stays a sampling limit. The HXMS peptide's
amides 142-144 are labels 141-143.
- Local ff2.1 on the fixed table, T 0.85, 12 seeds (`scratchpad/tm4_hdx_local/`; findings 5.5): TM4 137-148
  lose protein protection in at most 0.06% of frames at 298 K, against 10-25% in the pre-fix set, and 142-144
  never exchange (dG > 5.8, the frame-count bound). The pre-fix lability does not reproduce on current inputs.
- ff2.1 REMD on the fixed table (`glpG-RKRK-79HIS.ff_2.1_sc`): cancelled by the user 10-14 01:11 in block 2;
  files kept, HDX not run (remote_jobs.md §1 "Finished runs").
- Block-2 REMD readouts of ff_3.0_gdepth and ff_3.1_depth_data (findings 5.5; tools in
  `checks/hdx_ff30_20261014/`, recipe remote_jobs.md §1 "TM4 HDX on ff_3.1_depth_data: the recipe"): run on the
  misplaced-bead chains (findings 4.17), so replaced by the repaired chains below. Block 2's walkers do not cross
  the ladder (each rung's frames come from the few walkers near it), and a finite TM4 dG, as HXMS gives, needs
  ~3e5 independent frames per rung.
- The placement-fixed ff3.0 chains, blocks 2-3 (Phase 18; findings 5.6), in PDB numbers at 298 K:
  - The peptide's 142-144 are censored above 6.25 in 79ALA, the experiment's variant, and in 79ALA_S115T, in the
    10-22 figures and with the A61 accessibility shift undone (`hdx_accfix.sbatch`): the model agrees with
    HXMS's 6.1-7.4. In the 79HIS variants 142 resolves at 5.43 and 4.00 once the shift is undone.
  - The closing rule fails in every variant once the shift is undone: the 138-140 turn, Val145 (4.02 in 79ALA)
    and 148. No HXMS peptide isolates those residues (143-154 mixes nine amides).
  - The core is censored only through the membrane term (PS_protein: 79ALA's 142-144 at 1.8-4.6), on an engine
    without A83's ALA and GLY side chains.
- **Closed by the user 22:15 10-23** ("Ok. Close the question"), after being told (and after a first misreport
  of label 144 as a peptide amide) that TM4's dG is not expected to come down and on the measured peptide does
  not. The rule's fraying and the weak protein protection move to the Upside 3 glpG run (Phase 24 2b 4), which
  has A83 and the other hybrid fixes. No corrected figures are drawn: the 10-22 set carries the A61 shift, and
  its seven amides before a proline or the chain end would need their accessibility recomputed. HDX on ff3.1's
  repaired chains was dropped by the user 22:50 10-23 ("Yes. Drop it").

### Phase 17 - ff3.2: ff3.1's depth training started from ff3.0 (STOPPED 2026-10-18 in step 69; rerun in Phase 21)
The user's word (10-15): "same as ff3.1, except that it starts with ff3.0 instead of ff2.1". It was
ff31_depth_data's run (mode `data`; trainer byte-identical, ConDiv.py `df3e2497...`, rama_basin.py
`235a1bf7...`) with ff_3.0_gdepth's six released files as the start, in `$P/training/ff32_depth_data`
(`checks/ff32_deploy_20261015/`). The trainer takes its base library from `upside_input/rama.dat` and its Adam
groups from `init_param/`; the offsets start at zero on ff3.0's library and each depth's prior is centred there,
so the trainer's printed "NDRD" and the steric check's "against NDRD" mean the base library. Stopped by the user
08:33 10-18 in step 69 for the angular-gradient defect and rerun as `ff32_depth_data_angfix` (Phase 21); its
readouts to step 68 and panels f32_00-02 are the defect side. Retired with the old ff3.1 set
(`.bak_pre_merge_20261026`).

### Phase 18 - the hybrid side-chain placement defect: fix and rerun glpG (STARTED 2026-10-17, user)
Findings 4.17: the hybrid preparation stored `placement_data` gathered per row but `layer_index` as the absolute
library row, so the engine placed every side-chain bead at `rcf[li[li]]`. User 10-17: "fix the bug and rerun the
affected simulations". Decisions:
- **The layout is upside_config's (and master's):** `placement_data` is the full library table (86 rows; the
  1-body scalar likewise) and `layer_index` the absolute row. The engine is correct and does not change.
- **Repair, do not rebuild:** only the two `placement_data` datasets change, so membrane, packing, ions,
  coordinates and every other dataset stay bit-identical. The repair recomputes the full tables from the
  build-time `ff_2.1/sidechain.h5` by the preparation's own rule and refuses unless the seed's gathered rows
  equal them exactly.
- **The base seeds are repaired too** (`popepopg_REMD_mdw2/seeds/glpG-RKRK-<V>.up`, backups
  `.bak_pre_placefix_20261017`), so `validate_ff.sh` patches correct seeds. The deterministic seed (79HIS,
  ff_2.1_nobp) keeps one row per side-chain residue at its kept state, the side-chain/lipid and backbone/lipid
  couplings on and unchanged in form.

Done (`checks/placefix_20261017/`: `fix_placement.py`, `single_state.py`, `verify_seed.py`, `build.sh`,
`verify_all.sh`, `verify.log`): the fix in `py/martini_prepare_system_lib.py` (`load_sidechain_rotamer_payload`) and
its `$P` and `/beagle3` copies; placement at `rcf[li]` to 3.6e-6 A (old: up to 7.9 A), 1-body exact (old: up to 22
E_up off). All 10 chains were rerun (ff3.0 and ff3.1 in the 4 variants; 79HIS under ff_2.1_nobp and ff_2.1_ctrl38,
stopped 10-18 in Phase 21); the old chains are kept as `<V>.<ff>.placedefect_20261017`, their trajectories deleted
10-18 (user; HDX results kept). The 8 ff3.0/ff3.1 chains ran all five blocks: block 1's TM4 holds, with no glycine
at phi > 0 at T 0.70 (findings 1.35); ff3.0's HDX on blocks 2 and 3 is Phase 16's closing readout (findings 5.6:
the blocks disagree beyond their error bars, so a longer run is the user's decision).

7. [x] TM4 over blocks 1-5, read 11-04 on ff3.0's four chains (the old ff3.1's were deleted 11-02; findings 1.35):
   frays in the 79ALA backgrounds, holds in the 79HIS ones. ff3.0's trajectories deleted after (remote_jobs.md).

### Phase 19 - a second native-state panel from mdCATH (DONE 2026-10-17, user)
The Peng benchmark's 16 natives and the 44-domain panel (60 proteins) do not resolve ff_2.1_nobp from its
control; ~75-100 native-start proteins are needed. User 10-17: "find something new into the validation set";
time step 0.009, as trained. Decisions:
- **The 90 domains** of `ff3_selection/mdcath_candidates/` (selection rules in its README: mdCATH CATH 4.2.0
  S20 domains, 30 per length bin 50-100 / 101-150 / 151-250, fold classes balanced, unrelated to the training
  and panel proteins, folded at 320 K in all 5 all-atom replicas). 89 are used: 1gqoB00 is excluded
  (`excluded/README.txt`: mdCATH's stored starting structure has C-N 3.1 and 4.8 A at its C-terminus).
- **The all-atom reference is mdCATH's 320 K** (CHARMM22*, TIP3P), not Charron's amber99sb-ildn at 300 K:
  paired comparisons between Upside force fields hold; absolute folded fractions are not compared with the
  44-domain panel's.
- **A separate panel** in `ff3_selection/{domains,aa,runs}_mdcath/`. `panel.py` gains an optional time step on
  `run` (default 0.015, so every existing and automatic panel runs as before) and an `aa_mdcath` builder;
  `panel_mdcath.sbatch` beside `panel.sbatch`.

**Result (findings 1.21):** the deterministic side chain loses folding on held-out proteins. Over the 89
domains no-BP folds 0.264 against 0.365 (-0.101 [-0.155, -0.051], Wilcoxon p 7e-5; lower in 59 domains, higher
in 19), mean Q -0.035 (p 2e-7), CA-RMSD +0.21 A (p 0.004), in the same direction in every length bin and fold
class. ff_2.1 keeps the probabilistic side chain. Open: whether the 44-domain panel's -0.016 differs from this
because of its proteins or its time step 0.015; not needed for the decision.

### Phase 20 - audit the no-BP training protocol, then close the branch (STARTED 2026-10-18, user)
User 10-18: no-BP performs below ff2.1; check that the training protocol was correct and that its training has
ended. Read only (findings 1.21 "Audit"): the single-state library, configs and releases are correct; the
engine and trainer at n_rot 1 are correct, derivatives to float precision; the run matches ff21_ctrl_dt009 but
for the side chain; the mdCATH verdict recomputes exactly. Found on the way: the shared angular-gradient defect
(findings 3.12; Phase 21). The gate failed no-BP's `rot` group at step 38, so the verdict held at matched
length, not at no-BP's fixed point; answered by Phase 21, which trained both to their fixed points (Phase 25).

### Phase 21 - fix the angular-gradient defect and retrain what it touches (STARTED 2026-10-18, user)
User 10-18: "fix the issue and continue training if needed. also check if other on-going training also have the
defect, resolve them as well. re-run the validations if necessary." Decisions:
- **The fix is in the engine** (`src/interaction_graph.h`; findings 3.12): the derivative added to [t2,t1] is
  the edge's parameter derivative evaluated at p[t2,t1] with the elements exchanged, in place of the [t1,t2]
  derivative copied unchanged. Fully symmetric parameter sets get bit-identical derivatives; only the pair term's
  angular entries change. Energies and forces do not use this path.
- **The trainer's uniform 2x on pair coefficients stays.** It sets the regulariser's relative weight and has
  since ff2.1's own training (master's `rotamer_parameter_estimation.py` builds the same layout).
- **No-BP and its control continue from step 38 on the fixed engine, each to its own fixed point** (the gate at
  step 76, up to 6 epochs), then the mdCATH panel at dt 0.009. No new release names in `parameters/`.
- **ff3.2, the only training running, was stopped in step 69 and rerun from step 0; ff3.1 (`ff31_depth_data`)
  was rerun too**, each differing from its source only by the fix. ff3.0 and ff2.1 were trained with the defect
  and are not retrained unless the user asks.

Done: the fix built on midway2 (`checks/angfix_20261018/obj`, md5 libupside.so `3d9c6366...`, upside
`7136d99a...`), bitwise equal to the old build in energy, forces and every other derivative ConDiv reads
(`parity.log`, `parity_glpg.log`), each angular entry against tied central differences ratio 1.000 (`fd.log`),
and installed in `$P` 08:35 10-18 (backups `.bak_pre_angfix_20261018`). Reruns `ff31_depth_data_angfix` and
`ff32_depth_data_angfix` (`checks/angfix_20261018/build_rerun.sh`). Item 7, the validations: no-BP and the
control converged at step 76 (10-19, 10-20); `mdc_nobp_03` folds 0.241, `mdc_ctrl_03` 0.353, paired -0.111,
p 1e-5 (findings 1.21; Phase 25). The reruns' epoch-end 44-domain panels are in findings 1.37 (Phase 27 item 5).

### Phase 22 - why de novo folding fails on large proteins (STARTED 2026-10-19, user)
User 10-19: Upside folds large proteins de novo poorly; investigate, read the literature and other force fields'
designs, propose, run any test needed. The user's hypothesis: folding is co-translational (S. Fried, JHU: about
a third of E. coli proteins do not refold in vitro, more often the multi-domain ones), so a thermodynamic force
field cannot fold them from a denatured chain. The 456 training proteins are 50-150 residues, the 44-domain panel
50-75, the Peng de novo set at most 92. De novo failure is thermodynamic (the native basin is not Upside's
free-energy minimum, so no route reaches it) or kinetic (it is the minimum, but bulk search, or a missing
co-translational route, does not find it); each test separates the two.
- Items 1-7 (findings 11e): the mdCATH 151-250 bin's lower folded fraction (ff_2.1_ctrl38 0.24 against 0.42
  and 0.44) is the cutoff, not drift; on 2LZM, 2RN2, 1YOB and the 1ENH control at T 0.80 the compact non-native
  globules sit 150-180 E_up above native and vectorial growth equals bulk; no coil folds at T 0.88 or 0.95; the
  proposal went to the user 10-19 06:50.
- Item 8: `environment_coverage_hb`'s covering-residue types fixed (user 10-19: "Fix the bug"; findings 3.13):
  the writer stores one type per residue, as `environment_coverage_sc` does. A file written before the fix runs
  bitwise as master; a newly written file departs from master by intent. Repo and local `obj/` only: ff_2.1 and
  every training so far ran with the defect.
- Item 9: `upside3.md` holds the findings and the proposal.

**Review.** For these in-vitro refolders Upside's de novo failure is search into compact non-native globules
that it ranks well above the native; N-to-C growth does not help. Proposed, for the user's decision: compact
misfolds as negative examples in ConDiv, training past 150 residues, and replica exchange for large-protein
benchmarks. Not measured: the globules' free energy against the native's, and a size contrast under one protocol
(the 1ENH control does not fold at constant T in 60,000 tu); both need replica exchange.

### Phase 23 - plan Upside 3 (STARTED 2026-10-19, user)
User 10-19: Upside 3 starts after the side-chain verdict, with five goals: a debug of the whole repo, all bugs
found by us as Jumper did on master (this branch, master's code first); unit alignment; new parameters chosen
after the verdict; large proteins (Phase 22); GPU support with one force field for both backends. The plan, the
order and the decisions waiting on the user are in `upside3.md` (research reports in
`scratchpad/upside3_planning_20261019/`): defect inventory and method (2); units (3: E_up's kJ/mol value is FF2's
temperature calibration, 298 K / 0.85, so A+ is recommended: reduced units inside the engine, standard units at
every input and output); GPU (6). Decided 10-19: Upside 3 is a new release, so no list of fixes and their effects
on ff2.1 is published; each fix is proven correct, but none is measured against ff2.1 or master, and ff2.1 is not
re-validated on the audited engine (the user's 10-22 validation of the debug runs ff2.1 on it, Phase 24 2b).
Decided 10-29: the GPU "should produce the same results as the current ffs, without retraining or training another
version of ff". It runs the current force fields from the same files, to the audited CPU engine's results, and
computes energies and forces only (`upside3.md` 10.2); it no longer waits on the training changes or the side-chain
choice.
4. [ ] The user's decisions (`upside3.md` 11); the audit is Phase 24.
5. [ ] **Training objective terms for HDX, then large proteins, then IDRs** (user 10-28, in that priority; design and
   steps in `upside3.md` 8.2). Terms added to ConDiv's contrast from the per-frame dV/da the engine already returns; the
   energy's functional form stays spline tables.
   - [x] H1, the one measurement the midway2 data allow (10-28, findings 5.7): ubiquitin on the Peng-benchmark REMD.
     The pipeline reproduces Peng's FF2. A full benchmark REMD resolves 14-15 of the 24 measured amides at 273 K, and
     ff2.1 opens the globally exchanging amides within the folded state. ff3.0's ubiquitin melts about 0.06 T_up below
     ff2.1's.
   - [ ] H2-H5 (`upside3.md` 8.2; the large-protein and IDR baselines are 8.3.5 and 8.4.6). H2 reads ff3.2's ubiquitin
     arms of its Peng validation, whose replicas are kept until H2 is read (ff3.1's native arm's replicas were
     deleted 11-04 after scoring, so ff3.1 has none); H3 and H4 need new runs or data; H5 waits on H3, H4 and the training temperature
     mapping (`upside3.md` 6.3, 11 item 9).
6. [x] **`upside3.md` restructured into one ordered design** (user 10-29: "make it logical and ordered, and fill in
   the gap between steps"; the approved plan is `~/.claude/plans/mighty-prancing-candle.md`). Stages S1-S7 and track
   G with gates; cross-cutting standards, the temperature scale and the base force field each in one place; labels
   H1-H5 (HDX), L1-L6 (large proteins), B0-B2 and I1-I3 (IDRs); every gap filled from the records or written as a
   proposal for the user; cross-references in the other records updated; the channel-potential steps measured in
   dynamics.
   **Review (10-29).** Every number of the old text kept except the stale code line numbers, now cited by function.
   Gaps filled from the records: no Upside 3 training tree yet (every chain on the pre-audit engine), no base or
   side-chain choice, the trainer's loss (findings 1.37) against the Phase 29 rule, Peng's S14-S16 calibration and
   its overlap with the HDX test set, Ferrari's per-domain data, the debug validation's results. New proposals for the
   user (`upside3.md` 11): S1's gate, accumulation and acceptance rules, the training temperature mapping, S3's gate
   for co-training changes, the refusal design; the GPU's purpose, answered 10-29 (above). The lateral-pressure
   proposal Phase 24 item 5.1 pointed to was missing and is now in 8.5.

### Phase 24 - Upside 3 audit: find and fix every bug (STARTED 2026-10-19, user)
User 10-19: "the fix only goes into upside 3. You can modify on current branch (martini-dev) as this will be the
'upside 3'. Later in the stage we should make it refuse to run ff2.1 (not yet ...). You can start the debug progress
now, on this branch." Code concise, simple and logical, in master's style (never example 00's). Method in
`upside3.md` 4.2; the defect ledger is findings 3.14.

**Resume here ("continue upside3 progress"; state at 2026-10-22 19:10, Mac with the scratchpad; item 2b to 10-27).**
Working files (fixtures, baselines, build trees, reviewer outputs) are under the gitignored `scratchpad/upside3_audit/`.
Item 2b's runs are finished (11-02; remote_jobs.md §1 line 3); the Peng seed replicates are scored (11-02); owed: glpG's HDX.
1. [x] `obj/` holds every fix to A231 (earlier binaries kept as `obj/*.bak_a<N>`, one per install from `*.bak_a101` to
   `*.bak_a230`, and `*.bak_20261021`). Baselines `w56_obj.npz` on `fixtures_w21` (the writers' current set:
   fixtures_w20 plus A231's `sequence_id`; older sets are refused in part: pre-A130 anchors, pre-A172 channel configs,
   pre-A222 pulling configs) and `h53_obj.npz` on `fixtures_hybrid_w6` (fixtures_hybrid_w4 less A227's 13 items plus
   A231's 4 datasets; h53 is h52 moved at float32 rounding by A226's fast-math code generation). The w4 hybrids were
   regenerated 10-22 through the full prep with A189-A192: KE/1.5kT 6.6 1.008, burn-in 1.000, production 1.007;
   peptide C-N mean 1.319, max 1.82 A; CA Rg 12.3-13.2 A; PO4 hold end 102150.
2. [x] Pass 2: every unit has both reviews (`pass2/<unit>/report.md`, all 18 P1 and 7 P2 units). Fixed from pass 2:
   placement and writers A108-A121, coordinate nodes A127-A132, engine graph A133-A138, dry-MARTINI engine (P2-01, -03,
   -04) A139-A142, driver/REMD/thermostat (P1-10) A143-A149, local backbone (P1-08) A150-A155 (D3 decided for the user,
   item 5.5), Monte Carlo (P1-11) A156-A158, hybrid coupling (P2-02) A159-A163, membrane and surface (P1-09) A164-A173
   (its M2, M3 left for the membrane retraining: findings A173), shipped parameters (P1-18) A174 (its D1 and D7 for
   the user, item 5.6), I/O and C API (P1-12) A175-A182, hybrid preparation (P2-05) A183-A188 (its D1, D6, D8 and the
   membrane-thickness definition decided, item 5.6), ConDiv trainer (P2-07) A193-A196 (D1 decided by the user: the
   SI's SARW), parameter glue (P1-14) A197-A199 and A204 (D1 and D6 decided for the user 10-22: torch replaced by
   closed-form numpy gradients, so a process holds one OpenMP runtime; the 86 dead `rotscalar` coordinates removed,
   which stops older checkpoints resuming in a newer tree, remote_jobs.md 0b), hybrid analysis (P2-06) A200-A203,
   build (P1-17) A205-A213, analysis (P1-15) A214-A220, examples (P1-16) A221-A225.
2a. [x] The closing criteria (upside3.md 4.3) on the final engine, met 10-22 (`scratchpad/upside3_audit/closing/`):
   * invariance: every test passes (equipartition on chig and on the hybrid production state, pivot MC, REMD,
     thermostat v / mv fast / mv slow; `closing/inv.summary`).
   * node_fd: all 30 configs; every parameter derivative exact or "not implemented" (A29-A31); every flagged
     direction, all in the rotamer node, line-integrates to noise (`closing/{fd,li,li2}.summary`, `gfix/li_75.log`),
     the pair-table ones once `tests/line_integral.py` follows node_fd's tied direction (A230).
   * sanitizers: an ASan + UBSan build, `tests/regress.py dump` on all 30 configs through the library and the binary:
     every config exit 0, no report, 4260 arrays (`closing/san.summary`; the preload goes to Python.app, whose
     `python3.11` launcher drops DYLD_*).
   * registry (`registry_20261022.md`, on fixtures_w20 and fixtures_hybrid_w5): 44 of 50 node types matched; the 6
     written-but-never-read rows hold only items read outside the engine's construction, each named in
     `pass2/registry/report.md`; its fixes A226-A229.
   ff2.1 stays usable (user 10-22: "keep ff2.1 useable for now"); refusing it waits.
2b. [ ] Validation of the debug on midway2 (user 10-22: "fix D9 with (a). submit Peng 16 validation protein with ff2.1
   and glpG hybrid simulation with ff2.1, on midway2, as a validation for our entire debug"):
   1. [x] D9 by design (a), A231: a required per-element `sequence_id`, the exclusion key of every node that
      excludes by sequence distance, written by every writer; baselines w56 and h53 (item 1); registry 44 of 50
      matched as before.
   2. [x] The Upside 3 tree `/beagle3/trsosnic/yinhan/upside3` with `env_upside3.sh` (CLAUDE.md "Shared Upside
      Deployment"). It is separate from the shared deployment because running campaigns read that tree (the
      released-ff2.1 Peng arms are this validation's control) and CLAUDE.md forbids deploying under them. Code
      rsynced from the Mac working tree, built on the midway2 login node, its own `.venv` with mdtraj.
   3. [ ] Peng: the 16 proteins x {native, denovo}, the benchmark's own protocol (`ff3_benchmark/bench_run.py`:
      14-replica REMD on Table S2's ladders, dt 0.009, frame 100), ff_2.1, in
      `/beagle3/trsosnic/yinhan/upside3_validation/peng/`; compared with the released-ff2.1 arms on the pre-audit
      binary (`ff3_benchmark/runs/<prot>_<kind>_ff_2.1/`) and Peng's FF2 (Table S2, Fig S4). All 32 scored (11-02,
      findings 3.15): native level with Upside 2 (+0.012 [-0.001, +0.028]), de novo lower, -0.030 [-0.049, -0.014],
      p 0.004, in 14 of 16, with one run per arm per engine. Seed replicates (top7 and alpha3d de novo, seeds 1, 101
      and 202 on both engines; `peng_seeds/score/seed_spread.py`, 11-02): top7 -0.087 +- 0.031 (SE over seeds), beyond
      its seed spread; alpha3d -0.086 +- 0.064, within it (findings 3.15). The gap's cause is unproven.
   4. [ ] glpG-RKRK-79HIS hybrid with ff2.1 through `example/16.MARTINI/run.py`, the repo's job script
      rewritten for the Upside 3 core (A232; user 10-22: "the job scripts should also be updated to support it"):
      the preparation as a job, then 4 REMD blocks of 48 replicas on T 0.70-0.90, each stopped at `--time-limit` and
      queueing the next. `run_sim_1afo_full.sh` sources the Upside 3 tree too (A233). The campaign's `run_remd.py`
      stays as it is under the running campaigns (Upside 2; it also carries a NaN-rollback guard). At its block ends,
      besides KE, C-N and TM helicity: TM4's protein protection at 298 K (`pp_fail`, dG under PS_protein) with HDX
      through the fixed analysis code (A61 fixed), the follow-up of Phase 16's closed question. Block 1 read 10-26
      (findings 3.15): KE, C-N and five of six TM helices as Upside 2's hybrid, but TM4 is lost (DSSP alpha 0.44 at
      T 0.70 against 0.93-0.96), already in preparation stage 7.0. The loss is ff2.1's own (Known Errors). Done 11-02:
      ff2.1 4 blocks (TM4 0.44 to 0.05), ff3.0 3 blocks (TM4 0.85, 0.82, 0.85; findings 3.15). The HDX readout waits
      on the user's choice of force field.
   Next after it: stage S2 (`upside3.md` section 3), then the open items below; refusing ff2.1 is S7, at release
   (`upside3.md` 9.3).
2c. [x] The protein harness rerun on HEAD (Mac Studio, 10-29; the hybrid fixtures are MacBook-only, not rerun):
   regress 0 of 4064 arrays against the 10-27 records; invariance remd, thermostat, mc and chig equipartition PASS;
   index oracle 0 mismatches on all 28 fixtures after A237 (the 10-27 flags were the oracle's own); node_fd's 6
   flagged directions line-integrate to noise, 1rkl_channel's once the surface score's steps are removed (findings
   3.14, after the ledger; `upside3.md` 11 item 15).
   Open from the build review (`pass2/build/report.md`): the cluster binaries set flush-to-zero and
   denormals-are-zero in every process that loads the library (U1, `set_fast_math` found by `nm` in both trees),
   measured at the first midway2 build (remote_jobs.md 0b); `sse2neon.h` is an unpinned Homebrew header whose
   defaults set the Mac's `_mm_rcp_ps` accuracy (U2). Its B9 decided for the user 10-22: `.devcontainer/` deleted
   (A213); a container for Upside 3 is written when there is a runtime to test it on.
   Open, performance only: the g-JF round evaluates the forces once more than it needs (MARTINI review U4, about
   20% of a hybrid's cost); skipping it needs to know that nothing moved the coordinates since the last evaluation
   (barostat, rigid fit, MC, exchange), so it waits until correctness work is done. U6, a hot rigid protein in 6.6,
   is a transient.
3. [x] The glpG protocol validation (findings 8b): the rigid-stage gate failed one healthy prep in three, and every
   prep's flexible production closes the belt ends, so the gate measures production's steps 15,000-20,000 (A126,
   user 10-21: option (a)); `run.py` starts REMD from the gated production's last frame. Rerun end to end through the
   workflow by the 2b.4 preparation (10-23, gate passed: 0.0% of 39 belt sites bare; remote_jobs.md §1 line 3) and the
   ff3.1 and ff3.0 preparations of 10-27; first tested on recorded productions (`glpg_gate_20261021/gate_a/test.log`).
4. [x] The open fixes of 10-21, all in the ledger and `obj/`: D5 as A103, the four decisions of item 5.1-5.4 as
   A104-A107, placement D1-D3 as A108-A110, writers items 1-5, 7-12 and 15 as A111-A121 (6 and 13 are A106-A107), the
   `HBondEnergy` input checks as A122, BP D4 as A123, the dead-code batch as A124, the reader preconditions as A125;
   D9 as A231 (item 2b.1).
5. Decided 10-21 (user: "make decisions for me, do whatever test you need to"; the evidence is in each ledger row):
   1. `membranelateral_potential` deleted (A105): an xy force on surface CBs with no energy. Its forces do not sum
      to zero, as any energy of the xy-uniform membrane must (1qhj net force up to 0.096 E_up/A,
      `scratchpad/upside3_audit/lateral/net_force.log`), so logged energy, pivot MC and replica exchange could not
      see it. A lateral-pressure term, if wanted, is a new energy, the integral of P(z) A(z) dz over the protein's
      cross-section (upside3.md 8.5).
   2. The rotamer-library generators deleted (A104): nothing builds or calls them, no Upside 3 plan rebuilds
      libraries from Dunbrack, and master keeps them as the libraries' provenance (noted in up.md).
   3. Chain breaks (A106, writers patch 6 as written): every shipped workflow (tests/fixtures.py, examples 03-16)
      converts with `--record-chain-breaks` just before configuring, so only a reused single-chain input directory
      is affected, and it now fails loudly; re-run the converter.
   4. Insertion codes are refused by `PDB_to_initial_structure.py` with a request to renumber (A107): every consumer
      of the PDB numbering is integer-only, and none of the project's 46 PDBs has one.
   5. `--secstr-bias` (local backbone pass 2 D3, A155; user 10-22: "make decision for me"): its windows were
      non-periodic in phi, stepping by the full bias across +-180. Decided: one phi window for both basins, phi < 0
      continued past +-180 to +120 (the empty gap between beta and alpha_L), the psi windows on the same form, each
      the original pair of 10 degree edge sigmoids summed over its periodic images. A cosine-form window was
      rejected: on a 120 degree gap with 10 degree edges it leaves 3.5% of the bias on alpha_L.
   6. Decided by the user 10-22 09:05 (all four as recommended), with the backbone bead types by DSSP through
      martinize's ssClassification (user 10-22 09:20: "Use DSSP"; A191; it changes helical proteins too, glpG 84 of
      210 residues):
      * the side-chain/environment table rebuilt from each rotamer state's ideal-geometry heavy atoms at its library
        (chi1, chi2) rows, chi3/chi4 sampled, mapped by martinize's martini22 centres (A192; confirmed against 1146
        crystal side chains: the old table put the beads on Upside's trained rotamer point and vector), averaged
        azimuthally at the simulation temperature (0.8647 T_up), not 25 T_up; `martini.h5` md5 7496dd90, old one
        `martini.h5.bak_pre_martinize_beads_20261022`;
      * dead content deleted (A189): `parameters/ff_2.0/`, `parameters/packing/`, `dryMARTINI/DDM.pdb`, and the
        unread martini.h5 datasets;
      * the ion count from the P-P thickness, the prep's fallback measuring it from the PO4 beads (A190).
      The hybrid fixtures were re-prepared end to end on them (`fixtures_hybrid_w4/`, baseline h52, item 1). The
      cluster's live hybrid files are re-prepared before Upside 3 runs them (remote_jobs.md 0b).
   Still with the user: `upside3.md` 11.
6. Git is the user's. The commit should include `tests/` and `example/15.SpatialTransformation/calc_rmsd.py`.

**Revised decisions** (CLAUDE.md Development Rules reworded to match, 10-19):
* martini-dev is Upside 3. Fixes go only here; master is not changed.
* ff2.1 keeps running on this branch for development. Refusing it comes later.
* The harness lives in `tests/` (new). Its generated fixtures and outputs go to the gitignored
  `scratchpad/upside3_audit/`.

**Harness and pass 1 (10-19 to 10-21).**
* `tests/` (method `upside3.md` 4.2): `fixtures.py` (ff_2.1 configs for chig, 1dfn, 1ubq and membrane 1qhj, plus one
  config per writer option group, 44 of 44 node writers), `node_fd.py` (each potential leaf's forces and parameter
  derivatives against fourth-order differences on the realized float32 steps, with tie maps), `regress.py` (node
  outputs, forces and parameter derivatives on 3 frames plus a seeded REMD-and-pivot run, compared bitwise; a
  hybrid's REMD runs without pivot moves, which the engine refuses with the PO4 hold's fixed atoms),
  `line_integral.py` (settles a flagged direction whose stencil crosses a kink), `invariance.py` (pivot MC detailed
  balance, REMD on a harmonic ladder, equipartition with the SE from the pooled statistical inefficiency g, since
  lipid kinetic energy stays correlated for about 55 frames under the diffusion-calibrated friction, findings
  10.33), `index_oracle.py` (a numpy reconstruction of every writer's index arrays; 69 of 69 injected errors
  caught), `registry.py` (each node's C++ reads against its Python writers). Sanitizer build
  `scratchpad/upside3_audit/build_asan`. `get_param_size` was added to the C API and `upside_engine.py`.
* Pass 1, one reviewer per P1 unit, each claim with a reproducer: every unit reported; ledger A1-A101 fixed and
  verified (A2, a heap overflow, had made the first baseline irreproducible), plus the Python option batch (13
  options repaired). The user's 10-19 approvals are applied: exact derivatives or "not implemented" for every
  parameter, dead nodes and options deleted, and the angular cosine held to [-1,1]. Not changed on purpose: the
  H-bond 3.5 A cutoff (findings 3.14).
* Hybrid preparation (A70-A73, A102, A126): every prep MD stage at dry-MARTINI's own Langevin friction (mass / 4 ps;
  the overdamped overlap_settling friction deleted); per-lipid diffusion targets (user 10-20) in
  `LIPID_LATERAL_DIFFUSION_UM2_S`: DOPC 11.5, POPE 8.6, POPG 11.5 um^2/s, from ratios to DOPC (findings 2.5; no
  same-method set exists); ions keep the DOPC clock for now.
* The starting membrane (findings 8, 8b). GROMACS running CHARMM-GUI's own protocol on the same step5 bilayers
  expands the box as Upside's prep does, then recovers over about 12 ns as the overlapped leaflets separate; the
  prep covers about 5 ps and only sees the jump, so there is no engine defect. Decision: the prep starts from a
  GROMACS-equilibrated membrane of the same potential (no geometric correction; the leaflet-translation patch is
  dropped). Started at the GROMACS equilibrium (`scratchpad/bilayer_ref/upside_from_gmx_eq/`), Upside holds DOPC;
  POPE/POPG ends 1.3% smaller and 0.4 A thicker, still drifting at 0.73 ns. Open:
  [ ] a longer POPE/POPG run to settle the 1.3%; [ ] equilibrated membranes as the prep input.

### Phase 25 - are the finished force fields done training; the side-chain verdict (DONE 2026-10-20, user)
User 10-20: "check if the finished ffs needs no more trainings (converge, etc)"; "check if the deterministic /
probabilistic side chain have a conclusion". Read only; working files in
`/project/trsosnic/yinhan/checks/converge_20261020/`. The gate on `ff31_depth_data` and `ff31_depth_all` over steps
57-75: both converged (lowest dhb 0.033, 0.032); their depths show no pull (`depth_pull.py`). The defect left no
detectable pull: the control's first fixed-engine epoch and the ff3.2 rerun's epoch 0 show no angular pull, and
`mdc_ctrl_03` equals `mdc_ctrl38` (-0.012, p 0.67).

**Review (findings 1.21, 1.36).** Released and finished force fields need no further training by the gate's rule:
ff_3.0_gdepth, ff_3.1_depth_data, ff31_depth_all and the control; the defect's direct test for ff3.1 is its rerun
(step 76 on 10-21). Two exceptions, both reported to the user: ff_3.0_gdepth's glycine alpha_L depth drifts one way
(+0.015 nats a round), and no-BP's ~110 angular coefficients are not at a fixed point. The side-chain verdict is in:
at both fixed points the deterministic side chain folds 0.241 against 0.353 on 89 held-out domains. Choosing the
side-chain model for Upside 3 (upside3.md gate) is the user's.

### Phase 26 - the side-chain comparison on one metric, TM-score (DONE 2026-10-21, user)
User 10-21: "I prefer to use a single metric (TM score?) for all tests". TM-score as the Peng benchmark is scored
(`ff3_benchmark/scoring/tmscore.py`: fixed correspondence, superposition optimised, step 2, normalised by the full
chain): per domain, 4 seeds pooled, frames after T_EQUIL, mean TM and the fraction of frames at TM > 0.5. mdCATH at
the fixed points (`mdc_nobp_03`, `mdc_ctrl_03`) by `checks/converge_20261020/tm_panel.py`; the 44-domain panel rerun
at both fixed points at the training time step 0.009 (tags `p44_nobp_03`, `p44_ctrl_03` in `ff3_selection/runs/`,
params `params/mdc_<nobp|ctrl>_03`, `checks/converge_20261020/panel44_dt009.sbatch`); Peng's benchmark for the
38-step releases only (user: "44-domain panel only").

**Review (findings 1.21).** Mean TM to the native: mdCATH -0.034 (p 7e-8), the 44 domains -0.026 (p 0.007), Peng
natives -0.029 (p 0.18), de novo +0.010 (p 0.49); TM > 0.5 shares equal within error.

### Phase 27 - released ff2.1 on the same three tests (STARTED 2026-10-21, user)
User 10-21: "run validation with released ff2.1 so we know if we should keep ff2.1 or not". Both side-chain arms are
ConDiv retrainings of ff2.1 (`ff_2.1_ctrl38` shares none of ff2.1's trained files, only its Ramachandran library), and
both sit below Peng's published FF2 on his benchmark (findings 1.21), so neither answers whether ff2.1 itself should
stay. Decisions:
- **The released files** are `ff3_selection/params/ff21_released` (md5 equal to `/beagle3` `parameters/ff_2.1` and
  `common/rama.dat`); every run uses each test's own protocol, tool and binary, unchanged, so ff2.1 pairs with the
  retrained arms on the same domains, seeds and time step.
- **Time step 0.009 everywhere**, as the retrained arms' panels and Peng's protocol ran (the 10-01 44-domain
  `ff21_released` panel ran at 0.015 and kept no trajectories, so it was rerun).
- **Peng benchmark:** the 30 arms not yet run, with `bench.sbatch` as `validate_peng.sh` submits them; the two lambda
  arms of 10-01 to 10-04 (11d) are the same protocol and binary (`/beagle3` `obj/` unchanged since 09-30).
- **Analysis runs by itself** (no agent polls): the panels' TM jobs wait on their panels by dependency; a 1-CPU
  collector job wakes every 6 h, scores each finished arm with `score_arms_ff.py`, and when all 32 are scored writes
  the three-way table and stops. Files in `/project/trsosnic/yinhan/checks/ff21_released_20261021/`.
1. [x] mdCATH panel `mdc_ff21_released` (89 domains x 4 seeds) and the 44-domain panel `p44_ff21_released`.
2. [x] TM against the fixed points (`tm_panel.py`): on the 44 domains released ff2.1 0.681, the control 0.650, no-BP
   0.623 (findings 1.37); on mdCATH 0.727, 0.702, 0.667 (remote_jobs.md §1 line 1).
3. [ ] Peng's 30 remaining arms (scoring in progress, remote_jobs.md §1 line 1); compared against the step-37
   `ff_2.1_ctrl38` and `ff_2.1_nobp` and against Peng's published FF2 (0.55 native, 0.42 de novo).
4. [x] User 10-21: does the loss reach the other trained force fields? 44-domain panel at 0.015, every epoch end
   against `ff21_released` (`panel44_trajectory.py`): every run of this workflow folds 0.13-0.23 less, most of it in
   epoch 0 (findings 1.37).
5. [x] The same table with the owed epoch-end panels (`panel44_trajectory_final.txt`) and Peng's arms of
   ff_3.0_gdepth and ff_3.1_depth_data scored (`peng_trained.txt`): findings 1.37.
6. [ ] Review: keep ff2.1, or a retrained force field, by the three tests.

### Phase 28 - is retraining's folding loss a workflow defect (STARTED 2026-10-21, user)
User 10-21: "is this a defect in our training workflow? can you fix it?" Every run of this workflow loses 0.13-0.23
folded in epoch 0 (findings 1.37). Decisions: attribute before fixing, by the round-3 method on the round-4 control's
epoch 0 (`c9_00`, ff2.1's library, so every difference from ff2.1 is training), one file at a time put back to
ff2.1's; no workflow change without the user, since every candidate fix (rates, a prior toward ff2.1, freezing
groups, early stopping on a validation panel) departs from the FF2 workflow the project fixes (findings 10.17, memory
`condiv-must-be-ff21-workflow`), while a port defect against the SI is fixed directly. Found (findings 1.37;
`checks/ff21_released_20261021/resets_c9_00.txt`):
- No single file restores ff2.1's folding: the side-chain table (+0.089), the backbone burial (+0.081) and the H-bond
  (+0.068) each give back part and together about all of the 0.225 lost; environment (+0.028) and sheet (+0.035)
  are not resolved. hb and env are pulled one way from step 0, rot spreads as noise, and the training proteins'
  RMSD is flat over four epochs.
- No port defect: the ensemble definitions are Kleinmann's (Peng's code). The loss holds at the training time step
  0.009 (released ff2.1 0.603 against the control's fixed point 0.430). The SI-rate twin `ff30_gdepth_si` lost as
  much (`ds_02` 0.406), so the rates set speed, not the loss.
- Peng and Jumper stopped on a fixed budget short of convergence on purpose, Jumper having seen accuracy fall when
  training continued (findings 1.23). Our gate is our own addition, and every run started from ff2.1, so each is a
  second budget on top of FF2's.
- Why ff3.1 and ff3.2 end at their starting maps (user 10-21): each depth sits where its 0.2-nat prior, centred on
  the start, balances the data; the glycine pull is consistent only pooled over the 37 maps. From ff2.1 the data
  pull glycine alpha_L shallower (t +5.6), from ff3.0's map not at all (findings 1.34). A retraining that should move
  the maps needs a pooled prior, which 1.29 found to transfer best (Phase 29).

The fix options (item 5) were answered by the user's Phase 29 revised decision: every Adam group held at ff2.1.

### Phase 29 - ff3.1.1 and ff3.2.1: the depth training with a pooled prior (STARTED 2026-10-21, user)
User 10-21: "apply the fix and train ff3.1.1 and ff3.2.1", after Phase 28 found each depth held at its start by its
own 0.2-nat prior while the data's pull is clear only pooled. Decisions:
- **The model** (`training/ff31/rama_basin.py`, mode `pooled`): a map's offset on a basin is its central residue's
  shared offset plus the map's own deviation, c_mb = q_X(m),b + d_mb. The deviations are `data` mode's, at its width
  0.2 and on its 1,442 depths. The shared offsets pool the maps of one central residue only (rule 1.8: central GLY
  says nothing about central ALA; the coupling is estimation, as NDRD's hierarchical prior), each under its own
  Gaussian prior, trainable where the residue's pooled native sites reach `data`'s threshold of 20. Pre-proline maps
  (X|right|PRO) keep only their own deviations, as ff3.0's glycine pair left them out (findings 1.14). Same damped
  Newton step per epoch, solved jointly over shared offsets and deviations.
- **The shared offsets' width is cross-validated**, as the per-map widths were (findings 1.29): 10 folds over the
  proteins of each source run's epoch 0, the trainer's own step and held-out map-level gap, widths 0 (= `data`) to 2.
  Chosen: 0.5 (findings 1.34).
- **Verification before launch, no test jobs** (findings 10.16): with the shared offsets off, mode `pooled` gives
  `data`'s step exactly on recorded epoch-0 data (1e-17), and P and its transpose are adjoint (4e-16).

**Revised decision (10-21 13:40, user: "besides rama map not moving bug, try to resolve the accuracy lost bug as well
in this run").** The first chains (`ff311_depth_pooled`, `ff321_depth_pooled`, every group trained) were cancelled 3
minutes into step 0 and kept as `<run>.cancelled_20261021`. Reasons, from Phase 28: the loss is spread over the
co-trained side-chain table, backbone burial and H-bond (findings 1.37); holding some groups while the rest train
lost more folding before (the 10-07 frozen arms, findings 1.21); and Peng and Jumper stopped those terms on a fixed
budget that ff2.1 has already had (findings 1.23). So:
- **Every Adam group is held at ff2.1** (rate 0, `checks/ff311_20261021/freeze_groups.py` on the initial checkpoint;
  the trainer is unchanged). Only the depths train, under the pooled prior: ff3.1.1 = ff2.1 plus depths from NDRD's
  map, ff3.2.1 = ff2.1 plus depths from ff3.0's map. ff3.0's other files carry the loss (folded 0.464), so ff3.2.1
  takes ff2.1's; the two runs differ only in the starting map, and share ff3.1's rerun's minibatch order and
  held-out set (same `ConDiv.py` md5 `df3e2497...`, inputs and fixed engine; target 76, no gate).
- **Runs** `ff311_depths_only` and `ff321_depths_only` (`build_depths_only.sh`).

At step 76 (10-23; findings 1.34): ff3.1.1's glycine alpha_L +0.41 from ff2.1, ff3.2.1 at ff3.0's map; every map but
glycine's alpha_L agrees between the runs to their noise, and that gap went 0.517 to 0.179, still closing ~0.07 an
epoch. The per-run panels were superseded by Phase 30 (user 10-23): the two runs train on until they agree, and the
merged force field is what gets validated.

### Phase 30 - merge ff3.1.1 and ff3.2.1 into ff3.1 and validate it (STARTED 2026-10-23, user)
User 10-23: "train 3.1.1 and 3.2.1 further until they agree, then we can merge them into one (overwrite 3.1) and run
validation on it"; chosen: glpG hybrid included in the validation, merge by averaging. Decisions:
- **Agreement rule, fixed before the data:** at one round, (a) glycine's alpha_L shared-offset gap is no larger than
  the rms of the other trainable shared-offset gaps (the runs' noise, 0.092 at round 4), and (b) the glycine maps'
  gap rms over the data-reached depths is no larger than the non-glycine maps'. Train in epoch pairs: 114, then 133;
  at 152 stop and report.
- **Merge:** the mean of the two runs' offsets per (map, basin) in the common frame (NDRD plus offsets, ff3.2.1's
  including ff3.0's pair), written by `rama_basin.write_library`; the other five files are ff2.1's.
- **Validation:** released ff2.1's three soluble tests on the same binaries and protocols (44-domain and mdCATH
  panels at dt 0.009 with TM comparisons against `ff21_released`, Peng's 32 arms) and glpG on Upside 3 (TM4 at each
  block end, HDX on block 2).
1. [x] p311/p321 dropped from `submit_new.sh` (backup `.bak_pre_phase30_20261023`); both runs continued to step 114.
2. [x] At 114 (10-26): KE/1.5kT at most 1.019, 24 of 24 proteins at every step. Round 6: (a) met (0.072 against
   0.084), (b) not (0.092 against 0.084), short only by the constant +0.063 shift of ff3.0's fitted pair, which the
   data cannot move (findings 1.34). The user merged on that (10-26: "merge into ff3.1 now").
3. [x] Merged and released 10-26 as `ff_3.1` (`scratchpad/ff31_merge/`; `checks/ff31_merged_20261026/release.sh`):
   rama.dat `5b66c78a`, the other five ff_2.1's; glycine figure to ~/Downloads; the old ff3.1 and ff3.2 sets and
   their results moved aside as `.bak_pre_merge_20261026` (user: "from now on ff3.1 refers to the merged version").
   Not through `validate_ff.sh`: it extracts a checkpoint and patches Upside 2 glpG seeds, and the glpG validation
   runs on Upside 3 instead (user's choice).
4. [ ] Validation (`checks/ff31_merged_20261026/submit.sh`: panels mdc_ff31 and p44_ff31 with their TM jobs, Peng's 32
   arms and the collector; glpG on Upside 3, Phase 24 2b.4); read each as it finishes (remote_jobs.md §1 lines 2-3).
   Panels read 10-27: ff3.1 is level with released ff2.1 on both (mean TM +0.006 [-0.001, +0.012] on the 44 domains,
   +0.005 [-0.001, +0.012] on mdCATH's 89) and above the control on both (+0.037, +0.031). glpG on Upside 3: ff3.1
   loses TM4 in preparation as ff2.1 does, ff3.0 holds it (findings 3.15; Known Errors). Peng's arms are still running.
   Seven of them (proteinB, proteinL and top7 both kinds, proteinG de novo) stopped 10:31 10-28 when their chunks'
   self-resubmission met the Slurm controller outage; resubmitted 10-28 on the user's word with `bench.sbatch`'s own
   resubmission line, each resuming from its checkpoint (remote_jobs.md §1 line 2).

### Phase 31 - glpG in a DDM micelle with ff3.1 on Upside 3 (STARTED 2026-10-26, user)
User 10-26: "setup a DMM simulation for glpG, which we did experiment on ... DMM should be single layer instead of
bilayer ... run it with ff3.1"; "I will build a DDM micelle with charmm-Gui (without protein in it), but you will
provide me the radius and number of DDM molecules"; "you should implement the missing dihedral that is needed".
Decisions:
- **The single layer is a micelle.** A planar DDM layer cannot cover glpG's 28.2 A belt: one C12 tail under a
  maltose head gives a lamellar core of 11-14 A (findings 7.1, readme section 5), and the old DDM slab unfolded TM4.
- **DDM is the official Martini 2 topology**, cgmartini `martini_v2.0_DDM_01.itp` (Bjorling and Westenhoff; md5
  `0d587b39`): its three proper dihedrals (GROMACS type 1) are implemented in Upside 3, and its five head
  constraints are carried as the file's own `-DFLEXIBLE` harmonic bonds, Upside having no constraint solver.
- **Sizes from measurement.** Empty DDM micelle (Oliver et al. 2013, SAXS): aggregation number 135-155, oblate core
  semi-axes 15.7-16.1 and 27.9-28.3 A, head shell 5.4-5.8 A, monomer 698 A^3. glpG's belt (oriented CG, |z| <= 14.1
  A) has effective radius 15-18.5 A and reaches 26 A, so the empty micelle's core cannot hold it; a 16 A tail annulus
  around the belt holds 150 (rounded rim) to 200 (square rim) DDM. Told to the user: build N 140, head radius ~26 A.
1. [x] **The dihedral term** `dihedral_spring` (findings 3.14 A236, 10-26), beside `dist_spring` and `angle_spring`:
   V = k (1 + cos(n phi - phi_s)), phi GROMACS's, from one canonical spline of 1 + cos(u); the preparation writes it
   from the ITP's `[ dihedrals ]`, the ITPs read under `FLEXIBLE`; the DDM ITP is the official file. Equal to
   GROMACS's energies and forces on 600 DDM dihedrals, `node_fd` clean, `regress` bitwise unchanged, `index_oracle`
   clean. Built locally only (`obj_dihedral/`); the cluster's Upside 3 tree takes it at step 4, built on midway2.
   Still owed with step 3: DDM in the hybrid's own integration (KE/1.5kT, bond lengths, dihedral distributions
   against GROMACS on the same complex).
2. [ ] **The glpG-DDM complex, built locally by self-assembly** (user 10-26: "charmm-gui fails to build micelle this
   size. Let's build it locally"). The empty micelle is skipped: its core cannot hold glpG. GROMACS 2026.4 on the
   Mac, the potential Upside reproduces (findings 8: reaction field eps_r 15, potential-shifted LJ to 1.2 nm, `sd`
   at 303.15 K), NVT: the prep's oriented CG glpG held fixed, 200 official-topology DDM started in a shell around its
   belt, run until the number bound and the belt coverage stop changing. Free detergent is then removed and the
   complex is the environment the preparation takes. `scratchpad/ddm_micelle/`.
   The 200 ns run (`results.txt`, `contacts_profile.py`; findings 7.5), steady over 100-200 ns: all 200 DDM bound,
   no free micelle, tails within 0.94 nm of 95% of the belt's exposed beads, TM4 tail-covered (W146 tail 1.00). But
   maltose heads lie against the belt: a polar maltose bead within 0.6 nm of 56% of the exposed beads at |z| < 0.5 nm
   (tails 62%), and A37, V41, V122, A190, V194, A197, M198 and V201 are touched by heads only, where POPE/POPG around
   the same glpG puts no head at the belt core. Not a clean detergent ring, so not yet used. Dry Martini gives DDM
   alone the SAXS micelle (findings 7.5: 140 DDM form a micelle of 124, tail core semi-axes 14.5 / 25.3 / 29.3 A
   against 15.7-16.1 / 27.9-28.3, heads 99.5% outside), so the coat around glpG is the protein-detergent contact of
   GROMACS's Martini glpG or the excess detergent; the side-chain/lipid table covers DDM's bead types, so the
   arrangement to judge is the one Upside relaxes to (step 3).
3. [ ] **The preparation's micelle path.** It is bilayer-only today: `require_bilayer_lipids` refuses a single-tail
   detergent, and the xy barostat, tiling and ion slab assume a bilayer. A micelle runs NVT; DDM needs its own
   diffusion target for the friction calibration, from a measurement; the belt solvation gate applies as it is.
4. [ ] **glpG-DDM with ff3.1 on Upside 3**: preparation and REMD, read as the POPE/POPG run is (KE, C-N, DSSP per TM
   helix, TM4).

### Phase 32 - debug the ff3.1 training: does it recover ff3.0's glycine pair (STARTED 2026-10-27, user)
User 10-27: "I suspect there is a defect in our training ... assuming we did everything correct it should converge to
a correct rama map (ff3.0) for GLY", then "setup your debug for ff3.1 training". Findings 1.34 (10-27): the trainer
held ff3.0's glycine pair for four rounds in ff3.2.1; ff3.1 leaves ff3.0 through glycine's beta' basin, which ff3.0
never trained and which pulls in every epoch, and through rounds 5-6, where the two runs moved together because they
share the minibatch order. Decisions:
- **The test with a known answer.** ff3.1.1 as built (`build_depths_only.sh`: ff2.1, every Adam group held, depths
  in mode `pooled` from NDRD) with one change, glycine parameterised as ff3.0 is: only its shared alpha_R and alpha_L
  train; its four other shared offsets and every GLY map's deviations stay 0 (a mask on the initial checkpoint,
  `restrict_gly.py`; the trainer is unchanged). Every other map trains as in ff3.1.1. Expected if the training is
  correct: glycine's shared pair reaches ff3.0's (alpha_R -0.107, alpha_L +0.517 in the common frame), the pair
  ff3.2.1's data held at no pull for four rounds.
- **Two runs, two minibatch orders**, so their spread is the training's own noise: `ff311_glypair` (ff3.1.1's order,
  so it differs from ff3.1.1 only in the glycine mask) and `ff311_glypair_mb2` (the order `initialize` drew; the same
  456 proteins and held-out set, which a fixed seed draws from the sorted list).
- **Reading, fixed before the data:** per round, both runs' glycine pair, its per-epoch data pull (t over 19 steps,
  `pooled_pull.py`) and -ln(P_alpha_L / P_alpha_R) over the 37 GLY|X maps (`gly_LR.py`, ff3.0 -0.598). Passed when the
  pull on both is under |t| 2 for two epochs and the pair, averaged over both runs' last two rounds, is within 0.1 of
  ff3.0's on each basin (0.1: the other residues' shared-offset gap between ff3.1.1 and ff3.2.1, rms 0.084-0.116).
  Passed means the trainer is sound and ff3.1's glycine is the beta' basin and the co-training; failed means a defect
  in the glycine pair's path, then located.
- Target 76 (four rounds), stop for the user. Scripts and logs `/project/trsosnic/yinhan/checks/ff31_debug_20261027/`.
1. [x] Built 22:50 10-27 (`build.sh`, `build.log`): 1,517 trainable depths -> 1,403 (glycine's shared offsets 6 -> 2,
   its 40 maps' 110 deviations -> 0); both initial checkpoints equal ff3.1.1's in every other field, native site
   counts to 4e-12 (summation order), the held-out 46 the same, `_mb2`'s order differing in 452 of 456 places
   (`compare_init.txt`).
2. [x] Submitted 23:03 10-27 to 76: 49634971 (R from 23:05), 49634972 (R from ~23:30); remote_jobs.md
   §1 line 2.
3. [x] Read at each round end and at 76. Round 1 (read 22:10 10-28, `ff311_glypair` / `_mb2`): glycine's pair aR
   +0.021 / +0.009, aL +0.216 / +0.163 (ff3.0 -0.107 / +0.517); epoch-0 pull aL t +4.1 / +4.2, aR -1.7 / -1.5;
   -ln(P_aL / P_aR) -1.025 / -1.065 (ff2.1 -1.219, ff3.0 -0.598). alpha_L moves toward ff3.0's and is still pulled;
   the test reads rounds 3-4. Round 2 (read 09:30 10-29), `_mb2` only: aR -0.006, aL +0.301; epoch-1 pull aL t +3.0
   (summed +30.8), aR -1.5, other shared offsets |t| at most 2.1; -ln(P_aL / P_aR) -0.913. alpha_L keeps moving
   toward ff3.0's and is still pulled. `ff311_glypair` hung in step 34 from 22:41 10-28 (a worker's Pool child died
   of a glibc double free, cause not identified); cancelled 09:20 10-29 with 3g0m's files kept, its insurance link
   reruns step 34 (remote_jobs.md §1 line 2), so its round 2 comes later. Round 3 (read 12:40 10-29), `_mb2` only
   (`rama_round_03.dat` md5 `f65a593e`): aR +0.012, aL +0.370; epoch-2 pull aL t +1.6, aR -0.1, the first epoch under
   |t| 2, other shared offsets |t| median 0.6, max 2.9; -ln(P_aL / P_aR) -0.863. The pair is still 0.12 and 0.15
   from ff3.0's, so the rule waits on round 4 and on `ff311_glypair`. Margin constant at +0.192, 24 of 24 proteins
   every step, KE/1.5kT at most 1.013.
   **Verdict at 76, read 16:50 11-04 (Mac Studio; `pooled_pull.py <run> 2 3` and round 4 as recorded): failed.**
   Rounds 3 and 4, `ff311_glypair` / `_mb2`: aR +0.051, +0.069 / +0.012, +0.012; aL +0.375, +0.456 / +0.370, +0.388.
   Averaged over both runs, aR +0.036 and aL +0.397, 0.143 and 0.120 from ff3.0's (rule: within 0.1 on each). Pull:
   `_mb2` under |t| 2 in epochs 2 and 3 (aL +1.6, +0.6), `ff311_glypair` not (aL +1.5, then +2.7). On the old
   native target the trainer does not reach ff3.0's pair in four rounds: alpha_L still rises and is pulled, alpha_R
   sits at +0.01 to +0.07 with no pull (|t| at most 0.3), so the data does not separate it from ff3.0's -0.107. The
   cause is not located. Phase 33's native-target defect is a candidate, unproven: this test never ran on the fixed
   target. Not rerun (decided 11-04 at the user's delegation, "decide for me what needs to be run"): the target is
   retired, and on the fixed target ff3.2's shared glycine pair after round 3 is aR -0.095, aL +0.580, within 0.1 of
   ff3.0's (Phase 38 item 11), though with every term trained, which this test excludes. The same test on the fixed
   target is proposed only if ff3.2's glycine moves away from ff3.0's pair.

### Phase 33 - fix the training's native target; ff3.1 retrained on it, and ff3.2 with every term trained (STARTED 2026-10-28, user)
User 10-28: "try a version with everything trainable, start with current ff3.1, call it ff3.2. Fix ff3.1 according to
your plan, still call it ff3.1 ... start from either ff3.1.1 and ff3.1.2 or current ff3.1"; then: ff3.2 on the fixed
target, the Phase 32 runs kept.
- **The defect** (findings 1.34, 10-28): ConDiv's native side is a replica under the force field itself with a weak
  distance restraint (1/3^2, holding ~1 A RMSD, the original trainer's), which holds the fold but not the backbone
  dihedrals. At ff2.1 it keeps each residue in its native basin 0.98 (alpha_R) down to 0.49 (glycine beta'), 0.40
  (non-glycine pPII') and 0.004 (non-glycine beta'); where it leaves the native, the native side follows the map as
  the free side does, so a depth's gap need not close and its offset drifts (glycine beta', every epoch of both ff3.1
  runs, the same gap in the held-out proteins).
- **The fix:** the native replica alone also holds each residue's phi and psi near the native's: a periodic harmonic
  spring per dihedral (the engine's `Dihedral` and `Spring` nodes, as Upside's omega spring is written), width 20 deg
  at T0 (k = T0 / sigma^2), a folded residue's typical phi/psi spread, so it holds the basin and leaves the basin's own
  breadth. The distance restraint stays. The gradient is unchanged: every <dV/da> is still evaluated under the free
  Hamiltonian, so the springs change only which conformations the native side samples. Implemented in the runs' own
  `ConDiv.py` (`add_native_rama_restraint`); no other change to the workflow.
- **Gate at step 0, fixed before the data:** in the step's native replicas, the own-basin population of glycines
  native in beta' rises from 0.49 to at least 0.75 and alpha_R natives stay at 0.95 or more; the free ensemble's is
  unchanged by construction. Failed means the springs do not hold the native, and both chains stop for the user.
- **Start: the released ff3.1** (rama.dat `5b66c78a`; its merged shared offsets and deviations from
  `scratchpad/ff31_merge/merge_ff31.npz`), every other file ff2.1's; base library NDRD with priors centred on it, as
  ff3.1.1's; ff3.1.1's minibatch order and held-out set. The round-0 library must equal `parameters/ff_3.1/rama.dat`.
  - `ff31_nativefix` (ff3.1): ff3.1.1's design, every Adam group held at ff2.1, every depth trainable (all six basins
    of every central residue that reaches the data threshold), nothing else held.
  - `ff32_nativefix` (ff3.2): the same, every Adam group at the trainer's rates (everything the engine gives a
    gradient for), as the ff3.1/ff3.2 reruns of 10-18 trained.
- Target 76 (four rounds), stop for the user. Read per round: glycine's six shared offsets and their pull,
  -ln(P_aL / P_aR) on the GLY|X maps (ff3.0 -0.598, ff3.1 -0.825), native-basin fidelity by class, and for ff3.2 the
  H-bond margin and the epoch-end panels as the earlier all-term runs were read.
  Scripts and logs `/project/trsosnic/yinhan/checks/ff31_fix_20261028/`.
1. [x] `ConDiv.py` (md5 `54383004`) with `add_native_rama_restraint`; static check on 1bm8 and 3v0r
   (`check_restraint.py`): the engine's dihedral equals the spring centres to 2e-6 rad, the springs cost 1e-10 E_up at
   the native and 0.5 k d^2 to 3e-7 off it, the force finite.
2. [x] Built 02:20 10-28 (`build.sh`, `build.log`): both round-0 libraries md5 `5b66c78a`, the released ff_3.1's; the
   merge's few untrainable deviations (at most 1e-7, the basin fit's rounding of ff3.0's pair) kept as fixed parts of
   the map; rates zero (`ff31_nativefix`) and the trainer's (`ff32_nativefix`); every other field ff3.1.1's
   (`compare_start.py`). The gate script on ff3.1.1's own step 0 reads glycine beta' 0.356: FAIL, as the old target.
3. [x] Submitted 02:26 10-28 to 76: 49635325, 49635326 (remote_jobs.md §1 line 2).
4. [x] **The step-0 gate, read 10:50 10-28: FAIL as written, and the threshold is out of reach** (findings 1.34).
   Glycine beta' natives (n 21) 0.656 in ff3.2, 0.660 in ff3.1 (old target on the same proteins 0.356); alpha_R
   0.998. An ideal 20 deg Gaussian about each native gives these glycines 0.726 at most, since half sit within a
   spring width of a basin edge, so no 20 deg spring can reach 0.75. Over epoch 0 so far (ff3.2 steps 0-10, n 211)
   0.754 against that ceiling's 0.769. Both chains died 10:41 10-28 in a Slurm controller outage (ff3.2 in step 11,
   ff3.1 in step 16; no code error). User 10-28 "proceed with your best decision": both continue, since the springs
   hold each basin at a 20 deg spring's ceiling and the 0.75 was set above it; the threshold is not changed, the
   decision rests on the measured ceiling. Links 49635669 (ff3.2) and 49635327 (ff3.1, nice 10) released 10:58 10-28,
   ff3.2 first; the debug chains (Phase 32) at nice 10 behind ff3.2.
5. **ff3.2 has priority over ff3.1** (user 10-28, after Jumper's and Peng's statements that ConDiv trains all
   interaction parameters jointly: Jumper thesis p. 42-43, 65-66; Peng 2022 p. 551, SI p. 3). On every scheduling or
   reading conflict ff3.2 goes first.
6. [ ] Read at each round end. **Round 1, ff3.2** (read 22:10 10-28; ff3.1 is at step 16, its link pending): glycine's
   shared offsets aR +0.026 -> -0.131, aL +0.422 -> +0.517, beta' -0.253 -> -0.316 (beta, pPII, pPII' within 0.04),
   so its pair is ff3.0's (-0.107 / +0.517) after one round; epoch-0 pull aR t -3.2, aL +3.3, beta' -2.6, while the
   held-out proteins pull aL the other way (-3.9 summed, training +40.2); -ln(P_aL / P_aR) -0.567 (ff3.0 -0.598,
   ff3.1 -0.825). Native fidelity over epoch 0 (`native_fidelity.py`): glycine beta' natives 0.759 (n 366; ff3.1's
   steps 0-15 0.756), alpha_R 0.997, non-glycine pPII' 0.904, non-glycine beta' 0.096 (n 26). **H-bond margin +0.192
   -> +0.106 (step 18) -> +0.059 (step 31), under the +0.10 line**, E_beta now below E_alpha; at the same steps the
   all-term ff3.1 rerun `a31` read +0.160 / +0.119 and the control `c9` +0.135 / +0.130, so this run loses the margin
   about twice as fast. 24 of 24 proteins every step, KE/1.5kT at most 1.017. User 10-28 "proceed with your
   suggestion": the run keeps training under the +0.10 line, as the earlier all-term runs did.
   Glycine at 22:40 10-28: the pooled glycine map's alpha_L - alpha_R (`gly_rama_ff31_maps.py`'s cells) is -0.812
   in ff3.2's round 1, against ff3.0 -0.836, the merged ff3.1 -1.064 and ff2.1 -1.460; over epoch 1's first 14
   steps the pull is gone (alpha_R t -0.7, alpha_L +0.3, beta' -1.1), so the data holds it at ff3.0's balance.
   ff3.1 (steps 0-15, no round written yet) is pulled the same way: alpha_L t +3.3, alpha_R -1.4, beta' -3.0.
   **Round 2, ff3.2** (read 09:30 10-29): glycine's shared offsets aR -0.131 -> -0.196, aL +0.517 -> +0.536, beta'
   -0.316 -> -0.350 (beta +0.057, pPII +0.056, pPII' -0.052); the epoch-1 pull is under |t| 2 on every glycine basin
   (aR -1.3, aL +1.1, beta' -1.5; held-out aL -3.6 summed), the other shared offsets' |t| median 0.8, max 3.1.
   -ln(P_aL / P_aR) -0.567 -> -0.480, now 0.12 above ff3.0's -0.598 (less alpha_L than ff3.0). Native fidelity over
   epoch 1: glycine beta' natives 0.773 (n 366), alpha_R 0.997, non-glycine pPII' 0.906, non-glycine beta' 0.098
   (n 26). H-bond margin +0.080 at epoch 1's end, +0.036 at step 50, +0.051 at step 54 (dhb -0.489 -> -0.541); 24 of
   24 proteins every step, KE/1.5kT at most 1.016.
   **Round 3, ff3.2** (read 12:40 10-29; `rama_round_03.dat` md5 `67dd261c`): aR -0.196, aL +0.536 -> +0.567, beta'
   -0.350 -> -0.385 (beta +0.047, pPII +0.049, pPII' -0.056); the epoch-2 pull aR t -0.2, aL +1.3, beta' -1.3, the
   other shared offsets' |t| median 0.6, max 2.0. -ln(P_aL / P_aR) -0.451, a further 0.03 from ff3.0's -0.598 (less
   alpha_L still). Native fidelity over epoch 2: glycine beta' 0.777, alpha_R 0.989; non-glycine alpha_R 0.997, pPII'
   0.907, beta' 0.098 (n 26). The margin came back to +0.058 at epoch 2's end and is +0.061 to +0.067 over steps
   57-61; 24 of 24 proteins every step, KE/1.5kT at most 1.015.
   **Round 4, ff3.2, at the target** (read 17:10 10-30; step 76 done 11:58 10-30; `rama_round_04.dat` md5
   `96e265d2`): aR -0.196 -> -0.231, aL +0.567 -> +0.554, beta' -0.385 (beta +0.066, pPII +0.062, pPII' -0.053); the
   epoch-3 pull aR t -0.6, aL -0.0, beta' -0.4, the other shared offsets' |t| median 0.5, max 2.2, so the training data
   no longer moves glycine; the held-out proteins still pull aL down (-4.8 summed). -ln(P_aL / P_aR) -0.425, 0.17
   above ff3.0's -0.598 (rounds 1-4: -0.567, -0.480, -0.451, -0.425; less alpha_L each round, by less each round).
   Map step rms by round 0.104, 0.074, 0.051, 0.044; the held-out gap left after each round's step about 1.00, as in
   every round. Native fidelity over epoch 3: glycine beta' 0.782, alpha_R 0.989; non-glycine alpha_R 0.997, pPII'
   0.906, beta' 0.099 (n 26; free 0.002). Margin +0.042 (step 69) to +0.065 (step 75), dhb -0.575 -> -0.605, E_beta -
   E_alpha -0.036 -> -0.071 over steps 65-75; 24 of 24 proteins every step, KE/1.5kT at most 1.013. **Convergence
   gate at step 76** (`convergence_gate.py`, steps 58-76, run by hand 17:30 10-30;
   `checks/ff31_fix_20261028/gate_step76_ff32_nativefix.txt`): NOT CONVERGED, `dhb` alone pulled (p 0.0045 against the
   cut 0.0050; every other group p 0.82-1.00), the one-scalar H-bond term that fell -0.575 -> -0.605 over steps 65-75.
   **Round 5, ff3.2, at step 95** (read 14:05 10-31; step 95 done 03:10 10-31; `rama_round_05.dat` md5 `ad998fe5`):
   aR -0.231 -> -0.223, aL +0.554 -> +0.547; the epoch-4 pull aR t +0.2, aL -0.0, beta' -0.6, pPII' +0.8, the other
   shared offsets' |t| median 0.4, max 1.7; held-out aL still down (-6.2 summed). -ln(P_aL / P_aR) -0.439, the first
   round to move back toward ff3.0's -0.598 (by 0.014; rounds 1-5: -0.567, -0.480, -0.451, -0.425, -0.439), so the
   glycine balance has stopped moving at 0.16 less alpha_L than ff3.0 (P_aL / P_aR 1.55 against ff3.0's 1.82 and the
   merged ff3.1's 2.28). Map step rms 0.041. Native fidelity over epoch 4: glycine beta' 0.785, alpha_R 0.989;
   non-glycine alpha_R 0.998, pPII' 0.907, beta' 0.100 (n 26). Margin +0.064 (step 76) rising to +0.095 (step 94);
   dhb fell every step, -0.608 -> -0.647; E_beta - E_alpha -0.080 -> -0.076; 24 of 24 proteins every step, KE/1.5kT
   at most 1.013. **Convergence gate at step 95** (steps 77-95, by hand 14:00 10-31;
   `checks/ff31_fix_20261028/gate_step95_ff32_nativefix.txt`): CONVERGED, dhb p 0.027 against the cut 0.005, every
   other group 0.71-1.00. dhb still moved one way every step; the gate passes it under its family-wise cut (0.05 / 10
   groups), and a single-group 0.05 test would not.
   **ff3.1 (`ff31_nativefix`), rounds 1-4** (target 76 reached 11:20 10-31; read 10-31, user: "Check the basin depth
   difference for alpha l and r in GLY for ff3.1"; round 4 md5 `7693f10a`): glycine's shared aR +0.026 -> -0.052,
   -0.123, -0.155, -0.171, aL +0.422 -> +0.498, +0.496, +0.493, +0.495. Pooled glycine map (GLY|right|ALL plus the mean
   change of the 37 GLY|X maps, nonhelix zero) alpha_L - alpha_R by round -1.064 (released) -> -0.904, -0.831, -0.798,
   -0.781 (ff3.0 -0.836, ff2.1 -1.460, ff3.2 at 95 -0.675); -ln(P_aL / P_aR) -0.825 -> -0.666, -0.594, -0.561, -0.545
   (ff3.0 -0.598, ff3.2 -0.439). The epoch-3 pull is gone (aR t -0.3, aL +0.3, beta' -1.3; other shared offsets |t|
   median 0.5, max 1.9); the held-out proteins still pull aL down (-4.8 summed). Native fidelity over epoch 3: glycine
   beta' 0.771, alpha_R 0.989; non-glycine alpha_R 0.997, pPII' 0.906. No error line, KE/1.5kT at most 1.014. With
   every Adam group held, the depths settle 0.11 (L - R) from ff3.2's, which trains the other terms with them
   (`checks/ff32_fold_20261031/gly_depths.py` and `gly_LR_runs.py`, run on the Mac with the repo .venv).
7. [ ] ff3.2's epoch-end panels, the folding test of the fix (epoch 0 read 23:30 10-28: folds 0.449, -0.167 [-0.206,
   -0.129] against released ff2.1, level with `a32_00` and `a31_00`, so the fix does not cure epoch 0's loss; findings 1.37): tag `n32_EE` added to `submit_new.sh` (user 10-28), so
   each epoch end gets the 176-run panel the all-term runs had and reads against `a31_EE`, `a32_EE` and `c9_EE`
   (findings 1.37's table). `extract_ff.py` on the epoch-0 checkpoint writes a `rama.dat` md5-equal to the run's
   `rama_round_01.dat` (`0c9b7978`; `checks/ff31_fix_20261028/extract_test_n32_00/`). **Epoch 1** (`n32_01`, read
   17:10 10-30 by `paired_npz.py`): folds 0.414, -0.201 [-0.246, -0.160] p 1.8e-08 against released ff2.1 (40 of 44
   domains lower); -0.035 [-0.076, +0.004] p 0.22 against its own epoch 0 (mean Q -0.022, p 0.005); -0.066 [-0.110,
   -0.025] p 0.003 against the earlier all-term rerun's epoch 1 `a32_01` (0.481); -0.035 p 0.09 against `a31_01`
   (0.449); level with the control `c9_01` (0.405; +0.009, p 0.98). **Epochs 2 and 3** (read 14:05 10-31): `n32_02`
   folds 0.436, -0.180 [-0.227, -0.136] p 5.9e-08 against released ff2.1 (40 of 44 lower), level with `a32_02` (0.435)
   and `a31_02` (0.462; -0.027, p 0.13), above `c9_02` (0.375; +0.061, p 0.003); `n32_03` 0.431, -0.184 [-0.231,
   -0.140] p 1e-07 (39 of 44 lower), level with `a32_03` (0.432), `a31_03` (0.420) and `c9_03` (0.430), +0.017 p 0.25
   against its own epoch 1. **Epoch 4** (`n32_04` 49686492, 14:00-14:46 10-31, read 15:00): folds 0.447, -0.169
   [-0.215, -0.126] p 9.4e-08 against released ff2.1 (40 of 44 lower), mean Q 0.783 against 0.829; level with its
   own epoch 3 (+0.016, p 0.45) and with the epoch-3 references `a32_03` (0.432), `a31_03` (0.420) and `c9_03`
   (0.430), none of which has an epoch-4 panel.
8. [x] **ff3.2 extended one epoch, to step 95, and converged there** (user 10-30: "Go ahead and continue training ff3.2", after the step-76
   gate failed on `dhb` alone): the same chain (`train_chain.sbatch ff32_nativefix 95`, the run's own `ConDiv.py`
   `54383004`, resuming from step 76's checkpoint), so round 5 is written at step 95. The run has no
   `after_training.sbatch`: at 95 run `convergence_gate.py` by hand; not converged, one more epoch up to 114 as the
   gated runs had, then stop for the user. Read round 5 and `n32_04` as rounds 1-4. The rule fixed 10-30 before the
   panels: the full validation (as ff3.1's, plan.md Phase 30) only if the last epoch's panel folds level with
   released ff2.1 (paired 95% interval of the difference including 0); otherwise ff3.2 stops and findings 1.37
   records that the native-target fix does not cure the all-term loss. Step 95 reached 03:10 10-31; the gate there
   reads CONVERGED (item 6), so there is no extension to 114; round 5 is read (item 6). **Decided 15:00 10-31:
   `n32_04`'s interval excludes 0 (item 7), so ff3.2 stops without the full validation**; findings 1.37 records it,
   with what the loss is (item 9).
9. [x] **What is wrong in ff3.2's folding** (user 10-31: "Check what is wrong in ff3.2 that fails to fold protein";
   read only, findings 1.37 "What the loss is"; `checks/ff32_fold_20261031/`). The domains swell with helices and
   strands intact; ff3.2 makes every state away from the native cheaper relative to it, unfolded frames by 12.8 E_up
   on all 16 Peng proteins (released ff3.1 +1.2); about half through the two scalars the DSE target moves (backbone
   burial scale, H-bond energies), the rest through the side-chain tables both targets move. No code defect found;
   A193's SARW fix is missing from this run's trainer. Read as a temperature, ff3.2 at T acts as ff2.1 6% hotter for
   global stability (about 298 K at the panel's 280 K), but not term by term (packing 0.92, H-bond 0.97, backbone
   desolvation 1.34), so no temperature compensates it; ff3.2 is not deployed (item 8's rule).
   [ ] The causal test, for the user to decide (it departs from the SI's lambda only as a diagnostic): one epoch from
   ff2.1 at lambda 0 against the control's epoch 0 (`c9_00`, 0.390), judged by the 44-domain panel and `peng_full.py`'s
   unfolded shift (minutes per parameter set). If lambda 0 keeps ff2.1's folding, the DSE's weight is the cause.

### Phase 34 - disordered regions in Upside 3: how accurate, and what could improve them (STARTED 2026-10-28, user)
User 10-28: "for upside3, investigate the possibility of enhancing intrinsic disordered region simulation accuracy";
then "record everything in upside3.md. We will implement this later." The record is `upside3.md` 8.4 (what the model
gives a disordered chain, the literature, B0 and its reading rules, candidates I1-I3 and steps, with its "Resume here"
in section 12, for a cold session); this phase tracks only status. User 10-28 after B0: "record this ... I will work
on it later in another session".
Decided by the user 10-28: no explicit electrostatic term (the trained terms already carry electrostatics, so it
would count them twice); salt dependence is outside the model.
1. [x] Literature and the SAXS panel (`scratchpad/idp_20261028/lit/`).
2. [x] B0 on the Mac, 10-28 11:56-19:37, recorded in `upside3.md` 8.4. ff2.1 and ff3.1 make the same error. Both
   are 11-15% compact, both collapse ACTR and A1-LCD into globules, and both leave Cv78 too expanded. 10-11 of 16
   entries are unresolved. The routing is I1 with variants; I2 is not indicated and I3 is undecided.
3. [ ] B1 (ACTR at 298 K, about 30 min on the Mac) and B2 (longer sampling or REMD; the user's choice), as in 4.2's
   Steps.
4. [ ] I1-I3: after B0, in the user's order (`upside3.md` 8.1), on the user's word.

### Phase 35 - Upside 3 on the GPU (track G; STARTED 2026-10-29, COMPLETE 2026-11-06, user)
User 10-29: "start implementing GPU for upside3. The target is to be able to use GPU to run simulations with CPU
trained ffs, without any additional training or modifications. You can add GPU support for hybrid as well. The target
GPU to use is the GPU partition on midway2 and midway3, and on similar clusters"; then "connect to midway2 and midway3
to check their gpu available, and test on them when you implement the code"; tests on midway2 `gpu2` and midway3
`gpu` only, not beagle3. The standard is `upside3.md` 10.2; the approved design (route, architecture, milestones) is
`upside3.md` 10.3 and `~/.claude/plans/swift-roaming-floyd.md`.
Key decisions (reasons in `upside3.md` 10.3): the CPU engine stays the reference and host, a device engine is built
from its constructed nodes; one batch is every system of a process, flattened with pre-offset indices; the potentials'
math is shared through `UPSIDE_HD` and lane-type templates, the loops are per backend; the device files also compile
as C++ (host backend) so the device path is tested on the Mac; 64-bit fixed-point accumulation (deterministic, a
non-finite value reads back as NaN); BP one block per system on the CPU's schedule; one CUDA 11.5 binary from midway2
for K80 to A100; development tree `/beagle3/trsosnic/yinhan/upside3_gpu` (the Upside 3 tree is under its validation).
**Complete 2026-11-06** (user: "hybrid preparation can be on CPU. Also modify the scripts in examples so they can run
on GPU on midway2 or midway3. After that you can mark upside3 GPU as complete"). The hybrid's preparation stages and
NPT runs stay on the CPU, which closes item 7's last box and item 11's rest; the examples that the device runs have
an `engine` switch (`upside3.md` 10.5), checked by stub runs (below). **Resume here (state at 2026-11-06, Mac
Studio).** **Validated and closed 11-03 (user: "close what is
validated"): the GPU runs soluble proteins on ff2.1 and ff3.1 with the CPU's results**, by `upside3.md` 10.2's ladder
(findings 11g): V0 by construction; V1 of every node on the 14 standard configs, the 16 Peng proteins' 320 frames
(native, misfolded globules, both BP tolerances) and mixed batches, on the K80, V100, RTX6000 and A100 (items 10-12);
V2 of 1ubq by seeds and its energy drift (items 10-11); V3's melting at a tenth of the Peng duration and its Q and Rg
distributions (items 10-11). **Closed 11-04: mixed batches, by the mdCATH panel** (item 12 step 5: 89 domains x 4
runs per force field in batches of 88-92 systems; mean TM, folded fraction and Q equal to the CPU's). Its scope:
dynamics of one protein or a batch of different proteins per process; a config with an unported node (`Const3D`, `WallSpring_*`, membrane and
restraint terms) is refused by name; V3 cannot resolve a native-basin potential difference below 1-2 E_up. The
energy gain without the thermostat and KE/1.5kT 1.005 are the CPU engine's, which the GPU reproduces (Known Errors).
Open: the hybrid (1rkl 7.0's V1 and V2 and glpG's V1 and V2 pass, and glpG's TM4 block agrees with the CPU's;
stage 7.1, area per lipid, rigid groups and the barostat). Speed: one V100 with 4 ladders samples as fast as one 14-core broadwl
job on proteins, the K80 slower than 4 broadwl cores; glpG's 48 replicas run 850 tu an hour on a V100, 1,100 on the
A100, against 60 on a 28-core broadwl node; an mdCATH batch of 88 systems of 169-242 residues runs 1.8 times as
fast on an A100 as on a 28-core broadwl node after g9b (2.0 times), and batches of 57-166 residues 1.4-1.5 times
as fast on a V100. Build g9b (g9a with the `--gpu` host work parallel over systems, bitwise g9a's on a V100) is the
GPU tree's `obj` since 00:36 11-04.
- Running: nothing (remote_jobs.md line 5 has no job). The panel's trajectories were deleted 11-04 after step 5's
  read (its npz, logs and TM series kept; remote_jobs.md "Disk").
- The computer: the code is HEAD `be8f652` (g9b's `src/main.cpp`; committed and pushed 10:39 11-04 at the user's
  word, "git push all local changes"; g9a's mixed batches and `run.py`'s `remd_engine` in `3aa607d`); the working tree's `src`, `py` and `tests` equal `scratchpad/gpu_work/dev16_home`'s
  (`dev15_home` is g9a). The gitignored tools are in
  `scratchpad/gpu_work` (`g0/`, and `mdc/` for step 5); the bundle `/project/trsosnic/yinhan/checks/handoff_gpu_20261102/`
  holds their 11-02 state. Device tests on a Mac use a home built by `build_home.sh` (the repo's `obj/` lacks the
  device API); the CPU regress baseline is `dev14_record.npz`. The GPU tree has `obj` = `obj_g9b` (md5 `959ee66d`),
  `obj_g9a` (`e928f0fa`), `obj_g8a` and the earlier builds beside it. Before `sync.sh`, compare `scratchpad/gpu_work/configs` with
  the cluster's `checks/configs` by md5 (the cluster's hold the released ff3.1's and `glpg.stage_7.0.up`, which the
  Mac lacks); single files are sent by `scp`.
- The examples (11-06): `engine = "cpu"` or `"gpu"` in 01, 02, 03, 04, 05, 10, 11, 12, 13 and 15 (06-09 use nodes
  the device refuses); with `"gpu"` the self-submitting ones take one task, 4 cores and one GPU of `gpu_partition`
  (also the QOS: `gpu` on midway3, `gpu2` on midway2), and the in-place ones and the `0.run.sh` scripts add `--gpu`
  inside a GPU allocation. `example/16.MARTINI/run.py` gained `cluster` (midway2 or midway3), so its preparation and
  its REMD blocks run on one cluster: before, `remd_engine = "gpu"` submitted midway2's preparation and a midway3 GPU
  block joined by `afterok`, a dependency no scheduler can hold across clusters. Checked on the Mac by stub runs
  (`upside` and `sbatch` logging their arguments, the device-API home `dev16_home`'s `py`): with `"cpu"` every
  command line equals the one before the edit (13's two scripts but for their random seed), and with `"gpu"` each
  differs only by `--gpu` and the GPU job options; the hybrid's six cluster and engine variants submit as designed.
  No example has run on a GPU node since the edit.
- Left, a deployment rather than validation: the GPU code into the Upside 3 tree once its campaigns end (never while
  they run; ff3.2's glpG REMD runs on it now), then a rebuild on midway2 and `regress.py compare` against the tree's
  build.
- The estimate told the user 11-02 15:30: the mdCATH panel 11-05 to 11-07 and glpG TM4 on the GPU 11-06 to 11-09,
  both done 11-04; rigid groups and the barostat 11-07 to 11-10; the whole about 11-10.
- The hybrid's CPU validation (glpG on Upside 3; Phase 24 item 2b 4, Phase 37; remote_jobs.md line 3): ff2.1 (4 of 4
  blocks) and ff3.0 (3 of 3, TM4 held) read; the released ff3.1's block 1 read, TM4 DSSP 0.74 at T 0.70 (Upside 2's
  ff3.0 0.96). The TM4 HDX readout waits on the user's choice of force field. User choices so far: 11-01 (a), V1 of
  the released ff3.1, and "run both"; 11-02 "continue the validation", "Soluble remainder", "fix the bug you find".

Code in `src/device*.{h,cu}` and the ported nodes'
`make_device`; working files in the gitignored `scratchpad/gpu_work/` (`build_home.sh`/`regress_home.sh` build and
regress a copy of the tree; `dev12_record.npz` is the CPU baseline since A239, 4260 arrays, and `regress_hybrid.sh`
rebuilds a copy and regresses the hybrid fixtures alone; `hyb_dyn/pair.sh` runs one hybrid trajectory on both engines;
`sync.sh` sends the tree to `/beagle3/trsosnic/yinhan/upside3_gpu`, whose `build_gpu.sh TAG` builds it on the midway2
login node; `make_subset.py` writes the development configs in `configs/`; `g0/` holds the job scripts `v1.sbatch`,
`v2run.sbatch`, `v2cmp.sbatch`, `diag.sbatch`, and V3's `v3_run.py`, `v3.sbatch`, `v3_analyse.py`; `ff31_configs.py`
writes the ff3.1 configs; `deadzone_repro.py` is A239's reproducer). Jobs: remote_jobs.md line 5. Every CPU change
is checked with `regress_home.sh` against `dev12_record.npz` (A238 and A239 the only intended moves); every device
node with `tests/device_v1.py` on the host backend (the Mac) and on a GPU.
1. [x] G0 probe (findings 11g): K80 driver 460 (CUDA 11.2) and V100 driver 535 both run one CUDA 11.5 binary built on
   midway2 for sm_37 to sm_86; Threefry bits equal, fixed-point sums bitwise and repeatable. The A100 and RTX6000 jobs
   are still queued.
2. [x] G1 infrastructure: `device.h` (backend, `parallel_for`, block-per-system `per_system`/`block_for`/`block_scan`,
   fixed-point `DeviceAccum`), `device_engine` (graph flattening, the CPU's evaluation order, per-node potentials,
   Verlet and multi-step integration, OU thermostat on the same streams, an error flag for failed preconditions),
   `upside --gpu` (host work at frames, pivot moves and replica exchange on the CPU engines), the C API and
   `DeviceUpside`, `tests/device_v1.py`, CMake `UPSIDE_CUDA`. The host backend fills new memory with NaN.
3. [x] G2 bonded and backbone nodes (Distance3D, Angle, Dihedral, Spring, rama_coord, cat_pos, affine_alignment,
   backbone_pairs, infer_H_O, protein_hbond, hbond_energy, rama_map_pot), with the device pair list and interaction
   graph (`device_interaction_graph.h`) and the pair potentials shared in lane-templated headers.
4. [x] G3 placement, environment and coverage nodes (placement_*, weighted_pos, environment_coverage,
   hb_environment_coverage, hbbb_coverage, hbond_coverage, sigmoid_coupling, bb_sigmoid_coupling). V1 passes on the host
   backend for every node of 1ubq, chig and 1dfn without rotamer; on the K80 it passes but for `hbond_energy`'s
   exact-zero test (item 8).
5. [x] G4 rotamer BP (`device_rotamer.cu`): residue pairs by a dense lookup per pair class, numbered in (residue 1,
   residue 2) order; 1xN bead pairs folded into the other residue's energies and 11 pairs into the energy, as
   `fill_holders`; one thread block per system runs `solve_for_marginals`'s schedule (its first-update quirk included),
   each residue's messages multiplied in holder order 33, 36, 66; Bethe energy and envelope derivatives. V1 of the full
   ff2.1 graph passes on the host backend (chig, 1ubq, 1dfn) and on the K80 (chig; 1ubq, 1ubq_rotamer and 1dfn but for
   item 8, before its fix). Configs with nodes not yet ported (`Const3D`, `WallSpring_*`, membrane and restraint terms) are refused
   with their names; V1 reports them as skipped.
6. [x] V2 (`tests/device_dyn.py`, 14-replica REMD of 1ubq, CPU against `--gpu` from one start and seed) passes on
   the K80 and the V100, with and without rotamer (10-29): every replica's KE/1.5kT, potential and N-CA, CA-C, C-N
   means within 2.9 block standard errors (70 per run, one seed each), KE/1.5kT 1.00 +- 0.01 throughout. Speed is not there yet: per
   system and step the K80 takes 323 us (no rotamer) and 837 us (full) against 87 and 253 on 4 broadwl cores, the V100
   101 and 208 against 75 and 181 on 4 Cascade Lake cores. G5's speed work (launch count, host syncs) and V3 follow
   the hybrid.
7. [x] G6 the hybrid (`device_martini.cu`). Revised decision (10-29): `martini_potential` keeps the CPU's candidate
   table, the CSR adjacency over its non-zero pairs with the hybrid's static pair filters applied once on the host,
   in place of the class-by-class index of `upside3.md` 10.3 item 10: the same pairs by construction, and glpG's table
   is only 8.3 M pairs (33 MB per system). The skin cache is a slab per atom of its CSR row, and
   `martini_sc_table_1body`'s a slab per rotamer row of the environment atoms, so neither has a capacity to exceed.
   `martini_interaction.h` holds the pair term the CPU node and the device share; the side-chain table and the two
   springs keep their CPU code and the device has a copy, since sharing them moved the CPU's forces at rounding
   (findings 11g). Dynamics: masses, fixed and z-fixed atoms in both integrators, the g-JF step, the thermostat on the
   g-JF bath's complement, the PO4 hold's set changes uploaded after the host refresh; rigid groups and the barostat
   stay refused.
   - [x] shared pair term; the whole tree's CPU regress 0 of 4260 arrays against `dev6_record.npz` (`dev8_record.npz`)
   - [x] device nodes `martini_hybrid_position`, `martini_potential`, `martini_sc_table_1body`, `dist_spring`,
     `angle_spring`: V1 passes on the host backend for 1rkl 7.0 and 7.1 (the 1.0 A frames past the fixed-point range,
     findings 11g)
   - [x] `DeviceDynamics` and `upside --gpu` for the hybrid: one start and seed, 3.05e-5 A apart after one round
     (rounding), KE/1.5kT equal to 1e-6 over 100 rounds (host backend)
   - [x] build g5a on the GPUs: V1 passes every node of the 11 fixtures, the hybrid's included, on the K80 and V100
   - [x] V2 of the hybrid's GPU lipid excess found and fixed (A239, findings 11g): the g-JF position update's force
     term, below half an ulp of a large coordinate, rounded away where the compiler kept the source order (nvcc;
     clang for z); the update is now x + b dt/m (p + dt/2 f + beta/2) on the CPU and the device. The K80 with the
     fix matches the CPU without exchange (lipid KE +0.1 SE); regress moves only the hybrid dynamics (12 of 4260)
   - [x] V2 of the hybrid on build g7b passes (REMD of 1rkl 7.0, 16 CPU, 16 K80, 8 V100 seeds): Bonferroni 1 on both
     GPUs, the lipid KE of every replica within 1.8 SE (findings 11g)
   - [x] V2 of stage 7.1; rigid groups and the barostat on the device: not ported (user 11-06, "hybrid preparation
     can be on CPU"); the device refuses them by name and they run on the CPU
8. [x] `hbond_energy`'s exact-zero skip removed (user 10-29: "fix the issue"), ledger A238: its sensitivity is the
   energy per H-bond at every count, so backends no longer flip it.
9. [x] G5 speed, step 1: parallelism within an element. nvprof on the K80 (10-30, findings 11g): the time is in
   kernels, not launches or transfers, and in each slow kernel one thread loops serially over hundreds of partners:
   the pair lists' rebuild (1ubq without rotamer: 450 of 2000 us per evaluation) and cut (210 us), BP's incidence
   lists (900 us of full 1ubq's 10.9 ms; BP itself 3.2 ms), the hybrid's `martini_sc_table_1body` (5.2 of 11.4 ms)
   and `martini_potential` (2.1 ms). Decision: `device.h` gets `lanes_for`, one warp per element (one lane on the host
   backend), with an ordered compaction (ballot and prefix count, so the order is the serial one), an in-order float
   sum (the serial order, by shuffles) and `AccumLocal`, a thread's fixed-point partial sum flushed once (integer sums,
   so the same bits as one atomic per term). Every kernel changed keeps its results bitwise; the check is a `--gpu`
   trajectory of each build, compared bit for bit.
   - [x] `lanes_for`, `AccumLocal`, `slots_for` (a launch over a pair list's used slots) and a warp-shuffle
     `block_scan` in `device.h`; the pair-list rebuild, `martini_potential`, `martini_sc_table_1body`, BP's incidence
     lists and the 12 capacity-sized kernels use them (builds g6a, g6b)
   - [x] bitwise against g5b on the K80 (1ubq without and with rotamer, 1rkl 7.0) and on the host backend; per system
     step 360 to 245 us, 867 to 581 us, 14.5 to 6.5 ms (findings 11g)
   - [x] BP's solve (46% of full 1ubq): its passes fused, 3 synchronisations per update instead of 12 (g6c, bitwise
     against g5b on the K80; full 1ubq 581 to 526 us per system step, BP 3.55 to 2.78 s of the profile)
   - [x] g7b (g6c plus A239) installed as the tree's `obj/` (10-30 12:15)
10. [x] The validation ladder completed on g7b (user 10-30, "run them": V1 on ff3.1 and on every GPU type, the protein
    V2 by seeds, then V3).
    - [x] V1 on ff3.1 (chig, 1ubq, 1dfn written as `tests/fixtures.py` does, with ff_3.1 and its `rama.dat`; the same
      31 nodes as ff2.1) and the 11 fixtures: every quantity passes on the host backend, K80, V100, RTX6000 and A100,
      with the same worst errors on the four GPUs (findings 11g); the released ff3.1 (11-01, the three configs
      rewritten, only the Ramachandran tables differ) passes every quantity on the host backend and the V100
      (61986517), the same worst error per config as before
    - [x] V2 of the protein by seeds (`tests/device_dyn.py`, now with the backbone angles): 1ubq 14-replica REMD,
      2000 tu, ff2.1 without and with rotamer and ff3.1, 16 CPU, 16 K80, 8 V100 seeds; Bonferroni 0.087 to 1 over the
      112 comparisons of each config and GPU, and the per-seed means over replicas within 2.0 SE (findings 11g)
    - [x] V3 (REMD melting, CPU against GPU, ff2.1 and ff3.1). Scale decided by the user 10-30: a tenth of each Peng
      native duration, since one full arm takes about a week on a V100. The BP validation's 4 proteins (proteinG,
      homeodomain, WWdomain, NTL9), native start, the Table S2 ladders, dt 0.009, frames every 100 tu, exchanges every
      10; 4 ladders per engine and force field: a CPU job runs 1 (14 broadwl cores), a V100 job 4 (independent by their
      streams, exchanges only within a ladder). `g0/v3_run.py` (from `bp_validation/bp_run.py`, build `obj_g7b`,
      chunked and resubmitting), `v3.sbatch`, `v3_analyse.py` (Q, CA-RMSD, Rg, potential per temperature and T_mid,
      Welch's t over ladders, Bonferroni); a smoke test at 400 tu ran the resume path and the analysis (10-30). All 40
      arms complete in one chunk each (10-31 22:35): Bonferroni 0.12 to 1 per protein and force field, 0.93 over the
      456 comparisons; T_mid within 0.009 (3 K) and 2.2 SE in every pair (findings 11g)
11. [x] The rest of `upside3.md` 10.2's ladder. User 11-02 ("continue the validation", then "Soluble remainder"):
    the soluble part first, needing no device code; each reading fixed here before its data. The soluble part is
    closed 11-03 (Resume here); the hybrid's and the mdCATH panel's rungs run as item 12.
    - [x] V3's Q and Rg distributions, from item 10's 40 arms (no new runs). `g0/v3_dist.py`: per ladder,
      temperature and observable, the deciles 0.1-0.9 over the frames after the first 20% (v3_analyse.py's frames);
      each decile's CPU and GPU ladders (4 and 4) compared by Welch's t, as V3's means were. The distributions agree
      if the smallest p times the number of comparisons is above 0.05 per protein and force field and over all of
      them. Beside, not deciding: each temperature's two-sample KS distance between the engines' pooled frames.
      **Passes** (49746180, 11-02): Bonferroni 1 for each of the 8 protein and force field pairs and over all 1754
      comparisons (smallest p 0.013); the refactored `v3_analyse.py` reproduces `v3/analysis.txt` (empty diff).
    - [x] V1 at scale: the 16 Peng proteins with ff2.1 (Upside 3's Peng runs, `upside3_validation/peng/runs/`;
      the config is the native arm's replica 0 `/input`), 10 frames each, the last frame of replicas spread over
      the ladder: 4 from the native arm (0, n/3, 2n/3, n-1) and 6 from the de novo arm (0-3, n/2, n-1), the low
      ones misfolded globules where the protein does not fold; each frame's CA-RMSD to native and Rg are reported.
      Each config twice: BP at the production tolerance (1e-3) and at 1e-6 (max_iter 1000 both). `device_v1.py`
      gains `--frames` (frames per config from an npz) for it; `g0/v1x_prep.py` writes the configs and frames.
      Passes if every quantity of every frame passes V1's test at 1e-6. At 1e-3 a failing quantity is accepted
      only if the same frame passes at 1e-6 (the stop test fired at another iteration, a known BP property) and
      the frame's potential differs by at most 7.5e-4 E_up and its derivative by at most 8e-5 of its RMS, the BP
      stop-test fix's precedent (`upside3.md` 10.2); otherwise V1 fails.
      On g7b it **failed by this rule** (62076233, V100, 11-02): 30 of 32 configs pass every quantity at both
      tolerances (no BP failure at 1e-3; worst rotamer output error 6.9e-6). BBA and NTL9 fail one quantity each, at both
      tolerances: `hbond_energy`'s output at native replica n-1 (the same frame counted as frame 2's system 1 and
      frame 3's system 0), relative errors 1.67e-2 and 1.50e-3, rounding scale 0. Cause, measured on the CPU (Mac):
      those frames' H-bonds are weak, and both engines form an occupancy as `1.f-expf(-s)` (`hbond.cpp:185`,
      `device_hbond.cu:155`), which for s ~ 1e-6 keeps only whole steps of 2^-24 (BBA's two bonds at exactly 56
      and 1 steps, NTL9's at 571, 70 and 7). One ulp of `expf` between CUDA and libm is one step, 1/56 of BBA's
      energy; the absolute error is 2.2e-7 E_up in both, on frame potentials of 70 and 119 E_up. The CPU's value
      is flat under 4-ulp moves, so the test has no rounding allowance for it. Fixed as A240 (user 11-02: "fix the
      bug you find"): the steps come from two subtractions, `1-hb` per pair (`hbond_interaction.h`) and
      `1-expf(-s)`; both engines now use `log1pf` and `expm1f` (`expm1f` alone left every step in place). On the
      Mac the reproducer, node_fd and regress pass (findings 3.14 A240), and V1 on the host backend passes BBA and
      NTL9 at both tolerances and the 14 standard configs. **Passes on g8a** (V100, read 17:20 11-02): 62096441,
      the 32 configs with `--frames`, passes every quantity at both tolerances (1280 per config; BBA's worst error
      now 8.0e-6, NTL9's 1.4e-5, the largest of all 2.9e-5, WWdomain's `hbond_coverage`), and 62096444, the 14
      standard configs, passes every quantity with the same worst error per config as g7b's 61891693, except that
      the three 1dfn configs lose their worst, `hbond_energy` at 1.5e-4 (now 6.9e-6, 2.6e-7 and 6.9e-6). g8a is
      installed as the tree's `obj` (17:22 11-02; g7b kept as `obj_g7b`).
    - [x] V2's energy drift: `device_dyn.py nve`, N copies at T 0.85 without exchange, the thermostat thermalizing
      the start only (its interval beyond the run). Per copy over the whole run: the drift of the total energy
      (potential plus kinetic; least-squares slope, E_up per 1000 tu), its RMS about that line, the mean potential
      and KE/1.5kT. 1ubq ff2.1 with rotamer and 1ubq ff3.1 (`configs/1ubq.up`, `1ubq_ff31.up`, item 10's V2
      configs), 5000 tu, 4 copies per seed, 16 CPU seeds (broadwl) against 8 V100 seeds. They agree if `compare`'s
      Bonferroni is above 0.05 for each config. **Passes** (g7b; CPU 49746182-83, V100 62076234-35, read 16:55
      11-02): Bonferroni 0.65 (full, largest difference 2.2 SE) and 1 (ff3.1, 1.5 SE). Not decided by the rule:
      both engines gain about 90 E_up per 1000 tu per copy and KE/1.5kT averages 1.12-1.14, so the drift is the
      model's or the integrator's, not the device's (findings 11g); its cause is not measured. The V2 and V3
      results on g7b stand for g8a: A240 moves each occupancy by at most 7.6e-8 on both engines.
    The harness changes passed on the host backend first (MacBook Pro, `scratchpad/gpu_work/dev13_home`: BBA's two
    configs pass every quantity; the 1e-6 tolerance moves BBA's potentials by up to 8e-4 E_up; `nve` on chig keeps
    the total energy within 1 E_up RMS of its line, against 7 with the thermostat on). Submitted 13:40 11-02.
    Still open after these: V2's area per lipid and stage 7.1 (item 7). Closed by item 12: V1's hybrid frames beyond
    1rkl (glpG's), glpG TM4 and the mdCATH panel. Corrected
    11-02: glpG's TM4 REMD does not need item 7. Production sets no rigid group (`martini_hybrid.cpp`
    `rigid_protein_stage`), and glpG's `stage_7.0.up` (`upside3_validation/glpg_ff31`) has no `/input/barostat`;
    replica exchange refuses the barostat anyway. It needs (c) and glpG's own V1 and V2. Item 7's rigid groups and
    barostat serve the preparation stages and area per lipid, an NPT observable.
12. [x] (c) the production workflows on the GPU (user 11-02: "plan that, while still monitoring"). Approved 11-02 with
    the user's choices: glpG's TM4 on ff3.1, mdCATH on ff2.1 and ff3.1, the largest system's capacity, steps 1-2 now
    (steps 3-5 wait for the user's word). Goal: run the two workflows item 11 still owes, the glpG hybrid REMD and the mdCATH panel, on a
    GPU from the workflow scripts, each tested against the same workflow on the GPU tree's CPU engine.
    Measured for the plan (11-02):
    - glpG's `stage_7.0.up` (`upside3_validation/glpg_ff31`): 4712 particles, 25 nodes, all with a device port; no
      `/input/barostat`, no rigid group; `martini_potential`'s candidate table 7,491,136 pairs, which the device shares
      between replicas. Its CPU REMD (48 replicas on 28 broadwl cores, Upside 3 tree) runs 2,070-2,180 tu per replica
      in a 35.25 h block, 11 ms per system step; `run.py`'s comment of about 1,000 tu per block is stale.
    - The mdCATH panel (`ff3_selection/panel.py`): 89 domains x 4 seeds, 356 constant-temperature runs (T 0.8 after an
      anneal from t 96 to 400, 10,000 tu, frames every 10, pivot moves every 5, dt 0.009), one single-thread `upside`
      per run, 28 per node, on the ConDiv tree's Upside 2 binary; TM-score by `checks/converge_20261020/tm_panel.py`.
      A panel is 3.56 M system tu (6.3 h on one 28-core node). No Upside 3 run of it exists.
    - Different proteins cannot share one `--gpu` process (host backend, chig, 1ubq and 1dfn on ff3.1): the device's
      interaction graph requires equal `max_n_edge`, whose default is each system's full pair count, and
      `rama_map_pot` shares its maps, one 72x72 layer per residue. With both checks relaxed in a throwaway copy the
      batch runs, so no other node refuses; the interaction parameters are equal across the three. `backbone_pairs`
      already sizes every system's slab by the largest system's capacity.
    - Runners: the hybrid's `run.py` calls `upside` itself (one process of 48 replicas, sbatch `--wrap` from one
      `slurm_args_file`). The soluble runners are cluster files (`ff3_benchmark/bench_run.py`, its copies in
      `upside3_validation/peng/` and `bp_validation/`, and `g0/v3_run.py`, the only one with `--gpu`).
      `ff3_benchmark/bench_run.py` and `ff3_selection/panel.py` serve running campaigns (ff3.1's Peng arms, ff3.2's
      panels) and are not edited.
    Design (reasons):
    - The device engine's one change: batches of different proteins. Every system's pair slab takes the largest
      system's `max_n_edge`, as `backbone_pairs` does, and each system's own capacity is still tested, so the device
      refuses what the CPU refuses. `rama_map_pot` takes each system's maps, uploading identical maps once, so a batch
      of one protein keeps today's arrays. A batch of one protein stays bitwise (item 9's trajectory check).
    - The hybrid: `run.py` gains `remd_engine = "cpu"` or `"gpu"` in its stage 2 settings, with the GPU job's Slurm
      options and environment beside the CPU ones. Preparation stays on the CPU (the device refuses rigid groups and the
      barostat). A GPU block runs the same command with `--gpu`, on one GPU with 4 cores for the host work (frames,
      exchanges), under `env_upside3_gpu.sh` until (d); ladder, seeds and block restarts are the CPU's.
    - mdCATH: a runner in the GPU tree, `checks/mdc/mdc_run.py`, with `panel.py`'s configs, schedule, frames and seed
      count and an engine argument. A GPU job runs domains of similar size, 4 copies each, in one `upside --gpu`
      process; each run keeps its own output file, so `tm_panel.py` and `paired_npz.py` read it unchanged. The
      reference is the same runner on the GPU tree's build without `--gpu`: the same code, so the comparison isolates
      the device (the Upside 2 panels differ from it by the audit).
    - Not in (c): GPU options in the Peng runners (V3 covered Peng's REMD and no remaining test needs them), (b), (d).
    Steps, each test fixed before its data; every job is listed here for approval:
    1. [x] Mixed batches on the device, on the Mac's host backend first: the change above; CPU regress 0 of 4260
       arrays moved against `dev14_record.npz`; `tests/device_v1.py` gains a batch mode testing every system of one
       mixed process against its CPU engine, run on chig, 1ubq and 1dfn (ff2.1 and ff3.1). On a V100: a 1ubq
       trajectory bitwise equal to g8a's, and the mixed V1 (2 jobs, under 10 min each).
       Mac part done 18:12 11-02 (`scratchpad/gpu_work/dev15_home` against `dev14_home`, HEAD): V1 of the 14 standard
       configs prints the same lines as HEAD's; the mixed batches pass every quantity (ff2.1 worst 5.2e-6, ff3.1 the
       same); host-backend `--gpu` trajectories of 1ubq with exchange and of 1rkl 7.0 are bitwise HEAD's; CPU regress
       0 of 4260 arrays. Build g9a on the GPU tree (remote_jobs.md line 5); on the V100 19:25-19:30 11-02: 1ubq (14
       systems, 101 frames) and 1rkl 7.0 (4 systems, 41 frames) bitwise equal to g8a, at g8a's speed; the mixed V1 passes
       every quantity (ff2.1 and ff3.1, 192 per protein). g9a installed 21:54 11-02, after the `/beagle3` inode quota
       that blocked file creation 19:46-21:51 was reclaimed.
    2. [x] The hybrid's GPU option and glpG on the device. `run.py` as above, its stale comment corrected. glpG's V1
       frames: 10 production frames from `glpg_ff31` block 1 (replicas across the ladder, beginning and end), at BP
       tolerances 1e-3 and 1e-6, by item 11's V1 rule (1 V100 job, under 1 h). A probe of 48 glpG replicas with
       `--gpu` for 20 tu on a V100 and an A100: device memory and time per system step, which set the GPU type and a
       block's length (2 jobs, under 30 min).
       `run.py`'s option written; glpG's frames `checks/configs/v1g/` (`g0/v1g_prep.py`); V1 62105839 and the probes
       62105864 (V100) and 62105868 (A100, both by `g0/glpg_probe.sh`, through `run.py`'s own submission), 18:20.
       glpG's V1 passes all 1040 quantities at both tolerances (worst rotamer output 2.0e-6 and 2.4e-6, 4 times their
       rounding scale). The V100 probe ran 210 tu of 48 replicas in 903 s through `run.py`'s GPU path: 794.7
       us/system/step, 850 tu an hour against the CPU block's 60 (14 times one 28-core broadwl node), 4,382 of 16,384 MiB,
       KE/1.5kT 0.993-1.017; then `run.py` could not write `block_count` (the quota). A CPU-length block of 2,100 tu is
       about 2.5 h on a V100. The A100 probe (62105868 failed at its start on the same quota; rerun 62120973, 23:43):
       628.9 us/system/step, 1,100 tu an hour, 4,495 of 40,960 MiB, KE/1.5kT 0.994-1.014, so a 2,100 tu block is about
       1.9 h on an A100. glpG fits either GPU; the A100 is 1.27 times the V100. Steps 3-5 started 11-03 at the user's
       word, on the V100.
    3. [x] glpG V2 by seeds (`device_dyn.py`, as 1rkl's): REMD of stage 7.0 with 4 replicas over 200 tu, 16 CPU seeds
       on broadwl (about 1 h each on 14 cores) against 8 V100 seeds (`v2seq.sbatch`, about 1 h); KE/1.5kT, protein
       and lipid KE, potential, C-N and lipid diffusion per replica; passes if `compare`'s Bonferroni is above 0.05.
       Submitted 18:25 11-03 on g9a: the config `checks/configs/glpg.stage_7.0.up` is `glpg_ff31`'s stage 7.0
       checkpoint promoted as `run.py` does (its input `pos` and `mom` bitwise the CPU REMD's `remd_start.up`); CPU
       seeds 5001-5016 with 4 cores each (one OpenMP thread per system, so more cores add nothing; about 1.6 h by
       the CPU block's speed), V100 seeds 5201-5208 in one job; work `g0/v2s_glpg/`. **Passes** (read 22:10 11-03,
       `g0/v2cmp_glpg_20261103.txt`): Bonferroni 0.94 over the 48 comparisons, the largest difference 2.6 SE
       (replica 3's CA-C-N angle), lipid KE/1.5kT within 1.2 SE, lipid D within 0.8 SE and the potential within 1.5
       SE in every replica. Per system step the V100 took 4.2 ms at 4 replicas, the CPU 50 ms on 4 broadwl cores.
    4. [x] glpG TM4 REMD on the GPU: `run.py` with `remd_engine = "gpu"` from the CPU run's `stage_7.0` checkpoint
       in a new run directory, 48 replicas, blocks of the CPU's length (about 2,100 tu). Phase 24's readout (TM4 DSSP,
       KE and C-N at T 0.70) against the CPU's block. One trajectory per engine cannot test agreement sharply (V2
       does); a difference beyond the CPU block's half-to-half spread (ff3.1: 0.78 / 0.70) is investigated. Jobs: one
       GPU job per block, its length from step 2's probe.
       Submitted 18:30 11-03 by `g0/glpg_tm4.sh v100 04:00:00` (one block, `run.py`'s own submission; its settings
       those of `glpg_ff31`'s `run.py` but the engine): 3.5 h of REMD, about 2,975 tu at the probe's speed, so the
       block covers the CPU block's 2,150 tu (216 frames); both are read over t 0-2150 by
       `upside3_validation/glpg/glpg_u3_block1.py`, the GPU block's remainder beside.
       **Agrees** (read 23:40 11-03; `g0/glpg_tm4_v100_t0-2150.txt`, `_t2160-2790.txt`, `glpg_tm4_rungs_t0-2150.txt`):
       62200634 ran 2,790 tu (845 us per system step, 800 tu an hour), KE/1.5kT 0.998-1.003 over the 48 replicas.
       At T 0.70 KE/1.5kT 1.000 / 1.005 (protein / lipid; CPU 1.004 / 1.005), C-N mean 1.314 sd 0.120, max 1.80 A
       (CPU 1.312, 0.120, 1.79). TM4 DSSP 0.653 (halves 0.735 / 0.572) against the CPU's 0.741 (0.784 / 0.697): a
       difference of 0.087, the CPU halves' spread, so investigated by `g0/glpg_tm4_rungs.py`, every rung of both
       blocks. The GPU-CPU difference changes sign along the ladder (rungs 0-15 -0.046, 16-31 +0.042, 32-47 +0.055;
       19 of 48 below), and the ladder mean, which the exchanges cannot change, is 0.706 against 0.689 and falls
       alike by quarter (GPU 0.787, 0.729, 0.675, 0.632; CPU 0.766, 0.694, 0.659, 0.638). So T 0.70's difference is
       in which configurations occupy the cold rungs, not a bias of the engine. T 0.795 reads 0.717 against 0.705.
    5. [x] mdCATH on the GPU: `mdc_run.py`; the CPU reference (356 runs on broadwl, about 7 h on 28 cores per force
       field) and the GPU panel (356 systems in size-binned batches on V100s, its cost from a one-batch probe, which
       also times the host work: frames and pivot moves run serially over systems in `main.cpp`'s `--gpu` loop, and
       if they take more than 10% of the run that loop is parallelized as the CPU branch is, bitwise). The engines
       agree if, over the 89 domains, the paired bootstrap 95% interval of the mean TM-score difference (seeds pooled
       per domain, `tm_panel.py`) contains 0 and the paired Wilcoxon p is above 0.05; the folded fraction is reported
       beside, not deciding.
       `checks/mdc/mdc_run.py` (from `scratchpad/gpu_work/mdc/`): `panel.py`'s upside_config options and
       `run_upside` schedule, the tree's ff2.1 and ff3.1 files (md5 equal to the panels' `mdc_ff21_released` and
       `mdc_ff31nf`), the 89 domains sorted by length in 4 batches (57-95, 95-128, 129-166, 169-242 residues; 88-92
       systems), each batch one process of the same build on either engine, seed = batch number; each run's
       `<code>.s<k>.up` and `panel._analyze`'s npz (its phi/psi check against the engine included), tags
       `u3_ff21_<engine>` and `u3_ff31_<engine>`. `panel.py` imports `rama_basin`, which the GPU tree lacks; it is
       read from the ConDiv tree's `training/` (numpy only). Tested on the midway2 login node (2 domains, 20 tu, both
       force fields). Submitted 18:45 11-03: the CPU reference's 8 batches (28 broadwl cores each) and the probe on
       batch 4 of ff3.1 (500 tu, V100), whose result sets the GPU panel's 8 jobs. The CPU reference is complete
       (49777036-43, 1.6 to 3.4 h a batch; 57 to 124 us per system step on 28 cores, batch 1 to 4, equal for the two
       force fields): 356 runs per force field, each to t 9990, no non-finite potential or coordinate, KE/1.5kT after
       t 2000 1.006 on average (1.000 to 1.011 per run, the CPU dynamics' known 1.005), the last frames' longest C-N
       1.89 A (ff2.1) and 1.94 A (ff3.1).
       The V100 probe 62201376 failed 22:50 11-03 when upside started, its log lost with the probe's temporary
       directory (`mdc_run.py` now prints the end of upside's log on a failure; md5 `41f5e515`). On the Mac
       (`dev15_home`, same code) batch 4's 88 systems run on the host backend, and 22 of them for 12 tu (pivot moves
       included) give the CPU engine's KE/1.5kT (2.1-4.0 at T 0.04 before the anneal, on both engines) and no
       refusal. The host backend's peak footprint is 32.1 GB against the CPU engine's 12.0 GB: about 20 GB of device
       arrays, over the V100's 16. Each interaction graph's slab is the largest system's capacity, the full pair
       count n1 x n2 by default, for every system: 227 M slots at batch 4 (environment and H-bond coverage 177 M,
       rotamer 44 M), about 16.5 GB, plus about 3.5 GB of Ramachandran maps (22 domains' own, 2 nodes). The slab
       grows as the longest domain's length squared, so batches 1-3 take about 3, 6 and 10 GB and fit a V100; batch 4
       goes on an A100 (40 GB). The probe was resubmitted 23:49 11-03 on an A100 (62225254).
       The A100 probe (00:01-00:17 11-04): 19,755 MiB of device memory (the Mac's estimate 0.6 GB low); 500 tu in 335
       s, 68.5 us per system step, 1.8 times the CPU reference's 123.2 us on 28 broadwl cores; without pivot moves
       and frames 55.4 us, so the host work was 19.1%, over the 10% rule. The host loop now runs the systems in
       parallel through the CPU branch's own loop (`main.cpp` `for_each_system`, which both branches call; the
       systems' frames, re-centering and moves use only their own state and streams, and the logger's writes are
       already serialised). On the Mac (`dev16_home`): CPU regress 0 of 4260 arrays against `dev14_record.npz`; 4
       threads against `dev15_home` (`mc_bitwise.sh`, frames and moves every tu at T 0.8, pivots accepted) bitwise
       equal on the CPU engine (chig, 1ubq, 1dfn and two mdCATH domains, 2 copies) and the host backend (chig, 1ubq,
       1dfn, 3 copies; three mdCATH domains, 2 copies). Built g9b on midway2 (upside md5 `959ee66d`), installed as
       `obj` 00:36 11-04. Batch 3's device memory measured on the Mac: 10.4 GB (19.7 GB footprint against 9.2), so
       batches 1-3 go on V100s. Submitted 00:38 11-04: g9b's bitwise check against g9a on a V100 (62227671:
       `g0/bitwise.sbatch`, which now fails its job on a difference; 1ubq ff3.1 x 8, frames, swaps and pivot moves
       every tu, 4 threads), and behind it (`afterok`) the panel's 8 jobs 62227672-79 (ff2.1 then ff3.1, batches
       1-4; 1-3 V100 with 4 h, 4 A100 with 3 h, 4 cores each), tags `u3_ff21_gpu`, `u3_ff31_gpu`. Expected about 1.5
       h a batch. The check passed (03:46: 21 frames of 8 systems bitwise equal). The panel ran 04:57-08:29 (1.2-1.9 h
       a batch): 39.5, 51.5 and 65.4 us per system step for batches 1-3 on V100s (CPU 56.5, 76.9, 91.5; 1.43, 1.49,
       1.40 times) and 60.6 for batch 4 on the A100 (CPU 123.5; 2.04 times; 68.5 in the probe before g9b). Checks
       as for the CPU reference: 356 runs per force field to t 9990, no non-finite potential or coordinate,
       KE/1.5kT after t 2000 1.006 (1.001-1.011 per run), the last frames' longest C-N 1.85 A (ff2.1) and 1.84 A
       (ff3.1). Folded fraction, GPU - CPU over the 89 domains (`paired_npz.py`): ff2.1 -0.002 [-0.023, +0.019]
       p 0.91, ff3.1 -0.003 [-0.030, +0.023] p 0.92; mean Q +0.000 for both; every length bin's and fold stratum's
       interval contains 0. **The engines agree by the rule** (`tm_panel.sbatch`, 49803361-62, 08:33-10:36; logs
       `checks/mdc/logs/tm_mdc_u3_<ff>_gpu_cpu_<id>.out`): mean TM, GPU - CPU over the 89 domains, ff2.1 -0.000
       [-0.005, +0.005] p 0.86 (0.727 both; 43 domains lower, 46 higher), ff3.1 +0.000 [-0.005, +0.006] p 0.91 (0.731
       against 0.730; 44 and 45). TM > 0.5: -0.006 [-0.021, +0.010] p 0.33 and -0.006 [-0.022, +0.010] p 0.53. Of the
       28 length-bin and stratum rows (2 measures, 2 force fields) one is below 0.05: ff2.1's mixed stratum, TM > 0.5
       0.999 against 0.992, +0.007 [+0.000, +0.016] p 0.046, the GPU higher, within the 1.4 expected by chance.
    6. [x] Report to the user what runs on a GPU and how fast, and update `upside3.md` 10 and the readiness statement
       (10:50 11-04).
    Estimate (11-02): steps 1-2 11-03 to 11-04, step 3 11-04, steps 4 and 5 11-05 to 11-07.

### Phase 36 - side-chain pair tables coarse-grained from an all-atom force field (STARTED 2026-10-29, user)
User 10-29: derive the SC-SC spline table bottom-up from an all-atom force field (e.g. Amber), "like what we did in
hybrid model", train only "the correction for this CG spline table", and use the CG step to find which side chains
need more than the one-bead table; then "you need to make corresponding change in hybrid spline table as well". The
record is `upside3.md` 7.4 (the model today, routes, the F0-F5 diagnostic, the geometric screen, the hybrid's
corresponding change, open questions, BU1-BU6); this phase tracks only status.
1. [x] Investigated and recorded (10-29), with the geometric screen (`scratchpad/sc_bottomup_20261029/`). Found on
   the way: `rotamer_center_fixed` is a trained output of Jumper's fit, not the C++ tools' (up.md 2.1 corrected).
   Decided by the user 10-29: each side chain stays one CG bead, with more spline tables for orientation and the
   like, never a second bead (findings 10.43); the correction is a table of the CG table's dimension, zero if the CG
   is exact; the force field is for me to find (item 3).
2. [ ] BU1: the user's remaining choices (`upside3.md` 11, item 16): the route, the reference condition, the
   threshold, the stage.
3. [x] The all-atom force field (10-29): DES-Amber with TIP4P-D, ff19SB/OPC as the pilot's second. The pilot checks
   both against the salt-bridge association constants and, for buried pairs, against water-to-cyclohexane transfer
   and the N-methylacetamide dimer in CCl4 (`upside3.md` 7.4).
4. [x] The user's two points (10-29), measured and recorded in `upside3.md` 7.4 ("Tables against sites", "Buried
   and exposed pairs"). The 136 GB was the dense grid only. Most contacts are buried, so two limits mixed by burial
   are proposed (findings 10.44).
5. [x] The user's 10-30 points, measured and recorded in `upside3.md` 7.4 (`scratchpad/sc_bottomup_20261029/sites/`).
   * Sites only "if more than 1 bead can make it more accurate and faster". Sites run in today's engine (exact
     against ff2.1 when co-located) and cost +21% (ubiquitin) and +39% (CA II), against my 10-29 estimate of +1% and
     +21% (findings 10.45). On a steric proxy one bead with a frame (F4) is more accurate than sites for all five
     residues, at a cost of at least +17% and +21% on the pairs that involve them. One bead stays.
   * F4 still needs a fit that places no false well (MET-LEU) and a second direction in the placement node (BU3).
   * Implicit water: no model can be the reference, so W_water stays explicit. PB or an accurate GB is a candidate
     only for U_dry's electrostatics, after a check on buried charges.

### Phase 37 - deploy the fixed-target ff3.1 and validate it (STARTED 2026-10-31, user)
User 10-31: "Deploy ff3.1, run validations (Peng, glpG, and other tests we ran)", then, asked which ff3.1: the retrain on
the fixed native target replaces the released merge. From now on ff3.1 is `ff31_nativefix` at step 76 (Phase 33 item
6: glycine settled, -ln(P_aL/P_aR) -0.545); the 10-26 merge is retired with a backup, as the merge retired its
predecessors (Phase 30). Scripts and logs `/project/trsosnic/yinhan/checks/ff31_nativefix_release_20261031/`.
- **Files:** `extract_ff.py` on step 76's checkpoint; rama.dat must be md5-equal to the run's `rama_round_04.dat`
  (`7693f10a`) and the other five files ff_2.1's values (every Adam group was held), so they are released as ff_2.1's
  files byte for byte, as the merge was.
- **Release** (release.sh of 10-26, adapted): the four trees' `parameters/ff_3.1` (ConDiv `$P`, the shared Upside 2
  deployment, the Upside 3 tree, `ff3_selection/params`), md5-verified; the merge kept as `.bak_pre_nativefix_20261031`
  beside each. Running Peng arms never re-read their parameters (bench_run.py builds the config once), but the new arms
  reuse the run directory names, so the merge's 5 running Peng chunks and its collector are cancelled first and its
  runs, scores and Upside 3 glpG directory move aside with the same suffix. The GPU tree (`upside3_gpu`, line 5's own
  session) is not touched; the Mac repo's `parameters/ff_3.1` is updated with a backup.
- **Validation, the merge's tests** (Phase 30's submit.sh): the 44-domain and mdCATH panels at dt 0.009 with TM-score
  against released ff2.1 and the control; Peng's 32 arms on the pre-audit Upside 2 engine and their collector; glpG on
  Upside 3 (the merge's `glpg_ff31/run.py`, one REMD block, TM4 read at T 0.70 against ff3.0's 0.85). Tags `p44_ff31nf`,
  `mdc_ff31nf`; Peng and glpG under `ff_3.1`.
1. [x] Extracted 16:45 10-31 (`extract.log`): rama.dat md5 `7693f10a` (= `rama_round_04.dat`); bb_env, hbond,
   environment and sheet md5-equal to ff_2.1's, sidechain.h5 equal to 2e-16 (`compare_sidechain.py`), so released as
   ff_2.1's files.
2. [x] 16:47-16:49 10-31: the merge's 5 Peng chunks and collector cancelled (`release.log`); `release.sh apply`
   (`release_apply.log`): four trees md5-verified; the merge's ff_3.1 in `<tree>/backup/` and
   `params_superseded/bak_pre_nativefix_20261031/`, its 32 Peng runs and 223 logs in `runs_superseded/`,
   `logs_superseded/`, scores `score/ff_3.1.bak_pre_nativefix_20261031`, glpG `glpg_ff31.bak_pre_nativefix_20261031`;
   Mac repo `parameters/ff_3.1/rama.dat` replaced (old `rama.dat.bak_pre_nativefix_20261031`).
3. [x] 16:50 10-31 (`submit.sh`, ids in `submitted.log`, remote_jobs.md §1 line 2).
4. [ ] Read each as it ends. Panels read 23:40 10-31 (findings 1.34): 44 domains fold 0.631 against released
   ff2.1's 0.603 (+0.028 [-0.004, +0.059]), mdCATH 0.515 against 0.488 (+0.027 [+0.000, +0.053]), level with the
   merge; mean TM level with ff2.1 on both (+0.001 and +0.004, mdCATH read 01:30 11-01). glpG keeps TM4 through
   preparation (stage 7.0 DSSP mean 0.95; merge 0.36); REMD block 1 (read 07:30 11-02, findings 3.15) holds TM4 at
   DSSP 0.74 at T 0.70, between ff3.0's 0.85 and the merge's 0.42, falling within the block (halves 0.78 / 0.70) at
   both ends of the helix. Open: Peng's 32 arms (2 scored; the 30 chains broken by midway2's fault restarted 12:18
   11-02, remote_jobs.md §1 line 2).

### Phase 38 - fix the all-term training's folding loss, then continue ff3.2 (STARTED 2026-10-31, user)
User 10-31: "Fix the folding issue and continue working on ff3.2". The loss is measured (findings 1.37 "What the loss
is"): the trained terms make every state away from the native cheaper, about half through the two scalars the DSE
target moves. Fix the cause in the training, never by a temperature or by freezing terms (Phase 33 item 9).
- **Resume here (the MacBook Pro's handoff 16:50 11-02; basin arm read 03:05 11-03 by the Mac Studio).** The next
  steps, in order; times, ids and paths in remote_jobs.md §1 line 2 and its "Resume here" step 4:
  1. [x] Epoch 1 (ended 16:47 11-02), read 17:40 by item 6: `z32_01` level with `z32_00`, still under released
     ff2.1, and the unfolded chain compacting further (item 6). Epoch 2's first link ran on midway2-[0316-0327],
     whose busy cores sit at 1.2 GHz; at the user's word (18:10) they are excluded and the link was moved after step
     39 (remote_jobs.md §1 line 2); at 31-43 min a step since, step 57 about 05:30 and step 76 about 17:00 11-03.
  2. [x] The basin arm (`ff32_basin`, step 19 at 00:48 11-03; `w32_00` and its readouts by 02:56), read 03:05 by item
     10's rule: (a) and (b) both met (item 10, findings 1.37 "The fix"). Reported to the user, who decides.
  3. [ ] Adopted by the user 11-03 ("Ok. Proceed", to the recommendation of 03:10): **ff3.2 is `ff32_basin` trained
     on** from step 19 to 76 under the restraint, one epoch a link (item 11; 57 steps, 30-41 h at the present 31-43
     min a step). `ff32_dse0` trains on to 76 (about 21:05 11-03) as the comparison at the same length and the
     fallback; item 7 applies to both.
  4. The chosen force field is judged by item 10 (c) on mdCATH, Peng and glpG's TM4 against released ff2.1 (about 3
     days of jobs).
- **Step 1, the causal test.** Two arms built exactly as `ff32_nativefix` (Phase 33 build.sh: its trainer copy, start,
  minibatch order, engine), each changing one thing, one epoch (19 steps), against its epoch 0 `n32_00` (0.449):
  `ff32_dse0` with the DSE target off (dse_weight 0), and `ff32_sarwsi` with the SARW replica as the SI defines it
  (backbone burial scaled with the other terms, A193's fix). Each judged by the 44-domain panel at its epoch end and by
  `peng_full.py`'s unfolded-state shift (n32_00 -11.9 E_up, released ff3.1 +1.2). Read: an arm that keeps ff2.1's
  folding (paired interval against released ff2.1 including 0) and has no unfolded shift names the cause.
- **Step 2, the fix**, from step 1's answer, decided with the user where it departs from the FF2 workflow; then ff3.2
  trains every term (depths included) under it.
1. [x] Built 16:52-16:56 10-31 (`/project/trsosnic/yinhan/checks/ff32_fix_20261031/`: `make_trainers.py`, `build.sh`,
   `build.log`, `build_sarwsi.log`): each trainer differs from md5 `54383004` by its one change; each initial
   checkpoint equals `ff32_nativefix`'s but for run paths and round-off at most 7e-15 (`compare_init.py`). Static
   check (`sarw_check/`): on a 3g0m SARW frame the backbone burial term is 0 in `ff32_sarwsi`'s SARW config and
   -77 E_up in the old one.
2. [x] Submitted 16:57 10-31 (`ff32_dse0` 49692528, `ff32_sarwsi` 49692529); at 19 each `after_training.sbatch`
   submits its panel (`z32_00`, `s32_00`) and then `peng_full.sbatch` on it.
3. [x] Read 14:10 11-01 (`paired_npz.py`, `peng_full_<tag>.txt`; findings 1.37): DSE off (`z32_00`) folds 0.575,
   +0.127 [+0.091, +0.163] over `n32_00`, -0.040 [-0.076, -0.007] under released ff2.1, unfolded shift +3.80 E_up;
   the SI's SARW (`s32_00`) 0.446, level with `n32_00`, shift -9.93. The DSE target carries about three quarters of
   the loss; the SARW's backbone burial none of it.
4. [x] Step 2 decided by the user 11-01 ("Proceed to continue work on ff3.2. Double check it is training all
   parameters, not holding anything fixed"): **ff3.2 is `ff32_dse0` continued**, every term trained, the DSE target
   off (a departure from FF2's dual target, the user's). The run directory keeps its name; tags `z32_EE`. The run
   that kept the DSE (`ff32_nativefix`, `n32_EE`) is stopped and stays the comparison.
5. [x] Parameter check (`checks/ff32_fix_20261031/param_moves.py`, `ckpt_rates.py`, outputs `params_check.txt`):
   every group the engine gives a gradient for trains and moved in epoch 0: the side-chain pair, coverage and
   hydrophobe tables (`rot`, 31,420 coefficients, rate 0.0125), environment centre, scale, sharpness and weights,
   the backbone burial scale, H-bond entries 0-3, all 14 sheet values, and 1,517 Ramachandran depth offsets (pooled
   mode; 2,955 of 5,040 basin offsets moved after round 1). Held because the engine returns no derivative, as in
   FF2: the backbone burial centre, sharpness and H-bond weight (commented out in `BackboneSigmoidCoupling`), the
   environment coverage geometry (`param_deriv` "not implemented"), H-bond entries 4-11, the hydrophobe and rotamer
   bead positions. Not parameters of the model's trained set: rotamer probabilities, the map shapes inside each basin,
   springs, backbone sterics. The basin offsets without 20 native sites (in the map or its residue) stay at NDRD.
6. [ ] Training on from step 19 since 14:43 11-01, one epoch a link, to step 76 (the end of training above; user
   11-02, which replaced the convergence gate and its extension to 114). At each epoch end the run's
   `after_training.sbatch` (md5 `11bd82e9` since 12:33 11-02; the gated version `.bak_pre_budget_20261102`, older
   `.bak_pre_epochs_20261101`) submits the panel `z32_EE` (dt 0.015) and `peng_full.sbatch` on it, writes
   `rg_check_ff32_dse0.txt` (`rg_check.py`: free-replica Rg against `ff31_nativefix`, which ran ff2.1's files on the
   same proteins in the same order; epoch 0 at T 1.1 -2.2%, with the DSE +3.9%), then trains one more epoch; at step
   76 it submits `z32_03` and stops. Its run under stubs gave both branches (step 38 trains on to 57, step 76 stops).
   It is the only submitter of this run's panels (`submit_new.sh` does not list it, so a pass cannot take the panel
   without its shift). Read each epoch: folding against released ff2.1, the unfolded shift, the Rg ratio (unfolded
   chains must not keep compacting), glycine's depths (`gly_depths.py` on the round file).
   **Epoch 1, read 17:40 11-02:** `z32_01` folds 0.578, level with `z32_00` (+0.002 [-0.033, +0.039]) and under
   released ff2.1 (-0.038 [-0.072, -0.005], p 0.038); unfolded shift +5.61 E_up (epoch 0 +3.80); glycine's
   alpha_L - alpha_R -0.791 (round 1 -0.877, released ff3.1 -0.781). The Rg check fails its condition: at T 1.1 the
   free replicas are 0.957 [0.907, 0.995] of epoch 0's Rg on the same proteins (smaller in 77%), 0.939 of
   `ff31_nativefix`'s (epoch 0 0.978); T 0.8 level (1.004). Reported to the user; the run trains on.
   **Epoch 2, read 17:50 11-03:** `z32_02` folds 0.599, level with released ff2.1 (-0.016 [-0.059, +0.026], p 0.33)
   and with `z32_00` (+0.024 [-0.016, +0.065]), under `w32_00` at p 0.088 (-0.031 [-0.069, +0.007]); unfolded shift
   +6.38 E_up; Rg at T 1.1 0.923 of `ff31_nativefix`'s (epoch 1 0.939); round 3's glycine alpha_L - alpha_R -0.758.
7. [ ] At step 76: `paired_npz.py` of `z32_00` to `z32_02` against `z32_03`, and of `w32_00` to `w32_02` against
   `w32_03` for ff3.2 since 11-03 (item 11); step 76 unless one is better at p < 0.05. Then the chosen epoch against
   released ff2.1: level, the full validation (mdCATH, Peng, glpG); below, ff3.2 is not deployed (Phase 33's rule).
   **Epoch 0, read 12:45 11-02** (user: "check ff3.2 basin depth difference between GLY alpha L and R, and whether
   it keeps protein folded"; `gly_depths.py` and `gly_LR_runs.py` on the Mac, which reproduce every earlier row). The
   pooled glycine map's alpha_L - alpha_R moved from the start's -1.064 (the merge) to -0.877 in round 1, next to
   ff3.0 (-0.836); released ff3.1 is -0.781 and the run with the DSE -0.804 after the same epoch and -0.675 at step
   95; -ln(P_aL / P_aR) -0.825 to -0.640 (ff3.0 -0.598). Shared glycine offsets aR +0.026 to -0.062, aL +0.422 to
   +0.515. The 44-domain panel (`panel.py select`, 37 domains with 100 folded frames in every tag): `z32_00` folds
   0.575 (released ff2.1 0.615, -0.040 [-0.076, -0.007]; with the DSE 0.449), and its native-basin retention
   against all-atom is level with released ff2.1 in every class: helical glycines -0.072 +- 0.051 (ff2.1 -0.091),
   left-handed glycines -0.016 +- 0.027 (+0.029), helix -0.019, beta -0.017. In training (`check_step.py`, free
   against restrained), helical glycines in alpha_R 0.836 over epoch 0 and 0.863 over epoch 1's first 15 steps, on
   round 1's map; left-handed glycines in alpha_L 0.903 and 0.900; free CA-RMSD of matched minibatches unchanged
   (2.75 and 2.75 A, 2.14 and 2.30 A).
- **Step 3, the remaining loss** (user 11-02: "find the source of the remaining loss and try to resolve it"): `z32_00`
  folds -0.040 [-0.076, -0.007] under released ff2.1 with the DSE off, while its unfolded shift is +3.80 E_up (unfolded
  states raised against native), so the remaining loss is not the DSE's global one. Diagnosed as the DSE's share was;
  the fix is decided with the user where it departs from the FF2 workflow, and never by freezing terms or a temperature.
  Scripts, logs and outputs `/project/trsosnic/yinhan/checks/ff32_rest_20261102/` (`submitted.log`).
8. [x] Where and through which terms, read 13:25 11-02 (findings 1.37 "The remaining loss"): tertiary, through a
   directed change of the side-chain coverage and burial tables that favours the crystal over the thermal native
   and, as a side effect, compact non-native states that no replica of the target presents. `diagnose.sbatch`
   (`fail_mode.py`, `decomp.py`, `decomp2.py`), `decomp3.py` (direction or noise for the hydrophobe coverage and
   the burial; new), `crystal_view.py`.
9. [x] Causal, read 14:20 11-02 (`submit_resets.sh`; `resets_z32_00.txt`; findings 1.37): the side-chain file
   (+0.054 over `z32_00`) and the burial (+0.041) each bring it level with released ff2.1, the backbone burial's reset
   costs -0.096, the H-bond, sheet and library resets change nothing resolved: the signs item 8 predicted.
10. [ ] **The fix: the native-basin restraint** (user 11-02, chosen over persistent free replicas and no change). The
   native-restrained replica is held in the native basin as the panel defines it, not at the crystal: a one-sided
   wall (the engine's pair `WallSpring`, `advanced_config` `pair_wall`) at 1.2 d0 on every native CA contact
   (|i-j| >= 4, d0 < 8 A; `panel.py` `native_contacts`, `fraction_native`), zero inside and the old restraint's
   1/3^2 outside (both are 0.5 k x^2), replacing the random backbone elastic network at 1/3^2. The 20 deg phi/psi
   springs stay (a folded residue's spread; they keep the glycine target and the chirality the distances cannot).
   Inside the basin the replica samples the free Hamiltonian, so the NSE becomes native basin against free
   ensemble, the gradient of log P(folded). Departs from FF2's restraint.
   - [x] Trainer `ConDiv_basin.py` (md5 `0213050c`) = `ConDiv_dse0.py` (`a86a9062`) with that one change
     (`make_basin.py`, asserted replacements; the wall written with `upside_nodes`' `DistanceCoord` and
     `WallSpringPotential`).
   - [x] Smoke test, 13:52 11-02 (`smoke_basin.py`, training protein 3g0m, 138 residues, 263 native contacts): the
     configs differ only by the restraint nodes; the wall is 0 at the native and equals the analytic energy on frames
     to 4e-5; at T0 over 3000 tu (second half) the old restraint holds RMSD 0.92 A (0.79-1.05), Q 0.983, the wall
     1.33 A (1.10-1.56), Q 0.961 (min 0.909): out of the pin, inside the basin (the free T0 replica sits near 2 A).
   - [x] Run `ff32_basin` (tag `w32_00`) built 13:53 11-02 by `build_basin.sh` (initial checkpoint equal to
     `ff32_dse0`'s but for paths, largest float difference 8.6e-15; minibatches equal; round-0 library `5b66c78a`);
     one epoch (19 steps), first link 49746562; at 19 its `after_training.sbatch` submits panel `w32_00`,
     `peng_full` and the Rg check.
   - [x] Read against `z32_00` 03:05 11-03, by the rule below (findings 1.37 "The fix"): (a) met, the restrained
     replica at 1.23-1.53 A (old 0.85-0.91, free 2.1-3.0), the crystal's advantage -1.76 E_up and the loose frames'
     -0.40 on the 32 domains `crystal_view.py` and `decomp2.py` keep (`z32_00` on them -3.69 and -1.13), 3g0m's 3-7
     A dip gone; (b) met, folded 0.630 against 0.575 (+0.055 [+0.022, +0.091], p 0.004; level with released ff2.1,
     +0.015), mean TM 0.684 against 0.679 (+0.006 [-0.002, +0.013]). Unfolded shift +6.99 E_up (`z32_00` +3.80); Rg
     at T 1.1 0.943 of `ff31_nativefix`'s (`z32_00` 0.978, `z32_01` 0.939). The user decides.
   - **Decision rule, fixed 14:25 11-02 before any result of the arm** (user 11-02: "double check we are not forcing
     ff3.2 to produce the result we want. all edits should make sense, not just to twist parameter to generate
     results we want"). (a) Mechanism: the arm's change must do what the reason for it says, independent of folding:
     the restrained replica leaves the crystal pin and stays in the basin, the crystal's advantage over native-like
     frames (`crystal_view.py`; `z32_00` -3.79 E_up) and the loose frames' (`decomp2.py`; -1.08) both shrink. (b)
     Folding on the 44 not lower than `z32_00` in folded fraction and in mean TM. Only (a) and (b) together carry the
     arm on, as the user decides. No variant of the wall (its 1.2, its 8 A, its stiffness) is tried on any outcome:
     those are the panel's own definitions and the old restraint's stiffness, set before the arm, and a failure is
     recorded, not tuned. (c) Whatever ff3.2 becomes is judged on mdCATH TM, Peng's TM and glpG's TM4 against released
     ff2.1: proteins and metrics never used to choose. The 44-domain panel has chosen the DSE's removal, found this
     loss, and will choose this arm and the epoch, so its numbers for ff3.2 are optimistic and are not evidence of
     quality by themselves.
11. [ ] **ff3.2 under the restraint** (user 11-03; Resume here step 3): `ff32_basin` from its step-19 checkpoint to
    76, one epoch a link, nothing else changed (driver `0213050c`, the same minibatch order and solver). Its
    `after_training.sbatch` becomes `ff32_dse0`'s per-epoch form with tags `w32_EE` (the one-epoch version kept as
    `.bak_pre_epochs_20261103`), and adds `mech.sbatch` (`decomp2.py`, `crystal_view.py` on the epoch's panel,
    afterok it). Checked by `--test-only` and a stub run of both branches; no test job. Read each epoch as item 6
    (folding against released ff2.1 and `w32_00`, the unfolded shift, Rg, glycine's depths) and the crystal and
    loose gaps, which show whether the hot chain's compaction (Rg 0.943 at epoch 0) goes into the compact
    non-native states the restraint removed. Training keeps no frames (each step's divergence file holds the
    replicas' mean Rg and RMSD only), so the hot replicas' native contacts cannot be read from it. If the gaps grow
    back or folding falls, that is recorded and reported; the wall is not changed.
    **Step 31, read 17:45 11-03:** 30 min a step, 24 of 24 proteins, KE/1.5kT at most 1.015. The H-bond margin has
    sat at +0.098 to +0.106 since step 19, under the +0.10 line at steps 24-26 and 31 (`ff32_dse0` +0.131 to +0.159
    over epoch 3); the user was notified and the run trains on.
    **Epochs 0-2 read 15:30 11-04** (MacBook Pro; user: "check ff3.2 for (1) basin depth difference between alpha L
    and R in GLY, (2) whether it holds protein folded, as this failed in last version"). Step 70 at 15:08, 30 min a
    step, 76 about 18:10. (1) Pooled glycine alpha_L - alpha_R (`gly_depths.py`, dse0 round 1 reproducing -0.877):
    rounds 0-3 -1.064 / -0.915 / -0.829 / -0.774 (per map -ln(P_aL/P_aR) -0.825 / -0.677 / -0.591 / -0.537), round k
    the library of epoch k, so `w32_02` runs round 3; released ff3.1 -0.781, ff3.0 -0.836, ff2.1 -1.460, the previous
    ff3.2 -0.675. Shared glycine offsets after round 3 alpha_R -0.095, alpha_L +0.580 (ff3.0's pair -0.107, +0.517).
    Still moving, about 0.06 a round. (2) Folding: the panel folds `w32_00` / `w32_01` / `w32_02` 0.630 / 0.633 / 0.674
    against released ff2.1's 0.615 (+0.015 [-0.015, +0.045], +0.018 [-0.016, +0.051], +0.059 [+0.015, +0.103] p
    0.023; mean Q +0.014 p 0.029 at epoch 2), where the previous ff3.2 folded 0.447 (-0.169); unfolded-state shift
    +6.99 / +6.35 / +7.94 E_up, 0 of 16 negative (previous -11.94, 16 of 16); states beyond 15 A +19.1 / +23.0 / +28.3
    E_up above native (`crystal_view`). Against that, two trends: the loose-minus-native-like gap against released
    ff2.1 (`decomp2`) -0.40 / -0.93 / -1.20 E_up over 31-32 domains, through hydrophobe H-bond coverage and rotamer
    (backbone environment coupling +0.35 the other way), and the hot chain's Rg at T 1.1 0.943 / 0.896 / 0.880 of
    `ff31_nativefix`'s (T 0.8 level, 0.996-0.997). Folding has not fallen; by item 11 the growing gap is reported, and
    `w32_03` (about 19:30) shows whether folding follows it.
    **Why, read 16:40 11-04** (findings 1.37 "Why the two diagnostics grow"): both are states the contrast never
    contains (free replicas restart at the native each step and the free phase is the three coldest): the loose gap
    is the 11-02 side-chain packing mechanism regrowing (directed, half uniform hydrophobe coverage), the compaction
    is backbone burial favouring collapsed unfolded chains (-1.95 of -2.23 E_up at `w32_02`), which the DSE target
    held the other way (`n32_04` +3.08). Fixes, each a departure from FF2's workflow for the user to choose and to be
    tested as a one-epoch arm from step 19 by a rule fixed before its data: a denatured-state term with a set point
    (the unfolded ensemble held at released ff2.1's size instead of pushed toward the SARW), or non-native states in
    the free phase (persistent chains across a protein's visits, or warmer replicas in the free ensemble); or no
    change, with `w32_03`'s read and the epoch choice deciding.
    **Chosen 11-04 (user: "proceed with your suggestion, but make sure the mistakes can be corrected later"): no
    change.** `w32_03` and the epoch choice are read by a detached job (`checks/ff32_epoch3_20261104/`, job 49816524 ->
    `read.sbatch` -> `summary.txt`) by item 7's rule written into `choose.py` before the data (dry-run on `w32_02`):
    step 76 unless an earlier epoch folds better at p < 0.05; the chosen epoch below released ff2.1 (95% interval
    under 0) is not deployed, otherwise its full validation follows. If `w32_03` folds below `w32_02`, the epoch choice
    takes `w32_02` and the fixes above are tested as arms. To keep every step reversible: ff3.2's and `ff32_dse0`'s
    epoch-end checkpoints, their round libraries and the trainer are archived read-only with md5s
    (`checks/ff32_checkpoints_20261104/`, README with the arm recipe); nothing is released into `parameters/` before
    the validation passes; `ff31_nativefix`, `ff311_glypair`, `p44_ff21_released` and every `w32_*` / `z32_*` panel stay.
    **What follows each verdict, decided 11-04** (user: "decide for me what needs to be run"):
    - Validation: item 10 (c) as `checks/ff31_nativefix_release_20261031/submit.sh` ran it for ff3.1: the mdCATH panel
      and its TM against `mdc_ff21_released`, Peng's 32 arms and their collector, glpG on Upside 3 on the CPU
      (preparation and one REMD block, TM4 DSSP), like-for-like with ff2.1's, ff3.0's and ff3.1's glpG runs. The
      candidate runs under its own name and replaces no released force field. Its ubiquitin arms' replicas are kept
      for H2 (`upside3.md` 8.2). Script `/project/trsosnic/yinhan/checks/ff32_validate_20261104/submit.sh` (README;
      it takes the epoch from `summary.txt`, refuses a BELOW verdict, stages `ff_3.2_<epoch>`; dry run passed 11-04).
    - Fix arms, if the chosen epoch is below released ff2.1 or `w32_03` folds below `w32_02`: first persistent free
      chains across a protein's visits. The contrast then contains the loose and compact non-native states both
      diagnostics come from, and the objective stays the free ensemble at T0. The set-point term would impose
      released ff2.1's unfolded size as a target, and warmer free replicas change the objective. Built as a one-epoch
      arm from step 19, its rule fixed before its data.
    **Verdict, 19:13 11-04 (`checks/ff32_epoch3_20261104/summary.txt`): `w32_03` chosen, level with released ff2.1;
    its validation runs.** Rule 1: no earlier epoch folds better at p < 0.05 (`w32_02` +0.029 [-0.002, +0.061] p
    0.073 over `w32_03`, so not below by the rule). Rule 2: `w32_03` 0.645 against released ff2.1's 0.615, +0.029
    [-0.008, +0.066] p 0.065, LEVEL. Beside, deciding nothing, epoch 3 continues the trend of epochs 0-2: free-replica
    Rg at T 1.1 0.854 of ff3.1's (0.943, 0.896, 0.880 before), loose-minus-native-like gap -1.717 (32 of 32
    negative), states beyond 15 A +32.8 above native, unfolded shift +8.71, glycine L - R -0.770 (round 3 -0.774
    reproduced). The validation (`checks/ff32_validate_20261104/submit.sh apply`, 19:26 11-04, by the watch) runs
    `ff_3.2_w32_03`, staged beside the released force fields, not released; jobs in `remote_jobs.md` line 2. No fix
    arm runs. Its mdCATH panel `mdc_w32_03` (dt 0.009, 356 runs, COMPLETED 01:37 11-05, read 01:55 with
    `paired_npz.py`): 89 domains fold 0.566 against released ff2.1's 0.488 (+0.078 [+0.049, +0.107] p 9.7e-07, higher
    on 64, lower on 20) and ff3.1's 0.515 (+0.051 [+0.023, +0.080] p 5.6e-05); mean Q +0.017 and +0.014. Its TM files
    (49823717, 02:22 11-05, `checks/ff32_validate_20261104/tm_mdcath_*`): mean TM 0.739 against released ff2.1's
    0.727 (+0.012 [+0.006, +0.018] p 1.1e-04, higher on 62 of 89) and ff3.1's 0.731 (+0.008 [-0.000, +0.016] p
    0.023, higher on 53); TM > 0.5 level with both (0.967 against 0.972 and 0.974).
    **Peng interim, 12:20 11-06** (`checks/ff32_interim_20261106/peng_interim.txt`, job 49882604; the scorer's own
    method on every frame each ff3.2 arm had written in finished chunks, against released ff2.1's first frames up to
    the same count; it reproduces both BBA native jsons to every digit; not the verdict, which is `peng_ff32.txt`):
    level over the 32 arms, mean dTM -0.004 [-0.017, +0.009] (p 0.31, 12 higher, 20 lower; native -0.005, de novo
    -0.003). The losses sit in BBA (-0.056 native, -0.032 de novo), the WW domain (-0.065, -0.046) and cspA de novo
    (-0.079); the gains in ubiquitin (+0.077, +0.075) and NuG2 de novo (+0.074). BBA's mechanism is findings 1.39:
    the backbone burial scale fell from -0.301 to -0.136 over `w32_00`-`03` with nothing in the contrast to hold it.
12. [ ] **ff3.2 with persistent free chains** (user 11-06: "fix the ff3.2, train it and re-do the validations"; rule
    fixed 12:45 11-06, before any run). Why: every free replica restarts at the native each visit and the contrast
    takes its second half (4000 of 8000 time units at T 0.80-0.85), so the free ensemble never holds the non-native
    states the trained terms make cheaper, and nothing in the objective pushes them back up (11-04 "Why the two
    diagnostics grow"). Three readings grew every epoch: the hot chain's compaction, the loose-minus-native-like
    gap, and the backbone burial scale's fall. In Peng's long runs BBA turns into one helix and the WW domain loses
    TM. The fix chosen 11-04 for this case, first of its list: persistent contrastive divergence (Tieleman 2008).
    The objective stays the free ensemble at T0; only its sampling changes.
    - The change (`checks/ff32_pcd_20261106/make_pcd.py`, asserted replacements on `ConDiv_basin.py` `0213050c` ->
      `ConDiv_pcd.py` `abb27c42`): each protein's 12 free replica slots start from their own final frames of the
      protein's previous visit, one epoch earlier at the same minibatch (the minibatches are fixed and walked in
      order), written to that step's directory as `<code>.free_final.npy`. Epoch 0 starts at the native, as before,
      and a missing or mis-shaped file raises. The restrained and SARW replicas start at the native, as before.
      Nothing else changes: the basin wall, DSE off, the minibatch order, the solver, 76 steps.
    - The run `ff32_pcd` (tags `p32_EE`), built by `build_pcd.sh` exactly as `ff32_basin` was, from the same
      initial checkpoint (compared by `compare_init.py`), trained 0-76 one epoch a link. Epoch 0 has no previous
      visit, so it repeats ff32_basin's protocol at new random seeds, and the persistence acts in epochs 1-3. Its
      `after_training.sbatch` is ff32_basin's with the names changed (panel, unfolded shift, gaps, Rg check, next
      link). About 38 h of training at 30 min a step.
    - Checks before any run, no test job: the driver differs from `0213050c` only by the edits; a static test on
      the login node (`static_test.py`: the previous-visit paths, frames written to a copy of a real config read
      back equal and the engine's energy finite on them); the initial checkpoint equal to ff32_basin's but for
      paths; the after-training script's stub run of both branches; `--test-only` of the first link.
    - Read each epoch as item 11 (folding on the 44 against released ff2.1 and against `w32_EE` at the same epoch,
      the unfolded shift, the gaps, glycine's depths) and the backbone burial scale. The Rg check and the free
      replica's RMSD now read persisted chains, so they are reported but no longer comparable with ff32_basin's.
      Decision rule, fixed now: (a) mechanism, the persisted chains reach non-native states: in epoch 1 the median
      over proteins of free replica 1's mean CA-RMSD (`divergence.pkl` `rmsd`) above epoch 0's; without it the
      change did nothing and is reported as such. (b) Folding on the 44 at step 76 not below `w32_03` (paired 95%
      interval not wholly under 0). The epoch is chosen by item 7's rule (step 76 unless an earlier epoch folds
      better at p < 0.05), and a chosen epoch below released ff2.1 is not deployed. (c) The chosen epoch runs item
      11's full validation (mdCATH panel and TM, Peng's 32 arms and collector, glpG on Upside 3), against released
      ff2.1, with `w32_03` beside it. No variant of the change is tried on any outcome. The backbone burial scale
      has no target value; it is read, never set.
    - **Revised 13:25 11-06, before any of its data was read: 57 proteins a step, 80 steps** (Jumper's suggestion, more
      proteins in each round; the user chose it from the measurement, findings 1.40: ff32_basin's gradient noise scale
      is 127-2,221 proteins by group in epoch 0 and beyond the 456 after, so a step of 24 is 1-16% signal). The run
      at 24 (`ff32_pcd`, 49883975) was stopped after 1 step and is not read. `ff32_pcd57` (tags `p57_EE`):
      `ConDiv_pcd.py` with `minibatch_size` 57 (`make_pcd57.py` -> `ConDiv_pcd57.py`), 8 steps an epoch, 80 steps =
      10 epochs, Adam and every other setting unchanged; the depth offsets move once an epoch, so 10 rounds. Its
      minibatches are ff32_basin's shuffle of the same 456 proteins re-sliced into 8 (`match_order.py`). Each link
      runs 16 steps (`train_chain57.sbatch`: 57 tasks of 14 cores, 798 CPUs), and its `after_training.sbatch` reads
      epochs 1, 3, 5, 7 and 9 (steps 16-80). The persistent chains get 9 revisits. About 2.5 times the 76-step run's
      CPU (4,560 protein visits against 1,824). The rule above stands with these numbers: (a) at epoch 1 against
      epoch 0; (b) folding at step 80 not below `w32_03`; the epoch chosen by item 7's rule among the five read.

### Phase 39 - the most-probable-state side chain: BP picks one state per residue (STARTED 2026-11-03, user)
User 11-03, the PI's request: after BP gives each residue's state probabilities (e.g. 50 / 30 / 20%), keep only the
most probable state and ignore the others, then retrain. The user delegated the design ("decide them for me, and go
ahead and implement and deploy it") and ruled that **the engine change must not be git-committable**: it lives only
in `scratchpad/bpmax_20261103/` (git-ignored) and in its own cluster tree, never in this repo's `src/` or `py/`.
- **What it is.** BP runs as now. Each residue takes s* = argmax of its node marginal (ties to the lowest index);
  the energy is E = sum_i E_i(s_i*) + sum_ij E_ij(s_i*, s_j*) from the raw 1-body and pair values (no entropy, no
  mutual information), and every sensitivity is 1 for the chosen states and 0 otherwise, so forces and ConDiv's
  parameter derivatives (`ConDiv.py` `get_param_deriv`) are those of E at fixed states. A NaN marginal gives NaN
  weights, so it reaches the energy and forces as in BP mode. Differs from `ff21_nobp_dt009` (one fixed state per
  residue type, no BP, 1-body 0): here the state follows the environment step by step.
- **Decisions (delegated, 11-03).** Per-residue argmax, as the PI described (not the joint MAP of max-product BP).
  The ff2.1 retraining protocol of `ff21_ctrl_dt009` / `ff21_nobp_dt009` (init ff2.1, its minibatch order, dt 0.009,
  same trainer and rates) to step 76 and the same convergence gate; compared at the fixed points with the BP control
  and the fixed-state arm on the 89 mdCATH domains and the 44-domain panel, by TM-score. Run name `ff21_bpmax_dt009`,
  "the most-probable-state side chain".
- **The earlier arms' engine history** (user 11-03: "the program should be the same state as the other tests"). Both
  ran steps 0-37 on `$P`'s binary of 10-02 (before the angular-gradient fix, findings 3.12) and 38-76 on the fixed one
  of 10-18; the fix changes only the side-chain pair angular parameter gradient, not the simulations. This arm does the
  same: steps 0-37 on `obj_preangfix`, 38 on `obj_angfix`, both with the mode; the tree's `obj` symlink is switched
  at step 38 by the run's `setup/switch_engine.sh`. A first link on the fixed engine from step 0 was cancelled in step
  0 for this.
- **Where.** The patch is against `$P/src/rotamer.cpp` (the Upside 2 engine both earlier arms trained on), not Upside
  3's. A dedicated tree `/project/trsosnic/yinhan/upside2-md-bpmax` holds the patched sources, its two builds, `$P`'s `py/`
  with `write_rotamer` setting `most_probable_state = 1`, and the training scripts, so every config written there
  (training and its panels) carries the mode and runs on the patched binary. `$P` is not touched.
- **Physics risk, measured before training.** E jumps whenever a residue's argmax changes, and forces never see the
  jump, so the dynamics conserves no energy and need not sample exp(-E/kT). Measure on the patched build, BP against
  max mode on the same engine: argmax switches per step, the size of the jumps, KE/1.5kT, NVE drift. Large effects
  are reported before any training.
1. [x] Patch (`scratchpad/bpmax_20261103/rotamer_bpmax.patch`, 83 lines, `rotamer.cpp` only), stock and patched
   builds of `$P`'s source on the Mac (`build_stock`, `build_max`).
2. [x] Checks on the Mac (findings 1.38; scripts and results in `scratchpad/bpmax_20261103/tests`, `res`): BP mode
   bitwise equal to the stock build on 3 training proteins; max mode's energy equals the chosen states' energy
   rebuilt independently, its forces and parameter derivatives agree with central differences; the dynamics
   measurement: about one switch every second 3-step round, KE/1.5kT 0.8-2.4% above BP's under the thermostat over
   T 0.80-1.10 (at most 1.025), NVE drift 8-40 times BP's. Judged a property of the model, not a defect: proceed and
   report.
3. [x] Cluster tree `/project/trsosnic/yinhan/upside2-md-bpmax` built on midway2 with `$P/obj`'s compiler and flags,
   twice: each build in BP mode is bitwise the binary the earlier arms ran at the same steps (`obj_preangfix` against
   `upside.bak_pre_angfix_20261018`, `obj_angfix` against `$P/obj`) in energy, forces, parameter derivatives and a
   30-tu trajectory (1aba, 1a62); in max mode the two differ only in the rotamer angular parameter derivatives.
4. [x] Run directory as `ff21_nobp_dt009` was made: init_param md5-equal to released `ff_2.1`, every input, the
   trainer, the Python and the venv the earlier arms'; gate, engine switch and mdCATH panel in the new tree. Setup
   49766481 done; link 49766810 to step 38 running from 12:42 11-03, on `obj_preangfix` (remote_jobs.md §1 line 1).
4b. [x] Step 38 reached 10:20 11-04: `obj` switched to `obj_angfix` (`setup/engine_history.log`, `gate_49806985.out`)
   and link 49806998 ran to 76 with its workers on `obj_angfix/upside`.
4c. [x] Step 76 reached 07:08 11-05 (24 of 24 proteins every step, no error line, step 76's KE/1.5kT at most 1.031).
   The gate (`gate_49846050.out`, steps 58-76, 10 groups, pass if every p > 0.005) CONVERGED: smallest p 0.033
   (`dhb`, size 1), then 0.742 (`bbenve`); its mdCATH panel `mdc_bpmax_03` on epoch_03_minibatch_18 queued itself.
5. [x] At the fixed point: the mdCATH panel (automatic) and the 44-domain panel at dt 0.009 (in the new tree),
   TM-score against the fixed points of the BP control and the fixed-state arm, with released ff2.1 beside. Script
   `/project/trsosnic/yinhan/checks/bpmax_fixedpoint_20261104/submit.sh`: `plan` passed `--test-only` on all three
   jobs and the watch ran `apply` at 07:49 11-05 (its `submitted.log`); the panels and TM jobs are in remote_jobs.md §1
   line 1.
   - **44-domain panel `p44_bpmax_03`, read 09:50 11-05** (176 runs, KE/1.5kT 1.009-1.035, no non-finite output; files
     `tm_p44_p44_bpmax_03_vs_*.txt` in the check directory; mean TM over frames after 2000, 4 seeds pooled, Wilcoxon):
     level with the BP control `p44_ctrl_03`, 0.652 against 0.650 (+0.003 [-0.006, +0.011] p 0.51; higher on 25 of
     44); above the fixed-state arm `p44_nobp_03`'s 0.623 (+0.029 [+0.012, +0.046] p 0.0023; higher on 32); below
     released ff2.1's 0.681 (-0.028 [-0.041, -0.018] p 3.6e-06; lower on 35; TM > 0.5 0.839 against 0.886, p 0.0065),
     where the BP control also reads 0.650.
   - **mdCATH panel `mdc_bpmax_03`, read 14:45 11-05** (356 runs, KE/1.5kT 1.009-1.032, no non-finite output; files
     `tm_mdcath_mdc_bpmax_03_vs_*.txt`): mean TM 0.696, level with the BP control `mdc_ctrl_03`'s 0.702 (-0.006
     [-0.013, +0.000] p 0.25; lower on 45 of 89, higher on 44); above the fixed-state arm `mdc_nobp_03`'s 0.667 (+0.028
     [+0.017, +0.039] p 1.5e-06; higher on 65); below released ff2.1's 0.727 (-0.032 [-0.040, -0.023] p 7.1e-10; lower
     on 69). TM > 0.5 0.937 against 0.946, 0.920 and 0.972.
   - **Result, both panels:** keeping only BP's most probable state folds as well as full BP and better than one
     fixed state per residue type, so the probability weighting itself buys nothing measurable here, while letting the
     state follow the environment does. All three retrained arms fold below released ff2.1 by the same margin, so that
     gap belongs to the retraining, not to the side-chain treatment.

## Known Errors / Blockers

* **Upside 3 protein runs do not conserve energy without the thermostat** (11-02, Phase 35 item 11's V2 drift;
  findings 11g): 1ubq at T 0.85 with ff2.1 and with ff3.1, 4 copies x 5000 tu, gains about 90 E_up of total energy
  per 1000 tu per copy (KE/1.5kT 1.12-1.14 averaged over the run) on the CPU and on the V100 alike. Production runs
  carry the thermostat, which removes the heat. Not diagnosed; candidates are BP's early stop at tol 1e-3 (forces
  not the gradient of the energy), cutoff jumps such as the H-bond's (findings 3.14 below A239), and dt 0.009. A
  first test would repeat `device_dyn.py nve` on the CPU at BP tol 1e-6 and at a smaller dt.

* **Upside 3's glpG hybrid with ff2.1 loses TM4 in preparation stage 7.0** (findings 3.15, 10-26), answered 10-27:
  the loss is ff2.1's own, not Upside 3's. The same preparation and seeds, read at the end of stage 7.0: ff3.0 (A235,
  `upside3_validation/glpg_ff30`), whose placement-fixed Upside 2 chains hold TM4, holds it on Upside 3 too (TM4 DSSP
  alpha 0.88 at the last frame, last tenth 0.82); ff3.1 (A234, `glpg_ff31`), ff2.1 with only its map trained, loses it
  (0.59, last tenth 0.36) as ff2.1 does (0.24, 0.41). ff3.0's REMD block 1 at T 0.70 confirms it against Upside 2's
  0.96. Which force field the Upside 3 hybrid HDX uses is the user's call; no Upside 3 hybrid HDX result is read before.

* **Every pre-audit hybrid dry-MARTINI result carries these defects** (audit 10-20; fixed in the repo, findings 3.14
  A70-A76 and A83; Phase 24). The cluster trees running the campaign are unchanged.
  * The side-chain/lipid table treated the choline bead NC3 (+1) as neutral (charge from the bead type): ASP-NC3
    misses up to -2.15 E_up at 5 A, and LYS-NC3 at 6 A has the wrong sign. It also read its angle mirrored in the
    engine and placed the CB 0.57 A off (A76).
  * The NPT barostat never scaled a lipid, and the NPT stages total 4.6 ps against CHARMM-GUI's 10 ns, so
    production starts at an unrelaxed box.
  * The rotamer 1-body is not ff2.1's: the prep wrote -log of the unweighted phi/psi-grid average instead of the
    phi/psi-dependent term (relative rotamer energies in a helix off by up to 9.1 E_up, LEU). The live glpG REMD
    seeds (`popepopg_REMD_mdw2/seeds/`) carry `placement_fixed_scalar`, while master's ff2.1 workflows use
    `--dynamic-rotamer-1body`.
  * CHARMM-GUI's bead overlaps were never removed and the PO4 restraint's reference sat 53 A off in z, so only the
    hard z freeze held the heads (local preps heated to KE/1.5kT 26-33 in stage 6.3; the cluster glpG preps did not
    overheat but carry the offset). `scratchpad/relax_popepopg.sh` calls the deleted `prepare --mode bilayer`.
  * The ion count neutralised residue charges the protein does not carry: net +3 e (1rkl), -6 e (glpG).
  * Every ALA and GLY side chain was missing from the protein's side-chain node, its placement and the H-bond
    coverage (A83): glpG 42 of 210 residues, 1rkl 3 of 31.
  `parameters/ff_2.1/martini.h5` was rebuilt 10-20 (md5 `cc90a933...`; the old `2e612261...` kept as
  `martini.h5.bak_20261020`; ff_3.0 configs use it too) and again 10-22 (Phase 24 item 5.6). Every prepared hybrid
  `.up` must be re-prepared; the patched engine refuses the old table, and since A103 also any stage file without
  `combined_energy_end_slopes`.
* **Hybrid HDX lipid accessibility from the pre-audit `martini_hdx_membrane_accessibility.py` is shifted by one
  residue** (findings 3.14 A61, fixed in the repo 10-20): glpG, numbered from 1, gave each amide its N-terminal
  neighbour's accessibility. The shared deployment's copy (`/beagle3/trsosnic/yinhan/upside2-md/py`, md5 `ca115325`)
  is that version, so every hybrid HDX run through it carries the shift whatever its date: the 10-14 and 10-17 TM4
  readouts and the 10-22 ff3.0 block figures too (findings 5.6). The 10-22 ff3.0 sets (blocks 2, 3 and 2-3, four
  variants) were regenerated with the repo's code 11-01, then again with the core test
  (`checks/hdx_acc_flicker_20261101/a61/` and `core/`, findings 5.8); their work dirs still hold the shifted arrays.
  To regenerate the rest (on the cluster, by the campaign's drivers):
  * `/project/trsosnic/yinhan/popepopg_REMD/<V>/hdx*/results/<V>_popg_remd_<i>_ACC.npy`, and the same arrays in
    `popepopg_REMD_mdw2/<V>.<ff>*/hdx_sc600x12_*/results/` (for the 10-22 set, `acc_shift_corrected/` undoes the
    shift from the stored arrays for every amide whose next residue is a donor);
  * the combined `_PS.npy` (`_PS_protein.npy` stays valid);
  * every calc_hdx_ht output in those directories (`_dG_profiles.npz`, `.csv`, `_ranked.csv`, `_<T>.csv`, PNGs);
  * anything `scratchpad/local_popg_79HIS/acc.sh` produced; findings 5.3's `acc` values.
  Drivers that call the deleted `example/00` helpers need repointing: `hdx_cluster.sbatch` and three
  `scratchpad/np_slurm/*.sh` (`get_info` is now `py/get_info_from_upside_traj.py`; the `write_hybrid_energy` step is
  dropped, since the projection copies the coupled potential and the MBAR consumers subtract the mean). The
  `scratchpad/local_popg_79HIS/{hdx,acc}.sh` drivers are retired, not repointed: they call the deleted
  `write_hybrid_energy.py`, and their trajectories are gone.
* **Every ConDiv run before 2026-10-18 08:35 trained the side-chain pair term's off-diagonal angular coefficients on
  a mixed gradient** (findings 3.12; fixed in the engine then): ff2.1's retrainings to step 37, ff3.0, ff3.1, ff3.2.
  The corrected gradient shows no pull at ff3.0's or the control's fixed point (Phase 25); ff3.0 and ff2.1 are not
  retrained.
* **Every hybrid glpG result before 2026-10-17 ran with misplaced side-chain beads** (findings 4.17; Phase 18):
  TM1/TM4 helicity, glycine flips and HDX, under every force field, on the cluster and locally.
* **Older glpG TM4 sets.** Every count before 2026-10-07 11:22 is invalid (findings 3.11): those runs used the
  FF1-form ff_3.0's coverage tables (every local set of Phases 10-11 until then, the 09-30 ff_3.0 validation of
  Phases 3-5). Every hybrid run before 2026-10-10 08:39 used the under-sampled side-chain/lipid table (findings 4.10,
  4.15): their comparisons with each other stand, their absolute TM1 and TM4 values do not, and they are not compared
  with `_sc` sets. `patch_glpg.py` carries both fixes (md5 `97ca3294...`); another computer must update its local copy
  before use (remote_jobs.md "Setting up the next computer").
* **Total-potential jumps in local glpG runs** (findings 1.22): 13 of 105 pre-fix seeds, one blow-up. The likely
  cause is the old table's kicks (findings 4.10; none in the 24 fixed-table seeds of 4.13 and 4.14), unproven for the
  jumps themselves. Check the `_sc` sets for them; do not drop frames. One-frame backbone excursions below the 3000
  cut sit where the helices fail (findings 4.9).
* **The TM1 term is not adopted** (user, 10-10 11:20; findings 4.14). Its f was trained in Upside's implicit membrane
  and it is anchored at z = 0; any later test of it watches the protein's depth against the explicit bilayer
  (`analyze_hbmem.py` health block).
* **midway2-login2's `/project` hangs** (since 10-07 10:48). Route `/project` work through login1 (remote_jobs.md §8).
* **The cluster's `training/` keeps the pre-10-02 layout.** Do not sync the repo's `training/` over it.
* **`run_remd.py` rolls back NaN replicas and continues**, which the NO GUARDS rule forbids. It is kept unchanged so
  glpG stays comparable across campaigns. Whether to remove it is the user's call.
* **CPU work runs on midway2 broadwl** unless the user names midway3 caslake. No `amd`, `beagle3` or GPU partitions,
  and not `broadwl-lc` (its nodes cannot see `/project`). GPU work (Phase 35) runs on midway2 `gpu2` and midway3
  `gpu`, never `beagle3` (user 10-29).
* **HDX post-processing items await the user's decision** (user, 10-11; findings 5.5): three docstring/README
  statements that call TM helices non-exchanging, and the 7 A tail-contact criterion. Measured 11-01 (findings 5.8):
  around a backbone amide the tail shell peaks at 6.9 A, so the test flickers at the bilayer midplane, and every
  resolved TM4 amide exchanges only just past the radius. **Revised 11-01 (user: TM4 and TM5 still show finite
  points with censoring alone):** inside the hydrocarbon core (between the leaflets' ester planes) an amide is
  shielded when a tail sits inside its own tail shell, the first minimum of the core amide N-tail g(r), measured
  like the 7.00 A (8.12 A on glpG); outside the core the 7.00 A contact stays. Censored amides stay off scale. It is
  master's TM-slab AND lipid-facing conjunction, so a lipid-free pocket (S201) stays exposed. Done 11-01 for the four
  ff3.0 sets and slides 17-18: TM5 clears in all four, TM4 in 79HIS; three single-event TM4 points remain (5.8), and
  whether the resolution limit should count correlated frames (RT ln(ESS/g), ~2 kcal/mol lower) is the user's call.
  The Upside 3 tree still holds the pre-core `martini_hdx_membrane_accessibility.py` until the next sync. The
  PF = 1000 uptake cap was replaced in Upside 3 by the per-temperature ESS bound with censored amides (A64, 10-20,
  under the user's 10-19 direction to fix every bug in Upside 3; master unchanged), told to the user 10-20. Do not
  compare uptake with the HXMS data before the remaining items are settled.
* **lambda's helix 2 / helix 4 misorientation** is a known ff2.1 limitation with no single side-chain cause (findings
  11d). It is re-tested after ff3.0, not trained against.
