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

**Release by checkpoint selection** (user, 10-01/05). At convergence the gate validates its final
epoch-end checkpoint automatically (`validate_ff.sh`: Peng's 32 arms and glpG's 4 REMD variants); a run
that reaches 13 epochs unconverged stops for review. The selection panel (44 all-L CATH domains against
all-atom) chooses, Peng and glpG judge, and the release is the user's.

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

7. [ ] Owed: TM4 at the ends of blocks 2-5 of the 8 chains (remote_jobs.md).

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
4. [ ] The user's decisions (`upside3.md` 7); the audit is Phase 24.
5. [ ] **Training objective terms for HDX, then large proteins, then IDRs** (user 10-28, in that priority; design and
   steps in `upside3.md` 4.1). Terms added to ConDiv's contrast from the per-frame dV/da the engine already returns; the
   energy's functional form stays spline tables.
   - [x] M1, the one measurement the midway2 data allow (10-28, findings 5.7): ubiquitin on the Peng-benchmark REMD.
     The pipeline reproduces Peng's FF2. A full benchmark REMD resolves 14-15 of the 24 measured amides at 273 K, and
     ff2.1 opens the globally exchanging amides within the folded state. ff3.0's ubiquitin melts about 0.06 T_up below
     ff2.1's.
   - [ ] M2-M6 and T1 (`upside3.md` 4.1). M2 waits only on ff3.1's ubiquitin arms finishing; M3-M6 need new runs or
     data; T1 waits on M3 and the temperature mapping (`upside3.md` 7 item 5).

### Phase 24 - Upside 3 audit: find and fix every bug (STARTED 2026-10-19, user)
User 10-19: "the fix only goes into upside 3. You can modify on current branch (martini-dev) as this will be the
'upside 3'. Later in the stage we should make it refuse to run ff2.1 (not yet ...). You can start the debug progress
now, on this branch." Code concise, simple and logical, in master's style (never example 00's). Method in
`upside3.md` 2.2; the defect ledger is findings 3.14.

**Resume here ("continue upside3 progress"; state at 2026-10-22 19:10, Mac with the scratchpad; item 2b to 10-27).**
Working files (fixtures, baselines, build trees, reviewer outputs) are under the gitignored `scratchpad/upside3_audit/`.
Running on midway2 (item 2b; ids, logs and next actions in remote_jobs.md §1 line 3): the Peng arms and the glpG chain.
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
2a. [x] The closing criteria (upside3.md 2.2 item 5) on the final engine, met 10-22 (`scratchpad/upside3_audit/closing/`):
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
   3. [ ] (running) Peng: the 16 proteins x {native, denovo}, the benchmark's own protocol (`ff3_benchmark/bench_run.py`:
      14-replica REMD on Table S2's ladders, dt 0.009, frame 100), ff_2.1, in
      `/beagle3/trsosnic/yinhan/upside3_validation/peng/`; compared with the released-ff2.1 arms on the pre-audit
      binary (`ff3_benchmark/runs/<prot>_<kind>_ff_2.1/`) and Peng's FF2 (Table S2, Fig S4).
   4. [ ] (running) glpG-RKRK-79HIS hybrid with ff2.1 through `example/16.MARTINI/run.py`, the repo's job script
      rewritten for the Upside 3 core (A232; user 10-22: "the job scripts should also be updated to support it"):
      the preparation as a job, then 4 REMD blocks of 48 replicas on T 0.70-0.90, each stopped at `--time-limit` and
      queueing the next. `run_sim_1afo_full.sh` sources the Upside 3 tree too (A233). The campaign's `run_remd.py`
      stays as it is under the running campaigns (Upside 2; it also carries a NaN-rollback guard). At its block ends,
      besides KE, C-N and TM helicity: TM4's protein protection at 298 K (`pp_fail`, dG under PS_protein) with HDX
      through the fixed analysis code (A61 fixed), the follow-up of Phase 16's closed question. Block 1 read 10-26
      (findings 3.15): KE, C-N and five of six TM helices as Upside 2's hybrid, but TM4 is lost (DSSP alpha 0.44 at
      T 0.70 against 0.93-0.96), already in preparation stage 7.0. The loss is ff2.1's own (Known Errors); the HDX
      readout waits on the user's choice of force field.
   Next after it: refusing ff2.1 (the step `upside3.md` names after the audit), then the open items below.
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
      cross-section (upside3.md proposal).
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
   Still with the user: `upside3.md` 7.
6. Git is the user's. The commit should include `tests/` and `example/15.SpatialTransformation/calc_rmsd.py`.

**Revised decisions** (CLAUDE.md Development Rules reworded to match, 10-19):
* martini-dev is Upside 3. Fixes go only here; master is not changed.
* ff2.1 keeps running on this branch for development. Refusing it comes later.
* The harness lives in `tests/` (new). Its generated fixtures and outputs go to the gitignored
  `scratchpad/upside3_audit/`.

**Harness and pass 1 (10-19 to 10-21).**
* `tests/` (method `upside3.md` 2.2): `fixtures.py` (ff_2.1 configs for chig, 1dfn, 1ubq and membrane 1qhj, plus one
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
3. [ ] Read at each round end and at 76.

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
6. [ ] Read at each round end.

### Phase 34 - disordered regions in Upside 3: how accurate, and what could improve them (STARTED 2026-10-28, user)
User 10-28: "for upside3, investigate the possibility of enhancing intrinsic disordered region simulation accuracy";
then "record everything in upside3.md. We will implement this later." The record is `upside3.md` 4.2 (what the model
gives a disordered chain, the literature, B0 and its reading rules, candidates I1-I3); this phase tracks only status.
Decided by the user 10-28: no explicit electrostatic term (the trained terms already carry electrostatics, so it
would count them twice); salt dependence is outside the model.
1. [x] Literature and the SAXS panel (`scratchpad/idp_20261028/lit/`).
2. [ ] B0 on the Mac (4.1 M6): ff2.1 and the SARW as stage A, ff3.1 as stage B (`panel_run.py`, launched 10-28
   12:10, about 5.5 h on 11 cores; `results.txt` ends in DONE). Read with 4.2's rules and record the result there.
3. [ ] I1-I3: after B0, in 4.1's order, on the user's word.

## Known Errors / Blockers

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
  readouts and the 10-22 ff3.0 block figures too (findings 5.6). To regenerate (on the cluster, by the campaign's
  drivers):
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
  and not `broadwl-lc` (its nodes cannot see `/project`).
* **HDX post-processing items await the user's decision** (user, 10-11; findings 5.5): three docstring/README
  statements that call TM helices non-exchanging, and the unmeasured interface boundary of the 7 A tail-contact
  criterion. The PF = 1000 uptake cap was replaced in Upside 3 by the per-temperature ESS bound with censored amides
  (A64, 10-20, under the user's 10-19 direction to fix every bug in Upside 3; master unchanged), told to the user
  10-20. Do not compare uptake with the HXMS data before the remaining items are settled.
* **lambda's helix 2 / helix 4 misorientation** is a known ff2.1 limitation with no single side-chain cause (findings
  11d). It is re-tested after ff3.0, not trained against.
