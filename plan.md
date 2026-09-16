# ff3.0B: GG-only GLY Ramachandran symmetry, trained and tested on glpG

## Project Goal

ff3.0 symmetrizes the whole central-GLY row of the Ramachandran dimer library. Only the
`GLY|GLY` entries are forced symmetric by achirality; the 38 `GLY|X` entries describe a
chiral local unit and are symmetrized as a modelling choice. Build **ff3.0B**, which
symmetrizes `GLY|GLY` only, retrain the force field against it with the same ConDiv
protocol, and test it on glpG. The question it answers: is ff3.0's glpG improvement
attributable to the achirality-justified GG part, or to the XGX part?

## Architecture and Key Decisions

**The rama library is a training input, so this is a genuine retrain.**
`ConDiv.py:387` sets `rama_library = os.path.join(input_dir, 'rama.dat')`, where
`input_dir` is the run's own `upside_input/`. Every training protein's `.up` is built
against that library, so the fitted `rot` and `env` parameters are conditional on it.
ConDiv trains `rot` and `env` only; `hb` and `sheet` are held at their ff_2.1 init values
(`ConDiv.py:566-573`), which matches the file evidence: `bb_env.dat`, `sheet` and
`hbond.h5` are identical between `parameters/ff_2.1` and `parameters/ff_3.0`, and only
`environment.h5` and `sidechain.h5` differ.

**ff3.0B is defined as the minimal change from ff2.1.** `scratchpad/symmetrize_gly_gg_only.py`
applies the identical mirror and the identical arithmetic average as ff3.0's
`symmetrize_gly_rama.py`, but only to the two `GLY|GLY` entries per group. Verified: those
become symmetric to 0.000000, all 42 coil / 38 sheet `GLY|X` entries stay bit-identical to
`rama.dat`, all non-GLY central rows stay bit-identical, and the ALA chiral control is
untouched. The `GLY|ALL` marginal is deliberately left raw: it averages over mostly chiral
neighbours, and glpG runs in `mixture` mode which never reads it.

**Measured scope of the perturbation, which sets expectations.**

| | GG share | XGX share |
|---|---|---|
| training set dimer entries drawn by glycines (456 proteins, 3311 GLY) | 378 (**5.8%**) | 6145 (94.2%) |
| glpG glycine `dG(aR->aL)` shift from ff2.1 | +1.23 E_up (**6.2%**) | +18.59 E_up (93.8%) |

glpG contains **no GGG** and only two GG pairs (96-97 in TM2, 132-133 in the TM3/TM4 loop).
ff3.0B therefore changes **4 of 210** residue maps, and **zero** inside TM4: residues 136,
143 and 149 are `T-G-V`, `M-G-Y`, `R-G-E`, all bit-identical to ff2.1. TM4 core shift is
+0.000 E_up under ff3.0B against +3.127 under ff3.0.

**Consequence: the only route by which ff3.0B can change TM4 is the retrained `rot`/`env`
tables.** That makes this a sharp falsification test rather than a gradual interpolation:
if TM4 still frays under ff3.0B, ff3.0's benefit is XGX-driven; if TM4 recovers, the
benefit came through the global refit and not through the glycine maps at all.

**Two arms, run concurrently, because the cheap one answers the attribution question on its own.**
The rama-only arm holds `rot`/`env` fixed at `ff_3.0_trained` and varies *only* the rama
maps in an existing glpG seed, which isolates the map contribution with no training. The
training arm delivers the self-consistent ff3.0B the PI asked for.

## Execution Phases

### Phase 0 - scoping (done)
- [x] Census glycine contexts in the 456-protein training set and in glpG
- [x] Confirm `ConDiv.py` reads the rama library, so a retrain is required
- [x] Establish the reference protocol: 500 minibatches, 4 nodes / 12 tasks / 8 cpus,
      broadwl, ~565 s/step healthy (915 s/step degraded), 36 h chain links, 6.4 GB
- [x] midway2 control socket open; `/project` 679 G free, `/beagle3` 1.4 T free

### Phase 1 - build the library (done)
- [x] `scratchpad/symmetrize_gly_gg_only.py`, GG-only mirror average
- [x] Verify GG symmetric, GLY|X bit-identical, non-GLY rows bit-identical, ALA control intact
- [x] Quantify the effect on glpG's 23 glycines against ff2.1 and ff3.0

### Phase 2 - rama-only glpG arms (no training, ~12 h)
- [ ] Install `rama3B.dat` on midway2 under `parameters/common/`
- [ ] Three arms from one glpG seed, `rot`/`env` fixed at `ff_3.0_trained`, only
      `rama_map_pot/rama_pot` differing: `rama.dat` / `rama3B.dat` / `rama3.dat`
- [ ] Submit through the existing `run_arm_test.sbatch` path (28 replicas, 12 h)
- [ ] Score with `tm_health.py` on the corrected windows (TM4 core 135-151, TM1 core 30-48),
      not the raw 131-152 window whose GLY132/133 N-cap caps the ceiling at 0.818

### Phase 3 - train ff3.0B (~78-110 h, 4 chain links)
- [ ] `training/gg-only/` as a copy of `training/gly-sym/` with `upside_input/rama.dat` = rama3B
- [ ] `ConDiv.py initialize init_param upside_input pdb_list run_output`
- [ ] Chain `continue_mdw2.sbatch` to target 500 minibatches, matching the reference exactly
- [ ] `extract_ff.py` -> `parameters/ff_3.0B_trained/{environment.h5,sidechain.h5}`

### Phase 4 - test ff3.0B on glpG (~12 h)
- [ ] glpG arm test with ff3.0B `rot`/`env` + `rama3B.dat`, against ff2.1 and ff3.0 arms
- [ ] Report per-residue helix fraction for TM1 and TM4; decide the attribution

## Known Errors / Blockers

* **`glpG-RKRK-79HIS.up` is not in `popepopg_REMD_mdw2/seeds/`**; only
  `glpG-RKRK-79ALA_S115T.up` is there, while `run_arm_test.sbatch` defaults to the 79HIS
  variant. Pick the seed explicitly rather than relying on the default.
* **`/project` is 83% full (679 G free).** Put the new training run's `run_output` on
  `/beagle3` if it approaches the headroom; one run is ~6.4 GB, so this is a watch item
  rather than a blocker.
* **Five of the user's jobs are already running on midway2**, so a 4-node 36 h chain will
  queue behind them.
* **The predicted outcome of Phase 3/4 is a negative** (ff3.0B behaves like ff2.1),
  because the rama perturbation reaches 0.4% of training map lookups and 0% of TM4.
  That is the point of a control, but it is 4 days of 96 CPUs to confirm a negative.
* **The retrain is below the reproducibility floor of the trainer.** `findings.md:2277` measures two
  ConDiv runs from a *common ancestor* drifting 23 degrees apart in direction and 15-17% rms apart on
  every trained table, and `findings.md:1000` sets the rule: "Any future claim that a retrain
  'improved' the force field must clear the 21% reproducibility floor." ff3.0B's rama perturbation
  reaches 5.8% of the glycine dimer entries, about 0.4% of all residue-map lookups in the training
  set, and 0% of TM4. A single ff3.0B retrain therefore cannot be told apart from refit noise, and
  any TM4 difference it shows would be attributable to the stochastic refit rather than to the GG
  symmetrization. Interpreting it needs replicate retrains branched at step 0, which multiplies the
  cost. The rama-only arm of Phase 2 has no such problem: it is a single-variable change with no
  refit, so its result is interpretable on its own.
* **Build trap that would silently invalidate any arm:** `ConDiv.py` trained with
  `environment_potential_type = 0` (nonlinear) while `upside_config`'s CLI default is `1` (sigmoid).
  32 earlier benchmark jobs "never read the retrained table and measured ff_3.0's sidechain against
  ff_2.1's environment" (`remote_jobs.md:325-331`). Every arm must be built with
  `--environment-type nonlinear` / `environment_potential_type=0`.
* **The hybrid glpG builder has no force-field flag.** `py/martini_prepare_system.py:1921-1925` and
  `py/martini_prepare_system_lib.py:3639` hardcode `parameters/ff_2.1` and `rama.dat`. Arms are
  therefore made by post-build injection into a copy of the seed, which is what the cluster chain
  already does (`py/martini_inject_coverage.py --sidechain`, `scratchpad/rd1_env_test/inject_env.py
  --environment-library`). For the rama-only arm this is an advantage: overwriting
  `input/potential/rama_map_pot/rama_pot` changes nothing else.
