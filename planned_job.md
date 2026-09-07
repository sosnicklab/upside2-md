# Planned job: core-FF retraining → glpG production

Handoff note, written 2026-09-07. **Live job state is in `remote_jobs.md`, not here** — that file is
the single source of truth for job ids, logs and paths. Technical direction is in `plan.md`,
measured knowledge in `findings.md`. This file is only the short version plus what happens next.

---

## What was done

**Retraining the Upside core FF.** ConDiv contrastive-divergence training on midway2, refining
`parameters/ff_2.1` (verified byte-identical start: `init_param/*` md5-match ff_2.1). It trains the
rotamer tables (`pair_interaction`, `coverage_interaction`, `hydrophobe_interaction`) and the
environment energies; hbond and sheet are held fixed. Output goes to
`parameters/ff_3.0_trained/` — a new directory; there is no `parameters/ff_3.0`.

**The ladder was wrong and was fixed.** `ConDiv.py:54` `n_threads` sets *both* the OpenMP thread
count and the REMD replica count (`:397`). It had been set to 4 for node packing, truncating the
ladder from the reference 8 replicas (top T≈0.99) to 4 (top T=0.86) and weakening the CD negative
phase. Now 8 replicas on 4 nodes; thread density per node is unchanged so step time did not move.

**The install pipeline could not have worked, and was fixed.** `extract_ff.py` could not load any
checkpoint (`Target`/`Update` are `__main__` types) — under `set -e` that would have aborted the
install *after* the whole multi-day run. Also fixed: a nested-`sbatch` environment leak that made
every worker die instantly (see `remote_jobs.md`), a Slurm requeue loop that silently restarted
training from a stale checkpoint, and `STEPS_PER_JOB` overrunning the 36 h wall.

**The real cause of the TM4 instability was identified.** A local four-arm, three-replicate
experiment on the real POPE/POPG bilayer seed (143 POPE : 136 POPG, box 99.77² × 180 Å):

| arm | TM4 helix | diverged |
|---|---|---|
| ff_2.1, no coverage nodes | 0.441 [0.298-0.633] | 0/3 |
| trained pair only | 0.588 [0.495-0.668] | 1/3 |
| ff_2.1 + coverage nodes | 0.562 [0.318-0.766] | 0/3 |
| **trained pair + coverage nodes** | **0.782 [0.657-0.863]** | 0/3 |

Neither half works alone. glpG's hybrid build omits `hbond_coverage` and
`hbond_coverage_hydrophobe` — intra-protein terms that dry-MARTINI is structurally forbidden to
supply (`exclude_intra_protein_martini = 1`) — and those are 56% of the trained parameters. Full
detail and caveats in `findings.md` §4.

**Consequence: the 12 h A/B arm test was skipped.** The winner is fixed to ARM_B.

---

## What happens when training finishes — no human action needed

`MAX_STEPS=500`, `STEPS_PER_JOB=130`, `WINNER=ARM_B`. The chain is armed and self-driving; every
`sbatch` is issued from inside a running job, so a login-node outage does not affect it.

```
training reaches 500 → check_continue.sbatch
    extract_ff.py                       -> parameters/ff_3.0_trained/{sidechain,environment}.h5
    back up 4 glpG seeds                -> seeds/<V>.up.bak_pre_ff3_<stamp>
    submit decide_and_launch.sbatch
        WINNER=ARM_B (fixed)
        check_quota.py 120              -> abort if <120 GB headroom, nothing touched
        verify the backups exist        -> abort, nothing touched
        martini_inject_coverage.py x4   -> abort, restore from backup
        patch_seeds.py x4               -> abort, restore from backup
        check_hybrid_up.py --require    -> abort, nothing cancelled
        -- only now -- scancel mdw2_glpG*, delete replicas, block_count=0, submit_remd.sh x4
        remove armtest dirs, log quota
        push_progress.sh                (EXIT trap: runs on every path, including failures)
```

**The safety property:** every failure mode before the `scancel` leaves production running untouched
on the old force field. The single destructive stretch is fenced by `.production_relaunched`.

`push_progress.sh` commits and pushes `TRAINING_STATUS.md` plus the trained FF, so the outcome is
readable off-cluster without a login.

---

## First thing to check in a new session

```bash
ssh -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu   # if dead: expect scratchpad/mdw2_master.exp
D=/project/trsosnic/yinhan/upside2-md-mdw2/training/gly-sym
find $D/run_output -name checkpoint.pkl -path '*/epoch_*/checkpoint.pkl' | wc -l   # steps done / 500
cat $D/.ff_installed $D/.production_relaunched 2>/dev/null                          # has the install fired?
squeue -u yinhanw
```

**Never let an ssh command retry into password auth** — it triggers an RCC IP ban that took ~3 h to
clear and needed a VPN to work around. Check the socket with `-O check` first; if it is dead,
reconnect deliberately with the expect script (the user approves a Duo push).

**After the first production block completes**, measure TM4 against the pre-install baseline in
`popepopg_REMD_mdw2/BASELINE_TM_pre_ff3.txt`. That is the first real evidence about the final FF,
since skipping the arm test removed the pre-production check. Use
`$D/tm_health.py` on the cluster (its dihedral is the validated IUPAC convention — do not "fix" it).
Note it reads only `output_previous_*` groups, so for a fresh run point it at `/output`.

---

## Open items

* **TM4 is likely but not proven fixed.** ARM B measured 0.782 against a >0.8 pass criterion, at
  n=3, single temperature, and with a mid-training (step 269/500) FF. The final FF should do better.
* **Rg trends compact** in ARM B: 19.48 Å mean against a 20.4 Å crystal, one replicate at 18.52.
  This is the direction "over-burying the TM surface" would appear. Watch it in production.
* **The real deliverable is HDX ΔG, not helix fraction.** The recorded problem was helical-core
  CA-RMSD 4.15 Å (bilayer) and H-bond occupancy 0.844. Core RMSD was never measured in the local
  test. Adding it to the production check would speak to the fold, not just the helix.
* **NP is out of the chain**, being rebuilt from scratch instead of patched — its replicas carry 98
  accumulated output groups (246 GB) and its coordinates are old-FF unfolded. `NP-1AO6/build_np_ff3.py`
  and `verify_np_ff3.py` are written and dry-run validated; they need `--rama-library` pointing at
  `parameters/common/rama3.dat` (the default). Not started.
* **The trained environment table has no home in glpG** and never will — the environment node would
  double-count explicit dry-MARTINI lipid. That part of the training is unused by design.
* **Unexplained:** `rama_map_potential` std varied 360-1955 across otherwise identical local runs.
  Do not treat that column as a health metric until understood.
* **The GLY question is unresolved.** With symmetric GLY and ff_2.1, TM4 is still weak (0.441), so
  GLY symmetry was necessary but not sufficient. Whether GLY asymmetry ever contributed
  *independently* is untested — the experiment is a de-symmetrized negative control, a few hours
  locally.
* **The chain scripts exist only on the cluster** (`training/` is gitignored). Cluster `.bak_*`
  copies exist, but nothing is versioned. Same pattern that lost `build_all.py`.

---

## Two traps that cost real time

**There are two rama libraries.** `parameters/common/rama.dat` is the OLD force field and is
asymmetric in GLY *by design* — do not "fix" it (I did, and reverted). The GLY-symmetric library is
`parameters/common/rama3.dat`. Also: the correct GLY mirror is `m[::-1,::-1]` for `.up`
`rama_map_pot` maps but `roll(m[::-1,::-1], 1)` for the library's `dimer_pot`. Using the wrong one
fabricates ~3.3 E_up of asymmetry on perfectly good seeds. Always check against a chiral control
(ALA/SER/HIS must show ~10-11 E_up) and against `dG(aR->aL)`, which reads 0.000 on a correct map.

**Slurm snapshots a batch script at submission.** Editing a `.sbatch` does nothing to jobs already
queued — cancel and resubmit them, or the edit silently has no effect. This caused two separate bugs.
