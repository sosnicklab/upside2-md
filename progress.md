# Progress log

High-level execution diary. Job ids, states and log paths live in `remote_jobs.md` only;
technical findings live in `findings.md`; technical direction lives in `plan.md`. The full log up
to 10-07 12:00 is in git history.

---

## Summary to 2026-10-06

* **09-18 to 09-24, trainer.** The ConDiv port then in use trained FF1's form: spline burial, no
  backbone desolvation, no DSE (findings 9t, 9u). It was replaced by Kleinmann's port of Peng's FF2
  trainer, restored to the SI. From ff2.1, 8 of 9 groups held a fixed point (9v). AWH on 40 glycine
  dipeptides finished, and ff14SB agrees with ff99SB-ILDN to 0.045 nats (9r).
* **09-27 to 09-30, round 2 (basin offsets).** ff_3.0 was released 09-30 02:25 and withdrawn at
  20:34 (findings 1.14-1.15). Its glpG chains ran on seeds left with the FF1-form ff_3.0's coverage
  tables (3.11).
* **10-01 to 10-04, round 3 (`ff30_glyhb`).** AWH glycine map, glycine H-bond offsets in
  `HBondEnergy` (12-entry configs bitwise unchanged), side-chain rate / 10, the selection panel and
  release by checkpoint selection. The panels lost helix and the offsets went left-handed (1.17);
  stopped 10-04 21:58.
* **10-02, cleanup.** `training/` was cut to five files and made site-neutral. `findings.md` was
  condensed, and CLAUDE.md's compile and Slurm sections were rewritten (user).
* **10-02 to 10-05, rotamer-BP validation figure for Tobin.** Finished; it is in
  `/beagle3/.../bp_validation/static/`. The deck is not sent.
* **10-04/05, round 4 designed.**
  - Library: BioEmu-fitted glycine library (`rama31.dat`); push probe d -0.014 [-0.041, +0.012];
    non-glycine maps kept NDRD (1.19).
  - Runs: ff30_bio and ff30_gdepth started on caslake; automatic validation deployed.
  - Poly-Gly reference started (Phase 12).
* **10-06, restart at dt 0.009** (user). dt 0.015 destroyed free replicas, and both epoch-0
  panels lost folding and helix (1.21-1.22). The dt 0.009 runs and the matched control
  ff21_ctrl_dt009 started on midway2. TM4 moved to 12 seeds. Those counts were later found invalid
  (3.11).

## 2026-10-07

* **Frozen arm** (user). `ff30_bio_fz`, `ff30_gdepth_fz` and the control `ff21_ctrl_fz` run with
  `hb`, `dhb`, `hbg` and `sheet` at learning rate 0. Initial force fields byte-identical to their
  twins'. ff30_bio_fz started 10:05; at step 0 `hbond.h5` and `sheet` were md5-identical to
  ff2.1's.
* **How Jumper and Peng trained the H-bond** (user question; findings 1.23). Read Jumper's thesis,
  both 2018 papers, Peng's 2022 SI and every trainer.
  - No soluble trainer froze or staged `hb`. It trained jointly at a small rate, and both authors
    stopped early. Without DSE, H-bonds got stronger; Peng's DSE pulls them weaker.
  - Our `hb` rate is twice the SI's.
  - ff_2.0 to ff_2.1 kept `hbond.h5`, `sheet` and `bb_env.dat` byte-identical.
  - Corrected 9e and 1.17: ff2.1's H-bond energies came from FF2's ConDiv.
* **TM4 test defect found and fixed** (user approved; findings 3.11).
  - Defect: `patch_glpg.py` never wrote the two coverage tables. The live seeds carry the FF1-form
    ff_3.0's (rewritten so by the 09-30 deploy), so every TM4 count to date is invalid. At the seed
    frame the rotamer energy is -19.3 as tested against -121.7 with bio_start's own tables.
  - Fixed script deployed: local, `$P/training/` and `tm4_local/scripts/`, through login1 (login2's
    `/project` hangs). The pristine-seed gate passes as before. Notes for other computers are in
    remote_jobs.md and the cluster README.
  - Local inputs re-patched; pre-fix runs moved to `runs_precov_20261007/`.
  - The other Claude session (the watch) coordinated throughout.
* **No remote job needed cancelling.** Trainings, panels and Peng arms build their inputs with
  `upside_config` from the force field's own files. Only patched glpG inputs were affected.
* **TM4 reruns on fixed inputs.** bio_start ran 11:24-~13:10, then b9_01 to ~15:00. c9_00 was
  extracted 11:59 and b9_00 fetched.
  - The driver was replaced twice (a kill tested on a dummy first; running sets untouched). The
    last, `rerun_cov3_20261007.sh`, stops after b9_01 for the 16:00 shutdown. The rest of the queue
    is in remote_jobs.md "Handoff".
* **User correction: TM4 is judged only in the hybrid** (findings 10.19, memory
  `tm4-judged-in-hybrid-only`).
* **Handoff** (user: this computer shuts down at 16:00).
  - plan.md rewritten around the TM4 repair job (530 -> ~190 lines). remote_jobs.md header and
    "Resume here" became one handoff section, and the §1 TM4 block was replaced.
  - Pre-fix TM4 claims corrected in plan.md, remote_jobs.md and findings.md (1.14-1.22, 3.11, the new
    12c); architecture.md and GLY_sym.md moved from round 3 to round 4.
  - findings 1.13-1.22 and 9c-9w condensed from drafts by four agents, each checked and spliced
    (findings.md 4,704 -> 4,325 lines; the four files together 6,817 -> 5,645). Every number in
    each draft traces to its original. One substantive fix: the backbone term's engine derivative
    is `scale` only (9u now agrees with 9v).
* **13:00-13:20, SI-rate arm and the frozen runs** (user). `ff30_bio_si` and `ff30_gdepth_si` were
  submitted: the two designs with every group but rot at the SI's learning rate, otherwise their
  dt 0.009 twins. Built by `checks/si_init_20261007/` scripts, with initial force fields
  byte-identical to the twins'. The user left the frozen runs to my judgement:
  - ff30_bio_fz kept (running; with b9 and bs, the bio design gets three H-bond rates);
  - ff30_gdepth_fz and ff21_ctrl_fz cancelled before they started (queue start 10-08 10:44 and
    23:59).
  - `submit_new.sh` maps `bs_EE` and `ds_EE` (backup `.bak_pre_si_20261007`).
* **bio_start on fixed inputs (12 seeds, 11:24-13:15):** unwound 5, flipped 4, every seed with
  KE/1.5kT 1.001-1.011. The same force field on pre-fix inputs gave 4 and 1 (Fisher p 1.00, 0.32).
  b9_01 started 13:19:30 after a waiter deadlock (findings 10.14) and ends about 15:10. The watch
  records results in findings 1.24.
* **13:35 ff21_ctrl_fz resubmitted** (user): with c9 it is the cleanest test of the hypothesis behind
  the freeze, H-bond drift destabilizing TM4, with no glycine change. That hypothesis rested on
  pre-fix TM4 counts and an unresolved b9_00 difference (p 0.41), and is untested on fixed inputs.
* **14:00-14:10 poly-Gly production** (user approved). The collapse was read with gmx gyrate,
  polystat and mindist -pi: compact within 4 ns, reopening, image distance >= 5.2 nm. Built four
  7.5 nm dodecahedron starts (~28,950 atoms) and submitted the self-chaining `prod.sbatch` as
  49204428. The collapse ran at 14 ns/day, not the planned 55; 1 us may take about a month.
* **b9_01 on fixed inputs (12 seeds, 13:19-15:14):** unwound 6, flipped 5, against bio_start's 5 and
  4 (Fisher p 1.00). Mean last-block TM4 0.83 against 0.90, Mann-Whitney one-sided p 0.23.
  Unresolved; slightly toward unwound. The local queue ended for the shutdown; "continue jobs" on
  the next computer starts at remote_jobs.md "Resume here" (pointer in CLAUDE.md).
