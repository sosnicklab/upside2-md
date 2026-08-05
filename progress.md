# Progress log

## Current phase (2026-08-05) — NP full-ion rebuild; glpG micelle REMD half-lost to NaN

### NP-1AO6 rebuilt with 0.15 M bulk salt — DONE, submitted
`np_hybrid.build_ions` now takes `salt_molar` + the solvent-accessible volume and spends the neutralizing
excess FROM the salt budget (`salt_molar = 0` reproduces counterions-only, so no dead toggle). Per face:
**K+ 706 / Cl- 488, net +218**, 5174 atoms, solvent 7.82e6 Å³ in the 200 Å box — against the old
counterions-only 218 K+ / 0 Cl-. All six faces rebuilt and passed the build check (finite, no blow-up,
gold_shift 0.00, Rg 26.0–27.2 Å, drift ≤ 0.5 Å). Seeds uploaded and verified on the cluster; production is
**job 53080076**. Counterion builds and trajectories preserved (`counterion_backup/` locally 1.0 GB,
`prod_counterion/` on the cluster 3.2 GB) so the salt-vs-no-salt comparison stays available.

**Open decision:** spending the excess from the budget puts the achieved ionic strength at **0.127 M, not
0.15 M**, because the fixed 218-charge excess is a large fraction of a 706-pair budget. I = 0.15 M exactly
would need K+ 815 / Cl- 597.

### 6-face stability probe — RUNNING
A 2000-step build check cannot answer the stability question: the counterion failure appeared at ~step
97,000 of a 112,966-step chunk, and the glpG seed passed a 400-step smoke test then died at step 12,190.
Seed smoke tests bound loading and setup, NOT stability. So all six faces run a full chunk-equivalent
(112,966 steps, dt=0.005) concurrently at 2 threads each; ~2.3 h. `analyze_probe.py` then reports the
per-frame peptide C-N trace, first breach, and Rg per face.

### glpG micelle REMD — 2 of 4 lost to NaN
Started 00:43. After ~11.5 h: 79HIS NaN from step 12190, 79HIS_S115T NaN from step 6624 (all 48 replicas
each); 79ALA and 79ALA_S115T clean. Mechanism (findings 79): at step 12144 in 79HIS, 47 of 48 replicas were
healthy while `system 16` hit Rg 935,837 Å / potential 9.6e15 — **one replica diverged and replica exchange
carried the non-finite energy into every swap within one step**. Not temperature-driven (systems 16/17 sit
at T=0.77, mid-ladder). Transient positive potential spikes are normal and recover (the survivors peaked at
+36,450 and +41,012 with hbonds 195–202); the rare event is a spike that does not.
Stopped 53036661 and 53036664 with STOP + scancel so they could not resubmit onto corrupted state; left
53036667 and 53036669 running. **`run_remd.py` needs a finiteness guard before the two lost variants are
re-run** — the NP driver already does this correctly via `health()`, the REMD driver never got the
equivalent. The survivors carry the same exposure until it lands.

## Previous phase (2026-08-04) — glpG DDM rebuilt as a MICELLE, REMD re-run

The lamellar-slab glpG-DDM REMD results are void: the CHARMM-GUI DDM slab's tail core is 11–14 Å against
a 28–30 Å protein belt, which unfolded TM4 through its TM4:TM6 GxxxG interface. Diagnosis with numbers is
`findings.md` Update 76; experiment (pD 9 HXMS peptide 140–144) and the implicit-membrane run agree TM4 is
protected at ~10–13 kcal/mol, so the hybrid was the artifact, not the reference.

### Implementation — DONE (`py/martini_prepare_system{,_lib}.py`)
- Morphology derived from acyl-tail count in the lipid ITP (`derive_environment_morphology`): one tail →
  micelle, ≥2 → bilayer. DDM → micelle, DOPC/POPC/POPG → bilayer, detergent+lipid mixtures rejected. No
  name table and no CLI flag that can select the unphysical combination.
- `build_detergent_micelle` wraps a monolayer on the hydrophobic belt (tails in, heads out). Belt comes
  from the OPM reference's own thickness REMARK (14.1 Å half-thickness, matching glpG's measured 28.2 Å
  span). Random sequential adsorption over the shell volume, **innermost first**, spaced by the force
  field's contact distance. Count is derived: 186 DDM for glpG, in the experimental PDC range.
- Barostat forced off for micelle morphology at every stage — a finite aggregate has no lateral tension,
  so a tensionless xy barostat would squeeze the box onto it.
- Ions: `estimate_salt_pairs` takes an excluded volume (bit-identical arithmetic for the bilayer);
  `place_ions` takes one `reject(trial)` predicate — z-band for a slab, nearer-to-tail-than-head for an
  aggregate.
- Gates: **G1** at build time, environment tail span must cover the belt within one contact distance (the
  lamellar DDM slab fails outright); **G2** at the production handoff on equilibrated coordinates, every
  belt backbone site must have an environment bead within twice the contact distance. Build-time
  single-bead bald spots are reported, not fatal.

### Validation — DONE
Short local 79HIS run: compact micelle (env r95 43 Å in a 137 Å box), tail-z span 45 Å over the 28 Å belt,
first-shell distance tightening 7.04 → 5.76 Å, all six TM helices including TM4 at 0.97–1.00 helicity, box
fixed (barostat confirmed off).

Two builder rewrites were needed: convex-hull support-point seeding gave a sparse 32-molecule shell, and
random-order tail-tip seeding smeared tails across the shell (63 Å span, worse coverage). Also, the
CHARMM-GUI step5 file is PRE-minimization with 0.24 Å bead clashes, so it is not a usable packing-density
reference — the force field's contact distance is used instead.

### Cluster — in progress
- Cancelled the four slab REMD jobs (52970530–533) after writing STOP files so the self-resubmit could not
  re-chain. NP-1AO6 (53032366) left running.
- Freed 74 GB: removed the four invalid slab replica trajectory directories. Kept the slab seeds (886 MB,
  provenance) and `hdx_results`.
- New base `/home/yinhanw/project/glpG_DDM_micelle_REMD` (env.sh, run_remd.py, remd.sbatch,
  submit_remd.sh with BASE repointed, HDX driver). Updated prep code synced to the cluster `UPSIDE_HOME`.
- All four variants prepared (`scratchpad/glpg_ddm_micelle/prep_all4.sh`), 20,000-step burn-in each. Every
  one passed the stage-7.0 solvation gate: 0 bare belt sites, 0 beyond acyl-tail reach, local tail core
  45.1–47.5 Å against the 28.2 Å belt, nearest environment bead mean 5.59–5.66 Å.
- Seeds uploaded and md5-verified local vs cluster; four 48-replica T-REMD jobs submitted, all PENDING:
  **53036661** 79HIS, **53036664** 79HIS_S115T, **53036667** 79ALA, **53036669** 79ALA_S115T. Ladder and
  settings unchanged from the slab runs (`linspace(√0.70,√0.90,48)²`, self-resubmitting, 35:45:00). The
  micelle systems are cheaper: 2946 atoms vs the slab's 5011.
- Cancelling the running slab jobs forfeited their queue position, so these start from the back of the
  caslake queue.
- **Seed smoke-tested on the cluster binary before the queue wait** (400 steps, exact REMD flags, login
  node): KE/1.5kT = 1.010, hbonds 194.5 → 201.8, Rg 20.4–20.5 Å, total_potential −7908 → −7626, no NaN.
  The seeds' `/input/pos` and `/input/mom` are fully finite. Worth doing routinely — a seed that fails to
  load would otherwise die instantly after hours of queueing. (A 3-step run reports
  `avg_kinetic_energy/1.5kT -nan`; that is a degenerate-average artifact, not a fault.)
- Remaining: P6 re-run the HDX ΔG analysis once trajectories accumulate; P7 DOPC bilayer regression.

### NP footprinting — stopped itself, needs attention
Job 53032366 **COMPLETED** at 16:01 after 1:42:38 — not cancelled by this phase. Its own health check found
face 1's peptide C–N at 5.29 Å (> 2.5 Å) and it ended the chain rather than propagate a destroyed system
("no resubmit"). The other five faces were healthy (C–N 1.73–1.83 Å). This is the unresolved backbone-tearing
blocker already recorded for that track, now recurring at dt=0.005.

## Notes carried forward
- Git strictly read-only; all edits uncommitted. Master (`upside2-md-master`) parity preserved: the
  bilayer path is unchanged in behaviour.
- midway3 default shell is zsh (1-indexed arrays) — use explicit index mapping in upload loops.
- beagle3 binary needs `module load hdf5/1.14.3` at runtime (libhdf5.so.310).
- A micelle has no fixed normal (asphericity fluctuates 0.19–0.33, long axis drifts 6→15° off box z), so
  depth and tilt must be measured against the aggregate's instantaneous short principal axis, never box z.
  HDX ΔG itself is normal-independent.
