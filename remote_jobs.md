# Remote jobs on midway3 — status and handbook

Snapshot: **2026-08-13 CDT**. Written so a fresh session can pick up cold. Everything needed to
connect, check health correctly, and react to a failure is here. Job state below is live; superseded jobs
are not listed, only summarised in §8 where they carry a lesson.

---

## 0. Connect first (needs a Duo push on the user's phone)

Key-based auth is NOT enabled on midway3; password + Duo is the only method. The ControlMaster socket
expires roughly hourly, so expect to redo this most sessions.

```bash
ssh -S ~/.ssh/cm-mdw3.sock -O check yinhanw@midway3.rcc.uchicago.edu   # alive?
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw3_master.exp    # if not: opens socket, USER MUST APPROVE DUO
ssh -S ~/.ssh/cm-mdw3.sock yinhanw@midway3.rcc.uchicago.edu '<command>'
```

`~/project` is a symlink to `/project/trsosnic/yinhan/` (note: **yinhan**, not yinhanw).
Load the python env on the cluster with `source ~/project/NP-1AO6/env.sh` before any h5py work.

---

## 1. Current jobs

Snapshot **2026-08-17 15:30 CDT**.

| JobID | Name | Campaign | State | Notes |
|---|---|---|---|---|
| 53441123 | `popg_glpG-RKRK-79HIS` | **POPE/POPG** | PENDING | **relaunched 2026-08-17 15:25 from seeds patched to `current_stage = production`** (findings 116). Block 1 of 12, 48 replicas, T 0.70–0.90, dt 0.009, 36 h/block, self-resubmitting |
| 53441124 | `popg_glpG-RKRK-79HIS_S115T` | **POPE/POPG** | PENDING | same |
| 53441125 | `popg_glpG-RKRK-79ALA` | **POPE/POPG** | PENDING | same |
| 53441126 | `popg_glpG-RKRK-79ALA_S115T` | **POPE/POPG** | PENDING | same |
| 53411347 | `np_1AO6_prod` | **NP** | RUNNING 28 h | rolled over from 53372760 (COMPLETED 2026-08-16 10:41); self-resubmitting |

**Cancelled 2026-08-17 15:15: 53410263–66** (the rigid-protein run, findings 116) and their analysis jobs
53439860–63 / 53440383, which completed but produced unusable figures. Their data is archived at
`popepopg_REMD/<variant>_rigidbug/` — **184 G, deletable**, together with the older
`<variant>_pre_cbfix/` (102 G). Seeds backed up as `seeds/<V>.up.bak_production_handoff` before patching.

**MUST VERIFY on the new run, do not skip:** the protein's *internal* CA RMSD between distant frames of the
raw trajectory must be non-zero (local reference: RMSD 2.46 ± 0.563 at T = 0.70). Reading
`current_stage` is not sufficient — that is how the frozen run passed inspection for two days.

**`HDX_LIVE=1` is mandatory while the REMD jobs run.** `run_remd.py` renames `output` to the lowest free
`output_previous_<n>` at every chunk, so the live group can vanish between being listed and being read; the
flag makes `martini_remd_concat.py` skip it, at a cost of ~200 of 5043 frames.

### ⚠ All four POPE/POPG jobs are simulating a RIGID protein — findings 116

`/input/stage_parameters.current_stage` is **`production_handoff`** in the four cluster seeds (and in all
4 × 48 live replica files). `martini_hybrid.cpp:637-641` holds the protein as a rigid group whenever
`preprod_protein_mode == rigid_body` **and** the stage is not exactly `production`. Both are true here.
Measured in the raw trajectory: internal CA RMSD **0.000 Å** across whole chunks; projected Rg 17.60 ± 0.000,
RMSD 0.00 ± 0.003, H-bond 194.6 ± 0.01, identical at T = 0.70 and T = 0.90. The local seed has
`current_stage = production` and moves properly (RMSD 2.46 ± 0.563).

The 2026-08-15 verification note that these were "fine" checked `hybrid_interface_active_stage`
(lines 644-648), which *does* accept `production_handoff`. That is a different gate from the rigid one. One
string, two predicates, opposite accept sets.

Consequences, all downstream: the hole in the bilayer, the −2.5σ potential drift, lipid-shielded fraction
stuck at 0.50 vs the local 0.87, and HDX profiles with 188 of 203 amides off scale. **The cluster HDX figures
are unusable** — no amide can exchange when no H-bond ever breaks.

**Fix (one attribute, not applied):** `set_stage_label(<file>, "production")` on the four seeds and the
4 × 48 replica files, then restart. Costs the 28–30 h × 4 already spent, so it is the user's call.
MBAR itself is fine — ESS 13.4%, ladder overlap 0.25, better than the local run.

**Superseded 2026-08-16:** jobs 53349015/17/18/20 were cancelled at 16 chunks of block 1 and their data
moved to `popepopg_REMD/<variant>_pre_cbfix/` (102 G total, deletable). They carried the 0.568 A CB
placement defect of findings 102. Nothing usable was lost: the protection state has zero variance at that
point in a chain (48 replicas from one common seed need several blocks to decorrelate), so block 1 was not
analysable. The four seeds were patched in place and verified, then resubmitted as 53410263-66.
Verification of the correction before applying it: placing CB by each convention and comparing against the
**actual crystal CB atoms** over 187 residues gives mean deviation 0.576 A (raw, as shipped) versus
0.047 A (corrected) -- the corrected placement is right, independent of any code reading.

**Confirmed loading on the new binary 2026-08-15 00:53**: all four reached `[remd] block 1/12` and
`calibration 2000 steps x 48` with `n_atom 4949`, no `expected 1 arguments but got 2`. The seed migration
held.

**These four were verified rather than resubmitted (2026-08-15).** The BB proxy now reads its backbone O
from `infer_H_O`, so `martini_hybrid_position` takes two arguments and every config written by the older
prep fails to load. Checked, in order: `env.sh` points `UPSIDE_HOME` at the rebuilt
`~/beagle3/yinhan/upside2-md`; the four seeds were upgraded in place with
`py/martini_upgrade_hybrid_args.py` and one was proved to load and run on the new binary; the variant
directories are empty, so `run_remd.py` materialises replicas from those upgraded seeds and there are no
stale copies and no `STOP` files; the seeds carry an `/output`, which the first-block `reseed()` requires;
their stage is `production_handoff`, which `martini_hybrid.cpp:646-647` treats as active, so the SC-env
interface is on. The 36 h `REMD_WALL_SEC` equals the Slurm limit but the chunk loop guards with
`remaining() - MARGIN > last_wall*1.25` (35 min plus a chunk), so the chain still has room to reseed and
resubmit. **No resubmission needed.**

The NP job cannot be patched the same way — it holds its six replicas open for the whole block — so
`np_prod.sbatch` now runs the same upgrade in the gap before `upside` reopens them. **Confirmed working
2026-08-15 00:25**: the block rolled over to 53372760, all six reported `upgraded`, and the block then
loaded every config under the new binary. That line is a one-shot migration — delete it once a block
reports "already upgraded".

**Reading a `.up` that a running job owns.** `run_np_prod.py` reseeds by renaming `output` to the lowest
free `output_previous_<n>`, so on a live file that group can vanish between being listed and being read;
a first footprint attempt died on exactly that (`KeyError: object 'output' doesn't exist`). Rotated blocks
are never renamed again, so an analysis should read only `output_previous_*` and skip the live group —
that costs ~1.4% of frames and removes the race. Set `HDF5_USE_FILE_LOCKING=FALSE` to open at all, and
never let `upside_engine` touch a live file: it opens read-write.

Launcher `popepopg_REMD/{submit_remd.sh,remd.sbatch,run_remd.py,env.sh}`, seeds in `popepopg_REMD/seeds/`,
logs `popepopg_REMD/logs/remd.<jobid>.out`, data `popepopg_REMD/<variant>/`. Seeds are the stage-7.0
checkpoints from `popepopg_glpG/<variant>/checkpoints/`: 4709 atoms (1050 protein / 3393 lipid = 261×13 /
266 ions), box 99.869² × 123.697 Å, 7.48 M pairs, ions regenerated at 0.15 M with
`--membrane-thickness-angstrom 42.75`. This copy of `run_remd.py` uses a **per-chunk RNG seed**; the DDM
copy passes a constant one, which is what made its rollback re-run the identical failing chunk.

Prep history: 53325451–54 failed the stage-7.0 belt gate (34–46 of ~130 belt sites with no environment
bead within 13.92 Å) because the lipid deletion used the LJ *minimum*, `2^(1/6)·σ_max = 6.96 Å`, removing
41 lipids/leaflet where glpG's 822 Å² footprint justifies ~14. Resubmitted as 53332680–83 with
`--protein-lipid-cutoff 4.7` (bead σ, i.e. genuine overlap): gate passed, 261 lipids retained,
min protein–lipid distance 5.053 Å.

### NaN-hunt result (53324867 / 53324868, COMPLETED 9:15–9:18, 340 000 steps, exchange disabled)

**The trigger is real and measured. The LJ core table's floor is directly implicated.**

* Sub-Ångström approaches are **common, not rare**: 4742 and 5032 reported events below 1 Å per run.
* **~1440 events per run reach r < 0.6 Å** (inside the `r = max(r, 0.1*sig)` floored plateau), closest
  **0.0355 Å** and **0.0466 Å**.
* At 0.0355 Å the true dry-MARTINI LJ force is **5.8e29 E_up/Å**; the table's largest delivered force
  anywhere was **3.4e11** — roughly **18 orders of magnitude too weak**, because below 0.1σ the tabulated
  potential is a *constant* and therefore exerts **no force at all**.
* One event sat at **0.0804 Å with the whole box's maximum force only 7.28 E_up/Å** — conclusive proof the
  pair felt nothing.
* 93% of sub-0.6 Å events coincide with a >1e6 E_up/Å force elsewhere in the box (a neighbour being
  launched); 3% have no force above 1e3 anywhere.
* Offending pairs are **environment–environment** (LIPID–LIPID and LIPID–ION; PROTEIN = 0–1049,
  LIPID = 1050–2723, ION = 2724–3155), matching `lipid_kinetic` being what explodes.

Mechanism, with one link still unproven: ION and LIPID are integrated as **overdamped Brownian** (only
420 protein atoms are velocity-Verlet), and an overdamped step is proportional to the force, so a large
force produces a large displacement that can overshoot *through* a partner. Once inside the floored
plateau there is no restoring force to eject it, so it stays and launches its neighbours. **Entry via
Brownian overshoot is inferred, not yet measured**; the floor removing the exit is measured.

This **retracts** the earlier conclusion that the core was ~500 kT out of reach and therefore not the
trigger. It is reached thousands of times per run.

Consequence for the pending POPE/POPG jobs: they use the same table and the same integrator, so they are
expected to hit the same failure. Fix the table core before trusting their output.

## 2. THE TWO CAMPAIGNS ARE DIFFERENT SIMULATIONS

Conflating them caused several errors this session, including a threshold copied from NP that killed a
healthy 6 h glpG block. **Never transfer settings, thresholds, or analysis between them.**

| | **NP** (`np_1AO6_prod`) | **glpG** (`remd_glpG-*`) |
|---|---|---|
| method | regular MD, 6 independent trajectories, single T=0.8647, no exchange | **REMD**, 48 replicas, T ladder 0.70–0.90, configuration exchange |
| purpose | nanoparticle adsorption footprinting (K190 exposure) | **HDX** protection factors / ΔG |
| system | 1AO6 albumin 578 res + 5 nm MPA-AuNP, 4198 atoms, box 200 Å | glpG 210 res in a DDM micelle, 3156 atoms, box 137 Å |
| composition | PROTEIN 2890 + GOLD 887 + MPA 203 + ION 218 (counterions only) | PROTEIN 1050 + LIPID 1674 + ION 432 |
| integrator | **pure velocity-Verlet**, no `/input/brownian` | **MIXED**: 2736 atoms (ION+LIPID+630 protein) overdamped **Brownian**; 420 protein atoms velocity-Verlet |
| timestep | **0.001**, freely settable at runtime | **0.009, HARD-LOCKED** by `/input/brownian/numerical_time_step`; `martini_brownian.cpp:100` throws on mismatch. Friction is tuned against it for lipid D=11.5 µm²/s — **do not change it** |
| detection | non-finite positions OR ≥5 stretched bonds | non-finite potential (whole chunk) OR ≥5 stretched bonds |

---

## 3. NP campaign — `np_1AO6_prod`

**Unfolding is the expected result, not a failure.** 1AO6 albumin spreads on the MPA-AuNP surface —
that is the phenomenon the footprinting measures (K190 exposure). Protein Rg of 38–154 Å against a
native ~27–30 Å is therefore normal for this campaign and must **not** be reported as a blow-up or a
health problem. Judge NP health on non-finite frames, peptide C–N bonds and `avg_kinetic_energy/1.5kT`,
not on Rg or H-bond loss. (Contrast the glpG campaign, where the protein *must* stay folded and Rg
~19 Å is the health signal.)


**Dir** `~/project/NP-1AO6/` — `prod/` holds `np.run.{0..5}.up` + `np.<jobid>.out`, `block_count`.
**Driver** `run_np_prod.py` · **sbatch** `np_prod.sbatch` (sets `NP_DT=0.001`, see §8) · **submit** `submit_np.sh`
**Footprint analysis** `np_footprint.py` (uploaded 2026-08-15) → `np_footprint.npz`, `np_footprint.log`.
**Orientation map** (verify by checksum if ever re-uploading — a zsh 1-indexing bug shifted it once):

```
np.run.0 = 0-0-0   np.run.1 = 90-0-0   np.run.2 = 180-0-0
np.run.3 = 270-0-0 np.run.4 = 0-90-0   np.run.5 = 0-270-0
```

Self-resubmits up to `NP_MAX_BLOCKS=8`; chunks ~104 time units, ~1765 t.u. per 36 h block.
Local build source: `scratchpad/NP-footprinting/` (`build_all.py`, `np_hybrid.py`, six face dirs).

**Expected behaviour:** the protein progressively unfolds on the NP (Rg 26 → 84 Å and still climbing).
This is the intended science, confirmed by the user — do **not** treat rising Rg as a failure. Backbone
stays intact (0 broken bonds). Caveat worth revisiting: on `run.3` the unfolded protein now spans
246 Å in a 200 Å box, so it can interact with its own periodic image.

**Status check:**
```bash
f=$(ls -t ~/project/NP-1AO6/prod/np.*.out | head -1)   # newest block log
grep "^\[np\]" $f | tail -20          # driver trace + per-chunk bond counts
grep -ic nan $f                        # expect 0
```

---

## 4. glpG campaign — `remd_glpG-*` (HDX)

**Dir** `~/project/glpG_DDM_micelle_REMD/` — one subdir per variant, each with 48 `*.run.N.up`,
`remd.<jobid>.out`, `block_count`.
**NOTE the directory name**: `glpG_DDM_micelle_REMD`. The older `glpG_DDM_REMD` (lamellar, 72 G) is gone.
**Driver** `run_remd.py` · **sbatch** `remd.sbatch` · **submit** `submit_remd.sh <variant>`
**Variants:** `glpG-RKRK-79HIS`, `glpG-RKRK-79HIS_S115T`, `glpG-RKRK-79ALA`, `glpG-RKRK-79ALA_S115T`

Config: 48 replicas, T 0.70–0.90, `REMD_DT=0.009`, `--replica-interval 0.09`, `--exchange-criterion 0`,
swap sets A=(0-1,2-3,…) B=(1-2,3-4,…), 300 frames/chunk, `REMD_MAX_BLOCKS=12`.

### HDX analysis for the POPE/POPG campaign (staged 2026-08-15, not yet run)

`popepopg_REMD/hdx_cluster.sbatch`, one job per variant:

```bash
B=/home/yinhanw/project/popepopg_REMD
for V in glpG-RKRK-79HIS glpG-RKRK-79HIS_S115T glpG-RKRK-79ALA glpG-RKRK-79ALA_S115T; do
  sbatch --job-name=hdx_$V --output=$B/logs/hdx.%j.out --partition=caslake --account=pi-trsosnic \
    --nodes=1 --ntasks-per-node=1 --cpus-per-task=16 --time=04:00:00 \
    --export=ALL,PDB_ID=$V $B/hdx_cluster.sbatch
done
```

It runs `py/martini_remd_concat.py` first — **required**, because the chained driver rotates `output` to
`output_previous_<n>` every chunk, so `/output` alone is just the last ~300 frames. The concat joins the
chunks in restart order, renumbers `time`, and drops whole any chunk with a non-finite potential (a
rolled-back chunk is not a sample). Verified against a 3-chunk file: boundary frames in place, all 16
datasets carried, strided output identical to the naive slice. Then the same path as the local run —
`example/00.AnalysisScripts` + `write_hybrid_energy.py` + `plot_ref_style.py`.

**The cluster copy of the analysis pipeline drifts from the repo.** On 2026-08-15 the cluster's
`calc_hdx_ht.py` and `4.calc_D_uptake.py` still lacked the reference subtraction from findings 91, so an
HDX run there would have silently returned uniform MBAR weights. Re-upload the analysis scripts before
trusting any cluster-side result; only the C++ build is kept current by `install.sh`.

**CRITICAL — exchange recirculates a destroyed configuration.** `run_remd.py` reseeds each replica from
`output/pos[-1]`. If that frame is destroyed the next chunk starts destroyed, then an exchange swaps a
healthy configuration in and the bad one moves to another slot. Consequences:

- A single bad configuration appears in many slot files over time (observed walking 3 → 13 → 23 → 33).
- **Per-replica "this file is clean" is meaningless.** Count destroyed configurations *per frame across
  all 48 slots*; the count was 1, then 2, and grew.
- It never self-clears. That is why a detection gate exists.

**Status check:**
```bash
for v in glpG-RKRK-79HIS glpG-RKRK-79HIS_S115T glpG-RKRK-79ALA glpG-RKRK-79ALA_S115T; do
  d=~/project/glpG_DDM_micelle_REMD/$v; f=$(ls -t $d/remd.*.out | head -1)
  echo "$v block=$(cat $d/block_count) nan=$(grep -ic nan $f) DESTROYED=$(grep -c DESTROYED $f)"
done
```

---

## 5. How to check health CORRECTLY

**`isfinite` is not a health check.** At a real glpG failure the environment coordinates were
**±4.65e12 Å** — numerically finite, physically destroyed. And in a forced NP tear the protein reached
**431 broken bonds with the potential still finite at +3e5**, so no energy-based test fires at all.

Use the broken-bond **count** (healthy 0–2; torn 279–431 — a two-order-of-magnitude gap):

```python
import h5py, numpy as np
def n_broken(up, group="output", frame=-1):
    with h5py.File(up, "r") as h:
        pos = np.asarray(h[group]["pos"][frame, 0])
        if not np.isfinite(pos).all(): return -1          # -1 => non-finite
        nm  = np.array([x.decode() for x in h["/input/atom_names"][:]])
        pc  = np.array([x.decode() for x in h["/input/particle_class"][:]])
        rid = h["/input/residue_ids"][:]
    prot = np.where(pc == "PROTEIN")[0]
    C = {int(rid[a]): a for a in prot if nm[a] == "C"}
    N = {int(rid[a]): a for a in prot if nm[a] == "N"}
    res = [r for r in sorted(C) if (r+1) in N]
    cn = np.linalg.norm(pos[[C[r] for r in res]] - pos[[N[r+1] for r in res]], axis=1)
    return int((cn > 2.0).sum())
```

Cheap whole-chunk scan (glpG failures go fully NaN, so the scalar potential catches them):
```python
v = np.asarray(h[group]["potential"][:]).reshape(-1);  bad = int((~np.isfinite(v)).sum())
```
Also useful: protein Rg, and `avg_kinetic_energy/1.5kT` at the end of a log (healthy ≈ 1.0).

---

## 6. The detection gate and rollback (what "DESTROYED" / "ROLLBACK" in a log means)

- NP `health()`: non-finite positions OR `n_stretched >= NP_CN_COUNT` (5) → **ends the chain**
- glpG `destroyed()`: non-finite potential anywhere in the chunk OR `n_stretched >= REMD_CN_COUNT` (5) → **rolls back that replica** (as of 2026-08-12 driver)

**NP**: a gate trip looks like `[np] DESTROYED ...` followed by `no resubmit`; job exits rc=0/COMPLETED
short of its wall limit. A COMPLETED job that did not resubmit means the gate fired — check the log.

**glpG (updated driver)**: a rollback looks like `[remd] ROLLBACK #N filename: reason` followed by
`rolled back M/48 replicas; continuing chain`. The chain does NOT terminate. The NaN output is rotated
to `output_previous_N` as normal history; the rolled-back replica restarts the next chunk from its
pre-chunk positions. A replica that repeatedly blows up gets rolled back repeatedly; it does not get
dropped from the ladder. Watch for high rollback counts on the same file — that indicates a replica with
a persistent physics problem that won't self-correct.

**Rollback mechanism**: before each chunk `run_remd.py` snapshots `/input/pos` (3.6 MB total for 48
replicas). On NaN detection it overwrites the last `output/pos` frame and `output/potential[-1]` with
the pre-chunk values so that `reseed()` on the next iteration picks up the clean state.

**These driver scripts are NOT in git.** They live on the cluster at
`~/project/glpG_DDM_micelle_REMD/run_remd.py` and `~/project/NP-1AO6/run_np_prod.py`.
No version history exists for them. Edit directly on the cluster.
A running job keeps the version it loaded at start; edits take effect at the **next block**.

---

## 7. If the chain terminates — manual rollback procedure

**glpG (old driver, or gate fired before rollback logic):** patch last output frame of each NaN file
with the last finite frame from `output_previous_0` (end of block 1), then resubmit.

```python
import h5py, numpy as np
from pathlib import Path
run_dir = Path("~/project/glpG_DDM_micelle_REMD/<variant>").expanduser()
for fn in sorted(run_dir.glob("*.run.*.up")):
    with h5py.File(str(fn), "r") as h5:
        n_bad = int((~np.isfinite(np.asarray(h5["/output/potential"][:]))).sum())
        has_prev = "output_previous_0" in h5
    if n_bad:
        with h5py.File(str(fn), "r+") as h5:
            prev_pos = np.asarray(h5["/output_previous_0/pos"][-1, 0, :, :])
            last = h5["/output/pos"].shape[0] - 1
            h5["/output/pos"][last, 0, :, :] = prev_pos
            h5["/output/potential"][last, 0] = 0.0
# then: bash ~/project/glpG_DDM_micelle_REMD/submit_remd.sh <variant>
```

**NP**: `run_np_prod.py` `reseed()` is idempotent, so a config with no `/output` but a valid
`/input/mom` (`restart_valid=1`) restarts fine. To reset a chain: `echo 0 > prod/block_count`, `rm -f prod/STOP`.

---

## 8. History and known issues

**NP**
- Aug 2 (`52881141`): all 6 faces destroyed. Root cause **dt=0.005 exceeded the velocity-Verlet limit
  for the MARTINI LJ core**. Proven by A/B from an identical pre-tear frame, only dt differing:
  0.005 → 42 broken bonds, Epot +9130, avgKE/1.5kT 5.34; **0.001 → 0 broken, flat Epot, avgKE 1.006**.
  Fixed by `NP_DT=0.001`. Also rebuilt: fixed 200 Å complex-centred box (was 232–284, orientation-
  dependent) and counterions only (218 ions, was 2442–4324).
- Aug 5 (`53089047`): gate caught a tear on `np.run.3` and correctly stopped the chain.
- Aug 10 (`53234804`): `NP_DT` raised 0.001 → 0.005 based on a 50 t_u test. Runs 1 (90-0-0) and
  4 (0-90-0) destroyed at t≈250 (block 1); runs 0,2,3,5 healthy. Root cause confirmed: **backbone
  spring bonds (k=48 E_up/Å², ω=9.8 rad/t_u) become unstable during protein unfolding on the NP
  surface.** Protein unfolding drives backbone bonds to 3–7× thermal oscillation amplitude by t=250.
  At dt=0.005 the integration has too few force evaluations per half-period (≈1.2 frames/half-period)
  to track the rapid force reversal when bond amplitude is large, allowing a coincident force alignment
  at frame 935 to inject 12 E_up of kinetic energy into a single backbone C atom in 0.27 t_u.
  Note: the MARTINI BB-env LJ is NOT the culprit here (BB closest approach 4.1 Å throughout;
  blow-up 70+ Å from NP surface; residues 119-122). The limiting factor is backbone spring accuracy
  at large amplitude, not small-oscillation stability (ω×dt=0.049 << 2 is met throughout).
  Fixed: reverted `NP_DT=0.001`; deleted /output from destroyed runs; resubmitted as job 53251911.

**glpG**
- Aug 9 (`53088898`, 79ALA): real blow-up, rolled back; contaminated trajectories deleted, job log survives.
- Aug 11 (`53233848`, 79HIS block 2): 45/48 replicas DESTROYED at step 18476/18534 (99.7% through).
  Root cause: **replica 11 physics blow-up** (7 protein atoms went NaN at t=52.45 ps, frame 94/300 of
  the chunk). Replica 11 then spread full NaN (all 1050 protein atoms) to 44 other slots via REMD
  exchange over the next ~50 ps. Contributing factor: **ion escape from seed** — all four variant seeds
  already had ions at 24.9–185.6 Å from the protein centroid before any simulation. Ions beyond the
  12 Å dry-MARTINI cutoff experience zero force → free Brownian diffusion outward. This is a preparation
  design issue in `place_ions` (large virtual box, no PBC, no confinement); ions contribute zero
  energy/force once escaped, but ionic screening is absent. Silent for many blocks because protein and
  micelle remain physically correct. Did NOT cause the blow-up directly; underlying replica 11 failure
  cause still unidentified.
  Fixed: (1) `run_remd.py` updated with NaN rollback (automated per-chunk rollback rather than chain
  termination); (2) 45 NaN-infected `.up` files manually patched with block-1 final positions;
  (3) resubmitted as job 53294511 (block 3).
- **Mechanism RESOLVED Aug 10 (findings 88), and it was never the timestep.**
  `HybridPositionNode::propagate_deriv` discarded the O site's sensitivity entirely and dropped BB's
  fourth (O) weight, while all five backbone sites (N, CA, C, O, BB) sat in the MARTINI pair table. Net
  effect: **more force thrown away than delivered** (12537 vs 5909 E_up/A) and **3590 E_up/A of
  one-sided force per step acting on the environment**, since the env partner's sensitivity was passed
  through while the protein's was not. That is non-conservative, does net work regardless of step size,
  and explains a blow-up localised on one residue's backbone. Fixed by (a) restricting the protein side
  of the pair table to BB and (b) propagating O/BB sensitivities through the placement Jacobian.
  Two earlier dt claims are both withdrawn: the "8.3× over the limit" figure, and a later omega*dt
  analysis of mine that measured O-site contacts which exert no force at all.
  Still open: a **dt-independent +2–3% `avg_kinetic_energy/1.5kT` excess** (2.43/2.16/2.18% across a 4×
  dt range, so not discretisation error). Present in 1rkl/1AFO after the fix too. Cause unidentified.
  Do **not** use `--potential-deriv-agreement` as a gate on these systems: with a 1e-3 A step against a
  ~1e4 E_up float potential it is round-off dominated and returns 0.4196 both before and after the fix.
- Aug 9 (`53166153`, 79HIS): **false positive**, my fault. A `CN_MAX=2.5 Å` threshold borrowed from NP
  fired at 2.52 Å on a perfectly healthy chunk and cost a 6 h block. Healthy glpG reaches 2.659 Å.
  `CN_MAX` has since been removed from both drivers in favour of the broken-bond count.

**Environment gotchas**
- Midway3 **home** is at 28.6 G of a 30 G quota (`code/lammps` 8.4 G, `code/python` 7.9 G, `.cache`
  6.1 G, `.local` 5.9 G). Not simulation data, but jobs can fail oddly if home fills.
- Do not run scripts from `/tmp` on the login node: another user's `/tmp/inspect.py` shadows the stdlib
  and breaks numpy imports. Pipe via stdin or run from `~/project`.
- Cancelling a job is not always a no-op — a 9-second run was enough for `reseed()` to consume `/output`
  on all six NP configs and strand the next job.
- `/project` is 1.6 T of 3.9 T used, **2.3 T free** (2026-08-13). The four POPE/POPG chains are
  estimated at 278 G total, so disk is not a constraint.

**Scratch/snapshot dirs that can be deleted**: `~/project/glpG_snapshot/` (237 M, a copy of
`glpG-RKRK-79HIS.run.47.up` taken for local debugging).

---

## 9. glpG-DDM micelle campaign — closed 2026-08-13

All four variants ended at block 2–3 of 12, every one reporting `COMPLETED` with exit 0 while having
failed: 79HIS_S115T lost 48/48 replicas, 79ALA 38/48, 79ALA_S115T 1/48. A green exit code from this
driver means nothing.

`53294511` (79HIS, relaunched with the rollback driver) was **cancelled**: its ROLLBACK rounds #1/#2/#3
reported byte-identical non-finite frame counts across all 48 replicas, because `run_remd.py` passes a
constant `--seed` every chunk, so each rollback re-runs the identical failing chunk deterministically. It
burned 20.5 h producing only duplicates. Its `STOP` file is still in place (see §1).

HDX ΔG was recovered from the clean pre-blow-up window of every replica and delivered to
`~/Downloads/glpG_DDM_micelle_HDX_dG/`. That required fixing a silent MBAR defect: `calc_hdx_ht.py` and
`4.calc_D_uptake.py` built reduced potentials from raw energies, and a hybrid coupled potential
(~-7.6e3 E_up) overflows `exp(-u)` so the solver never left `f_k = 0` — returning a flat average over the
whole ladder while looking like a reweighted ensemble. See findings 90–91.

The blow-up trigger is **still unidentified**; jobs 53324867/53324868 in §1 are the instrumented hunt for
it. Ruled out by measurement: the backbone over-count, stale pair lists, minimum image, the table build
formula, timestep margin in the sampled range, and thermal access to the LJ core.
