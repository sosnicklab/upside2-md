# Remote jobs on midway3 — status and handbook

Snapshot: **2026-08-09 23:20 CDT**. Written so a fresh session can pick up cold. Everything needed to
connect, check health correctly, and react to a failure is here.

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

## 1. Current jobs (all RUNNING as of the snapshot)

| JobID | Name | Campaign | Started | Elapsed | Node | Block |
|---|---|---|---|---|---|---|
| 53137729 | `np_1AO6_prod` | **NP** | 08-08 17:44 | 1-05:36 | midway3-0080 | 1/8 |
| 53166591 | `remd_glpG-RKRK-79ALA_S115T` | **glpG** | 08-09 14:44 | 8:36 | midway3-0066 | 2/12 |
| 53167136 | `remd_glpG-RKRK-79HIS_S115T` | **glpG** | 08-09 14:44 | 8:36 | midway3-0087 | 2/12 |
| 53172920 | `remd_glpG-RKRK-79ALA` | **glpG** | 08-09 14:44 | 8:36 | midway3-0097 | 2/12 |
| 53195377 | `remd_glpG-RKRK-79HIS` | **glpG** | 08-09 22:45 | 0:35 | midway3-0154 | 3/12 |

All time limits `1-11:45:00` (35 h 45 m). Health at snapshot: **0 NaN, 0 gate trips, all clean.**
Disk: NP-1AO6 40 G, glpG_DDM_micelle_REMD 48 G (4 x 48 replicas, ~12 G each).

```bash
squeue -u yinhanw -o "%.10i %.28j %.9T %.11M %.11l %R"
sacct -u yinhanw -S 2026-08-08 -X -o "JobID%12,JobName%28,State%22,Elapsed,Timelimit,Start%17,End%17,ExitCode"
```

---

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

**Dir** `~/project/NP-1AO6/` — `prod/` holds `np.run.{0..5}.up` + `np.<jobid>.out`, `block_count`.
**Driver** `run_np_prod.py` · **sbatch** `np_prod.sbatch` (sets `NP_DT=0.001`) · **submit** `submit_np.sh`
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
f=~/project/NP-1AO6/prod/np.53137729.out
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

## 6. The detection gate (what "DESTROYED" in a log means)

Both drivers check after every chunk and **end the chain instead of resubmitting** if a system is
destroyed. It does not prevent blow-ups; it stops silent contamination. It never touches the physics.

- NP `health()`: non-finite positions OR `n_stretched >= NP_CN_COUNT` (5)
- glpG `destroyed()`: non-finite potential anywhere in the chunk OR `n_stretched >= REMD_CN_COUNT` (5)

A gate trip looks like `[np] DESTROYED ...` / `[remd] DESTROYED ...` followed by `no resubmit`, and the
job exits **rc=0 / COMPLETED**. So a COMPLETED job that ran far short of its time limit and did not
resubmit means the gate fired — check the log before assuming success.

**These driver scripts are NOT in git.** They live on the cluster with a gitignored local copy at
`scratchpad/np_slurm/run_np_prod.py` and `scratchpad/np_slurm/run_remd_micelle.py` (the latter deploys
to `run_remd.py`). No version history exists for them. Edit the local copy, then `scp` to deploy.
A running job keeps the version it loaded at start; edits take effect at the **next block**.

---

## 7. If the gate fires — rollback procedure (glpG)

Done once already for 79ALA. Verify the target chunk is clean **before** deleting anything.

```python
# 1. find the last fully clean chunk across ALL 48 replicas (cheap: scan scalar potential)
# 2. then, per replica:
with h5py.File(f, "r+") as h:
    v = np.asarray(h["output_previous_3"]["potential"][:]).reshape(-1)
    if (~np.isfinite(v)).any(): raise SystemExit("target not clean, ABORT")
    for g in ("output","output_previous_4","output_previous_5","output_previous_6","output_previous_7"):
        if g in h: del h[g]
    h.move("output_previous_3", "output")
# 3. resubmit:  bash ~/project/glpG_DDM_micelle_REMD/submit_remd.sh <variant>
```
NP equivalent: `run_np_prod.py` `reseed()` is idempotent, so a config with no `/output` but a valid
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

**glpG**
- Aug 9 (`53088898`, 79ALA): real blow-up. First explosion `output_previous_4` frame 220, **slot 47 =
  the hottest replica (T=0.90)**, t=91.08. Rolled back to `output_previous_3` (~5 chunks / 625 t.u.
  discarded); contaminated trajectories were deleted, so the raw evidence is gone. The job log
  `remd.53088898.out` survives.
- **Mechanism is OPEN.** An earlier claim of "same dt root cause as NP, 8.3× over the limit" was wrong
  twice over: it applied NP's velocity-Verlet analysis to a Brownian-majority system, and the number
  itself was miscomputed (per-pair it is ~0.56× of the harmonic limit — marginal at worst).
  A local single-replica stress at **T=1.30 for t=1080 at fixed dt=0.009 stayed completely stable**, so
  temperature alone does not reproduce it. Leading untested hypothesis: **replica exchange itself**
  (glpG stores no momentum; velocities are re-thermalized each chunk, so how they are handled when a
  configuration moves between temperature slots is a candidate). Next step would be a minimal 2–4
  replica exchange test locally.
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
- `/project` group quota 1.56 T of 3.49 T — comfortable.

**Scratch/snapshot dirs that can be deleted**: `~/project/glpG_snapshot/` (237 M, a copy of
`glpG-RKRK-79HIS.run.47.up` taken for local debugging).
