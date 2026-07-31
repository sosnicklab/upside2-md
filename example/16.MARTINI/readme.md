# Dry-MARTINI / Upside hybrid HDX simulations

This example builds and runs a membrane protein as a **hybrid** system and turns the trajectories
into per-residue hydrogen–deuterium exchange (HDX) protection and uptake:

- the **protein** is on the Upside trained force field — a foldable N/CA/C/O backbone with
  reconstructed side-chain rotamers, so it opens and closes; and
- the **environment** (a lipid bilayer or a detergent such as DDM) is **dry-MARTINI**: coarse-grained
  and *implicit-solvent* (no water beads; ions are explicit).

Protein and environment couple through the Upside side-chain–environment (SC-env) and
backbone–environment (BB-env) interface potentials. Every non-bonded interaction is read from
pre-built spline tables in `parameters/ff_2.1/martini.h5`; the C++ engine does no unit conversion
(dry-MARTINI → Upside units are applied once, at h5-build time).

HDX free energies come from **equilibrium open probabilities**, not simulated labeling time. A single
trajectory rarely samples the transient openings that expose buried amides, so the simulation is run
as **temperature replica exchange (T-REMD)** and reweighted to the target temperature with **MBAR**.
The full workflow is therefore three stages:

```
  Stage 1  prepare  →  one equilibrated production seed  (.up)   per protein variant
  Stage 2  REMD     →  48 replica trajectories at a temperature ladder
  Stage 3  analyse  →  MBAR ΔG_open → protection factors → uptake
```

---

## 0. Inputs you provide

| Input | Purpose |
|---|---|
| all-atom protein PDB | the protein to fold/simulate (any starting orientation) |
| OPM reference PDB | a membrane-oriented homolog from [OPM](https://opm.phar.umich.edu) (midplane at z=0). The protein is superposed onto it to set insertion orientation and depth. |
| CHARMM-GUI Martini membrane | the bilayer/detergent from CHARMM-GUI Martini Maker (`gromacs/step5_charmm2gmx.pdb` + `gromacs/system.top`). Either a *dry* (`dry_martini_v2.1`) or *wet* job works — water is stripped either way. |

Membrane bead names must match a dry-MARTINI topology under `dryMARTINI_itp/` (e.g. `DDM`, `DOPC`);
no atom remapping is done.

## 1. Environment and force field

```bash
source ../../.venv/bin/activate
source ../../source.sh
```

Build the consolidated dry-MARTINI force-field file once (and whenever its inputs change; writes
`parameters/ff_2.1/martini.h5`). `run.py` does this automatically if the file is missing:

```bash
python ../../py/martini_gen_params.py --upside-home ../../
```

---

## 2. Stage 1 — Prepare the equilibrated production seed

`martini_prepare_system.py run-hybrid-workflow` inserts the protein into the membrane and equilibrates
it into a single production-ready `.up`. This is done **once per protein variant**; the result is the
seed that Stage 2 replicates across the temperature ladder.

The preparation is a standard membrane-protein insertion:

1. **Coarse-grain the protein** with `martinize.py` (MARTINI 2.2) into a BB+SC bead envelope. This
   envelope is a **prep-time geometric aid only** — it is never simulated; the runtime protein is the
   Upside backbone + rotamers.
2. **Orient** the protein by superposing its Cα onto the OPM reference (sequence alignment of the
   shared fold, then Kabsch with iterative core refinement). The membrane normal becomes z and the
   transmembrane region is centered on the midplane. The *same* rigid transform is applied to the
   runtime backbone and to the CG envelope, so they stay consistent.
3. **Read the CHARMM-GUI membrane**, stripping water and discarding its ions (both regenerated below),
   and **tile/crop** it in-plane to the protein footprint plus a PBC-safe belt.
4. **Insert** the protein and **remove overlapping lipids** — any lipid clashing with the CG envelope
   (backbone *and* side chains) is deleted, so lipids never pack into the space the rotamers occupy.
5. **Size the box**: in-plane = protein span + PBC margin (rounded to whole membrane tiles);
   z = protein span + PBC margin. Only as large as PBC safety requires.
6. **Add ions** at `--salt-molar` (default 0.15 M) over the *solvent-accessible* volume, placed in the
   solvent slabs outside the membrane (see §5 — ions).
7. **Equilibrate (NPT)**: a minimization then a graded release of the lipid-head restraint (stages
   6.0 → 6.6) under a **semi-isotropic, xy-only, tensionless Monte-Carlo barostat** (see §5 —
   pressure). This relaxes the membrane area; the CHARMM-GUI box is *not* trusted as-is.
8. **Produce (NVT)** at the equilibrated box (stage 7.0). The whole protein is one rigid body through
   the 7.0 minimization + burn-in, then relabeled `production` for flexible dynamics.

Edit the inputs/settings at the top of `run.py`, then:

```bash
python run.py
```

or call the workflow directly:

```bash
python ../../py/martini_prepare_system.py run-hybrid-workflow \
  --pdb-id glpG-RKRK-79HIS \
  --protein-aa-pdb ~/Downloads/glpG-RKRK-79HIS.pdb \
  --lipid-name DDM \
  --charmm-gui-dir ~/Downloads/charmm-gui-<jobid> \
  --protein-orientation-mode opm --opm-reference ~/Downloads/2nr9.pdb \
  --membrane-thickness-angstrom 48.8 \
  --run-dir outputs/glpG-RKRK-79HIS
```

Common knobs (CLI flags or the matching `UPPER_CASE` environment variables): `--salt-molar`,
`--protein-pbc-margin`, `--eq-6*-nsteps`, `--prod-70-nsteps`, `--temperature`, `--prep-seed`, `--seed`.
`--protein-orientation-mode` may be `opm` (recommended for membrane proteins), `input` (trust the PDB),
or `lay-flat` (PCA).

Outputs land under `outputs/<run>/`: `hybrid_prep/` (packed PDB + mapping), `checkpoints/` (per-stage
`.up`), `logs/`, and per-stage `.vtf`. The **stage-7.0 production `.up` is the seed** for Stage 2.

---

## 3. Stage 2 — Temperature-REMD production

REMD samples the opening/closing events HDX depends on and gives MBAR the multi-temperature data it
needs. One OpenMP `upside` process advances **all replicas together**, attempting temperature swaps
between neighbours.

**Temperature ladder.** 48 replicas spanning 0.70–0.90 T_up (≈ 245–316 K; the physiological target
0.8647 T_up = 303 K sits near the top). Spacing is uniform in √T so the exchange acceptance is roughly
constant along the ladder:

```python
import numpy as np
N, T_lo, T_hi = 48, 0.70, 0.90
ladder = np.linspace(np.sqrt(T_lo), np.sqrt(T_hi), N) ** 2     # the --temperature list
SET_A  = ",".join("%d-%d" % (i, i+1) for i in range(0, N-1, 2))  # even pairs 0-1,2-3,...
SET_B  = ",".join("%d-%d" % (i, i+1) for i in range(1, N-1, 2))  # odd  pairs 1-2,3-4,...
```

**Run.** Copy the Stage-1 seed to `<id>.run.{0..47}.up`, then advance a chunk with:

```bash
upside <id>.run.0.up ... <id>.run.47.up \
  --temperature <ladder>              \  # the 48 comma-separated T values
  --replica-interval 0.09 --exchange-criterion 0 \  # attempt T-swaps; 0 = temperature exchange
  --swap-set <SET_A> --swap-set <SET_B> \           # alternate even/odd neighbour pairs
  --time-step 0.009 --integrator v \
  --thermostat-timescale 5.0 --thermostat-interval -1 --seed 1 \
  --duration-steps <n> --frame-interval <dt> \
  --disable-recentering --record-momentum --restart-using-momentum
```

**Continuation.** A wall-guarded driver runs this in chunks and resubmits itself before the wall. Each
chunk: reseed every replica's `/input` pos+momentum from its own last `/output` frame
(`restart_valid`), rename `/output` → `output_previous_K` so per-replica history accumulates in one
file, then run the next chunk. It calibrates ms/step from a short first chunk, sizes chunks to ~2 h,
and stops at a `STOP` file or a max-block cap. Because replicas restart from momenta, temperature
exchange and restarts are seamless.

**Slurm (RCC midway3/beagle3).** One node, one task, `--cpus-per-task=48`, `OMP_NUM_THREADS=48`,
`--partition=caslake --account=<acct>`, `--time=35:45:00` (under the 36 h QOS cap). Run from
`/project` space and load `hdf5/1.14.3` at runtime (the binary needs `libhdf5.so.310`); see the Slurm
setup in the repo `CLAUDE.md`. Submit one job per protein variant.

The result is 48 per-replica `.up` files per variant, each with the accumulated
`output_previous_*` history.

---

## 4. Stage 3 — HDX analysis

The hybrid run has no separate HDX estimator. It projects the trajectory into the standard
protein-only contract and reuses the maintained workflow in `example/00.AnalysisScripts`, which
extends the REMD/MBAR method in `example/04.HDX`:

```bash
bash run_hdx_analysis.sh
```

The wrapper builds the FASTA, initial coordinates, chain-break file, and heavy-atom-coverage
`<id>-HDX.up`; projects each production `.up` into an HDX-view file (keeping the **full hybrid
potential and T_up** so MBAR can reweight, while exposing only mapped N/CA/C coordinates); and runs
the trajectory, protection-state, uptake, stability, and summary scripts. `T_up` is read from replica 0
and used as the uptake/stability target; set `HDX_T_UP` only to override.

For the 48-replica REMD output:

```bash
N_REP=48 HYBRID_SOURCE_DIR=/path/to/replicas \
HYBRID_SOURCE_PATTERN='<id>.run.{replica}.up' \
HDX_WORK_DIR=/path/to/hdx_analysis bash run_hdx_analysis.sh
```

`results/*_PS_protein.npy` is the Upside protein protection state; `results/*_PS.npy` is identical
unless calibrated water-accessibility arrays are supplied via `WATER_ACCESSIBILITY_DIR` /
`WATER_ACCESSIBILITY_PATTERN`.

---

## 5. Physical model and key parameters

**Pressure control (waterless membrane).** dry-MARTINI has no explicit solvent, so nothing sets the
box normal — the barostat must be **xy-only**. Equilibration uses Upside's **Monte-Carlo barostat**
(`box.cpp`): it scales molecule centers of mass (the rigid protein is never distorted), accepts on the
exact NPT weight (no virial — correct for a vacuum-containing box), freezes z (`mc_dmax_z = 0`), and
targets **zero lateral tension** (`target_p_xy = 0`), the correct state for a dry membrane.
Compressibility is 3×10⁻⁴ bar⁻¹. Production is NVT at the equilibrated box. This mirrors the
dry-MARTINI CHARMM-GUI protocol (semi-isotropic, `compressibility 3e-4 0.0`, `ref_p 0 0`, then
`Pcoupl = no` for production).

**Ion concentration.** Ions represent 0.15 M ionic strength in the *solvent* region, so the count is
`molarity × N_A × [A_xy·(L_z − t_membrane) − V_protein]`, not the whole box. `t_membrane` is the
**equilibrated** dry-MARTINI thickness: run a short pilot equilibration once, measure the relaxed lipid
z-extent, and pass it as `--membrane-thickness-angstrom` (for DDM ≈ 49 Å). Ions are placed only in the
solvent slabs, never the hydrophobic core.

**Reviewing a prepared structure.** Upside stores unwrapped coordinates, so to inspect a stage
`.up`/`.vtf` in VMD, wrap coordinates **per molecule** with the minimum image relative to the membrane
center — otherwise molecules that crossed the periodic boundary appear scattered outside the box.
(Display-only; the engine uses minimum image internally.)

## 6. Caveats

- **Membrane accessibility.** dry-MARTINI has no explicit water, so a membrane-accessibility correction
  (calibrated against explicit-water MARTINI, atomistic simulation, or experiment) must be applied
  after the stock protection calculation — otherwise every amide inside the membrane reads as fully
  protected.
- **Timescale.** HDX uses equilibrium open probabilities, not simulated time, as the labeling clock;
  the dry-MARTINI lipid-diffusion / protein-timescale mismatch does not multiply the uptake time but
  can still block convergence.
- **Convergence.** Treat a single trajectory as qualitative until timestep and block convergence,
  independent replica / REMD overlap, and peptide-level experimental validation are established.
