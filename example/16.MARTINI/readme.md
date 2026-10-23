# Dry-MARTINI / Upside hybrid HDX simulations

This example builds and runs a membrane protein as a **hybrid** system and turns the trajectories
into per-residue hydrogen–deuterium exchange (HDX) protection and uptake:

- the **protein** is on the Upside trained force field — a foldable N/CA/C/O backbone with
  reconstructed side-chain rotamers, so it opens and closes; and
- the **environment** (a lipid bilayer) is **dry-MARTINI**: coarse-grained
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
| CHARMM-GUI Martini membrane | the bilayer from CHARMM-GUI Martini Maker (`gromacs/step5_charmm2gmx.pdb` + `gromacs/system.top`). Either a *dry* (`dry_martini_v2.1`) or *wet* job works — water is stripped either way. |

Membrane bead names must match a dry-MARTINI topology under `dryMARTINI_itp/` (e.g. `DOPC`, `POPE,POPG`);
bonds and angles are mapped from the ITP by bead name, and a molecule whose bead names differ from the
ITP's is refused.

## 1. Environment and force field

```bash
source ../../.venv/bin/activate
source ../../source.sh
```

Build the consolidated dry-MARTINI force-field file once (and whenever its inputs change; writes
`parameters/ff_2.1/martini.h5`). `run.py` does this automatically if the file is missing:

```bash
python ../../py/martini_build_tables.py --upside-home ../../
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
3. **Require a bilayer-forming lipid.** `require_bilayer_lipids` counts acyl chains as connected
   components of the apolar (`C1`–`C5`) bond subgraph in the lipid ITP and refuses a single-tailed
   amphiphile such as DDM: one tail with a bulky head forms micelles, not a lamellar phase that can solvate
   a TM domain (see §5 — morphology).
4. **Build the environment**: **tile/crop** the CHARMM-GUI membrane in-plane to the protein footprint plus
   a PBC-safe belt, then **insert** the protein and **remove overlapping lipids** — any lipid clashing with
   the CG envelope (backbone *and* side chains) is deleted, so lipids never pack into the rotamer space.
5. **Size the box**: in-plane = the protein's larger span + PBC margin on both axes (it turns freely in
   the membrane plane), rounded up to whole membrane tiles along each axis; z = protein span + PBC margin.
6. **Add ions** at `--salt-molar` (default 0.15 M) over the *solvent-accessible* volume — the box minus
   the membrane slab and the protein outside it — and never inside the hydrophobic interior (see §5 —
   ions). The counterions neutralise the charges the topology carries (the charged chain ends and the
   anionic lipids; side chains carry no explicit charge), and every stage file is checked to be neutral.
7. **Equilibrate** in CHARMM-GUI's order. 6.0 minimises on the soft-core potential (GROMACS soft core,
   lambda 0.01, sc-alpha 4), which pushes the assembly's overlapping beads apart; 6.1 minimises on the
   published potential; 6.2 → 6.6 run NPT dynamics under a harmonic phosphate-height restraint released
   200 → 100 → 50 → 20 → 10 kJ/mol/nm², 500 steps each and 60,000 for 6.6, where the box converges within
   30,000 steps and the protein's lipid face within 60,000 (`findings.md` 8b). The barostat is
   **semi-isotropic, xy-only, tensionless** Monte-Carlo (see §5 — pressure), relaxing the membrane area; the
   CHARMM-GUI box is *not* trusted as-is.
   Before every MD stage the workflow refuses a configuration with a pair inside the repulsive core (pair
   energy above kT ln(1e9 N_pairs)).
8. **Produce (NVT)** at the equilibrated box (stage 7.0). The whole protein is one rigid body through
   the 7.0 minimization + burn-in, then relabeled `production` for flexible dynamics. Production runs at
   least 20,000 steps, and `assert_environment_solvation` takes the belt sites that face lipid (fewer than 60
   protein backbone atoms within 10 Å) over steps 15,000-20,000 and fails the build if more than 10% of them
   have no environment bead within 8 Å; it warns at the same fraction for acyl-tail beads (see §5 —
   morphology; `findings.md` 8b). It measures production because the prep's rigid protein leaves the belt
   ends 7-16% bare whatever production then does, and the flexible protein closes them. Every stage draws its
   own noise (a seed hashed from `--seed` and the stage name), and the momentum and its round pass through
   every handoff, so a production continued with `CONTINUE_STAGE_70_FROM` is the continuous run.

`run.py` drives both stages as Slurm jobs. Edit the inputs/settings at its top and run, on a login node:

```bash
python run.py
```

Without a production seed it submits the Stage-1 preparation and, after it (`afterok`), the first REMD block
(Stage 2); with the seed it submits the next REMD block, until `n_blocks` are done. Each job sources
`job_env_file` and runs `run.py` again with `HYBRID_JOB_STAGE` naming its stage, so the work runs on compute
nodes. To run **only** the Stage-1 preparation, e.g. locally to inspect the packed system before committing to
production, call the workflow directly:

```bash
python ../../py/martini_prepare_system.py run-hybrid-workflow \
  --pdb-id glpG-RKRK-79HIS \
  --protein-aa-pdb ~/Downloads/glpG-RKRK-79HIS.pdb \
  --lipid-name POPE,POPG \
  --charmm-gui-dir ~/Downloads/charmm-gui-<jobid> \
  --protein-orientation-mode opm --opm-reference ~/Downloads/2nr9.pdb \
  --run-dir outputs/glpG-RKRK-79HIS
```

Without `--opm-reference` the stage-7.0 solvation gate cannot run (the belt is undefined) and says so.

Common knobs (CLI flags or the matching `UPPER_CASE` environment variables): `--salt-molar`,
`--protein-pbc-margin`, `--eq-6*-nsteps`, `--prod-70-nsteps`, `--temperature`, `--prep-seed`, `--seed`.
`--protein-orientation-mode` may be `opm` (recommended for membrane proteins), `input` (trust the PDB),
or `lay-flat` (PCA).

Outputs land under `outputs/<run>/`: `hybrid_prep/` (packed PDB + mapping), `checkpoints/` (per-stage
`.up`), `logs/`, and per-stage `.vtf`. `run.py` runs stage 7.0's production for the 20,000 steps the solvation
gate measures; the stage-7.0 `.up` that passed it is the **seed** for Stage 2, which starts from its last frame.

---

## 3. Stage 2 — Temperature-REMD production

REMD samples the opening/closing events HDX depends on and gives MBAR the multi-temperature data it
needs. It runs in **blocks** of one job each: one OpenMP `upside` process advances **all replicas together**,
attempting temperature swaps between neighbours (the standard Upside REMD driver, see `example/04.HDX`), and
stops itself at `--time-limit`, `time_margin` before the wall. The first block copies the seed into `n_rep`
replicas and starts each from the seed's last frame; every later block restarts each slot from its own last
frame with the walker label it held, archives the block's `/output` as `output_previous_N`, and draws its own
seed (`stage_seed(randomseed, "remd.segment<N>")`). A block queues the next one until `n_blocks` blocks are done;
raise `n_blocks` and re-run `run.py` to extend. The chain stops at an engine error; nothing is rolled back.

**Temperature ladder.** 48 replicas spanning 0.70–0.90 T_up (≈ 245–316 K; the physiological target
0.8647 T_up = 303 K sits near the top). Spacing is uniform in √T so the exchange acceptance is roughly
constant along the ladder — `run.py` builds it with `np.linspace(√T_low, √T_high, n_rep)**2` and the
neighbour swap sets with `ru.swap_table2d`.

**Speed and frames.** glpG's hybrid advances ~49 time units per wall hour with 28 replicas on a 28-core broadwl
node and ~28 with 48 (measured on the glpG campaign, 23 ms per system step), so a block covers ~1,000 time
units; `frame_interval = 10` gives ~100 frames per replica per block. **Momentum is not recorded** and blocks
restart from positions, re-thermalizing velocities, so a replica grows by one frame of positions per frame. Do
**not** add `--record-momentum` or shrink `frame_interval` toward per-step recording.

A block's command is, in effect:

```bash
OMP_NUM_THREADS=28 upside <id>.run.0.up ... <id>.run.47.up \
  --duration 1e6 --time-limit <wall - time_margin> --frame-interval 10 --temperature <ladder> \
  --time-step 0.009 --integrator v --thermostat-timescale 5.0 --thermostat-interval -1 \
  --seed <block seed> --disable-recentering \
  --replica-interval 10 --swap-set <even pairs> --swap-set <odd pairs>
```

**Slurm (RCC midway2).** Every job takes the partition and node exclusions from `slurm_args_file`
(`/project/trsosnic/yinhan/slurm/midway2.args`, CLAUDE.md), `--account`, `--time=run_time` (35:45:00, under the
36 h cap), one node with `--cpus-per-task=min(n_rep, cpus_per_node)`, and sources `job_env_file`, the Upside 3
tree's environment (`/beagle3/trsosnic/yinhan/upside3/env_upside3.sh`). Submit from midway2; one chain per
protein variant.

The result is 48 per-replica `.up` files per variant, each with the accumulated `output_previous_*`
history.

---

## 4. Stage 3 — HDX analysis

The hybrid run has no separate HDX estimator. It projects the trajectory into the standard
protein-only contract and reuses the maintained workflow in `example/00.AnalysisScripts`, which
extends the REMD/MBAR method in `example/04.HDX`:

```bash
bash run_hdx_analysis.sh
```

The wrapper builds the FASTA, initial coordinates, chain-break file, and heavy-atom-coverage
`<id>-HDX.up`; projects the production into an HDX-view file (keeping the **full hybrid potential and
T_up** so MBAR can reweight, while exposing only mapped N/CA/C coordinates); and runs the trajectory,
protection-state, uptake, stability, and summary scripts. A single production run is its
`<id>.stage_7.0.up`, `stage_7.1.up`, ... continuations, joined in order into one trajectory (the default
`HYBRID_SOURCE_PATTERN='<id>.stage_7.{segment}.up'`, one replica). A single run's `T_up` is read from its file
and used as the uptake, stability and plotted temperature. A ladder has no one temperature, so a REMD analysis
names its target with `HDX_T_UP` (it must lie on the ladder); `PLOT_T_TARGETS` sets the plotted ones (default `HDX_T_UP`).

For the 48-replica REMD output:

```bash
N_REP=48 HDX_T_UP=0.8647 HYBRID_SOURCE_DIR=/path/to/replicas \
HYBRID_SOURCE_PATTERN='<id>.run.{replica}.up' \
HDX_WORK_DIR=/path/to/hdx_analysis bash run_hdx_analysis.sh
```

`results/*_PS_protein.npy` is the stock Upside protein-only protection state (H-bonded **or** buried by
other protein atoms). `results/*_PS.npy` is that folded with the **membrane-accessibility correction**
(below, `results/*_ACC.npy`) via `combine_hdx_protection.py --water-accessibility`. An external
accessibility array can instead be supplied via `WATER_ACCESSIBILITY_DIR` / `WATER_ACCESSIBILITY_PATTERN`.
Per-residue outputs carry the PDB residue numbering.

**Membrane-accessibility correction (required for membrane proteins).** The protein-only protection has
no lipid term, so a transmembrane amide facing the bilayer — shielded from water by lipid tails but not
buried by protein — drops out of the protected state whenever its backbone H-bond momentarily breaks,
fragmenting the continuous +∞ ΔG that a non-exchanging TM helix must show. `martini_hdx_membrane_accessibility.py`
recovers the missing term from the **full hybrid trajectory** (which still has the lipid): per frame, an
amide backbone N is marked water-**inaccessible** if any lipid hydrophobic-tail bead lies inside the tails'
own first coordination shell. Those amides then read as protected, so lipid-facing TM residues rejoin the
continuous +∞ while genuinely solvent-exposed loops stay exchangeable.

The criterion is the bilayer boundary and it carries **no free parameter**. The radius is the flat first
minimum of the intermolecular tail–tail g(r) measured from the trajectories' own lipids (7.00 Å on the glpG
POPE/POPG patch). It is measured once, from the g(r) pooled over every replica of the analysis, applied to all
of them so the protection state MBAR averages is one function of the coordinates on every rung, and written
to `results/<id>_<sim>_all_ACC_radius.txt`. One contact is the boundary because at that radius a
phosphate bead — the conventional edge of the bilayer — has a median of **zero** intermolecular tail
neighbours, against 2 for an ester bead, 6 for the first tail bead and 11 mid-core. The test is local rather
than a slab between the leaflet planes: a slab also protects every amide lining the protein's own polar
interior, which for glpG merges the interfacial loops into the helices (contiguous +∞ runs of 40–50
residues). Master ANDs its slab with a lipid-facing `surface` test for the same reason.

---

## 5. Physical model and key parameters

**Aggregate morphology (why a single-tail detergent is refused).** *DDM is retired as an environment
(2026-09-09); glpG runs in a POPE/POPG bilayer. This section is the physics behind the tail-count check
that stops any single-tail detergent from being built as a slab.*

A lamellar phase's hydrophobic thickness is
`2·V_tail/APL`, and the area per molecule is set by the head. DDM has one C12 tail and a bulky maltose
head (APL ≈ 40 Å²), giving a tail core of only **≈ 11–14 Å** — while glpG's hydrophobic belt is **28.2 Å**
(OPM). Reaching 28 Å would need APL ≈ 23 Å², which a maltose head cannot approach, so **no DDM bilayer can
ever fit a TM domain**; that is exactly why DDM micellizes. Curvature resolves it: on a micelle the heads
spread over a *larger outer* surface while the tails fill an inner core set by the chain length, and the
coverage the protein feels is the belt extent those tails wrap — set by the protein, not by APL.

Building DDM as a slab anyway put half of glpG's belt against maltose headgroups and unfolded TM4 by
letting detergent into its TM4:TM6 GxxxG interface (`findings.md` Update 76). Hence a single-tailed
amphiphile is refused by its tail count, not by name. `assert_environment_solvation` hard-fails when the
lipid face meets vacuum and warns on the mismatch signature, lipid-facing sites beyond acyl-tail reach. Note that a
*packed-state* span is not a usable test: CHARMM-GUI templates are pre-minimization and compressed, so
every lipid reads ≈ 20 Å (DOPC 20.4, POPE/POPG 20.6) regardless.

**Pressure control (waterless membrane).** dry-MARTINI has no explicit solvent, so nothing sets the
box normal — for a **bilayer** the barostat must be **xy-only**. Equilibration uses Upside's
**Monte-Carlo barostat** (`box.cpp`): it scales molecule centers of mass (the rigid protein is never
distorted), accepts on the exact NPT weight (no virial — correct for a vacuum-containing box), freezes z
(`mc_dmax_z = 0`), and targets **zero lateral tension** (`target_p_xy = 0`), the correct state for a dry
membrane. Compressibility is 3×10⁻⁴ bar⁻¹. Production is NVT at the equilibrated box. This mirrors the
dry-MARTINI CHARMM-GUI protocol (semi-isotropic, `compressibility 3e-4 0.0`, `ref_p 0 0`, then
`Pcoupl = no` for production).

**Ion concentration.** Ions represent 0.15 M ionic strength in the *solvent* region, so the count is
`molarity × N_A × [V_box − V_excluded]`, not the whole box:
`V_excluded = A_xy·t_membrane + V_protein,outside`, where `t_membrane` is the **equilibrated** dry-MARTINI
P-P thickness, between the two leaflets' mean PO4 heights (POPE/POPG 41.8 A, DOPC 43.8 A at the GROMACS
equilibrium with 0.15 M salt, findings 8; pass it as `--membrane-thickness-angstrom`, else the packed lipids'
own P-P is used) and `V_protein,outside` counts only the residues outside that slab, since the slab already
excludes the transmembrane part. Ions never land between the headgroup planes.

**Reviewing a prepared structure.** Upside stores unwrapped coordinates, so to inspect a stage
`.up`/`.vtf` in VMD, wrap coordinates **per molecule** with the minimum image relative to the membrane
center — otherwise molecules that crossed the periodic boundary appear scattered outside the box.
(Display-only; the engine uses minimum image internally.)

## 6. Caveats

- **Membrane accessibility.** dry-MARTINI has no explicit water, and the stock protection is protein-only,
  so without the correction in §4 a lipid-facing TM amide (water-shielded by lipid but not protein-buried)
  reads as *exposed* whenever its H-bond flickers, breaking the +∞ continuity of non-exchanging TM helices.
  The lipid-tail-contact correction fixes this; its `--cutoff`/`--min-contacts` are a first-order physical
  definition and are best re-calibrated against explicit-water MARTINI, atomistics, or experiment.
- **Timescale.** HDX uses equilibrium open probabilities, not simulated time, as the labeling clock;
  the dry-MARTINI lipid-diffusion / protein-timescale mismatch does not multiply the uptake time but
  can still block convergence.
- **Convergence.** Treat a single trajectory as qualitative until timestep and block convergence,
  independent replica / REMD overlap, and peptide-level experimental validation are established.
