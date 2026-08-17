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
3. **Pick the aggregate morphology** from the amphiphile's own topology. `derive_environment_morphology`
   counts acyl chains as connected components of the apolar (`C1`–`C5`) bond subgraph in the lipid ITP:
   **one tail → micelle**, **two or more → bilayer**. So `DDM` builds a micelle and `DOPC`/`POPC`/`POPG`
   build a bilayer; mixing a detergent with a bilayer-former is rejected as ambiguous. This is not a
   user-selectable flag — a single-tailed detergent with a bulky head cannot form a lamellar phase that
   solvates a TM domain (see §5 — morphology), so the unphysical combination must be unreachable.
4. **Build the environment.**
   * *Bilayer*: **tile/crop** the CHARMM-GUI membrane in-plane to the protein footprint plus a PBC-safe
     belt, then **insert** the protein and **remove overlapping lipids** — any lipid clashing with the CG
     envelope (backbone *and* side chains) is deleted, so lipids never pack into the rotamer space.
   * *Micelle*: **wrap** a monolayer onto the protein's hydrophobic belt, tails inward and heads out. The
     belt comes from the OPM reference's own `1/2 of bilayer thickness` REMARK. Monomers fill the shell
     volume by random sequential adsorption, **innermost first**, spaced by the force field's contact
     distance; the count follows from the geometry (186 DDM for glpG). Nothing is deleted, so there is no
     first-shell void to heal.
5. **Size the box**: bilayer — in-plane = protein span + PBC margin (rounded to whole membrane tiles),
   z = protein span + PBC margin. Micelle — a cube around protein + aggregate plus vacuum padding, so the
   aggregate can reorient without its periodic image changing shape.
6. **Add ions** at `--salt-molar` (default 0.15 M) over the *solvent-accessible* volume — the box minus
   what the aggregate and protein exclude — and never inside the hydrophobic interior (see §5 — ions).
7. **Equilibrate**: a minimization then a graded release of the lipid-head restraint (stages 6.0 → 6.6).
   For a **bilayer** this runs under a **semi-isotropic, xy-only, tensionless Monte-Carlo barostat**
   (see §5 — pressure), relaxing the membrane area; the CHARMM-GUI box is *not* trusted as-is. For a
   **micelle** the barostat is **off at every stage**: a finite aggregate has no lateral tension to
   relax, so box scaling is not a degree of freedom and a tensionless xy barostat would simply squeeze
   the box onto the micelle.
8. **Produce (NVT)** at the equilibrated box (stage 7.0). The whole protein is one rigid body through
   the 7.0 minimization + burn-in, then relabeled `production` for flexible dynamics. Immediately before
   production, `assert_environment_solvation` verifies on the **equilibrated** coordinates that every
   hydrophobic-belt site has an environment bead within twice the contact distance, and reports acyl-tail
   reach and the local tail-core thickness against the belt (see §5 — morphology).

`run.py` drives both stages. Edit the inputs/settings at its top and run:

```bash
python run.py
```

With `continue_sim = False` it builds the seed (Stage 1) and submits the REMD job (Stage 2); re-run
with `continue_sim = True` to extend production. To run **only** the Stage-1 preparation — e.g. to
inspect the packed system before committing to production — call the workflow directly:

```bash
python ../../py/martini_prepare_system.py run-hybrid-workflow \
  --pdb-id glpG-RKRK-79HIS \
  --protein-aa-pdb ~/Downloads/glpG-RKRK-79HIS.pdb \
  --lipid-name DDM \
  --charmm-gui-dir ~/Downloads/charmm-gui-<jobid> \
  --protein-orientation-mode opm --opm-reference ~/Downloads/2nr9.pdb \
  --run-dir outputs/glpG-RKRK-79HIS
```

`--opm-reference` is **required** for a micelle: the detergent belt is placed on the protein's published
hydrophobic band. `--membrane-thickness-angstrom` is a lamellar ion-counting input only and is rejected in
micelle mode, where the aggregate's own volume is measured instead.

Common knobs (CLI flags or the matching `UPPER_CASE` environment variables): `--salt-molar`,
`--protein-pbc-margin`, `--eq-6*-nsteps`, `--prod-70-nsteps`, `--temperature`, `--prep-seed`, `--seed`.
`--protein-orientation-mode` may be `opm` (recommended for membrane proteins), `input` (trust the PDB),
or `lay-flat` (PCA).

Outputs land under `outputs/<run>/`: `hybrid_prep/` (packed PDB + mapping), `checkpoints/` (per-stage
`.up`), `logs/`, and per-stage `.vtf`. The **stage-7.0 production `.up` is the seed** for Stage 2.

---

## 3. Stage 2 — Temperature-REMD production

REMD samples the opening/closing events HDX depends on and gives MBAR the multi-temperature data it
needs. `run.py` copies the Stage-1 seed into `n_rep` replicas and submits one OpenMP `upside` process
that advances **all replicas together**, attempting temperature swaps between neighbours. This follows
the standard Upside REMD driver (see `example/04.HDX`).

**Temperature ladder.** 48 replicas spanning 0.70–0.90 T_up (≈ 245–316 K; the physiological target
0.8647 T_up = 303 K sits near the top). Spacing is uniform in √T so the exchange acceptance is roughly
constant along the ladder — `run.py` builds it with `np.linspace(√T_low, √T_high, n_rep)**2` and the
neighbour swap sets with `ru.swap_table2d`.

**Frame recording (keeps the disk bounded).** Output is written at `frame_interval = 100` upside time
units (~2000 frames per replica over `duration = 200000`) — the Upside standard — and **momentum is
not recorded**; restarts (`continue_sim = True`) reseed **position only** and re-thermalize velocities.
Each replica therefore stays ~120 MB (all four variants ≈ 23 GB), instead of the multi-hundred-GB
growth that an over-dense, momentum-recording, self-accumulating output would produce. Do **not**
add `--record-momentum` or shrink `frame_interval` toward per-step recording.

The submitted command is, in effect:

```bash
OMP_NUM_THREADS=<n_rep> upside <id>.run.0.up ... <id>.run.47.up \
  --duration 200000 --frame-interval 100 --temperature <ladder> \
  --time-step 0.009 --integrator v --thermostat-timescale 5.0 --thermostat-interval -1 \
  --seed 1 --disable-recentering \
  --replica-interval 10 --swap-set <even pairs> --swap-set <odd pairs>
```

**Slurm (RCC midway3/beagle3).** `run.py` submits one node, `--cpus-per-task=n_rep`,
`OMP_NUM_THREADS=n_rep`, `--partition=caslake --account=<acct>`, `--time=35:45:00` (under the 36 h QOS
cap). Run from `/project` space and load `hdf5/1.14.3` at runtime (the binary needs `libhdf5.so.310`);
see the Slurm setup in the repo `CLAUDE.md`. To extend a run, re-run `run.py` with `continue_sim = True`
(it archives each replica's `/output` to `output_previous_N` and restarts from the last frame). Submit
one job per protein variant.

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

`results/*_PS_protein.npy` is the stock Upside protein-only protection state (H-bonded **or** buried by
other protein atoms). `results/*_PS.npy` is that folded with the **membrane-accessibility correction**
(below) via `combine_hdx_protection.py --water-accessibility`. A calibrated external accessibility array
can instead be supplied via `WATER_ACCESSIBILITY_DIR` / `WATER_ACCESSIBILITY_PATTERN`.

**Membrane-accessibility correction (required for membrane proteins).** The protein-only protection has
no lipid term, so a transmembrane amide facing the bilayer — shielded from water by lipid tails but not
buried by protein — drops out of the protected state whenever its backbone H-bond momentarily breaks,
fragmenting the continuous +∞ ΔG that a non-exchanging TM helix must show. `martini_hdx_membrane_accessibility.py`
recovers the missing term from the **full hybrid trajectory** (which still has the lipid): per frame, an
amide backbone N is marked water-**inaccessible** if any lipid hydrophobic-tail bead lies inside the tails'
own first coordination shell. Those amides then read as protected, so lipid-facing TM residues rejoin the
continuous +∞ while genuinely solvent-exposed loops stay exchangeable.

The criterion is the bilayer boundary and it carries **no free parameter**. The radius is the flat first
minimum of the intermolecular tail–tail g(r) measured from the trajectory's own lipids (7.00 Å on the glpG
POPE/POPG patch, reproducible across all 16 replicas); one contact is the boundary because at that radius a
phosphate bead — the conventional edge of the bilayer — has a median of **zero** intermolecular tail
neighbours, against 2 for an ester bead, 6 for the first tail bead and 11 mid-core. The test is local rather
than a slab between the leaflet planes: a slab also protects every amide lining the protein's own polar
interior, which for glpG merges the interfacial loops into the helices (contiguous +∞ runs of 40–50
residues). Master ANDs its slab with a lipid-facing `surface` test for the same reason.

---

## 5. Physical model and key parameters

**Aggregate morphology (why DDM is a micelle).** A lamellar phase's hydrophobic thickness is
`2·V_tail/APL`, and the area per molecule is set by the head. DDM has one C12 tail and a bulky maltose
head (APL ≈ 40 Å²), giving a tail core of only **≈ 11–14 Å** — while glpG's hydrophobic belt is **28.2 Å**
(OPM). Reaching 28 Å would need APL ≈ 23 Å², which a maltose head cannot approach, so **no DDM bilayer can
ever fit a TM domain**; that is exactly why DDM micellizes. Curvature resolves it: on a micelle the heads
spread over a *larger outer* surface while the tails fill an inner core set by the chain length, and the
coverage the protein feels is the belt extent those tails wrap — set by the protein, not by APL.

Building DDM as a slab anyway put half of glpG's belt against maltose headgroups and unfolded TM4 by
letting detergent into its TM4:TM6 GxxxG interface (`findings.md` Update 76). Hence the morphology is
derived from tail count, not chosen. `assert_environment_solvation` hard-fails on any belt site facing
vacuum and warns on the mismatch signature — measured contrast: micelle **45.1 Å** local tail core with 0
sites beyond tail reach, versus **13.9 Å** and 2 sites for the old slab, against a 28.2 Å belt. Note that a
*packed-state* span is not a usable test: CHARMM-GUI templates are pre-minimization and compressed, so
every lipid reads ≈ 20 Å (DOPC 20.4, POPE/POPG 20.6) regardless.

**Pressure control (waterless membrane).** dry-MARTINI has no explicit solvent, so nothing sets the
box normal — for a **bilayer** the barostat must be **xy-only**. Equilibration uses Upside's
**Monte-Carlo barostat** (`box.cpp`): it scales molecule centers of mass (the rigid protein is never
distorted), accepts on the exact NPT weight (no virial — correct for a vacuum-containing box), freezes z
(`mc_dmax_z = 0`), and targets **zero lateral tension** (`target_p_xy = 0`), the correct state for a dry
membrane. Compressibility is 3×10⁻⁴ bar⁻¹. Production is NVT at the equilibrated box. This mirrors the
dry-MARTINI CHARMM-GUI protocol (semi-isotropic, `compressibility 3e-4 0.0`, `ref_p 0 0`, then
`Pcoupl = no` for production). For a **micelle** the barostat is off at every stage — a finite aggregate
has no lateral tension, so there is nothing for box scaling to relax.

**Ion concentration.** Ions represent 0.15 M ionic strength in the *solvent* region, so the count is
`molarity × N_A × [V_box − V_excluded]`, not the whole box. For a bilayer
`V_excluded = A_xy·t_membrane + V_protein`, where `t_membrane` is the **equilibrated** dry-MARTINI
thickness: run a short pilot equilibration once, measure the relaxed lipid z-extent, and pass it as
`--membrane-thickness-angstrom` (for DDM-as-bilayer ≈ 49 Å). For a micelle the aggregate's own measured
volume is used and that flag is rejected. Ions never land in the hydrophobic interior: outside a slab's
headgroup planes, and for an aggregate never closer to a tail bead than to an interfacial head bead.

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
- **A micelle has no fixed normal.** Unlike a periodic slab (pinned by PBC: asphericity 0.21, normal within
  3° of box *z*), a finite aggregate fluctuates in shape and orientation — measured 0.19–0.33 asphericity
  with its long axis drifting 6° → 15° off *z*. Any depth or tilt analysis must therefore use the
  aggregate's **instantaneous short principal axis**, never box *z*, or environment wobble reads as protein
  bending. HDX ΔG is normal-independent, so the wobble does not affect that deliverable.
- **Timescale.** HDX uses equilibrium open probabilities, not simulated time, as the labeling clock;
  the dry-MARTINI lipid-diffusion / protein-timescale mismatch does not multiply the uptake time but
  can still block convergence.
- **Convergence.** Treat a single trajectory as qualitative until timestep and block convergence,
  independent replica / REMD overlap, and peptide-level experimental validation are established.
