# Hybrid Upside / dry-MARTINI membrane simulations

This example prepares and runs a membrane protein as a **hybrid** system:

- the **protein** stays on the Upside trained force field — a foldable N/CA/C/O backbone with
  reconstructed sidechain rotamers, so it can change conformation; and
- the **environment** (a lipid bilayer or a detergent such as DDM) is **dry-MARTINI** — a
  coarse-grained, *implicit-solvent* force field (no explicit water beads; ions are explicit).

The protein and environment interact through the Upside side-chain–environment (SC-env) and
backbone–environment (BB-env) interface potentials. All non-bonded interactions are read from
pre-built spline tables in `parameters/ff_2.1/martini.h5`; the C++ engine performs no unit
conversion (dry-MARTINI → Upside units are applied once, at h5-build time).

---

## 1. Inputs you provide

| Input | Purpose |
|---|---|
| all-atom protein PDB | the protein to fold/simulate (any orientation) |
| OPM reference PDB | a membrane-oriented homolog/structure from [OPM](https://opm.phar.umich.edu); its membrane midplane is at z=0. The protein is superposed onto it to set the insertion orientation and depth. |
| CHARMM-GUI Martini membrane | the bilayer/detergent, built with CHARMM-GUI Martini Maker (`gromacs/step5_charmm2gmx.pdb` + `gromacs/system.top`). Either a *dry* (`dry_martini_v2.1`) or a *wet* Martini job works — water is stripped either way. |

The membrane bead names must match a dry-MARTINI topology under `dryMARTINI_itp/` (e.g. `DDM`,
`DOPC`); no atom remapping is done.

## 2. Environment setup

```bash
source ../../.venv/bin/activate
source ../../source.sh
```

Build the consolidated dry-MARTINI force-field file once (and whenever the force-field inputs
change; writes `parameters/ff_2.1/martini.h5`):

```bash
python ../../py/martini_gen_params.py --upside-home ../../
```

## 3. What the preparation does

`martini_prepare_system.py run-hybrid-workflow` performs a standard membrane-protein insertion:

1. **Coarse-grain the protein** with `martinize.py` (MARTINI 2.2) into a BB+SC bead envelope. This
   envelope is a **prep-time geometric aid only** — it is never simulated; the runtime protein is
   the Upside backbone + rotamers.
2. **Orient** the protein by superposing its Cα onto the OPM reference (sequence alignment of the
   shared fold, then Kabsch with iterative core refinement). The membrane normal becomes z and the
   transmembrane region is centered on the membrane midplane. The *same* rigid transform is applied
   to the runtime backbone and the CG envelope, so they stay consistent.
3. **Read the CHARMM-GUI membrane**, stripping water and discarding its ions (both are regenerated
   below), and **tile/crop** it in-plane to the protein footprint plus a PBC-safe belt.
4. **Insert** the protein and **remove overlapping lipids** — any lipid clashing with the CG
   envelope (backbone *and* sidechains) is deleted, so lipids never pack into the space the
   rotamers will occupy.
5. **Size the box**: in-plane = protein span + PBC margin (rounded to whole membrane tiles);
   z = protein span + PBC margin. The box is only as large as PBC safety requires.
6. **Add ions** at `--salt-molar` (default 0.15 M) over the *solvent-accessible* volume
   (box − membrane slab − protein), placed in the solvent slabs outside the membrane. Because the
   dry-MARTINI membrane relaxes to a thickness different from the CHARMM-GUI start, the ion count
   uses the **equilibrated** membrane thickness (`--membrane-thickness-angstrom`) so the production
   system realizes the target concentration (see §5).
7. **Equilibrate (NPT)**: a minimization then a graded release of the lipid-head restraint
   (stages 6.0 → 6.6) under a **semi-isotropic, xy-only, tensionless Monte-Carlo barostat**
   (see §5). This relaxes the membrane area — the CHARMM-GUI box is *not* trusted as-is.
8. **Produce (NVT)** at the equilibrated box (stage 7.0). The whole protein is one rigid body
   through the 7.0 minimization + burn-in, then relabeled `production` for flexible dynamics.

Outputs land under `outputs/<run>/`: `hybrid_prep/` (packed PDB + mapping), `checkpoints/` (per-stage
`.up`), `logs/`, and per-stage `.vtf` trajectories.

## 4. Running

Edit the inputs/settings near the top of `run.py`, then:

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
`--protein-pbc-margin`, `--eq-6*-nsteps`, `--prod-70-nsteps`, `--temperature`, `--prep-seed`,
`--seed`. `--protein-orientation-mode` may be `opm` (recommended for membrane proteins), `input`
(trust the PDB), or `lay-flat` (PCA).

## 5. Physical model and key parameters

**Pressure control (waterless membrane).** dry-MARTINI has no explicit solvent, so nothing sets the
box normal — the barostat must be **xy-only**. Equilibration uses Upside's **Monte-Carlo barostat**
(`box.cpp`): it scales molecule centers of mass (the rigid protein is never distorted), accepts on
the exact NPT weight (no virial needed — correct for a vacuum-containing box), freezes z
(`mc_dmax_z = 0`), and targets **zero lateral tension** (`target_p_xy = 0`), the correct state for a
dry membrane. Compressibility is 3×10⁻⁴ bar⁻¹. Production is NVT at the equilibrated box. This
mirrors the dry-MARTINI CHARMM-GUI protocol (semi-isotropic, `compressibility 3e-4 0.0`,
`ref_p 0 0`, then `Pcoupl = no` for production).

**Ion concentration.** Ions represent bulk 0.15 M ionic strength in the *solvent* region, so the
count is `molarity × N_A × [A_xy·(L_z − t_membrane) − V_protein]`, not the whole box. `t_membrane` is
the **equilibrated** dry-MARTINI thickness: run a short pilot equilibration once, measure the relaxed
lipid z-extent, and pass it as `--membrane-thickness-angstrom` (for DDM this is ≈ 49 Å). Ions are
placed only in the solvent slabs (never the hydrophobic core).

**Reviewing a prepared structure.** Upside stores unwrapped coordinates, so to inspect a stage
`.up`/`.vtf` in VMD, wrap coordinates **per molecule** with the minimum image relative to the
membrane center — otherwise molecules that crossed the periodic boundary appear scattered outside
the box. (This is display-only; the engine uses minimum-image internally.)

---

## 6. HDX analysis

The hybrid run has no separate HDX estimator. It projects the trajectory into the standard
protein-only contract and reuses the maintained workflow in `example/00.AnalysisScripts`, which
extends the REMD/MBAR method in `example/04.HDX`:

```bash
bash run_hdx_analysis.sh
```

The wrapper builds the FASTA, initial coordinates, chain-break file, and heavy-atom-coverage
`<id>-HDX.up`, projects the production `.up` into an HDX-view file (keeping the full hybrid potential
and `T_up` for MBAR while exposing only mapped N/CA/C coordinates), and runs the trajectory,
protection-state, uptake, stability, and summary scripts. `T_up` is read from replica 0 and used as
the uptake/stability target; set `HDX_T_UP` only to override. Multiple replicas use the same
contract:

```bash
N_REP=16 HYBRID_SOURCE_DIR=/path/to/replicas \
HYBRID_SOURCE_PATTERN='<id>.run.{replica}.up' \
HDX_WORK_DIR=/path/to/hdx_analysis bash run_hdx_analysis.sh
```

`results/*_PS_protein.npy` is the Upside protein protection state; `results/*_PS.npy` is identical
unless calibrated water-accessibility arrays are supplied via `WATER_ACCESSIBILITY_DIR` /
`WATER_ACCESSIBILITY_PATTERN`.

**Caveats.** Dry-MARTINI has no explicit water, so a membrane-accessibility correction (calibrated
against explicit-water MARTINI, atomistic simulation, or experiment) must be applied after the stock
protection calculation — otherwise every amide inside the membrane reads as fully protected. HDX
uses equilibrium open probabilities, not simulated time, as the labeling clock; the dry-MARTINI
lipid diffusion / protein-timescale mismatch does not multiply the uptake time but can still block
convergence. Treat a single trajectory as qualitative until timestep and block convergence,
independent replicas/REMD overlap, and peptide-level experimental validation are established.
