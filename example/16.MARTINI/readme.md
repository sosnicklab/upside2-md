This example runs hybrid Upside/dry-MARTINI membrane systems with
full-resolution 14-bead DOPC lipids.

Setup:

```bash
source ../../.venv/bin/activate
source ../../source.sh
```

Build the consolidated dry-MARTINI force-field file when the force-field inputs
change (writes `parameters/ff_2.1/martini.h5`):

```bash
bash build_martini_h5_m1.sh
```

Run a full-resolution 1RKL system:

```bash
bash run_sim_1rkl_full.sh
```

The shell wrappers call `run_sim_hybrid.sh`. Common overrides are environment
variables:

```bash
RUN_DIR=outputs/test_1rkl \
RUNTIME_PDB_ID=test_1rkl \
EQ_60_NSTEPS=500 \
PROD_70_NSTEPS=10000 \
bash run_sim_1rkl_full.sh
```

The Python wrapper exposes the same workflow through editable settings near the
top of `run.py`:

```bash
python run.py
```

Outputs are written under `outputs/`. Each run contains stage checkpoints, logs,
and VTF files.

At the stage-7 handoff, the production BB--environment and SC--environment
interactions are active in a dedicated `production_handoff` stage. The complete
protein remains one rigid body through both minimization and burn-in: it can
translate and rotate, but its internal geometry cannot absorb interface clashes.
The workflow relabels the system `production` only after burn-in, which removes
the rigid-body constraint and begins flexible production. No finite positional
spring or gradual release is used.

## HDX analysis

The hybrid simulation does not implement a separate HDX estimator. It projects
the full trajectory into the standard protein-only trajectory contract and then
uses the maintained workflow in `example/00.AnalysisScripts`, which extends the
REMD/MBAR method in `example/04.HDX`.

Analyze the current single-replica 1RKL production trajectory with:

```bash
bash run_hdx_analysis.sh
```

The default analysis work directory is:

```text
outputs/martini_1rkl_hybrid_full/hdx_analysis
```

The wrapper generates the ordinary FASTA, initial coordinates, chain-break file,
and heavy-atom-coverage `1rkl-HDX.up`. It then projects
`checkpoints/1rkl.stage_7.0.up` into an HDX-view file and runs the existing
trajectory, protection-state, uptake, stability, and summary scripts. The
projection keeps the full hybrid potential and `T_up` for MBAR while exposing
only mapped N/CA/C coordinates to the protein analysis engine. The wrapper also
reads `T_up` from replica 0 and uses it as the uptake and stability target, so a
single-temperature run is not silently reweighted to an unsampled condition.
Set `HDX_T_UP` only when deliberately overriding that detected value.

Multiple hybrid replicas use the same wrapper contract:

```bash
N_REP=16 \
HYBRID_SOURCE_DIR=/path/to/hybrid/replicas \
HYBRID_SOURCE_PATTERN='1rkl.run.{replica}.up' \
HDX_WORK_DIR=/path/to/hdx_analysis \
bash run_hdx_analysis.sh
```

`results/*_PS_protein.npy` is the unchanged Upside protein protection state.
`results/*_PS.npy` is the downstream array and is identical unless calibrated
water-accessibility arrays are supplied through `WATER_ACCESSIBILITY_DIR` and
`WATER_ACCESSIBILITY_PATTERN`.

Dry MARTINI contains no explicit water. The workflow therefore does not treat
every amide between global phosphate planes as quantitatively protected; that
rule saturates a transmembrane helix. A membrane accessibility correction must
be calibrated separately against explicit-water MARTINI, atomistic simulation,
or experiment and is applied transparently after the stock protein protection
calculation.

The HDX calculation uses equilibrium open probabilities rather than simulated
time as the labeling clock. The current mismatch between DOPC molecular lateral
diffusion and the empirical protein timescale does not directly multiply the
HDX uptake time. It can still prevent convergence. Do not treat the current
single-trajectory result as quantitative without timestep convergence, block
convergence of protection probabilities, independent replicas or REMD overlap,
and peptide-level experimental validation.

For the supplied 1RKL stage-7 trajectory, molecular DOPC diffusion is about
`0.0152 um^2/s` versus the factor-four-corrected target `11.5 um^2/s`. After the
first 100 frames are discarded, first- and second-half protection probabilities
differ by `0.078` on average across donors and by `0.621` for the worst donor.
This trajectory therefore fails the transport and block-convergence gates. Its
protein-only protection profile is suitable for qualitative inspection, not a
quantitative membrane HDX prediction.
