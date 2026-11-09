# A modified Boltz-2 to replace Upside's sampling for HDX

**Status (2026-11-06).** Phase 0 is done (10-21). Phase 1's HDX half is read out: Boltz-2 ranks EHEE's amides at 10^5
samples (censored Spearman 0.74; Upside FF2 0.61) but not HEEH_rd4_0097's at 10^4 (0.04; Upside 0.32) or Ub psWT's at
10^5 (-0.21, and -0.06 with a query-only MSA; Upside 0.43). By experiment's urea dependence EHEE's measured amides
exchange through subglobal openings, HEEH's helices by global unfolding and Ub psWT's amides below the bound locally,
so Gate 1 ends in its "wrong states" outcome; its mdCATH half is not run. Scored as well by Peng's second metric, the
m-value (11-03, the PI), Boltz-2 opens every measured amide locally: its m does not follow HEEH's (Spearman -0.18;
Upside FF2 0.475) and is unresolved on EHEE and Ub psWT at 10^5, so the verdict stands. In the side test BioEmu2's m
follows HEEH's as well as FF2's (0.465), but it opens EHEE only globally (0.18; FF2 0.40). Phase 2 (user 10-27: no
retraining of Boltz, small trained add-ons allowed) is a trained add-on that tilts frozen Boltz-2's ensemble; no
fine-tune is planned. Phase 2's E1 is held (10-28): above kappa 1 the Feynman-Kac `log_z` falls under its provable
floor, and the diagnosis job (10-29) did not reproduce the tuning's runs at the same seed, so its test was not read;
a reproducibility job (10-29) showed that no steered run reproduces its seed across processes, even on one GPU, so a
repeat test must be statistical. Read statistically (10-29), the diagnosis pointed to a defect in the weights at low
power; a CPU toy with a known normalizing constant (10-29) then found `sampler.py`'s weights and `log_z` right on
Boltz's sampler, which points I61's shortfall to the estimator's variance with Boltz-2's network. More particles
(1,024 and 4,096, read 11-03) do not reduce it, and E1 is closed (user 11-03). Phase 3 (11-03) keeps Boltz-2 frozen
and replaces the steered sampler with replica exchange whose move is Boltz-2's own renoise-denoise, so no path weight
is carried. Its toy test T0 passes (11-03): MBAR recovers the toy's rare-event truth (P(d >= 7.5) 4.15e-4 against
3.70e-4). On Boltz-2 (T1-T1d, 11-03 to 11-06) every move that remembers its state leaves the chains more
protected than Boltz-2's own samples, by 2.4-3.9 standard errors in mean N, and the moves that do not are fresh
draws. At the user's word T2b ran the ladder on T1d's move anyway, its calibration measuring whether that bias
matters at HDX's precision (Decision 7). Read 11-07, it fails: I61 reads 0.47 of plain (the limit 0.5) and the ladder
opens none of V70, I13, W45 and L67 in 20 states, so by its rule Phase 3 closes at plain resolution, which is the
user's call. The user chose to go on (11-07): T2c tilts each of the four amides' own opening excess, with I61 as
calibration, and reads Boltz-2's open fraction there at a decade's precision (Decision 8). Read 11-08, T2c fails its
calibration: I61's ladders are not stationary (+3.7 standard errors against the limit 3), so no amide is read; beside,
the four targets' estimates sit 3.4-7.3 decades under experiment, on tails not calibrated. The user chose to go on
(11-08): T2d bounds the four amides with 5 x 10^5 more of Boltz-2's exact plain samples, at a decade, with no move
to bias it (Decision 9). Read 11-09: V70 opens once in 5 x 10^5, too stable by more than a decade; I13, W45 and L67
open 6, 5 and 2 times, 6-15 times under experiment and unresolved at a decade; the set is read alone, because D52
opens 4% less often than in the 10^5 (p 1.4e-4). Reported. The user asked for the next job (11-09): E2a tests
Phase 2's premise, Boltz-2's own order, on Ferrari et al.'s 11 new HDX-NMR domains (Decision 10). Jobs are recorded in
`remote_jobs.md` §1. Phase 0-1's
run data on `/beagle3` and BioEmu2's samples were deleted 11-03 at the user's word (logs and scripts kept; the pooled
sets the readouts use are in the Mac Studio's `scratchpad/boltz_hdx/flags/samples/`; `remote_jobs.md` disk section).
Phase 2's `.venv` and `cache/mols` were packed 11-03 for `/beagle3`'s file quota: unpack both before running E1 again
(the two `tar -xf` commands in that section). Every number below comes from reading code, docs and papers on
2026-10-21, or from the dated measurements; numbers marked *estimate* are unmeasured.

## 1. Goal

Build a modified Boltz-2 that samples protein conformations for HDX in place of Upside (user, 2026-10-21). Judge it by
the per-residue HDX dG it gives through the repo's estimator, against experiment and against Upside, and judge its
dynamics against all-atom MD on mdCATH. Find out whether a route with no training is good enough before spending GPU
time on training.

**Built on Boltz, not BioEmu** (user, 2026-10-21): Jumper does not accept BioEmu because its data are selected. No
BioEmu weights, samples, training data or benchmarks are used (findings 10.34). BioEmu2's code may be used (user,
2026-10-21): it is MIT-licensed, and its layer, denoising and steering code is itself adapted from Boltz (its
`NOTICES.md`), so it ports onto Boltz's modules directly. The one exception is section 7's side test of BioEmu2's
weights for HDX (user, 2026-10-24), kept separate from this work.

## 2. Background

### 2.1 What the HDX estimator needs from a sampler

From findings 5.1 and 5.2:
* The estimator is equilibrium. Each frame gets a binary protected flag per amide (backbone H-bond score, an Asp/Glu
  side-chain contact proxy, burial); MBAR gives the protected fraction p, and dG = RT ln(p/(1-p)) under EX2. Frame
  order is not used. Dynamics matter for decorrelation, and for whether EX2 holds.
* The flags need only N, CA and C per residue. The hybrid already feeds the estimator that way (a projection onto
  `3*n_res` coordinates), so Boltz frames go through the same code.
* With binary flags the largest resolvable dG is about RT ln(N_eff): 4.1 kcal/mol at 298 K from 1,000 independent
  samples, 6.8 from 10^5. Protected amides in stable proteins sit at 5 to 12 kcal/mol, open fractions of 10^-4 to
  10^-9.
* Upside reaches its present limit (5.5 kcal/mol on glpG) with REMD and MBAR, which need an energy and a temperature.
  Boltz-2 has neither: `temperature_bins` is defined in `src/boltz/data/const.py:468`, but no model code reads it. Boltz
  sample frequencies are therefore the only population estimate, and nothing in its training calibrates them as
  Boltzmann weights.
* Boltz sees no membrane, so glpG and the TM4 question stay with Upside. Tests use soluble proteins.

### 2.2 Boltz-2 at `b1ebfc4`

Source `github.com/jwohlwend/boltz` (version 2.2.1, Python 3.10 to 3.12). Paths are relative to `src/boltz/`.
* Prediction runs in `bf16-mixed` (`main.py:1262`). cuEquivariance kernels are used unless `--no_kernels`.
* `--method` picks a tag from `data/const.py:440` (`md` is 0). The tag is embedded and added to the trunk's single
  representation (`model/modules/trunkv2.py:199`), so it reaches the diffusion module only through the trunk.
* MD training data: poses from MISATO, ATLAS and mdCATH, used in only 8k of 97k training steps (stages 2 and 3). The
  preprint evaluated MD conditioning on held-out mdCATH and ATLAS clusters and reports MD-conditioned ensembles less
  diverse than MD.
* The trunk (MSA module and PairFormer) is computed once per input and reused for every diffusion sample.
* Diffusion module, `model/modules/diffusionv2.py`: `AtomDiffusion.sample` (line 295), `preconditioned_network_forward`
  (251), `sample_schedule` (276). The reverse step's noise is off below `gamma_min = 1.0` (lines 191, 336). Training
  noise is `sigma = 16 A * exp(-1.2 + 1.5 N(0,1))` (lines 538-541). The default `step_scale` is 1.5 for Boltz-2
  (`main.py:1231`), which sharpens samples toward one structure; 1.0 is the unscaled sampler.
* Boltz-2 training code is not released ("coming soon" in its README and `docs/training.md`). Only Boltz-1 can be
  retrained with official code.

### 2.3 Why time-lagged data is needed for dynamics

An equilibrium ensemble fixes the energy landscape, not the timescale. Upside gets dynamics by assuming Langevin motion
on an energy trained with contrastive divergence; its times are in natural units (`README.md:102`), and the hybrid's
clock is calibrated from a measured lipid diffusion (`up.md:188`). A diffusion model has no energy to integrate, and a
lag input means something only if training used pairs of MD frames separated by that lag (an implicit transfer
operator, Schreiner et al. 2023). This was the basis of the 10-21 plan's lag-conditioned fine-tune, now dropped
(Decision 4).

### 2.4 Why the Arts et al. 2023 score-as-force method is not used on Boltz

Arts et al. used one Calpha bead per residue, proteins up to 56 residues and equilibrium MD samples. Their score is the
gradient of a scalar energy, which makes the force conservative; a direct-vector score was about 10x worse on
chignolin (TIC JS). Friction and temperature were copied from the MD, and they still report no mapping from
coarse-grained to real time. Boltz's denoiser outputs per-atom vectors (not conservative) and is not equivariant, and
of its training noise draws 6.5% fall below 0.5 A and 0.5% below 0.1 A, so the score at bond-vibration scale is poorly
trained. Plainer et al. 2025 found small-noise scores inconsistent (Fokker-Planck violations) even in purpose-built
models, and fixed it with a training-time regularizer.

Onsager-Machlup action minimisation (Raja et al., ICML 2025, OM-TPS; read 11-06, user: "check if Onsager-Machlup
action optimization can be used in Boltz job") is the same score-as-force idea applied to paths. The learned score at
a small nonzero diffusion time is taken as the drift of an overdamped SDE, and a discretised path between two given
endpoints is moved by gradient descent to minimise the OM action: a smoothness term, the drift's norm and its
divergence, the last by Hutchinson's estimator. It returns one most-likely path per optimisation, not the path
ensemble, which the authors state. Their free energies and rates came from umbrella sampling, committors or MD with
a physical force field around the paths: alanine dipeptide's barrier of about 6 kcal/mol; a Muller-Brown rate of
1.3e-5 against the true 5.4e-5. Their models were trained on equilibrium samples of the system itself (a 2-D
potential, alanine dipeptide, Arts et al.'s per-protein C-alpha models of the fast folders, a flow model on
tetrapeptide MD), which is what makes the score a force. It is not used for Boltz-2's HDX, for three reasons:
* HDX needs the populations of opened states, which are endpoint free energies. An action is not one without a
  potential, and Boltz-2's non-conservative denoiser defines none.
* The opened endpoints (V70, I13, W45, L67) are the states Boltz-2 never produces, so there is nothing to connect the
  native to.
* The drift would be Boltz-2's score at small noise, the regime where its moves leave the chains more protected (T1b,
  T1d) and where the reasons above apply.
If T2b's support fails, the same ladder can tilt on a coordinate aimed at one amide (its burial, as E1's potential was)
instead of N. MBAR's reweighting stays exact for any coordinate, and it needs no path or endpoint.

### 2.5 Cost reference points

Gathered for the dropped fine-tune (Decision 4). Boltz-1: 68k steps, batch 128, crops of 384 to 512 tokens, no GPU
count given; AlphaFold3 about 150k steps, batch 256, about 4x Boltz-1's compute (Boltz-1 report). Boltz-2: 97k steps
in 4 stages, affinity models on 128 A100s (preprint). AlphaFlow's MD fine-tune of full AF2: 11 h on one 8x A100
machine, about 90 A100-hours; its PDB stage (which made AF2 generative, and which Boltz does not need) 267 h
(AlphaFlow, Table 2).

## 3. Decisions

1. **Boltz-2 at `b1ebfc4` is the base model.** No BioEmu weights, samples, data or benchmarks; BioEmu2's code may be
   ported (section 1).
2. **Phase 1: independent `md`-tagged Boltz-2 ensembles through the HDX estimator, with no training.** This is the
   reference every modification has to beat.
3. **Phase 2: a small trained add-on tilts frozen Boltz-2's ensemble** (user 2026-10-27), sampled through Boltz-2's
   own steering and trained on experimental HDX-MS. The errors are site-specific inside the folded ensemble, so a
   global term cannot fix them, and experiment is the only teacher whose ceiling is above Upside (Phase 2).
4. **No fine-tune of Boltz-2** (user 2026-10-27: no compute for it). The renoise-denoise chain and the lag-conditioned
   fine-tune of the 10-21 plan are dropped; both targeted time order, not populations (assessment 2026-10-24).
5. **No Arts-style score-as-force dynamics on Boltz** (section 2.4).
6. **Phase 3: replica exchange on frozen Boltz-2, its own renoise-denoise as the Metropolis move** (2026-11-03,
   after E1 closed). It replaces Phase 2's Feynman-Kac sampler: a move that leaves Boltz-2's distribution unchanged
   needs only the tilt in its acceptance, so no path weight can collapse. Decision 4 dropped renoise-denoise as a
   model of time order; here it is an equilibrium move, and nothing is read from the chain's order.
7. **The ladder runs on T1d's move, its bias measured by T2's calibration** (2026-11-06). Every move that remembers
   its state failed T1's invariance rules toward protection (T1b's moves at 3.95 A and below, T1d's from the 1
   A-noised state), and the moves that keep the distribution (T1c's, 7.67 and 15.95 A) are fresh draws (T2a). Closing
   Phase 3 at Boltz-2's plain resolution was proposed at 11:40; the user chose to go on ("move to next Boltz job and
   submit it for me"). T1d's failure stays on the record: its chains sit 0.10 NH more protected, a 0.17% shift that
   T1's rules resolve. T2's calibration, a factor of 2 on five amides Boltz-2 opens in 12-80 of 10^5 plain samples,
   is the test of whether that bias matters for HDX.
8. **After T2b, each ladder tilts one amide's opening excess (T2c)** (2026-11-07). T2b failed calibration (I61 0.47
   of plain) and support; its rule closed Phase 3 at plain resolution, and the user chose to go on ("Propose next boltz
   job for me and run it"). T2b's top rung opened about 7.5 NHs, the easiest ones, so a tilt on N does not reach the
   core amides; section 2.4's remedy for that, written before T2b's data, is a coordinate aimed at one amide. T2c
   asks E1's question at a decade's precision, the scale of the gap it measures; T2b's failure at a factor of 2 stands.
9. **After T2c, Boltz-2's exact plain samples bound the four core amides (T2d)** (2026-11-08). T2c failed calibration
   (I61 not stationary) and the user chose to go on ("Propose next boltz job and submit it for me"). Its diagnosis:
   the target ladders open for the first time at sweeps 617-1766 (seed 7 never opens V70, W45 or L67), so the read half
   holds each ladder's one relaxation; I61's tilted rungs rise between the halves at all 14 rungs in all four jobs, and
   its lambda 0 reads 0.44 of plain as T2b's read 0.47, the move's drift toward protection that T1-T1d measured. A
   move that walks carries that bias and the moves that keep the distribution are fresh draws (T2a), so the decade
   question goes to Boltz-2's exact samples, in the direction they can settle: more samples bound an amide that never
   opens. T2c's failure stands; no ladder variant is run.
10. **After T2d, Boltz-2's own order on the add-on's domain class (E2a)** (2026-11-09). T2d read the four core
   amides alone (V70 too stable by more than a decade; I13, W45 and L67 under experiment, unresolved at a decade),
   and the user asked for the next job ("propose and run next BOltz job"). Every sampler past plain has failed
   (E1's weights, T1-T2c's moves), so the add-on would work on plain samples, keeping Boltz-2's order inside each
   element. E2a tests that order on Ferrari et al.'s 11 new site-resolved domains before any library or training.

## 4. Phases

### Phase 0: setup and first test on midway3

Cluster rules are in `~/.claude/midway3.md`. Slurm flags go on the command line:
`--partition=gpu --qos=gpu --account=pi-trsosnic --gres=gpu:N`, with `--time` under the QOS's 36 h cap.

- [x] midway2 cannot host it (2026-10-21): of its GPU partitions `pi-trsosnic` reaches only `gpu2`, five K80 nodes
      (Kepler, compute capability 3.7, which CUDA 12 does not support) and one V100 node, midway2-0664,
      `DOWN+NOT_RESPONDING` since 2024-06-07; CentOS 7.9, glibc 2.17.
- [x] midway3 `gpu` partition (2026-10-21): 4 GPUs per node on 5 V100 nodes, 4 RTX 6000 nodes (one down) and one A100
      node (midway3-0294). QOS `gpu`: at most 36 h, 16 GPUs and 4 nodes per user. RHEL 8.10, glibc 2.28.
- [x] Working directory `/beagle3/trsosnic/yinhan/boltz_hdx` (`$B` below), not `/project`: a venv is tens of
      thousands of files and `/project` metadata is 60x slower (findings 10.14). The project path, if needed, is
      `/project/trsosnic/yinhan/` (`~/.claude/midway3.md`'s `/project/trsosnic/yinhanw/` exists but is empty).
- [x] Boltz cloned at `b1ebfc4` into `$B/boltz`; `$B/.venv` from module `python/3.11.9` (`$B/install.sh`, login node)
      with `torch 2.14.1+cu126` from the PyTorch CUDA 12.6 index, then `pip install -e "./boltz[cuda]"`. The PyPI torch
      is a CUDA 13 build without V100 support (findings 10.14). Weights and CCD in `BOLTZ_CACHE=$B/cache`, downloaded
      on the login node by `$B/prep.py`.
- [x] MSAs from the ColabFold server, fetched on the login node by `$B/prep.py` through `boltz.main.compute_msa` (the
      call `--use_msa_server` makes): ubiquitin 16,384 sequences; EHEE_rd2_0005, de novo, its own sequence only, so
      it runs effectively single-sequence. Compute nodes reach RCSB and Hugging Face (HTTP 200); ColabFold's root URL
      gave 404, which shows the host answers but not that MSA requests work, so MSAs stay a login-node step.
- [x] First test, 2026-10-21: one 1-GPU job per GPU type (output `$B/runs/boltz_t0_<gpu>_<jobid>/`, logs `$B/logs/`,
      script `$B/first_test.sbatch`). Peng 2022's two HDX proteins, ubiquitin (1UBQ, 76 residues) and EHEE_rd2_0005
      (40 residues, de novo), each as a default prediction (5 samples, `step_scale` 1.5, no tag) and an `md`-tagged
      ensemble (100 samples, `step_scale` 1.0). `$B/check_samples.py` measures against the 1UBQ crystal and the repo's
      `example/04.HDX/pdb/EHEE_rd2_0005.pdb` and saves each sample's N/CA/C to `<run>.npz`.
  * Hardware: A100-PCIE 40 GB (compute capability 8.0), V100-PCIE 16 GB (7.0), Quadro RTX 6000 24 GB (7.5), driver
    535.216.03. All three ran bf16 with the cuEquivariance kernels. Peak GPU memory 2.4 to 5.0 GB. With one seed the
    default samples agree across GPU types within 0.01 A in every CA metric.
  * Speed (prediction time, trunk included, 10 samples in parallel): 100 `md` samples of ubiquitin took 93 s on the
    A100, 117 s on the V100 and 122 s on the RTX 6000, about 1 s per sample; EHEE 54 to 71 s. So 10^3 ubiquitin
    samples cost about 20 V100-minutes and 10^5 about 33 V100-hours.
  * Default samples are near-native with clean geometry: CA-RMSD to the reference median 0.78 A (ubiquitin) and
    0.93 A (EHEE); peptide C-N 1.327 +- 0.008 A; at most one heavy-atom pair under 2.5 A per sample.
  * `md` samples are wider, with worse local geometry (bond spreads in the table below): ubiquitin pairwise CA-RMSD
    mean 2.0 A (max 4.7), CA-RMSD to 1UBQ median 1.6 A (max 3.0), CA RMSF up to 7.7 A at Gly76. Peptide C-N sd 0.054
    to 0.059 A against 0.014 in crystal structures (Engh and Huber 1991); about 7% of C-N bonds outside 1.25 to 1.45 A,
    along the whole chain; on the V100, 2 of 100 ubiquitin samples have a C-N bond outside 1.05 to 1.65 A. 86 to 97 of
    every 100 samples have a heavy-atom pair under 2.5 A, up to 10 pairs.
- [x] Tag by step-scale test, 2026-10-21 (`$B/tag_stepscale_test.sbatch`, output
      `$B/runs/boltz_tag_stepscale_v100_61231966/`): {no tag, `md`} x {`step_scale` 1.5, 1.0}, 100 samples each at
      seed 2, both proteins; close contacts classified by atom pair locally.

  | cell | C-N sd (A), ubq / EHEE | N-CA sd | pairwise CA-RMSD (A), ubq / EHEE | samples with a contact < 2.0 A, ubq / EHEE |
  | --- | --- | --- | --- | --- |
  | no tag, 1.5 | 0.008 / 0.009 | 0.004 | 1.28 / 0.26 | 3 / 2 |
  | `md`, 1.5 | 0.008 / 0.009 | 0.004 | 1.57 / 0.45 | 1 / 5 |
  | no tag, 1.0 | 0.031 / 0.036 | 0.018 | 1.61 / 0.52 | 15 / 17 |
  | `md`, 1.0 | 0.055 / 0.054 | 0.035 | 1.99 / 0.88 | 25 / 27 |

  * `step_scale` 1.0 causes the bond-length scatter: at 1.5 both tags give clean bonds, at 1.0 the `md` tag nearly
    doubles it. The `md` cell at 1.0 reproduces the first test with a new seed.
  * The `md` tag widens the ensemble at both step scales; `md` at 1.5 is almost as wide as no tag at 1.0, with clean
    bonds.
  * Most contacts under 2.5 A are polar pairs (Lys NZ or Arg NH with Asp/Glu carboxylate O, and backbone N-O): salt
    bridges and H-bonds drawn 0.3 A or more too short (N-O is about 2.8 A or longer). Neither reference structure has
    a contact under 2.5 A. One steric overlap stands out: in no-tag, 1.5 ubiquitin model 33 the Arg42 and Arg72
    guanidinium groups sit on top of each other (CZ-CZ 0.69 A, 15 contacts); the same pair comes within 1.75 A in
    no-tag, 1.0 model 41.
  * The conflict: the unscaled sampler (1.0) is the one whose samples follow the model's distribution, and it has
    broken local geometry; `step_scale` 1.5 has clean geometry but sharpens the distribution toward its main state,
    which biases the populations HDX reads toward protection.
- [x] Phase 1 sampler setting: `md` at `step_scale` 1.0, because its bond-length scatter does not reach the protection
      flags (2026-10-21, local). Every set from the tag by step-scale test went through `py/get_protection_state.py`
      unchanged, with configs from `example/04.HDX/1.config.py`'s recipe (ff_2.1, heavy-atom coverage) built from the
      two reference structures. The engine comes from a scratch tree (py/ copied and obj/ built at repo `6d295b0`),
      because the repo's `obj/` was older than its `py/` and lacked `get_param_size`. Files in
      `scratchpad/boltz_hdx/flags/` (`run_flags.sh`, `samples_to_traj.py`, `compare_flags.py`). The estimator reads
      only N, CA and C (Upside places its own side-chain beads from the backbone), so Boltz's side chains, with their
      short salt bridges and the Arg42/Arg72 overlap, do not enter.
  * Geometry noise alone: the `md`, 1.5 samples plus isotropic noise per atom (0.038 A, which gives C-N sd 0.054 A,
    the 1.0 spread) flip 0.2% (ubiquitin) and 1.3% (EHEE) of the combined flags on the same frames; the H-bond flag
    alone flips 2.8% and 0.8%, mostly in amides that burial already protects. Per-residue protected fractions move by
    at most 0.03, mean 0.002, below the binomial error of 100 samples, 0.02 to 0.05.
  * `md` at 1.0 against `md` at 1.5, the whole difference: mean |dp| 0.019 (ubiquitin) and 0.016 (EHEE), max 0.16 and
    0.12. At 1.0 fewer amides are never seen open (p = 1): 45 against 55 of 72 for ubiquitin and 24 against 29 of 39
    for EHEE.
  * Every set stays near the native flags: the reference structures protect 59 of 72 and 32 of 39 amides, and 96 to
    98% of sample flags agree with them. With 100 samples, every amide protected in the native state is censored at
    RT ln(100) = 2.7 kcal/mol.

### Phase 1: independent Boltz-2 ensembles, no training

HDX:
- [ ] `md`-tagged ensembles at `step_scale` 1.0 (Phase 0) for ubiquitin and EHEE_rd2_0005, then other soluble
      proteins with HDX data and Upside results. Done: 10^3 samples per protein at seed 3 (`$B/phase1_hdx.sbatch`,
      `$B/phase1_ubq_rerun.sbatch`), and the sets below.
  * Chunking: Boltz split the samples into `multiplicity % max_parallel_samples + 1` chunks (`diffusionv2.py:383`),
    one chunk for 1,000 and for 10, so `--max_parallel_samples` had no effect and 1,000 ubiquitin samples ran out of
    memory on the 16 GB card. `$B/boltz` now splits into chunks of `max_parallel_samples` (uncommitted edit in our
    clone). Chunking leaves the samples unchanged: the noise and random rotations are drawn for every sample before
    the chunk loop; of EHEE's 100 samples at seed 3 two unchunked runs agree exactly at PDB precision, and chunks of 10
    move no backbone atom by more than 0.001 A, with the same confidence order (`$B/compare_chunks.py`). Ubiquitin's
    1,000 in chunks of 100 took 20 V100-minutes, peak 8.2 GB.
  * The 1,000-sample sets match Phase 0's 100. EHEE: CA-RMSD median 0.95 A, Rg 9.22 A, C-N 1.337 +- 0.053 A, 6.7% of
    C-N outside 1.25-1.45 A (7.2%). Ubiquitin: CA-RMSD to 1UBQ median 1.49 A (1.58 and 1.62), pairwise 1.88 A (1.97,
    1.99), Rg 11.35 +- 0.16 A (reference 11.49), C-N 1.343 +- 0.055 A with 6.8% outside 1.25-1.45 A (7.1%, 6.9%), CA
    RMSF up to 7.3 A at Gly76. 17 of 1,000 ubiquitin samples have a C-N bond outside 1.05 to 1.65 A (Phase 0: 2 and 1
    of 100), and 2 have one over 2 A, a broken chain, both in the C-terminal tail (70-71 at 2.05 A, 75-76 at 3.71 A;
    Phase 0 had 74-75 at 2.89 A). The samples stay unfiltered with uniform weights.
- [x] Protection flags from every sample's N/CA/C through the repo's protection-state code, uniform weights; dG per
      residue against experiment and against Upside. An amide with p = 1 is censored at RT ln(N), not a value.
      Pipeline `scratchpad/boltz_hdx/flags/` (`run_flags.sh`, `dg_phase1.py`, one `readout()` for every set); each
      rerun has reproduced the earlier flags, dG tables and profiles byte for byte. A 10^5 set is its 10 job sets
      copied (md5-equal) and pooled by `flags/pool_seeds.py` (the 10-25 sets in seed order, later ones in ledger
      order) into `samples/<set>_n100000.npz`. One figure per system, `dg/<name>_<set>_dG_vs_residue.png`, from
      `example/00.AnalysisScripts/helpers/plot_ref_style.py` (findings 10.36) on `dg/<name>_<set>_dG_profiles.npz`,
      the pipeline's layout with the experiment's temperature in K as its one temperature. Comparisons with experiment
      and Upside are in `dg/dg_phase1.log` and below, not in the figures.
  * References: Peng 2022 SI Tables S5 and S7 (NMR HX, 0 M urea) and S8 and S10 (Upside FF2), in
    `flags/reference/`; findings 11b. RT is the experiment's, 298 K for EHEE and 273 K for ubiquitin, so the bound
    RT ln(1000) is 4.09 and 3.75 kcal/mol. The ubiquitin data are for Ub psWT (F45W/H68N).
  * 10^3 readout (2026-10-22). Protected fractions agree with Phase 0's independent 100 samples (Fisher's exact test,
    smallest p 0.011 over 72 ubiquitin donors, 0.127 over 39 EHEE donors). Boltz-2 never opens 34 of 72 ubiquitin and
    21 of 39 EHEE donors. Every measured amide but EHEE's R40 (4.02) lies above the bound (4.9 to 8.5 kcal/mol), so
    these sets show only whether Boltz opens amides experiment shows are stable. With experiment and Upside censored at
    the same bound, Boltz-2 falls below it at 2 of 39 measured amides (ubiquitin R54 and Y59, with 1 and 2 open
    samples of 1,000), and Upside FF2 at 11 (EHEE F6 Y24 V37 R40; ubiquitin L50 R54 T55 L56 Y59 L67 V70).
  * Rank order cannot be tested at 10^3: Upside's correlation with experiment (r 0.63 and 0.46) comes from values
    above 4 kcal/mol, where Boltz resolves nothing. Reaching 8.5 kcal/mol needs about 1.7 x 10^6 samples at 298 K and
    6.4 x 10^6 at 273 K; 10^4 samples reach 5.5 and 5.0.
  * Over-protection cannot be tested: Boltz never opens 6 EHEE amides (Y4 D8 E11 E18 A20 H27, all with assigned NMR
    peaks) and 12 ubiquitin amides that have no reported dG, and the paper does not say why they are missing.
- [x] Ubiquitin psWT (F45W/H68N) ensemble, so the comparison is like for like: 1,000 `md` samples at seed 3 with its
      own ColabFold MSA. Reference: 1UBQ with residues 45 and 68 renamed and cut to backbone and CB (no psWT structure
      is used anywhere in Peng 2022, whose Upside runs also start from 1UBQ). CA-RMSD to the reference median 1.56 A,
      pairwise 1.90 A, Rg 11.34 +- 0.16 A, C-N 1.343 +- 0.055 A, the wild type's statistics (1.49, 1.88, 11.35,
      1.343 +- 0.055). 36 of 72 donors never open (wild type 34), 2 never protected; of the 24 measured amides Boltz-2
      falls below min(experiment, bound) only at L56, by less than 0.001 kcal/mol (one open sample of 1,000), where
      the wild type fell below at R54 and Y59; Upside FF2 falls below at 7. The mutations at 45 and 68 do not change
      the picture.
- [x] Rank order on HEEH_rd4_0097, the third Peng 2022 protein (decided 2026-10-22, delegated by the user). Its 20
      measured amides lie at 2.2 to 4.3 kcal/mol (Table S6; HX at 278 K, pD_read 4.6), so 10^4 samples (bound
      5.09 kcal/mol) resolve all of them, where ubiquitin and EHEE at 10^6 would cost hundreds of V100-hours for
      8.5 kcal/mol. Upside FF2's agreement is published (Table S14, 0 M: Pearson 0.185, Spearman 0.317, RMSE
      1.06 kcal/mol; its per-residue Table S9 is a copy error). Sequence: the 43-residue design that Table S6 numbers
      from D1; reference the design model Peng simulated, `/project2/trsosnic/share/HEEH_rd4_0097.pdb`. A run's sample
      count is bounded by GPU memory, because Boltz's confidence step keeps every sample's PAE and PDE logits on the GPU
      and concatenates them (about 2 x N x L^2 x 256 bytes plus the copy), so the 10^4 are 10 runs of 1,000 at seeds
      3-12 sharing one input, pooled into one set (`$B/phase1_heeh_seeds.sbatch`; 719-6,991 s per seed). The pooled
      10^4: CA-RMSD to the design model median 1.70 A (max 10.46), pairwise 1.10 A, Rg 9.47 +- 0.12 A (model 9.27),
      C-N 1.340 +- 0.048 A, the other `md` sets' spread. Readout 2026-10-24 at 278 K, 40 donors (`dg/dg_phase1.log`):
      * Sampling is not the limit: the two halves of the seeds agree donor by donor (smallest Fisher p 0.125, mean
        |dp| 0.001), and their dG at the measured amides differ by 0.40 kcal/mol rms (Pearson 0.88 between halves).
      * 16 of the 20 measured amides resolve, at 1.49-5.09 kcal/mol against experiment's 2.16-4.28. The other 4 are
        never open: V12, L13, V35 and L39 at >= 5.09, where experiment has 3.75, 3.15, 3.84 and 3.68.
      * No rank order: Spearman 0.042 over all 20 (censored ranked at the bound; p 0.86; bootstrap over the 10 seeds
        -0.055 to 0.166), 0.092 over the 16 resolved (Pearson 0.214, RMSE 0.95, mean offset +0.14 kcal/mol); RMSE over
        all 20 at least 1.09. At 20 amides neither Boltz-2's Spearman nor Upside FF2's 0.317 differs from zero or from
        the other.
      * Where it goes wrong: the hairpin's V19 and W21, experiment's least stable amides (2.16, 2.34), read 3.94 and
        4.32; T24 (2.54) reads 1.49 and E4 (4.18) 2.68. Boltz-2's protected hairpin amides (V19, W21, Y26, D28: 3.9 to
        >= 5.09) sit as high as the helices', where experiment puts the hairpin lowest.
      * No unfolded state: no sample has more than 35% of its donors open (median 17.5%, the loops), and 6 of 10^4 lie
        beyond 8 A of the design model. Open fractions against experiment's exp(-dG/RT): the hairpin's V19 and W21
        0.0008 and 0.0004 against 0.020 and 0.014, T24 0.063 against 0.010, the helices 0 to 0.008 against 0.0004 to
        0.0047. The helices' amides sit close together at 3.5-4.3 kcal/mol, as exchange through whole-protein unfolding
        would put them (confirmed by the urea item below).
      * For Gate 1 this is the "too narrow" outcome on HDX (censored dG where experiment is finite), with some amides
        too open as well.
- [x] Local exchange (user 2026-10-24, after HEEH): does Boltz-2 rank amides whose exchange comes from fluctuations of
      the folded state? EHEE_rd2_0005 and Ub psWT at 10^5 `md` samples each (bound 6.82 kcal/mol at 298 K and 6.25 at
      273 K), which resolves 9 of EHEE's 15 and 12 of Ub psWT's 24 measured amides (4.0-6.8 kcal/mol): 100 runs of
      1,000 per protein at seeds 13-112, 10 V100 jobs each (`$B/phase1_seeds.sbatch`), about 20 and 35 V100-hours.
      Readout fixed before the data: Spearman against experiment over the resolved amides and with censored ones ranked
      at the bound, the halves of the seeds against each other, and the open fractions against exp(-dG/RT). If Boltz-2
      ranks them, the no-training route holds where exchange is local; if not, it is closed and the training question
      is next.
      * Inodes: Boltz writes 5 files per sample, and the first submission filled `/beagle3`'s group inode quota
        (findings 10.39). `pack_predictions.sh` packs each `predictions/<name>` directory into `<name>.tar` and removes
        it once the tar lists every file; `check_samples.py` reads the models from the tars (on Phase 1's packed EHEE
        set it writes the old reader's npz byte for byte); `phase1_seeds.sbatch` packs each seed as it ends and keeps a
        seed whose tar holds 1,000 models. Submission goes through the ledger `$B/submitted_1e5.txt`.
      * Readout 2026-10-25, sets `samples/{ehee,ubqps}_md_s10_n100000.npz`, figures
        `dg/{ehee,ubqps}_md_s10_n100000_dG_vs_residue.png`.
        - The pooled sets match the 10^3 ones: CA-RMSD to the reference median 0.94 and 1.52 A, Rg 9.22 and 11.34 A,
          C-N 1.338 and 1.343 +- 0.05 A; no sample repeats; 47 and 131 of 10^5 have a C-N bond over 2 A.
        - Sampling is not the limit. The halves (seeds 13-62 against 63-112) agree donor by donor (smallest Fisher p
          0.019 and 0.064). The measured amides open in few samples, EHEE's in 0-63 of 10^5 and Ub psWT's in 0-80, most
          in 1-10, so a value carries counting noise of about RT/sqrt(k) (0.4 kcal/mol at k = 2); EHEE's halves differ
          by 0.46 kcal/mol rms at its measured amides, Ub psWT's by 0.12.
        - The openings are not broken chains: of the samples that open a measured amide, 4 of 134 (EHEE) and 0 of 109
          (Ub psWT) have a C-N bond over 2 A. Many lie 4-7 A from the reference, and Ub psWT's I3, F4 and V5 open
          together in the same 2 samples.
        - EHEE_rd2_0005 (298 K, bound 6.82 kcal/mol) is ranked. 14 of 15 measured amides resolve, at 4.36-6.82 against
          experiment's 4.02-8.41 (A15 never opens). Spearman 0.735 (p 0.003) over the 14 and 0.741 over all 15 with
          experiment, Boltz-2 and Upside censored at the bound (seed bootstrap 95% 0.51-0.86); Upside FF2 on the same
          comparison 0.609 (p 0.02), below Boltz-2 in 84% of the resamples, so the two do not differ. Pearson 0.736;
          Boltz-2 sits 0.91 kcal/mol below experiment on average, its open fractions mostly 2-36x experiment's (W19
          equal, R40 4x lower).
        - Ub psWT (273 K, bound 6.25) is not ranked: Spearman 0.173 (p 0.57) over the 13 resolved, -0.206 over all 24
          censored (bootstrap -0.39 to 0.10); Upside FF2 0.566 and 0.432 (p 0.03). Boltz-2 never opens 11 measured
          amides, among them four of experiment's six least stable, V70, W45, I13 and L67 (4.88-5.27 kcal/mol, open
          fractions 6e-5 to 1.2e-4), and opens D21 (7.73) at 1.2e-4, 180x experiment, and strand 1's I3, F4 and V5
          (about 8) at 2e-5, 50x. With the halves this close, the misplaced openings are the model's ensemble.
        - For Gate 1: the no-training route ranks one of the two proteins, as well as Upside FF2, and fails on the
          other, where it opens the wrong parts.
- [x] MSA control on Ub psWT (user 2026-10-25, "propose and apply fix", chose the controlled 10^5 tests). Diagnosis of
      the 10^5 sets first, from their raw flags:
      * Protection is burial-dominated. At the measured amides the backbone H-bond breaks in 1-95% of samples (Ub
        psWT's N25 95%, T55 all), but Upside's burial term (BurialLevel > 5; medians 7-22) keeps the amide protected,
        so an amide opens only in the rare samples where it is both unbonded and exposed, the 4-7 A excursions above.
        The dG therefore measure how often a model samples partially unfolded states.
      * The seeds are exchangeable: the per-seed dispersion index of each donor's H-bond and protection counts is
        0.95-1.06 (median; at most 1.4, none significant at p 0.001) for both proteins, so the 1,024 MSA rows Boltz
        subsamples per run do not move the ensemble.
      * The errors are in the model's distribution, and HEEH's missing unfolded state cannot be fixed by a setting;
        Boltz-2 training code is not released and mdCATH holds no unfolding. One input differed in kind between the
        protein Boltz-2 ranks and the one it fails: EHEE's ColabFold MSA is the query alone (HEEH's three sequences),
        Ub psWT's 16,383 sequences.
      The test: Ub psWT at 10^5 in EHEE's regime, `inputs/ubqps_nomsa.yaml` with `inputs/msa/ubqps_nomsa_0.csv`
      (EHEE's form, the header and the query twice), at the MSA set's seeds 13-112, 10 jobs, about 34 V100-hours.
      Gate before the other nine jobs, fixed in advance: job 1's samples must be a folded ubiquitin (CA-RMSD to the
      reference median below 4 A, Rg within 15% of 11.49 A), or the answer is that Boltz-2 needs the MSA to fold it.
      Passed 2026-10-26 (median 1.54 A, Rg 11.34 A, the MSA set's 1.52 A and 11.34 A). Readout fixed before the data:
      if the censored Spearman's block-bootstrap 95% interval excludes the MSA set's -0.21 and reaches Upside FF2's
      0.43, MSA depth is the defect on Ub psWT and Phase 1 takes the query-only regime for every protein; otherwise MSA
      depth is ruled out and the no-training route stays closed on Ub psWT. HEEH's missing unfolded state stands either
      way.
      * Readout 2026-10-27, set `samples/ubqps_nomsa_md_s10_n100000.npz` (EHEE's and Ub psWT's MSA sets rebuilt from
        their jobs the same way, and the scratch engine rebuilt at 6d295b0, reproduce the whole 10-25 `dg_phase1.log`,
        tables and profiles byte for byte), figure `dg/ubqps_nomsa_md_s10_n100000_dG_vs_residue.png`.
        - The set is the MSA set's folded ubiquitin: CA-RMSD to the reference median 1.54 A (MSA 1.52), 7 of 10^5
          beyond 8 A (none), Rg 11.34 +- 0.17 A, C-N 1.342 +- 0.054 A, 124 samples with a C-N bond over 2 A (131); no
          sample repeats. The halves agree (smallest Fisher p 0.057; 0.32 kcal/mol rms at 15 measured amides).
        - Without the MSA, Boltz-2 opens Ub psWT more often. 4 donors never open (MSA 14), 21 of the 24 measured amides
          resolve (13), and the measured amides open 201 times against 119, a third of them in samples beyond 4 A of
          the reference (8% with the MSA). No sample has more than a third of its donors open, so these are partial
          excursions, not an unfolded state.
        - The added openings are not in experiment's order: Spearman -0.196 (p 0.40) over the 21 resolved, Pearson
          -0.109, 1.01 kcal/mol below experiment on average; over all 24 censored at the bound -0.062 (bootstrap 95%
          -0.245 to 0.089), against Upside FF2's 0.432. Experiment's least stable V70, W45 and I13 open at 1-2e-5
          against 0.7-1.2e-4, while strand 1's I3, F4 and V5 open at 1-1.6e-4, about 300x experiment, D21 at 1.4e-4
          (200x) and Y59 at 8.5e-4 (20x).
        - By the rule fixed before the data: the interval contains the MSA set's -0.21 and does not reach 0.43, so MSA
          depth is ruled out as the defect on Ub psWT and the no-training route stays closed there.
- [x] Which opening each measured amide exchanges through (2026-10-27): experiment's dG against urea (Peng 2022 SI
      Tables S5-S7), each amide's slope as a fraction of Table S4's global m-value (EHEE 0.607, HEEH 0.871, Ub 1.270
      kcal/mol/M, fitted on the most stable amides), near zero for a local opening and near one for global unfolding
      (Peng 2022, p. 554). `flags/urea_classes.py`, output `dg/urea_classes.log`; the slope is a straight line over
      the reported concentrations, where the paper's per-residue values come from a two-state fit.
      * Ub psWT exchanges locally below the bound. V70, I13, L69, Y59, W45 and L67 (4.9-5.5 kcal/mol) have 0.01-0.23 of
        m_global; the amides at 7.7-8.5 (I3, F4, V5, L15, D21, I30, K27) 0.70-0.92, global (L15 and K27 are Table S4's
        pair). Boltz-2 never opens four of the local exchangers (V70, I13, W45, L67) and opens the global ones in a few
        samples (I3, F4, V5 and L15 in 2 of 10^5, D21 in 12). Its error grows the more local the opening (Spearman of
        error against m fraction -0.43 over 23 amides; L50's slope rests on two concentrations and is not read):
        Boltz-2's folded ubiquitin lacks the real one's local fluctuations, and its rare excursions open strands that
        experiment opens only on global unfolding.
      * No measured EHEE amide is local. The 9 below the bound have 0.29-0.62 of m_global, subglobal openings; A12,
        R13, E36 and T29 0.94-1.00, global (A12 and R13 are Table S4's). Boltz-2's ranking of EHEE is a ranking of its
        subglobal openings, and its error does not follow the class (+0.20).
      * HEEH exchanges by global unfolding in its helices (16 of 17 helix amides at 0.72-1.10 of m_global; Table S4
        fits m_global on 8 of them) and by a subglobal opening in its hairpin (V19 and W21 0.63, T24 0.46). Boltz-2,
        with no unfolded state, spreads the helix amides from 2.7 to >= 5.09 kcal/mol where global unfolding holds
        experiment's at 3.2-4.3, and leaves V19 and W21 1.8-2.0 kcal/mol too protected (error against class -0.12).
      * The 10^5 test was framed as local exchange on both proteins, but EHEE's is subglobal: Boltz-2 ranks one
        protein's subglobal openings and fails on another's local ones and a third's global and subglobal ones.
- [x] m-values, Peng 2022's second metric (user 2026-11-03, the PI: "it is not only delta G we care about. We also care
      about the slopes"). Peng's Figs 5-7 and 10 and Table S14 score every amide's dG and its m-value, the slope of dG
      against urea, at 0 M and at high urea; the readouts above scored dG only. Peng's denaturant model weights a
      conformation with N protected NHs by exp(-s N [den] / RT) (main text p. E; `example/04.HDX/4.calc_HDX.py`), so
      at 0 M m_i = s dN_i, where dN_i is the number of protected NHs lost in the samples that open amide i: the size of
      the opening. Pearson and Spearman of m at 0 M do not depend on s. Experiment's m is the 0 M slope of eq. S18
      fitted to Tables S5-S7 with Table S4's global channel; scoring FF2's Tables S8 and S10 against it reproduces
      Table S14's FF2 m-value Pearson and Spearman to within 0.03 (EHEE at 0 and 4 M, Ub at 0 M), though not Ub's RMSE
      (0.35 against 0.52). Code `flags/urea/` (`parse_si.py` writes `reference/*_urea.csv` from the SI text,
      `check_s14.py`, `mvalue.py`); log `urea/mvalue.log`.
      * Boltz-2 opens every measured amide locally. Its openings break 0.9-6.1 of a sample's 31-59 protected NHs (one
        Ub psWT sample 10.6), including the amides that experiment's slopes place at 0.85-1.0 of m_global (EHEE's A12,
        R13, T29 and E36; HEEH's helices). Its EHEE dG ranking (Spearman 0.74) comes from openings of the wrong size:
        experiment exchanges A12 and R13 by global unfolding, Boltz-2 through openings of 5.4 and 0.9 NHs.
      * HEEH_rd4_0097: Boltz-2's m does not follow experiment (Spearman -0.18 over 16, bootstrap over the 10 seeds
        -0.36 to 0.09); Upside FF2 0.475 (Table S14).
      * EHEE and Ub psWT are not resolved at 10^5. An amide's m is the mean over the samples that open it, and most
        measured amides open in 1-9 (Ub psWT median 1). EHEE Spearman 0.37 over 14 (-0.34 to 0.72; FF2 0.41 on the
        same 14); Ub psWT 0.33 over 12 (-0.03 to 0.67; FF2 0.65 on the same 12), where I3, F4, V5, L15 and T55 share
        one value from the same 2 samples.
      * The s-dependent half of Table S14 (m RMSE; dG and m at 4, 2.13 and 2.4 M) is not computed. Peng calibrates s
        on the model's own global unfolding (Table S4), and Boltz-2 samples none: the fewest protected NHs in any
        sample are 21 of 39 (EHEE), 26 of 40 (HEEH) and 48 of 72 (Ub psWT).

mdCATH 320 K (CHARMM22*/TIP3P, 5 replicas of about 500 frames at 1 ns), equilibrium only; not run:
- [ ] Test set: 8 to 12 domains from the 89-domain panel (`ff3_selection/` on the cluster; `plan.md` Phase 19), over
      the 50-100, 101-150 and 151-250 residue bins, plus 3 to 5 domains that switch states at 320 K (multimodal Calpha-
      RMSD or TICA distribution across replicas). Every domain must come from Boltz-2's held-out mdCATH clusters, since
      its training included mdCATH; if that split cannot be found, mark the results as possibly optimistic.
- [ ] The panel's local copies keep backbone atoms only. Fetch heavy atoms for the chosen domains, reusing the
      range-request approach of `mdcath_fetch.py`.
- [ ] Metrics on Calpha or heavy atoms against the 5 MD replicas: per-residue RMSF; JS divergence of pairwise-distance
      and Calpha-RMSD distributions; TICA fitted on MD with Boltz samples projected onto it, and the free-energy surface
      on TIC1 and TIC2.

Gate 1:
* dG tracks experiment at least as well as Upside, and the mdCATH distributions match: equilibrium needs no training.
* Ensembles too narrow (censored dG where experiment is finite, RMSF below MD): Phase 2.
* Ensembles in the wrong states, not only too narrow: reconsider the approach before training.
* **Verdict on the HDX half (2026-10-27): wrong states.** Boltz-2 ranks EHEE, but HEEH lacks its global and subglobal
  openings (too narrow) and Ub psWT its local fluctuations, while opening strands that exchange only globally (the
  urea item above). By this rule the approach was reconsidered before any training, and the mdCATH half waits on that.
  The user's choice, 2026-10-27: no retraining of Boltz-2, a small trained add-on on the frozen model (Phase 2).

### Phase 2: a trained add-on that tilts frozen Boltz-2's ensemble (2026-10-27)

User 2026-10-27: "You are not allowed to retrain boltz, as we do not have the computation resources for that. But you
can train small models as addons to make it work." Plan approved the same day.

Why this form. Mixing Boltz-2's ensembles with experiment's own global unfolding population (Table S4's dG_global, the
best any global add-on could do) leaves every Spearman as it was (EHEE 0.74, HEEH 0.04, Ub psWT -0.19): an unfolded
state adds the same open fraction to every amide. The errors are site-specific inside the folded ensemble (Phase 1's
urea item), so the add-on has to change where Boltz-2 opens, during sampling. What-if substitutions (rescored in
memory, about +-0.05): Ub psWT reaches about 0.48 if V70, I13, W45 and L67 open at experiment's values and about 0.73
if strand 1 and D21 are also suppressed; HEEH 0.32 with its helix amides fixed and 0.70 with the hairpin too.

* Target p(x) proportional to p_Boltz(x) exp(-U_phi(x; T)), sampled by Boltz-2's own steering with its weights frozen:
  `model/potentials/potentials.py` (`Potential`, `get_potentials`) and `diffusionv2.py` `AtomDiffusion.sample`
  (gradient guidance on the denoised x0; Feynman-Kac resampling with the guided/unguided correction `ll_difference`).
* U_phi = -sum_s mu_s O_s(x) - sum_i lambda_i o_i(x). o_i is a differentiable per-amide open score from N, CA and C, a
  small surrogate trained to reproduce `get_protection_state.py`'s H-bond and burial flags on non-test frames; O_s is
  its mean over secondary-structure element s. mu and lambda come from a small network on the frozen `s_trunk` and the
  temperature, lambda_i shrunk toward zero so the order inside an element stays Boltz-2's, which is what ranks EHEE.
  This is the maximum-entropy form: the smallest change to Boltz-2's distribution that moves the open fractions.
* A global term: a small network for dG_unf, entering as p_open,i = (1 - f_U) p_open,i(folded) + f_U. It sets absolute
  values and the crossover to global exchange, never the order.
* Two sampler changes, in our wrapper (subclass or patch at import; nothing in Boltz's tree or a venv): every final
  particle kept with its weight, since the sampler draws one per group at the last step (`torch.multinomial(..., 1)`),
  which biases events near 1e-4; and guidance only where the reverse step is stochastic, since noise is off below
  `gamma_min` and `ll_difference` does not apply there.
* Readout: weighted flags, p a weighted mean and the bound RT ln(N_eff). `dg_phase1.readout()` gains weights; uniform
  weights must reproduce its present output byte for byte.
* Training signal: Ferrari et al. (Rocklin lab; bioRxiv 2025.03.20.644235, PMC11974690, published in Nature),
  multiplexed HDX-MS of 5,778 domains of 28-64 residues (3,590 stable) at 20 C, pH 6 and 9, each with an unassigned set
  of opening free energies, 13 also with site-resolved HDX-NMR. Data through the authors' form
  (https://forms.gle/RwJwvfw6WN4gjXaD9); raw spectra PRIDE PXD061702. Experiment is the only teacher whose ceiling is
  above Upside: Upside REMD caps the add-on at Upside, and mdCATH's 1 ns frames cannot resolve 1e-4 openings.
* Loss per domain: the 1-D Wasserstein distance between the estimator's dG set and the measured set, plus dG_unf
  regression, on a fixed library of 10^3 unguided and 10^3 element-steered samples per domain with path weights and
  precomputed flags, so no gradient passes through Boltz. Held out: the three test proteins, the EHEE and HEEH
  topologies and beta-grasp domains; the 13 NMR domains validate the residue assignment.
* Code in `scratchpad/boltz_hdx/addon/` (git-ignored) and `$B/addon/` on midway3.

- [x] E1, the falsifier (about 7 V100-hours, no training; closed 11-03, below): untrained Boltz-2's own rare open fractions by particle
      resampling, one amide at a time, with an opening potential whose bias is undone in the weights; 1,024 particles
      per amide, two seeds.
      * First, the wrapper with the add-on off reproduces plain Boltz-2's samples for the same seed, coordinates equal.
      * Calibration: on amides the unguided sets open at least 10 times (Ub psWT D21 at 12 of 10^5, HEEH T24, the
        others by count) the estimate lies within a factor of 2 of the unguided frequency, with N_eff >= 50.
      * Support: V70, I13, W45 and L67 get open samples with no C-N bond over 2 A and the rest within 3 A CA-RMSD, and
        finite estimates (+-50%).
      * Either failing stops the route: the local class cannot be fixed by a Boltz add-on, and Ub psWT tops out near
        0.20, below Upside's 0.43.
      * Local CPU sandbox, 2026-10-27 (`scratchpad/boltz_hdx/addon/sandbox/`: the midway3 Boltz tree, torch 2.14.1,
        the same weights): the parity check passes (8 samples, seed 7, coordinates equal), and the steered path runs
        end to end. The first opening potential (Boltz's own N and O atoms) left Upside's H-bond at V70 in place
        (0.43-0.98): the estimator infers H and O from N, CA and C, so the potential has to act through that geometry
        and through the coverage burial.
      * The estimator's geometry in torch, 2026-10-28 (user: "run E1"). Every quantity the estimator thresholds is a
        function of N, CA and C alone: the side-chain beads are fixed in each residue's frame and weighted by
        exp(-E(phi, psi)) from a spline table, with no belief propagation. `protection.py` repeats the engine's nodes
        (infer_H_O, protein_hbond, affine_alignment, the bead placement and weights, hbbb_coverage,
        environment_coverage_hb) from the protein's `-HDX.up`, exported to an npz since the midway3 venv has neither
        PyTables nor h5py. `check_protection.py` on 400 frames each of Ub psWT's 10^5 and HEEH's 10^4 sets, half with
        some amide's burial within 2 of the criterion: every node within float32 rounding of the engine (burial within
        4.4e-5), no protection flag differs (0 of 28,800 and 0 of 16,000), and the autograd gradients match central
        differences to 4e-8.
      * The potential (`opening.py`): kappa times the target amide's excess over the three criteria in burial units,
        (burial - 5)+ + 4 (H-bond - 0.01)+ + 4.6 (ASP/GLU coverage - 0.05)+, zero exactly where the estimator calls
        the amide open. Burial carries it because burial is what keeps amides closed: in the 10^5 set an amide whose
        burial is at most 5 is open 91-100% of the time, except where ASP/GLU coverage holds it (Q2, D21, K27, A28).
      * Two sampler changes (`sampler.py`, 2026-10-28; parity unaffected, as both act only when steering). The
        guidance shift enters the weight correction centered, since the step's centering of the coordinates removes
        its mean before the noise is drawn. And before each resampling the log of the particles' mean weight is
        accumulated into log_z, the Feynman-Kac estimate of the tilt's normalizing constant, so the open fraction
        has a second estimate, exp(log_z) times the mean of w times the open flag, that needs no closed particle at
        the end; the self-normalized estimate needs the tilted particles to keep Boltz-2's closed states, which no
        tilt strong enough to open V70 (about 1e-7) can do with 1,024 particles. Fixed before the data: both checks
        read the Feynman-Kac estimate, the same estimator on calibration and support, with N_eff over the open
        particles; "+-50%" means each seed's estimate within 50% of the two seeds' mean. The self-normalized estimate
        is reported alongside.
      * CPU test of the steered path (V70, 4 particles, kappa 1): it runs; the stored final energy equals the one
        recomputed from the final backbone; one particle reaches V70 burial 16.2, beyond the unguided 1% quantile
        (15.7).
      * What the support check is up against, from the 10^5 unguided flags: an amide is open only with H-bond <= 0.01
        and burial <= 5, and V70, I13 and L67 never come below burial 10.2 (0.001% quantile; 0-1 of 10^5 below 9), W45
        below 6.4. Their tails fall a decade per 1.5-2.5 burial units, so Boltz-2's own open fraction there extrapolates
        to about 1e-7 to 1e-8 before the H-bond condition, against experiment's 6e-5 to 1.2e-4. Amides that do open have
        heavy tails (Y59 88 samples below 5, D21 12). Extrapolated, not measured: E1 measures it.
      * GPU pilot 61663101, 2026-10-28 (`remote_jobs.md` §1 line 4): the parity check passes on the V100 (8 samples,
        seed 7, coordinates equal), and 64 unguided samples through the wrapper match the pooled 10^5 flags (852
        openings over the 41 amides the pool opens at least 10 times, against 854.5 expected; no amide off by a
        two-sided binomial test at p < 0.05; longest C-N at most 1.56 A). Readout in
        `scratchpad/boltz_hdx/addon/runs/boltz_addon_e1_pilot_ubqps_v100_61663101/readout.txt`.
      * Tuning 61701305 (`remote_jobs.md` §1 line 4): 256 particles, seed 1, kappa 0.5-4 and guidance weight 0.01
        and 0.05 on V70 and I61. Rule fixed before the data (`addon/e1_tune_ubqps.txt`): E1 uses the setting
        that gives V70 the most open particles among those whose I61 run has a Feynman-Kac N_eff of at least 25;
        agreement with the unguided frequency is not used.
      * E1 proper, at that setting, 1,024 particles, seeds 1 and 2. Calibration: Ub psWT D21 (12 of 10^5 open
        unguided), E34 (35), Q41 (12), Y59 (80) and I61 (38); HEEH T24 (629 of 10^4) and Q5 (36).
        Support: Ub psWT I13, W45, L67, V70.
      * Tuning readout, 2026-10-28 on the Mac Studio (`e1_readout.sh` on the 12 npz; its scratch engine first
        reproduced the pilot's 852 openings against 854.5 expected and C-N at most 1.56 A; output
        `scratchpad/boltz_hdx/addon/runs/boltz_addon_e1_tune_ubqps_v100_61701305/readout.txt`). V70 opens only at
        kappa 0.5 (3 of 256) and kappa 2 with guidance 0.05 (93 of 256); kappa 1 and 4 open none. Every open sample
        has C-N under 2 A and the rest within 2.3-3.7 A CA-RMSD. I61's Feynman-Kac estimate against the unguided
        3.8e-4: 2.1e-5 at kappa 0.5 (2 open), 6.4e-4 and 4.1e-4 at kappa 1 (N_eff 187, 176), then 3.4e-15 and
        5.5e-14 at kappa 2 and 3.0e-37 at kappa 4. The rule picks kappa 2, guidance 0.05.
      * **The Feynman-Kac estimate fails above kappa 1, so E1 proper is held** (user 10-28 had chosen "read E1
        tuning, run E1"; held on this readout). The potential is zero on open states, so the tilt's normalizing
        constant Z(kappa) = E[exp(-U)] is at least the open fraction (ln 3.8e-4 = -7.9 for I61), and ln Z(kappa)
        must be convex in kappa. The runs' `log_z` is -6.3, -7.3, -33.3 and -84.1 at kappa 0.5, 1, 2 and 4 for I61
        (V70 -9.2, -18.9, -39.7, -88.6): concave, and at kappa 2-4 22-76 nats under that floor. The final step adds
        nothing (mean log-weight about 0), and the collapse follows kappa, not the guidance weight (at kappa 2,
        -33.3 and -30.5). The rule's N_eff over the open particles (187-243) does not see it. Cause not identified:
        weight degeneracy at 256 particles, or a defect in how the resampling increments enter `log_z`.
        `sampler.py`'s bookkeeping reads correctly (the increments telescope along each lineage; the guidance
        ratio is the Gaussian one). E1 proper at kappa 2 would fail calibration on I61 for the estimator's sake,
        and the route would be stopped for the wrong reason.
      * The diagnosis job (user 10-28: "go ahead on boltz job"), its reading fixed here before the data.
        `sampler.py` records a `trace` at every resampling: the step, t_hat, the particles' mean and least energy,
        the mean potential increment, the mean guidance log-ratio since the last resampling, the log mean weight
        (`log_z`'s increment, except at the last step) and the effective sample size. Nothing it computes changes.
        Runs, all on I61 at guidance 0.05: kappa 1 and 2 at 256 particles with seed 1 (the tuning's two runs again),
        kappa 2 at 256 with seed 2, and kappa 1 and 2 at 1,024 with seed 1, kappa 2 at 1,024 with seed 2.
        1. The two repeated runs reproduce the tuning's `log_z` (-7.93, -30.49) and open counts (177, 241 of 256);
           otherwise the instrumented sampler differs, and nothing else is read.
        2. Degeneracy if both 1,024-particle kappa-2 runs have `log_z` above both 256-particle kappa-2 runs, closer
           to the floor -7.9; otherwise the shortfall does not shrink with the particle count, which points to a
           defect in the weights.
        3. The trace locates the steps where `log_z` loses its nats, and whether the effective sample size
           collapses there.
        Submitted 23:41 10-28 as **61755382** (`boltz_addon_e1_diag_ubqps_v100`, V100, 3 h; `remote_jobs.md` §1 line
        4), task file `$B/addon/e1_diag_ubqps.txt`. Code in `$B/addon/`: `sampler.py` `83cdb782` and `e1_rare.py`
        `7f0e1ee5` (the tuning's as `.bak_pre_trace_20261028`), ahead of the MacBook Pro's `addon/`. On stub Boltz
        modules (a toy denoiser and potential, CPU) the instrumented sampler's coordinates, weights, energies and
        `log_z` equal the old one's bitwise, and the trace's increments sum to `log_z`.
      * **Diagnosis readout, 08:20 10-29: rule 1 fails, so nothing else in the job is read.** 61755382 ran 23:43-01:04
        10-29 on midway3-0277, all six runs written, no FAILED line. The repeated runs give `log_z` -6.33 (kappa 1)
        and -28.35 (kappa 2) against the tuning's -7.93 and -30.49, with the same seed, particles and settings in
        both npz, and their samples differ wholesale (backbone up to 48 A apart). Between the two jobs only the
        trace bookkeeping changed (no added random draw) and the node (tuning 0280, diagnosis 0277) under bf16
        autocast. The pilot's parity check showed one process reproducing itself at a seed, but no run has been
        repeated across jobs before, so whether Boltz-2 reproduces a seed across jobs on these GPUs is not known.
        Cause not identified; the log_z failure above kappa 1 stays open and E1 proper stays held. Readout
        `scratchpad/boltz_hdx/addon/runs/boltz_addon_e1_diag_ubqps_v100_61755382/readout.txt` (npz and log copied
        there).
      * The reproducibility job (user 10-29: "fix them all"), its reading fixed here before the data. The tuning's
        I61 kappa 1 guidance 0.05 run (256 particles, seed 1) six times on one GPU in one job, each in its own
        process (`addon/e1_repro.sbatch`): twice with the current code (`new`), twice with the tuning's
        (`old`, the `.bak_pre_trace_20261028` files in `$B/addon/repro_old/`), and twice with the current code under
        `torch.use_deterministic_algorithms(True, warn_only=True)` and `CUBLAS_WORKSPACE_CONFIG=:4096:8` (`det`),
        whose warnings name the nondeterministic operations reached. `addon/repro_compare.py` reads it. Bitwise
        equality of backbones, log-weights and energies decides each pair:
        1. `new_a`/`new_b` or `old_a`/`old_b` unequal: one GPU does not reproduce a seed across processes, so the
           diagnosis's rule 1 cannot be met by equality whatever the code, and the trace's first differing row says
           where the runs part. Then `det_a`/`det_b` equal puts the cause in operations PyTorch's deterministic mode
           controls (named by the warnings), and a rerun of the diagnosis can be deterministic; unequal puts it
           outside them (Boltz's CUDA kernels among the candidates), and any repeat test must compare `log_z` over
           seeds rather than demand equality. The choice of test is the user's.
        2. Both pairs equal and `old_a`/`new_a` unequal: the trace change altered the GPU computation; locate it.
        3. Both pairs and `old_a`/`new_a` equal, `old_a`/tuning unequal: runs reproduce on one node and not across
           nodes (the tuning ran on 0280).
        4. Everything equal: this node reproduces the tuning, and the diagnosis's difference is node 0277's.
        The six `log_z` beside the tuning's -7.93 and the diagnosis's -6.33 are a first look at the run-to-run
        spread at kappa 1, read only as that.
        Submitted 09:24 10-29 as **61786565** (`boltz_addon_e1_repro_ubqps_v100`, V100, 1:30; `remote_jobs.md` §1
        line 4); it runs on midway3-0277, the diagnosis's node, so `new_a` against the diagnosis's run is also a
        repeat across jobs on one node. Old code md5 `20bdd6b7` (`e1_rare.py`) and `9ae0eb7e` (`sampler.py`).
      * **Reproducibility readout, 10:15 10-29: rule 1, a seed is not reproduced across processes on one GPU.**
        61786565 ran 09:29-10:08 on midway3-0277, all six runs written. Every pair differs wholesale (backbone up to
        45 A apart), the current code's and the tuning's alike, so the diagnosis's failed rule 1 says nothing about
        the trace change. The runs part at the first trace row, step 0, before any resampling (mean energy 13.569
        against 13.557), and resampling on those weights picks different particles. Deterministic mode does not
        close it: no deterministic-mode warning in the log, the step-0 difference shrinks from about 1e-3 to about
        1e-6 of the mean energy and remains, so its source is outside what PyTorch's deterministic mode governs
        (cuequivariance's kernels the candidates; not identified). A repeat test has to compare `log_z` over repeats
        or seeds, not demand equality; which test, and whether the diagnosis's rules 2-3 are read under it, is the
        user's. First look at the spread: the eight `log_z` at this setting (six here, the tuning's, the diagnosis's)
        have mean -8.26, s.d. 1.58, range -10.32 to -6.24, and the diagnosis's kappa 1 run sits 1.6 from the
        tuning's, inside it. Readout `scratchpad/boltz_hdx/addon/runs/boltz_addon_e1_repro_ubqps_v100_61786565/readout.txt`
        (npz and log copied there; readout also on `/beagle3`).
      * The diagnosis's rules 2-3 read statistically (user 10-29: "go with your recommendation", the reading rather
        than more seeds first), fixed here before its kappa 2 runs are looked at. The noise allowance is the kappa 1
        spread just measured: s.d. 1.58 per run, so a difference between two two-run means has s.d. 1.58, and 3.2
        nats is twice that. The spread at kappa 2 is not measured; the two seeds at each particle count give a
        first look at it, reported beside and not used to decide.
        2. Degeneracy if the mean `log_z` of the two 1,024-particle kappa 2 runs exceeds that of the diagnosis's two
           256-particle runs by more than 3.2 nats and both 1,024 runs lie above both 256 runs. A difference under
           3.2 nats means the shortfall under the floor -7.9 does not shrink detectably from 256 to 1,024 particles,
           which points to a defect in the weights. Anything else is undecided. The tuning's 256-particle kappa 2 run
           is reported beside, not used.
        3. The trace, on each of the four kappa 2 runs: the resampling rows carrying the most `log_z` loss and the
           effective sample size there. Since the runs share no particles, a step is located only if the same rows
           carry the loss in all four.
        `addon/diag_rules23.py` reads it.
      * **Rules 2-3 readout, 10:40 10-29: by rule 2 the shortfall does not shrink with the particle count, which
        points to a defect in the weights.** Kappa 2 `log_z` at 256 particles -28.35 and -33.63 (seeds 1, 2), at
        1,024 -30.78 and -57.33: the 1,024 mean is 13.1 nats lower, not higher (the tuning's 256 run -30.49 lies
        among the 256 runs). The test's power is low, which the reading has to carry: the two 1,024 seeds differ by
        26.5 nats (the 256 seeds by 5.3), far beyond the kappa 1 spread the 3.2-nat allowance came from, so a
        few-nat gain from more particles would not show. Rule 3: row 0 (step 0, t_hat 4608) carries -25.6 to -26.6
        nats in all four runs, the potential on the first denoised prediction (mean energy 27.1-27.3 in every
        run, the closed state); the rest of the loss falls on different rows in each run, some with the effective
        sample size collapsed (ESS/N 0.001-0.04) and some not (0.45-0.65), so no later step is located. Beside the
        rules: 86-90% of the final particles are open (energy 0), so along any lineage the potential increments sum
        to about zero and the step-0 charge should be repaid in later rows; the per-row means of the potential
        increment instead sum to -195 to -331, meaning the energies of the particles each resampling selects rise
        again by the next one, while the guidance log-ratio sums to only -0.11 to -0.14 per run. The loss sits in
        the potential increments, not in the guidance correction. Whether that is a bookkeeping defect or the
        estimator's own variance at this rarity, these runs cannot separate. Readout
        `scratchpad/boltz_hdx/addon/runs/boltz_addon_e1_diag_ubqps_v100_61755382/rules23.txt` (also on `/beagle3`).
        E1 proper stays held.
      * The toy test (user 10-29: "yes" to it), its reading fixed here before the data. `addon/toy.py` runs
        `sampler.py` unchanged on Boltz's own `AtomDiffusion` at b1ebfc4 (the midway3 tree's Python files, copied to
        `addon/sandbox/`) with E1's settings, on the Mac Studio's CPU in fp32 (Boltz2.forward runs the sampler with
        autocast off). Only the network is replaced, by a deterministic, rotation- and translation-equivariant denoiser
        of a Gaussian ubiquitin: Ub psWT's CA trace with Gaussian-network fluctuations (contacts within 7.3 A, mean
        RMSF 1 A). It aligns the reference onto the noisy coordinates and shrinks the residual in each normal mode by
        lambda / (lambda + sigma^2), so that, like Boltz-2's network, it predicts a folded structure at every noise
        level. The event is the toy's I61: residues 61 and 56 (CA 5.6 A apart) are open when their distance is at
        least 7.5 A, the 3.8e-4 upper quantile in 200,000 calibration samples (I61's unguided open fraction). The
        potential is kappa 7.09 (7.5 - d)+, which is 13.5 on the reference at kappa 1, I61's step-0 energy. The parity
        check passes (8 toy samples, seed 7, coordinates equal). The toy's own distribution is what its plain sampler
        draws, so the truth is measured from 10^6 fresh plain samples: Z(kappa) = E[exp(-kappa 7.09 (7.5 - d)+)] and
        the open fraction, with standard errors. The estimate is ln Z-hat, `log_z` plus the log of the final
        particles' mean weight (in I61 the final step added about 0, so `log_z` alone was read there), and the
        open-fraction estimate is E1's Feynman-Kac one. Runs: kappa 1 and 2, guidance weight 0 (no guidance) and 0.05;
        256 and 1,024 particles with 64 seeds each, 4,096 with 16, 16,384 with 4.
        1. Per setting, the weights are right if the mean of Z-hat / Z over the seeds lies within 3 standard errors of
           1 (the seeds' spread combined with the truth's) at both 1,024 and 4,096 particles, and wrong if it lies more
           than 3 off at both; otherwise undecided. Heavy tails make the seeds' spread understate the error, which is
           why one particle count alone does not decide. The same for the open-fraction estimate, reported beside.
        2. If all four settings pass, `sampler.py`'s weights and `log_z` are right on Boltz's sampler, and I61's
           shortfall is the estimator's variance with Boltz-2's network, or comes from something the toy lacks. If
           only the guided settings fail, the guidance ratio is wrong, and they are rerun with `alignment_reverse_diff`
           off: the ratio leaves out that the noisy coordinates are aligned onto the guided prediction. If an unguided
           setting fails, the Feynman-Kac bookkeeping itself is wrong.
        3. Beside, not deciding: the median of ln Z-hat - ln Z at each particle count, and the trace's step-0 charge
           and the sums of the mean potential increment and of the guidance log-ratio, against I61's.
        A toy run depends on rounding: the same seed at another batch size gives another run (batched and
        one-at-a-time denoiser outputs differ by about 1e-6 relative, the trace differs from step 0, and resampling
        amplifies it), so the toy, like the GPU, is read over seeds.
      * **Toy readout, 12:55 10-29: rule 1 passes in all four settings, so `sampler.py`'s weights and `log_z` are
        right on Boltz's sampler.** Truth from 10^6 plain samples: open fraction 3.70e-4 +- 0.19e-4, ln Z(1) -6.470
        +- 0.017, ln Z(2) -7.314 +- 0.033 (the open states carry 24% and 56% of Z). Mean Z-hat / Z at 1,024 and
        4,096 particles: kappa 1 without guidance 1.031 +- 0.082 and 0.977 +- 0.060, with guidance 1.046 +- 0.038 and
        1.018 +- 0.039; kappa 2 without guidance 0.749 +- 0.132 and 0.927 +- 0.185, with guidance 1.115 +- 0.150 and
        1.011 +- 0.077 (all within 2 standard errors). The open-fraction estimate agrees likewise (0.70-1.22 of the
        truth, within 1.7 standard errors). Beside the rule: the median of ln Z-hat - ln Z rises toward 0 with the
        particles in every setting; it is largest at kappa 2 without guidance, -1.91, -0.93, -0.49 and -0.21 at 256,
        1,024, 4,096 and 16,384 particles, where at 256 the mean Z-hat / Z is 0.505 +- 0.095, the estimator's
        variance at few particles. The toy's step-0 charge is I61's size by the choice of SCALE (-13.5 and -27.0
        nats), and it is repaid through the run, mostly at the last step (its log mean weight +5.7 to +17.2), since
        the toy's prediction keeps the pair near the reference until the noise falls below about 1 A. `log_z` alone
        is -12.9 and -24.8, under the floor by as much as I61's, so in the toy only ln Z-hat, the final step
        included, reads true; in I61 the last step added about 0 and 86-90% of the particles ended open, so its
        `log_z` reading stands. The sum of
        the mean potential increments is -10.6 to -23.3 in the toy, near minus the step-0 energy, against I61's -195
        to -331. In the toy the particles a resampling selects keep their energies, while in I61 their energies rise
        again by the next resampling. The guidance log-ratio sums to -0.6 to -1.1 here (I61 -0.1). So the difference
        is in what the network predicts between resamplings: Boltz-2's denoised structures move the target amide's
        energy from one resampling to the next far more than an equivariant shrinking denoiser does, each
        resampling selects on fluctuations that do not last, and that costs nats every time. This points to the
        estimator's variance with Boltz-2's network, which 1,024 particles did not reduce in I61, or to something
        else the toy lacks; the mechanism is read from the trace means, not measured per particle. Readout `scratchpad/boltz_hdx/addon/runs/toy/readout.txt`
        (`toy_readout.py`; truth and runs beside it, md5 of `toy.py` and `sampler.py` in `md5.txt`). E1 proper stays
        held; what to change (the tilt's schedule over the noise levels, or the kappa used) is the user's.
      * Resampling at low noise only (user 10-29: "run both", on the proposal), its reading fixed here before the
        data. `sampler.py` gains `resampling_t_max`: resampling happens only at steps whose t_hat is at most it, and
        the first resampling, wherever it falls, charges the whole energy, so the weights stay unbiased for any limit;
        guidance is unchanged. `e1_rare.py` and `toy.py` set it with `--resample-below T`; by default there is no
        limit, and three toy grid runs rerun with the changed code equal the old ones bitwise. The limit is
        sigma_data, 16 A: below it Boltz's c_skip exceeds 1/2 and the prediction follows the noisy coordinates more
        than the network's guess. That leaves resampling at t_hat 13.8 down to 2.0 (9 steps) and the last step. The
        old I61 runs at kappa 2 support it: their resamplings above 16 A cost -25.5 to -60.7 nats, those below -8.2
        to +5.7. At t_hat 13.8 those runs still averaged energy 23.5-26.1, so the first resampling will charge most
        particles' energy at once and only the persistence below 16 A can repay it.
        1. On the toy, the same rule as the toy test at the limit (four settings; 256 and 1,024 particles with 64
           seeds, 4,096 with 16): passing in all four shows the changed bookkeeping right; a failure stops the GPU
           job.
        2. On the GPU, I61 at kappa 2, guidance 0.05, 256 particles, seeds 1 and 2 (`addon/e1_lownoise_ubqps.txt`).
           ln Z-hat (`log_z` plus the log of the final particles' mean weight) has the floor -7.9. The collapse is
           removed if both seeds give ln Z-hat of at least -11.1 (the floor less the 3.2-nat allowance from kappa
           1), not removed if both lie below -11.1 (the old schedule gave -28.4 and -33.7 at these seeds), and
           otherwise undecided.
        3. Beside, not deciding: E1's calibration criterion on I61 (the Feynman-Kac estimate within a factor of 2
           of 3.8e-4, N_eff over the open particles at least 50), the open count, the trace's first resampling (its
           charge and effective sample size) and the self-normalized estimate. If the collapse is removed, the next
           steps (V70 at this setting, the tuning rule rerun under the limit, E1 proper) are the user's. If not, the
           limit does not rescue the estimator at kappa 2.
        Toy readout, 13:18 10-29 (`runs/toy_low/readout.txt`, the toy test's truth): rule 1 passes in all four
        settings, so the GPU job goes out. Mean Z-hat / Z at 1,024 and 4,096 particles: kappa 1 without guidance
        0.892 +- 0.063 and 0.937 +- 0.059, with guidance 0.965 +- 0.035 and 1.007 +- 0.030; kappa 2 without guidance
        0.812 +- 0.170 and 1.677 +- 0.488, with guidance 0.994 +- 0.068 and 0.916 +- 0.070 (all within 1.7 standard
        errors); the median of ln Z-hat - ln Z shrinks with the particles as before. The trace now has 10 rows, the
        first at t_hat 13.8, which carries the whole charge (-13.5 and -26.9). In the toy the limit changes little,
        as expected: its loss between resamplings was small to begin with. Submitted 13:19 10-29 as **61805263**
        (`boltz_addon_e1_lownoise_ubqps_v100`, V100, 1 h; `remote_jobs.md` §1 line 4) with `sampler.py` `21df0b41`
        and `e1_rare.py` `e4314698` (the previous ones as `.bak_pre_lownoise_20261029`, md5 `83cdb782` and
        `7f0e1ee5`).
        GPU readout, 14:40 10-29 (`addon/lownoise_read.py` on the run's two npz, run from midway2): **rule 2, the
        collapse is not removed.** ln Z-hat -15.05 (seed 1) and -20.26 (seed 2), both below -11.1, against the old
        schedule's -28.4 and -33.7 at the same seeds: the limit recovers 13 nats and leaves 7-12 short of the floor.
        The first resampling, at t_hat 13.8, charges -25.3 and -22.7 nats (ESS/N 0.085 and 0.005), and the steps
        below it do not repay that. Beside: open 197 and 244 of 256; Feynman-Kac open fraction 2.7e-7 and 1.6e-9
        against the unguided 3.8e-4, so the calibration criterion is not met; self-normalized 0.28 (N_eff 12.6) and
        0.95 (N_eff 255). By the reading fixed above, the limit does not rescue the estimator at kappa 2; what
        follows is the user's. `e1_readout.sh` (flags and C-N) has not been run on these samples.
      * Particle count at low noise (user 11-02: "Particle scaling, full"), its reading fixed here before the data.
        The 256-particle runs' traces show the weights collapsing onto a few particles at most resamplings (ESS/N
        0.00-0.03 at t_hat 11.1-5.6 in seed 1, 0.005-0.01 at 13.8-8.9 in seed 2), with single rows repaying +8.2 to
        +8.4 nats: one lineage's luck, the regime where the toy's deficit shrank with the particle count. The test is
        whether it shrinks here. No code changes: `e1.sbatch` with `e1_rare.py` `e4314698` and `sampler.py`
        `21df0b41`, I61 at kappa 2, guidance 0.05, `--resample-below 16`, at 1,024 and 4,096 particles with seeds 1-4
        and 16,384 with seeds 1-2 (task files `addon/scale/`), beside 61805263's two runs at 256 (-15.05, -20.26).
        `addon/scale_read.py` reads them; ln Z-hat, the floor -7.9 and the line -11.1 as in the rule above.
        1. Removed if at some N of at least 1,024 every seed gives ln Z-hat of at least -11.1 and meets the
           calibration criterion on I61 (the Feynman-Kac open fraction within a factor of 2 of 3.8e-4, N_eff over the
           open particles at least 50); the smallest such N is E1's working particle count. Then the next steps (V70
           at that N, the tuning rule rerun under the limit, E1 proper) are the user's.
        2. Not removed if every 16,384 seed lies below -11.1 and their mean rises less than 3.2 nats over the
           256-particle mean (-17.66): more particles do not rescue kappa 2 at any count E1 can afford, and what to
           change in the sampler is the user's.
        3. Otherwise undecided; the trend of the mean with N is reported.
        Beside, not deciding: open counts, the self-normalized estimate, each run's first-resampling charge and
        ESS/N, its least ESS/N over the rows, and the seeds' spread at each N.
        Submitted 13:09 11-02 from `$B/addon` (`remote_jobs.md` §1 line 4), seven V100 jobs: **62073874** (1,024,
        seeds 1-4 in turn), **62073875-62073878** (4,096, one seed each), **62073879-62073880** (16,384, one seed each);
        about 16 V100-hours, from 256 and 1,024 particles' 375 s and 1,225 s per run. The readout **62073881**
        (`scale/read.sbatch`, caslake, after all seven) writes `addon/runs/scale_readout.txt`. `scale_read.py`
        `5f959e3d` reproduced 61805263's and 61755382's numbers before submission.
        Outcome: the 1,024 and 4,096 arms finished 11-02; the two 16,384 arms ran to their end and lost their results
        at save to `/beagle3`'s inode quota, and the readout job failed with them.
        **Readout, 11-03 on the Mac Studio (the eight finished runs pulled from `$B/addon/runs/`): removal fails at
        1,024 and at 4,096 particles, and the rule's other branch is defined at 16,384, so formally undecided.** The
        256-particle runs reproduce 61805263 (-15.05, -20.26). ln Z-hat at 1,024: -12.29, -19.15, -11.13, -16.79
        (mean -14.84); at 4,096: -17.61, -16.94, -14.54, -16.52 (mean -16.40, s.d. 1.32), every seed below the line
        -11.1 and 7-10 nats under the floor -7.9. The mean rises 1.26 nats from 256 to 4,096 and falls from 1,024 to
        4,096, where the toy's deficit shrank at every step. No seed meets calibration: the Feynman-Kac open fraction is
        2e-8 to 1.3e-5 against the unguided 3.8e-4. The first resampling charges -16 to -24 nats at ESS/N 0.000-0.004
        at every particle count. Kappa 1 estimated I61 right (6.4e-4 and 4.1e-4) but opens no V70 in 256, so the tilt
        this estimator survives does not reach the states E1 needs, and the tilt that reaches them breaks it.
        Readout `scratchpad/boltz_hdx/addon/runs/scale_readout.txt` (also at `$B/addon/runs/`); `scale_read.py`'s
        last-but-one line takes 4,096 as the largest N and is not the rule's verdict. Rerunning the two 16,384 arms
        (about 10 V100-hours, after unpacking Boltz's `.venv` and `cache/mols`) would decide it formally.
        **E1 closed, 11-03** (user: "proceed with your recommendation. but don't give up on the big goal, try another
        potential solution"): the 16,384 arms are not rerun and the Feynman-Kac sampler is dropped. The add-on's
        design (U_phi, the data, E2-E4) carries over to Phase 3's sampler.
- [ ] E2a, Boltz-2's own order on the add-on's domain class (user 11-09: "propose and run next BOltz job"; Decision 10;
      rule fixed 08:37 11-09, before any of its data). The add-on keeps Boltz-2's order inside each element (lambda_i
      shrunk toward zero) because that order is what ranks EHEE; Phase 1 tested it on three proteins (EHEE 0.74, HEEH
      0.04, Ub psWT -0.21). Ferrari et al.'s HDX-NMR set (Zenodo 10.5281/zenodo.14983481, `Dataset_4_HDXNMR.json`,
      CC BY 4.0, md5 `207a51dc`) gives per-position opening free energies for 13 domains of the library's class, 11 of
      them new here: 5 three-helix and 3 EEHEE designs (two of them point mutants of the others) and 3 natural LysM
      domains. If Boltz-2 does not order these, an add-on that keeps its order cannot reach Upside, and one that
      changes it needs site labels that only these 13 domains carry.
      Design. Each domain as its measured construct (No-His where measured, else 6xHis; `addon/nmr_domains.py` writes
      `addon/inputs/nmr_domains.tsv` and `flags/reference/<name>_experiment.csv` from the rows of averaged
      `free_energy_nc_corrected`), MSAs from the ColabFold server by `addon/prep_nmr.py` (the LysM domains 16,384
      sequences; the No-His designs 2-3 rows, their own sequence and at most one other design; the 6xHis designs
      210-225, other tagged sequences matched on the tag). 5 x 10^4 plain `md` samples a domain, 50 runs of 1,000 at
      seeds 13-62, one MSA trunk each, through `e1_rare.py --unguided` after its `--check` (`addon/nmr.sbatch`,
      `t2d.sbatch` with the input named, 25 seeds a job): 22 jobs, about 125 GPU-hours. The bound RT ln N at
      293 K, Ferrari's temperature, is 6.30 kcal/mol, under which lie 44-100% of a domain's measured amides (all of
      them in 8 of the 11). Flags on the Mac by Phase 1's pipeline (`run_flags.sh`'s steps), each domain's reference
      structure the predict step's `model_0`.
      Rules (`addon/nmr_read.py`, written before the data):
      1. Every sample counts; a seed with fewer than 1,000 samples or a job whose `--check` fails is not read.
      2. Per domain, Phase 1's censored Spearman rho_d against the averaged per-position values over all measured
         amides, Boltz-2 and experiment censored at the bound; dG = RT ln(p / (1 - p)) at 293 K. A domain where
         Boltz-2 leaves every measured amide at the bound has no order and counts as 0 (added 09:05, before data).
      3. Across domains, rho-bar is the mean over 9 families (a point mutant and its parent averaged as one family);
         its 95% interval by bootstrap over families, each domain's seeds resampled inside it as Phase 1's blocks.
      Verdict. Order holds: rho-bar at least 0.45, Upside FF2's mean Spearman over Peng's three (0.61, 0.43, 0.32), and
      its interval above 0; E2's library and E3 follow on this order. Order fails: the interval under 0.45; an add-on
      that keeps Boltz-2's order cannot reach Upside, and Phase 2's training premise fails on its own domain class;
      reported, and whether the Boltz line goes on is the user's call. Otherwise unresolved, reported. No variant
      (more seeds, other constructs) is run on any outcome. The bar is Upside's on other proteins, since Upside has not
      run these; Boltz-2's own reference and experiment are the same domains.
      Beside, not deciding: each rho_d with its interval and the Spearman over resolved amides; the two mutant pairs'
      per-site dG differences, Boltz-2's against experiment's; Phase 1's EHEE and HEEH sets against Dataset_4's values
      for them, a check of Peng's data against Ferrari's; the per-condition rows at 293 K; donors open per sample, C-N
      over 2 A and Rg.
      * Checks before GPU time: the code is T2d's (`e1_rare.py` `e4314698`, `sampler.py` `21df0b41`, `opening.py`
        `2518d645`, equal on the cluster); `nmr.sbatch` differs from `t2d.sbatch` only in the input's name and in
        keeping the predict step's PDB; `--test-only` passed on all 22. `nmr_read.py` (`cc8a0e77`) on sets of known
        answer (`addon/runs/nmr_check/synth.py`, 50,000 independent samples a domain, output `nmr_check_out.txt`):
        open fractions at experiment's values read +0.929 [+0.824, +0.949], order holds; permuted -0.017, reversed
        -0.932, both fail; experiment plus noise of SD 1.2, 2.0 and 2.5 kcal/mol +0.643 and +0.499, holds, and
        +0.444 [+0.330, +0.546], unresolved. The input steps run on hhh0012's predicted PDB (60 donors). Its beside
        on Phase 1's existing sets: Boltz-2's EHEE 10^5 against Dataset_4's EHEE 0.773 (against Peng's 0.74), HEEH
        10^4 0.009 (0.04); Dataset_4 against Peng's experiment, Spearman 0.975 over 15 EHEE amides, 0.661 over 19 HEEH. The V100s were all taken at 08:40 with 25
        jobs ahead, so E2a runs on the RTX 6000s, one GPU type for every domain: Phase 0 measured the RTX 6000 at the
        V100's speed with the same bf16 kernels, and each job's `--check` runs on its own GPU. Submitted 08:43 11-09 from `$B/addon` with `sbatch
        --job-name=boltz_plain_nmr_<name>_rtx6000 --partition=gpu --qos=gpu --account=pi-trsosnic --gres=gpu:1
        --constraint=rtx6000 --time=12:00:00 nmr.sbatch NAME FIRST LAST` (`nmr.sbatch` `4a8cbbb5`;
        `addon/submitted.log`): 62687144 and **62687506-62687516**, the 12 blocks of hhh0012, hhh0062, hhh0557,
        eehee0642, eehee0871 and eehee0871m; the other 10 go in as the `gpu` QOS's 12-job limit frees (the watch).
        A first submission at 16 GB (`dc7e2216`, from T2d's 5.7 GB peak) ran out of memory on midway3-0283, whose
        `$SLURM_TMPDIR` is in RAM (findings 10.55); its jobs were cancelled before any seed and resubmitted at 40 GB.
- [ ] E2, data (the data are public on Zenodo, findings 5.7, so no form is needed): a library pilot of about 200 domains (about 85 V100-hours; the
      full 600 about 250), N, CA and C only for `/beagle3`'s inodes; family leakage checked against the library's
      domain list (its Table S1).
- [ ] E3, training: the surrogate (per-flag accuracy on held-out frames reported), then the tilt and the global term
      (CPU or under 10 GPU-hours), selected on the 13 NMR domains.
- [ ] E4, the test, its rule fixed here before the data: steered sampling at experiment's temperature for EHEE, HEEH
      and Ub psWT, two seeds, about 40 V100-hours. Censored Spearman >= Upside FF2's on all three (0.61, 0.43, 0.32)
      with a block bootstrap, and Ub psWT's correlation of error with urea class (now -0.43) within +-0.2. Ablations: no
      tilt, no global term, element terms only, lambda_i unshrunk.

### Phase 3: replica exchange on frozen Boltz-2, its own denoiser as the move (2026-11-03)

User 2026-11-03: "proceed with your recommendation. but don't give up on the big goal, try another potential
solution", after E1 closed (Decision 6).

Why this form. Phase 2 steered the reverse diffusion and carried path weights, and at the tilt that opens Boltz-2's
rare amides those weights collapse (E1). The weights are needed only because steering changes the path. A Metropolis
chain whose move leaves Boltz-2's distribution unchanged needs none: a tilted target p_Boltz(x) exp(-U(x)) is sampled
with acceptance min(1, exp(-(U' - U))), and Boltz-2's density, which no sampler gives, never enters.
* The move, renoise-denoise: from a sample x, draw x_k = x + sigma_k eps at step k of Boltz-2's schedule and run its
  own sampler from step k to the end. If the sampler returned the exact posterior of its own noised distribution, this
  would be the two-block Gibbs sampler on (x, x_k), reversible with respect to p_Boltz, so the Metropolis ratio would
  need only U. Boltz-2 is not exactly that: its steps are discrete, the network is not a conservative score (Plainer
  et al. 2025 found small-noise inconsistency), and below gamma_min = 1 A its steps are deterministic. So the chain's
  stationary distribution at U = 0 can drift from plain Boltz-2's; T1 measures it, and it is the route's main risk.
  sigma_k is drawn at random from a fixed set at every sweep (a mixture of reversible moves is reversible): large ones
  give near-independent proposals, small ones local moves.
* The ladder: rungs with U_r = lambda_r N(x), N the number of protected NHs from the estimator's own flags
  (`protection.py`, equal to the engine's): Peng's denaturant model, so each rung is a urea concentration up to the
  scale s. Neighbouring rungs swap with acceptance min(1, exp((lambda_i - lambda_j)(N_i - N_j))). No gradient is
  needed, so no surrogate.
* The readout: MBAR over every rung with reduced potentials lambda_r N (p_Boltz cancels), so each amide's open fraction
  at lambda = 0 resolves well below RT ln(samples), and its m-value comes off the same ladder. The pooled library with
  its MBAR weights also reweights exactly to p_Boltz exp(-U_phi) for any add-on, the training base Phase 2 needed.
* What it cannot change: the target is still Boltz-2's distribution, and Phase 1 shows it wrong where it is resolved
  (Ub psWT's D21 and strand 1 opened 50-200x too often, no unfolded HEEH, openings of 1-6 NHs). This sampler measures
  Boltz-2's tails and m-values; correcting them stays the add-on's job, now on weights that do not collapse.
* Code: `scratchpad/boltz_hdx/addon/remc.py` (git-ignored, beside `sampler.py`).

- [x] T0, the code on Phase 2's toy (CPU, Mac Studio; `toy.py`'s Gaussian ubiquitin and its 10^6 plain samples as
      truth: P(d >= 7.5 A) 3.70e-4 +- 1.9e-5, P(d >= 7.75 A) 5.2e-5 +- 7.2e-6, ln Z(kappa 1) -6.470 +- 0.017 and
      ln Z(kappa 2) -7.314 +- 0.033 for E1's pair potential). Tilt U_r = -lambda_r d on the pair, 8 rungs from 0 to
      8 /A, 64 independent ladders, moves at sigma 1, 2, 4 or 8 A, 400 sweeps after 50 discarded. Rules, fixed before
      the data, with standard errors over the 64 ladders combined with the truth's:
      1. Parity: the renoise-denoise loop started at step 0 from the plain sampler's initial noise reproduces
         `AtomDiffusion.sample` bitwise (8 samples, seed 7). Failing stops T0.
      2. Invariance: rung 0 alone (lambda 0, every move accepted) gives the truth's mean and s.d. of d and its
         P(d >= 6.5) within 3 standard errors.
      3. Ladder: MBAR's P(d >= 7.5), P(d >= 7.75) and ln Z(kappa 1, 2) at lambda 0 within 3 standard errors.
      Passing 1-3 means the code is right on the toy, and T1 is proposed; failing 2 or 3 is a code defect to find
      before any GPU time. The toy cannot test Boltz-2's own consistency, which is T1's question.
      * Readout, 17:35 11-03 (`addon/t0_read.py`, `runs/t0/readout.txt`): rules 1 and 2 pass, rule 3 fails. Parity
        equal at seeds 7 and 11. Invariance: mean d -2.1 se, P(d >= 6.5) -0.3 se, s.d. +1.2 se, no trend, the same
        after moves at every sigma. Ladder: P(d >= 7.5) 1.84e-4 (-6.2 se), P(d >= 7.75) 1.45e-5 (-4.9 se), ln Z -6.72
        and -7.79 (-4.7, -5.6 se). Diagnosis: the top rungs had not equilibrated. Rungs 0-2 match the truth reweighted
        to their tilt (rung 1 mean d 5.982 against 5.994), while rungs 3-7 still rise between the halves of the run
        (rung 5 7.035, then 7.210, against 7.302) at move acceptance 2-6%. So the run length was short, not the code
        wrong; 400 sweeps had been fixed before the data.
      * Rerun, its rule fixed after that readout and before the rerun's data: the ladder run alone, same settings and
        seed 22, 3,000 sweeps, reading sweeps 1,501-3,000. (a) Its first 450 sweeps equal `ladder.npz` bitwise (the
        `remc.run` refactor that takes the move as a function changes no draw). (b) Stationarity: every rung's mean d
        agrees between the two halves of the read range within 3 standard errors over the ladders; if not, still
        not equilibrated, and the rule is not read. (c) Rule 3 as written. Beside, not deciding: each rung's mean d and
        P(d >= 7) against the truth reweighted to its tilt where that has an effective sample size of at least 1,000
        (rungs 0-4), a test of the tilted chain's stationary distribution with no MBAR in it, so of the move's
        reversibility.
      * **Rerun readout, 18:10 11-03 (`runs/t0/readout_long.txt`): T0 passes, 8 of 8.** (a) The first 450 sweeps
        equal `ladder.npz` bitwise. (b) Every rung stationary (-0.4 to +0.5 se between halves). Rule 3: P(d >= 7.5)
        4.15e-4 (+1.3 se), P(d >= 7.75) 6.8e-5 (+1.7 se), ln Z -6.435 and -7.249 (+0.9, +1.0 se). Beside: rung means
        on the reweighted truth where it is well sampled (rung 1 5.994 and P(d >= 7) 0.0331, against 5.994 and
        0.0333), so no sign of a reversibility bias on the toy at this precision. The top rungs accept 1-4% of moves
        and needed about 1,500 sweeps to equilibrate from plain starts: the cost that sets T2's budget.
- [ ] T1, invariance on Boltz-2 (GPU, about 3.5 V100-hours; user 11-03: "proceed"). One V100 job on Ub psWT,
      `remc_boltz.py` (code beside `remc.py`). The 11-03 cleanup deleted the processed input, so the job first
      reruns `boltz predict` on `inputs/ubqps.yaml` with its stored MSA. Boltz's `.venv` and `cache/mols` stay packed:
      the job extracts both tars to the node's local disk (the venv's Boltz is an editable install of
      `$B/boltz/src`), so `/beagle3`'s file quota (1,163,491 of 1,268,000) is not touched. Then: the parity check
      (`--check`, seed 7); 1,024 plain samples from the same trunk (`--plain`); 256 chains at lambda 0 started from
      that call's own plain samples, 100 sweeps, moves at sigma 1, 2, 4, 8 or 16 A (45-82 of the 200 steps, about 0.31
      of a sample, about 0.38 s per move on a V100 from Phase 1's 20 minutes per 1,000 samples). Rules, fixed here
      before the data, on sweeps 21-100 against the same-trunk reference (the 256 starts and the 1,024):
      1. Parity equal, or nothing else is read.
      2. Every one of the 72 donors' protected fraction within 3.4 standard errors (Bonferroni over 72 at 0.05; the
         chains' standard error over chains combined with the reference's binomial one).
      3. Mean N and mean CA-RMSD within 3 standard errors; the share of states with a C-N bond over 2 A not above the
         reference's by more than 3.
      4. No trend: mean N over sweeps 21-60 and 61-100 within 3 standard errors, paired by chain.
      Passing 2-4 means Boltz-2's own move keeps its distribution, and T2 is proposed. A failure is located by the
      states after moves at each sigma; the candidate remedies (drop those sigma, or run the chain on the 1 A-noised
      state with the deterministic tail as a read-out map) get their own rule before any rerun. Beside, not deciding:
      the same comparisons against the pooled plain 10^5 set (other MSA subsamples), and the time per move. T1 tests
      invariance only (every move is accepted at lambda 0); reversibility under a tilt is T2's calibration.
      * Rule 2's statistics corrected 18:35 11-03, before any T1 data: run on synthetic chains drawn from the pooled
        10^5 flags, which must pass, the readout failed. A donor the reference never opens has a zero binomial error
        at its own proportion, so 0.23% in the chains read z +6.9. The reference's error is now binomial at the
        proportion pooled over chains and reference (the chains counted by their effective number of states), the
        two-sample test under the null; the C-N share likewise. The synthetic pass case then passes (largest |z|
        2.69) and a case with one donor shifted fails (-34.7). `addon/t1_read.py`.
      * 62199777 (submitted 18:13 11-03) failed rule 1 at 19:31: `denoise_from` at step 0 against Boltz-2's sampler
        DIFFERENT, max |d| 36.9 A; the sampler against itself equal. Cause: `remc.schedule` reads the module's
        `num_sampling_steps`, which Boltz2 leaves at AtomDiffusion's default of 5 (`Boltz2DiffusionParams` has no
        such field; the 200 steps are passed per call), so the check ran a 5-step path, and T1's moves would have
        started at steps 3-4 of 5 instead of 118-155 of 200. The toy sets its own count, so T0 could not see it.
        Fixed 11-04 in `remc_boltz.py` (the module's count set to `SAMPLING_STEPS`; `remc.py` `e5e9dfb2` with a
        docstring only; backups `.bak_pre_nsteps_20261104`). The rerun 62262317 passed rule 1 (10:59 11-04: `equal`,
        the sampler against itself equal) and then failed at the first observation: `structure.residues["atom_idx"]`
        is a field of a 43-byte structured array, and torch refuses `start[1:]`'s stride as an index. Fixed by
        copying it to a contiguous int64 array (`remc_boltz.py` md5 `7556f0ae`, backup `.bak_pre_stride_20261104`).
        The code after the check had not run before, so it was run end to end on the Mac first, through the real
        `remc.py`, `protection.py` and `ubqps_protection.npz` with Boltz stubbed by a network that denoises to the
        1UBQ reference (session scratchpad, not kept): the check equal, `t1.npz` written with every array, 72 donors,
        CA-RMSD 0 and the reference's C-N; the harness reproduces both failures on the earlier files. Rerun
        **62263853**, same `t1.sbatch`, rules unchanged (`remote_jobs.md` §1 line 4).
      * **Readout, 18:30 11-04 (`addon/t1_read.py`, `runs/boltz_remc_t1_ubqps_v100_62263853/readout.txt`): T1 fails
        by its rules, 2 of 3; nothing submitted.** 1: equal. 2 fails: largest |z| 4.12, donor rows 53 (+4.1) and 31
        (+3.8) against the limit 3.4. 3 passes: mean N 58.721 against 58.598 (+2.7 se), CA-RMSD 1.594 against 1.560 A
        (+2.3 se), C-N share 0.15% against 0.08% (+0.6 se). 4 passes: 58.705 then 58.737 (+1.0 se).
        Located, not deciding. The two failing donors are open in 1 and 2 of the 1,280 reference samples and in 0 and 2
        of the chains' 20,480 states. A Fisher exact test with every chain state counted as independent, which favours
        finding a difference, gives p 0.059 and 0.019, against 7e-4 for 0.05 over 72 donors: at one or two counts the
        rule's normal approximation overstates the difference, so these two z alone do not show drift. The signs agree,
        though: every listed deviation has the chains more protected (rows 50, 70, 30 at +2.3 to +2.7, mean N +2.7 se),
        and states with 18 or more open NHs are 51 of 20,480 against 8 of 1,280 (p 0.02; at 17 or more, p 0.9). The
        shift is there by sweeps 1-10 (58.73 against the starts' 58.56), flat after, and the same after moves at every
        sigma (58.71-58.73), so it locates to no sigma. The reference's 1,280 samples limit every one of these reads.
        Time per move 0.34 s (87.2 s per sweep of 256 chains), so T2 as estimated is about 18 V100-hours.
- [ ] T1b, T1 again with an 8x reference (user 11-04: "proceed with your recommendation on boltz"). T1's rule 2 failed
      on two donors its 1,280-sample reference opens once or twice, where the reference's binomial error rests on one
      or two counts, and rule 3's mean N read +2.7 se with the reference's error the larger part. One V100 job,
      `t1b.sbatch`: `t1.sbatch` with the run's seed 2 (T1's 1) and `--plain 8192` (T1's 1,024), so the reference is
      8,448 same-trunk samples; the parity check, the 256 chains, 100 sweeps and the sigma set are T1's. For the larger
      reference `remc_boltz.py` observes 1,024 replicas per call (`OBSERVE_BATCH`; the burial coverage would not fit a
      16 GB V100 at 8,448): at T1's sizes that is one call, the code T1 ran; on T1's 256 final frames, chunks of 50
      give flags and C-N bonds equal bitwise and CA-RMSD within 2e-5 A, and the frames' flags equal T1's sweep-100
      record (`addon/runs/t1b_check/observe_test.py`, CPU). About 5.5 V100-hours.
      Rules, fixed here before its data: T1's rules 1-4 as written, read by `t1_read.py` unchanged on T1b's own data.
      T1b is a new test, not pooled with T1, and T1's failure stays on the record. Pass: Boltz-2's move keeps its
      distribution at this precision, and T2 is built with its rule written before its data, checked with
      `--test-only` and submitted. Fail: the move drifts; it is located by sigma as for T1, and a remedy gets its own
      rule before any run; nothing is submitted. Beside, not deciding: T1's rows 53 and 31, the sign of every
      deviation over 2 se, and the share of states with 18 or more open NHs.
      * **Readout, 06:55 11-05 (62293143, COMPLETED 06:03; `addon/t1_read.py`, `runs/boltz_remc_t1b_ubqps_v100_62293143/`
        `readout.txt`, beside reads `beside.txt` from `beside.py`): T1b fails by its rules, 2 of 3; nothing
        submitted.** 1: equal. 2 fails: donor row 0 (residue 2) reads +4.1 against the limit 3.4, protected 0.9989 in
        the chains and 0.9968 in the 8,448-sample reference (open in 23 of 20,480 chain states, 22 of 256 chains, and
        27 of 8,448). These counts are not the one or two that limited T1. 3 passes: mean N 58.688 against 58.632
        (+2.4 se), CA-RMSD 1.582 against 1.567 A (+1.6 se), C-N share 0.13% against 0.13% (+0.0 se). 4 passes: 58.679
        then 58.696 (+0.6 se). 86.8 s per sweep.
        Beside, not deciding. T1's rows 53 and 31 now read +0.7 and -1.1, so T1's two failures were small-count noise,
        as located. Of the 7 rows over 2 se, 6 have the chains more protected: residues 2 (+4.1), 75 (+3.1), 46 (+2.7),
        21 (+2.6), 65 (+2.4) and 66 (+2.4); residue 32 is the exception (-2.0). Over all 72 rows, 35 read positive and
        18 negative. States with 18 or more open NHs: 0.40% against 0.43% (p 0.76); with 17 or more, 1.78% against
        2.13% (p 0.046). Against the pooled plain 10^5, the largest |z| is 3.80. The same sign as T1 (+2.7 se on
        mean N), at about half T1's size: Boltz-2's move leaves the chains slightly more protected.
        Located by sigma. Mean N of the states after a move: +2.7 and +2.6 se at 0.92 and 1.90 A, +0.9, +1.4 and
        +0.6 se at 3.95, 7.67 and 15.95 A. Residue 2 reads +3.3 to +4.0 at 0.92-3.95 A and +1.3 to +1.5 at the two
        largest. The change in N across one move, which is zero at every sigma for a move that keeps the
        distribution, is +0.016, +0.034 and +0.030 at the three small sigma and -0.042 and -0.064 at 7.67 and 15.95 A
        (each within 2 se). The small moves add protection and the large ones, which are close to fresh samples,
        take it back. That fits a bias carried by the moves at 4 A and below, but the chain mixes all five, so the
        states before each move are the chains' own biased ones and this read cannot isolate a sigma. The
        mechanism, the denoiser's late steps sharpening a renoised state more than the sampler's own path does, is
        unproven. For HDX the size matters: residue 2's open fraction 0.11% against 0.32% is a factor of 2.9, wider
        than T2's factor-of-2 target.
- [x] T1c, T1b with only the two large moves (user 11-05: "decide for me the next move for Boltz"). Why this test: in
      a chain sitting above the reference, a move that keeps Boltz-2's distribution can only pull N back toward it.
      In T1b the moves at 7.67 and 15.95 A did (-0.053 and -0.032 per move from the chains' +0.056, about what a fresh
      sample would give), and the three at 3.95 A and below pushed it further up (+0.018, +0.028, +0.028). The
      one-move test from fresh plain samples proposed at 06:55 was dropped before any run: T1b's scatter of N across
      one move (s.d. 1.26 at 0.92 A to 2.28 at 15.95 A) gives a standard error of 0.020-0.036 per sigma at 4,096
      pairs, against effects of 0.02-0.05. One V100 job, `t1c.sbatch`: `t1b.sbatch` with the run seed 3 and
      `--sigmas 8,16` (the schedule's 7.67 and 15.95 A); the parity check, 256 chains, 100 sweeps and 8,192 same-trunk
      plain samples as T1b. About 6 V100-hours (only the two longest moves of T1b's 45-82 steps).
      Rules, fixed here before its data: T1's rules 1-4 as written, read by `t1_read.py` unchanged on T1c's own data,
      not pooled with T1 or T1b. Pass: the large-move set keeps Boltz-2's distribution at this precision, and T2 is
      built on it with its rule written before its data; lambda 0 cannot measure acceptance under the tilt, which
      near-independent moves make lower, so T2's design starts there. Fail: the drift does not come from the small
      moves alone; the second remedy, the chain on the 1 A-noised state with the deterministic tail as a read-out
      map, gets its own rule before any run, and nothing is submitted. Beside, not deciding: the change in N across
      one move at each sigma, residue 2's row, the sign of every row over 2 se, and the time per sweep.
      * **Readout, 00:50 11-06 (62369218, COMPLETED 22:38 11-05; `t1.npz` md5 `1c7a8385` on the Mac and the cluster;
        `runs/boltz_remc_t1c_ubqps_v100_62369218/t1_read.txt`): T1c passes, rules 1-4.** 1: equal. 2: largest |z|
        3.21 against 3.4 (row 21, protected 1.0000 in the chains and 0.9996 in the 8,448-sample reference), then
        +2.4, +2.2, -2.1, +2.0. 3: mean N 58.628 against 58.613 (+0.8 se), CA-RMSD -0.5 se, C-N share -0.1 se.
        4: 58.634 then 58.622 (-0.5 se). Beside from `t1_read.py`: mean N after a move 58.640 at 7.67 A and 58.619 at
        15.95 A; against the pooled plain 10^5 the largest |z| is 7.66, the same-trunk reference itself sitting below
        it (58.613 against 58.648; other MSA subsamples, as in T1b); 109.6 s per sweep of 256 chains, 0.43 s a move.
        So the large-move set keeps Boltz-2's distribution at this precision, and T2 follows, from acceptance under
        the tilt (T2a). The other beside reads are in T2a's output.
- [x] T1d, the chain on the 1 A-noised state (user 11-06: "Proceed with your recommendation", option (b) of T2a's
      readout; rule fixed 01:55 11-06 before any code ran). Why: T2a found the large moves near-independent, so only a
      local move can walk a ladder into the tails, and T1b's local moves (3.95 A and below) carried a drift toward
      protection. Below gamma_min = 1 A Boltz-2's steps are deterministic but for the random rotation, the regime
      where the late steps may sharpen a renoised state more than the sampler's own path does (unproven). Here the
      chain's state is the sampler's own state at the last level above gamma_min (step j, sigma_j), carried with its
      read-out y, the sample the remaining steps make of it. A move adds independent noise of variance
      sigma_k^2 - sigma_j^2, reaching step k's sigma_k, and runs the sampler's stochastic steps from k back to j; the
      read-out is drawn again for the proposal. N and every observable are taken on the read-out. With (x, y) as the
      state the acceptance stays exp(-lambda (N' - N)) though the read-out draws random rotations, and y's marginal is
      Boltz-2's own if the move keeps the sampler's distribution at level j. Code: `remc.py` (`run`'s `level`,
      `denoise_from`'s `end`, `noised_level`), `remc_boltz.py --noised`, `toy_remc.py --noised --sigmas`.
      Checks before any GPU time, on the Mac; failing one is a code defect to find, not a result:
      (i) the refactor changes no draw of the old chain: `toy_remc.py chain` with T0's ladder settings (seed 22, 8
      rungs to lambda 8, 64 ladders, 450 sweeps) equals `runs/t0/ladder.npz` bitwise;
      (ii) the new path on the toy: 512 chains at lambda 0, moves at 2 and 4 A from the noised state, 450 sweeps,
      seed 23, read by T0's rule 2 (mean and s.d. of d and P(d >= 6.5) within 3 standard errors of the truth);
      (iii) `--test-only` of the job.
      The job, `t1d.sbatch`: `t1c.sbatch` with the run seed 4 and `--sigmas 2,4 --noised` (moves at 1.90 and 3.95 A
      from the state at sigma_j); the parity check adds a split check (the sampler run to step j, then on to the end,
      equals the sampler bitwise); 256 chains started from new plain states at step j, 100 sweeps, 8,192 same-trunk
      plain samples. About 5 V100-hours (moves of about 55-65 steps; plain samples at 1.1 s).
      Rules, fixed here before its data: T1's rules 1-4 as written, read by `t1_read.py` unchanged on T1d's own data,
      rule 1 also requiring the split check equal. Pass: the move from the noised state keeps Boltz-2's distribution
      at this precision, and T2a's memory read (`t2a_read.py`) on T1d's data decides as T2a's rule does: an upper
      bound at 0.2 or more at either sigma, and the ladder is built on that move with its own rule first. Fail:
      located by T1b's reads (mean N by sigma, the change in N across a move); the move's remedies are then spent,
      nothing is submitted, and closing Phase 3 at Boltz-2's plain resolution goes to the user. Beside, not deciding:
      T1c's beside list and the memory read.
      * **Checks, 05:40 11-06 (`runs/t1d_toy/`): (i)-(iii) pass; submitted.** (i) The toy ladder rerun equals
        `ladder.npz` bitwise in every array over 450 sweeps (`ladder_regress.log`). (ii) Rule 2 on the noised chain
        (`invariance_noised_read.txt`): mean d 5.650 against 5.653 (-2.4 se), P(d >= 6.5) 0.0593 against 0.0606
        (-2.3 se), s.d. 0.5437 against 0.5457 (-2.2 se), all within 3 but all on the same side, a shift of 0.003 A in
        mean d; first and last quarter equal (5.648, 5.649). (iii) `--test-only` passed. Code `remc.py` `72397622`,
        `remc_boltz.py` `4fc53964`, `t1d.sbatch` `ba884737` (backups `.bak_pre_noised_20261106` on the Mac and the
        cluster; `t0_read.py` gained `--invariance`, its rule-2 output on T0's file unchanged). **62419935**
        `boltz_remc_t1d_ubqps_v100`, submitted 05:42 11-06 (midway3 `gpu`, V100, 8 h).
      * **Readout, 11:35 11-06 (62419935, COMPLETED 11:09, wall 17,123 s; `t1.npz` md5 `d54b7d33` on the Mac and the
        cluster; `runs/boltz_remc_t1d_ubqps_v100_62419935/t1_read.txt`, `t2a_read.txt`): T1d fails rules 2 and 3;
        nothing submitted.** 1: both `[check]` lines equal (step 0, and the split at step 154, sigma 1.016 A). 2 fails:
        row 68 (residue 73) reads +3.7 against 3.4, protected 0.2576 in the chains and 0.2205 in the 8,448-sample
        reference, then residue 72 +3.1 and residues 18 and 69 +2.6. 3 fails: mean N 58.711 against 58.609 (+3.9 se);
        CA-RMSD 1.590 against 1.573 A (+1.5 se), C-N share +0.2 se. 4 passes: 58.698 then 58.725 (+0.7 se). 79.5 s per
        sweep of 256 chains, 0.31 s a move.
        Located by T1b's reads. Mean N after a move is 58.699 at 1.90 A and 58.719 at 3.95 A, both about 0.1 above
        the reference. The change in N across one move is +0.003 (+0.2 se) at 1.90 A and -0.003 (-0.3 se) at 3.95 A.
        The mixed chain's stationarity forces only the average of the two to zero, so each move on its own leaves the
        chains' raised level where it is. A move that kept Boltz-2's distribution would pull chains 0.10 above it back
        by about (1 - r) times that, -0.07 at 1.90 A and -0.09 at 3.95 A at the memory read's r, if N's regression
        across a move is linear (assumed, not measured), as T1b's large moves did (-0.053, -0.032). So both moves from
        the noised state carry the bias, with the sign of T1b's small moves and nearly twice T1b's size (+0.102 against
        +0.056 in mean N). Putting the chain below gamma_min did not remove it, and its mechanism remains unproven. By
        T1d's rule the move's remedies are spent: nothing is submitted, and closing Phase 3 at Boltz-2's plain
        resolution goes to the user.
        Beside, not deciding. Over all 72 rows 40 read positive and 15 negative, and the 4 rows over 2 se all have the
        chains more protected. Residue 2, T1b's failure, reads -0.2 (open in 31 of 20,480 chain states and 12 of
        8,448). Against the pooled plain 10^5 the largest |z| is 3.49. The memory read: r +0.330 [+0.305, +0.353] at
        1.90 A and +0.160 [+0.141, +0.179] at 3.95 A (per-donor medians +0.101 and +0.005), so the 1.90 A move would
        carry a ladder by T2a's rule, but T1d's rule reads memory only on a pass, and a ladder on a move that does not
        keep the distribution would sample the wrong tails. A move costs 0.29 of a plain sample. The toy check (ii)
        had its three reads on one side at about -2.3 se, inside its rule; whether that was this bias at the toy's
        precision is not known.
- [ ] T2, calibration and support (GPU): the N ladder on Ub psWT; I61, D21, E34, Q41 and Y59 within a factor of 2 of
      their plain frequencies, V70, I13, W45 and L67 resolved, as E1 asked. *Estimate*, from the toy's equilibration:
      8 ladders x 12 rungs x 2,000 sweeps, about 20 V100-hours; fixed after T1 measures the time per move.
      * **T2a, acceptance under the tilt, from data in hand (CPU, the Mac; rule fixed 00:55 11-06, before T1c's moves
        were read).** Why first: with only the large moves, a proposal that forgets its starting state is a plain
        Boltz-2 sample whatever the rung. Every rung's chain then filters plain samples, and the ladder samples the
        tilted tails no better than plain samples reweighted by exp(-lambda N) at the same number of proposals (the
        independence sampler's bound); its one gain would be a move's cost against a full sample. How much a large
        move remembers therefore decides T2's design. Read on T1c, sweeps 21-100, where every move is accepted:
        1. The correlation of N before and after one move, per sigma (7.67 and 15.95 A), with a 95% bootstrap interval
           over chains.
        Decision: if the upper bound of 1 is under 0.2 at both sigma, the large-move ladder is not built and nothing
        is submitted. That is reported for the user's decision, with what could sample the tails instead: plain
        samples drawn by large moves at their cost, reweighted; a local move that passes T1's rules; or the chain on
        the 1 A-noised state. If it is 0.2 or more at either sigma, the ladder is built on that move: rungs spaced
        from the pooled plain N distribution reweighted to each lambda for a swap acceptance of at least 0.2, move
        acceptance per rung predicted from the measured memory, its own rule written before its data, `--test-only`,
        then submitted. Beside, not deciding: each donor's agreement before and after a move against independent
        draws' (p^2 + (1 - p)^2); T1c's own beside reads (the change in N across one move, residue 2's row, the
        sign of every row over 2 se); the cost of one move against one plain sample; and the plain open counts of
        the nine target amides in the 8,448 same-trunk and pooled 10^5 samples.
        * **Readout, 01:00 11-06 (`addon/t2a_read.py` `910d5efe`, `runs/boltz_remc_t1c_ubqps_v100_62369218/`
          `t2a_read.txt`): the large moves are near-independent, so the large-move ladder is not built and nothing
          is submitted; the next step is the user's.** 1: r +0.029 [+0.009, +0.048] at 7.67 A (9,216 pairs) and
          -0.011 [-0.029, +0.007] at 15.95 A (11,264), both far under 0.2. Beside: per-donor memory across a move,
          on the 23 donors protected between 1% and 99%, median +0.003 (range -0.012 to +0.065). Change in N across
          one move +0.017 (+1.0 se) at 7.67 A and -0.011 (-0.7 se) at 15.95 A, both zero within error as for a move
          that keeps the distribution. Residue 2, T1b's failure, reads +0.4 (open in 28 of 20,480 chain states and
          13 of 8,448). Six rows over 2 se, four more protected in the chains (residues 24 +3.2, 57 +2.4, 53 +2.2,
          73 +2.0; 8 -2.1, 59 -2.0), and 32 positive against 23 negative over all 72: no sign of T1b's one-sided
          drift. One move costs 0.43 s against at most 1.11 s for a plain sample (0.39). Plain open counts, same-trunk
          8,448 and pooled 10^5: I61 1 and 38, D21 1 and 12, E34 3 and 35, Q41 1 and 12, Y59 0 and 80, L67 1 and 0;
          V70, I13 and W45 0 and 0.
          What it means. This move set samples Boltz-2's distribution correctly, but as fresh draws at 0.39 of a plain
          sample's cost, so a tilted ladder on it gains nothing on the tails over those draws reweighted. Resolving
          an amide open at frequency f takes about 10 / f draws: at 1e-5, 10^6 moves, about 120 V100-hours; V70, I13
          and W45, never open in 10^5, may sit lower. The tails need a move that walks: the moves at 3.95 A and below
          remember the state but carried T1b's drift, and the remaining remedy, the chain on the 1 A-noised state with
          the deterministic tail as a read-out map, is untested. Options for the user: (a) stop Phase 3's ladder and
          use large-move draws as a 2.6x cheaper plain sampler for what plain sampling resolves; (b) test the 1
          A-noised-state chain, with T1's invariance rules on a new rule before any run (about 4-6 V100-hours); (c)
          close Phase 3 at Boltz-2's plain resolution.
      * **T2b, the N ladder on T1d's move (user 11-06: "move to next Boltz job and submit it for me"; Decision 7;
        rule fixed 12:25 11-06, before any ladder data).** Why this move: of the moves tested, only the ones that
        remember their state can walk a ladder into the tails (T2a), and T1d's at 1.90 A remembers enough by T2a's
        rule (r +0.330 [+0.305, +0.353]). Its bias is known (T1d); T2b's calibration measures whether it matters.
        Design. Chains on T1d's state (`--noised`: the sampler's state at step 154) with moves from 1.90, 3.95 and
        7.67 A (`--sigmas 2,4,8`; revised 13:05 from T1d's 1.90 and 3.95 A after check (i), below), 12 rungs at lambda 0-3 per protected NH (spacing 0.273). The rungs come from the
        pooled plain N (mean 58.65, s.d. 1.57; its lower tail falls 1.6-1.9 in ln P per NH from N 57 to 53), so the
        top rung flattens the tail, and a spacing times the s.d. of 0.4-0.8 predicts swap acceptance well above 0.2.
        Two jobs, seeds 5 and 6 (two MSA trunks, as the plain set pools 100), each 4 ladders (48 replicas) for 2,000
        sweeps, the second half read; about 10-20 V100-hours each at T1d's 0.31 s a move, 36 h allowed (the `gpu`
        QOS limit). If the first 10 sweeps' time projects past 34 h, the watch cancels and resubmits with the
        sweeps that fit, the rule unchanged.
        Checks before GPU time, a failure being a code defect to find, not a result: (i) the noised-state ladder on
        the toy (`toy_remc.py chain --noised --sigmas 2,4`, seed 24, T0's long ladder settings: 8 rungs to lambda 8,
        64 ladders, 3,000 sweeps, 1,500 burnt) read by T0's rule 3 with its stationarity check (`t0_read.py
        --stationary`); (ii) `t2_read.py` on synthetic ladders resampled from the pooled plain set at lambda 0-0.6
        recovers the plain open fractions, and its histogram MBAR equals pymbar's weights; (iii) `--test-only`.
        Rules (`addon/t2_read.py`, written before the data): 0. every rung's mean N stationary between the halves of
        the read range within 3 se over the 8 ladders, or nothing is read; 1. both `[check]` lines equal in both
        jobs; 2. calibration: MBAR's lambda 0 open fraction of I61, D21, E34, Q41 and Y59 (the two jobs' mean)
        within a factor of 2 of the pooled plain 10^5's; 3. support: V70, I13, W45 and L67 each open, with no C-N
        bond over 2 A, in at least 20 read states in each job. Pass: T3 reads Ub psWT's dG at the ladder's
        resolution, against experiment and Upside FF2, with T1d's bias stated beside it. Calibration fails: the bias
        matters at HDX's precision, and Phase 3 closes at plain resolution. Support alone fails: the N tilt does not
        reach those amides; reported, and the ladder's lambda range is not widened on that outcome. Beside, not
        deciding: each amide's estimate with a 95% bootstrap interval over ladders, every donor against plain, MBAR's
        mean N at lambda 0 against plain, per rung the move and swap acceptance, mean N, CA-RMSD and the C-N share,
        acceptance per sigma, round trips per ladder.
        * **Checks, 12:50 11-06: (ii) and (iii) pass, (i) fails; nothing submitted.** (ii) On synthetic ladders
          resampled from the pooled plain set (lambda 0-0.6, 2 jobs x 4 ladders x 2,000 sweeps) `t2_read.py` gives
          ratios 0.78-1.20 on the five calibration amides and mean N 58.64 and 58.65 against plain's 58.65, support
          fails as it must (no V70 open in plain), and its histogram MBAR equals pymbar's weights to 1e-16.
          (iii) `--test-only` passed for both seeds (`t2b.sbatch` `96f5f5d6`). (i) The noised-state ladder on the toy
          (`runs/t2_toy/ladder_noised.npz`, 44 min; `ladder_noised_read.txt`) is stationary (every rung within 1.3 se
          between halves) but fails T0's rule 3: P(d >= 7.5) 2.75e-4 against the truth's 3.70e-4 (-3.0 se),
          P(d >= 7.75) 2.85e-5 against 5.2e-5 (-2.7 se), ln Z(kappa 2) -7.531 against -7.314 (-3.2 se). Every tilted
          rung sits below the truth reweighted to it (rung 3, lambda 3.43, ESS 34,935: mean d 6.647 against 6.676),
          where T0's sample-space ladder with moves at 0.92-7.67 A sat at or above it and passed. By the rule above a
          failed check is a defect to find, so T2b waits: two toy ladders at lambda 0 and 3.43 with moves at 1.90 and
          3.95 A, one in sample space and one on the noised state, tell the noised path from the move set.
        * **Diagnosis, 13:02 11-06 (`runs/t2_toy/diag_read.txt`, `t2_toy_diag.py`; 256 ladders x 2 rungs, 1,000
          sweeps, the second half read): the move set, not the noised path.** At lambda 3.43 the sample-space chain
          reads mean d 6.652 against the reweighted truth's 6.676 (-2.3 se) and P(d >= 7) 0.262 against 0.279
          (-1.9 se); the noised chain 6.637 (-4.0 se) and 0.252 (-3.0 se). The two differ by 1 se, and both are
          right at lambda 0 (5.651 and 5.655 against 5.654). So moves at 3.95 A and below under-reach the tilted
          tail in either space, which T0's ladder, with a 7.67 A move beside them, did not (rung 3 at 6.690). No
          code defect. **Revised 13:05, before any ladder data on Boltz-2: T2b's moves are 1.90, 3.95 and 7.67 A**
          (`--sigmas 2,4,8`; the large move is T1c's, which keeps Boltz-2's distribution, and is a fresh draw that the
          tilt mostly rejects). Check (i) is rerun on it (`runs/t2_toy/ladder_noised_248.npz`, seed 27, T0's long
          settings) by T0's rule 3; the rest of T2b's rule is unchanged.
        * **Check (i) rerun, 13:36 11-06 (`runs/t2_toy/ladder_noised_248_read.txt`, 33 min): passes T0's rule 3, 7 of
          7; submitted.** Stationary (every rung within 1.5 se between halves); P(d >= 7.5) 3.13e-4 against 3.70e-4
          (-1.7 se), P(d >= 7.75) 4.60e-5 against 5.2e-5 (-0.6 se), ln Z(kappa 1) and (kappa 2) -2.1 se each. Beside,
          and the limit of this pass: the tilted rungs read as low as without the 7.67 A move (rung 3, lambda 3.43:
          6.648 against the reweighted truth's 6.676; the failed run 6.647), so the better tail estimates may owe
          part of their margin to the seed, and the toy's tails still come out about 0.85 of the truth. T2b's
          calibration on Boltz-2 is the test of whether a bias of that kind matters at HDX's precision. Acceptance
          0.30 at each of the three sigmas. **62453020** (seed 5) and **62453022** (seed 6)
          `boltz_remc_t2b_ubqps_v100`, submitted 13:37 11-06 from `$B/addon` (`t2b.sbatch` `e5a4dc5c`, midway3 `gpu`,
          V100, 36 h; `--test-only` passed first; `addon/submitted.log`).
        * **Read 00:55 11-07: T2b fails, calibration (I61) and support; nothing submitted.** Both jobs COMPLETED
          (seed 5 22:52, seed 6 23:57 11-06; 2,000 sweeps in 32,575 and 32,799 s), pulled to
          `addon/runs/boltz_remc_t2b_ubqps_v100_<id>/` (`t2.npz` md5 `a6f55034` and `9196d762`, equal on the cluster;
          logs, `predict.log`), read by `t2_read.py` into `addon/runs/t2b_read.txt`. Rule 0: stationary, every rung
          within 3 se between the halves of sweeps 1001-2000 (the largest rung 10, -2.9 se; rungs 7-11 all lower in
          the second half, -0.6 to -2.9 se). Rule 1: both `[check]` lines equal in both jobs. Rule 2, at lambda 0
          against the pooled plain 10^5, the jobs' mean: I61 1.77e-4 against 3.80e-4 (38 open), ratio 0.47, FAIL; D21
          1.10, E34 0.69, Q41 0.73, Y59 0.60 pass, so four of the five read below plain. Rule 3: open with no C-N bond
          over 2 A in [seed 5, seed 6] read states, V70 [0, 0], I13 [5, 0], W45 [0, 0], L67 [0, 0], FAIL. By the rule
          the bias matters at HDX's precision and Phase 3 closes at plain resolution; closing goes to the user.
          Beside: MBAR's mean N at lambda 0 58.764 and 58.657 against plain's 58.648; over the 36 donors open at
          least 50 times in plain, ln(ladder / plain) median -0.00, largest +0.56 (residue 51). Per rung, lambda 0 to
          3: move acceptance 1.00, 0.79, 0.60, 0.42, 0.26, 0.16, 0.08, 0.04, 0.02, 0.01, 0.00, 0.00; swap acceptance
          0.75-0.80 inside the ladder (0.39 and 0.37 at its two ends); mean N 58.73 down to 51.12; median CA-RMSD
          1.39-1.60 A; C-N over 2 A at most 0.0016. Acceptance per sigma 0.30, 0.28, 0.27 at 1.90, 3.95 and 7.67 A.
          Round trips per ladder 36, 55, 45, 39 and 19, 76, 62, 43. So the ladder mixes in lambda and its top rung
          holds chains about 7.5 amides less protected than plain, but above lambda 1.9 almost no move is accepted,
          and those chains do not include the four amides plain never opens.
- [x] T2c, one ladder per amide (user 11-07: "Propose next boltz job for me and run it"; Decision 8; rule fixed 16:30
      11-07, before any of its data). Each ladder tilts one amide's opening excess E_a, `opening.py`'s potential at
      kappa 1: burial above 5, plus 4 times the backbone H-bond score above 0.01, plus 4.6 times the ASP/GLU coverage
      above 0.05. E_a is zero exactly where the estimator calls the amide open, and rung r samples p_Boltz
      exp(-lambda_r E_a), so MBAR reweights to lambda 0 as it does for N. The question is E1's, which E1 could not
      answer: how often Boltz-2 itself opens experiment's four least stable amides. Their open fractions at 273 K are
      V70 1.24e-4, W45 8.9e-5, I13 6.9e-5 and L67 6.0e-5 (Peng 2022 SI); the plain tails extrapolate to about 1e-7
      to 1e-8 (E1). The answer is read at a decade's precision. T2b measured the move's bias on five resolved amides
      at 0.47-1.10 of plain, well inside a decade. A pass also leaves a weighted library holding these amides' open
      states, the support Phase 2's add-on would be trained on.
      Design. T2b's chain and moves (the 1 A-noised state; moves at 1.90, 3.95 and 7.67 A). 14 rungs at lambda 0,
      0.25, 0.5, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.5, 1.7 and 2.0 per burial unit, set from the pooled
      plain excesses: their lower tails fall 0.85-1.10 in ln P per unit, so tilts near 1 flatten them, and the
      spacing is 0.1 there, where the tilted spread is widest; plain reweighted to lambda 0-0.75 predicts swap
      acceptance 0.5-0.8 on every target. Five ladders a job, one per amide: I61 for calibration (open in 38 of the
      pooled 10^5; its excess, mean 12.8 and s.d. 2.5, is the nearest in shape to the targets' 13.8-18.4), then V70,
      I13, W45 and L67. Four jobs, seeds 7-10 (four MSA trunks), 70 replicas each, 2,000 sweeps, the second half read:
      about 14 V100-hours each at T2b's 0.34 s per replica and move, 56 in all. Code: `remc_boltz.py --targets
      --lambdas`; `remc.py` unchanged.
      Checks before GPU time, a failure being a code defect to find: (i) `remc_boltz.py`'s excess on the pooled 10^5
      backbones equals the one computed from the engine's own BurialLevel, Hbond and CTerHbond arrays, and is zero
      exactly where the stored flags are open, and each replica reads its own ladder's amide; (ii) `t2c_read.py` (MBAR
      by pymbar 4.0.3) on synthetic ladders (the pooled plain set resampled at its low rungs) recovers plain's I61,
      support fails as it must, and rule 0 passes stationary data; (iii) `--test-only` of the four jobs.
      Rules (`addon/t2c_read.py`, written before the data), per amide over its four ladders:
      0. Stationary: every rung's mean excess agrees between the halves of the read range within 3 standard errors,
         the error the larger of the spread over the ladders and the batch-means error within them (blocks of 100
         sweeps); an amide that is not stationary is not read. (Revised 16:15, before any data, by check (ii): the
         spread over four ladders alone follows Student's t with 3 degrees of freedom, so a limit of 3 flagged half of
         the stationary synthetic amides.)
      1. Parity: both `[check]` lines equal in every job, or nothing is read.
      2. Calibration: I61's lambda 0 open fraction (MBAR on its four ladders) within a factor of 10 of the pooled
         plain 10^5's (3.8e-4); otherwise the tails are not trusted at the verdict's precision and no amide is read.
      3. Support: the amide open, with no C-N bond over 2 A, in at least 20 read states of every job.
      4. Precision: the 95% bootstrap interval of its lambda 0 open fraction over the four ladders spans at most a
         factor of 10.
      Verdict for each target amide passing 0-4, its estimate against experiment's open fraction: more than 10 times
      lower, Boltz-2 is too stable there by more than a decade; within a factor of 10, within a decade; more than 10
      times higher, too open. Every target too stable: Phase 1's missing local openings are measured on Boltz-2's own
      distribution, correcting them is the add-on's job (Phase 2), and the ladders' states with their MBAR weights are
      its support; reported, nothing submitted. Any within a decade: T3 reads that amide's dG at the ladder's
      resolution, T2b's bias beside. An amide failing 3 or 4 is reported unresolved. No variant (rungs, moves,
      sweeps) is run on any outcome. Beside, not deciding: each job's estimate per amide; per rung the move and swap
      acceptance, mean excess and open share; round trips per ladder; the open states' CA-RMSD and C-N; in each
      amide's ladders the other donors at lambda 0 against plain, and MBAR's mean N.
      * **Checks, 16:20 11-07 (`addon/runs/t2c_check/`): (i)-(iii) pass; submitted.** (i) `excess_check.py`, the
        functions compiled from `remc_boltz.py`'s own source: on 2,978 of the pooled backbones (995 of them the
        lowest-excess frames of a target or I61) the coordinate equals the engine's excess within 3.9e-5, is zero on
        exactly the open flags (I61's 38 open frames among them), and every replica of the 70 reads its own ladder's
        amide. (ii) `synth.py`, four synthetic jobs drawn from the pooled plain set at 0.3 times the grid (lambda
        0-0.6), read by `t2c_read.py` (`synth_read.txt`): I61 3.79e-4 against plain's 3.80e-4 (ratio 1.00; jobs
        3.6-4.0e-4), the other donors at plain (mean N 58.65), support 0 open for the four targets as it must be. Rule
        0 as first written flagged four of the five stationary amides (largest -8.5 se) and was revised (above);
        revised, it flags 3.5% and 4.5% of 400 stationary synthetic sets per amide (`drift_calibration.py`; 3.7% for 14
        independent rungs), and V70's +3.1 in the four-job set is that rate. The reader takes 49 s. (iii)
        `--test-only` passed for seeds 7-10. Code `remc_boltz.py` `e912eb54` (backup `.bak_pre_targets_20261107` on
        the Mac and the cluster, `4fc53964`, T2b's), `t2c.sbatch` `647359c3` (`t2b.sbatch` with the T2c arguments),
        `t2c_read.py` `6527a503`. **62551232-62551235** `boltz_remc_t2c_ubqps_v100` (seeds 7-10), submitted 16:26
        11-07 from `$B/addon` (midway3 `gpu`, V100, 35:30 h; `addon/submitted.log`); projected start 05:11 11-08.
        Without the V100 constraint an A100 would start at 02:34, which would mix GPU types against T2b's measured
        bias, so all four stay on V100 as T1-T2b did.
      * **Read 07:00 11-08: T2c fails, calibration (I61 not stationary); no amide is read; nothing submitted.** All
        four COMPLETED 06:05-06:31 11-08 (2,000 sweeps in 45,886-46,992 s), pulled to
        `addon/runs/boltz_remc_t2c_ubqps_v100_<id>/` (`t2c.npz` md5 `1c71727b`, `e07d6f2f`, `49aeeafa`, `7ceb5626`,
        logs too, equal on the cluster; `predict.log`), read by `t2c_read.py` `6527a503` into `addon/runs/t2c_read.txt`
        (exit 1, its failure code). Rule 0: V70, W45 and L67 stationary (largest +2.4, -2.5, -1.9 se); I61 +3.7 se at
        rung 6 and I13 -3.0 se at rung 10, not stationary, not read (the limit flags 3.5-4.5% of stationary synthetic
        sets per amide). Rule 1: both `[check]` lines equal in all four jobs. Rule 2: I61 at lambda 0 1.65e-4 [1.10e-4,
        2.48e-4] against plain's 3.80e-4, ratio 0.44 (jobs 0.97-2.77e-4), inside the factor of 10, but not stationary,
        so calibration fails and no amide is read. No variant is run; what follows is the user's call.
        Beside, deciding nothing. Rule 3, open with no C-N bond over 2 A in the read states of seeds 7-10: V70 [0,
        1335, 1875, 103], I13 [2000, 980, 1686, 778], W45 [0, 920, 234, 1565], L67 [0, 461, 485, 699]; seed 7 never
        opens V70, W45 or L67. Rule 4's intervals span factors of 15 (I13) to 1,800 (V70). The lambda 0 estimates,
        V70 3.1e-10, I13 4.8e-10, W45 3.4e-8 and L67 3.1e-12, lie 3.4 (W45) to 7.3 (L67) decades under experiment's
        6.0e-5 to 1.24e-4, and every interval's upper end at least 2.7 decades under (W45 1.6e-7 against 8.9e-5); on
        tails the calibration did not pass, they are not a measurement. MBAR's mean N at lambda 0 58.67-58.73 against
        plain's 58.648; the other donors at lambda 0 against plain, ln ratio median -0.00 in every amide's ladders
        (range -1.16 to +0.73, I61's down to -5.60). Swap acceptance 0.55-0.93 inside the ladders, 0.25-0.44 at the
        two ends; round trips 13-91 per ladder; move acceptance 0.21-0.23 at each sigma, 0.00-0.06 above lambda 1.2.
        The open states' median CA-RMSD 1.44-2.04 A (V70's 4.68 A), none with a C-N bond over 2 A.
      * Diagnosis, 07:30 11-08 (`addon/runs/t2c_diag/diag.py`, read-only on the four `t2c.npz`): drift between the
        halves per rung, I61 +0.6 to +3.7 se and positive at all 14 rungs (V70 6, I13 2, W45 1, L67 0 of 14); in
        blocks of 200 sweeps each target ladder's top rung falls from an excess of 9-13 to 0-2 once, at sweeps
        617-1766 (first open state), and seed 7's V70, W45 and L67 ladders never open (lowest excess 5.18, 1.36,
        0.74). Not trunk sampling: in the pooled 10^5 the per-trunk open counts of I61, D21, E34, Q41 and Y59 are not
        overdispersed (chi-square 82-105 on 99 degrees of freedom); four random trunks read I61 at 0.47 of the pool or
        less with probability 0.22, but T2b and T2c, eight trunks, both read 0.44-0.47.
- [x] T2d, Boltz-2's exact plain samples bound the four core amides (user 11-08: "Propose next boltz job and submit it
      for me"; Decision 9; rule fixed 07:55 11-08, before any of its data). The pooled plain 10^5 never opens V70, I13,
      W45 or L67, whose experimental open fractions at 273 K are 1.24e-4, 6.9e-5, 8.9e-5 and 6.0e-5: its one-sided 95%
      bound, 3.0e-5, sits only a factor of 2-4 under them. 5 x 10^5 more samples put the bound of an amide that never
      opens at 5.0e-6 pooled (6.0e-6 alone), more than 10 times under every target. Exact sampling has no move to
      bias it, so no calibration amide is needed; if the amides open, the same samples measure them.
      Design. 500 runs of 1,000 plain `md` samples at seeds 113-612, one MSA trunk each (Phase 1's 10^5 used 13-112),
      through `e1_rare.py --unguided` with Phase 1's settings (its 64 unguided samples matched the pooled flags in
      61663101); `addon/t2d.sbatch` runs `e1_rare.py --check` first, then 25 seeds a job, each seed's backbones in
      `s<seed>.npz`, resumable. 20 V100 jobs (seeds in blocks of 25), about 7-8 h each, about 150 V100-hours.
      Flags on the Mac by Phase 1's pipeline (`flags/run_flags.sh`'s recipe), as the pooled 10^5's.
      Rules (`addon/t2d_read.py`, written before the data):
      0. Pooling: on every donor the union opens at least 10 times, the two-sided Fisher exact p of the new set's
         count against the 10^5's above 0.05 divided by the number of those donors; otherwise the new set is read
         alone, the failure reported.
      1. Every sample counts (Phase 1's sets are unfiltered); a seed with fewer than 1,000 samples or a job whose
         `--check` fails is not read.
      2. Per target, k open samples of N: too stable by more than a decade if the one-sided 95% upper bound (exact
         Poisson) is under experiment's open fraction / 10; too open by more than a decade if the one-sided 95% lower
         bound is over 10 times experiment's; within a decade if the two-sided 95% interval (Garwood) lies inside
         [experiment / 10, 10 x experiment]; otherwise unresolved at a decade.
      Verdict. All four too stable: Boltz-2's own distribution, sampled exactly, lacks Ub psWT's four least stable
      local openings by more than a decade; correcting them is Phase 2's add-on's job; reported, nothing submitted.
      Any within a decade: T3 reads that amide's dG at plain resolution. Unresolved or too open: reported. No variant
      (more seeds) is run on any outcome. Beside, not deciding: the calibration amides' frequencies new against old,
      mean protected NHs, Phase 1's Ub psWT readout at the pooled 6 x 10^5 (resolved measured amides, the bound
      RT ln N), C-N over 2 A and Rg, and T2c's beside estimates against the bounds.
      * Checks before GPU time: (i) the code is unchanged (`e1_rare.py` `e4314698`, `sampler.py` `21df0b41`,
        `opening.py` `2518d645`, equal on the cluster); each job's `--check` must pass. Boltz-2 subsamples the MSA
        inside the trunk (`trunkv2.py`, `torch.randperm`) after `seed_everything(seed)`, as in `boltz predict`, so
        each seed is its own trunk; the featurizer's own random state is fixed at 42 in both. (ii) `t2d_read.py`
        on sets of known answer (`addon/runs/t2d_check/synth.py`): 5 x 10^5 Bernoulli draws at the pooled rates
        pool (smallest p 0.40 against the limit 8.8e-4) and read all four too stable; with I61 tripled the set is
        not pooled and read alone; L67 drawn at experiment's rate reads within a decade, V70 at 30 times
        experiment's too open; the pooled 10^5's halves (seeds 13-62 against 63-112) pool (smallest p 0.064 against
        1.2e-3). (iii) `--test-only` passed for the 20 jobs. Submitted 08:05 11-08 from `$B/addon` with `sbatch
        --job-name=boltz_plain_t2d_ubqps_v100 --partition=gpu --qos=gpu --account=pi-trsosnic --gres=gpu:1
        --constraint=v100 --time=12:00:00 t2d.sbatch FIRST LAST` (`t2d.sbatch` `a74f2da6`, `t2d_read.py` `6803192b`;
        `addon/submitted.log`): **62607649-62607660**, seeds 113-412. The `gpu` QOS takes 12 submitted jobs a user
        (`QOSMaxSubmitJobPerUserLimit`), so the last 8 blocks, seeds 413-612 by 25, go in as slots free (the watch).
        A job that stops early is continued with its directory as OUT; finished seeds are kept. Four jobs stopped at
        18:25-18:36 11-08 with their seeds written whole: midway3's epilog deletes every /tmp entry of the user when any
        of the user's jobs on a node ends, so each finished job removed its neighbours' environments (`findings.md`
        10.55). `t2d.sbatch` (`65da611b`) now extracts to the job's own `$SLURM_TMPDIR`; the code each seed runs is
        unchanged. The four went back into their directories, with seeds 413-487, as 62625768-62625774, and the last
        five blocks (seeds 488-612) as 62627802-62627806 at 19:40. 62607660, the last on the old script, stopped the same
        way at 21:09 when 62607650 ended beside it, and went back as 62632705. All 20 directories COMPLETED with `[done]` by
        04:51 11-09, every `--check` equal: 500 seed files, pulled md5-equal to `addon/runs/t2d/`, each 1,000 finite
        unsteered samples at its own seed, pooled in seed order as `flags/samples/ubqps_md_s10_t2d_n500000.npz`; flags
        by `run_flags.sh`'s two steps for that set (500,000 x 72 donors).
      * **Read 05:55 11-09** (`t2d_read.py` `6803192b`, `addon/runs/t2d_read.txt`):
        0. Not pooled: 55 donors open at least 10 times in the union; the smallest Fisher p is D52's, 1.4e-4 (50,231
           of 500,000 against 10,445 of 10^5, open 0.1005 against 0.1045) against the limit 9.1e-4, so the new set
           is read alone. The gap is steady within each set (D52's halves 0.1003 and 0.1006 new, 0.1042 and 0.1047
           old; the new set's halves pool, smallest p 2.1e-3 against 9.1e-4), and D52's per-seed spread is 1.23 times
           the binomial's (0.94-1.32 over the 26 donors open 0.5-99.5% of the time), too little to carry it, so the
           two sets differ by a little at D52; the cause is not measured (`e1_rare.py --unguided` here, `boltz predict`
           and `check_samples.py` for the 10^5).
        1. N 500,000, every sample.
        2. V70: k 1, one-sided 95% upper 9.5e-6 under experiment's tenth 1.24e-5 -> **too stable by more than a
           decade**. I13: k 6, 1.2e-5 [4.4e-6, 2.6e-5] against 6.87e-5 -> unresolved at a decade. W45: k 5, 1.0e-5
           [3.3e-6, 2.3e-5] against 8.90e-5 -> unresolved. L67: k 2, 4.0e-6 [4.8e-7, 1.4e-5] against 6.04e-5 ->
           unresolved.
        Verdict: not all four too stable and none within a decade, so by the rule it is reported and nothing is
        submitted; no variant is run. The three unresolved amides open 6-15 times less often than experiment, their
        95% intervals wholly under it and reaching under its tenth; V70 lacks its opening by more than a decade.
        Beside, not deciding: the calibration amides new against old, I61 2.98e-4 (149) and 3.80e-4 (38), D21 8.4e-5
        and 1.2e-4, E34 3.3e-4 and 3.5e-4, Q41 1.6e-4 and 1.2e-4, Y59 7.1e-4 and 8.0e-4 (Fisher p 0.2-0.7); mean
        protected NHs 58.651 and 58.648; measured amides resolved in the union 22 of 24, the bound RT ln N 7.22
        kcal/mol (10^5: 6.25); C-N over 2 A in 628 of 5 x 10^5 samples (10^5: 131), Rg 11.34 +- 0.16 A (10^5: 11.34).
        T2c's lambda-0 estimates, V70 3.1e-10, I13 4.8e-10, W45 3.4e-8 and L67 3.1e-12, sit 2.0-5.2 decades under the
        lower ends of these exact intervals, as T2c's failed calibration warned.
- [ ] T3, the readout: dG and m-values for all three proteins against experiment and Upside FF2 (Table S14), and
      whether Boltz-2's own distribution holds an unfolded HEEH at all. At plain resolution these are Phase 1's
      readouts (Status): dG Spearman 0.74 on EHEE (FF2 0.61), 0.04 on HEEH (0.32), -0.21 on Ub psWT (0.43); m against
      HEEH's -0.18 (FF2 0.475), unresolved on EHEE and Ub psWT at 10^5; no unfolded HEEH in 10^4. At the ladder's
      resolution, Ub psWT's would have followed a T2b pass; T2b failed (11-07), and T2c failed its calibration (11-08).

## 5. Known risks and blockers

* No energy and no temperature input: Boltz-2's frequencies are not Boltzmann weights, and dG resolution is RT ln(N)
  at one condition (2.1). Phase 2's tilt supplies the correction, learned from experiment at 20 C, so its temperature
  and pD dependence (the tests are at 273-298 K, pD 4.6-7.6) is extrapolated, by 0.5-1 kcal/mol at worst.
* Global and subglobal openings (HEEH's helices and hairpin, ubiquitin's strand 1) are rare events on
  millisecond-to-second timescales. mdCATH's 320 K runs of up to about 500 ns do not contain them, so that data cannot
  teach their populations. Local openings (Ub psWT's V70, I13, W45, L67) are fluctuations of the folded state at
  populations near 1e-4; whether mdCATH-length MD samples them is not measured.
* Boltz-2 trained on mdCATH: Phase 1 domains must come from its held-out clusters. The 89-domain panel holds stable
  native states, a weak test of transitions (Phase 1 adds switching domains for this).
* `step_scale` 1.0 gives bond-length scatter 2 to 4x the crystal spread. It does not change the protection flags
  (Phase 0), but any use beyond the HDX estimator, such as comparisons of local geometry with MD, has to deal with it.
  Boltz-2 also draws salt bridges and H-bonds 0.3 A or more too short at every setting.
* Phase 2's sampling correction `ll_difference` treats the previous guidance shift as noise and ignores the update's
  Jacobian; if that bias matters, every steered estimate is biased (E1 tests it).
* Phase 2's labels are unassigned dG sets: the right values can land on the wrong residues, and only the 13 NMR
  domains check the assignment. The library holds 28-64 residue domains; Ub psWT has 76.
* With 15-24 measured amides per protein, Spearman intervals are about +-0.3, so a tie with Upside may stay undecided.
  Family leakage through the HDX-MS library is checked against its domain list.
* Phase 2's compute is about 300 V100-hours in all; E1 alone, about 7, commits nothing else.

## 6. References

* Boltz: https://github.com/jwohlwend/boltz ; Boltz-1 report https://gcorso.github.io/assets/boltz1.pdf ; Boltz-2
  preprint https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1.full
* Peng et al. 2022, Upside HDX (ubiquitin, EHEE_rd2_0005); findings 5.1
* Arts et al. 2023, Two for One: https://arxiv.org/abs/2302.00600
* Plainer et al. 2025, Consistent Sampling and Simulation: https://arxiv.org/abs/2506.17139 (code: noegroup/ScoreMD)
* Jing et al. 2024, AlphaFlow: https://arxiv.org/abs/2402.04845
* Schreiner et al. 2023, Implicit Transfer Operator Learning (NeurIPS 2023)
* BioEmu2 (MIT): https://github.com/microsoft/bioemu2 at `f99e0d8` (code only) and `06be369` (section 7)
* Engh and Huber 1991, Acta Cryst. A47, 392 (peptide bond geometry)

## 7. Side test: BioEmu2's weights for HDX (separate record)

User 2026-10-24: "Set up another job to test whether we can use bioemu2 ... to get hdx results"; asked about the rule
of section 1, "Lift it, side test only". Nothing here is a baseline or comparator for the Boltz work; it is judged
against experiment and Upside's published values only, through the same estimator as Phase 1.

BioEmu2 at `06be369` (2026-10-23; the install, checkpoint and CLI code are unchanged since `f99e0d8`):
* Python >= 3.12; `.[cuda12]` with the cu126 index (torch 2.13.0+cu126, the CUDA 12 build that keeps sm_70). Weights
  `microsoft/bioemu2` on Hugging Face, `checkpoints/bioemu2-preview/` (2.0 GB, MIT, not gated).
* Input: a system YAML with the sequence, and always an MSA, from ColabFold through BioEmu2's `get_msa`, cached in
  `BIOEMU2_MSA_CACHE` and reused offline. No temperature, denaturant or pH input; the equilibrium mode is the `md`
  method token. Sampler `edm` (unsteered), 200 steps.
* Output: 5 files per run whatever N (`samples.xtc` of independent samples, `topology.cif`, all heavy atoms).
* Its model card says sample frequencies are not calibrated Boltzmann probabilities, and no HDX evaluation exists.
* **Training data overlap.** MD set `MSR_megasim_merge` is 170 of MegaSim's 271 systems, MD of Tsuboyama 2023 wild
  types folded and unfolded. The 271 include Rocklin designs (EHEE_rd1-4, HEEH_rd3_0223, HEEH_rd4_0094, ...) and four
  ubiquitin mutants, but not EHEE_rd2_0005, HEEH_rd4_0097 or wild-type ubiquitin; which 170 were used is not
  published. Close relatives of all three test proteins, unfolded states included, may be in training, so even an
  agreement would not be a held-out test.

Two settings depart from the defaults, for the V100: `trimul_variant: full` in a copy of the checkpoint's
`config.yaml`, the variant BioEmu2's own size tiers give up to 200 tokens (`triangle_multiplication.py`; the default
`cueq` is a Triton kernel untested on sm_70, and Boltz turns its cuEquivariance kernels off below sm_80); and
`--precision 32-true`, which affects only the trunk (denoising runs in fp32 regardless; bf16 is emulated on Volta).

- [x] Local sandbox on the Mac (CPU, 2026-10-24): a 4-sample EHEE run with the edited config (fp32, 3.6 min) writes the
      5 files, all 4 passing BioEmu2's physicality check (CA-RMSD to the design 0.61-0.73 A, C-N 1.337 +- 0.029 A).
      `bioemu2_samples.py` turns a run into `<name>_n<N>.npz` with `backbone` (N, CA, C, A) and the Phase 1 checks;
      its batched Kabsch equals `check_samples.py`'s exactly.
- [x] midway3, `/project/trsosnic/yinhan/bioemu2_hdx` (`$E`; `/beagle3`'s inodes are kept for Boltz): uv's managed
      CPython 3.13.17, because the cluster's `python/3.13.11` module has no `_ctypes` and torch cannot import; torch
      2.13.0+cu126 with sm_70. Weights and MSAs fetched on the login node (Ub psWT 21,123 sequences, EHEE and HEEH 2
      and 4, as Boltz's). `--test-only` waits for a free gpu slot, because the 12-job limit refuses it too;
      `submit.sh` runs it before submitting.
- [x] One V100 job: EHEE_rd2_0005 10^3, Ub psWT 10^3, then HEEH_rd4_0097 10^4 (its measured amides lie below the
      10^4 bound, so it is the direct test), 1,880, 3,204 and 19,614 s (data paths in `remote_jobs.md` "Finished
      runs"). Every sample passes BioEmu2's physicality check. CA-RMSD to the reference median 0.74, 2.22 and 1.75 A;
      Rg 9.25, 11.31 and 9.57 A (references 9.02, 11.49, 9.27); C-N 1.337-1.340 +- 0.03 A. 3 of EHEE's 1,000 and 8 of
      HEEH's 10^4 lie beyond 8 A of the reference, extended chains with Rg up to 20.4 and 25.9 A; none of Ub psWT's.
- [x] Readout 2026-10-25: flags and dG as Phase 1, against experiment and Upside's published values; one
      `plot_ref_style.py` figure per protein. In `scratchpad/bioemu2_hdx/readout/` (`readout.sh`: Phase 1's
      protection-state code, `write_traj` and `readout()` on BioEmu2's own samples, which stay in that directory;
      `dg/dg_bioemu2.log`, `dg/<name>_bioemu2_n<N>_dG_vs_residue.png`). Each set's halves are its independent draw.
      * HEEH_rd4_0097 (278 K, bound 5.09): all 20 measured amides resolve (0.97-5.09 kcal/mol). Against experiment
        Spearman 0.239 (p 0.31; bootstrap over 10 blocks of 1,000, 0.136 to 0.352), Pearson 0.258, RMSE 1.04, mean
        offset +0.20; Upside FF2's published 0.317, 0.185, 1.06, inside that interval. Neither Spearman differs from
        zero at 20 amides. The halves agree (Pearson 0.98 at 18 amides).
      * HEEH's unfolded state is sampled, rarely: 4 of 10^4 samples have more than half their donors open, and the
        helix amides open at mostly 3-8e-4, that population, against experiment's 4e-4 to 3e-3. The errors are
        elsewhere: the hairpin's V19 and W21 open at 7e-4 and 6e-4 where experiment has 0.020 and 0.014, T24 at 0.15
        against 0.010, E4 at 0.028 against 5e-4.
      * EHEE (298 K, bound 4.09) and Ub psWT (273 K, bound 3.75) at 10^3 test only over-opening: BioEmu2 falls below
        min(experiment, bound) at 9 of 15 EHEE amides (Upside FF2 at 4) and at none of 24 Ub psWT amides (Upside 7).
      * Ub psWT's halves differ at T9 (protected 0.36 against 0.25, Fisher p 1.5e-4, below 0.05/72; blocks of 100 read
        0.21-0.43): the sample order is not fully exchangeable there. No other donor differs.
- [x] EHEE and Ub psWT at 10^5 (user 2026-10-25, the controlled tests), the resolving test, at BioEmu2's default
      settings and MSAs: 10 jobs of one run of 10^4 per protein, seeds 1000, 1100, ..., 1900. BioEmu2 seeds batch i
      of a run with SEED + i, so runs of 100 batches need seeds 100 apart; these are independent of the first job's
      seed 1. `$E/sample.sbatch` takes `NAME:N:SEED`, submitted through `$E/submitted_1e5.txt` and `$E/submit_1e5.sh`;
      about 52 and 89 V100-hours. Readout fixed before the data: the Boltz-2 10^5 statistics (resolved Pearson and
      Spearman, censored Spearman with a block bootstrap, Upside FF2 on the same comparison, halves, open fractions,
      shortfall below min(experiment, bound)); with HEEH this decides whether BioEmu2 ranks these proteins' HDX.
      Sets `samples/{ehee,ubqps}_bioemu2_n100000.npz`, pooled in ledger order (halves: seeds 1000-1400 against
      1500-1900; `readout.sh` reproduces every earlier set byte for byte); figures
      `dg/{ehee,ubqps}_bioemu2_n100000_dG_vs_residue.png`.
      * EHEE (2026-10-27): CA-RMSD to the design median 0.74 A, Rg 9.25 +- 0.52 A, C-N 1.337 +- 0.028 A with none over
        2 A, no sample repeats; 265 of 10^5 lie beyond 8 A (Rg up to 34.9 A).
        - BioEmu2 opens EHEE far too often. All 15 measured amides resolve, but at 3.52-4.70 kcal/mol against
          experiment's 4.02-8.41: every one falls below min(experiment, bound) (mean shortfall 2.43 kcal/mol; Upside
          FF2 at 11 of 15, 1.01), with open fractions of 4e-4 to 2.6e-3 against experiment's 7e-7 to 1.1e-3.
        - The cause is a globally unfolded population. 184 of 10^5 samples (0.18%) have more than half their donors
          open, and they give 77% of the measured amides' openings (61-94% per amide; 97% come from samples beyond 4 A
          of the design). An unfolded population near 0.2% puts every amide near RT ln(1/0.002) = 3.7 kcal/mol, where
          experiment's most stable amides (8.4 kcal/mol) put it near 1e-6.
        - The order is partly kept: Spearman 0.561 (p 0.03) and Pearson 0.648 over the 15; censored at the bound 0.457
          (bootstrap 95% 0.356 to 0.603), above Upside FF2's 0.609 in 1.8% of the resamples. The halves agree to
          0.08 kcal/mol rms (smallest Fisher p 0.162), so this is the model's ensemble, not counting noise.
      * Ub psWT (2026-10-27): seed 1700's `physicality.json` rejects one frame ("isolated_atoms"), which BioEmu2 still
        writes to `samples.xtc`; its N, CA and C, all the flags read, are ordinary (CA-RMSD 2.13 A, C-N at most
        1.40 A), and it is kept, as every sample is. CA-RMSD to the reference median 2.22 A, at most 4.17 A; Rg 11.31
        +- 0.16 A; C-N 1.338 +- 0.029 A, one sample over 2 A (bond 74-75, the C-terminal tail); no sample repeats.
        - BioEmu2 holds ubiquitin in its native basin. No sample lies beyond 4.2 A of the reference or has more than
          29% of its donors open, so the measured amides open only by local fluctuations of the folded state: 7 of the
          24 ever open (Y59 41 times, L69 14, L50 2, K27, A28, I30 and L56 once each), and the samples that open Y59 or
          L69 lie no farther from the reference than the rest (median 2.2-2.4 A).
        - Experiment's least stable amides open through fluctuations BioEmu2 does not sample. V70, W45, I13 and L67
          (open fractions 6e-5 to 1.2e-4) would open 6-12 times in 10^5 samples and never do (Poisson p 4e-6 to
          2.4e-3), while Y59 opens ten times as often as experiment's 4e-5.
        - The order is not kept: Spearman 0.788 (p 0.04) over the 7 resolved (Upside FF2 0.679 at the same 7), but
          censored at the bound over all 24, 0.063 (p 0.77; bootstrap 95% -0.068 to 0.425), above Upside FF2's 0.432
          in none of 2,000 resamples. The shortfall test sees only over-opening: BioEmu2 falls below min(experiment,
          bound) at 6 of 24 (mean 0.07 kcal/mol; Upside FF2 at 21, 1.65), 4 of them a single opening at RT ln(99999),
          1e-5 kcal/mol under the bound; the real ones are Y59 (4.23 against 5.45) and L69 (4.81 against 5.16). No
          donor differs between the halves (smallest Fisher p 0.011).
      * The test's answer, with HEEH: at its default settings and MSAs BioEmu2 does not rank these proteins' HDX as well
        as Upside FF2. It loses the order on Ub psWT (censored Spearman 0.06 against 0.43), keeps part of it on EHEE
        (0.46 against 0.61, above Upside in 1.8% of resamples), and on HEEH neither model's Spearman differs from zero.
        Its stability errors run both ways: EHEE unfolds globally in 0.18% of samples, about 2.4 kcal/mol too often,
        and ubiquitin never leaves its native basin. The training overlap above still holds.
- [x] m-values (2026-11-03; Phase 1's m-value item, same code and rules; EHEE and Ub psWT on the 10^5 sets, at
      `$E/readout/` since the MacBook Pro's upload of 16:25 11-03, pulled md5-equal and reproducing the 10-27 dG
      readout). BioEmu2 gets the size of the openings right on HEEH only.
      * HEEH_rd4_0097 (10^4): m Spearman 0.465 over all 20 (bootstrap over the 10 blocks 0.14 to 0.57; Pearson 0.521),
        FF2 0.475 (Table S14), where its dG gave 0.239. Its helix amides open in its few unfolded samples, breaking
        11-22 protected NHs (3-8 samples each), and experiment's slopes put them at 0.76-1.09 of m_global; Boltz-2's
        break 0.9-3.2.
      * EHEE_rd2_0005 (10^5): every measured amide opens through openings of 14-23 of 32 protected NHs, its 0.27%
        unfolded population, where experiment's slopes put 11 of the 15 at 0.22-0.67 of m_global. m Spearman 0.179
        (bootstrap 0.02 to 0.31), below FF2's 0.400 on the same 15. Boltz-2 errs the other way, at 0.9-5.4 NHs.
      * Ub psWT (10^5): 6 of 23 open, all at 1.1-2.2 NHs (K27, which experiment exchanges globally, in one sample);
        too few to score (Spearman 0.27, bootstrap -0.32 to 1.00). It samples no unfolded state, so no s.
- [x] The s-dependent half of Table S14 (2026-11-03, `flags/urea/high_urea.py`, log `urea/high_urea.log`). s by
      Peng's rule, the model's global stability given experiment's slope: s = m_global / (<N>_folded - <N>_unfolded),
      recalibrated in every bootstrap resample. Unfolded is CA-RMSD beyond 8 A, with more than half the donors open as
      the check. The model's global line has experiment's slope at 0 M (0.8712 and 0.6070), and finite-difference
      slopes equal the analytic ones to 4e-4 kcal/mol/M.
      * HEEH, 0 M: the whole row matches Upside FF2. dG Pearson 0.258, Spearman 0.239, RMSE 1.04 (FF2 0.185, 0.317,
        1.06); m 0.521, 0.465, RMSE 0.40 or 0.34 by the two definitions (FF2 0.541, 0.475, 0.404). Global stability
        3.94 or 4.32 kcal/mol against experiment's 4.32. Unfolded is 8 samples beyond a gap from 4.53 to 11.64 A, or 4.
      * HEEH, 2.13 M, not determined: m Spearman 0.37 or 0.36 (bootstrap 0.16 to 0.48; FF2 0.439), dG Spearman 0.03 or
        -0.01 (FF2 0.036), m RMSE 0.68 or 0.27 (bootstrap 0.37 to 1.39; FF2 0.216), global stability 1.25 or 2.36
        (experiment 2.46). These rest on the few unfolded samples. Under the CA-RMSD definition, the one sample with 8
        protected NHs holds 6% of all weight at 2.13 M and sets 18 of the 20 amides at 1.2-1.5 kcal/mol.
      * EHEE is not comparable with Peng's FF2 at high urea. Peng reads FF2 at the temperature where its global
        stability matches experiment's, and BioEmu2, with no temperature, is 3.51 or 3.73 kcal/mol against 8.57 (265
        or 184 unfolded samples; no gap in CA-RMSD, 7.96 to 8.20 A). At 4 M every amide then reads 0.4-1.0 kcal/mol
        against experiment's 3.3-6.1 (dG RMSE 4.8 or 4.2; FF2 1.14), its m all 0.84-0.93 against 0.19-0.61 (RMSE 0.50
        or 0.36; FF2 0.175). The dG order holds (Spearman 0.72 or 0.69, bootstrap 0.49 to 0.76; FF2 0.689), and the m
        order follows the definition (-0.22, bootstrap -0.50 to -0.02, or 0.57; FF2 0.600).
