# A modified Boltz-2 to replace Upside's sampling for HDX

**Status (2026-10-27).** Phase 0 is done (10-21). Phase 1's HDX half is read out: Boltz-2 ranks EHEE's amides at
10^5 samples (censored Spearman 0.74; Upside FF2 0.61) but not HEEH_rd4_0097's at 10^4 (0.04; Upside 0.32) or Ub
psWT's at 10^5 (-0.21, and -0.06 with a query-only MSA; Upside 0.43). By experiment's urea dependence EHEE's measured
amides exchange through subglobal openings, HEEH's helices by global unfolding and Ub psWT's amides below the bound
locally, so Gate 1 ends in its "wrong states" outcome; its mdCATH half is not run. Phase 2 (user 10-27: no retraining
of Boltz, small trained add-ons allowed) is a trained add-on that tilts frozen Boltz-2's ensemble; no fine-tune is
planned. Phase 2's E1 GPU pilot 61663101 (sampler parity and 64 unguided samples) is queued on midway3; the next step
after it is the user's (E1). Jobs are recorded in `remote_jobs.md` §1. Every number below comes from reading code,
docs and papers on 2026-10-21, or from the dated measurements; numbers marked *estimate* are unmeasured.

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

- [ ] E1, the falsifier (about 7 V100-hours, no training): untrained Boltz-2's own rare open fractions by particle
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
        end to end. The first opening potential (`opening.py`, Boltz's own N and O atoms) left Upside's H-bond at V70 in
        place (0.43-0.98): the estimator infers H and O from N, CA and C, so the potential has to act through that
        geometry and through the coverage burial, which the surrogate of E3 needs anyway.
      * What the support check is up against, from the 10^5 unguided flags: an amide is open only with H-bond <= 0.01
        and burial <= 5, and V70, I13 and L67 never come below burial 10.2 (0.001% quantile; 0-1 of 10^5 below 9), W45
        below 6.4. Their tails fall a decade per 1.5-2.5 burial units, so Boltz-2's own open fraction there extrapolates
        to about 1e-7 to 1e-8 before the H-bond condition, against experiment's 6e-5 to 1.2e-4. Amides that do open have
        heavy tails (Y59 88 samples below 5, D21 12). Extrapolated, not measured: E1 measures it.
      * Status: the GPU pilot 61663101 (the parity check on the GPU and 64 unguided samples) is queued on midway3
        (`remote_jobs.md` §1 line 4). Next is the user's: an opening potential through the estimator's own geometry,
        then E1's calibration and support runs (recommended), or the route closed.
- [ ] E2, data (the user requests it through the form): a library pilot of about 200 domains (about 85 V100-hours; the
      full 600 about 250), N, CA and C only for `/beagle3`'s inodes; family leakage checked against the library's
      domain list (its Table S1).
- [ ] E3, training: the surrogate (per-flag accuracy on held-out frames reported), then the tilt and the global term
      (CPU or under 10 GPU-hours), selected on the 13 NMR domains.
- [ ] E4, the test, its rule fixed here before the data: steered sampling at experiment's temperature for EHEE, HEEH
      and Ub psWT, two seeds, about 40 V100-hours. Censored Spearman >= Upside FF2's on all three (0.61, 0.43, 0.32)
      with a block bootstrap, and Ub psWT's correlation of error with urea class (now -0.43) within +-0.2. Ablations: no
      tilt, no global term, element terms only, lambda_i unshrunk.

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
