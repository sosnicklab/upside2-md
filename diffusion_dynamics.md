# A modified Boltz-2 to replace Upside's sampling for HDX

**Status: Phase 0 done (2026-10-21); Phase 1's first job is queued.** `step_scale` 1.0, the unscaled sampler, gives
bond-length scatter 2 to 4x the crystal spread, but the scatter leaves the HDX protection flags unchanged, so Phase 1
uses `md` at 1.0 (Phase 0). midway2 has no GPU that can run Boltz-2 (Phase 0). Jobs are recorded in `remote_jobs.md` §1.
Every number below comes from reading code, docs and papers on 2026-10-21, or from the dated measurements in
Phase 0; numbers marked *estimate* are unmeasured.

## 1. Goal

Build a modified Boltz-2 that samples protein conformations for HDX in place of Upside (user, 2026-10-21). Judge it
by the per-residue HDX dG it gives through the repo's estimator, against experiment and against Upside, and judge its
dynamics against all-atom MD on mdCATH. Find out whether a route with no training is good enough before spending GPU
time on a fine-tune.

**Built on Boltz, not BioEmu** (user, 2026-10-21): Jumper does not accept BioEmu because its data are selected. No
BioEmu weights, samples, training data or benchmarks are used (findings 10.34). BioEmu2's code may be used (user,
2026-10-21): it is MIT-licensed, and its layer, denoising and steering code is itself adapted from Boltz (its
`NOTICES.md`), so it ports onto Boltz's modules directly. Phase 3 takes its lag conditioning from there.

## 2. Background

### 2.1 What the HDX estimator needs from a sampler

From findings 5.1 and 5.2:
* The estimator is equilibrium. Each frame gets a binary protected flag per amide (backbone H-bond score, an Asp/Glu
  side-chain contact proxy, burial); MBAR gives the protected fraction p, and dG = RT ln(p/(1-p)) under EX2. Frame
  order is not used. Dynamics matter for decorrelation, and for whether EX2 holds.
* The flags need only N, CA and C per residue. The hybrid already feeds the estimator that way (a projection onto
  `3*n_res` coordinates), so Boltz frames can go through the same code.
* With binary flags the largest resolvable dG is about RT ln(N_eff): 4.1 kcal/mol at 298 K from 1,000 independent
  samples, 6.8 from 10^5. Protected amides in stable proteins sit at 5 to 12 kcal/mol, open fractions of 10^-4 to
  10^-9.
* Upside reaches its present limit (5.5 kcal/mol on glpG) with REMD and MBAR, which need an energy and a temperature.
  Boltz-2 has neither: `temperature_bins` is defined in `src/boltz/data/const.py:468`, but no model code reads it.
  Boltz sample frequencies are therefore the only population estimate, and nothing in its training calibrates them
  as Boltzmann weights.
* Boltz sees no membrane, so glpG and the TM4 question stay with Upside. Tests use soluble proteins.

### 2.2 Boltz-2 at `b1ebfc4`

Source `github.com/jwohlwend/boltz` (version 2.2.1, Python 3.10 to 3.12). Paths are relative to `src/boltz/`.
* Prediction runs in `bf16-mixed` (`main.py:1262`). cuEquivariance kernels are used unless `--no_kernels`.
* `--method` picks a tag from `data/const.py:440` (`md` is 0). The tag is embedded and added to the trunk's single
  representation (`model/modules/trunkv2.py:199`), so it reaches the diffusion module only through the trunk.
* MD training data: poses from MISATO, ATLAS and mdCATH, used in only 8k of 97k training steps (stages 2 and 3).
  The preprint evaluated MD conditioning on held-out mdCATH and ATLAS clusters and reports MD-conditioned ensembles
  less diverse than MD.
* The trunk (MSA module and PairFormer) is computed once per input and reused for every diffusion sample.
* Diffusion module, `model/modules/diffusionv2.py`: `AtomDiffusion.sample` (line 295),
  `preconditioned_network_forward` (251), `sample_schedule` (276). The reverse step's noise is off below
  `gamma_min = 1.0` (lines 191, 336). Training noise is `sigma = 16 A * exp(-1.2 + 1.5 N(0,1))` (lines 538-541).
  The default `step_scale` is 1.5 for Boltz-2 (`main.py:1231`), which sharpens samples toward one structure; 1.0 is
  the unscaled sampler.
* Boltz-2 training code is not released ("coming soon" in its README and `docs/training.md`). Only Boltz-1 can be
  retrained with official code.

### 2.3 Why time-lagged data is needed for dynamics

An equilibrium ensemble fixes the energy landscape, not the timescale. Upside gets dynamics by assuming Langevin
motion on an energy trained with contrastive divergence. Its times are in natural units (`README.md:102`), and the
hybrid model's clock is calibrated from a measured lipid diffusion (`up.md:188`). A diffusion model has no energy to
integrate. A lag input means something only if training used pairs of MD frames separated by that lag (an implicit
transfer operator, Schreiner et al. 2023).

### 2.4 Why the Arts et al. 2023 score-as-force method is not used on Boltz

* Arts et al. used one Calpha bead per residue, proteins up to 56 residues and equilibrium MD samples. Their score is
  the gradient of a scalar energy, which makes the force conservative. In their ablation, a direct-vector score was
  about 10x worse on chignolin simulations (TIC JS). Friction and temperature were copied from the MD, and they
  still report no mapping from coarse-grained to real time.
* Boltz's diffusion module fails these requirements:
  * the denoiser outputs per-atom vectors (not conservative) and is not equivariant;
  * of the training noise draws, 6.5% fall below 0.5 A and 0.5% below 0.1 A, so the score at bond-vibration scale is
    poorly trained.
* Plainer et al. 2025 found small-noise scores inconsistent (Fokker-Planck violations) even in purpose-built models,
  and fixed it with a training-time regularizer.

### 2.5 Cost reference points

| Run | Scale | Source |
| --- | --- | --- |
| Boltz-1 structure model | 68k steps, batch 128, crops of 384 to 512 tokens; no GPU count given | Boltz-1 report |
| AlphaFold3 | about 150k steps, batch 256, about 4x Boltz-1's compute | Boltz-1 report |
| Boltz-2 | 97k steps in 4 stages; affinity models on 128 A100s | Boltz-2 preprint |
| AlphaFlow MD fine-tune (full AF2) | 11 h on one 8x A100 machine, about 90 A100-hours | AlphaFlow, Table 2 |
| AlphaFlow PDB stage (made AF2 generative) | 267 h on the same machine; Boltz does not need this stage | AlphaFlow, Table 2 |

## 3. Decisions

1. **Boltz-2 at `b1ebfc4` is the base model.** No BioEmu weights, samples, data or benchmarks; BioEmu2's code may
   be ported (section 1).
2. **Phase 1: independent `md`-tagged Boltz-2 ensembles through the HDX estimator, with no training.** This is the
   reference every modification has to beat.
3. **Phase 2: a renoise-denoise chain on Boltz-2, with no training,** for time-ordered samples with a clock that is
   not calibrated.
4. **Phase 3: a frozen-trunk, lag-conditioned fine-tune of Boltz-2,** only if Phases 1 and 2 fail in a way that MD
   lag data could fix.
5. **No Arts-style score-as-force dynamics on Boltz** (section 2.4).

## 4. Phases

### Phase 0: setup and first test on midway3

Cluster rules are in `~/.claude/midway3.md`. Slurm flags go on the command line:
`--partition=gpu --qos=gpu --account=pi-trsosnic --gres=gpu:N`, with `--time` under the QOS's 36 h cap.

- [x] midway2 cannot host it (checked 2026-10-21 from midway2-login1). Of midway2's GPU partitions `pi-trsosnic`
      reaches only `gpu2`: five K80 nodes (Kepler, compute capability 3.7, which CUDA 12 does not support) and one
      V100 node, midway2-0664, `DOWN+NOT_RESPONDING` since 2024-06-07. The OS is CentOS 7.9 with glibc 2.17.
- [x] midway3 `gpu` partition (sinfo, 2026-10-21): 4 GPUs per node on 5 V100 nodes, 4 RTX 6000 nodes (one down) and
      one A100 node (midway3-0294). QOS `gpu`: at most 36 h, 16 GPUs and 4 nodes per user. RHEL 8.10, glibc 2.28.
      GPU memory and compute capability come from the first jobs' `nvidia-smi`.
- [x] Working directory `/beagle3/trsosnic/yinhan/boltz_hdx` (`$B` below), not `/project`: a venv is tens of
      thousands of files and `/project` metadata is 60x slower (findings 10.14). The project path, should it be
      needed, is `/project/trsosnic/yinhan/`; `~/.claude/midway3.md`'s `/project/trsosnic/yinhanw/` exists but is
      empty.
- [x] Boltz cloned at `b1ebfc4` into `$B/boltz`; `$B/.venv` from module `python/3.11.9` (`$B/install.sh`, login
      node) with `torch 2.14.1+cu126` from the PyTorch CUDA 12.6 index, then `pip install -e "./boltz[cuda]"`. The
      PyPI torch is a CUDA 13 build without V100 support (findings 10.14).
- [x] Weights and CCD in `BOLTZ_CACHE=$B/cache`, downloaded on the login node by `$B/prep.py`.
- [x] MSAs from the ColabFold server, fetched on the login node by `$B/prep.py` through `boltz.main.compute_msa`
      (the call `--use_msa_server` makes): ubiquitin 16,384 sequences; EHEE_rd2_0005, de novo, its own sequence
      only, so it runs effectively single-sequence. The job records whether compute nodes reach the server.
- [x] First test, 2026-10-21: one 1-GPU job per GPU type, about 7 min each (61216915 A100 on midway3-0294,
      61216916 V100 on 0279, 61216917 RTX 6000 on 0283; output `$B/runs/boltz_t0_<gpu>_<jobid>/`, logs
      `$B/logs/`). Peng 2022's two HDX proteins, ubiquitin (1UBQ, 76 residues) and EHEE_rd2_0005 (40 residues,
      de novo), each as a default prediction (5 samples, `step_scale` 1.5, no tag) and an `md`-tagged ensemble
      (100 samples, `step_scale` 1.0). Job script `$B/first_test.sbatch`; `$B/check_samples.py` measures against
      the 1UBQ crystal and the repo's `example/04.HDX/pdb/EHEE_rd2_0005.pdb` and saves each sample's N/CA/C to
      `<run>.npz`.
  * Hardware: A100-PCIE 40 GB (compute capability 8.0), V100-PCIE 16 GB (7.0), Quadro RTX 6000 24 GB (7.5), driver
    535.216.03. All three ran bf16 with the cuEquivariance kernels, so no `--no_kernels` repeat was needed. Peak
    GPU memory 2.4 to 5.0 GB. With one seed the default samples agree across GPU types within 0.01 A in every CA
    metric.
  * Speed (Boltz's prediction time, trunk included, 10 samples in parallel): 100 `md` samples of ubiquitin took
    93 s on the A100, 117 s on the V100 and 122 s on the RTX 6000, about 1 s per sample; EHEE took 54 to 71 s.
    So 10^3 ubiquitin samples cost about 20 V100-minutes and 10^5 about 33 V100-hours.
  * Network: compute nodes reach RCSB and Hugging Face (HTTP 200). ColabFold's root URL gave 404, which shows the
    host answers but not that MSA requests work; MSAs stay a login-node step.
  * Default samples are near-native with clean geometry: CA-RMSD to the reference median 0.78 A (ubiquitin) and
    0.93 A (EHEE); peptide C-N 1.327 +- 0.008 A; at most one heavy-atom pair under 2.5 A per sample.
  * `md` samples are wider: ubiquitin pairwise CA-RMSD mean 2.0 A (max 4.7), CA-RMSD to 1UBQ median 1.6 A (max
    3.0), CA RMSF up to 7.7 A at Gly76; EHEE pairwise mean 0.9 A. Their local geometry is worse. Peptide C-N sd
    is 0.054 to 0.059 A, against 0.008 in the default runs and 0.014 in crystal structures (Engh and Huber 1991).
    N-CA and CA-C sd is 0.035 A, against 0.004. About 7% of C-N bonds fall outside 1.25 to 1.45 A, spread along
    the whole chain. On the V100, 2 of 100 ubiquitin samples have a C-N bond outside 1.05 to 1.65 A.
    86 to 97 of every 100 samples have a heavy-atom pair under 2.5 A, up to 10 pairs; most are over-short polar
    pairs (next test).
  * The `md` runs changed the tag and `step_scale` (1.5 to 1.0) together, so the cause of the scatter is unknown.
    Whether it changes protection flags is unmeasured.
- [x] Tag by step-scale test, 2026-10-21 (61231966, V100 on midway3-0280, 18 min; `$B/tag_stepscale_test.sbatch`,
      output `$B/runs/boltz_tag_stepscale_v100_61231966/`): {no tag, `md`} x {`step_scale` 1.5, 1.0}, 100 samples
      each at seed 2, both proteins. Close contacts were classified by atom pair with a local script on the copied
      samples.

  | cell | C-N sd (A), ubq / EHEE | N-CA sd | pairwise CA-RMSD (A), ubq / EHEE | samples with a contact < 2.0 A, ubq / EHEE |
  | --- | --- | --- | --- | --- |
  | no tag, 1.5 | 0.008 / 0.009 | 0.004 | 1.28 / 0.26 | 3 / 2 |
  | `md`, 1.5 | 0.008 / 0.009 | 0.004 | 1.57 / 0.45 | 1 / 5 |
  | no tag, 1.0 | 0.031 / 0.036 | 0.018 | 1.61 / 0.52 | 15 / 17 |
  | `md`, 1.0 | 0.055 / 0.054 | 0.035 | 1.99 / 0.88 | 25 / 27 |

  * `step_scale` 1.0 causes the bond-length scatter. At 1.5 both tags give clean bonds; at 1.0 the `md` tag
    nearly doubles it. The `md` cell at 1.0 reproduces the first test with a new seed.
  * The `md` tag widens the ensemble at both step scales. `md` at 1.5 is almost as wide as no tag at 1.0, with
    clean bonds.
  * Most contacts under 2.5 A are polar pairs: Lys NZ or Arg NH with Asp/Glu carboxylate O, and backbone N-O.
    These are salt bridges and H-bonds drawn 0.3 A or more too short (N-O is about 2.8 A or longer). Neither reference structure has a contact
    under 2.5 A. One steric overlap stands out: in no-tag, 1.5 ubiquitin model 33 the Arg42 and Arg72
    guanidinium groups sit on top of each other (CZ-CZ 0.69 A, 15 contacts); the same pair comes within 1.75 A
    in no-tag, 1.0 model 41.
  * This is the conflict the test was set up to find. The unscaled sampler (1.0) is the one whose samples follow
    the model's distribution, and it has broken local geometry. `step_scale` 1.5 has clean geometry but
    sharpens the distribution toward its main state, which biases the populations HDX reads toward protection.
- [x] Phase 1 sampler setting: `md` at `step_scale` 1.0, because its bond-length scatter does not reach the
      protection flags (2026-10-21, local, no GPU). Every set from the tag by step-scale test went through
      `py/get_protection_state.py` unchanged, with configs from `example/04.HDX/1.config.py`'s recipe (ff_2.1,
      heavy-atom coverage) built from the two reference structures. The engine had to come from a scratch tree
      (py/ copied and obj/ built at repo `6d295b0`), because the repo's `obj/` is older than its `py/` and lacks
      `get_param_size`. Files in `scratchpad/boltz_hdx/flags/` (`run_flags.sh`, `samples_to_traj.py`,
      `compare_flags.py`). The estimator reads only N, CA and C: Upside places its own side-chain beads from the
      backbone, so Boltz's side chains, with their short salt bridges and the Arg42/Arg72 overlap, do not enter.
  * Geometry noise alone: the `md`, 1.5 samples plus isotropic noise per atom (0.038 A, which gives C-N sd
    0.054 A, the 1.0 spread). On the same frames the noise flips 0.2% (ubiquitin) and 1.3% (EHEE) of the
    combined flags; the H-bond flag alone flips 2.8% and 0.8%, mostly in amides that burial already protects.
    Per-residue protected fractions move by at most 0.03, mean 0.002. That is below the binomial error of 100
    samples, 0.02 to 0.05.
  * `md` at 1.0 against `md` at 1.5, the whole difference: mean |dp| 0.019 (ubiquitin) and 0.016 (EHEE), max 0.16
    and 0.12. At 1.0, fewer amides are never seen open (p = 1): 45 against 55 of 72 for ubiquitin and 24 against
    29 of 39 for EHEE.
  * Every set stays near the native flags: the reference structures protect 59 of 72 and 32 of 39 amides, and
    96 to 98% of sample flags agree with them. With 100 samples, every amide protected in the native state is
    censored at RT ln(100) = 2.7 kcal/mol.

### Phase 1: independent Boltz-2 ensembles, no training

HDX:
- [ ] `md`-tagged ensembles of 10^3 to 10^4 samples per protein (20 V100-minutes to 3.3 V100-hours for ubiquitin,
      Phase 0), at `step_scale` 1.0 (Phase 0), for ubiquitin and EHEE_rd2_0005, then other soluble proteins with
      HDX data and Upside results.
  * First job queued 2026-10-21 (61254153, one V100, `$B/phase1_hdx.sbatch`; `remote_jobs.md` §1): 1,000 samples
    per protein at seed 3.
- [ ] Protection flags from every sample's N/CA/C through the repo's protection-state code, uniform weights; dG per
      residue against experiment and against Upside. An amide with p = 1 is censored at RT ln(N), not a value.

mdCATH 320 K (CHARMM22*/TIP3P, 5 replicas of about 500 frames at 1 ns), equilibrium only:
- [ ] Test set: 8 to 12 domains from the 89-domain panel (`ff3_selection/` on the cluster; `plan.md` Phase 19), over
      the 50-100, 101-150 and 151-250 residue bins, plus 3 to 5 domains that switch states at 320 K (multimodal
      Calpha-RMSD or TICA distribution across replicas). Every domain must come from Boltz-2's held-out mdCATH
      clusters, since its training included mdCATH; if that split cannot be found, mark the results as possibly
      optimistic.
- [ ] The panel's local copies keep backbone atoms only. Fetch heavy atoms for the chosen domains, reusing the
      range-request approach of `mdcath_fetch.py`.
- [ ] Metrics on Calpha or heavy atoms against the 5 MD replicas: per-residue RMSF; JS divergence of pairwise-distance
      and Calpha-RMSD distributions; TICA fitted on MD with Boltz samples projected onto it, and the free-energy
      surface on TIC1 and TIC2.

Gate 1:
* dG tracks experiment at least as well as Upside, and the mdCATH distributions match: equilibrium needs no
  training. Phase 2 only if decorrelation or the EX1/EX2 question needs time order.
* Ensembles too narrow (censored dG where experiment is finite, RMSF below MD): Phase 2, then Phase 3.
* Ensembles in the wrong states, not only too narrow: reconsider the approach before training.

### Phase 2 (optional): renoise-denoise chain on Boltz-2, no training

- [ ] New sampler loop, each step:
      `y = x_k + sigma_j * eps`; run the reverse sampler from `sigma_j` to 0 with the `md` tag and
      `step_scale=1.0`; Kabsch-align the result to `x_k`. Compute the trunk once and reuse it.
      Code hooks: `AtomDiffusion.sample`, `preconditioned_network_forward`, `sample_schedule` (section 2.2).
- [ ] Keep the reverse step stochastic below 1 A. Boltz switches its noise off below `gamma_min=1.0`, and the
      stochastic step is what gives detailed balance:
      `p(x) K(x'|x) = integral p_sigma(y) p(x|y) p(x'|y) dy` is symmetric in x and x'. A deterministic step
      keeps the distribution stationary but is not guaranteed reversible.
- [ ] Sweep `sigma_j` (0.5, 1, 2, 4, 8 A) and the number of reverse steps. Check that the chain's distribution
      matches Phase 1's independent samples.
- [ ] Dynamics against mdCATH with a fitted clock: a Markov state model on MD TICA clusters, per-row JS divergence
      of the transition matrices at matched lag (Arts' metric), implied timescales up to one scale factor, a
      Chapman-Kolmogorov check, and the physicality pass rate along the chain.
- Limit: the clock is arbitrary, and the chain only visits states Boltz-2 already samples.

### Phase 3 (conditional on Gate 1): frozen-trunk, lag-conditioned fine-tune of Boltz-2

- [ ] Blocker: no Boltz-2 training code. Either wait for its release or write a training loop for the diffusion
      module only. The EDM loss weight and noise distribution are already in `diffusionv2.py` lines 532 to 542.
- [ ] Lag conditioning (Schreiner et al. 2023), ported from BioEmu2 at `f99e0d8` (`github.com/microsoft/bioemu2`,
      MIT; code only, section 1). Paths relative to its `src/bioemu2/`, read 2026-10-21:
  * `x_0` projection into the atom queries: `cond_to_q_trans = LinearNoBias(3, atom_s)`
    (`biomoa/modules/encoders.py:160`, applied at 314). Initialize it to zero here, so the fine-tune starts from
    Boltz-2 unchanged;
  * a gate on `x_0` that decreases with lag, so long lags fall back to the equilibrium distribution:
    `g(tau) = sum_k softplus(w_k) exp(-tau/tau_k)`, `MultiExponentialTimeEmbedding` (`biomoa/modules/time.py`),
    weights initialized to 1/n (`encoders.py:163-166`, 318-320);
  * the lag embedding in the single conditioning: `DampedHarmonicTimeEmbedding` (`time.py`; 32 log-spaced
    timescales from 1 ps to 1 ms, each with 4 oscillation frequencies), projected and added at
    `biomoa/modules/diffusion.py:417-418`, 450-455;
  * a per-sample flag `is_dynamics_conditioned` that zeroes both inputs, so lagged and unconditioned examples train
    together (`encoders.py:322-324`, `diffusion.py:452`);
  * the same random rotation and translation applied to `x_0` and to the noisy coordinates at every denoising step
    (`denoising/edm.py:136-149`), and endpoints aligned to `x_0` (`propagation.py:471`).

  BioEmu2 releases no training code either, so the blocker above stands.
- [ ] Data: mdCATH at 320 K (5,398 domains, 5 replicas of up to about 500 ns, 1 ns frames) gives lags from 1 ns to
      a few hundred ns, and the learned kinetics inherit CHARMM22*/TIP3P. Map frames to Boltz's heavy-atom order.
      Hold out the Phase 1 domains and their CATH superfamilies.
- [ ] MSAs for all domains (ColabFold server; days of CPU or wall time, no GPU).
- [ ] Precompute trunk embeddings, one pass per domain: tens of GPU-hours (*estimate*). Storage at 128 bf16 pair
      channels is about 6 MB at 150 residues and 23 MB at 300, 50 to 150 GB in total.
- [ ] Pilot on about 200 domains: one 4-GPU node, one 35.5 h job, about 150 GPU-hours. Check that the loss depends
      on `x_0` and the lag (ablate each at evaluation) and measure seconds per step.
- [ ] Full run: 100 to 1,000 GPU-hours (*estimate*), that is 1 to 7 chained 36 h jobs on one node, with
      checkpoint/restart chaining like `upside2-md/training/train_chain.sbatch`. Re-estimate after the pilot.
- [ ] Evaluate with the Phase 1 and Phase 2 metrics, against Phase 1's untrained Boltz-2 on the same domains.

## 5. Known risks and blockers

* No energy and no temperature input: Boltz populations cannot be reweighted, and dG resolution is RT ln(N) at one
  condition (section 2.1).
* HDX in stable proteins comes from rare openings on millisecond-to-second timescales. mdCATH's 320 K runs of up to
  about 500 ns do not contain them, so a Phase 3 fine-tune on that data cannot learn their populations.
* Boltz-2 trained on mdCATH: Phase 1 domains must come from its held-out clusters (Phase 1).
* The 89-domain panel holds stable native states, which makes it a weak test of transitions (Phase 1 adds switching
  domains for this).
* `step_scale` 1.0 gives bond-length scatter 2 to 4x the crystal spread. It does not change the protection flags
  (Phase 0), but any use beyond the HDX estimator, such as Phase 2's chain or comparisons of local geometry with MD,
  has to deal with it. Boltz-2 also draws salt bridges and H-bonds 0.3 A or more too short at every setting.
* Compute nodes reach RCSB and Hugging Face; whether they can query ColabFold is untested, so MSAs are fetched on
  the login node.
* Phase 3 is blocked on Boltz-2 training code, or on writing a diffusion-only training loop.

## 6. References

* Boltz: https://github.com/jwohlwend/boltz ; Boltz-1 report https://gcorso.github.io/assets/boltz1.pdf ;
  Boltz-2 preprint https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1.full
* Peng et al. 2022, Upside HDX (ubiquitin, EHEE_rd2_0005); findings 5.1
* Arts et al. 2023, Two for One: https://arxiv.org/abs/2302.00600
* Plainer et al. 2025, Consistent Sampling and Simulation: https://arxiv.org/abs/2506.17139 (code: noegroup/ScoreMD)
* Jing et al. 2024, AlphaFlow: https://arxiv.org/abs/2402.04845
* Schreiner et al. 2023, Implicit Transfer Operator Learning (NeurIPS 2023)
* BioEmu2 code (MIT, used for code only): https://github.com/microsoft/bioemu2 at `f99e0d8`
* Engh and Huber 1991, Acta Cryst. A47, 392 (peptide bond geometry)
