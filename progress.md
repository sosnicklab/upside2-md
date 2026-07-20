# Progress

## Production-freeze correction (2026-07-20)

- Reproduced the user-visible freeze in both stage-7 VTF files. DOPC molecular COM motion was only
  `0.0081 A` per saved frame and `0.17--0.24 A` net RMS over 1,001 frames.
- Identified the cause: production `tau_up=.0036` produced mass-6 friction `1666.7` and
  `alpha*dt/(2m)=1.25`, suppressing coordinates despite an apparently normal kinetic temperature.
- Rejected native-damping, old low-friction, and direct protein bare-friction controls because they moved
  DOPC but overheated the interface or reduced 1RKL helicity.
- Implemented the factor-four-corrected bare-particle mobility clock (`alpha_bead=0.1691804`) at the shared
  `.009` timestep, with no lipid substeps. Protein N/CA/C carriers receive additive friction from DOPC beads
  inside the existing 12 A spline cutoff; zero-contact carriers retain the Upside OU bath.
- Recompute static contact counts from actual coordinates after stage handoff, minimization promotion, and
  continuation. SC-env, BB-env, hard cores, virtual-site force projection, and spline tables are unchanged.
- Kept preparation separate from the kinetic claim: softened stages use native damping and full-core
  stages use strong FDT overlap-settling damping; production alone uses the particle clock.
- Two 50,000-step regressions completed. DOPC COM net RMS was `3.43 A` (1AFO) and `3.72 A` (1RKL), molecular
  D was `0.0132` and `0.0152 um^2/s` (still below target), and total kinetic energy was within 1% of target.
  1AFO retained 54/54 helical residues; the 1RKL transmembrane core at residues 10--29 remained helical.
- Physically valid fresh 1AFO and 1RKL workflows completed minimization, equilibration, burn-in, production,
  and VTF extraction with finite moving coordinates and the expected phase/contact friction metadata.
- Updated the workflow description and manuscript to state the verified particle-level observable and the
  unmatched molecular diffusion explicitly.

Files modified: `py/martini_prepare_system.py`, `py/martini_prepare_system_lib.py`,
`example/16.MARTINI/run_sim_hybrid.sh`, `example/16.MARTINI/drymartini_upside_interface.tex`, `plan.md`,
`findings.md`, and `progress.md`.

## HDX workflow audit (2026-07-20)

- Traced the full structural-protection-to-uptake path. It MBAR-averages binary protected states and uses an
  EX2-like analytic rate model; simulation time is not used as the uptake clock.
- Confirmed that the standard protein-only extractor cannot load hybrid trajectories (coordinate-count mismatch)
  and that the available hybrid-aware extractor is not wired into the analysis driver.
- Verified saturation of the current hybrid membrane rule on stage-6.6: all 29 1RKL donors are always protected,
  while 65.5% are H-bond-open; downstream uptake would therefore be dominated by the `k_chem/1000` floor.
- Found no hybrid `PS.npy`, percentD, or HDX artifacts in `example/16.MARTINI`; there is no existing calculated
  hybrid HDX result to validate.
- Reassessment: lipid/protein clock mismatch is indirect for this equilibrium estimator, but frozen or short,
  poorly decorrelated trajectories remain unsuitable. Current outputs support qualitative structural inspection,
  not quantitative HDX prediction.

Files modified for this audit: `plan.md`, `findings.md`, and `progress.md` only.

## 1RKL temperature and hybrid-HDX reassessment (2026-07-20)

- Confirmed production `T_up=0.8647` equals 303.15 K and is passed consistently to g-JF and OU baths. Last-window
  protein/lipid kinetic temperatures are within 1.1%/0.3% of target, excluding a temperature-conversion error.
- Recomputed DSSP over all 1,001 stage-7 frames. The 10--29 core averages 84.5% helix occupancy, briefly drops to
  four helical residues near 1.09 us, and ends at 19/20; core CA RMSD stays below 1.58 A.
- Corrected the prior physical-duration arithmetic: 50,000 numerical steps at the declared 40 ps/step are 2 us,
  not 18 ns.
- Identified timestep convergence—not integrator-specific temperature scaling—as the unresolved numerical check;
  prior evidence favors `.00225` over `.009` for the coupled hard interface.
- Defined the HDX adapter: map hybrid N/CA/C into a protein-only HDX analysis engine, preserve full-system energy
  and temperature for ensemble weighting, and treat membrane water accessibility as a separately calibrated term.
- Found a workflow documentation hazard: `T.npy` must contain Upside `kT`, not Kelvin as stated in the README.

Files modified for this reassessment: `plan.md`, `findings.md`, and `progress.md` only.

## HDX reuse architecture correction (2026-07-20)

- Recast the hybrid path as one H5 projection adapter rather than a separate HDX implementation.
- The adapter will map hybrid N/CA/C into a standard protein-only trajectory view while retaining full-system
  potential and Upside temperature for the existing MBAR analysis.
- `example/04.HDX` remains the method reference; `example/00.AnalysisScripts` remains the maintained uptake and
  experiment-comparison pipeline. Existing steps run unchanged after projection.
- Any membrane accessibility correction remains a transparent postprocessing layer with the stock PS retained.

## Hybrid HDX implementation and verification (2026-07-20)

- Added a hybrid H5 projector and optional calibrated water-accessibility combiner; stock protein PS is retained.
- Wired hybrid configuration through `00.AnalysisScripts` and added `16.MARTINI/run_hdx_analysis.sh`.
- Kept steps 1--6 and the `04.HDX` MBAR/EX2 method as the analysis path; fixed scalar-temperature handling and
  numerically stable log-weight normalization reached by the single-replica hybrid case.
- Verified exact mapped coordinates and exact full potential/H-bond/temperature/time preservation for 1,001
  frames. The full temperature-matched 1RKL workflow completed uptake, stability, and summary generation.
- Quantitative trust gate failed: first/second-half PS differs by up to `0.621`, minimum effective sample count is
  `3.74`, DOPC diffusion is `0.0152` versus `11.5 um^2/s`, and no calibrated membrane water accessibility exists.

Files added: `py/martini_hdx_project.py`, `py/combine_hdx_protection.py`, and
`example/16.MARTINI/run_hdx_analysis.sh`. Updated the `00.AnalysisScripts` driver, steps 2--4, stability helper,
both READMEs, `plan.md`, `findings.md`, and `progress.md`.
## 2026-07-20 — 1RKL temperature/BB-force audit paused for rerun

- Verified that the cited stage-7 trajectory ran at `T_up=.8647`, not `.80`; the current script's `.80` default
  was applied after that output was generated.
- Compared particle--BB tables, type wiring, pair admission, spline derivative, and virtual-site force propagation
  with `b1041bb`; found no raw pair-force regression and identified the then-current full-Jacobian route and
  historical timestep as material differences. The later correction restores the historical partial route.
- Made no physics or workflow-code changes. Per user direction, the next gate is a fresh `.80` simulation and
  residue-wise DSSP/provenance check.

## 2026-07-20 — Friction calibration and trust manuscript rewrite

- Rewrote drymartini_upside_interface.tex around one code-matched calibration derivation, including the
  factor-four mapping, numerical bead mobility/friction, contact-local protein drag, H5 audit fields, and
  preparation-versus-production distinction.
- Replaced broad thermodynamic/HDX trust claims with an evidence table and explicit implementation,
  configurational, and kinetic validation layers. Documented the molecular-diffusion failure, unvalidated protein
  friction, timestep/cutoff concerns, and current HDX convergence/accessibility failures.
- Corrected the Coulomb description from shifted to abruptly truncated and identified its quantitative validation
  consequence.
- Verification passed: exact arithmetic assertions, manuscript consistency checks, two-pass warning-free
  pdflatex, chktex, and git diff --check.
- Detected and re-audited freshly replaced stage-7 outputs. Updated the manuscript to the unified $T_\mathrm{up}=.8$
  artifacts and their new kinetic-temperature, DOPC COM diffusion, and DSSP diagnostics.

## 2026-07-20 — Unified hybrid temperature configuration

- Changed `run_sim_hybrid.sh` to define one `TEMPERATURE` and use it for both the Upside thermostat and DOPC
  friction calibration; removed the independent bilayer reference-temperature default/override.
- Left particle--BB forces, timestep, conservative tables, and all interface interactions unchanged pending the
  fresh `.80` simulation.
- Verification passed: `bash -n`, `git diff --check`, a static single-source contract assertion, and an execution
  trace showing `TEMPERATURE=.77` overrides a conflicting reference setting and exports the reference as `.77`.

## 2026-07-20 — 1RKL BB-force and stage-7 handoff correction

- Restored the historical BB reverse route required by Upside's regenerated-O cycle: N/CA/C receive
  14/54, 12/54, and 12/54 of the BB gradient, while sensitivities on regenerated O and BB are discarded.
- Localized helix loss to the first stage-7 minimization: stages 6.0--6.6 retain 23 helical residues, while the
  existing stage-7 trajectory starts at 13.
- Replayed the handoff from identical coordinates. The unrestrained minimization left 9 helical residues, an
  early spring-10 restraint left 17, and hard-fixing protein coordinates retained all 23.
- Replaced the intermediate fixed-plus-spring protocol with one rigid `production_handoff`: SC-env and BB-env
  are active while the complete protein remains a rigid body through minimization and burn-in. Flexible dynamics
  start only after the explicit `production` relabel.
- A reduced end-to-end workflow completed the new sequence. The maximum internal pair-distance change across all
  93 persistent N/CA/C carriers through handoff minimization and burn-in was `4.58e-5 A`. The shortened membrane
  equilibration makes this a wiring/invariant test only; a full fresh 1RKL run remains the structural gate.
- Final verification passed: exact H5 weight assertions, removed stage-7 spring/release controls, Python and shell
  syntax, full C++ build, two-pass warning-free TeX compilation, and `git diff --check`.

Files modified for this correction: `src/martini_hybrid.cpp`, `py/martini_prepare_system.py`,
`example/16.MARTINI/run_sim_hybrid.sh`, `example/16.MARTINI/readme.md`,
`example/16.MARTINI/drymartini_upside_interface.tex`, `plan.md`, `findings.md`, and `progress.md`.
