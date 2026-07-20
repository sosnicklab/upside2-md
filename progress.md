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
