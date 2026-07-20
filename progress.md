# Progress

## Standard-step dry-MARTINI friction clock (2026-07-20)

- Verified `example/02.ReplicaExchangeSimulation` uses the standard `.009` numerical timestep and set both
  protein and bilayer to that same step.
- Derived the factor-four molecular target: 40 ps physical protein time equals 10 ps raw MARTINI time and
  requires a whole-DOPC lateral MSD slope of `0.184 A^2/step` for `D=11.5 um^2/s`.
- Scanned five one-step frictions. Molecular DOPC diffusion saturated near `0.175 um^2/s`, 66x below target,
  while lower friction heated the interface; rejected molecular-diffusion matching on this path.
- Implemented the accepted fallback: native 4 ps dry-MARTINI relaxation becomes 16 ps physical, hence
  `tau_up=.0036` and `alpha_i=m_i/tau_up`. Production applies this FDT bath to environment particles and
  protein N/CA/C carriers initially inside the existing 12 A lipid interaction range.
- Removed lipid substeps, the overdamped branch, displacement caps, empirical friction multipliers, and dead
  flags. Positive-g-JF-friction particles are excluded from OU; noncontact carriers retain standard OU.
- Kept SC-env and BB-env active and spline-only. O/BB are differentiable virtual sites with zero momentum and
  BB-env gradients mapped through the complete N/CA/C Jacobian.
- Added explicit H5 clock/contact metadata and `analyze_dopc_diffusion.py`, which reconstructs mass-weighted
  14-bead COMs, unwraps PBC, removes bilayer drift, and reports molecular D beside the target.
- Corrected preparation staging and constraint handling: ordinary damping before production; fixed and
  z-fixed particles are honored inside g-JF. Targeted PO4 z-hold and O/BB momentum smokes pass exactly.
- Fresh final-code 5000-step workflows passed for 1RKL and 1AFO. DSSP helix counts stayed 7 and 27; mean
  protein/lipid kinetic energies were within 0.8% of `3kT/2`. Measured molecular D was about `7e-5 um^2/s`.
- Separate 50,000-step production continuations kept both helix counts and target kinetic temperature.
- Final verification passed: C++ build, Python compile, shell syntax, diff whitespace, synthetic diffusion
  regression, H5 metadata audit, fresh workflows, and LaTeX compilation.

Files modified: `src/martini_brownian.{h,cpp}`, `src/deriv_engine.cpp`, `src/thermostat.{h,cpp}`,
`src/martini_hybrid.cpp`, `src/martini_internal.h`, `src/martini_potential.cpp`, `src/main.cpp`, `src/martini.h`,
`py/martini_prepare_system.py`, `py/martini_prepare_system_lib.py`, `example/16.MARTINI/run.py`,
`example/16.MARTINI/run_sim_hybrid.sh`, `example/16.MARTINI/analyze_dopc_diffusion.py`, and the manuscript.
