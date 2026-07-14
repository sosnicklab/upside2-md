# Project Goal

Preserve the committed dryMARTINI force field and stable protein dynamics while
correcting CGL transport to the Upside clock: one integration stage is
`dt=0.004 T_up = 40 ps`, CGL and CGL-SC forces are evaluated on that same
stage, and the measured CGL lateral diffusion is interpreted with the
`14 * 4 = 56` coarse-graining correction.

# Architecture & Key Decisions

- Use commit `52f637e` as the behavioral baseline. Its CGL-CGL table, tempered
  Boltzmann projection, hybrid interaction tables, GLE implementation, and
  protein behavior remain unchanged.
- Make transport-only changes: set every production integration stage to
  `0.004 T_up` and calibrate only the permitted dryMARTINI mass/GLE transport
  against raw `D≈0.25 A^2/T_up`. Keep protein carrier masses at exactly
  `1 m_up`.
- Retain the committed two-mode CGL GLE (`tau=0.2,2.0`, coupling
  `0.30375,0.2205`). It acts only on CGL/environment particles; the Upside
  protein thermostat and force field are not slowed or rescaled.
- Preserve all CGL-CGL, SC-CGL, BB-CGL, and environment interactions. They are
  evaluated from the committed H5 spline tables on every global integration
  stage.
- Keep the required `25 T_up` tempered Boltzmann projection. Removing it
  causes CGL clustering.
- Do not add orientational restraints, force caps, exclusions, fitted membrane
  terms, explicit fork states, or new force-field reconstruction machinery.
- Start fresh after a mass/timestep change. Continuations must validate the
  stored transport signature and must not carry momentum across a mass change.
- Leave installed H5 files unchanged unless a force-field defect is separately
  demonstrated after the corrected transport baseline passes or fails.

# Execution Phases

- [x] Restore active Python/C++ force-field and hybrid code to commit
  `52f637e`, preserving only the minimal transport configuration.
- [x] Build and run focused static checks: unchanged table hashes, protein mass
  `1`, accepted CGL mass/GLE metadata, and stage `dt=0.004`.
- [x] Run a fresh protein-free CGL membrane and measure lateral diffusion,
  bilayer structure, clustering, and finite energies.
- [x] Run fresh short 1rkl and 1afo trajectories and check secondary structure,
  hbonds, RMSD, protein-bilayer angle, bilayer integrity, and finite dynamics.
- [x] Extend accepted runs through stages 7.0-7.3 and inspect trajectories and
  logs for ejection, angle drift, leaflet-gap bending, or explosion.
- [x] Update `example/16.MARTINI/cg_lipid_potentials.tex` only after the model
  and simulations are accepted.

# Known Errors / Blockers

- No current blocker.

## Revised Decisions

- Commit `9fb81b5` is the historical exact-step transport source, and its GLE
  kernel/schedule are already present in `52f637e`. Its `0.012` mass cannot
  be transplanted unchanged to the later pair-relaxed spline: a fresh 404-CGL
  run gives raw `D=0.910 A^2/T_up`, leaflet separation
  `16.4 -> 24.9 A`, and 34% maximum leaflet reassignment.
- Promote the stable `0.05` branch as the active default. On the fresh
  one-step 404-CGL trajectory, multi-origin x56 diffusion converges to
  `11.47/11.45 um^2/s` over the `10-40/20-60 T_up` windows with
  `R2 > 0.998`. This matches the reported `11.5 um^2/s` DOPC value at 303 K
  while the bilayer remains finite, closed, ordered, and laterally uniform.
- Reject `0.075`: its matched early probe remains too fast
  (x56 `30.9 um^2/s`). Promote `0.05` as the workflow default.

# Review

- The accepted implementation is the committed force field plus the transport
  delta: one global `0.004 T_up` Verlet step, CGL mass scale `0.05`, and the
  committed two-mode GLE. Protein mass and protein dynamics are unchanged.
- Fresh 404-CGL, 1RKL, and 1AFO runs pass diffusion, finite-energy, membrane,
  secondary-structure, and protein-angle gates. Stage 7.2 is not vertical and
  stage 7.3 does not explode.
- The methods paper was rewritten around the accepted model and fresh results;
  it compiles without warnings and all seven rendered pages were inspected.
