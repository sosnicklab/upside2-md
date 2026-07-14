# Findings

## Accepted Baseline

- Commit `52f637e` is the stable protein/force-field baseline. Its recorded
  stage-7 runs retained secondary structure (1rkl `0.968`, 1afo `0.986`),
  stayed finite through stage 7.3, and kept the 1rkl protein-bilayer angle near
  its initial value.
- The committed implementation already contains idempotent CGL mass scaling,
  the accepted two-mode CGL GLE, zero coarse stage-7 restraints/burn-in, the
  pair-relaxed CGL-CGL table, and the required 25-temperature tempered
  Boltzmann projection.
- The installed committed CGL table hash is
  `ba7049317e69a371812c69da8876edaa4ce81c21bef4c75a82a761176862b46c`.
- Commit `9fb81b5` is the historical one-step exact-GLE calibration. It uses
  CGL mass scale `0.012`, rotational tau `0.008`, GLE taus
  `0.2,2.0`, and couplings `0.30375,0.2205`. Its documented 288-CGL
  target is raw `D≈0.25 A^2/T_up`, or approximately
  `14 um^2/s` after the required factor of 56.
- That historical mass is not transferable to the later pair-relaxed
  `52f637e` conservative table. At `dt=0.004`, a fresh 404-CGL
  `mass_scale=0.012` trajectory gives raw `D=0.910 A^2/T_up`, opens the
  leaflets from `16.4` to `24.9 A`, and reaches 34% maximum leaflet
  reassignment. The GLE architecture remains valid; the mass value is rejected.
- The accepted one-step `mass_scale=0.05` 404-CGL run records exactly
  `120 T_up` for 30,000 steps at `dt=0.004`. It is finite, keeps leaflet
  separation near `16.5 A`, has `p05(|n_z|)=0.977`, and has no meaningful
  leaflet mixing. Multi-origin x56 diffusion is `12.40`, `11.47`, and
  `11.45 um^2/s` over `5-20`, `10-40`, and `20-60 T_up`, respectively;
  the latter two fits have `R2 > 0.998`.
- Physical comparison supports that range: a pure-DOPC study reports about
  `9.32 um^2/s` at 298 K (DOI `10.1039/D0CP05111J` cites the experimental
  value), while Lipid14 reports DOPC at 303 K as `6.48 um^2/s` in NPT,
  `9.49 um^2/s` in NVE, and cites experimental `11.5 um^2/s`
  (DOI `10.1021/ct4010307`). The accepted CGL trajectory brackets these
  values after the required x56 correction.
- Fresh 320-ns hybrid trajectories pass through stage 7.3 without nonfinite
  coordinates or energy. 1RKL retains mean/final DSSP `0.925/0.903`, mean
  angle `35.9 deg` versus `34.4 deg` initially, and a closed `15.19 A`
  bilayer. 1AFO retains mean/final DSSP `0.989/1.000`, both chain angles within
  `12.8-18.9 deg`, and a closed `15.31 A` bilayer.
- The installed `dopc.h5`, `particle.h5`, and `sidechain.h5` hashes are
  unchanged. There is no active C++ or force-field-builder diff.

## Physical Interpretation

- The timestep belongs to the global Upside integration stage. CGL-CGL and
  CGL-SC conservative forces therefore use the same `40 ps` stage as the
  protein backbone; there is no separate interaction timescale.
- CGL transport may be changed with dryMARTINI mass/GLE parameters, but protein
  mass, protein force constants, protein thermostat, and Upside time mapping
  must remain unchanged.
- The `14 * 4` factor corrects the measured one-CGL lateral diffusion for the
  molecular coarse-graining represented by one CGL particle. It must be
  applied to the measured diffusion, not to the simulation clock.

## Lessons

- When the committed version already passes the protein and force-field
  behavior gates, preserve it exactly and isolate the transport delta. Do not
  redesign the force field without a failing diagnostic obtained from that
  committed model under the corrected transport settings.
- Treat visual reports of a frozen bilayer as a transport-measurement problem
  first. Verify stage duration, mass scaling, GLE scope, and lateral MSD before
  changing conservative interactions.
- Never slow the Upside core to match dryMARTINI. Transport tuning applies only
  to dryMARTINI particles, while CGL-SC forces remain evaluated every global
  stage.
- Do not infer a transport parameter from a trajectory until the executable's
  requested-step to integrated-step mapping is verified. The pre-reset
  `0.05` estimate was entangled with a three-step scheduler; only the fresh
  one-step `0.05` validation is accepted.
- Historical transport parameters must be revalidated after a conservative
  table changes. Preserve the later stable table, then calibrate only the
  explicitly permitted dryMARTINI transport variable against diffusion and
  structure gates.
