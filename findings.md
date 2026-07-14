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
- The accepted coarse transport settings do not by themselves validate the
  full 14-bead workflow. A fresh full 1AFO run with protein mass `1` and global
  production `dt=0.004` remains stable through burn-in step 30,450, then shows
  an abrupt protein-energy spike above `+10,000` and sustained hbond collapse.
  That trajectory is rejected and must not be continued.
- The former stage-7 continuation path was not an exact restart: after promoting
  final position and momentum, it ran minimization, moved positions by up to
  `0.83 A`, and invalidated momentum. Repeated segments drove 1RKL's membrane
  angle from about `34 deg` to `6 deg`. Production continuation must go directly
  from state promotion/interface refresh to MD.
- Removing continuation minimization restores exact position/momentum and
  compaction boundaries. It does not explain the current 1RKL rotation: that
  checkpoint becomes vertical within its original stage 7.0 and is therefore
  independently rejected.
- A rigid whole-protein tilt scan of the complete stage-7 Hamiltonian has its
  minimum at the initial 1RKL tilt (`34.15 deg`); rotating toward vertical
  raises the potential by about `61 E_up`. The observed rotation is therefore
  not an equilibrium preference in the initial conservative potential.
- The installed `dopc.h5`, `particle.h5`, and `sidechain.h5` hashes are
  unchanged. There is no active C++ or force-field-builder diff.
- Reconstruction audit of `dopc.h5` confirms that CGL-CGL uses two isolated
  references, a `140 x 9 x 9` tensor and pair-tail relaxation; SC-CGL uses
  eight references, 29 placement-bead rows, `25 x 11 x 11` controls and
  positive-overlap tail relief; CGL-target uses twelve references and
  `91 x 321` controls with no compaction overlay. The CGL-CGL force-match,
  bilayer-PMF, and IBI counters are zero.

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

- A valid completed stage-7 checkpoint takes precedence over rerunning hybrid
  packing and equilibration. Before restarting a named production workflow,
  inspect its latest exact `PDB.stage_7.N.up` file and enable continuation
  explicitly; never infer continuation merely from an existing output folder.
- Do not force a legacy checkpoint through current validation. If particle
  masses or the global timestep changed, preserve its files for provenance but
  regenerate a compatible production stage before enabling exact momentum
  continuation.

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
- A methods paper must make the force field independently reconstructible.
  Naming tables and reporting validation is insufficient; document source
  ensembles, coordinates, statistical projection, relaxation optimization,
  spline fitting, unit conversion, and runtime force/torque evaluation.
