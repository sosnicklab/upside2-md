# Findings

## External / Technical Findings
- 2026-06-22: Coarse hybrid salt and CGL z-placement fix.
  - Missing salt ions in `run_sim_1afo.sh` and `run_sim_1rkl.sh` came from
    the coarse branch of `run_sim_hybrid.sh` overriding the preparer's default
    with `EXPLICIT_IONS=0`. Full-resolution workflows did not apply that
    coarse-only override.
  - CGL z-coordinate mismatch came from post-COM initial conditioning:
    `build_cg_lipid_array` placed each CGL at the geometric COM of its source
    DOPC beads, but the later default leaflet z-separation step shifted upper
    and lower leaflet CGL positions apart. Default z-separation conditioning is
    now disabled by setting its default target to `0.0`; explicit
    `UPSIDE_CG_LIPID_MIN_LEAFLET_Z_SEP` overrides still work.
  - Targeted verification showed explicit-ion preparation adds salts for both
    reported systems: 1RKL `48` NA / `45` CL and 1AFO `47` NA / `51` CL. Direct
    CGL collapse matches source DOPC COM z with max absolute difference
    `1.4e-14 A` for 1RKL and `2.8e-14 A` for 1AFO.
- 2026-06-20: Retained CGL output geometry and speed.
  - Current retained-output measurement confirms the user observation:
    1RKL CGL lipid top/bottom p95/p05 extend `18.15 A/10.75 A` beyond the
    protein p95/p05, compared with `6.93 A/7.11 A` in full resolution. 1AFO
    CGL extends `12.17 A/7.78 A`, compared with `3.26 A/4.21 A` full.
  - VTF rendering was a major visualization amplifier: each physical CGL was
    expanded into a synthetic head/tail display rod with median span
    `~25.6 A`. This is not force-field physics. User correction: do not hide
    the rods; the rods should remain visible and should use full-resolution
    DOPC bead-derived display offsets. New coarse workflows now keep rod
    display as the default and preserve the per-lipid offsets computed while
    collapsing full-resolution DOPC into CGL.
  - User correction: the request was to check whether CGL runs are slow
    because H5 files are being regenerated, not to disable the missing-H5
    fallback. The retained 1RKL/1AFO CGL logs show no `martini_gen_params.py`
    or `Generating...` messages, all runs report `MARTINI parameter libraries
    found`, and the parameter-H5 mtimes predate the retained checkpoints.
    The over-broad coarse missing-H5 fail-fast guard was removed.
  - CGL runtime slowness is partly kernel overhead, not only preparation.
    Skipping potential-value spline accumulation in derivative-only CGL nodes
    improves a copied 1AFO 1000-step benchmark from `3.6 s` to `2.9 s`
    without changing the logged initial/checkpoint energies, but full
    resolution is still `2.6 s`. Further speedup should target CGL
    orientation/table kernels and pair evaluation count, not physical
    parameter shortcuts.
  - Retained stage-7 timing confirms CGL is slower inside MD: 1RKL CGL
    production is `3889 us/step` versus `3436 us/step` full, and 1AFO CGL is
    `3709 us/step` versus `2488 us/step` full. This is consistent with the
    active one-particle CGL kernels: `cg_lipid_pair`, `cg_lipid_rotamer_sc`,
    and especially `cg_lipid_target` use orientation-resolved spline tables
    despite the lower particle count.
  - Lesson: distinguish physical CGL centers from visualization rods when
    interpreting bilayer top/bottom in VMD, but do not remove rods to avoid
    the issue. Fix rod length/placement from the corresponding full-resolution
    DOPC geometry, and keep physical CGL-center alignment as a separate model
    validation question.
  - Lesson: when the user asks whether a suspected cause is responsible,
    first prove or falsify that cause from logs/files. Do not replace that
    diagnostic with a guard or behavior change unless the evidence shows the
    guard is needed.
- 2026-06-20: Strict isolated-conformer CGL result.
  - Semi-isotropic NPT alone does not fix the single-conformer strict-H5
    Stokes-Einstein break. Under NPT, the old strict table gives stable
    geometry and a target-like point but still drops from `D40=0.241` at
    scale `1.1` to `0.106` at scale `1.2`, with mobility-line `R2=0.52`.
  - Implemented a strict-transferable isolated DOPC bonded-conformer ensemble:
    conformers are sampled by Metropolis moves under the dry-MARTINI DOPC ITP
    bond/angle energy only, at `T=0.8647`, with no bilayer trajectory,
    bilayer PMF, force matching, IBI, thickness target, capping, or
    orientation potential.
  - Promoted `parameters/dryMARTINI/dopc.h5` now records
    `dopc_reference_source=isolated_dopc_itp_bonded_mc_ensemble`,
    `dopc_reference_conformer_count=2`, `correction_layer=none`,
    `force_match_enabled=0`, `ibi_enabled=0`, `fluid_pmf_pair_sample_count=0`,
    and no diagnostic IBI/force-match datasets.
  - CGL-only 288-lipid no-ion NPT validation at converted experimental
    temperature passes at the target point:
    `mass_scale=0.012`, `rotational_tau=0.008`, `memory_taus=0.2,2.0`,
    `couplings=0.30375,0.2205` gives `D40=0.257 um^2/s`, fit `R2=0.999`,
    no leaflet crossings, `min |n_z|=0.961`, and same-leaflet minimum
    `6.20 A`.
  - Stokes-Einstein status is improved but not perfect: local target-side
    scales `0.7,0.75,0.8` give `D40=0.290,0.257,0.240` with line
    `R2=0.984`; broader `0.7-1.2` gives `R2=0.950`, with a noisy
    non-monotonic `0.9/1.0` pair. Replicates are still needed for a final
    robust transport claim.
  - Direct 1RKL and 1AFO coarse hybrid workflows with the new default H5 and
    GLE settings complete through stage 7 with all hybrid interactions active,
    exact `T=0.8647`, no explicit ions, and NPT enabled. Protein proxies are
    stable over the 10k production windows: 1RKL hbond sum `26.6 -> 27.3`
    with max `36.7`, Rg `12.7 -> 13.7`; 1AFO hbond sum `86.0 -> 78.3`,
    Rg `15.8 -> 15.5`.
  - Follow-up geometry audit: the protein-hybrid CGL geometry failure was a
    validator artifact. The old median-z leaflet labels misassigned
    near-midplane/ambiguous CGL centers; same-leaflet minima near `3.1-3.3 A`
    were actually opposite-leaflet pairs separated by `~29 A` in z and with
    opposite orientation signs. The corrected orientation-sign leaflet gate
    passes current 1RKL and 1AFO: no orientation sign flips, same-leaflet
    minimum spacing `>=6.25 A`, leaflet separation `~30.1 A`, and
    `min |n_z| >= 0.941`. Remaining risk is replicate robustness, not a known
    physical leaflet failure in these outputs.
  - Independent protein-hybrid replicate validation now passes. 1RKL rep2
    keeps hbond sum `27.70 -> 31.05` within `25.57-35.37`, Rg
    `12.68 -> 13.12 A`, exact `T=0.8647`, no orientation-sign leaflet
    crossings, same-leaflet spacing `6.28 A`, leaflet separation `30.07 A`,
    and `min |n_z|=0.962`. 1AFO rep2 keeps hbond sum `84.86 -> 82.35`
    within `74.28-88.43`, Rg `15.82 -> 15.37 A`, exact `T=0.8647`, no
    orientation-sign leaflet crossings, same-leaflet spacing `6.29 A`,
    leaflet separation `30.29 A`, and `min |n_z|=0.956`.
  - Current conclusion: the strict isolated-conformer one-particle CGL model
    is physical under the stated constraints and passes CGL-only plus
    replicated 1RKL/1AFO hybrid stability gates. The remaining limitation is
    transport statistics, not a known protein-bilayer failure: local
    CGL-only Stokes-Einstein behavior near the target point is acceptable,
    while broader/replicate diffusion-line validation is still recommended.
  - H5 method audit for the installed artifacts: CGL-CGL, CGL-SC,
    CGL-particle, and SC-particle all record the same base generation
    contract, `explicit_dry_martini_constituent_projection`, using
    `_compute_pair_energy_and_gradient`, dry-MARTINI LJ/Coulomb terms,
    `1.2 nm` cutoff, `1e-6 nm` distance guard, native dry-MARTINI to Upside
    unit metadata, no correction layer, and `Tavg=25.0`. The source/target
    ensembles differ by family as expected. The representation is not
    identical: SC-particle stores `direct_rotamer_free_energy_kj_mol` on a
    full radial/angular grid, while CGL-CGL, CGL-SC, and CGL-particle store
    log1p-reduced CGL spline/table representations.
  - Corrected plot interpretation: diffusion versus raw GLE coupling scale is
    not the Stokes-Einstein plot. The SE check should use diffusion versus
    relative mobility, `1 / coupling_scale`, for fixed-shape GLE scans. The
    final isolated-conformer data are locally acceptable near the selected
    point (`0.70,0.75,0.80` scales give mobility-fit `R2=0.979`), but the
    broad scan is not a robust SE validation (`R2=0.947` and the `0.90/1.00`
    points are non-monotonic). Do not claim broad Stokes-Einstein linearity
    from the current single-seed scan.
- 2026-06-20: Strict CGL force-field provenance correction.
  - User constraint: production CGL force fields must not contain bilayer
    target information. The bilayer must emerge from a transferable
    dry-MARTINI molecular force-field projection, not from fitting CGL
    parameters to bilayer thickness, leaflet separation, or protein-bilayer
    alignment.
  - Consequence: prior PMF, force-matched, and IBI CGL-CGL tables are
    diagnostic-only even if they pass short stability and clock gates. Their
    metadata (`correction_layer`, `force_match_enabled`, `ibi_enabled`, or
    nonzero `fluid_pmf_pair_sample_count`) must be rejected by production
    injection.
  - Additional provenance issue: a direct CGL projection can still inherit
    bilayer information if its default reference coordinates are copied from
    a bilayer PDB template. Production DOPC CGL tables should instead use an
    isolated DOPC reference generated from ITP bonded geometry unless an
    explicitly transferable isolated-conformer ensemble is provided.
  - Lesson: do not describe bilayer-observable table training as physical
    production FF under this requirement. It is a diagnostic optimizer, not a
    transferable force field.
- 2026-06-20: Strict-transferable CGL-only clock check.
  - Rebuilt `parameters/dryMARTINI/dopc.h5` from isolated
    `isolated_dopc_itp_bond_angle_tree` DOPC geometry with direct
    dry-MARTINI constituent projection. The installed CGL-CGL metadata has
    `correction_layer=none`, `force_match_enabled=0`, `ibi_enabled=0`,
    `fluid_pmf_pair_sample_count=0`, and no force-match/IBI diagnostic
    datasets.
  - Production validation passes the clean installed H5 at
    `UPSIDE_MARTINI_TEMPERED_AVERAGE_TEMP_UPSIDE=25.0` and rejects a retained
    IBI diagnostic H5 with the new bilayer-trained metadata error.
  - 72-lipid strict-H5 checks are stable but noisy near the target. The
    larger 288-lipid ion-free bilayer gives a clean physical drag line at
    converted experimental temperature `T=0.8647`: GLE coupling scales
    `1.0, 1.5, 1.8` give `D40=0.946, 0.254, 0.130 um^2/s`, all stable, with
    mobility-line `R2=0.998`.
  - Current strict-production CGL clock candidate:
    `mass_scale=0.012`, `rotational_tau=0.008`,
    `memory_taus=0.2,2.0`, `couplings=0.6075,0.441`. This matches the
    `~0.25 um^2/s` 40-ps-step target in the 288-lipid CGL-only bilayer
    without bilayer-derived FF information.
  - Short no-burn-in protein smokes do not yet prove hybrid stability. 1RKL
    keeps hbond stable (`30.98 -> 32.58`) and CGL orientations upright, but
    the old median-plane leaflet gate flags six near-midplane CGL COM
    relabelings. 1AFO keeps compact log Rg and upright CGL orientations, but
    hbond drops `92.74 -> 73.14` over 10k steps. Before changing parameters,
    rerun with less abbreviated/default equilibration; if the 1AFO loss
    persists, inspect transferable CGL-SC/CGL-target projection and handoff.
  - Longer 100k 288-lipid bilayer-only validation changes the conclusion:
    the earlier `1.5x` GLE point cages at long time (`D40=0.089 um^2/s`),
    so it is not a correct-clock point. Lower drag remains geometrically
    stable with upright CGL orientations and no leaflet crossings. Long-run
    `D40` values for coupling scales `0.8, 1.0, 1.1, 1.2, 1.3, 1.5` are
    `0.300, 0.298, 0.278, 0.099, 0.083, 0.089 um^2/s`.
  - The best strict-H5 long-run timescale point is currently
    `memory_taus=0.2,2.0`, `couplings=0.4455,0.3234`, `mass_scale=0.012`,
    `rotational_tau=0.008`, giving `D40=0.278 um^2/s`, stable geometry,
    `min |n_z|=0.939`, and no leaflet crossings. It is only provisional:
    the long-run Stokes-Einstein line fails (`R2=0.64` over all six points)
    because diffusion drops sharply between coupling scales `1.1` and `1.2`.
- 2026-06-20: Retained CGL outputs still fail bilayer/protein alignment.
  - The user-observed VMD issue is real in retained outputs. In
    `outputs/martini_1rkl_hybrid`, CGL display endpoints extend about
    `+16.6 A/-12.6 A` beyond the protein 95th/5th z extents, while
    full-scale 1RKL is only about `+0.8 A/-1.5 A`.
  - Retained 1RKL CGL also shows a secondary-structure proxy loss:
    `output/hbond` total drops `27.75 -> 21.93`; retained full-scale 1RKL
    stays `34.15 -> 33.18`. This supports the visual concern that the
    protein is starting to destabilize in the CGL bilayer-gap region.
  - Current promoted CGL parameters improve short 1RKL hbond stability
    (`30.47 -> 32.19`) but do not fully fix bilayer/protein z alignment.
    Full-DOPC mapped lipid COM 95th percentile is `11.35 A`, while current
    CGL center 95th percentile is `18.01 A`.
  - A second sampled-bin IBI update can pull CGL centers closer to the
    full-DOPC COM target (`center p95=13.94 A`, display top excess `+5.8 A`)
    but destabilizes protein in a 10k smoke (`hbond 30.39 -> 22.33`,
    rising protein potential, kinetic ratio `1.600`). Do not promote IBI-2.
  - Timing finding: retained and current CGL workflows are slower in
    wall-clock because stage preparation/injection repeats expensive CGL
    table/node work for stage 6 and stage 7. The current 1RKL CGL smoke spent
    about `190 s` in two prepare/inject phases versus `24.5 s` in 10k
    stage-7 MD.
  - Rule: do not repair the alignment by changing only VTF display offsets.
    The physical issue is CGL-center axial distribution relative to mapped
    full-DOPC COMs and protein-interface stability. The next table-training
    objective must include both bilayer COM/thickness and protein hbond
    stability gates.
- 2026-06-20: Physical CGL bilayer-gap fix candidate.
  - A sampled-bin IBI update to only the CGL-CGL pair table can remove the
    stage-7 protein-bilayer leaflet gap without adding runtime restraints or
    orientation potentials. The source candidate
    `example/16.MARTINI/outputs/dopc_tavg25_1rkl_full_target_vs_cgl_thickness_ibi1.h5`
    has been promoted to `parameters/dryMARTINI/dopc.h5`.
  - The IBI candidate uses the existing physical policy: sampled
    full-DOPC/model distribution ratios update sampled CGL-CGL bins, while
    unsampled bins and the direct dry-MARTINI overlap core are preserved.
    CGL-SC, CGL-target, SC-particle, and generic MARTINI tables remain active
    from the base H5.
  - The candidate table changes transport; the old GLE couplings
    `0.22,0.16` become too fast (`D40=1.426 um^2/s`). Increasing only the
    FDT GLE drag to `couplings=0.528,0.384` with
    `memory_taus=0.2,2.0`, `mass_scale=0.012`, and
    `rotational_tau=0.008` gives `D40=0.2539 um^2/s` at the `40 ps/step`
    clock.
  - The Stokes-Einstein check is acceptable on this fixed table over a wide
    physical drag range: coupling scales `1.0,1.6,2.4,3.0` give
    `D40=1.426,0.516,0.254,0.200 um^2/s`, monotonic in mobility with
    linear-fit `R2=0.997`.
  - Default stage-7 protein-bilayer gates pass with this promoted table:
    1RKL final center separation/tail gap `30.13 A/-2.83 A`, 1AFO
    `29.96 A/-2.68 A`, both with final CGL orientation
    `|n_z| min/p05/mean` around `0.94/0.966/0.988` and stable protein log
    proxies. Direct no-override 1RKL smoke also completed through stage 7 and
    injected GLE couplings `[0.528,0.384]`.
- 2026-06-20: Direct coarse hybrid wrapper failure.
  - `example/16.MARTINI/1rkl.out` and `example/16.MARTINI/1afo.out` failed
    before MD because the direct wrappers did not export the accepted
    cleaned-H5 runtime defaults; the validator therefore compared installed
    `Tavg=25.0` CGL metadata against fallback `Tavg=10.0`.
  - The H5 files were not stale. Direct coarse hybrid runs must default to
    `UPSIDE_MARTINI_TEMPERED_AVERAGE_TEMP_UPSIDE=25.0`, no explicit ions,
    `TEMPERATURE=0.8647`, `CG_LIPID_MASS_SCALE=0.012`,
    `CG_LIPID_ROTATIONAL_THERMOSTAT_TIMESCALE=0.008`, and the selected
    two-mode FDT GLE (`taus=0.2,2.0`, `couplings=0.22,0.16`,
    `replace_markovian=1`) unless the caller intentionally overrides them.
  - Follow-on correction: full-resolution wrappers failed immediately because
    macOS bash with `set -u` treats `"${empty_array[@]}"` as an unbound
    variable. Coarse-only wrapper arguments must be appended conditionally to
    a non-empty command array, and both coarse and full branches must be
    smoke-tested after wrapper changes.
- 2026-06-20: CGL bilayer gap in retained protein-bilayer outputs.
  - `outputs/martini_1rkl_hybrid` has a real bilayer-thickness problem, not
    only a VMD rendering artifact. In stage 7, CGL center leaflet separation
    is `~42.4 A` and the displayed hydrophobic tail endpoint gap is
    `~10.2 A`.
  - The gap is not caused solely by the 40k stage-7 burn-in: a no-burn-in
    10k stage-7 control still expands from `~33.1 A` to `~42.5 A`.
  - `outputs/martini_1afo_hybrid` shows the same global geometry
    (`~41.4 A` center separation, `~10.3 A` displayed hydrophobic gap), so
    this is not a 1RKL-local protein artifact.
  - Bead-resolved `parameters/dryMARTINI/DOPC.pdb` has PO4-PO4 separation
    `~36.5 A` and tail means nearly meeting. The current CGL protein-bilayer
    outputs must therefore fail a bilayer-thickness/tail-contact gate even
    though their CGL orientation vectors remain mostly normal.
  - Rule: do not hide this by extending VTF display rods or adding a
    thickness restraint. The physical fix should retrain/rebuild the CGL-CGL
    conservative table against mapped full-DOPC cross-leaflet/thickness
    observables and then rerun CGL-only, 1RKL, and 1AFO gates.
- 2026-06-19: Cleaned H5 promotion and CGL-only retest.
  - Rebuilt a protein-compatible DOPC H5 at `Tavg=25.0`, but the resulting
    CGL-CGL numeric table was not the accepted clock table. It was stable in
    the 288-lipid ion-free CGL-only gate but too slow: production-only
    `D40=0.151 um^2/s` after 25% discard.
  - Installed final `parameters/dryMARTINI/dopc.h5` is a merged clean H5:
    accepted CGL-only CGL-CGL/CGL-target tables plus rebuilt current
    CGL-SC/protein-interface tables. Installed
    `parameters/dryMARTINI/sidechain.h5` is rebuilt with the same
    `Tavg=25.0` setting.
  - Metadata now agrees across CGL-CGL, CGL-SC, CGL-particle, and
    SC-particle tables for the shared projection family, unit conversion,
    cutoff, numerical guard, and tempered averaging temperature. The installed
    CGL-CGL table is bitwise identical to
    `example/16.MARTINI/outputs/dopc_npt100k_forcematch_radial_tavg25_cglonly.h5`.
  - Cleaned-H5 CGL-only retest at converted experimental temperature
    `T=0.8647`, `mass_scale=0.012`, `memory_taus=0.2,2.0`,
    `couplings=0.22,0.16`, and `rotational_tau=0.008` passed: full-run
    `D40=0.262 um^2/s`, production-only `D40=0.246 um^2/s` after 25%
    discard, stable geometry, no leaflet crossing, and orientation gate pass.
- 2026-06-19: Cleaned-H5 retained protein validation.
  - Runtime validation with the cleaned `Tavg=25.0` H5 requires exporting
    `UPSIDE_MARTINI_TEMPERED_AVERAGE_TEMP_UPSIDE=25.0`; leaving the runtime
    default at `10.0` correctly rejects the cleaned H5 during schema
    validation.
  - 1RKL CGL/no-ion 20k production passes with all hybrid interface nodes
    active and no production `restraint_position`: finite potential
    `-1759.6..-1629.5`, exact `T=0.8647`, Rg `11.72..12.59 A`, final
    H-bond sum `30.84`, final CGL `min abs(n_z)=0.871`, no CGL below
    `0.7`, no orientation flips, and final same-leaflet spacing `7.54 A`.
  - 1AFO CGL/no-ion 20k production also passes: finite potential
    `-1155.5..-1056.9`, exact `T=0.8647`, Rg `15.37..15.85 A`, final
    H-bond sum `76.87`, final CGL `min abs(n_z)=0.944`, no CGL below
    `0.7`, no orientation flips, and final same-leaflet spacing `7.44 A`.
  - 1RKL full-resolution/no-ion 20k production passes as a full-DOPC
    control: finite potential `-22237.4..-21858.5`, exact `T=0.8647`, Rg
    `12.53..13.02 A`, final H-bond sum `34.24`, active
    `martini_potential` and `martini_sc_table_1body`, and no production
    restraint.
  - 1AFO full-resolution/no-ion fails during the extended full-DOPC handoff
    before a valid production run: stage-7 burn-in loses all H-bonds, Rg
    reaches `3.0e6 A` with final logged Rg `1.54e5 A`, and total potential is
    `~1e16..1e19`. This is a full-resolution packing/barostat/handoff
    blocker, not a CGL H5 failure.
- 2026-06-19: Physical/H5/stability verification audit.
  - Active one-particle CGL production paths have no twist coordinate, no
    standalone orientation potential/restraint, and no active runtime
    force/energy cap. The dynamic CGL orientation is a rotational state with
    FDT-consistent damping/noise; stale capped CGL table metadata is rejected
    during injection.
  - SC/env and BB/env interactions are not disabled: generic pairs that would
    double-count are replaced by `cg_lipid_rotamer_sc`, `cg_lipid_target`,
    `martini_potential`, `martini_sc_table_1body`, and BB proxy projection
    nodes. The remaining caveat is short hybrid startup hold/ramp controls
    for SC/backbone and PO4 handoff; these are not CGL orientation
    potentials, but they are a production-startup control.
  - Current H5 generation is method-consistent by contract for CGL-CGL,
    CGL-SC, CGL-particle, and SC-particle: explicit dry-MARTINI constituent
    bead/rotamer projections with tempered Boltzmann averaging and runtime
    spline tables. The strict consistency mismatch is parameterization, not
    algorithm: the accepted CGL-only timescale H5 uses tempered averaging
    temperature `25.0`, while the installed/full protein-compatible H5 and
    SC-particle tables use `10.0`.
  - Retained HDF5 evidence proves the 288-lipid CGL-only bilayer and selected
    1RKL CGL 20k checkpoint. 1RKL full-resolution, 1AFO full-resolution, and
    1AFO CGL currently have transcript evidence only; their referenced
    checkpoint/VTF outputs are missing locally, so independent bilayer/CGL
    orientation and secondary-structure rechecks require regeneration.
- 2026-06-18: No-ion one-particle CGL GLE retune outcome.
  - All tests used converted experimental temperature `T=0.8647`, one CGL
    particle per DOPC, unchanged conservative spline tables, no explicit ions,
    and `mass_scale=0.012` so each lipid CGL is approximately one Upside mass
    unit.
  - `coupling=0.248`, `memory_tau=1.0`, `rotational_tau=0.008` is the best
    CGL-only timescale point found: 100k CGL-only `D40=0.252 um^2/s`, fit
    `R2=0.99975`, final `min abs(n_z)=0.911`, and stable bilayer geometry.
    It is not promotable because the no-ion 1RKL smoke final frame has one
    CGL below `abs(n_z)=0.7`.
  - `coupling=0.250`, `rotational_tau=0.008` is the best current 1RKL
    stability point: 20k no-ion 1RKL keeps finite energy, exact temperature,
    protein Rg `12.63..13.27 A`, final H-bond sum `27.44`, no leaflet
    crossing, final `min abs(n_z)=0.721`, and same-leaflet 3D spacing above
    `7 A`. It is not a correct-clock point because CGL-only `D40=0.209
    um^2/s`, implying about `48 ps/step` against the `0.25 um^2/s` target.
  - Stronger rotational damping (`rotational_tau=0.006`) did not repair the
    physical transport gate. Single-seed 100k CGL-only
    `coupling=0.236,0.240,0.242,0.244` gave
    `D40=0.211,0.228,0.305,0.289`; the target-window response is not
    Stokes-Einstein monotonic in `1/(coupling^2*memory_tau)`.
  - Current conclusion: mass scaling plus FDT-consistent translational GLE and
    rotational Langevin damping alone has not produced a one-particle CGL that
    simultaneously matches the 40 ps/step clock, passes 1RKL bilayer stability,
    and passes the Stokes-Einstein linearity check.
- 2026-06-18: Multi-exponential one-particle CGL GLE.
  - Implemented H5 schema `cgl_exponential_memory_gle_v2` as a finite sum of
    exponential memory kernels with positive `memory_tau[]` and `coupling[]`,
    plus one auxiliary momentum per CGL per mode. The runtime update remains
    FDT-consistent, preserves one CGL particle per lipid, and does not change
    conservative spline interactions or add orientation potentials.
  - Smoke validation confirmed `input/cgl_gle/aux_momentum` shape
    `(n_mode,n_cgl,3)`, output shape `(frames,n_mode,n_cgl,3)`, and exact
    restart promotion for all modes.
  - Best tested kernel:
    `memory_taus=0.2,2.0`, `couplings=0.22,0.16`,
    `rotational_tau=0.008`, `mass_scale=0.012`, no explicit ions. CGL-only
    100k gives `D40=0.234 um^2/s`, which is closer to the target than the
    previous stable single-mode point (`0.209`) but still slow.
  - The same kernel passes a 20k no-ion 1RKL stability gate: final
    `min abs(n_z)=0.731`, no CGL below `0.7`, no leaflet crossing, final
    same-leaflet 3D spacing `6.90..7.00 A`, exact temperature, and compact
    protein metrics.
  - It still fails the Stokes-Einstein requirement. Fixed-shape coupling
    scaling `0.90..1.10` produced `D40=0.267,0.297,0.234,0.241,0.194`
    versus mobility coordinate `1/sum(c_i^2 tau_i)`, with line `R2=0.65`.
    Therefore multi-exponential GLE is not yet a fully physical clock match.
- 2026-06-18: Reference-fit GLE from full-DOPC COM dynamics.
  - The full-resolution DOPC NPT reference trajectory has COM diffusion
    `3.277 um^2/s` in physical time and `13.108 um^2/s` after the standard
    Martini `x4` time correction, so the bead-resolved dry-MARTINI reference
    is not the source of the CGL clock mismatch.
  - Fitting the mapped lateral COM VACF with `memory_taus=0.2,2.0` produced a
    good normalized VACF fit (`R2=0.983`) but collapsed to a native-clock
    short-memory kernel: `couplings=2.595,2.5e-7`.
  - Using the fitted short-memory shape with target-scaled couplings did not
    recover the required CGL transport. CGL-only 100k checks at
    `coupling=0.45,0.55,0.65`, `rotational_tau=0.008`, and
    `mass_scale=0.012` were geometrically stable but gave
    `D40=0.189,0.141,0.199`, with no Stokes-Einstein linearity.
  - Conclusion: within the current one-particle CGL representation, physical
    transport-model changes alone have been exhausted. The next acceptable
    physical change should alter the trained conservative one-particle CGL
    bilayer target/model class, not add capping, force scaling, hidden
    particles, DPD, or an orientation restraint.
- 2026-06-19: Pair-only IBI correction on the radial force-matched CGL-CGL
  table.
  - Added a caging/order diagnostic for mapped full-DOPC versus one-particle
    CGL trajectories. The long `c0p30` CGL model has a structural first-shell
    mismatch: coordination is higher than target by `0.55`, pair probability
    is excessive near `7.875 A`, and first-shell orientation-dot is higher by
    `0.111`.
  - Built a physical sampled-bin IBI candidate on top of the radial
    force-matched full H5:
    `example/16.MARTINI/outputs/dopc_npt100k_forcematch_radial_ibi_c0p30long_tavg25_full.h5`.
    Only `cg_lipid_table/cg_lipid_pair` is replaced. All CGL-target,
    CGL-SC, CGL-particle, and SC-particle metadata are copied from the full
    base H5. The direct overlap core and unsampled bins are retained.
  - The IBI update touched `1944/5880` sampled non-core tensor bins with
    `deltaU=-17.541..16.583 kJ/mol`.
  - The new table is stable but not promotable. Short 20k runs are mobile and
    target-like at `scale=2.55` (`D40=0.241`), but 100k runs still cage:
    `scale=1.6 D40=0.226`, `scale=2.0 D40=0.177`,
    `scale=2.55 D40=0.141`. The best full-run point has block diffusion
    `0.771,0.191,0.113,0.096`, so apparent target diffusion is dominated by
    an initial transient.
  - Conclusion: pair-distribution IBI on the current one-particle pair table
    does not remove long-time caging. The next physical conservative change
    needs a model class that directly represents bilayer packing cooperativity
    or area/undulation response, not another pair-only sampled-bin correction.
- 2026-06-19: Stokes-Einstein recheck of the accepted two-mode GLE candidate.
  - User clarified that `D40=0.234 um^2/s` is acceptable against the
    `0.25 um^2/s` clock target, so the active question is the transport
    relation rather than exact target matching.
  - Repeating the common fixed-shape coupling scales with new seeds kept the
    CGL-only bilayer geometrically stable but produced much slower diffusion:
    `scale=0.95 D40=0.085`, `scale=1.00 D40=0.127`,
    `scale=1.05 D40=0.062 um^2/s`.
  - Combining old and new seeds at the common scales gives seed means
    `0.191`, `0.180`, and `0.151 um^2/s` for scales `0.95`, `1.00`, and
    `1.05`, with a monotonic mean line versus
    `1/sum(c_i^2 tau_i)` (`R2=0.905`). However, the seed standard deviations
    are `0.150`, `0.076`, and `0.126 um^2/s`, too large to validate a
    robust Stokes-Einstein relation.
  - Block diffusion shows the failure mode is real trajectory aging/slowing,
    not only a fitting artifact. Old near-target runs stayed mobile across
    blocks, e.g. old `scale=1.00` blocks
    `0.238,0.166,0.319,0.213`. New repeats commonly start mobile then slow:
    new `scale=1.00` blocks `0.335,0.115,0.059,0.081`; new `scale=1.05`
    blocks `0.178,0.111,0.040,0.015`.
  - Conclusion: the original two-mode point is acceptable on absolute clock
    for a favorable seed and remains stable, but it is not yet physically
    validated because diffusion is not seed-robust and the Stokes-Einstein
    gate is dominated by slow-state/caging variability.
- 2026-06-19: 288-lipid ion-free CGL-only Stokes-Einstein validation.
  - Added a calibration-runner `--cgl-bilayer-pdb` input and a reproducible
    2x2 DOPC PDB tiler. The generated template
    `example/16.MARTINI/outputs/dopc_2x2_ionfree.pdb` contains 288 DOPC
    lipids, 4032 atoms, box `100.184 x 100.184 x 85.000 A`, and no template
    ions.
  - Added production-only block analysis by allowing
    `analyze_cgl_diffusion_blocks.py` to discard an initial trajectory
    fraction before the MSD fit and block split.
  - The full protein-compatible force-matched H5 is rejected by the current
    strict CGL-target validator before MD. For the strict ion-free CGL-only
    bilayer gate there are no target particles, so the correct artifact is
    the CGL-only CGL-CGL H5:
    `example/16.MARTINI/outputs/dopc_npt100k_forcematch_radial_tavg25_cglonly.h5`.
    A full protein-interface H5 still needs current CGL-target metadata before
    final protein validation.
  - Smoke test at `memory_taus=0.2,2.0`, `couplings=0.22,0.16`,
    `mass_scale=0.012`, `rotational_tau=0.008`, and `T=0.8647` passed with
    288 CGL particles only, stable geometry, and short-run `D40=0.224`.
  - Six 100k production runs over scales `0.95,1.00,1.05` with two seeds per
    scale were all geometrically stable and target-like after discarding the
    first 25% of frames. Production-only `D40` means:
    `0.95 -> 0.241 +/- 0.019`, `1.00 -> 0.265 +/- 0.026`,
    `1.05 -> 0.255 +/- 0.003 um^2/s`.
  - The narrow target window is too flat/noisy for a standalone SE validation:
    production-only means over `0.95,1.00,1.05` fit versus
    `1/sum(c_i^2 tau_i)` with negative slope and `R2=0.413`.
  - Wider physical endpoints show the expected friction response without
    geometry failure: `0.80 -> D40=0.327`, `1.20 -> D40=0.175` after 25%
    discard. The production-only wide subset `0.80,1.00,1.20` gives a
    positive SE fit with `R2=0.926`; all five scale means give a positive fit
    with `R2=0.832`. Full-run means give `R2=0.996` for the wide subset and
    `R2=0.947` across all five means.
  - Conclusion: the severe SE failure in the 72-lipid runs was largely a
    finite-size/template-ion/initial-state artifact. In the 288-lipid ion-free
    bilayer the two-mode GLE is physically friction-responsive over a wider
    mobility span and has acceptable target-window diffusion. The practical
    parameter set remains `couplings=0.22,0.16` because it has prior no-ion
    1RKL stability evidence and gives production-only 288-lipid
    `D40=0.265 +/- 0.026 um^2/s`.

## Lessons / User Corrections
- User correction: `D40=0.234 um^2/s` for the original two-mode CGL GLE is
  acceptable against the `0.25 um^2/s` target. Do not over-optimize the
  absolute diffusion value once it is in that window. The active decision
  criterion is whether the physical transport model passes the
  Stokes-Einstein relation robustly.

- User correction: keep exactly one CGL particle per lipid. Do not recommend
  or start a multi-site CGL representation unless the user explicitly reopens
  that option. Within the one-particle constraint, physical next attempts
  must be conservative model-class changes such as trained density/area terms,
  not hidden particles, orientation restraints, capping, or friction-only
  tuning.

- 2026-06-18: One-particle conservative density correction tests.
  - Added a table-driven `cg_lipid_density` runtime term that preserves exactly
    one CGL particle per lipid. It evaluates `sum_i F(rho_i)` where
    `rho_i=sum_j w(r_ij)` over CGL centers. It applies center forces only,
    has no orientation derivative, no CGLD marker, no cap, no force scaling,
    and no additional orientation potential.
  - Built
    `example/16.MARTINI/outputs/dopc_npt100k_forcematch_radial_densityre_tavg25_cglonly.h5`
    from the radial force-matched pair table plus a relative-entropy
    local-density embedding. The density kernel cutoff was derived from the
    full-DOPC same-leaflet CGL-center RDF first minimum (`9.875 A`), and the
    embedding used `F=kBT ln(P_model(rho)/P_target(rho))` from the full-DOPC
    target and long caged CGL model.
  - Relative-entropy density was geometrically stable but worsened transport.
    Short NVT `tau=14,17,20` gave `D40=0.101,0.175,0.171` with
    Stokes-Einstein `R2=0.712`; `tau=20,23,26` gave
    `D40=0.174,0.147,0.146` with negative slope and `R2=0.765`.
    Distribution inspection showed the full-DOPC target had a broader
    high-density tail than the caged CGL model, so the valid
    relative-entropy update introduced attractive high-density wells and
    increased local caging.
  - Built
    `example/16.MARTINI/outputs/dopc_npt100k_forcematch_radial_densityharmonic_tavg25_cglonly.h5`
    using the same density kernel but a target-fluctuation harmonic embedding:
    `F=0.5*(kBT/var_target_rho)*(rho-mean_target_rho)^2`, with
    `mean rho=0.250`, `var rho=0.03596`, and
    `k=70.1 kJ/mol/rho^2` from the full-DOPC distribution.
  - Harmonic density was also stable but failed the transport relation:
    `tau=14,17,20` gave `D40=0.242,0.166,0.200`, with negative fitted
    Stokes-Einstein slope and `R2=0.308`. The target-like `tau=14` point is
    not promotable because the required linearity check fails.
  - Current conclusion: within one-particle Markovian Langevin CGL,
    conservative local-density many-body corrections can stabilize geometry
    but still do not produce the requested physical clock plus
    Stokes-Einstein behavior. The next physical one-particle option should be
    an FDT-consistent generalized Langevin / memory-kernel implicit-solvent
    model fitted from mapped full-DOPC COM dynamics, not more conservative
    fitting layers.
  - Implemented the one-particle CGL GLE option as an FDT-consistent
    exponential memory kernel with an auxiliary momentum Markovian embedding.
    It is activated only through `/input/cgl_gle`, preserves exactly one CGL
    particle per lipid, leaves the conservative CGL spline tables unchanged,
    and can replace the ordinary Markovian CGL thermostat for translational
    transport.
  - CGL-only GLE scans used the radial force-matched conservative table
    `example/16.MARTINI/outputs/dopc_npt100k_forcematch_radial_tavg25_cglonly.h5`
    at production temperature `T=0.8647`. The useful physical mass setting is
    `mass_scale=0.012`, close to one Upside mass unit for each 14-bead DOPC
    CGL. `mass_scale=0.02` stayed stable but was too slow and nonlinear.
  - Long 100k GLE transport points with `mass_scale=0.012`,
    `memory_tau=1.0`, and `rotational_tau=0.03` were geometrically stable:
    `c=0.65,0.50,0.40,0.28` gave `D40=0.140,0.134,0.193,0.277` in single
    long runs, with MSD-fit `R2>=0.993`. Adding a second `c=0.28` seed gave
    `D40=0.403`, so the weak-coupling edge is fast and seed-sensitive.
  - The 100k transport line is now qualitatively Stokes-Einstein consistent
    when using the correct GLE coordinate `1/(coupling^2*memory_tau)`:
    the long-point mean line for `c=0.65,0.50,0.40,0.28` gives `R2=0.968`.
    However, the production coupling is not yet seed-robust. Direct 100k
    target-window tests gave `c=0.295 -> D40=0.225` and
    `c=0.303 -> D40=0.230`, while `c=0.28` gave `0.277/0.403`. Current
    recommendation is to treat `coupling=0.29..0.30` as the practical window
    and run more seeds before protein-bilayer validation.
  - Follow-up target-window robustness selected `coupling=0.30` as the
    conservative stable CGL-only candidate. Two 100k seeds at `c=0.29`
    were reproducibly fast (`D40=0.298,0.289`), but the 500k check failed
    same-leaflet geometry (`same_leaflet_min_xy=0.184 A`), so `c=0.29` is
    rejected. Two 100k seeds at `c=0.30` gave `D40=0.192,0.245`, and the
    500k check stayed stable. In the 500k stable run, the short-lag linear
    regime gives `D40=0.249..0.252` for max lag 100 native-time, while the
    full long-lag fit to 500 native-time gives `D40=0.206`. Five 200-time
    blocks remain mobile and target-like (`D40=0.238,0.254,0.277,0.231,0.323`),
    so there is no long-time caging signal. Use `c=0.30` for the first
    protein-bilayer stability test, but report both short-lag and long-lag
    diffusion estimates.
  - Protein-bilayer workflow integration requires a full H5, not the
    CGL-only screening H5. The CGL-only radial force-matched table lacks the
    CGL-target/protein-interface metadata needed by `inject_cg_lipid_nodes`.
    The full compatible file built for 1RKL is
    `example/16.MARTINI/outputs/dopc_npt100k_forcematch_radial_tavg25_full.h5`,
    with CGL-CGL, CGL-SC, CGL-particle, and SC/protein-interface tables using
    the same force-matched CGL-CGL target and table-generation contract.
  - `py/martini_prepare_system.py` now injects CGL transport metadata into
    every coarse hybrid stage, including mass scaling, replacement of the
    ordinary Markovian CGL thermostat, and `/input/cgl_gle` restart state.
    GLE auxiliary momentum must be part of restart state; finalization now
    promotes `/output/cgl_gle_aux_momentum[-1]` into
    `/input/cgl_gle/aux_momentum`.
  - CGL dynamic orientation is also restart state. A continuation must promote
    `/output/cgl_orientation[-1]` into
    `/input/potential/cgl_orientation_state/direction` and
    `/output/cgl_orientation_mom[-1]` into `/input/cgl_orientation_mom`.
    Without this, a production restart combines final coordinates with stale
    orientation vectors and can produce a false instability. The bad 100k
    continuation had initial potential `14940`; after restart promotion the
    same source starts at `-2968.68`.
  - First 1RKL validation with the selected candidate completed a smoke run
    and a 20k-step continuation with all hybrid interactions active. The
    20k continuation stayed finite at `T=0.8647`, protein Rg stayed
    `12.65..13.23 A`, and the bilayer did not exchange leaflets under a
    two-cluster z assignment. The minimum same-monolayer CGL 3D spacing was
    `6.35 A`, final same-monolayer 3D minimum was `6.89 A`, and leaflet
    center separation was `32.9 A` at the final frame.
  - Lesson: do not use the median z value alone to classify leaflets when the
    leaflet counts are unequal or the bilayer/protein placement shifts the z
    distribution. It can mark opposite-leaflet particles as same leaflet and
    create false XY-overlap alarms. Use sign-based or two-cluster z leaflet
    assignment and report full 3D separation plus z separation for close XY
    pairs.
  - Restart-fixed 100k 1RKL continuation with explicit ions is numerically and
    protein-stable but not yet bilayer-orientation stable. It keeps finite
    energy, exact `T=0.8647`, protein Rg `12.10..13.99 A`, final H-bond sum
    `27.3`, no leaflet exchange, and minimum same-monolayer CGL 3D spacing
    `6.03 A`; however final CGL orientation has `min abs(n_z)=0.00035`,
    mean `0.659`, and `120/275` lipids below the CGL-only `0.70` threshold.
  - The tilted CGLs are strongly associated with explicit ions in the
    implicit-water hybrid setup: at the final restart-fixed 100k frame,
    `15/15` lipids with `abs(n_z)<0.1` and `109/120` lipids with
    `abs(n_z)<0.7` are nearest to ION targets. This points to mobile naked
    ions as the current physical-model problem. The correct next test is an
    implicit-salt/no-explicit-ion preparation, not friction retuning, table
    capping, or an added orientation potential.
  - `EXPLICIT_IONS=0` now exists for hybrid preparation. It is not validated
    yet: first smoke attempts generated larger non-comparable bilayers
    (`773` and `565` CGL lipids) and were interrupted in the CGL conditioning
    loop before MD. A controlled no-ion packing path with comparable bilayer
    size/geometry remains the next step.

- 2026-06-17: Dynamic one-particle vector CGL long-time gate.
  - Dynamic CGL orientation is implemented without fixed vectors or CGLD
    atoms. `compose_vector6d` consumes `pos` and `cgl_orientation_state`, and
    `/output/cgl_orientation` stays unit-normalized.
  - Short 20k scans can overestimate lateral diffusion because the bilayer has
    an initially mobile relaxation period. Long 100k CGL-only validation is
    now required before any protein-bilayer test.
  - Best short-run dynamic candidate tested: `Tavg=25`, production
    `T=0.8647`, `mass_scale=0.02`, translational tau `9`, rotational tau `1`,
    native transverse inertia. It gave 4/4 stable 20k seeds with mean
    `D40=0.244 +/- 0.074 um^2/s`, but a 100k run gave only
    `D40=0.034 um^2/s`.
  - Higher CGL table tempering did not fix the long-time cage: `Tavg=30`
    remained slow and failed orientation; `Tavg=50` produced severe
    orientation failure and unphysical diffusion artifacts.
  - Current conclusion: mass scaling, translational Langevin tau, rotational
    Langevin tau, and simple higher-tempered H5 tables are insufficient to
    match the Upside-core clock while preserving long-time CGL bilayer
    stability. The next physical change should be conservative-table
    retraining/rebuilding for the dynamic vector CGL against a fluid bilayer
    target, not more thermostat tuning.
  - Correction: do not promote an ensemble CGL table built with capped or
    subsampled conformer-pair combinations unless convergence is explicitly
    demonstrated. The table builder should average all conformer pairs for the
    selected physical conformer ensemble.
  - Implemented a clean ensemble-table path that samples DOPC conformers from
    a full-resolution Upside trajectory and averages all conformer pairs for
    the selected ensemble. The temporary conformer-pair cap was removed before
    validation.
  - The first ensemble table derived scalar orientation properties from the
    mean conformer geometry, which underestimated transverse inertia. The
    table builder now derives scalar properties per conformer and averages
    those physical scalars.
  - Corrected table:
    `example/16.MARTINI/outputs/dopc_ensemble2_allpairs_tavg25_scalaravg.h5`.
    It uses two full-DOPC trajectory conformers, all four conformer-pair
    combinations, `Tavg=25`, runtime Boltzmann temperature `0.8647`,
    transverse inertia `9381.587 g mol^-1 A^2`, orientation length
    `12.496 A`, and orientation mass `60.084 g/mol`.
  - Short CGL-only scans with the corrected table and rotational tau `0.1`
    found a stable narrow window at translational `tau=5.5,6.0,6.5`, but it
    is not promotable: `D40=0.173,0.176,0.212 um^2/s`, below the
    `0.25 um^2/s` target, and Stokes-Einstein `R2=0.815`. Higher tau values
    fail the orientation gate. This remains a conservative table/training
    definition problem, not something to fix by adding an orientation
    potential or by pushing rotational friction further.
  - Improved trajectory conformer sampling to select evenly across the
    frame-lipid pool after relaxation instead of tying conformer index to one
    frame and one lipid index. CGL-CGL bead-energy sampling was vectorized;
    a scalar/vectorized check matched to `1.9e-6 kJ/mol` sorted max
    difference.
  - Built CGL-only screening table
    `example/16.MARTINI/outputs/dopc_ensemble4_allpairs_tavg25_framelipid_cglonly.h5`
    from four full-DOPC trajectory conformers, all 16 conformer-pair
    combinations, `Tavg=25`, runtime Boltzmann temperature `0.8647`, no
    conservative cap attributes, no CGLD, and no SC-CGL table. This H5 is for
    bilayer screening only; protein-bilayer validation still requires a full
    SC-CGL table.
  - Four-conformer short screens: rotational tau `0.1` failed geometry at
    translational `tau=4.0,4.5,5.0` and `5.5,6.0,6.5` due one nearly
    membrane-parallel lipid. Stronger physical rotational damping
    (`rot_tau=0.03`, FDT-consistent) passed at `tau=4.0,4.5,5.0` with
    `D40=0.243,0.398,0.702 um^2/s` and Stokes-Einstein `R2=0.966`.
  - Long 100k validation rejects the four-conformer screening table:
    `tau=4.0`, `rot_tau=0.03` remained geometrically stable but slowed to
    `D40=0.054 um^2/s`. Four equal-time blocks gave
    `D40=0.212,0.041,0.045,0.085`, confirming long-time caging rather than a
    short-fit artifact.
  - Current conclusion remains: direct pair-energy CGL tables, even with a
    broader trajectory conformer ensemble and physical rotational damping,
    do not maintain target long-time diffusion. The next physical table
    change should be PMF/relative-entropy/force-matching against an
    equilibrated fluid-bilayer CGL pair distribution, not isolated two-lipid
    direct-energy averaging.
  - Implemented a fluid-bilayer PMF CGL-CGL table option. The final tested
    table,
    `example/16.MARTINI/outputs/dopc_fluidpmf4_oriented_contactcore_tavg25_cglonly.h5`,
    uses a full-DOPC orientation-resolved PMF in runtime coordinates
    `(r, -n1 dot n12, n2 dot n12)`, reverse-pair symmetrization, direct
    dry-MARTINI bead-energy overlap core below the DOPC-derived CGL contact
    distance, `Tavg=25`, production `T=0.8647`, no cap attributes, no CGLD,
    and no SC-CGL table for CGL-only screening.
  - Radial-only PMF was rejected: it produced target-like short diffusion but
    caused leaflet crossing because it removed orientation/leaflet
    information. The orientation-resolved PMF fixed geometry and kept all
    tested high-tau points stable.
  - Orientation-resolved PMF short scans: `tau=8,11,14` passed geometry with
    `D40=0.137,0.188,0.357` and Stokes-Einstein `R2=0.913`. A higher
    `tau=20,30,40` scan was stable and `tau=20` was target-like
    (`D40=0.236`), but the three-point line check failed due mobility
    saturation at high tau.
  - Long 100k validation rejects the current PMF table. `tau=14` gave
    `D40=0.033`; `tau=20` gave `D40=0.097`, with block values
    `0.281,0.065,0.021,0.012`. Geometry remained stable, so the blocker is
    long-time caging, not bilayer breakup or orientation failure.
  - Current conclusion: the PMF idea is physically cleaner and stabilizes
    bilayer geometry, but the available short, small, NVT full-DOPC reference
    is not a sufficient fluid-bilayer target for long-time CGL mobility. The
    next physical requirement is a better equilibrated full-DOPC reference
    trajectory, preferably larger and semi-isotropic/NPT at the target
    temperature/area, before another PMF or relative-entropy table rebuild.
  - Added an NPT-capable full-DOPC reference runner using the existing
    semi-isotropic barostat metadata. The 20k full-DOPC NPT reference ran at
    `T=0.8647` (`303.15 K`), `1 bar`, area per lipid `69.70 A^2`, finite
    energy/box, and no NC3 leaflet sign changes.
  - Built
    `example/16.MARTINI/outputs/dopc_npt20k_fluidpmf4_oriented_contactcore_tavg25_cglonly.h5`
    from that NPT reference: four conformers, all 16 conformer-pair direct
    core combinations, 1,027,512 ordered fluid-PMF pair samples, 1,569
    occupied tensor bins, `Tavg=25`, production `T=0.8647`, no cap attrs,
    no CGLD, and no SC-CGL table for CGL-only screening.
  - The NPT-derived PMF table passed a low-tau short Stokes-Einstein screen:
    `tau=8,11,14` gave `D40=0.132,0.160,0.198` with `R2=0.993` and all
    geometry gates passing. The target-window short scan at `tau=17,20,23`
    was geometrically stable and `tau=20` was target-like (`D40=0.247`), but
    the three-point high-tau line check failed (`R2=0.376`).
  - Long 100k validation rejects the 20k NPT-derived PMF table. `tau=20`
    gave full-run `D40=0.109`; block `D40` values were
    `0.279,0.129,0.049,0.116`, and 19 CGL orientations became near-parallel
    to the membrane. This is still not promotable to 1RKL. The next physical
    step is a longer and/or larger full-DOPC reference before PMF rebuild,
    not additional CGL thermostat tuning or orientation potentials.
  - Extended the full-DOPC NPT reference to 100k at `T=0.8647` (`303.15 K`)
    and 1 bar. It remained stable with area per lipid `69.70 A^2`, finite
    energy/box, no NC3 leaflet sign changes, and full-DOPC lateral diffusion
    `3.28 um^2/s` in physical units (`13.11 um^2/s` with the usual MARTINI
    `x4` acceleration).
  - Built
    `example/16.MARTINI/outputs/dopc_npt100k_fluidpmf4_oriented_contactcore_tavg25_cglonly.h5`
    from the 100k NPT reference. PMF coverage improved to 2,561,112 ordered
    pair samples, 121 occupied radial bins, and 2,986 occupied tensor bins.
    The table kept the direct dry-MARTINI overlap core, had `Tavg=25`,
    production `T=0.8647`, no cap attrs, no CGLD, and no SC-CGL table for
    CGL-only screening.
  - The 100k-NPT PMF table improved short-run orientation geometry but failed
    the Stokes-Einstein requirement. In NVT, `tau=8,11,14` gave
    `D40=0.137,0.161,0.138` with `R2≈0`; `tau=17,20,23` gave
    `D40=0.247,0.212,0.358` with `R2=0.532`. Enabling the same
    semi-isotropic NPT ensemble for CGL kept geometry stable but still failed
    linearity: `tau=17,20,23` gave `D40=0.195,0.119,0.315`, `R2=0.371`.
  - Current conclusion: one-shot pair PMF from the small 72-lipid reference,
    even with longer NPT sampling and CGL NPT screening, does not produce a
    physical linear transport regime. The next physical change should be an
    iterative relative-entropy/IBI or force-matching refinement of the
    dynamic-vector CGL-CGL table against full-DOPC CGL-center/orientation
    distributions, with pressure/area included in the target. Do not proceed
    by thermostat tuning, capping, orientation potentials, or 1RKL testing.
  - Implemented center/orientation pair-distribution extraction for both
    full-DOPC trajectories and dynamic one-particle CGL trajectories, then
    added sampled-bin relative-entropy/IBI table updates. The update rule is
    `U_new = U_old + kBT ln(P_model/P_target)`; unsampled bins retain the
    existing conservative table value, and the direct overlap core is
    preserved. No caps, force scaling, CGLD marker, fixed orientation, or
    additional orientation potential were introduced.
  - Built
    `example/16.MARTINI/outputs/dopc_npt100k_fluidpmf4_ibi1_tau17model_tavg25_cglonly.h5`.
    Metadata: target/model samples `2,561,112/1,027,512`, sampled tensor bins
    `1058/7007`, correction range `-15.079..9.360 kJ/mol`, production
    `T=0.8647`, `Tavg=25`, no cap attrs, and no SC-CGL table for screening.
  - IBI-1 short NVT gates improved: `tau=8,11,14` gave
    `D40=0.071,0.107,0.184` with `R2=0.958`; `tau=14,15.5,17` gave
    `D40=0.137,0.187,0.243` with `R2=0.999`. Geometry passed throughout.
    However, the 100k `tau=17` gate failed badly: full-run `D40=0.014`;
    four 50-native-time blocks gave `D40=0.165,0.049,0.030,0.019`.
    IBI-1 is not promotable to 1RKL.
  - Added a previous-table base-energy path so IBI-2 is truly iterative:
    `U_2 = U_1 + kBT ln(P_1/P_target)`, using the IBI-1 raw CGL-CGL energy
    grid as the conservative base. Built
    `example/16.MARTINI/outputs/dopc_npt100k_fluidpmf4_ibi2_tau17longmodel_tavg25_cglonly.h5`
    from the caged 100k IBI-1 model. Metadata: target/model samples
    `2,561,112/3,839,112`, sampled bins `1149/7007`, correction range
    `-19.244..9.339 kJ/mol`, production `T=0.8647`, `Tavg=25`, no cap attrs.
  - IBI-2 stayed geometrically stable but did not solve transport:
    `tau=14,15.5,17` gave `D40=0.086,0.137,0.174` with `R2=0.993`, and
    `tau=20,23,26` gave `D40=0.126,0.100,0.136` with `R2=0.064`.
    The current evidence says pairwise one-particle vector CGL distribution
    matching is insufficient for the requested clock unless the target/model
    representation is changed more fundamentally.
  - H5 generation consistency update: SC-particle, CGL-particle, CGL-SC, and
    CGL-CGL now share a documented base-table contract in generated metadata:
    explicit dry-MARTINI constituent LJ/Coulomb projection over the relevant
    physical rotamer/conformer/orientation ensemble, with the same unit
    contract, nonbonded cutoff, and numerical distance guard recorded. The
    SC-particle path now uses the same `_compute_pair_energy_and_gradient`
    helper as the CGL paths instead of its own ad hoc LJ/Coulomb loop.
  - Important limitation: CGL-CGL fluid-PMF/IBI remains a recorded correction
    layer on top of the common base table, not the universal method for all
    interactions. Applying those corrections to CGL-SC, CGL-particle, or
    SC-particle would require a corresponding physical protein/interface
    target ensemble; doing it without that target would be an unphysical
    parameter twist.
  - Implemented first CGL-CGL force-matching target:
    radial generalized force projection from the 100k full-DOPC NPT
    reference. The full-DOPC H5 does not store forces, so the implementation
    recomputes bead-level dry-MARTINI pair forces from trajectory coordinates,
    projects each lipid-lipid pair onto
    `(r, -n1 dot n12, n2 dot n12)`, integrates `dU/dr` along radial tensor
    bins, preserves the direct overlap core, and leaves unsampled bins at the
    existing physical base value.
  - Built
    `example/16.MARTINI/outputs/dopc_npt100k_forcematch_radial_tavg25_cglonly.h5`.
    Metadata: 100 full-DOPC frames, 511,200 ordered pair-force samples,
    2,482/7,007 tensor bins updated, production `T=0.8647`, `Tavg=25`, no
    cap attrs, no SC-CGL table for screening.
  - Radial force matching improved the CGL-only bilayer relative to IBI but
    did not pass the hard timescale gate. Short `tau=14,17,20` was stable and
    line-like (`D40=0.118,0.187,0.219`, `R2=0.958`), while `tau=20,23,26`
    reached target-like short diffusion but failed linearity
    (`D40=0.291,0.295,0.359`, `R2=0.796`). Long 100k validations still slowed:
    `tau=20` full `D40=0.213`, block `0.288,0.198,0.160,0.151`; `tau=26`
    full `D40=0.195`, block `0.247,0.232,0.193,0.176`.
  - Current conclusion: radial projected force matching is a substantial
    improvement but still leaves long-time slowing. The next physical
    conservative-table change should include angular generalized
    forces/torques in the CGL-CGL fit, or move to a more expressive physical
    CGL orientation representation. Do not try to force this table through by
    changing only Langevin friction.
  - Implemented generalized CGL-CGL force matching by projecting recomputed
    dry-MARTINI bead forces and torques onto `dU/dr`,
    `dU/d(-n1 dot n12)`, and `dU/d(n2 dot n12)`, then least-squares fitting
    a conservative tensor energy grid. The implementation preserves the
    direct overlap core, retains unsampled bins, adds no runtime potential,
    and uses no capping or force scaling.
  - Built
    `example/16.MARTINI/outputs/dopc_npt100k_forcematch_generalized_tavg25_cglonly.h5`.
    Metadata: 100 full-DOPC frames, 511,200 ordered pair-force samples,
    2,516/7,007 tensor bins updated, 6,043 LS equations, LS residual RMS
    `18.99 kJ/mol`, production `T=0.8647`, `Tavg=25`, no cap attrs, no
    SC-CGL table for screening.
  - First generalized short screens reject this table before long validation:
    `tau=14,17,20` is geometrically stable but nonlinear
    (`D40=0.168,0.211,0.185`, `R2=0.148`); `tau=20,23,26` is also stable
    but nonlinear (`D40=0.196,0.195,0.317`, `R2=0.744`). Because individual
    MSD fits are clean and geometry passes, the blocker is the transport
    relation, not bilayer breakup. The next physical check is a
    stricter-sampling generalized rebuild with more reference frames and
    higher per-bin count threshold before moving to a different
    representation.
  - Built stricter generalized table
    `example/16.MARTINI/outputs/dopc_npt100k_forcematch_generalized_frames300_mincount12_tavg25_cglonly.h5`
    with 300 full-DOPC frames, 1,533,600 ordered force samples, 2,510 updated
    tensor bins, min-count 12, production `T=0.8647`, and no cap attrs. The
    generalized LS residual worsened to `25.37 kJ/mol`, indicating the
    pairwise one-particle vector coordinates cannot cleanly represent the
    projected full-DOPC force/torque field.
  - The stricter generalized table failed short CGL-only validation:
    `tau=14` and `tau=20` failed orientation geometry; only `tau=17` passed
    geometry and it was too slow (`D40=0.125`). Current conclusion:
    generalized pair force matching is exhausted for the one-particle vector
    CGL representation. The next physical direction should be a more
    expressive lipid representation, preferably explicit physical multi-site
    CGL sites such as head/body-tail or head/two-tail sites with spline
    tables trained consistently from dry-MARTINI/full-DOPC targets. A
    one-particle many-body density/area-compressibility conservative term is
    a secondary option only if it is trained to physical bilayer targets and
    not used as a diffusion-tuning knob.
- 2026-06-17: Vector-only CGL reimplementation after user removed the need
  for a CGLD particle; fixed orientation later rejected by user.
  - Coarse CGL setup now writes one particle per DOPC lipid. It no longer
    appends CGLD orientation sites, no longer creates CGL-CGLD bonds, no
    longer registers CGLD nonbonded exclusions, and no longer writes
    `compose_vector6d/orientation_index`.
  - Runtime `compose_vector6d` already supports the no-`orientation_index`
    path: it uses the stored `direction` vector and propagates translational
    derivatives to the CGL center only. The generated DOPC-only CGL bilayer has
    88 atoms: 72 CGL plus 8 NA and 8 CL.
  - The scan used the existing `Tavg=25` CGL table and production temperature
    `T=0.8647`, so the bilayer was simulated at the experimental converted
    temperature. No conservative spline interactions were flattened, capped,
    or disabled.
  - Smoke H5 validation confirmed no `CGLD` atom types, no
    `compose_vector6d/orientation_index`, no orientation bonds, and zero
    CGLD thermostatting.
  - Stokes-Einstein checks pass for vector-only CGL at `mass_scale=0.02`.
    `tau=2.5,3.5,5.0` gave stable `D40=0.0586,0.0923,0.1397 um^2/s` with
    `R2=0.9997`. `tau=7.0,8.5,10.0` gave stable
    `D40=0.2277,0.3613,0.6105 um^2/s` with `R2=0.9705`.
  - Best fixed-vector diagnostic candidate for the `56x` target is `mass_scale=0.02`,
    `tau=6.0`: 4/4 stable 20k seeds with `D40=0.239 +/- 0.054 um^2/s`.
    `tau=7.0` was also stable in 5/5 seeds but faster/noisier:
    `0.308 +/- 0.100 um^2/s`.
  - User correction: CGL orientation must not be fixed. These fixed-vector
    scans are only a diagnostic baseline proving that removing CGLD removes
    marker instability; they are not a promotable physical CGL model.
  - Lesson: when the requested model is a vector particle, confirm whether the
    vector is a dynamic state before treating stored direction metadata as an
    acceptable implementation.
- 2026-06-17: Physical high-tempering continuation after user disallowed
  twisting/capping/orientation potentials.
  - Built separate copied DOPC CGL H5 files at table averaging temperatures
    `Tavg=25.0` and `Tavg=30.0` under `example/16.MARTINI/outputs/`.
    Metadata verified that production Boltzmann temperature stayed
    `0.8647`; only `azimuthal_average_temperature_upside` changed.
  - Rejected CGLD constraint-style marker fixes. Fixed-center projection,
    mass-weighted constraints, and removing independent CGLD bath coupling
    produced unphysical diffusion or orientation tumbling. The rejected runtime
    and CLI paths were removed after testing.
  - `Tavg=30`, `cgl_mass_scale=0.02`, native CGLD mass gave target-like
    diffusion but failed the orientation geometry gate for all
    `tau=18,20,22` single-seed points.
  - `Tavg=25`, `cgl_mass_scale=0.02`, native CGLD mass, shared tau passed the
    broad single-seed scan: `tau=18,20,22` had
    `D40=0.240,0.286,0.391 um^2/s`, all geometry-stable, with
    Stokes-Einstein `R2=0.950`.
  - Seed replicates remain the blocker. `Tavg=25`, `mass=0.02`, `tau=18.0`
    gave 3/4 stable seeds with all-seed mean `D40=0.250 +/- 0.038` and
    stable-seed mean `0.238 +/- 0.035`. `tau=18.5` gave 3/4 stable seeds but
    was too fast/noisy: stable-seed mean `0.301 +/- 0.092`.
  - Refining mass/friction found isolated target-like stable points, e.g.
    `mass=0.015,tau=16` (`D40=0.257` in one scan), but did not produce a
    three-point stable Stokes-Einstein window. Do not promote a protein-bilayer
    validation candidate yet.
  - A marker-projection internal-orientation implementation was tested and
    rejected, then removed. At `Tavg=25`, `CGL mass=0.02`, `cgl_tau=18`,
    `cgld_tau=18,5,2,1,0.5,0.2,0.1`, all short scans failed orientation
    geometry and many produced unphysical `D40` values from several to hundreds
    of `um^2/s`.
  - Increasing only CGLD inertia is stabilizing but too slow. With
    `Tavg=25`, `CGL mass=0.02`, `tau=18`, `CGLD mass scale=5` gave stable
    `D40=0.0677`; scale `10` gave stable `0.0344`; scale `2` still failed
    orientation. Raising CGL tau for `CGLD mass scale=5` reached only
    `D40=0.130` at stable `tau=30`; `tau=35` failed geometry. Lowering CGL
    mass to `0.01` with `CGLD mass scale=5` remained too slow and failed at
    `tau=30`.
  - Current conclusion: the hidden-marker representation is the blocker. The
    next physical implementation must be a real rotational state variable for
    CGL orientation, not a projected atom or another conservative orientation
    parameter.
- 2026-06-17: CGL `56x` target scan after rejecting DPD.
  - DPD was rejected because this MARTINI hybrid model uses implicit water.
    Do not add pairwise DPD transport to this workflow.
  - With a `56x` coarse-graining factor against the estimated pure-DOPC
    reference, the working CGL target under the Upside-core `40 ps/step` clock
    is `14 / 56 = 0.25 um^2/s`.
  - Added CGL/CGLD decoupled thermostat support to
    `example/16.MARTINI/cgl_timescale_calibration.py`.
  - Rebuilt a separate DOPC CGL table at tempered average temperature
    `20.0`: `example/16.MARTINI/outputs/dopc_tempered_20.h5`. Production
    scans still used the experimental converted temperature `T=0.8647`.
  - Best broad scan with `cgl_mass_scale=0.02`, native CGLD mass, shared tau,
    and `Tavg=20`: stable `tau=18,20,22` gave `D40=0.204,0.301,0.447` and
    passed the Stokes-Einstein line check (`R2=0.987`).
  - Seed replicates near the target are promising but not production-ready:
    `tau=18.5` had 3/4 stable 20k seeds with stable-seed mean
    `D40=0.259 +/- 0.024 um^2/s`; `tau=19` had 3/4 stable seeds with mean
    `0.293 +/- 0.036`; `tau=18` had 3/4 stable seeds with mean `0.247` but
    high scatter.
  - Decoupling CGLD tau to `15` stabilized one scan but broke the
    Stokes-Einstein trend. Scaling the orientation carrier bond by `4x` also
    did not help; lower tau points failed geometry and the stable points were
    insufficient for linearity.
  - Current recommendation: do not promote a protein-bilayer parameter yet.
    The physical H5-tempering + mass-scaling path makes the target plausible,
    but the remaining CGLD orientation-marker failures need a cleaner rigid
    orientation treatment or longer robustness testing before protein
    validation.
- 2026-06-17: Full-resolution dry-MARTINI DOPC mismatch audit.
  - Added a full-resolution DOPC COM diffusion analyzer and runner:
    `example/16.MARTINI/analyze_dopc_diffusion.py` and
    `example/16.MARTINI/full_dopc_diffusion_reference.py`.
  - Fixed `inject_particles_table` for full-resolution MARTINI: the runtime
    expects `coefficient_indices` to reference the local `coefficients` rows,
    so injected spline grids must be ordered in that same local row order.
  - The canonical `parameters/dryMARTINI/DOPC.pdb` is not directly stable for
    full-resolution production Upside MD. It has hard-core inter-lipid
    contacts that give initial production energy around `6.7e12 E_up` and
    immediate coordinate blow-up. A GROMACS steepest-descent minimization using
    the dry-MARTINI topology was needed to produce
    `example/16.MARTINI/outputs/full_dopc_gmx_relax/em.upside.pdb`.
  - Stable full-DOPC 20k-step Upside runs at `T=0.8647` (`303.15 K`),
    `time_step=0.002`, `inner_step=1` gave:
    `tau=5 D40=0.004993 um^2/s`, `tau=10 D40=0.005190 um^2/s`,
    `tau=20 D40=0.004871 um^2/s`.
  - The physical Upside time unit derived from project units is
    `0.20289664298287868 ps`; therefore `time_step=0.002` is
    `0.00040579328596575736 ps`, not a Martini-like `40 ps` step.
  - The same trajectories convert to raw physical diffusion
    `4.92, 5.12, 4.80 um^2/s` for `tau=5,10,20`; applying the common Martini
    `4x` effective-time convention gives `19.69, 20.47, 19.21 um^2/s`.
    This is reasonably close to the Dry Martini SI DOPC-equivalent value
    `24.5 um^2/s`, given the tested run is at `303.15 K` while the SI row is
    at `310 K`, and the local trajectory covers only about `8.12 ps` raw
    physical time.
  - The `tau=5` full-DOPC run was geometrically stable: final nonbonded
    minimum distance `4.07 A`, no pairs below `3 A`, area/lipid `69.70 A^2`,
    PO4-PO4 thickness `39.78 A`.
  - Stokes-Einstein check failed for full-DOPC in this range:
    `D40` versus `tau` slope `-1.15e-5`, `R^2=0.298`.
  - Source audit found no mass-aware Langevin/integrator bug. The `mv`
    integrator updates coordinates with `dt / mass` when MARTINI masses are
    registered, and the OU thermostat scales random momentum kicks with
    `sqrt(mass)`.
  - A same-topology GROMACS NVT diagnostic from the minimized 72-lipid setup
    was also slow and poorly linear over 1 ns. This points at the local small,
    minimized, NVT setup/analysis horizon as a poor paper reproduction, not an
    Upside-integrator-only explanation.
  - Current interpretation: the apparent paper mismatch was mainly an
    accounting error from using the empirical Upside-core `40 ps/step` mapping
    as the physical dry-MARTINI clock. Keep those clocks separate. A rigorous
    paper-level comparison still requires a larger equilibrated DOPC bilayer,
    semi-isotropic/NPT conditions, `310 K`, and ns-us production analysis.
- 2026-06-17: Stokes-Einstein diagnostic fix.
  - The previous diffusion analyzers used a single-time-origin MSD
    (`r(t)-r(0)` only). That is too noisy for the short 72-lipid trajectories
    used here and can make a three-point `D` versus `tau` test fail even when
    the individual MSD fits look linear.
  - Updated both CGL and full-DOPC analyzers to use multi-time-origin lateral
    MSD up to `max_lag_fraction=0.5`, reporting the lag-time fit window and
    sample counts. The CGL calibration runner now passes this setting through.
  - Reanalyzing the existing full-DOPC `tau=5,10,20` trajectories with the
    multi-origin estimator changes the trend from negative to positive, but it
    still fails the three-point line check because `tau=20` is slightly below
    `tau=10` on this short trajectory (`R2=0.476`). Treat this as a poor
    diagnostic window, not as proof of an integrator defect.
  - Native-mass CGL at `tau=2.5,5,10`, `T=0.8647`, `time_step=0.002`, 20k
    steps passes the Stokes-Einstein line-fit gate:
    `D_native=0.02307,0.04006,0.08175 A^2/native_time`, `R2=0.9977`, all
    three rows geometry-stable. Artifacts:
    `example/16.MARTINI/outputs/cgl_stokes_fix_native_tau2p5_5_10/`.
  - Matched full-DOPC from the GROMACS-minimized PDB at `tau=2.5,5,10`,
    `T=0.8647`, `time_step=0.002`, 20k steps also passes the current
    line-fit gate: `D_native=0.005835,0.009000,0.011219 A^2/native_time`,
    `R2=0.9175`. Artifacts:
    `example/16.MARTINI/outputs/full_dopc_stokes_fix_tau2p5_5_10/`.
  - Rule for future calibration scans: choose tau windows inside the stable,
    monotonic friction-controlled regime. Do not include high-tau points that
    fail geometry or enter saturated/noisy dense-bilayer mobility unless the
    purpose is explicitly to show breakdown of Stokes-Einstein scaling.
- 2026-06-17: User correction on physical CGL timescale matching.
  - Do not twist CGL potentials, flatten mobility barriers, weaken spline
    tables, or otherwise change physical interactions solely to make CGL match
    the Upside-core `40 ps/step` empirical clock.
  - Any CGL table rebuild must be justified by a physical inconsistency in the
    current coarse-graining input, ensemble, temperature, structural target, or
    transport model. A diffusion target alone is not enough.
  - Valid outcomes include concluding that the single-particle CGL
    representation has a separate physical clock, or adopting a physically
    motivated transport model such as momentum-conserving thermostatting, but
    not arbitrary parameter tuning.
  - User clarified that mass scaling is acceptable because the Upside core uses
    unit masses. Exact core-like CGL mass checks were therefore run:
    `CGL=1.0 m_up` with native `CGLD`, and `CGL=1.0 m_up` with
    `CGLD=1.0 m_up`.
  - Exact unit-mass results still do not reach the `40 ps/step` target. Best
    stable row is `CGL=1.0 m_up`, `CGLD=1.0 m_up`, `tau=5`, with
    `D40=0.1331 um^2/s`; `tau=10` and `20` fail orientation/packing. This is
    still about `13x` below a very lenient DOPC/8 target, `26x` below DOPC/4,
    and `105x` below `14 um^2/s`.
- 2026-06-17: Reported dry-MARTINI DOPC lateral diffusion target.
  - The ACS Figshare SI for Arnarez et al., Dry Martini
    (`10.1021/ct500477k.s001`) reports the DOPC-equivalent row as `dry`,
    `310 K`, `small`, `PC`, tail entry `DV` with `CCDC/CCDC` tails. Its
    lateral diffusion is `2.45e-2` in units of `1e-5 cm^2/s`, i.e.
    `24.5 um^2/s` with `0.5 um^2/s` reported error.
  - The matching wet Martini row is `5.89e-2 x 1e-5 cm^2/s`, i.e.
    `58.9 um^2/s`, so dry Martini is about `2.4x` slower than wet Martini for
    this DOPC-equivalent bilayer in the SI table.
  - This creates a real target-definition issue: comparing CGL directly to
    experimental DOPC values around `8-14 um^2/s` is not the same as matching
    the dry-MARTINI model's own reported dynamics. If the goal is to preserve
    dry-MARTINI kinetics under the Upside runtime, the calibration reference
    should be the dry-MARTINI DOPC SI value, with explicit handling of Martini
    effective-time scaling.
  - Collapsing 14-bead DOPC into CGL can justify a separate effective-time
    factor, but the measured stable CGL scan would require a very large one:
    `24.5 / 0.046216 = 530x` to make the best stable row match raw
    dry-MARTINI DOPC, or `14 / 0.046216 = 303x` to match the pure-DOPC
    experimental reference. Ordinary Martini-style factors near `4-10x` do not
    make the current stable CGL diffusion plausible at `40 ps/step`.
- 2026-06-16: Higher CGL table tempering test.
  - The current `dopc.h5` has two temperature concepts. Runtime
    `boltzmann_temperature_upside=0.8647` reconstructs physical spline
    energies from the log1p-reduced table. Rebuild-time
    `azimuthal_average_temperature_upside=10.0` controls the tempered PMF
    averaging used to generate the table.
  - Changing `boltzmann_temperature_upside` in an already-built input is not a
    valid test of higher PMF tempering; it rescales stored log1p-reduced energy
    deviations. A meaningful test requires rebuilding `dopc.h5`.
  - Added `UPSIDE_MARTINI_TEMPERED_AVERAGE_TEMP_UPSIDE` as a table-generation
    and validation setting. The default remains `10.0`, preserving current
    behavior. Added `--dopc-h5` and `--tempered-average-temp` to
    `example/16.MARTINI/cgl_timescale_calibration.py` for copied-table scans.
  - Built copied high-temper tables:
    `example/16.MARTINI/outputs/cgl_tempered_tables/temp20/dopc.h5` and
    `example/16.MARTINI/outputs/cgl_tempered_tables/temp50/dopc.h5`.
  - Native-mass `T_avg=20` scan (`tau=5,10,20,50,100`, 20k steps): stable
    rows were only `tau=5` and `tau=10`, with `D40=0.019804` and
    `0.022046 um^2/s`. Higher tau rows failed geometry.
  - Native-mass `T_avg=50` scan: stable rows were again only `tau=5` and
    `tau=10`, with `D40=0.013194` and `0.019805 um^2/s`. Higher tau rows
    failed geometry.
  - Combined `T_avg=50` with `mass_scale=0.5,0.3,0.2,0.1` and
    `tau=5,10,20`. Only `mass_scale=0.5,tau=5` passed geometry
    (`D40=0.019653 um^2/s`). The fastest row,
    `mass_scale=0.1,tau=20`, reached `D40=0.503283 um^2/s` but failed with
    `8` flips, `13` bad-orientation rods, and same-leaflet minimum `0.443 A`.
  - Conclusion: higher tempered table generation at `T_avg=20` or `50` does
    not resolve the timescale/stability problem. It can change diffusion, but
    the rows that move toward the target still fail CGL orientation/geometry.
  - Updated next direction: decouple translational and rotational CGL dynamics
    before further table tempering. Specifically, scan lower CGL center mass
    while keeping CGLD orientation mass and/or orientation damping/stiffness
    near native values, then validate whether diffusion can increase without
    rod flips.
- 2026-06-16: CGL mass-scale/dissipation scan result.
  - Ran a 40-point CGL-only bilayer scan under
    `example/16.MARTINI/outputs/cgl_timescale_scan_20k/`:
    `mass_scale=1.0,0.8,0.6,0.5,0.4,0.3,0.2,0.1`, `tau=5,10,20,50,100`,
    `nsteps=20000`, `time_step=0.002`.
  - Stable-gated points remained far below the DOPC target under a `40 ps/step`
    mapping. The best stable point was `mass_scale=0.2,tau=5`, with
    `D40=0.046216 um^2/s`; compared with the `14 um^2/s` pure-DOPC reference,
    this is about `303x` too slow. Native-mass stable points were
    `D40=0.013215` at `tau=5` and `0.030222` at `tau=10`.
  - Rows closest to the target were unstable. The closest row,
    `mass_scale=0.2,tau=100`, reached `D40=7.003 um^2/s` but failed with
    leaflet crossings `1/1`, `14` flips, `34` bad-orientation rods, and low
    diffusion-fit quality `R^2=0.663`.
  - Dominant failure mode across unstable rows was CGL orientation failure:
    `31` rows had bad orientation alignment, `22` had flips, `6` had
    same-leaflet overlaps, and `1` had leaflet crossing.
  - No mass scale had at least three stable tau points, so no candidate passed
    the Stokes-Einstein requirement. Native mass and `mass_scale=0.8` had only
    two stable tau points; all lower mass scales had one or zero.
  - Conclusion: uniform CGL/CGLD mass scaling plus scalar Langevin dissipation
    cannot defensibly match the Upside core `~40 ps/step` timescale while
    preserving bilayer stability. The failure happens before protein-bilayer
    validation, so protein-bilayer checks are not needed to reject this
    parameter family.
  - Next directions that preserve hybrid physical interactions:
    1. Decouple CGL translational and rotational dynamics: lower only the CGL
       center mass while keeping CGLD orientation mass/spring/rotational
       damping stable, then rescan diffusion and orientation.
    2. If dynamics-only decoupling is insufficient, regenerate the CGL-CGL
       spline table with a mobility-calibrated lateral barrier/corrugation
       while preserving area, thickness, orientation, and DOPC structural
       observables.
    3. Consider an implicit-water-compatible transport model instead of
       independent OU Langevin friction. It must not require explicit water
       or disable physical interactions.
- 2026-06-16: Rebuilt CGL-only timescale calibration workflow.
  - Current checkout did not contain the previously documented calibration
    source files; only stale `__pycache__` entries were present. Recreated the
    workflow as `example/16.MARTINI/analyze_lipid_diffusion.py` and
    `example/16.MARTINI/cgl_timescale_calibration.py`.
  - The CGL-only bilayer input can be built from
    `parameters/dryMARTINI/DOPC.pdb` through the existing preparation API, but
    it must also inject `particle.h5` and `dopc.h5`; conversion alone leaves
    runtime nodes without the required spline-table data.
  - Mass scaling is implemented as a multiplier on existing `CGL` and `CGLD`
    masses in copied HDF5 inputs, preserving their relative inertia. CGL
    friction is controlled through `/input/thermostat_timescale`, while
    non-CGL atoms keep the global thermostat timescale.
  - The analyzer computes lateral diffusion from unwrapped `xy` CGL MSD with
    `D_xy = slope(MSD_xy) / 4` and reports the effective MD step implied by
    DOPC references from Scientific Reports 9, 1508 (2019),
    `https://www.nature.com/articles/s41598-018-37814-x`: simulated
    unlabeled DOPC `8.4 um^2/s`, dilute FCS tracer about `12 um^2/s`, and
    estimated pure DOPC about `14 um^2/s`.
  - Stokes-Einstein validation is now explicit: for each fixed mass scale,
    stable-gated points are fit as `D_xy` versus `tau = 1/gamma`. At least
    three stable tau points are required for a meaningful line-fit decision;
    two points are reported but marked `insufficient_stable_points`.
  - A short verification grid under
    `example/16.MARTINI/outputs/cgl_timescale_grid_2k/` ran
    `mass_scale=1.0,0.5,0.1` by `tau=5,20` for `2000` steps. Native mass
    passed bilayer geometry at both tau values. `mass_scale=0.5` and `0.1`
    failed orientation geometry, even though lower mass increased apparent
    diffusion.
  - Current short-grid values under a `40 ps/step` mapping: native mass gives
    `0.0283-0.0412 um^2/s`; failed `mass_scale=0.1` gives
    `0.332-0.503 um^2/s`. The latter cannot be used for calibration because
    it fails the stability gate.
  - No production parameter is recommended from the short grid. A valid
    recommendation requires longer trajectories, at least three stable tau
    points per candidate mass scale, and a supplied protein-bilayer stability
    check through `--protein-input`.
- 2026-06-16: 1RKL MARTINI secondary-structure regression.
  - User-reported bad VTFs in both `martini_1rkl_hybrid` and
    `martini_1rkl_hybrid_full` share the same generated runtime setting:
    `--temperature 1.2`.
  - Root cause is the uncommitted wrapper change in
    `example/16.MARTINI/run_sim_hybrid.sh` that changed the default
    `TEMPERATURE` from committed `0.8647` to `1.2`.
  - Current bad-output metrics: CGL stage-7 hbond sum `23.92 -> 14.28`, CA Rg
    `12.73 -> 14.02 A`; full-lipid stage-7 hbond sum `26.38 -> 15.45`, CA Rg
    `12.69 -> 13.73 A`.
  - Controlled copied-checkpoint test from the same CGL stage-7 start:
    `T=0.8647` gives hbond final `33.48`, last-20 mean `32.15`, CA RMSD max
    `1.03 A`; `T=1.2` gives hbond final `22.30`, last-20 mean `21.55`, CA RMSD
    max `1.34 A`.
  - Rule: do not change wrapper default temperature while debugging MARTINI
    timescale/integrator behavior. Temperature changes affect protein
    secondary-structure stability and must be explicit experimental inputs.
  - Hybrid verbose progress logs should print actual MD step counts. Do not put
    rounded elapsed time in the leading progress field, because short frame
    intervals can repeat the same integer and look like duplicate steps.
- 2026-06-15: User correction on CGL diffusion target.
  - Do not use full-resolution DOPC diffusion as the calibration target for
    this task. The requested target is the CGL diffusion expected under the
    Upside core time mapping of `40 ps` per integration step.
  - The required deliverable is a dissipation sweep: dissipation constant on
    the x-axis, measured CGL diffusion on the y-axis, and a horizontal target
    line for the expected diffusion rate under the 40 ps/step mapping.
  - Rule: when the user asks for a timescale calibration, define the target
    conversion explicitly before comparing against any convenient reference
    trajectory.
- 2026-06-16: User correction on final CGL timescale interpretation.
  - Keep Upside core and dry-MARTINI runtime settings uniform when possible,
    including the same simulation timestep.
  - Hybrid MARTINI workflow time accounting must be explicit. Do not hide the
    old `3x` inner-step multiplier in frame intervals or restarts. The MARTINI
    workflow should not expose a separate inner-step knob; use one internal
    `MARTINI_MD_INNER_STEP=1` constant for frame intervals, restart public
    time, and the `upside --inner-step` argument.
  - Short native-mass checks support that default: explicit `--integrator mv
    --inner-step 1` completed 20k steps for both the CGL-only bilayer and the
    1RKL protein-bilayer checkpoint. Both outputs ended at `time_final=40.0`
    for `nsteps=20000`, `time_step=0.002`; bilayer geometry passed in both,
    and the protein-bilayer run passed protein stability with final CA RMSD
    `1.65 A`.
  - The primary dry-MARTINI timescale estimate should come from measured
    CGL-only bilayer lateral diffusion compared with a real DOPC diffusion
    reference. The Upside `40 ps/step` mapping is a comparison point, not a
    fitted target to force by unstable mass/friction settings.
  - Scientific Reports 9, 1508 (2019), DOI
    `10.1038/s41598-018-37814-x`, reports DOPC diffusion values useful as
    references: simulated unlabeled DOPC `8.4 +/- 0.4 um^2/s`, dilute FCS
    tracer about `12 um^2/s`, and an estimated pure-DOPC value about
    `14 um^2/s` after correcting for dye hydrodynamic drag.
  - Do not use 1RKL/protein-bilayer trajectories to measure the dry-MARTINI
    lipid timescale; protein crowding and protein-lipid perturbation make that
    diffusion value system-specific. Use protein-bilayer runs as stability
    checks only.
  - The observable is lateral diffusion, computed from `xy` MSD with
    `D_xy = slope(MSD_xy) / 4`; do not use 3D diffusion for bilayer lipid
    timescale calibration.
  - Native-mass CGL-only bilayer measurement at `tau=5` and
    `time_step=0.002` gives lateral diffusion `0.017319996 A^2/native_time`,
    or `0.008659998 um^2/s` under the Upside `40 ps/step` mapping. Compared
    with real-DOPC references, this implies effective dry-MARTINI steps of
    `0.041238 ps` for `8.4 um^2/s`, `0.028867 ps` for `12.0 um^2/s`, and
    `0.024743 ps` for `14.0 um^2/s`.
  - Rule: report CGL-only bilayer lateral diffusion in native units and infer
    the dry-MARTINI ps/step explicitly from the chosen real-DOPC reference. Do
    not present low-mass unstable target crossings or protein-bilayer
    diffusion as valid dry-MARTINI lipid timescale measurements.
- 2026-06-15: CGL-specific diffusion calibration path.
  - The runtime OU thermostat previously accepted only one global
    `--thermostat-timescale`; MARTINI mass-aware noise was already present,
    but no per-particle dissipation control existed.
  - Added optional `/input/thermostat_timescale` semantics: absent dataset means
    existing global behavior, present dataset gives one positive OU timescale
    per atom.
  - Coarse hybrid stage generation now writes CGL/CGLD-specific thermostat
    timescales only when `CG_LIPID_THERMOSTAT_TIMESCALE` (or
    `--cg-lipid-thermostat-timescale`) is positive. Non-CGL atoms remain at
    the global `THERMOSTAT_TIMESCALE`.
  - Short 1RKL smoke with `CG_LIPID_THERMOSTAT_TIMESCALE=8.0` completed through
    stage 7.0 and produced `248` atoms at global timescale `5.0` plus `550`
    CGL/CGLD atoms at timescale `8.0`.
  - Added a local lateral-diffusion analyzer. For the corrected target, use
    `40 ps` per actual integration step, not full DOPC comparison. With
    production `time_step=0.002`, `1 native time = 20 ns`.
  - The dissipation sweep over `tau=5,10,20,50,100,200,500` did not reach the
    40 ps/step target candidate of `0.5 um^2/s` (`1.0 A^2/native_time`).
    Measured CGL diffusion stayed around `0.012-0.016 um^2/s`, with the best
    point at `tau=20`, `gamma=0.05`, `D=0.015631 um^2/s`.
  - Current conclusion: no defensible production CGL dissipation constant can
    be selected from dissipation-only tuning in this tested hybrid setup. The
    best point is about `32x` below the target line.
  - Native active coarse 1RKL CGL masses are much heavier than protein carrier
    atoms: `CGL=84 m_up` (`1008 g/mol`), `CGLD=6.42074 m_up` (`77.05 g/mol`),
    while `N/CA/C/O=6 m_up` (`72 g/mol`).
  - A copied-input diagnostic with `CGL/CGLD mass=1.0 m_up` was numerically
    stable for 5000-step stage-7 runs across `tau=5,10,20,50,100,200,500`.
    The best measured point was `tau=10`, `gamma=0.1`, `D=0.101894 um^2/s`,
    which is still about `4.9x` below the `0.5 um^2/s` target line.
  - Working rule: do not infer CGL kinetic timescale from geometric bead size
    alone. Check the effective mass and dense-bilayer caging/friction response
    directly, and keep any mass retuning explicit as dynamics calibration
    rather than a hidden force-field change.
  - User correction: timestep setup must remain uniform across the hybrid
    simulation and consistent with other Upside examples. Keep `--time-step
    0.002` for the whole system during CGL/dry-MARTINI sweeps; do not introduce
    a dry-MARTINI-only timestep or CGL-only subcycling to solve calibration.
  - Multi-mass sweeps with the uniform `time_step=0.002` completed without
    short-run numerical failures for `CGL/CGLD mass=0.5,0.25,0.1,0.05,0.02,
    0.015,0.01,0.005 m_up`.
  - The first direct crossing of the `0.5 um^2/s` target was at
    `mass=0.015 m_up`, `tau=5`, `gamma=0.2`, with measured
    `D=0.500349 um^2/s`. The log-interpolated recommendation for that mass is
    `tau=5.0635`, `gamma=0.19749`.
  - Lighter masses also cross the target but overshoot more strongly in the
    tested grid: `mass=0.01 m_up` recommends `tau=18.64`, `gamma=0.05365`;
    `mass=0.005 m_up` recommends `tau=18.06`, `gamma=0.05536`. Prefer the
    least aggressive mass reduction that crosses the target, pending longer
    validation.
  - Longer validation rejects `mass=0.015 m_up` at the uniform production
    timestep. The 20k-step bilayer-only copied-input run failed with CGL
    leaflet crossings `14/35`, `bad_parallel=26`, `bad_flip=45`, and
    same-leaflet nearest-neighbor minimum `0.040 A`. The native-mass
    bilayer-only control passed the same 20k-step run, so this is a mass
    override problem rather than a generic test-input failure.
  - The 20k-step 1RKL protein-bilayer validation with `mass=0.015 m_up` and
    `tau=5.0635` also failed CGL geometry: leaflet crossings `30/4`,
    `bad_parallel=8`, `bad_flip=2`, and same-leaflet nearest-neighbor minimum
    `0.120 A`. Protein coordinates stayed finite, but CA Rg drifted
    `12.71 -> 14.33 A`, final CA RMSD was `5.52 A`, and hbond last-20 mean
    fell to `24.19`.
  - The `mass=0.015` gamma/diffusion data do not show the expected
    Stokes-Einstein `D proportional 1/gamma` behavior. Linear fit quality was
    `R^2=0.507` for `D` versus `tau=1/gamma`, with a negative slope; `D`
    versus `gamma` had `R^2=0.715`, which is not the expected friction
    scaling. Treat the short-run diffusion crossing as a noisy unstable-regime
    artifact, not a valid kinetic calibration.
  - Added a stability-gated mass/tau workflow that runs standalone bilayer
    validation first, then protein-bilayer validation, then CGL diffusion only
    for points that pass both geometry gates. The workflow keeps the same
    uniform `time_step=0.002` for all atoms and writes `gamma` versus diffusion
    plots with failed points marked separately.
  - Stable-gated sweeps over `mass=0.02,0.03,0.05,0.1,0.2,0.5,1.0 m_up` and
    `tau=5,10,20,50` found no valid calibration point. All `mass<=0.5` rows
    failed standalone bilayer geometry; `mass=1.0,tau=5` passed standalone
    bilayer and protein stability, but the protein-bilayer CGL geometry failed.
    Therefore no production CGL dissipation constant is recommended from the
    tested low-to-moderate mass range.
- 2026-06-15: Cross-system CGL cutoff/performance check.
  - User correction: 1RKL and 1AFO must use the same CGL cutoff rules; do not
    validate or discuss the cutoff change using only one of the two protein
    examples.
  - Active 1RKL and 1AFO coarse outputs use identical CGL runtime cutoff
    metadata: CGL-CGL `41.3 A`, CGL-target `26.6 A`, SC-CGL `26.6 A`,
    bead-level nonbonded cutoff `1.2 nm`, DOPC axis radius `14.55848 A`, and
    DOPC perpendicular radius `4.114344 A`.
  - Current active stage-7 performance satisfies the expected ordering for both
    examples: 1RKL CGL/full `1707.99/3328.29 us/systems/step`; 1AFO CGL/full
    `1923.47/2522.99 us/systems/step`.
- 2026-06-15: Current post-directional-pairlist validation.
  - The old active `martini_1rkl_hybrid` stage-7 output contained real CGL
    orientation defects rather than only a VTF rendering issue: the final frame
    had three `aligned-z < 0.70` rods, and two were already bad at the
    production handoff.
  - Regenerating the default 1RKL coarse hybrid output with the current
    directional CGL pairlist produced clean final geometry: aligned-z
    min/p05/mean `0.770/0.849/0.949`, `bad_parallel=0`, `bad_flip=0`, and no
    leaflet crossings.
  - Current measured performance does not support another cutoff change:
    active 1RKL coarse CGL hybrid is `1707.99 us/systems/step` versus
    full-resolution lipid `3328.29 us/systems/step`; same-size DOPC-only CGL
    bilayer is `350.29 us/systems/step` versus full-resolution lipid
    `797.32 us/systems/step`.
  - Next performance work should start from a fresh profile. SC-CGL should not
    be directionally narrowed unless side-chain extent bounds make the filter
    conservative for all rotamers.
- 2026-06-15: CGL directional pairlist correction.
  - User correction: because `dopc.h5` CGL tables are trained from two DOPC
    molecules, a center-distance cutoff that admits a third DOPC-sized gap is
    not a defensible runtime candidate rule. Pairlist filtering should use the
    represented DOPC bead envelope plus the dry-MARTINI bead nonbonded cutoff;
    the wider two-body radial table support should not drive candidate-pair
    work when the two represented molecules cannot have bead-level contact.
- 2026-06-15: Hybrid cleanup style rule.
  - User correction: cleanup is not only deletion of stale code. The hybrid
    interface should match the human-written master branch style across C++,
    Python, shell, and Markdown. Avoid generated-sounding comments, debugging
    chronology, patch-history explanations, and compatibility wrappers in the
    active-development files.
- 2026-06-14: Exact CGL pairlist acceleration.
  - It is physically safe to skip CGL-specific pair evaluations only when the
    pair is outside the table cutoff plus a rebuild buffer, because the runtime
    spline evaluator still enforces the original table cutoff and taper for
    every active pair. This is a neighbor-cache optimization, not a shorter
    force-field cutoff.
  - Rebuilding when any endpoint moves by at least half the buffer preserves
    the standard pairlist safety condition: a pair that was outside
    `cutoff + buffer` cannot enter the true `cutoff` before the list is
    rebuilt.
  - The first 1RKL coarse timing check after this change ran at
    `3335.10 us/systems/step` for a 5000-step stage-7 trajectory, bringing the
    CGL workflow into the same range as the current full-lipid timing instead
    of being dominated by dense CGL loops.
- 2026-06-14: Phase 14 cleanup verification lesson.
  - When removing an H5 reader fallback, audit every writer that can create
    that node, including transient workflow helpers. The cleaned
    `restraint_position` runtime correctly requires `spring_const_xyz`, but
    stage-7 burn-in still wrote the old scalar `spring_const` until the short
    hybrid smoke exercised that path.
  - Wrapper scripts are part of the schema surface. Removing generator options
    such as fit-relax controls also requires updating build wrappers before
    regenerated-H5 validation can be trusted.
- 2026-06-14: Current active CGL gap/orientation diagnostic.
  - The reported central CGL bilayer gap in active coarse outputs is not a
    physical void in the HDF5 trajectory. Using stored CGL direction signs for
    leaflet assignment, active 1AFO and 1RKL coarse stage-7 finals have no CGL
    flips, no leaflet crossings, and negative tail-centerline gaps
    (`-4.897 A` and `-3.807 A` respectively).
  - Median-z leaflet splitting can falsely report flips in wrapped/tiled CGL
    systems. For CGL orientation diagnostics, use
    `input/potential/compose_vector6d/direction` when present.
  - The VTF files already carry `radius 4.114` on synthetic CGL display atoms.
    A viewer/rendering mode that ignores atom radius or emphasizes only sparse
    rods/points can still show an apparent center gap even though the DOPC bead
    envelope represented by the CGL vector overlaps across the midplane.
- 2026-06-14: Physical-integrity audit of uniform tempered PMF.
  - The active uniform tempered-PMF implementation satisfies the requested
    physical-design constraints in the audited artifacts: no CGL twist
    coordinate, no standalone CGL orientation potential, no active force cap,
    no arbitrary interaction scaling, no hidden-bead relaxation, and no
    excluded-area/nonnegative projection.
  - Installed `dopc.h5` CGL-CGL, SC-CGL, and CGL-particle tables all use
    `tempered_boltzmann_free_energy` at `tau=10.0`, `fit_relax_steps=0`,
    `sample_dist_min_nm=1e-6`, and no cap/scale/twist attrs.
  - Installed `sidechain.h5` SC-particle table uses
    `tempered_boltzmann_free_energy` at `tau=10.0`, `fit_relax_steps=0`,
    `relaxation=rigid_rotated_geometry`, and no cap/scale attrs.
  - Runtime Phase 11 `.up` audit passed: no CGL twist/orientation-potential
    node, `martini_potential/force_cap=0`, `protein_env_interface_scale=1`,
    SC-env force caps zero, no production restraint node, and duplicate
    residue/rotamer row multiplicity is one for both SC-CGL and SC-particle.
  - A stale unread `nonprotein_hs_force_cap=100.0` metadata default was found.
    It was not read by the active C++ runtime, but it was changed to `0.0` in
    the generator default and Phase 11 checkpoint attrs to avoid misleading
    future audits.
- 2026-06-14: User decision after production-PMF failure.
  - The model requirement is now one hidden-state reduction method for all
    four interaction classes. Since production-temperature PMF failed on the
    first coarse workflow, the next active experiment is uniform tempered PMF
    with `tau=10.0` for CGL-CGL, SC-CGL, CGL-particle, and SC-particle.
  - Dynamics caveat: a tempered PMF is a configurational coarse-graining over
    unresolved rigid orientations. It can improve structural transferability,
    but it is not expected to preserve microscopic kinetics or diffusion
    constants automatically. Dynamics claims require separate calibration or
    comparison to reference trajectories.
  - Focused 1RKL diagnostic with uniform tempered PMF passed structural
    metrics after correcting a manual restart stage-label mistake: no CGL
    flips/crossings, aligned-z p05 `0.845`, same-leaflet NN p05 at least
    `6.886 A`, protein hbond last20 `34.73`, and Rg last20 `12.72 A`.
  - The dynamics mismatch remains: compared with active full-resolution DOPC
    COM motion, CGL lateral MSD is about `3x` larger at lag `15.0` time units
    (`1.400` vs `0.470 A^2`). The practical resolution is an effective-time
    calibration for lipid lateral dynamics, or a separate dissipative/friction
    model; changing the PMF averaging method alone cannot guarantee kinetic
    matching.
  - Restart lesson: `stage_label` used for files/logs is not the same as the
    runtime `/input/stage_parameters/current_stage`. Manual production
    restarts must leave `current_stage=production`; setting it to a numeric
    label such as `7.1` can activate the wrong runtime stage semantics and
    produce invalid energies.
- 2026-06-14: Universal production-temperature PMF validation result.
  - Applying production-temperature PMF hidden-state averaging uniformly to
    CGL-CGL, SC-CGL, CGL-particle, and SC-particle fails on the first coarse
    workflow tested (`run_sim_1rkl.sh`).
  - The failure is not a workflow crash. The run completed, but the final CGL
    bilayer collapsed: aligned-z min/p05/mean `-0.063/0.448/0.821`,
    `bad_parallel=7`, `bad_flip=1`, leaflet crossings `4/4`, and same-leaflet
    nearest-neighbor min/p05 `2.285/2.886 A`.
  - Protein observables failed at the same time: hbond first/final/min/last20
    `28.74/6.93/0.56/6.43` and Rg first/final/last20
    `12.56/10.40/10.37 A`.
  - This supports the earlier concern that production-temperature pair PMFs are
    too attractive for dense CGL-CGL packing because the hidden axial states are
    independently optimized by the pair PMF. A single production-PMF rule across
    all pair classes is therefore not acceptable without a new physical
    representation.
  - SC-CGL implementation/cost check: the generated SC-CGL table has 18
    explicit sidechain types and `18 x 1 x 2541` spline parameters
    (`21 x 11 x 11` per residue). It does include CGL hidden frame sampling,
    but `sidechain_bead_frame_count=1` in the active table and ALA/GLY are not
    explicit sidechain types. Therefore SC-CGL need not dominate wall time over
    CGL-CGL, whose grid is `120 x 9 x 9` and whose hidden samples evaluate two
    14-bead DOPC geometries.
- 2026-06-14: Universal tempered-PMF validation result.
  - Applying the same `tau=10.0` tempered-PMF hidden-state reduction to
    CGL-CGL, SC-CGL, CGL-particle, and SC-particle is not a valid production
    model for the current Phase 9 validation matrix.
  - All four requested workflows completed, and the proteins/full-resolution
    bilayers remained stable by hbond/Rg and DOPC leaflet metrics. The failure
    is specific to coarse CGL orientation: final `1AFO` coarse had
    `bad_parallel=2`, `bad_flip=2`; final `1RKL` coarse had
    `bad_parallel=4`, `bad_flip=4`.
  - The bad CGL rods were present by the first recorded production frame or
    appeared during the stage-7 burn-in, while same-leaflet spacing and leaflet
    crossing checks stayed acceptable. This points to loss of orientation
    transferability when the high-temperature PMF reduction is applied to the
    hybrid CGL-target/SC-CGL surface, not to a global bilayer collapse.
  - Current lesson: keep the dense CGL-CGL tempered-PMF choice separate from
    the less densely packed SC-CGL, SC-particle, and CGL-particle reductions
    unless a new physical representation is validated. A single hidden-state
    averaging rule across all pair classes is too blunt for this model.
  - User correction: do not overinterpret the Phase 9 flipped CGLs as proof
    that tempered CGL-CGL itself is wrong; CGL-CGL was already using tempered
    PMF before this experiment, and rare flips can be seed sensitive. The
    cleaner preferred experiment is production-temperature PMF, without
    tempering, applied uniformly to CGL-CGL, SC-CGL, SC-particle, and
    CGL-particle.
  - Phase 10 implementation check: SC-CGL generation is one task per explicit
    sidechain type, not one task per residue instance in a simulation. The
    active MARTINI sidechain set is 18 types; ALA and GLY have no explicit
    sidechain bead geometry in the orientation library/mapping and are skipped.
    CGL-CGL can still appear more expensive because each unresolved sample
    evaluates a 14-by-14 DOPC bead pair set, while SC-CGL evaluates 14 times
    the sidechain bead count.
- 2026-06-11: 1AFO full-lipid stage-7 secondary-structure fix.
  - VTF distortion in `outputs/martini_1afo_hybrid_full/1afo.stage_7.0.vtf`
    was real trajectory geometry and originated during stage-7
    production-Hamiltonian burn-in. The damaged production frame 0 was the
    promoted burn-in endpoint, not a VTF export artifact.
  - Full-resolution hybrid protein BB proxy types must be structure-derived
    when no protein ITP secondary structure is supplied. The old fallback
    treated the 1AFO helices as coil and typed most BB proxies as `P5`; the
    corrected 1AFO mapping has 61 structure-geometry residues, 7 coil fallback
    residues, and 4 terminal charge overrides.
  - Stage-6 protein position restraints must select `protein_membership >= 0`.
    The previous selector restrained environment atoms (`protein_membership ==
    -1`) and left the protein unrestrained.
  - Structure-derived BB typing alone was not sufficient: the unrestrained
    stage-7 burn-in still promoted a damaged input. Virtual BB proxies stayed
    at the N/CA/C/O COM, so proxy drift was not the cause.
  - Accepted workflow fix: stage-7 burn-in may restrain protein positions as an
    equilibration protocol while all SC-env/BB-env interactions remain
    full-strength. The restraints must be removed before recorded production;
    the accepted active 1AFO full-lipid stage-7 production file has no
    `/input/potential/restraint_position` group.
- 2026-06-11: CGL-CGL training path snapshot.
  - The current CGL-CGL force-field table is generated from two rigid,
    canonical 14-bead DOPC geometries using direct dry-MARTINI LJ+Coulomb bead
    energies. The active path forbids hidden-bead relaxation
    (`fit_relax_steps=0`) and averages unresolved azimuthal/bead-frame samples
    as a tempered two-body PMF.
  - The current SC-CGL table also integrates out hidden SC/CGL bead-frame
    states, but it uses weighted energy expectation rather than a PMF or
    tempered PMF. This distinction should be preserved in methods text and
    reports.
  - Installed `parameters/dryMARTINI/dopc.h5` records
    `azimuthal_count=4`, `cgl_bead_frame_count=8`, `n_radial=120`,
    `n_angular=9`, `sample_dist_min_source=numerical_zero_guard_only`, and
    `energy_transform=log1p_reduced_tempered_pmf` for
    `cg_lipid_table/cg_lipid_pair`.
- 2026-06-10: Correction on CGL gap diagnosis.
  - User correction: a visible CGL gap with a horizontal particle must be mapped
    from VTF atom/rod ordinal back to the HDF5 CGL/CGLD indices before calling
    it a visualization-only issue. Aggregate tail-gap or percentile alignment
    metrics can miss a single visually obvious bad lipid.
  - Rule: for reported CGL visual defects, parse the VTF, identify bad displayed
    rod ordinals, map them to `compose_vector6d/elem_index` and
    `orientation_index`, then compare prepared, minimized, burn-in, and
    production HDF5 geometry.
  - The active `martini_1afo_hybrid` horizontal particle was a real trajectory
    geometry issue: VTF ordinal `104` / residue `130` mapped to HDF5 CGL/CGLD
    `464/736`. It started upright in the prepared file, but had an opposite
    leaflet CGL only `2.542 A` away laterally. Initial CGL conditioning must
    include opposite-leaflet lateral de-overlap derived from the DOPC
    perpendicular bead envelope, not only same-leaflet spacing.
- 2026-06-10: CGL VTF display envelope.
  - CGL centerline tail-gap metrics alone can be visually misleading in VMD
    because `martini_extract_vtf.py` emits synthetic head/tail rod atoms, while
    the physical CGL table represents a resolved DOPC bead envelope with
    nonzero perpendicular radius. If the VTF omits CGL display radii, a viewer
    can show an apparent bilayer gap even when the centerline tails overlap and
    the resolved bead envelope overlaps by several Angstroms.
  - Rule: CGL visualization should carry `max_perp_radius_ang` as display
    radius metadata in VTF output. This is a rendering attribute only; it must
    not be used as a simulation cap, interaction scale, or added potential.
- 2026-06-10: Phase 3 artifact root causes.
  - CGL leaflet classification in hybrid preparation should use the CGL
    direction sign when both leaflet signs are available. Median-COM splitting
    can misclassify wrapped or tiled lower-leaflet CGLs whose COMs sit near the
    bilayer center, creating apparent midplane/horizontal trapped lipids before
    production dynamics starts.
  - CGL initial z conditioning should use the DOPC-derived display-tail
    zero-gap scale (`2 * tail_projection_ang`) as the default. The older
    `2 * tail_projection_ang + contact` spacing created a visible central
    display gap, while removing z conditioning entirely left over-interdigitated
    CGL COM planes that could relax into horizontal rods.
  - Hybrid protein-lipid packing clearance should be derived from the active
    dry-MARTINI DOPC maximum nonbonded contact distance. The old `4.5 A`
    default was below the derived `6.959 A` DOPC contact, so full-lipid 1AFO
    could enter stage-7 burn-in with protein-lipid contacts inside the physical
    dry-MARTINI core and lose secondary structure.
- 2026-06-10: Method-section documentation rule.
  - User correction: paper methods should present the accepted physical model
    as one logical derivation, not as accumulated patches or debugging history.
    When updating `cg_lipid_potentials.tex`, remove superseded mechanisms
    instead of explaining why each old workaround was abandoned.
  - User correction: do not imply the CGL potential is tied to a global bilayer
    direction. Bilayer-normal alignment is a validation criterion for the
    tested membrane systems; the actual CGL tables use local pair coordinates
    and should be described as rotationally invariant.
- 2026-06-10: Phase 2C physical acceptance findings.
  - Rule: active dry-MARTINI table builds must not use finite physical-distance
    floors such as `0.10 nm`. A near-zero `1e-6 nm` guard is acceptable only as
    a numerical singularity guard for exact bead overlap and must be recorded
    as metadata (`sample_dist_min_source=numerical_zero_guard_only`).
  - Accepted SC-CGL representation: extended-support full tensor over the three
    runtime coordinates `(r, sidechain direction angle, CGL direction angle)`,
    with an invertible transformed control. This replaces the unsafe 4-mode
    separable SC-CGL representation without adding normalization, interaction
    scaling, hidden relaxation, capping, or a standalone CGL orientation
    potential.
  - Accepted SC-particle representation: direct full radial-by-angular rotamer
    table consumed by the runtime node. The older separable
    `rotamer_angular_energy_kj_mol` artifact remains only as legacy H5 data and
    is not used when `rotamer_full_energy_kj_mol` is present.
  - Accepted CGL-CGL representation: three runtime coordinates, not four
    runtime parameters. CGL-CGL uses an explicit tempered two-body PMF over
    unresolved axial/azimuthal bead-frame coordinates; this is a documented
    two-body representability choice, not many-neighbor normalization or an
    added orientation potential.
  - Fresh no-floor validation completed for 1RKL and 1AFO in both CGL and
    full-resolution lipid modes. Production logs were finite and stable, CGL
    lipids had no leaflet crossings, full DOPC leaflet separations remained
    stable, and protein hbond/Rg observables stayed bounded through production.
- 2026-05-26: Stage-7 1AFO/1RKL diagnostic root causes.
  - The protein coordinates in `*.stage_7.0.prepared.up` match the HDF5
    reference/PDB mapping for 1AFO and 1RKL; the apparent VTF frame-0 damage is
    the promoted endpoint of stage-7 burn-in, not a raw PDB import mismatch.
  - Full-resolution lipid mode exposed an SC-particle table factorization bug:
    the old SVD residual created runtime reconstructions as low as
    `-4.7e11 kJ/mol`, which is an unphysical lipid-sidechain well and explains
    helix bending/unfolding during burn-in.  The corrected decomposition keeps
    the orientation-average radial term and signed leading SVD mode in the
    resolved range, shifts rows upward only to prevent falling below the sampled
    dry-MARTINI minimum, and makes unresolved hard-core rows finite radial
    barriers to avoid float-storage cancellation wells.
  - User correction: do not add a standalone CGL orientation potential.  CGL
    orientation behavior must come from the CGL--CGL spline in `dopc.h5`.
    The CGL--CGL table should retain the resolved dry-MARTINI lipid--lipid
    attractions instead of stripping a radial background and replacing it with
    a separate correction.
  - User correction: do not use hidden-bead relaxation when tabulating CGL--CGL
    or SC--CGL direction-vector tables.  The table entries should be direct
    dry-MARTINI nonbonded energies from rotated full-resolution lipid and
    sidechain bead geometries at the sampled direction vectors.
  - HDF5 parameter generation should be atomic: write to a sibling temporary
    file and replace the production `.h5` only after the build succeeds, so an
    interrupted expensive table build cannot leave `dopc.h5` as a placeholder.
  - In sandboxed environments process-pool semaphores can be unavailable.  The
    dry-MARTINI table builder should keep process workers as the default for
    M1/Slurm but fall back to a thread pool with the same worker count instead
    of serializing the build.
- 2026-05-26: User correction on direction-vector spline variables and
  table-build speed.
  - Rule: describe and implement the CG potentials in terms of the runtime
    direction vectors. Table generation must rotate the resolved SC and CGL
    bead models over those sampled direction vectors: one orientable object for
    SC-particle and CGL-particle, both orientable objects for SC-CGL and
    CGL-CGL.
  - Rule: do not present around-vector bead-frame averaging as the physical
    requirement. If an implementation samples unresolved bead-frame rotations
    around a direction vector, document it only as secondary quadrature for
    resolved bead placement and keep the direction-vector spline coordinates
    primary.
  - Parallel table construction should follow Slurm CPU allocation conventions:
    prefer an explicit `UPSIDE_MARTINI_TABLE_WORKERS`, otherwise use
    `SLURM_CPUS_PER_TASK`/`SLURM_CPUS_ON_NODE`, and only then local CPU count.
    This matches the sweep infrastructure pattern where `--cpus-per-task` is
    set by generated Slurm scripts.
- 2026-06-09: User correction on Phase 2C scope.
  - Rule: the accepted CGL-CGL method may use three runtime table parameters,
    but the implementation priority remains a transferable two-body model. Do
    not make many-neighbor normalization the first-choice fix, because the same
    method must extend to SC-CGL where that normalization has no clear physical
    analogue.
  - Rule: SC-CGL, SC-particle, and CGL-particle must be rebuilt and tested with
    the same physical direct-geometry philosophy as CGL-CGL: no capping, no
    hidden relaxation, no empirical interaction disabling, and no added
    standalone CGL orientation potential.
- 2026-06-09: Phase 2C direct-table runtime findings.
  - CGL-target Boltzmann-weight controls are numerically fragile for charged
    target hard cores. In float32, hard-core rows underflow to a flat reduced
    free-energy plateau with weak restoring force. An invertible
    `log1p((E - E_ref) / kBT)` reduced-PMF control preserves the same physical
    two-body PMF while keeping the spline representable.
  - CGL runtime spline evaluation previously clamped radial coordinates below
    one knot spacing with zero derivative. That is an implicit force cap even
    when the H5 table contains a physical repulsive core. Low-end radial
    evaluation must preserve a restoring slope.
  - The latest 1RKL coarse probe shows radial extrapolation alone is not
    sufficient: stage-6 minimization still drives same-leaflet CGL centers into
    sub-Angstrom contacts. The next diagnosis must compare raw direct CGL-CGL
    PMF samples, fitted controls, and runtime derivatives at the collapsed
    geometries before changing the model.
  - CGL-CGL direct PMF sampling must include the overlap core. Starting the raw
    grid at `5 A` missed a physical `10^3-10^4 E_up` wall at `0.35-2 A` for
    the collapsed angular sector and left sub-`5 A` behavior to extrapolation.
  - Raw-energy CGL-CGL fitting is also numerically unsafe once the physical core
    is included: true extreme-angle hard cores at long COM range can bleed into
    ordinary angular sectors. The same invertible `log1p` reduced-PMF control
    used for CGL-target keeps the two-body physical PMF while avoiding that
    spline artifact.
  - Table H5 metadata is not sufficient validation for transformed spline
    controls. The generated runtime `.up` node must also carry
    `log1p_reduced_transform`, `boltzmann_temperature_upside`,
    `energy_transform`, `spline_control_quantity`, and the
    `pair_interaction/reference_energy_eup` dataset. If injection drops those
    fields, Upside silently interprets log controls as raw energies and weakens
    the physical hard core by orders of magnitude.
  - Production-temperature Boltzmann CGL-CGL PMF still collapsed the focused
    dense CGL-only bilayer after runtime transform injection was fixed, while
    direct energy expectation was too stiff and launched the bilayer. The
    current passing CGL-only table is an explicit tempered two-body PMF over
    unresolved axial/azimuthal coordinates. This is not many-neighbor
    normalization, hidden relaxation, capping, or a standalone orientation
    potential, but it is a real modeling parameter and must remain visible in
    H5 metadata and reports.
  - Superseded by the 2026-06-10 full-tensor/no-floor validation: the
    4-mode separable SC-CGL representation was not safe for extended support,
    but the extended full tensor with transformed control now passes the fresh
    1RKL/1AFO CGL/full validation matrix.
- 2026-05-25: CGL-only bilayer validation of updated CGL-CGL runtime path.
  - The validation harness should exercise installed production spline tables
    (`particle.h5` plus `dopc.h5`) by default; local table refitting is a
    separate rebuild check and can be much slower.
  - A 72-lipid CGL-only DOPC bilayer with only CGL-CGL and CGL-CGLD runtime
    terms passed 200-step and 2000-step NVT checks after removing the separate
    leaflet-normal orientation spring. The 2000-step check showed no leaflet
    flips/crossings and kept CGL-CGLD lengths close to the canonical
    `11.139 A` orientation length.
- 2026-05-25: User correction on CG table orientation sampling.
  - Rule: CGL-CGL, CGL-SC, CGL-target, and SC-particle spline tables must not be
    derived from one default lipid or sidechain orientation. They must be fit
    from resolved dry-MARTINI bead energies after rotating the full-resolution
    lipid and sidechain bead models over the sampled CG orientation variables.
  - Runtime CG lipid geometry must not be re-derived from the first packed
    lipid in a simulation box. The packed distribution may provide initial
    positions and leaflet directions, but the CG model geometry must be
    canonical and consistent with the resolved model used to build the tables.
  - Working rule: when diagnosing full-resolution versus single-particle
    lipid mismatches, first distinguish table-generation physics from runtime
    particle instantiation. A table can be orientation-resolved while runtime
    `compose_vector6d` metadata is still physically inconsistent.
- 2026-05-25: User correction on CGL orientation ownership.
  - Rule: do not add a separate CGL leaflet-normal orientation potential when
    the CGL-CGL spline table already has lipid orientation coordinates, unless
    a distinct many-body coarse-grained orientation term is explicitly derived
    and documented. Otherwise the orientation physics is double-counted.
  - Audit finding: `parameters/ff_2.1/sidechain.h5` stores one sidechain
    placement bead per rotamer, while MARTINI sidechain definitions can have
    multiple beads. A resolved sidechain table cannot be made physical by
    inventing missing MARTINI bead positions; it needs a real geometry source
    or an explicit one-placement-bead model statement.
  - Resolution: use the Upside rotamer center/vector as the sidechain pose and
    expand it with MARTINI sidechain bead offsets derived from `martinize.py`
    bond lengths/connectivity. SC-env and CGL-SC table generation then evaluate
    all resolved MARTINI sidechain beads for each rotamer.
- 2026-05-25: User correction on sidechain table ownership.
  - There is no dry-MARTINI SC-SC interaction table in this workflow. The
    SC-particle table is shared in both lipid modes for non-CGL environment
    particles. Full-resolution lipid mode also uses SC-particle for explicit
    lipid beads, while CGLipid mode replaces only the explicit-lipid part with
    SC-CGL. Do not describe the ordinary Upside rotamer `pair_interaction`
    table as the dry-MARTINI SC-SC path for this task.
- 2026-05-25: 1AFO first-frame coarse output after carrier routing fix.
  - The user-observed first-frame bend in
    `example/16.MARTINI/outputs/martini_1afo_hybrid` was not caused by the
    minimizer. On a copied current-binary minimization-only check, fragment
    bends changed only `40.50/42.89 deg -> 40.07/42.19 deg`.
  - The first production frame is the promoted endpoint of the 40k-step stage-7
    burn-in, not the raw minimized structure. This matters for diagnosis:
    frame 0 can reflect burn-in dynamics even when the minimizer is harmless.
  - The existing generated coarse output was stale relative to the current C++
    force-routing binary. Its old stage-7 minimization log started at
    `162.96 E_up`, while the same prepared file evaluated with the current
    binary starts at `2867.25 E_up` and minimizes to a different state.
  - Regenerating the coarse stage-7 handoff with the current binary kept the
    physical interface active and produced production-frame-0 fragment bends
    `28.98/21.38 deg`, first/last hbond `82.37/84.70`. The ignored generated
    checkpoint and VTF under `outputs/martini_1afo_hybrid` were replaced with
    this current-binary result.
- 2026-05-25: User correction on N/CA/C/O force ownership.
  - Rule: do not infer "BB proxy only" CGL targeting as excluding force
    contributions from backbone carrier atoms. N/CA/C/O carriers must
    accumulate sidechain/rotamer and backbone-environment contributions
    additively.
  - Corrected diagnosis: mapped N/CA/C/O carrier atoms have `ROLE_BB` and a BB
    map index, so the C++ derivative path was treating direct carrier
    sensitivities as projectable virtual BB-proxy sensitivities. Direct carrier
    gradients must be added directly to the carrier coordinate; only true
    virtual BB proxies are projected through the BB map.
  - Rejected model: adding all N/CA/C/O carriers as independent CGL targets
    double-counts the hidden backbone representation. A copied 1RKL probe with
    carrier CGL targets collapsed to final hbond `6.5` over 40k steps; after
    the C++ projection fix, a 10k carrier-target probe still fell to hbond
    `9.6`.
  - Accepted model: CGL target sites include BB proxies and non-protein point
    targets, not independent N/CA/C/O carrier copies. Carrier atoms still
    receive sidechain/rotamer forces and projected BB-proxy forces. If another
    node ever applies a direct carrier gradient, it stays direct and is not
    re-projected.
  - Verification: the copied 1RKL proxy-target/additive-force probe held hbond
    `31.43` at the 40k-equivalent endpoint with kinetic ratio `0.989`, versus
    the stale coarse output final hbond `23.71`. The 1AFO current injection
    produced `182 CGL x 170` targets (`72` BB proxies plus `98` ions), with
    both fragments retaining `Qd/+1` and `Qa/-1` endpoints.
- 2026-05-25: 1RKL single-particle lipid bend after charged-endpoint fix.
  - Terminal endpoint typing was already correct in the current generated HDF5:
    1RKL has `Qd/+1` and `Qa/-1` endpoint BB proxies in both coarse and full
    modes; 1AFO has the same charged endpoints on both fragments
    (`0-35` and `36-71`) in both modes.
  - Superseded assumption: full/coarse mismatch should not be fixed by making
    N/CA/C/O carrier atoms environment-inactive. The physical invariant is
    additive force accumulation with correct direct-vs-projected C++ routing.
- 2026-05-25: 1AFO single-particle lipid helix bend after the CGL-SC rotamer
  fix.
  - Current `example/16.MARTINI/outputs/martini_1afo_hybrid` artifacts are not
    stale: stage-7 contains `cg_lipid_rotamer_sc`, so the remaining 1AFO bend is
    a second issue.
  - Geometry audit showed the coarse single-particle chain-1 helix is already
    bent at stage-7 production frame 0 (`~89.7 deg` half-chain angle), while the
    full-lipid chain remains near `45 deg`.  Stage-6 output and stage-7
    prepared input are not bent; the kink is introduced during the 40k-step
    stage-7 burn-in, not by handoff minimization.
  - The driver is the charged terminal BB target projection: explicit-DOPC
    `Qd`/`Qa` directional PMF rows are non-transferable as an additive
    single-particle CGL attraction after lipid headgroup, solvent, and ion
    response have been integrated out.
  - Physical fix: keep `Qd`/`Qa` charged in the dry-MARTINI protein definition
    and generic interactions, but split charged BB proxy targets into
    `cg_lipid_target_charged_bb_excluded_volume` with nonnegative CGL controls.
    Ordinary uncharged BB/protein targets still use the full CGL-target table;
    ions continue to use their excluded-volume target node.
  - Verification on a regenerated copied 1AFO stage-7 file produced 68 ordinary
    CGL targets, 4 charged-BB excluded-volume targets, and 98 ion
    excluded-volume targets.  A full 40k-step burn-in from the regenerated
    minimized handoff ended with chain angles `32.3 deg` and `36.9 deg`, versus
    the old coarse output's `37.3 deg` and `85.7 deg`.
- 2026-05-24: Secondary-structure divergence between full lipid and
  single-particle lipid workflows.
  - Full-resolution lipid mode feeds explicit DOPC sidechain environment
    contacts into the sidechain rotamer one-body path through
    `martini_sc_table_1body`; in the checked 1RKL output this node saw 4041
    environment atoms.
  - Single-particle lipid mode had `martini_sc_table_1body` seeing only ions
    (93 atoms in the checked 1RKL output), while the dry-MARTINI-derived
    CGL-SC table was evaluated as a standalone `cg_lipid_sc` CB potential.
    This kept CGL-SC active but bypassed rotamer probabilities, so the coarse
    workflow used a different sidechain/lipid coupling path from full mode.
  - Physical fix: keep the same CGL-SC spline active but inject it as the
    rotamer one-body coordinate node `cg_lipid_rotamer_sc`, appended to the
    `rotamer` arguments. The generated standalone `cg_lipid_sc` node is removed
    on reinjection to avoid double-counting.
  - Runtime lesson: Upside node type names are prefix-matched, so new registered
    node names must not start with an existing registered type such as
    `cg_lipid_sc`.
- 2026-05-24: MARTINI cleanup verification.
  - `example/16.MARTINI/test_cg_lipid/run_test.py` is a practical short
    executable check for the cleanup because it covers MARTINI table generation,
    stage conversion, CG-lipid node injection, and minimization without running
    full production trajectories.
  - The minimal DOPC+GLY fixture should pass `active_residue_names=("GLY",)` to
    table generation. Building every sidechain table is unnecessary for that
    fixture and turns a focused smoke test into a long parameter-generation run.
- 2026-05-24: User correction on `py/martinize.py`.
  - `py/martinize.py` is third-party vendored code and must not be edited,
    formatted, annotated, deleted, or otherwise touched during MARTINI cleanup.
  - Working rule: isolate any dependency on martinize behind project-owned
    wrappers, but leave the vendored source byte-unchanged.
- 2026-05-24: User correction on cleanup depth.
  - Removing debug switches is not enough for this task. The required bar is a
    branch-vs-master cleanup of all MARTINI Python/C++ additions: straighten the
    workflow logic, match local master style, reduce duplicate implementations,
    and remove generated-looking structure.
  - Working rule: after visible debug cleanup, inspect the actual diff against
    `/Users/yinhan/Documents/upside2-md-master` and normalize new code against
    nearby master patterns before calling the refactor complete.
- 2026-05-24: User correction on `run_sim_1rkl.sh` environment bootstrap.
  - `source.sh` appends variables such as `PYTHONPATH`, `CPLUS_INCLUDE_PATH`,
    `LIBRARY_PATH`, and `LD_LIBRARY_PATH` without default guards.
  - Working rule: shell launchers that source legacy project bootstrap scripts
    must not enable `set -u` until after bootstrap has completed, unless every
    sourced script has been audited for unset-variable safety.
- 2026-05-22: User correction on stage-7.4 in-plane lipid orientations.
  - The CGLD-CGL vectors lying in the x-y plane are a real physical model defect, not a VTF-only issue. The previous symmetric VTF display style is preferred and was restored.
  - Current 1RKL stage-7.4 output before the fix had `34` lipids with `|n_z|<0.5` and `14` with `|n_z|<0.25`; the count grew across stage-7 chunks.
  - Root cause: the CGL-CGL table clips all negative tensor controls to zero. That avoids a known additive many-neighbor collapse but also removes enough bilayer orientational cohesion that sparse CGL particles can become near-free rotors.
  - Directly restoring the stored radial attractive background is not acceptable: on a copied stage-7.4 artifact, `cg_lipid_pair` accumulated from about `-1453` to `-3953 E_up` over one 10k-step chunk.
  - Superseded fix: do not add `cg_lipid_leaflet_orientation`. The accepted
    direction is to retain the full resolved dry-MARTINI lipid--lipid
    attraction in the CGL--CGL spline stored in `dopc.h5` and avoid any
    separate one-body or long-range orientation correction.
- 2026-05-21: 1RKL stage-7.4 edge-lipid VTF geometry check.
  - Do not run `obj/upside` diagnostics directly on original `.up` artifacts, even with `--duration-steps 0`; the program can rewrite the `/output` group. Use a copied HDF5 file for derivative checks.
  - The reported edge-lipid oddity in `outputs/martini_test_1rkl_hybrid/1rkl.stage_7.4.vtf` is not a special physical box-edge interaction. HDF5 CGL-CGLD minimum-image lengths are normal (`11.13-12.04 A`), edge and non-edge z distributions are comparable, and no displayed lipid z coordinate leaves the half-box.
  - VTF edge artifacts come from periodic imaging: a continuous single-particle lipid display rod whose center is near an x/y periodic boundary can put one endpoint outside the primary image. The extractor must keep that rod continuous rather than independently wrapping endpoints.
  - `py/martini_extract_vtf.py` now uses the stored per-lipid `display_head_offset_ang` and `display_tail_offset_ang` metadata instead of replacing the display geometry with an artificial symmetric rod.
  - The trajectory still shows a real slow increase in CGL z spread and tilted CGLD-CGL directions across stage-7 chunks. The C++ derivative path is present: CGL spline nodes write 6D sensitivities and `compose_vector6d` propagates orientation sensitivity to CGL/CGLD coordinates.
- 2026-05-21: User correction on latest 1RKL production drift and ion ownership.
  - Rule: when parsing workflow logs, isolate the actual production stage section. `1rkl.0.log` also contains stage-6 and minimization output; using every `total_potential` line overstates the production baseline.
  - Exact production sections in `example/16.MARTINI/1rkl.0-4.log` still are not stationary: logged production moves from `111.75` to `-2.21 E_up` overall, and a copied continuation from `1rkl.stage_7.4.up` has first/last fifth total-potential means `-6.37 -> -41.05 E_up`.
  - The residual drift is not the old ion/CGL-target sink: over the five logged production chunks, `cg_lipid_target` is approximately flat (`-31.51 -> -29.01 E_up`). The main remaining downward terms are `cg_lipid_pair` (`350.80 -> 294.28 E_up`) and `cg_lipid_sc` (`-85.23 -> -151.68 E_up`).
  - Ion ownership check on current 1RKL stage-7 checkpoints: ions have `14,415` generic MARTINI ion-protein pairs, including `11,532` ion-backbone-carrier and `2,883` ion-BB-proxy pairs, plus `4,278` ion-ion pairs. There are zero generic ion-CGL/CGLD pairs because every CGL/CGLD pair is owned by dedicated spline nodes; ion-CGL is represented by `cg_lipid_target_ion_excluded_volume` with all 93 ions as targets and nonnegative controls.
  - Current workflow design is physically motivated but not yet production-stationary. The fix should target additional physical equilibration or the CGL-CGL/CGL-SC relaxation source, not ion-protein exclusions or parameter scaling.
- 2026-05-21: Root cause of the latest 1RKL `-800 -> -1100 E_up` drift.
  - Component parsing confirmed the current drift is dominated by `cg_lipid_target` (for example `1rkl.1.log` target `-1007.93 -> -1280.76 E_up`) while `cg_lipid_pair` is not the main sink.
  - The dominant residual target sink came from terminal charged BB proxy targets (`Qd`/`Qa`) using an explicit-DOPC charged-point attractive well in the additive single-particle CGL projection.
  - User correction accepted: do not simply remove charge from the CGL-target
    interaction while leaving the dry-MARTINI protein definition charged.
    Superseded resolution: keep charged terminal BB types in the dry-MARTINI
    protein definition and route only their direct CGL single-particle
    projection through an excluded-volume target node, because the full
    explicit-DOPC charged-target PMF is not transferable as an additive CGL
    attraction.
  - Validation: the first copied 10k-step continuation re-equilibrated old charged-model coordinates (`85.55 -> -13.30 E_up`, first/last fifth means `59.35 -> 8.89 E_up`). The next continuation had total potential `-13.30 -> -26.54 E_up`, range `[-73.93, 23.57]`, first/last fifth means `-1.60 -> -35.71 E_up`; a third continuation reversed the mean change (`-22.91 -> -12.66 E_up`, range `[-58.95, 63.38]`, kinetic ratio `1.017`). This removes the monotonic multi-hundred-`E_up` target sink and leaves fluctuation-scale wandering after re-equilibration.
- 2026-05-21: User correction after explicit burn-in.
  - Rule: do not classify residual production drift as solved just because copied continuations from older endpoints fluctuate.  Re-run or inspect the exact fresh logs generated after the workflow change and confirm the intended burn-in was actually applied.
  - Fresh 1RKL logs reportedly still move from about `-800 E_up` to about `-1100 E_up` over four trajectories after the burn-in change, so the residual source must be reopened.
- 2026-05-21: User correction on 1AFO multi-chunk drift.
  - Rule: do not judge residual equilibration from a single production chunk when later restart chunks exist.  Compare absolute component baselines across all chunks, including `1afo.1.log`, `1afo.2.log`, and `1afo.3.log`, before concluding a system has plateaued.
  - The previous statement that 1AFO was flat was only true for `1afo.0.log` first-to-last and is insufficient for the user's reported multi-chunk accumulation from about `-700 E_up` to below `-1000 E_up`.
  - Reanalysis showed restart bookkeeping is consistent: momentum is valid, transition starts advance by `10000` steps per chunk, and chunk handoff potentials are continuous.
  - The shared residual is production-Hamiltonian interface relaxation, mainly `cg_lipid_sc` and `cg_lipid_target`, not a continuing ion sink.  Copied continuations from current 1AFO and 1RKL stage-7.3 endpoints fluctuate around the already-relaxed basin.
  - Accepted correction: run an explicit stage-7 burn-in under the same production Hamiltonian before the named production segment.  Promote final burn-in positions and momenta to `/input`, advance `sc_env_transition_step_start`, clear `/output`, and then start logged production.  This preserves physics and reclassifies equilibration rather than hiding it as production.
- 2026-05-21: Follow-up after ion-target split shows residual 1RKL relaxation, not a continuing ion sink.
  - Current `1rkl.0-2.log` still contain substantial relaxation from the corrected stage-7 launch, mostly in `cg_lipid_target` with smaller `cg_lipid_sc` relaxation, but the latest chunk is much closer to stationary behavior: `1rkl.3.log` total potential changes only `-28 E_up` first-to-last.
  - A copied continuation from current `1rkl.stage_7.3.up` final coordinates and momenta ran 10k more steps over 60 public time units with total potential `-1045.92 -> -1062.16 E_up`, range `[-1138.67, -981.46]`, and first/last fifth mean delta `-34.37 E_up`; this is fluctuation-scale behavior compared with the previous thousand-`E_up` ion adsorption sink.
  - Current `1afo.0.log` total potential is flat first-to-last (`-730.55 -> -730.09 E_up`), with compensation among protein, CGL-pair, CGL-SC, and CGL-target components.
  - Ion ownership audit: 1RKL and 1AFO stage-7 HDF5 files have zero generic CGL-ion pairs, while ions remain in generic dry-MARTINI ion/protein and ion/ion pairs.  CGL-ion interactions are owned by the excluded-volume-only target node.
  - Ion spatial audit: final ion-CGL minima are `16.76 A` for 1RKL and `16.83 A` for 1AFO, and ion-z distributions do not collapse toward the membrane.  This supports trusting ion behavior as generic dry-MARTINI ion behavior plus a CGL excluded-volume boundary, not an attractive DOPC-ion PMF.
- 2026-05-21: Residual 1RKL/1AFO energy drift is a CGL-target ion adsorption artifact.
  - The progress log's old `protein_potential` bucket included `cg_lipid_target`, so it made the protein appear to be the large energy sink even when internal Upside protein terms were small and stable.
  - Per-node diagnostics on copied final-frame HDF5 files showed the real total-potential drift is dominated by `cg_lipid_target`, not CGL-CGL or CGL-SC.
  - Splitting target particles by class showed mobile ions dominate the long-term sink: 1RKL stage 7.3 was about `-1558 E_up` from BB targets and `-3401 E_up` from ion targets, and one copied continuation lowered the ion target further to about `-4170 E_up`.
  - Rule: the directional explicit-DOPC CGL-target table is appropriate for BB/protein targets, but mobile ions must not reuse the attractive DOPC-ion energy well as an additive many-ion potential in the single-particle lipid model. The missing hydrated/solvent ion PMF makes this a non-transferable ion adsorption well.
  - Accepted correction: split mobile ions into a separate CGL-ion node derived from the current `cg_lipid_target` table by retaining only nonnegative controls. This preserves excluded volume and leaves BB/protein CGL-target interactions active.
  - Rule: MARTINI progress logs must report `cg_lipid_target` separately from `protein_potential`; otherwise future diagnostics will misattribute environmental target relaxation to the protein.
- 2026-05-20: Predescu integrator convention and MARTINI-specific opt-in dynamics.
  - Predescu et al. 2012 define HOH integrator coefficients so each coefficient set sums to the number of stages `q`; with this convention, `dt` is the average time per force evaluation and one outer integrator application advances `q*dt`.
  - The shared `DerivEngine::integration_cycle(VecArray,float,float,IntegratorType)` coefficient pattern `(1.5-3a, 1.5-3a, 6a)` and `(3b, 3-6b, 3b)` is consistent with that convention for `q=3`; do not change it for MARTINI-only behavior.
  - Rule: do not alter shared C++ integrator semantics used by ordinary Upside/master workflows to fix MARTINI-only dynamics.
  - Rule: when adding a generic CLI switch for MARTINI needs, the default behavior must remain semantic parity for non-MARTINI workflows unless the user explicitly approves broader changes.
  - The hybrid MARTINI issue is not that Predescu coefficients are mathematically wrong; it is that this workflow combines physical dry-MARTINI masses, stiff CGL-CGLD orientation carriers, active tabulated interface forces, and a rigid stage-6 to flexible production release.
  - Removing `/input/mass` from a copied stage-7 artifact failed: total potential `-505.9 -> -24239.58 E_up`, Rg `12.7 -> 14.9 A`, kinetic ratio `5.70`.
  - Scaling pair energies by mass is not physical because it changes the scalar potential and thus the configurational distribution; per-particle mass-weighted forces for a pair would also not be a single conservative pair potential.
  - Superseded timing conclusion: do not restore MARTINI workflow public-time scheduling to `inner_step=3`. The active requirement is fixed 1:1 slow/fast stepping for MARTINI, implemented with explicit `mv` and internal `MARTINI_MD_INNER_STEP=1`, while leaving shared C++ integrator semantics unchanged.
- 2026-05-20: 1RKL stage-7 continuation failure from CGL-CGL many-neighbor attraction.
  - The exact output uses a fixed square stage-7 box; the apparent XY ellipse is membrane morphology/coordinate drift, not anisotropic box-vector drift.
  - The total potential descent is dominated by CGL-CGL: CGL-CGL interacting pairs increase from about `2,794` at stage-7 start to over `6,200` by the end of `7.2`.
  - Rigid-protein copied-HDF5 relaxation still drives `cg_lipid_pair` from about `-90.9k` to `-116.7k E_up`, ruling out protein release timing as the root cause.
  - Rule: explicit DOPC-derived CG-lipid projections must use the same bead-level dry-MARTINI nonbonded cutoff as the generic runtime MARTINI potential; otherwise hidden bead LJ/Coulomb tails become nonphysical long-range CGL attractions.
  - Rule: isolated DOPC-DOPC attractive controls are not automatically transferable to an additive many-neighbor single-particle membrane. CGL-CGL should keep excluded-volume and orientation-compatibility controls; non-transferable attractive CGL-CGL controls act as an energy sink and must be rejected by table metadata.
  - Copied-HDF5 validation with attractive CGL-CGL controls removed kept final CGL nearest-neighbor min/p05/mean `6.278/6.729/7.821 Å`, XY aspect `1.08`, hbond sum `26.75`, CA RMSD `4.16 Å`, and Rg `12.47 Å` over `10k` total stage-7 steps.
- 2026-05-20: User correction after residual 1RKL validation.
  - Rule: do not close a stability fix based only on a sibling validation directory when the user is inspecting a specific output path; audit or regenerate the exact named artifact.
  - Rule: protein secondary-structure stability needs an explicit retention metric such as hbond count plus aligned CA RMSD, not only protein Rg and lipid contact distances.
  - Rule: if an exact output looks unchanged after a code edit, first distinguish stale generated artifacts from a failed physics correction before adding another model change.
  - The exact regenerated 1RKL output used fresh DOPC charge and CGL-target excluded-area attrs, but CGL-SC still became a deep sink (`-1319 E_up`) while hbond retention fell to `1.22` and final CA RMSD reached `5.86 Å`.
  - Patching copied HDF5 CGL-SC controls inside the DOPC contact shell to be non-attractive improved final hbond retention to `17.96`, final CA RMSD to `3.66 Å`, and final CGL-SC energy to `-73.6 E_up` without collapsing CGL-CGL spacing.
  - Rule: for CGL-SC, sub-contact rows are unresolved overlap with the hidden DOPC molecule; preserving attraction there is not physical even if the relaxed explicit-bead sample is attractive.
- 2026-05-20: Residual 1RKL protein instability after the first minimization fix.
  - Stage-6 minimization was being written only to `/output`; without promoting minimized `/output/pos[-1]` and `/output/box[-1]` to `/input`, the following MD stage started from the original unrelaxed coordinates.
  - DOPC CG table fitting must use molecule charges from `dry_martini_v2.1_lipids.itp`: NC3 is type `Q0` but charge `+1`, PO4 is `-1`, and the remaining beads are neutral. Type-inferred charges are insufficient for lipid tables.
  - CGL-target charged rows can contain attractive spline-control overshoot inside the DOPC contact excluded-area domain. Constraining only those excluded-area controls to nonnegative values is a table-fitting correction to preserve dry-MARTINI excluded volume, not an interaction cap or scale.
  - Direct copied-HDF5 validation on 1RKL improved final CA RMSD from `5.203 Å` to `3.902 Å` and kept Rg near `12.861 Å` over 10k stage-7 steps.
  - Rule: when validating a minimization stage, compare the initial potential of the subsequent MD to the minimizer final potential; otherwise a successful minimization may still be a no-op for production.
  - Rule: do not infer molecule bead charges solely from MARTINI bead type when topology records explicit per-particle charges.
- 2026-05-19: Protein stability regression after CGL-CGL table changes.
  - The reported 1RKL/1AFO protein instability is downstream of lipid collapse, not evidence that protein-env interactions should be disabled or weakened.
  - 1RKL bad stage 7 reaches CGL-CGL nearest-neighbor p05 `0.215 Å` and CGL-protein min `0.267 Å` by the final frame.
  - A source-equivalent CGL-CGL unresolved-core patch plus minimized stage-6 handoff keeps final CGL-CGL min/p05 `5.968/6.200 Å` and CGL-protein min/p05 `5.197/9.729 Å`.
  - Rule: an unresolved spline core can be force-field-derived and still insufficient if it copies a weak resolved boundary row. For CGL-CGL, unsampled inner rows must be tied to a dry-MARTINI sampled short-range energy expectation, with attrs that stale guards can verify.
  - Rule: a short minimization under the full active Hamiltonian is a physical handoff relaxation. It is not equivalent to a cap, restraint, interaction exclusion, or force scaling.
  - Rule: for long MARTINI workflows, validate physics fixes by patching copied stage HDF5 artifacts and comparing geometry metrics against the bad output before spending time on full reruns.
  - Rule: always compare lipid nearest-neighbor and lipid-protein distances along with protein Rg/RMSD; protein metrics alone can hide that the primary failure is membrane collapse.
- 2026-05-18: CGL-CGL excluded-area validation.
  - A WCA excluded-area projection is physical only if the fitted spline does not try to represent the divergent zero-distance limit; sampling it to `1.4 Å` created enormous coefficients and immediate bilayer blow-up.
  - The accepted CGL-CGL table samples the compressed CGL-center region from `5.0 Å`, adds a DOPC-contact WCA term with `epsilon=kBT` and `sigma=r_contact/2^(1/6)`, and constrains only the DOPC excluded-area spline rows to nonnegative controls.
  - Reduced CGL-CGL angular grids are not trustworthy for bilayer morphology: the coarse `8r x 5^2 theta` table still allowed same-leaflet overlap. CGL-CGL must keep the full `16r x 7^2 theta x 2^2 phi` sampling grid even when other table parts use reduced settings.
  - Bilayer-only validation passed for 5000 steps after this change with final aligned-z min/p05/mean `0.431/0.695/0.898`, no flips/crossings, same-leaflet nearest-neighbor min/p05 `4.937/5.769 Å`, and CGL-CGLD RMSD `0.275 Å`.
- 2026-05-18: User correction after PMF cap-removal attempt.
  - Rule: a pair PMF over unresolved lipid azimuths can be physical for an isolated two-body system but still nonphysical when reused as an additive pair potential, because it can overcount many-body rotational relaxation and drive lipid aggregation.
  - Rule: if protein remains stable while the user expects protein-env stress, verify that SC-env and BB-env interactions are actually present and nonzero before interpreting protein stability as a model success.
  - Rule: do not accept protein stability as evidence when lipid/env interactions are wrong; inspect component interactions and generated output morphology.
- 2026-05-18: User correction on CG-lipid energy caps.
  - Rule: do not introduce fixed numeric caps as model physics. If a spline table needs finite behavior in an unresolved region, derive that boundary from the underlying dry-MARTINI model or document the physical statistical-mechanics rationale.
  - Rule: successful stability is not enough if it depends on an arbitrary cap; replace it with a force-field-derived boundary or remove it.
- 2026-05-17: Focused hybrid validation after preserving all protein-env paths identified two implementation-level physics errors.
  - `cg_lipid_target` should include the virtual MARTINI `BB` proxy site, but not also the `N/CA/C/O` backbone carrier atoms as independent dry-MARTINI targets. Including both double-counts the same backbone physical site.
  - Forces from CGL-target interactions on a virtual `BB` proxy must be projected through the same hybrid BB proxy-to-carrier map as the standard Martini BB-env potential.
  - A standalone 50k-step bilayer-only run passed, while the hybrid copied artifact collapsed only when target forces could drive lipids through the old soft CGL-CGL core; this isolates the remaining instability to the unresolved bead-overlap core, not to missing or excluded protein-env interactions.
  - The initially tested hard finite CGL-CGL unresolved core (`5000 kJ/mol`) stabilized the copied hybrid run, but the user correctly rejected it as a nonphysical fixed cap.
  - Corrected rule: use dry-MARTINI-derived energy-expectation samples and first-resolved energy-expectation boundary rows for unresolved spline regions; reject stale tables lacking those attrs rather than reusing cap-based tables.
- 2026-05-17: User correction on CGL-target inclusion during protein-instability debugging.
  - Rule: do not exclude `BB` or other protein/environment targets from `cg_lipid_target` to stabilize the protein.
  - Rule: if a newly added CG-lipid interface destabilizes a previously working protein-env path, correct the physics, table construction, force routing, or restart/handoff semantics while keeping the full interface present.
- 2026-05-17: User correction after short copied-HDF5 validation: the fresh `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/1rkl.stage_7.0.vtf` still shows protein instability.
  - Rule: copied-artifact short runs are useful triage, but final acceptance must inspect the exact fresh output path the user names.
  - Rule: when a user reports persistent instability after a proposed physical fix, re-open the diagnosis and verify the exact VTF/HDF5 trajectory before adding more model changes.
- 2026-05-17: User reminder on physical-model integrity during 1RKL bilayer/protein stabilization.
  - Rule: do not interpret a passing short trajectory as acceptable if it was achieved by weakening, scaling, disabling, or restraining interactions.
  - Rule: acceptable fixes must be physical table-construction or handoff-semantics fixes, such as charging bonded deformation during explicit-bead relaxation, clamping unresolved spline cores to excluded-volume repulsion, and using canonical force-field reference geometry for potential tables.
  - Rule: when diagnosing bilayer orientation near a protein, distinguish actual parallel rods (`abs(n_z)` small) from simple global-z leaflet sign mismatches; report absolute bilayer-normal alignment and 3D CGL-CGL nearest-neighbor distances.
- 2026-05-16: User correction on CG-lipid orientation physics.
  - Rule: do not repair bilayer orientation with an orientation spring, z-axis pin, empirical pair scale, or disabled interactions. Treat CGL-CGL many-neighbor normalization as a fallback only after transferable two-body table construction is exhausted, because the same CGL method must extend to CGL-SC where normalization has no clean analogue.
  - Rule: CGLD should be treated as the endpoint of the physical orientation vector; pair interactions must apply derivatives to the vector through the CGLD coordinate, not through a separate orientation parameter.
  - Rule: radial or orientation-averaged CGL-target interactions cannot torque the single-particle lipid correctly. CGL-target interactions need directional tables over `(r, n_CGL.rhat)` built from explicit DOPC bead-target dry-MARTINI energies.
  - Rule: generic Martini pair lists should exclude `CGL` and `CGLD`; CGL-CGL, CGL-SC, and CGL-target interactions should be owned by dedicated spline nodes to avoid double-counting and hidden-site nonbonded artifacts.
  - Verification of the current bilayer-only probe: `1000` steps with the stage-6.0 NPT contract kept all lipids physically oriented (`bad_parallel=0`, `bad_flip=0`) and preserved CGL-CGLD length near the DOPC-derived `11.139 Å`.
- 2026-05-15: VTF lipid coloring needs endpoint-specific metadata, not only endpoint coordinates.
  - Current extracted VTF has `282` `LIPH/HYDROPHILIC/UNK` endpoints and `282` `LIPT/HYDROPHOBIC/DOPC` endpoints.
  - Existing stage metadata already stores DOPC-derived display offsets with mean endpoint span `24.270477 Å`; the visual length problem is metadata/coloring, not coordinate span.
  - Superseded by the 2026-05-16 user correction: do not use companion `.vmd` scripts or `ResType` mapping for this workflow.
- 2026-05-15: Regenerated 1RKL VTF endpoint coloring metadata is correct.
  - `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/1rkl.stage_7.0.vtf` now contains `282` `LIPH/HYDROPHILIC/LIPH/7` endpoints and `282` `LIPT/HYDROPHOBIC/LIPT/6` endpoints.
  - First-frame lipid endpoint bond span averages `24.270496 Å`, so the visible bar length remains at the DOPC-derived physical display length rather than growing into an artificial long same-color segment.
- 2026-05-16: User correction on VTF lipid coloring.
  - Rule: do not create companion `.vmd` files for this workflow; keep coloring support in the `.vtf` atom metadata only.
  - Rule: drop `ResType` support when it requires external VMD-side mapping.
  - A single LIPH-LIPT bond can still look like one same-colored rod in VMD bond renderings; use two same-colored half-bonds so `Name`, `Type`, `Element`, and `ResName` coloring can split hydrophilic and hydrophobic sides directly.
  - The full NC3-to-tail span from the starting PDB can look too long for the single-particle display rod; use the DOPC-derived orientation length for the visual rod span.
  - Regenerated 1RKL VTF verification: total visual span is `11.631147 Å`, split into two `5.8156 Å` side bonds, with no companion `.vmd` file.
- 2026-05-16: User correction: `LIPH` and `LIPT` still appeared single-colored under VMD `Name` coloring.
  - Rule: when relying on VMD default `Name` coloring without a companion script, use built-in atom names with known default categories.
  - Use `name N` for the hydrophilic side and `name C` for the hydrophobic side; keep `type=HYDROPHILIC/HYDROPHOBIC`, `resname=LIPH/LIPT`, and `atomicnumber=7/6` for semantic metadata.
- 2026-05-16: User correction: use the MARTINI coloring scheme for the single-particle lipid sides.
  - Rule: display atoms for CGL sides should reuse actual DOPC MARTINI bead categories, not generic element names.
  - Use `PO4/Qa` for the hydrophilic side and `C1A/C1` for the hydrophobic side while keeping side-specific `resname=LIPH/LIPT`.
  - Verified regenerated `1rkl.stage_6.0.vtf` and `1rkl.stage_7.0.vtf` both carry `564` `PO4/Qa/LIPH/15` display atoms and `564` `C1A/C1/LIPT/6` display atoms.
- 2026-05-15: Latest `example/16.MARTINI/outputs/martini_test_1rkl_hybrid` failure is a preproduction physics omission.
  - `stage_6.0.up` contains `cg_lipid_pair`, `compose_vector6d`, `dist_spring`, and `martini_potential`, but not `martini_sc_table_1body` or `cg_lipid_sc`.
  - `stage_7.0.up` contains `martini_sc_table_1body` and `cg_lipid_sc`, but the handoff from `6.0` already has a CGL-protein clash.
  - Trajectory inspection: `6.0` input has min CGL-protein distance `7.608 Å` and same-leaflet CGL minima near `6.959 Å`; `6.0` final has min CGL-protein distance `2.275 Å` and upper-leaflet same-leaflet minimum `0.588 Å`.
  - Rule: rigid preproduction may keep protein geometry fixed, but it must not omit active SC-env/BB-env/CGL-SC interface physics when relaxing lipids around the protein.
- 2026-05-15: Verification of the preproduction-interface fix shows the old same-leaflet XY diagnostic can look low when median-z leaflet assignment mixes tilted or crossing lipids; use 3D CGL-CGL nearest distance as the clash check.
  - Fixed reused-table `6.0` output: min CGL-protein distance `5.357 Å`, 3D CGL-CGL nearest distance `7.126 Å`, and CGL-CGL p05 nearest distance `7.675 Å`.
  - Fixed 10000-step production output: Rg range `12.8-14.1 Å`, final Rg `13.2 Å`, and `cg_lipid_sc` range `[-14.70, 18.88] E_up`.
- 2026-05-15: Audit of current 1RKL CG-lipid stabilizers found the largest non-ITP controls in the single-particle CGL representation.
  - DOPC-derived geometry from `parameters/dryMARTINI/DOPC.pdb` gives first-lipid head-to-tail span about `25.70 Å`, COM-to-tail projection about `11.14 Å`, max perpendicular bead radius about `4.11 Å`, and transverse rotational inertia about `9560 g/mol Å^2`.
  - Dry-MARTINI DOPC nonbonded sigmas in `dry_martini_v2.1.itp` are `0.47-0.62 nm`, so a physically derived CGL contact threshold is `2^(1/6) * 0.62 nm = 6.96 Å`, matching the old empirical `7.0 Å` spacing closely.
  - The old CGLD length `5.0 Å`, orientation spring `50 E_up`, CGL-vs-X sigma cap `0.9 nm`, and CGL same-leaflet spacing `7.0 Å` should be replaced or documented as derived from these force-field/geometry quantities.
- 2026-05-15: Implementation verification showed stage conversion derives CGLD defaults from the actual packed 1RKL lipid reference, producing `length=11.631 Å`, `mass=85.427 g/mol`, `bond_fc=36.314 E_up/Å²`, and same-leaflet spacing target `6.959 Å`.
- 2026-05-12: User clarification for the restored pre-production branch: preserve it for any real bonded dry-MARTINI environment particles, not just bonded lipids.
  - Rule: workflow routing must distinguish physical environment bonds from synthetic CG lipid orientation helpers such as `CGL-CGLD`.
  - Consequence: the always-on `6.0` stage covers the current CG-lipid path, while `6.1-6.6` remain available only when genuine bonded environment topology is present.
- 2026-05-07: The latest 1RKL protein collapse is dominated by overlong `cg_lipid_sc`, not by ion bonds or fixed lipid vectors.
  - Existing output has `cg_lipid_sc.fit_r_max_nm=0.7` but runtime `cutoff_ang=16.8`, extending the fitted sidechain-lipid attraction far beyond its sampled `7.0 Å` support.
  - Component diagnostics showed `cg_lipid_sc` growing from about `-86 E_up` at frame 0 to about `-2488 E_up` by frame 200, while protein Rg collapsed from about `12.8 Å` to `8.6 Å`.
  - Static and full-run tests with `cutoff_ang=8.4` (`7.0 Å` fit support plus `1.4 Å` taper) kept protein Rg near `12.6-12.8 Å` over 10000 steps while preserving active SC-env and BB-env interactions.
  - The production progress label `prot_potential` was misleading because it included `cg_lipid_pair` and `cg_lipid_sc`; those terms need separate reporting for future diagnostics.
- 2026-05-07: The strange first VTF frame in the latest 1RKL output had a visualization bug and a separate real dynamics bug.
  - HDF5 frame 0 CGL positions are initially slab-like and `dist_spring` has zero ion bonds; the apparent ion/protein attachments came from `py/martini_extract_vtf.py` using original HDF5 CGL indices after mode-2 output remapping.
  - Correct VTF extraction must remap original CGL indices to output atom indices before adding CGH direction-marker bonds; verified regenerated VTF has 282 `CGL-CGH` bonds and zero CGH bonds sourced from ions/protein.
  - The real launch instability was dominated by CGL-ion effective LJ: old CGL-ion sigma was about `17.9 Å`, the closest ion was `6.366 Å`, and the old generated `martini_potential` had CGL-ion LJ max `97404.794 E_up` / sum `145929.948 E_up`.
  - The correct fix is not to disable hybrid protein/environment physics; cap only the isotropic CGL-vs-non-CGL effective LJ radius and place ions with a CGL-aware cutoff. With a `0.9 nm` sigma cap, the same old geometry estimates only `22.167 E_up` max CGL-ion LJ.
  - Debug summaries should always include CGL-ion/protein nearest-distance and LJ metrics when diagnosing single-particle lipid failures.
- 2026-05-07: DOPC-sidechain was not using the same directional CG table method as single-particle DOPC-DOPC.
  - DOPC-DOPC used `cg_lipid_quadspline_v3`, full multimode params, and `eval_multimode_pair()`.
  - DOPC-sidechain still used fixed 54-parameter `cg_lipid_sc_quadspline_v1` and runtime `eval_quadspline()`.
  - Both paths used the same angular sign convention, `ang1=-n1_dot_n12;ang2=n2_dot_n12`, but the radial/mode schema and evaluator were different.
  - Resolution: use `cg_lipid_sc_quadspline_v2` with full multimode params and the same evaluator style as DOPC-DOPC, while retaining SC-specific rotamer averaging and lipid-against-fixed-SC relaxation.
- 2026-05-06: User correction: `run_sim_1rkl.sh` initial-structure debugging must not require the expensive MARTINI table or pair-list generation path.
  - Rule: use `INITIAL_DEBUG_ONLY=1` for visual initial-structure debugging. This runs real packing and CG conversion, writes PDB/TSV/JSON diagnostics, and exits before `martini.h5` generation/use, HDF5 nonbonded pair-list generation, and dynamics.
  - Normal fresh workflow runs must also emit the same single-particle lipid PDB immediately after packing, before table generation, so users do not need to opt into debug-only mode just to get the artifact.
  - Verified 1RKL debug output under `example/16.MARTINI/outputs/debug_1rkl_initial_only/debug/`: the bilayer-only PDB has 282 physical CGL particles, no CGLD, and no ions.
  - Verified normal `run_sim_1rkl.sh` output under `example/16.MARTINI/outputs/debug_1rkl_normal_generates_pdb/debug/`: `1rkl.stage_7.0.input_debug.bilayer.pdb` is generated before table building and has 282 physical CGL particles, no CGLD, and no ions.
  - The generated 1RKL topology report has 282 true bonds, all CGL-CGLD orientation bonds; `ion_bond_count=0`.
  - The generated 1RKL CGL leaflet stats after conditioning have same-leaflet nearest-neighbor min/p05 of about 7.0 Å, and CGL-protein minimum distance is about 7.60 Å.
- 2026-05-06: User correction: debug PDB export must include an explicit single-particle lipid bilayer PDB, not just all-particle and visible whole-system views.
  - Rule: when debugging the CG lipid bilayer, always emit a CGL-only PDB (`*.input_debug.bilayer.pdb`) that excludes protein, ions, water, and hidden CGLD orientation sites.
  - The bilayer-only test artifact now contains exactly the 72 physical CGL particles, making the single-particle lipid slab directly inspectable.
- 2026-05-06: Initial-structure debug export separates actual topology from visualization artifacts.
  - The current existing 1RKL stage-7 HDF5 input has `282` true `dist_spring` bonds and all are CGL-CGLD orientation bonds; `ion_bond_count=0`.
  - The ion-bond appearance in a viewer is therefore not coming from `/input/potential/dist_spring/id`; it is likely viewer bond inference or a stale visualization format.
  - The same existing 1RKL stage predates the CGLD molecule-id fix: all `282` CGLD sites have molecule IDs different from their parent CGL, which can make hidden orientation sites look like a separate chunk in molecule-based visualization.
  - Newly generated bilayer-only input has corrected CGLD molecule IDs (`72` molecules for `72` CGL/CGLD pairs) and no debug warnings.
- 2026-05-06: The current 1RKL single-particle lipid failure has two separable causes.
  - The inspected 1RKL stage has bad initial CGL placement: upper-leaflet CGL z extends below the midplane and same-leaflet XY spacing reaches about `2.245 Å`, while the trajectory explodes by the first saved frame.
  - A clean 72-DOPC bilayer-only system has balanced 36/36 leaflets, expected direction signs (`lower dir_z ~= +0.995`, `upper dir_z ~= -0.993`), and no close same-leaflet XY overlaps below about `5.35 Å`.
  - The CG lipid quadspline table was also over-amplified: radial fitting did not match Upside's clamped deBoor runtime basis, angular basis knots were incorrectly energy-converted despite being dimensionless, and refitting `V_angular(r)` by averaging `residual/(Ang1*Ang2)` amplified small angular products into huge coefficients.
  - Runtime-matched radial fitting plus least-squares projection/clipping of `V_angular(r)` reduced the bilayer-only `V_angular` max from order `10^3 E_up` to about `19.5 E_up` and removed the immediate numerical blow-up in the isolated test.
  - The isolated single-particle DOPC model still spreads in z over a longer NVT run, so numerical stability does not yet prove physical bilayer morphology stability.
- 2026-05-05: The current CG lipid bilayer instability comes from mismatched training/runtime semantics, not from a missing generic damping term.
  - `cg_lipid_pair` was fitting a full isotropic radial energy while the stage still retained CGL-CGL LJ in `martini_potential`, so the isotropic lipid-lipid term was represented twice.
  - The CG-CG table path enabled short relaxation even though its own docstring says relaxation is invalid for CG-CG fitting; this can remove the repulsive wall that should remain in the effective interaction.
  - CG-CG angular samples were generated relative to the global z-axis while the runtime quadspline interprets angles relative to the pair displacement vector.
  - The generic `InteractionGraph<PosQuadSplineInteraction>` runtime path uses direct coordinate differences, so CG lipid pair and SC interactions miss minimum-image PBC even though the base MARTINI potential is PBC-aware.
  - The shared quadspline interaction signs do not directly encode the requested lipid convention `Ang1(-n1*n12) Ang2(n2*n12)`, so the CG lipid runtime should own this convention explicitly instead of changing shared bead interaction behavior.
  - Direct explicit-bead CG-CG sampling contains very large overlap outliers even outside the nominal CGL-CGL core. If those outliers are fit into the angular residual, they produce enormous spline coefficients and destabilize the one-particle lipid model; the isotropic CGL-CGL LJ core should own excluded volume while the angular residual remains bounded.
- 2026-04-30: User correction: AA-direct hybrid mapping must not silently alias missing BB proxies to CA coordinates.
  - If `bb_atom_index=-1` is a sentinel, runtime must keep it as missing-proxy state and derive BB interaction sites from mapped carriers.
  - Writing BB COM coordinates into physical CA coordinates changes the protein model and breaks helicity.
- 2026-04-30: AA-direct rigid-body behavior requires geometric projection, not just force/momentum filtering.
  - A rigid-body force redistribution alone allows internal drift under finite-step integration.
  - The stage-6 implementation must project coordinates back to a rigid manifold each step while preserving net translation/rotation.
- 2026-04-30: In AA-direct mode, proxy-fix routing must exclude physical carrier atoms.
  - Logic that fixed `ROLE_BB` atoms was valid only when those atoms were virtual MARTINI proxies.
  - With carrier-backed BB sites (`CA` role), that same logic incorrectly freezes real protein coordinates and destabilizes backbone geometry.
- 2026-04-30: BB mapping parity with martinize includes termini fragment overrides.
  - Extracted BB typing is incomplete unless charged fragment termini are explicitly set to `Qd`/`Qa` and charges follow BB type.
  - Missing this parity shifts interface electrostatics and can destabilize early-stage behavior.
- 2026-04-29: User correction: the previous production restart patch is incomplete.
  - Do not assume recording/copying momentum plus a transition counter is sufficient.
  - Re-check whether the restart source is the true final integrator state or only the last logged frame.
  - Continuation should not silently proceed from incomplete restart state.
- 2026-04-29: Production restart root cause after re-audit.
  - A restart that reuses saved momentum must not refresh hybrid reference carrier positions; doing so creates a position/momentum mismatch unique to split production runs.
  - Normal MD output must include the final integrator state. If the last saved frame is only the last periodic pre-integration sample, continuation advances the hybrid transition counter from a stale time and restarts from stale coordinates.
  - Older production files without `output/mom`, or with momentum but no validated final-state marker, are not exact restart sources; the workflow should fail before MD rather than rethermalize silently or use a stale frame.
- 2026-04-28: User correction: production restart parity must be verified explicitly.
  - Stage-7 continuation may look syntactically valid while missing runtime state such as momentum, box, or hybrid ramp/hold counters.
  - Verification must compare split production (`7.0 -> 7.1`) against an uninterrupted run of matching length, not just check that continuation executes.
- 2026-04-28: Production continuation was missing two state channels.
  - Without `--record-momentum`, stage files do not contain `output/mom`, so continuation rethermalizes instead of preserving the velocity state.
  - Without a persisted SC-env transition counter, every `7.x` restart begins the force-cap/backbone-hold schedule again from step zero.
  - Correct restart handoff needs to copy saved momentum into `input/mom`, pass `--restart-using-momentum`, and set `hybrid_control.sc_env_transition_step_start` from saved production time.
- 2026-04-28: User correction after the compact-reference fix: a clean initial stage file is not enough; stage-7 must be replayed for at least the early production window because the current run explodes by step `50`.
  - Step-0 energies can look plausible while force propagation is still catastrophically wrong.
  - Verification for this workflow must include a short production replay with Rg/energy inspection, not only HDF5 structure checks and one-frame extraction.
  - The fix must preserve active SC-env and BB-env interactions rather than hiding the explosion behind disabled physics.
- 2026-04-28: The compact runtime AA reference mapping requires C++ runtime consumers to use `hybrid_bb_map/atom_indices` directly.
  - Recomputing `reference_index_offset + reference_atom_indices` is invalid once raw AA PDB indices are compacted.
  - This specifically corrupted the backbone O refresh path and can destabilize production immediately.
- 2026-04-28: The active SC/environment table node needs the migrated SC force cap too.
  - `martini_potential` already capped startup pair forces, but `martini_sc_table_1body` was applying uncapped table gradients to environment atoms and CB feedback.
  - Applying `sc_env_lj_force_cap` to SC-table point/vector gradients preserves the SC-env interaction while preventing launch-force events.
- 2026-04-28: User correction after the hybrid workflow refactor: syntax and reduced-run smoke tests are insufficient when the first VTF frame depends on exact hybrid runtime atom indexing.
  - The Python stage injector must not preserve sparse AA PDB atom-index gaps as appended runtime particles.
  - When converting AA backbone references into stage-runtime carriers, append a compact set of the actually referenced N/CA/C/O atoms and map raw reference indices through that compact lookup.
  - This prevents padded inert reference rows from changing the production stage geometry surface while keeping SC-env and BB-env interface potentials active.
- 2026-04-28: `example/16.MARTINI/run_sim_1rkl.sh` must source the repo bootstrap before enabling `set -u`.
  - `source.sh` can reference `MY_PYTHON` before it is set in the caller environment;
  - enabling `set -u` before sourcing it aborts the workflow immediately;
  - the launcher should establish the project environment first, then enable strict shell mode.
- 2026-04-24: The temp checkout’s dry-MARTINI speedup is concentrated in three hot-path changes rather than in its unrelated feature diffs.
  - `MartiniPotential` adds a cached Verlet-style pairlist (`cache_buffer`, cached coordinates/box, active pair indices) so the full `pairs` table is not rescanned every force evaluation.
  - It also compacts coefficient rows into a unique parameter table and stores direct spline pointers per unique parameter row, removing repeated per-pair coefficient unpacking and spline-map lookup from the inner loop.
  - The active SC/environment runtime in this checkout is `martini_sc_table_1body`, and the temp checkout accelerates that path with a cached active-contact list over `(row, env)` contacts inside `cutoff + cache_buffer`.
- 2026-04-24: The active production SC/environment node in this repo is `martini_sc_table_1body`, not `martini_sc_table_potential`.
  - `py/martini_prepare_system_lib.py` injects `martini_sc_table_1body`;
  - therefore the temp-file active-contact cache belongs on `martini_sc_table_1body` in this checkout, while `martini_sc_table_potential` should remain as a compatibility surface.
- 2026-04-24: The checked-in `obj/` build tree is stale after the repo move.
  - `obj/CMakeCache.txt` still references `/Users/yinhan/Documents/upside2-md-martini`;
  - reliable verification currently requires a fresh configure/build outside that stale tree unless the user wants `obj/` reinitialized.
- 2026-04-14: The reported `ModuleNotFoundError: h5py` from `example/16.MARTINI/run_sim_1rkl.sh` was caused by a moved virtualenv, not by a missing `h5py` wheel in the repo `.venv`.
  - `./.venv/bin/python` imports `h5py` successfully;
  - the failing interpreter was Homebrew `python3.14`, reached because `.venv/bin/activate` still hardcoded `VIRTUAL_ENV=/Users/yinhan/Documents/upside2-md-martini/.venv`;
  - stale absolute paths also remained in console-script shebangs under `.venv/bin` such as `pip`, so the robust installer fix needs to repair both activation scripts and launcher shebangs in place.
- 2026-04-14: The current Python installer is Linux-biased in two ways that matter on local macOS.
  - it defaults to `PYTHON_BIN=python3.11`, but the host default `python3` on this Darwin arm64 machine is `/opt/homebrew/opt/python@3.14/bin/python3.14`, so the script should discover and validate a Python `3.11` interpreter instead of assuming one fixed executable path;
  - it installs required and optional packages in one `pip install` invocation under `set -e`, so any macOS-specific failure in optional extras aborts the full environment bootstrap.
- 2026-04-14: The current Darwin arm64 machine already exposes the Homebrew-native libraries most relevant to HDF5/PyTables source builds.
  - `brew --prefix` resolves successfully for `hdf5`, `lzo`, `bzip2`, and `c-blosc2`;
  - these prefixes can be exported as build hints without affecting Linux behavior.
- 2026-04-14: The repository-local virtual environment can remain reusable even when its `activate` script becomes stale after a repo move.
  - on this machine, `.venv/bin/activate` still hardcodes `VIRTUAL_ENV=/Users/yinhan/Documents/upside2-md-martini/.venv`, which breaks `source .venv/bin/activate` even though `.venv/bin/python` itself is valid and executable;
  - consequence:
    - the installer should drive the environment with `$VENV_DIR/bin/python` directly instead of depending on activation side effects.
- 2026-04-13: The current `hybrid-interface-sweep` stage-7 outputs require schema-aware analysis rather than assumptions copied from the bilayer workflow.
  - Fresh reduced `1rkl.stage_7.0.up` output from `/tmp/hybrid_interface_sweep_analysis_smoke/tasks/scale0p85_r01/run/checkpoints/1rkl.stage_7.0.up` shows:
    - `particle_class` does not provide a usable `LIPID` partition for this analyzer,
    - `/output/box` may be absent,
    - protein carrier selection should come from `input/hybrid_bb_map`, not proxy MARTINI atom classes.
  - The implemented analyzer therefore:
    - falls back to `input/potential/martini_potential.{x_len,y_len,z_len}` when `/output/box` is absent,
    - identifies lipid molecules from non-protein `PO4` beads plus `molecule_ids`,
    - identifies protein `CA` carrier atoms from `hybrid_bb_map.atom_indices[:,1]` where `atom_mask[:,1] != 0`.
- 2026-04-13: The new hybrid sweep analysis uses diffusion-style observables as the practical fluidity proxy.
  - Per completed task, it computes:
    - protein lateral COM diffusion relative to the bilayer COM,
    - lipid `PO4` lateral diffusion relative to the bilayer COM as a guardrail.
  - The implementation unwraps XY coordinates, removes bilayer COM drift, drops the first `20%` of frames as burn-in, and fits `MSD = 4Dt + b` over a fixed lag window.
- 2026-04-02: The active stage-7 SC/dry-MARTINI simulation path is not probabilistic and does not use a separate sidechain face-vector DOF.
  - `example/16.MARTINI/inject_sc_table_stage7.py` deletes legacy `rotamer`, `placement_fixed_scalar`, and `placement_fixed_point_vector_only` nodes, then injects only:
    - `affine_alignment`,
    - `placement_fixed_point_only_CB`,
    - `martini_sc_table_potential`;
  - the injected `martini_sc_table_potential` takes arguments `[pos, placement_fixed_point_only_CB]` and stores only:
    - `cb_index`,
    - `residue_table_index`,
    - `env_atom_index`,
    - `env_target_index`,
    - `grid_nm`,
    - `energy_kj_mol`,
    - residue/target order metadata;
  - `src/martini.cpp` then evaluates that node as a purely radial CB-to-environment table over `dist = |CB - env|`, with no rotamer probabilities, no angular term, and no sidechain orientation/face vector input;
  - the current generated stage artifact `/tmp/legacy_cleanup_short/checkpoints/1rkl.stage_7.0.prepared.up` confirms the active potential graph contains only `affine_alignment`, `placement_fixed_point_only_CB`, and `martini_sc_table_potential` among the relevant SC nodes.
- 2026-04-02: Upside still has a general vector-aware CB placement path, but the current stage-7 MARTINI integration does not use it.
  - `py/upside_config.py` contains a `placement_fixed_point_vector_only_CB` path with explicit stored direction vectors;
  - the active stage-7 injector instead uses `placement_fixed_point_only_CB`, so that vector/orientation capability is currently bypassed for the dry-MARTINI SC coupling.
- 2026-04-02: The active protein-CG generation path uses fresh `martini22` backbone bead assignments with charged termini, but the checked-in fallback `example/16.MARTINI/pdb/1rkl_proa.itp` is stale relative to the current `martinize.py`.
  - `example/16.MARTINI/run_sim_1rkl.sh` defaults `MARTINIZE_ENABLE=1` and regenerates the protein CG PDB/ITP through `martinize.py` before hybrid preparation;
  - `example/16.MARTINI/martinize.py` sets charged termini by default because `-nt` defaults to `False`, logs that chain termini will be charged, and explicitly overrides the first and last backbone bead types to `Qd` / `Qa`;
  - fresh generated output from the active workflow at `/tmp/legacy_cleanup_short/hybrid_prep/1rkl_hybrid_proa.itp` confirms:
    - first backbone bead: `Qd`, charge `+1.0000`,
    - last backbone bead: `Qa`, charge `-1.0000`,
    - representative internal backbone bead: residue 2 (`GLU`) uses neutral `P5`;
  - the checked-in fallback `example/16.MARTINI/pdb/1rkl_proa.itp` does not match that current contract:
    - first backbone bead: `Q5`, charge `+1`,
    - last backbone bead: `Q5`, charge `-1`,
    - representative internal backbone bead: residue 5 (`GLU`) uses `P2`;
  - consequence:
    - the active default workflow is correct as long as `MARTINIZE_ENABLE=1`,
    - the checked-in fallback ITP should not be treated as equivalent to fresh current `martini22` output.
- 2026-04-02: The remaining legacy surface after the direct-Upside rewrite was prep-only, not runtime.
  - `example/16.MARTINI/prepare_system_lib.py` was still exporting `hybrid_sc_map` and `example/16.MARTINI/validate_hybrid_mapping.py` was still requiring it, even though the active stage-7 workflow no longer consumed that group;
  - after removing that dead path, a fresh prep artifact in `/tmp/hybrid_mapping_bb_only/hybrid_mapping.h5` validates successfully and now contains only:
    - `/input/hybrid_bb_map`,
    - `/input/hybrid_control`,
    - `/input/hybrid_env_topology`.
- 2026-04-02: `src/box.cpp` and `src/box.h` no longer carry the retired hybrid rotamer/proxy subsystem.
  - the remaining code there is active NPT box propagation and compatibility plumbing;
  - no stage-7 legacy hybrid runtime behavior remained to remove from those files in this pass.
- 2026-04-02: The current production energy-drift regression came from two active stage-7 paths that bypassed the intended direct-Upside semantics.
  - direct `BB -> CA` MARTINI forces in `src/martini.cpp` were applying full protein-side feedback from step `0`, bypassing the existing `sc_env_backbone_hold_steps` ramp that older hybrid projection paths used;
  - `martini_sc_table_potential` likewise fed full `CB` force back into Upside from step `0`;
  - stage-7 still allowed protein-internal MARTINI proxy-proxy terms through `allow_protein_pair_by_rule(...)`, which kept legacy proxy bonded/nonbonded bookkeeping energy alive even though those protein proxies are no longer the physical protein model in direct-Upside production.
- 2026-04-02: Reinstating the startup protein-feedback ramp on the direct BB/CB paths helps substantially but is not sufficient by itself.
  - replay from the saved production handoff with only the direct-feedback-ramp restoration reduced the archived step-`2050` drift from `martini_potential 26387.78 / total 26216.96` to about `10259.78 / 10097.96`;
  - the remaining positive drift was then traced to leftover stage-7 protein proxy-proxy terms.
- 2026-04-02: Removing the remaining active protein-internal MARTINI proxy-proxy terms resolves the reproduced production drift on the saved handoff.
  - replay from `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up` with both fixes reached step `2050` at:
    - `martini_potential -25985.45`,
    - `total -26147.28`;
  - the archived current run at the same step is:
    - `martini_potential 26387.78`,
    - `total 26216.96`;
  - the patched replay remained negative through all checked saved lines `0 -> 2050`, so the reproduced long-horizon energy build-up no longer occurs in that stage-7 artifact.
- 2026-04-02: Reduced simulation masses are now restored as the active workflow state.
  - `example/16.MARTINI/prepare_system_lib.py` again writes `ff_mass / 12` into `/input/mass`;
  - `example/16.MARTINI/run_sim_1rkl.sh` and `example/16.MARTINI/test_prod_run_sim_1rkl.sh` again append reduced `/12` reference masses;
  - fresh shortened workflow output `/tmp/martini_mass_reduced/checkpoints/1rkl.stage_6.0.up` confirms the reduced values are present:
    - `4082` particles at `6.0`,
    - `8` particles at `3.75`;
  - the same run again starts stage `6.0` at `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`.
- 2026-04-02: The current huge `martini_potential` / `total` values in `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/logs/stage_7.0.log` are real runtime values, not a logging bug.
  - `src/main.cpp` prints the direct result of `split_engine_potential_terms(...)` and `sys.engine.potential`;
  - the active stage-7 log shows `martini_potential` on the order of `1.9e21` at step `0`, decaying to about `1.27e21` by step `1200`;
  - consequence:
    - formatting can be improved for readability,
    - but a print-path edit does not address the physical/force-field cause of those giant energies.
- 2026-04-02: Restoring native dry-MARTINI masses alone does not reproduce the earlier minimization blow-up.
  - the weighted 4-carrier BB mapping remained restored while only the mass path was switched back to native force-field values;
  - a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run in `/tmp/martini_mass_restore` again started stage `6.0` at:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`;
  - direct HDF5 inspection of `/tmp/martini_mass_restore/checkpoints/1rkl.stage_6.0.up` confirms native dry-MARTINI mass values are present in `/input/mass`:
    - `4082` particles at `72.0`,
    - `8` particles at `45.0`.
- 2026-04-01: The last two attempted cleanup changes were not simulation-safe and have been reverted.
  - the attempted switch to native dry-MARTINI masses in `/input/mass` caused immediate stage-6.0 minimization blow-up in the real workflow;
  - the attempted CA-only active `hybrid_bb_map` remap was also reverted together with that mass change, per user request;
  - after reverting both, a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run returned to the pre-regression startup:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`.
- 2026-04-01: Dry-MARTINI particle masses now stay in native dry-MARTINI profile units throughout the simulation workflow.
  - `example/16.MARTINI/prepare_system_lib.py` no longer divides force-field masses by `12` when writing `/input/mass`;
  - `example/16.MARTINI/run_sim_1rkl.sh` and `example/16.MARTINI/test_prod_run_sim_1rkl.sh` no longer append reduced-unit AA reference masses during stage-file augmentation;
  - verified on fresh shortened workflow outputs:
    - `/tmp/martini_mass_native/checkpoints/1rkl.stage_6.0.up`,
    - `/tmp/martini_mass_native/checkpoints/1rkl.stage_7.0.prepared.up`,
    - `/tmp/martini_mass_native/checkpoints/1rkl.stage_7.0.up`.
  - for the physical dry-MARTINI particles in the current 1RKL system, `/input/mass` now contains native force-field values only:
    - `4082` particles at `72.0`,
    - `8` particles at `45.0`.
  - the remaining appended AA reference rows carry unconverted bookkeeping masses (`18.6667`, `16.0`, `21.3333`, plus `1.0` placeholders for unused reference slots), not reduced `/12` values.
- 2026-04-01: The old weighted `BB -> N/CA/C/O` active carrier map is now fully removed from the live hybrid force path.
  - `example/16.MARTINI/prepare_system_lib.py` now exports CA-only active BB rows (`atom_mask=[0,1,0,0]`, `weights=[0,1,0,0]`) while keeping separate `reference_atom_indices` / `reference_atom_coords` for stage-7 reference geometry;
  - `example/16.MARTINI/run_sim_1rkl.sh` and `example/16.MARTINI/test_prod_run_sim_1rkl.sh` now inject CA-only runtime BB rows into every prepared stage file, rather than remapping all four backbone reference particles into active carriers;
  - `src/martini.cpp` now rejects any active `N/C/O` BB target rows at load time and uses separate stored reference-runtime indices for stage-7 `O` reconstruction, so reference geometry no longer implies four-way active force projection;
  - verified on fresh generated files:
    - `/tmp/martini_ca_only/hybrid_prep/hybrid_mapping.h5`,
    - `/tmp/martini_ca_only/checkpoints/1rkl.stage_6.0.up`,
    - `/tmp/martini_ca_only/checkpoints/1rkl.stage_7.0.prepared.up`,
    - `/tmp/martini_ca_only/checkpoints/1rkl.stage_7.0.up`,
    all show only `atom_mask=[0,1,0,0]` and `weights=[0,1,0,0]` for `/input/hybrid_bb_map`.
- 2026-04-01: The direct-Upside stage-7 architecture is now verified end to end in the real workflow.
  - `src/martini.cpp` active production now skips protein `SC` proxy MARTINI nonbonded pairs and treats the Upside `CA` carrier as the dry-MARTINI `BB` interaction site for protein-backbone/environment pairs;
  - `example/16.MARTINI/run_sim_1rkl.sh` and `example/16.MARTINI/test_prod_run_sim_1rkl.sh` no longer inject `hybrid_sc_map` into runtime stage files, so active stage logs now report `n_sc=0` and the legacy probabilistic SC proxy path is inactive;
  - direct stage-file inspection on the shortened full workflow confirms:
    - `6.6` has no `hybrid_sc_map`, no `martini_sc_table_potential`, and no `placement_fixed_point_only_CB`;
    - `7.0.prepared` has `martini_sc_table_potential` and `placement_fixed_point_only_CB`, with no `hybrid_sc_map` or `rotamer`.
- 2026-04-01: Fresh production-stage prepared files from the current `prepare_system.py` path may omit `/input/sequence`.
  - the stage-7 SC-table injector therefore cannot assume `sequence` exists in the prepared `.up` file;
  - the correct fallback source is the effective protein ITP passed by the workflow, using residue names from `[ atoms ]` rows with `BB` roles;
  - after adding that fallback, the injector now also writes the normalized sequence back into the stage file for downstream consistency.
- 2026-04-01: Stage gating for the new SC table is now verified against real generated stage files.
  - `example/16.MARTINI/run_sim_1rkl.sh` only injects the new SC table when `stage_label=production`;
  - observed generated files confirm the intended split:
    - `6.6` has no `martini_sc_table_potential`, no `placement_fixed_point_only_CB`, and still retains `/input/hybrid_sc_map`;
    - `7.0.prepared` has `martini_sc_table_potential`, `placement_fixed_point_only_CB`, no legacy `rotamer`, and no `/input/hybrid_sc_map`.
- 2026-04-01: With the repaired repository `.venv`, both the focused production helper and the real top-level workflow now execute the new stage-7 SC-table path successfully.
  - `example/16.MARTINI/test_prod_run_sim_1rkl.sh` completed a fresh stage-7 preparation/injection plus a short `5`-step production replay from the saved `6.6` handoff;
  - `example/16.MARTINI/run_sim_1rkl.sh` also reached stage 7 and completed a short production replay under reduced verification settings, proving that the actual workflow consumes the new stage-7 injector path.
- 2026-04-01: Shortened 6.x relaxation remains a poor physical benchmark even though the stage-7 wiring is correct.
  - a reduced full-workflow verification run with `20`-step stages reached stage 7 successfully, but stages `6.4-6.6` still developed enormous MARTINI energies under those intentionally shortened pre-production settings;
  - consequence:
    - use the short run only as an integration/smoke check for stage gating and file wiring;
    - use the intended relaxation horizon for any stability or thermodynamic assessment of the new force field.
- 2026-04-01: Completed external training results are available at `/Users/yinhan/Documents/SC-training/runs/default/results/assembled/`.
  - assembled summary confirms `18` non-empty residue sidechain types, `38` dry-MARTINI target particle types, and `684` residue-target tasks;
  - assembled `sc_table.json` provides a uniform `96`-point radial grid per residue-target table and is sufficient to build a compact native-unit `martini.h5` without rerunning training.
- 2026-04-01: The correct stage-7 replacement scope is narrower than “remove the old augmenter”.
  - removing `rotamer` / `placement_fixed_scalar` / `placement_fixed_point_vector_only` is correct for the new direct SC table path;
  - removing the full stage-7 augmenter would also remove the existing Upside backbone force nodes (`Distance3D`, `Angle`, `Dihedral_*`, rama, hbond, backbone pairs), which the user explicitly still wants summed with the new SC/dry-MARTINI force.
  - consequence:
    - stage 7 must still inject the Upside backbone nodes;
    - the replacement only swaps the sidechain-specific production path for `affine_alignment`, `placement_fixed_point_only_CB`, and the new `martini_sc_table_potential`.
- 2026-04-01: Local environment limitation for workflow verification remains active.
  - default `/opt/homebrew/opt/python@3.14/bin/python3.14` still lacks both `h5py` and `tables`;
  - the repository `.venv` points at a missing Homebrew Python 3.10 path in this environment, so the actual HDF5-mutating stage-7 injector and full `run_sim_1rkl.sh` execution cannot be run end to end here without repairing the Python environment first.
  - workaround used in this pass:
    - `martini.h5` generation is implemented through `python3 + h5import`, which works locally;
    - injector/runtime integration was validated through syntax checks, file-structure checks, and a full C++ rebuild.
- 2026-03-31: User clarified the intended dry-MARTINI integration target for the new sidechain work:
  - drop sidechain back-mapping as a design goal for this force-field effort;
  - replace it with a direct dry-MARTINI sidechain-type to dry-MARTINI particle-type interaction table integrated into the current Upside hybrid framework.
  - evaluation of `plan.md` should therefore focus on whether it provides a clean table-driven replacement for the existing SC back-mapping/coupling path, not on preserving that older path.
- 2026-03-31: User further clarified that the replacement must also drop the current sidechain-relaxation workflow path in `example/16.MARTINI/run_sim_1rkl.sh`, and that benchmarking should be run on that real workflow.
  - current workflow evidence:
    - `example/16.MARTINI/run_sim_1rkl.sh` still defines and writes production SC relaxation/control attrs under `/input/hybrid_control` (`sc_env_lj_force_cap`, `sc_env_coul_force_cap`, `sc_env_relax_steps`, `sc_env_backbone_hold_steps`, `sc_env_po4_z_hold_steps`, `sc_env_relax_dt`, `sc_env_restraint_k`, `sc_env_max_displacement`, PO4 clamp controls, SC energy dump controls);
    - the same script still injects stage-7 `rotamer`, `placement_fixed_scalar`, `placement_fixed_point_vector_only`, and `affine_alignment` nodes through `augment_production_rotamer_nodes(...)`;
    - `example/16.MARTINI/test_prod_run_sim_1rkl.sh` mirrors the same SC-relaxation control surface, so it is only a smoke helper for the current design, not the authoritative benchmark target for the replacement architecture.
  - consequence for plan review:
    - if the new force field removes SC back-mapping/relaxation entirely, both the runtime SC path and the workflow-side stage-7 augmentation/control wiring must be removed or bypassed together;
    - performance/stability benchmarking should use the full `run_sim_1rkl.sh` workflow (real handoff, real stage-7 preparation, real horizon), not only short helper probes.
- 2026-03-31: New `SC-training/` workflow implementation findings:
  - `example/16.MARTINI/martinize.py` can be imported safely as a Python module for forcefield metadata without triggering CLI execution; the `martini22` class provides reusable canonical residue sidechain bead definitions.
  - `example/16.MARTINI/ff_dry/dry_martini_v2.1.itp` provides a complete `atomtypes` + `nonbond_params` source sufficient to generate first-pass residue-target training tasks directly from repository data.
  - the implemented first-pass workflow currently produces `18` non-empty canonical residue sidechain types against `38` dry-MARTINI target particle types, for `684` residue-target tasks total.
  - the workflow now samples target positions over spherical shells around the sidechain center, but still assigns isotropic shell energies because explicit residue-sidechain geometry reconstruction is not yet part of the baseline model (`sum_beadwise_colocated_spherical_shells`).
  - benchmark orchestration now targets `example/16.MARTINI/run_sim_1rkl.sh`, but actual consumption of the assembled SC table still requires separate runtime integration in `src/martini.cpp`.
- 2026-03-31: User corrected the default SC-training target scope back to the full bundled dry-MARTINI particle-type list:
  - default manifest generation should include all `38` dry-MARTINI particle types present in `SC-training/data/dry_martini_v2.1.itp`, including the `S*` and `AC*` types, not just a filtered non-ring subset;
  - for the current sidechain definitions that means the complete table remains `18 x 38 = 684` residue-target tasks.
- 2026-03-31: Slurm parallelization check for `SC-training/`:
  - the actual training execution path is a Slurm array job, not a serial loop:
    - generated `train_array.sbatch` uses `#SBATCH --array=0-683` for the current manifest;
    - each array element executes exactly one task through `workflow.py run-array-task --round-manifest ... --task-id "$SLURM_ARRAY_TASK_ID"`;
    - result assembly is staged as a separate dependent collector job/script via `assemble-results`.
  - verification also exposed a launcher bug:
    - both `SC-training/run_local.sh` and `SC-training/submit_remote_round.sh` sourced repo-root `source.sh` under `set -u`, which could abort immediately when `MY_PYTHON` was unset in the parent environment;
    - the wrappers now seed `MY_PYTHON` from the chosen interpreter path when needed and temporarily relax `nounset` while sourcing `source.sh`.
- 2026-03-31: User-reported Midway launcher failure exposed a second Slurm wrapper bug:
  - when `submit_remote_round.sh` itself is executed under `sbatch`, `dirname "$0"` can resolve to Slurm's spool copy under `/var/spool/slurm/...` rather than the real `SC-training` directory;
  - that made the default `BASE_DIR` fall back under the spool tree (`.../runs/default`), which fails with permission errors and also breaks relative lookup of `workflow.py`;
  - the wrappers now resolve `SCRIPT_DIR` by searching for a directory that actually contains `workflow.py`, preferring the real submit location from `SLURM_SUBMIT_DIR` and `SLURM_SUBMIT_DIR/SC-training` when present;
  - simulated spool-copy execution verified both common submission styles:
    - `sbatch submit_remote_round.sh` from inside `SC-training/`;
    - `sbatch SC-training/submit_remote_round.sh` from the repository root.
- 2026-03-31: Double-check of the current MARTINI runtime against the new SC-table requirements:
  - force symmetry / two-way force is already the relevant runtime pattern:
    - direct pair evaluation uses `gi = -force`, `gj = force` before accumulation;
    - the probabilistic SC-env path computes equal-and-opposite `proxy_grad` and `env_grad`, applies the environment contribution directly to the dry/environment atom, and later projects the SC-side gradient back through `hybrid_sc_map`;
    - `project_sc_gradient_if_active(...)` and `project_bb_gradient_if_active(...)` both use additive gradient updates, so MARTINI feedback is summed with existing Upside forces rather than replacing them.
  - current protein-pair policy in `src/martini.cpp` excludes `BB-BB` and `SC-SC` MARTINI nonbonded interactions and only allows `BB-SC` for the same residue.
  - consequence for the replacement architecture:
    - the new SC/dry table should target non-protein dry-MARTINI particles and project its feedback onto protein carriers;
    - backbone dry-`BB`/`CA` to surrounding dry-particle coupling should remain a separate backbone/environment path;
    - direct protein internal SC-backbone interactions should not be reintroduced through the new table.
- 2026-03-31: Dry-MARTINI unit-contract finding for training vs simulation:
  - user corrected the desired design after that review:
    - training should stay in native dry-MARTINI units only;
    - simulation should use explicit unit-conversion parameters instead of hardcoded conversion numbers.
  - current corrected state:
    - `SC-training/workflow.py` records native-unit policy only and no longer emits baked Upside conversion values;
    - `example/16.MARTINI/prepare_system_lib.py` now requires explicit simulation env vars `UPSIDE_MARTINI_ENERGY_CONVERSION` and `UPSIDE_MARTINI_LENGTH_CONVERSION`, then writes the corresponding attrs (`energy_conversion_kj_per_eup`, `length_conversion_angstrom_per_nm`, `coulomb_constant_native_kj_mol_nm_e2`) into `martini_potential`;
    - `example/16.MARTINI/run_sim_1rkl.sh` and `example/16.MARTINI/test_prod_run_sim_1rkl.sh` now fail early if those required conversion parameters are not provided;
    - `src/martini.cpp` now derives runtime Coulomb scaling from those attrs instead of hardcoding `31.775347952181`.
- 2026-03-31: Portability check for uploading only `SC-training/`:
  - before this pass, the workflow was not actually self-contained:
    - defaults still pointed at `example/16.MARTINI/martinize.py` and `example/16.MARTINI/ff_dry/dry_martini_v2.1.itp`;
    - `run_local.sh`, `submit_remote_round.sh`, and generated Slurm scripts still required repo-root `.venv` and `source.sh`.
  - current state after the portability patch:
    - `SC-training/data/dry_martini_v2.1.itp` now bundles the dry-MARTINI nonbond parameter source used by training;
    - `SC-training/data/martini22_sidechains.json` now bundles the martini22 canonical residue sidechain bead definitions needed for task generation;
    - `workflow.py` now defaults to those bundled data files and can load sidechain definitions from JSON directly;
    - local and Slurm launchers now work without parent-repo environment files for training, using `python3` by default and optional local `.venv` activation if present.
  - residual limitation:
    - the optional benchmark hook still points to `example/16.MARTINI/run_sim_1rkl.sh`, so uploading only `SC-training/` is sufficient for training and Slurm staging but not for benchmark execution.
- Delayed production instability is reproducible from the saved baseline stage-7 log:
  - file: `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/logs/stage_7.0.log`
  - reported values match the user's failure window exactly:
    - `5000 / 100000 ... potential   595.84, martini_potential -26362.14, total -25766.30`
    - `10000 / 100000 ... potential 14036.12, martini_potential   516.82, total 14552.94`
  - the next root-cause pass must therefore focus on a delayed instability in the hybrid runtime, not only on the already-fixed carrier-export corruption.
- Saved baseline production artifact:
  - file: `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`
  - stage-7 output contains actual protein AA carrier atoms with roles `N/CA/C/O` on runtime indices `4090..4322`.
- Hybrid virtual-site runtime defect:
  - active hybrid `BB` proxies are refreshed from the live AA carrier map in `src/martini.cpp::refresh_bb_positions_if_active(...)`, so they are virtual coordinates rather than independent dynamic particles during stage-7 force evaluation.
  - before the current fix, those same `BB` proxies were still being thermalized and integrated as ordinary atoms:
    - a 20-step recorded-momentum replay from the exact `6.6 -> 7.0` workflow handoff showed frame-0 nonzero momentum on all 31 active `BB` proxies,
    - frame-0 `BB` momentum norm mean/max were about `4.55e-2 / 8.74e-2`,
    - frame-0 `BB` momentum also disagreed with the carrier-derived mapped momentum (mean/max mismatch about `4.80e-2 / 9.90e-2`).
  - this means the thermostat was injecting ghost kinetic energy into coordinates that the hybrid runtime immediately overwrites before force evaluation, which is a non-Hamiltonian energy source and a plausible delayed 5k -> 10k instability mechanism.
- Actual carrier geometry from that artifact:
  - restrained `N-CA` and `CA-C` bonds remain near expected values over saved frames.
  - `C-O` is not part of the injected `Distance3D` / `Spring_bond` backbone restraint set and drifts badly in the saved output.
  - measured `C-O` carrier distance changes from about `1.234 Å` at frame 0 to a mean of about `17.90 Å` and a max of about `30.30 Å` at the last saved frame.
- Source-path finding:
  - `example/16.MARTINI/run_sim_1rkl.sh` injects four active BB carriers with weights `[N, CA, C, O]`.
  - `augment_production_rotamer_nodes()` builds production backbone nodes using only `N/CA/C` (`ref_n_atom = 3 * len(residue_ids)`), so `O` receives no direct backbone restraints.
  - `src/martini.cpp::refresh_bb_positions_if_active(...)` still uses the four-carrier weighted COM, so drifting `O` contaminates the refreshed MARTINI `BB` proxy.
- Export-path finding:
  - `example/16.MARTINI/extract_martini_vtf.py::reconstruct_backbone_aa(...)` does not use the actual runtime AA carriers when they exist.
  - It backmaps each residue by adding a pure translation from the live `BB` coordinate to the stored reference `N/CA/C/O` geometry, with no residue rotation.
  - In the saved stage-7 artifact, that reconstructed backbone diverges strongly from the actual runtime carriers by the last saved frame:
    - mean carrier-position mismatch about `7.12 Å`,
    - max carrier-position mismatch about `20.75 Å`,
    - mean peptide `C(i)-N(i+1)` bond in the reconstruction about `7.78 Å` versus about `1.37 Å` in the actual runtime carriers.
- Implemented fix:
  - `src/martini.cpp` now reads `hybrid_bb_map/reference_atom_coords` into runtime state and reconstructs the active `O` carrier from the current `N/CA/C` local frame before BB COM refresh.
  - `example/16.MARTINI/extract_martini_vtf.py` now exports actual runtime carrier coordinates when `hybrid_bb_map/atom_indices` resolve to runtime roles `N/CA/C/O`.
- Implemented follow-up fix for delayed instability:
  - `src/martini.cpp` now installs active hybrid `BB` proxy atoms into the dynamic fixed-atom mask whenever hybrid production is active, so those proxy beads no longer receive thermostat kicks or independent integrator updates.
  - `src/main.cpp` now applies `martini_fix_rigid::apply_fix_rigid_md(...)` immediately after startup momentum initialization/thermalization, so any active virtual `BB` proxies start from zero momentum before the first saved production frame.
- Reopened source-path finding after the user's rerun:
  - the active-`BB` momentum fix is real but not sufficient, because the user reran the actual workflow and reproduced the same production log through step `10000`.
  - `src/martini.cpp::align_active_protein_coordinates(...)` was still applying a rigid-body transform directly to the live integrated protein coordinates and momenta whenever `integration_rmsd_align_enable=1`.
  - project-local design notes in `example/16.MARTINI/task_plan.md` state that rigid-body alignment is intended for coupling coordinates only and that saved trajectories should preserve the raw integrated state.
  - this means the production integrator was being externally rewritten every alignment cycle, which is another non-Hamiltonian hybrid-energy source independent of the already-fixed `BB` virtual-site momentum bug.
- Implemented current source fix:
  - `src/martini.cpp::align_active_protein_coordinates(...)` now performs BB refresh and RMSD/reference bookkeeping only; it no longer rotates/translates the real integrated protein coordinates or momenta.
  - `src/deriv_engine.cpp` comments were updated to match the corrected semantics at the call sites.
- Focused verification:
  - rebuilt successfully with `cmake --build obj`;
  - strict-copy replay from the broken last saved frame confirmed runtime repair:
    - pre-run `/input/pos` `C-O` mean/max `17.90/30.30 Å`,
    - first saved `/output/pos` frame after replay `C-O` mean/max `1.32/1.57 Å`;
  - exporter verification on that replay file:
    - `use_runtime_carriers=True`,
    - exported AA backbone coordinates matched actual runtime carriers exactly on the checked frame (`0.0 Å` mean/max mismatch).
- Virtual-BB verification:
  - 5-step per-frame recorded-momentum replay from the exact workflow handoff (`6.6 -> 7.0`) after the new fix:
    - active `BB` proxy momentum remained exactly zero on every checked saved frame (`5/5` frames, max `0.0`);
    - active `BB` proxy positions still matched the carrier-weighted COM to numerical precision (max mismatch about `2.15e-6 Å`);
    - early MARTINI energy remained well behaved over the checked frames (`-23613.40 -> -23801.97` from step `0 -> 4`).
- Alignment-semantic verification:
  - rebuilt successfully with `cmake --build obj` after removing the live-state rigid-body rewrite.
  - exact stage-7 handoff A/B replay (`integration_rmsd_align_enable=1` vs `0`) over 5 saved frames now agrees to numerical noise only:
    - `/output/pos` max absolute difference about `3.36e-5 Å`,
    - `/output/potential` max absolute difference about `3.83e-5`,
    - protein-only coordinate max norm difference about `4.38e-5 Å`.
  - this confirms the stage-7 integration-alignment flag is now effectively coupling-side bookkeeping rather than a raw-trajectory mutation.
- Delayed-horizon verification status:
  - the exact fixed-code production replay from the real `6.6 -> 7.0` handoff is healthy through its `5000`-step checkpoint (`time 10.0`) and no longer matches the broken baseline values.
  - broken baseline at step `5000`: `potential 595.84`, `martini_potential -26362.14`, `total -25766.30`.
  - current replay at step `5000` / `time 10.0`: `potential 768.92`, `martini_potential -26414.40`, `total -25645.48`.
  - the replay was stopped after that verified midpoint because reaching the full `10000`-step horizon is hour-scale in this environment, but the original failure signature is already absent at the midpoint.
- Residual hybrid-virtual-site finding:
  - short momentum-enabled replays show that reconstructed active backbone `O` carriers still carry independent momentum even after the `BB`-proxy fix.
  - measured `O` momentum norm mean/max over the 5-step saved replay grew from about `4.78e-2 / 9.50e-2` to about `1.15e-1 / 1.94e-1`, while active `BB` proxy momentum stayed exactly `0.0`.
  - because `O` contributes about `0.2963` of each `BB` COM weight, that is a real correctness issue, but fixing it cleanly requires projecting `O` feedback through the `N/CA/C -> O` reconstruction map rather than merely freezing `O`.

## Lessons
- 2026-04-02: A subsystem replacement is not complete while prep/export helpers still emit the retired schema.
  - Working rule: after removing a runtime path, search the preparation and validation scripts for the same groups, attrs, and helper names. If the active workflow no longer consumes them, delete that export/validation code instead of leaving inert compatibility baggage behind.
- 2026-04-02: A direct-carrier rewrite is not complete until startup gating and legacy proxy bookkeeping are both ported or removed.
  - Working rule: when replacing a projected hybrid force path with a direct-Upside carrier path, audit more than the pair-force location. Explicitly check whether startup feedback ramps, hold windows, and proxy bonded/nonbonded terms still exist in parallel and will continue injecting energy or stale bookkeeping into production.
- 2026-04-02: Treat simulation-stability feedback from the user as the priority over a cleaner unit contract.
  - Working rule: if the user reports the box explodes after a unit change, restore the previously stable runtime contract first, verify the real workflow startup, and only then revisit the desired unit model with targeted diagnostics.
- 2026-04-02: Distinguish logging readability fixes from simulation-state fixes.
  - Working rule: if a log line shows absurd magnitudes, first verify whether the underlying value is real. Only then decide whether to patch formatting, physics, or both.
- 2026-04-02: When two changes were reverted together after a regression, do not keep attributing the failure to one change without re-isolating it.
  - Working rule: if the user later asks to restore one half of a reverted pair, re-test that change in isolation before repeating the earlier causal claim.
- 2026-04-01: Do not treat a unit or carrier-mapping cleanup as valid until the first real minimization stage reproduces the old energy scale.
  - Working rule: for hybrid workflow changes that touch masses or BB carriers, run the real `6.0` preparation/minimization path immediately and compare the initial MARTINI potential before making broader claims about correctness.
- 2026-04-01: When the user asks for native-unit preservation, inspect the staged arrays directly rather than trusting comments or historical unit conventions.
  - Working rule: for mass/unit corrections, verify the actual HDF5 payloads and separate physical dry-MARTINI particles from appended bookkeeping/reference rows before declaring the unit contract fixed.
- 2026-04-01: When the user says a proxy carrier mapping must collapse to one carrier, audit the generated HDF5 payloads as well as the force code.
  - Working rule: for hybrid mapping changes, verify the raw exported map, an early runtime stage file, and the final active stage file. It is not enough to change the runtime pair-force logic if the staged `hybrid_bb_map` still carries legacy active slots.
- 2026-04-01: In shell-driven workflows, a default assignment is not enough when downstream Python reads `os.environ`.
  - Working rule: if a child process depends on a variable, set workflow defaults with `export`, not a plain shell assignment, or the subprocess will still see the variable as missing.
- 2026-04-01: When making workflow unit defaults match the documented contract, trace every required conversion gate, not just the one the user mentioned first.
  - Working rule: if a runtime path requires both energy and length conversions, do not stop after defaulting one of them; check the startup guards and defaults together so the workflow is actually runnable with the documented unit contract.
- 2026-04-01: Do not call a direct-Upside replacement complete while legacy runtime stage artifacts still exist.
  - Working rule: after replacing a proxy-based hybrid force path, verify both the runtime math and the staged HDF5 payloads/logs. If stage files still carry old groups like `hybrid_sc_map` or runtime logs still report active SC proxy rows, the old path has not been fully removed from the actual workflow.
- 2026-04-01: When replacing a subsystem inside a stage-specific workflow, separate “sidechain path” from “everything the old helper happened to inject”.
  - Working rule: before deleting or bypassing a stage helper, enumerate which injected nodes belong to the user’s replacement target and which ones provide still-required physics; otherwise it is easy to delete the old sidechain path and accidentally delete the backbone force field with it.
- 2026-03-31: When the user corrects the target architecture, re-evaluate the plan against the corrected end state immediately.
  - Working rule: if the user says an existing subsystem is being removed or replaced, stop judging the proposal by compatibility with that subsystem and instead separate reusable infrastructure from obsolete design assumptions.
- 2026-03-31: When the user says a feature is being dropped, include the workflow scaffolding and benchmark path in that scope check.
  - Working rule: verify not just the C++ mechanism but also the driver scripts, injected HDF5 attrs/nodes, and the benchmark entrypoint, otherwise a removed subsystem may still be exercised by the nominal workflow.
- 2026-03-31: When the user tightens force semantics, verify both mechanics and units end to end.
  - Working rule: for hybrid interaction changes, check Newton-pair symmetry, the exact accumulation point where forces re-enter the protein coordinates, the exclusion rules for internal protein pairs, and the training-to-runtime unit contract in the same pass.
- 2026-03-31: When the user rejects a baked numeric contract, remove it from both metadata and runtime.
  - Working rule: if the user says a conversion must be parameterized, do not leave the old numbers lingering in manifests, docs, or runtime defaults and call the design fixed; trace the value through generation, stored attrs, and evaluation code.
- 2026-03-31: Parameterized does not mean "parameter with a baked fallback".
  - Working rule: if the user requires a value to be supplied as a parameter, do not hide the old constant behind a default; require the parameter explicitly or make the fallback behavior clearly legacy-only.
- 2026-03-31: A workflow is not portable just because it lives in one folder.
  - Working rule: when a user intends to upload a subtree independently, verify it from a standalone copied folder and trace defaults, bundled data, launcher scripts, and generated batch scripts for hidden parent-repo dependencies before calling it self-contained.
- 2026-03-31: Do not narrow a domain table default with an ad hoc subtype filter when the authoritative forcefield source is already available.
  - Working rule: if the workflow is supposed to build the complete forcefield table, default to the full bundled type list unless the user specifies an explicit subset, and treat any user-provided type list as the contract to validate against.
- 2026-03-31: A wrapper-level verification is not complete until the canonical entrypoint survives its own environment setup.
  - Working rule: when checking a staged workflow, run the top-level launcher itself at least once instead of only the underlying Python function, because shell-layer `set -u` / activation / environment sourcing bugs can invalidate an otherwise-correct execution model.
- 2026-03-31: Slurm spool copies make `$0` an unsafe source of workspace paths.
  - Working rule: for batch launchers, do not derive persistent workflow directories from `dirname "$0"` alone; prefer explicit submit-directory hints such as `SLURM_SUBMIT_DIR` and validate candidates by checking for expected workflow files.
- 2026-03-20: A production-stage fix is not validated by a short replay alone when the user reports a delayed failure at a later step count.
  - Working rule: if the user gives a later failure horizon (for example `10000` steps after a `5000`-step smoke passed), extend verification to that horizon or a targeted replay that reaches the same instability before calling the issue fixed.
- 2026-03-20: A mechanistic source-level fix is still not validated if the user reruns the actual workflow and reproduces the same long-horizon log.
  - Working rule: when the user reports “same output” after a rerun, immediately reopen the diagnosis, discard the previous fix as non-dominant, and validate the next candidate on the exact `stage handoff + seed + horizon` combination rather than on shorter proxy checks.
- 2026-03-20: A hybrid "alignment" feature can silently violate design intent if it mutates live integrator state instead of coupling/reference coordinates.
  - Working rule: when a control path sounds like diagnostic or coupling-side alignment, verify in code that it does not rewrite saved raw coordinates or momenta unless the design explicitly requires that behavior.
- 2026-04-13: When the user uses a bilayer-only calibration as a proxy for a hybrid control, separate the calibration Hamiltonian from the production Hamiltonian explicitly.
  - Working rule: if a user says a factor should be measured in a bilayer-only run but applied only at the Upside/dry-MARTINI interface, do not reuse the calibration scope as the production scope; implement those two surfaces independently and verify both.
- 2026-04-13: Before planning against an older compatibility node, re-check which runtime node is actually live in the current checkout.
  - Working rule: when multiple stage-7 paths coexist, inspect the active workflow artifact and codepath first; do not assume the older compatibility node is still the production path.
- 2026-04-13: When asked to build a new sweep workflow for an existing scientific run, orchestrate the canonical entrypoint instead of cloning its stage logic.
  - Working rule: if the repo already has one shell workflow that owns preparation, stage handoff, and activation logic, make the new sweep workflow drive that entrypoint with task-local env overrides rather than reimplementing the stage ladder in Python.
- 2026-04-13: When adding analysis to an existing scientific sweep, inspect the real output schema first and derive selectors from the generated files rather than from a nearby workflow.
  - Working rule: before implementing the analyzer, inspect a real stage file for atom selectors, box storage, and protein-carrier metadata; do not assume that `particle_class`, box datasets, or atom-role conventions match a different workflow.
- 2026-04-13: When the user explicitly asks for Slurm-array parity with a reference workflow, verify the live submission path against that reference and state the array granularity clearly in code/docs.
- 2026-04-13: When the user corrects the scientific scope of a workflow, do not preserve the old workflow surface if it still answers the wrong Hamiltonian.
  - Working rule: if a workflow folder is supposed to calibrate bilayer-only softening for later interface use, remove protein-relative controls and observables from that workflow rather than trying to reinterpret them in place.
- 2026-04-13: "Softening factor" and "softened potential" are not interchangeable controls.
  - Working rule: when the user asks for a simple percentage-like factor on LJ+Coulomb interaction strength, implement and analyze a scalar multiplier on the evaluated pair energy/force or its equivalent coefficient scaling; do not substitute soft-core LJ or Slater-softened Coulomb shape parameters.
  - Working rule: if a reference workflow is provided, inspect its actual `submit-*` wrapper and `workflow.py` array-submission code, then align empty-manifest guards, one-task-per-array-element behavior, and user-facing docs rather than only assuming the current implementation is “close enough”.
- 2026-04-13: Do not expand a workflow correction into the reference workflow when the user explicitly scoped the fix to the sweep/orchestration layer.
  - Working rule: if the user says the canonical example workflow is already correct and only the sweep wrapper is wrong, keep code changes inside the sweep workflow and treat the example path as read-only reference context.

## 2026-04-13 (Interface-Only Hybrid Scaling Calibration)
- Active stage-7 SC production path in this checkout:
  - `martini_sc_table_1body` is the live production path,
  - `martini_sc_table_potential` remains only as a backward-compatibility surface.
- The new shared interface factor can be applied in the direct pair runtime without regenerating MARTINI spline tables:
  - direct LJ contributions are linear in `epsilon`,
  - direct Coulomb contributions are linear in the charge product `q_i q_j`,
  - therefore a single multiplicative `interaction_scale` on the evaluated pair energy/force reproduces `epsilon *= scale` and `q -> sqrt(scale) * q`.
- Bilayer calibration implementation was verified directly on staged HDF5 artifacts:
  - for a `pair_scale = 0.85` smoke run,
  - `stage_6.0.prepared.up` had `epsilon` ratios of exactly `0.85`,
  - charge ratios were exactly `0.921954`, matching `sqrt(0.85)`.

## 2026-04-13 (Root-Level Hybrid Interface Sweep Workflow)
- The correct simulation surface for a production interface-scale sweep is `example/16.MARTINI/run_sim_1rkl.sh`:
  - it already owns hybrid packing,
  - stage-file generation,
  - stage `7.0` activation,
  - `PROTEIN_ENV_INTERFACE_SCALE`,
  - and the expected output layout.
- A thin orchestration workflow is enough for this sweep:
  - manifest,
  - one task-local `RUN_DIR` per scale/replicate,
  - Slurm array submission,
  - compact assembled status summary.
- The reduced smoke run verified that the new workflow can reach a real task-local `stage_7.0.up` checkpoint through the canonical hybrid entrypoint.
- 2026-04-29: User correction on C++ cleanup output policy.
  - MARTINI-hybrid progress output must permanently include protein structural/energy signals, not just total potential.
  - Required default line content for hybrid runs: `hbonds`, protein `Rg`, `prot_potential` (protein-side/non-MARTINI potential), and `total_potential`.
  - This must be auto-enabled for hybrid only, with no additional CLI switch, and must not alter non-hybrid simulation output behavior.
- 2026-04-29: User correction on `prot_potential` semantics.
  - `prot_potential` for MARTINI-hybrid progress lines must exclude dry-MARTINI bulk terms.
  - It must include hybrid interface terms:
    - SC-env table interface (`martini_sc_table_1body` / `martini_sc_table_potential`),
    - BB-env interface contribution inside `martini_potential`.
  - A simple `non-martini node sum` is insufficient because BB-env interface energy is embedded in `martini_potential`.

- 2026-04-29: User correction on cleanup scope: stage-7 continuation is required and must be retained.
  - Do not over-prune by deleting all continuation behavior when asked for baseline cleanup.
  - Keep explicit production continuation (`--continue-stage-70-from`) as part of active workflow requirements.
  - Prefer pruning auto-discovery/redundant continuation branches rather than removing continuation capability itself.
- 2026-04-29: User correction on continuation UX after cleanup.
  - Removing Python autodiscovery requires equivalent continuation detection in workflow bash scripts.
  - `run_sim_1rkl.sh` must auto-detect latest `stage_7.N.up` and pass explicit `--continue-stage-70-from`.
  - When passing explicit continuation source, bash should neutralize legacy auto flags to avoid parser/behavior conflicts.

## 2026-04-29 (User Correction: Hybrid Protein Runtime Representation)
- The required runtime representation for this refactor is strict:
  - remove all CG protein particles from stage files;
  - keep only protein backbone carriers (`N/CA/C/O`) as protein particles;
  - compute BB-env from per-residue COM of `N/CA/C/O` and distribute BB forces back to those carriers.
- Stage policy required by user:
  - hybrid interface active in stage 6.x and 7.x;
  - stage 6.x protein is rigid-body (free translation/rotation, no absolute pinning);
  - rigid-body constraint must be removed for stage 7.x.
- Mapping policy required by user:
  - use `pydssp` for secondary-structure classification and backbone-type mapping.
- Lesson:
  - when the user specifies a physics-level representation change, do not keep a “nearly equivalent” proxy-particle compromise; implement the exact particle set and force-routing semantics requested.
- 2026-06-20: Half-step IBI candidate status.
  - Added explicit sampled-bin IBI under-relaxation for CGL-CGL table training.
    This is recorded as `ibi_step_size` plus raw/applied correction metadata in
    the H5. It is a table-optimizer line-search parameter, not a runtime force
    cap, interaction scale, or orientation restraint.
  - Full-step IBI2 remains rejected: it improves apparent bilayer thickness
    but drives 1RKL hbond sum to `22.33` in the short smoke and shows upward
    potential drift.
  - 0.75-step IBI is also rejected: 1RKL `avg_kinetic_energy/1.5kT=1.682`
    and potential drift are unacceptable.
  - Half-step IBI is the best current CGL-CGL table candidate. With
    `mass_scale=0.012`, `rotational_tau=0.008`, two-mode GLE
    `memory_taus=0.2,2.0`, and `couplings=0.405,0.294`, the 20k CGL-only
    point gives `D40=0.254 um^2/s`, fit `R2=0.994`, stable geometry, and
    `min |n_z|=0.957`. The paired 1RKL 10k no-burnin smoke keeps final hbond
    sum `30.00`, exact temperature, and final `min |n_z|=0.958`.
  - This is not promotable yet. The short coupling line
    `0.38,0.405,0.42,0.46` gives `D40=0.346,0.254,0.213,0.205`, which is
    monotonic but has Stokes-Einstein line `R2=0.82`; longer replicated
    CGL-only checks are required before changing defaults or installed H5s.
  - The VTF-reported large CGL bilayer span is partly a visualization
    representation problem. The extractor emits only two endpoint atoms per
    CGL with fixed display span `25.698 A`, so endpoint atoms dominate visual
    percentiles compared with full DOPC's 14 beads per lipid. A physical
    display fix should reconstruct a bead cloud from stored DOPC reference
    geometry or table metadata; do not shorten display rods by an arbitrary
    factor.

- 2026-04-29: Local runtime verification for the AA-direct MARTINI workflow is currently blocked by an environment dependency gap.
  - `python3 py/martini_prepare_system.py run-hybrid-workflow ...` fails at import time with `ModuleNotFoundError: No module named 'tables'` from `py/martini_prepare_system_lib.py`.
  - Static verification (`py_compile`) and C++ build still pass, but end-to-end runtime checks require installing `tables` in the active `.venv`.
- 2026-04-29: The rigid-body fix-rigid extension in `src/martini.cpp` requires declaration order care.
  - `register_fix_rigid_for_engine` calls `rebuild_rigid_groups`; a forward declaration is required before that call site.
- 2026-04-30: User correction on rigid AA geometry integrity before stage 7.0.
  - A force-only rigid treatment is not sufficient in minimization if trial/final line-search coordinates are not projected to the rigid manifold.
  - Working rule: when protein is declared rigid in stage 6.x, project rigid group coordinates on every minimizer state update (initial/trial/accepted/final), not only via force filtering.
  - Working rule: stage handoff for AA-direct mode must be a strict coordinate copy path; do not apply per-residue carrier refresh/alignment transforms during 6.6 -> 7.0.
  - Working rule: hybrid interface setup must never mutate AA backbone carrier coordinates at runtime; interface math can remap interaction sites and gradients, but cannot overwrite state vectors.
- 2026-04-30: User correction on helix recognition at stage 7.0 start.
  - Working rule: when users report “not helix in VMD,” validate both HDF5 checkpoint geometry and exported VTF first-frame geometry with an independent DSSP check before attributing the issue to a specific stage.
  - Working rule: if minimization is not part of the last known-good workflow, keep stage-6 minimization disabled by default and require explicit opt-in.
- 2026-04-30: Stage-7 start structure comparisons must use backbone mapping order, not raw protein-membership order.
  - Working rule: compare stage `7.0` start against PDB using `hybrid_bb_map/atom_indices` (`N/CA/C/O` per residue) to avoid index-order ambiguity.
  - Working rule: when old output directories already contain drifted artifacts, validate on a fresh run directory before concluding the current code path still regresses.
- 2026-04-30: User correction on BB typing source and fallback behavior in AA runtime mode.
  - Working rule: when parity is requested against the last working martinize behavior, do not introduce DSSP-dependent BB typing in the active AA runtime path.
  - Working rule: for no-ITP AA runtime prep, use martinize fallback backbone typing (`ss="C"` table behavior) plus termini/chain-break `Qd/Qa` overrides.
- 2026-04-30: User correction on interaction parity scope.
  - Working rule: if the user asks for last-working BB-env/SC-env interaction logic, preserve AA virtual-BB COM coordinate/gradient routing but avoid adding new cross-interface role filters that change pair inclusion semantics.
- 2026-04-30: User correction on BB force routing scope in current AA-runtime framework.
  - Working rule: once BB is represented as an explicit virtual proxy, BB->carrier force distribution must remain enabled in production path as well; do not gate projection logic behind startup-only ramps.
  - Working rule: `ROLE_BB` must remain bound to explicit `BB` particles only (never reinterpret `CA` as BB in this framework), while `N/CA/C/O` serve only as carriers for COM/force projection.
- 2026-04-30: Stage-7-only instability root cause can be hidden in unconstrained backbone O carriers when BB COM includes O.
  - Working rule: if BB COM uses N/CA/C/O while stage-7 backbone constraints are N/CA/C-centric, restore the runtime O reconstruction path from reference N/CA/C frame before BB COM update.
- 2026-04-30: Multi-chain AA runtime prep must not key protein residues by shared segid.
  - Working rule: when generating `/input/sequence` and protein chain grouping from runtime PDB, prefer `chain_id` over `segid`; shared segid labels (e.g., `PROA`) can collapse distinct chains and break stage-7 sequence/mapping consistency.
- 2026-04-30: User correction on 1AFO topology semantics (two peptide segments, not one continuous chain).
  - Working rule: for multi-segment peptide systems, keep residue identity keyed by true chain IDs and ensure chain-break adjacency residues are excluded from inferred H-bond donor/acceptor construction (`i-1` and `i` at each chain break `i`).
  - Working rule: sequence length parity with `hybrid_bb_map` is necessary but not sufficient; stage-7 backbone/hbond node generation must also receive chain-break exclusions to avoid cross-segment backbone chemistry artifacts.
- 2026-04-30: Hybrid MARTINI nonbonded path should fail fast if table coverage is incomplete.
  - Working rule: when policy requires spline-only LJ/Coulomb (except Ewald), remove analytic fallback branches from runtime pair evaluation and throw explicit runtime errors if a nonzero interaction lacks a corresponding spline.

## 2026-05-06 (Dynamic Single-Particle Lipid Orientation)
- The prior single-particle CGL lipid orientation was static in practice:
  - `cg_lipid_pair` produced direction derivatives,
  - but `compose_vector6d` discarded them because lipid directions were copied from a fixed HDF5 `direction` dataset.
- A minimal rotatable-vector representation can be added without creating extra interacting lipid beads:
  - append one hidden `CGLD` site per CGL lipid,
  - constrain CGL-CGLD distance with a harmonic bond,
  - compute the direction from normalized CGL-to-CGLD displacement,
  - propagate the normalized-vector derivative to both sites.
- The dynamic-vector bilayer-only smoke confirms vectors now rotate:
  - direction norms remain exactly unit-length in diagnostics,
  - 5000-step max/mean angular changes reached `149.876/65.677 deg`.
- The bilayer still spreads in `z` over 5000 steps, though less than the static-vector run:
  - static-vector final `z_std=13.384 Å`,
  - dynamic-vector final `z_std=10.390 Å`.
- Lesson:
  - when a directional potential is intended to rotate coarse particles, verify the full derivative path from pair evaluator through the coordinate node to physical integrator coordinates; accumulating a derivative on a virtual direction component is not sufficient if the upstream node drops it.

## 2026-05-06 (Full Directional CG Lipid Bilayer Stabilization)
- The unstable bilayer was not caused only by frozen vectors:
  - the residual-only directional table left the large isotropic CGL-CGL LJ core active;
  - that core strongly repelled same-leaflet neighbors at the observed initial spacing;
  - the residual directional term did not provide enough cohesive tail-tail or same-leaflet packing attraction.
- The stable local bilayer-only configuration uses:
  - `cg_lipid_quadspline_v3` full multimode CGL-CGL table;
  - zero CGL-CGL `martini_potential` coefficients when `cg_lipid_pair` is active;
  - explicit tail-tail and side-by-side cohesion modes;
  - `UPSIDE_CG_LIPID_PAIR_SCALE=0.02` default, applied only to radial energy channels;
  - initial same-leaflet CGL spacing conditioning to min/p05 `7.000/7.000 Å`.
- Verification result:
  - 72-DOPC bilayer-only `50000`-step NVT run stayed finite and slab-like;
  - `z_std` went from `6.512 Å` initially to `5.968 Å` at time `100`.
- Lesson:
  - when scaling spline tables with dimensionless angular factors, scale only energy-valued radial channels; scaling angular channels as well changes the effective energy by extra powers of the scale factor.
- 2026-05-12: Historical commit `3a98be1a4b4ebbbdf645fd1db1dcb84efa86af1e` kept the hybrid stage-6 branch available, with `6.0` prepared under NPT-enabled preproduction settings and later stage-6 MD/barostat relaxations on the extended bonded-environment route.
  - Working rule: when restoring the current CG-lipid path, preserve the historical stage-6 barostat contract and recover an explicit active `6.0` MD relaxation without deleting the bonded-environment branch.
- 2026-05-12: The fresh bilayer-only probe using the recovered `6.0` zero-target semi-isotropic NPT contract did not show meaningful short-horizon shrink: `50.091999 -> 50.100193 Å` in `x/y` over `500` steps, with `z=85 Å` unchanged.
  - Working rule: describe this verification as a contract and trend check, not as proof of a historical shrink magnitude unless an archived comparison artifact is available.
- 2026-05-12: The user supplied `/Users/yinhan/Documents/upside2-md-martini` as the bead-resolved dry-MARTINI reference repo for the required bilayer-only NPT comparison.
  - Working rule: for this box-relaxation question, compare the current single-particle model against a fresh dry-MARTINI bilayer-only run from that repo rather than treating historical workflow shape as a sufficient proxy.
- 2026-05-12: Direct matched bilayer-only NPT comparison shows opposite XY box trends over `500` steps:
  - bead-resolved dry-MARTINI: `50.091999 -> 50.071121 Å`, area `-0.083341%`;
  - current single-particle model: `50.091999 -> 50.100193 Å`, area `+0.032718%`.
  - Working rule: do not describe the current single-particle NPT box response as dry-MARTINI-consistent until its XY relaxation sign and scale agree in this direct bilayer-only comparison.
- 2026-05-28: B-spline angular underdetermination causes spurious CG-CG energy features.
  - Rule: the number of angular B-spline control points must not greatly exceed the number of angular data samples. With 15×15=225 angular controls fitting 7×7=49 data points per radial distance (4.6× underdetermination), Tikhonov regularization with λ=0.01 is too weak to prevent oscillations between sample points. These spurious features create unphysical orientational preferences in simulation.
  - Fix: reduce n_knot_angular from 15 to min(cos_theta_count + 2, 15) = 9, increase smooth from 0.01 to 0.1. This reduces angular std by 62-70% and eliminates most out-of-bounds interpolated values.
  - The same pattern exists in the SC-CGL table (9²=81 data → 15²=225 controls, 2.8× ratio) but is less severe.
  - When fitting tensor-product B-splines to underdetermined data, prefer matching knot resolution to data resolution over weak regularization. Weak regularization with excess DOFs creates interpolation artifacts that look like physical structure but are not.
