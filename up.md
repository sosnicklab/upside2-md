# Upside File Format Reference

All Upside data files are HDF5 archives. The `.up` extension is used for both simulation
input/state files and dynamic output files. Parameter force-field files use `.h5`. A few
scalar parameter files use plain text (`bb_env.dat`, `sheet`).

Units throughout are Upside-native unless marked otherwise: energy in E_up
(1 E_up = 2.914952774 kJ/mol), length in Angstrom, mass in m_up (12 g/mol),
temperature in T_up (1.0 T_up = 350.59 K). The conversion from dry-MARTINI
native units (nm, kJ/mol, e) to Upside units happens once at h5-build time in
Python; the C++ engine reads pre-converted tables and performs no further conversion.

---

## 1. Simulation Input/State File (`.up`)

Created by `py/upside_config.py` for pure-protein runs and by
`py/martini_prepare_system.py` for MARTINI and hybrid runs.
The same file accumulates output trajectories during simulation.

### Top-level groups

```
/input/      system topology, force-field nodes, initial state, dynamics params
/output/     trajectory (appended by the engine; absent before first run)
```

### 1.1 `/input` datasets (all system types)

| Dataset | Shape | dtype | Notes |
|---|---|---|---|
| `pos` | `(N, 3, 1)` | float32 | Initial coordinates in Angstrom; last dim is frame placeholder |
| `mom` | `(N, 3, 1)` | float32 | Momenta to restart from, read only with `--restart-using-momentum`; `@round` is the round they were recorded at, from which the engine resumes the time and every random counter; `@replica_index` is the replica the system held (replica exchange), read with or without the momentum (a position-only continuation keeps its walkers) and checked to be a permutation over the systems. Copied from `/output/mom`, `/output/round` and `/output/replica_index` of the frame restarted from (`generate_restart_config.py`, the examples' `continue_sim`) |
| `mass` | `(N,)` | float32 | Per-atom mass in m_up |
| `atom_names` | `(N,)` | `S4` | Atom names matching PDB convention (N, CA, C, PO4, NC3, ...) |
| `atom_roles` | `(N,)` | `S4` | The atom names again (N, CA, C, O, BB, PO4, ...); the engine reads the carrier and phosphate roles from it |
| `particle_class` | `(N,)` | `S10` | One of: PROTEIN, LIPID, ION (hybrid) |
| `charges` | `(N,)` | float64 | Partial charge in e |
| `molecule_ids` | `(N,)` | int32 | Per-atom molecule index (0-based, contiguous) |

Additional datasets present in protein or hybrid runs:

| Dataset | Shape | dtype | Notes |
|---|---|---|---|
| `sequence` | `(n_res,)` | `S3` | Three-letter AA codes (ASP, GLU, ...) |
| `type` | `(N,)` | `S8` | MARTINI bead type (Qd, P5, Na, ...) or atom element |
| `thermostat_timescale` | `(N,)` | float32 | Per-atom Langevin friction timescale; protein ~5, frozen beads ~1e8 |
| `residue_ids` | `(N,)` | int32 | 1-based residue index per atom |
| `residue_names` | `(N,)` | `S4` | PDB residue name per atom (hybrid); names each lipid bead's molecule for the per-lipid friction |

### 1.2 `/input/args` group (protein-only runs)

Stores the full command-line argument set as HDF5 attributes. Key attributes:

| Attribute | Type | Notes |
|---|---|---|
| `fasta` | bytes | Path to FASTA file |
| `initial_structure` | bytes | Path to `.npy` initial coordinate file |
| `invocation` | bytes | Full command line used to create this file |
| `bond_stiffness` | float64 | k_bond in E_up/Ang^2; default 48 |
| `angle_stiffness` | float64 | k_angle in E_up; default 175 |
| `omega_stiffness` | float64 | k_omega in E_up; default 30 |
| `rama_library` | bytes | Path to `rama.dat` (HDF5) |
| `rama_sheet_mixing_energy` | bytes | Path to `sheet` (plain text) |
| `reference_state_rama` | bytes | Path to `rama_reference.pkl` |
| `hbond_energy` | bytes | Path to `hbond.h5` |
| `rotamer_placement` | bytes | Path to `sidechain.h5` |
| `rotamer_interaction` | bytes | Path to `sidechain.h5` |
| `environment_potential` | bytes | Path to `environment.h5` |
| `bb_environment_potential` | bytes | Path to `bb_env.dat` |
| `membrane_potential` | bytes | Path to `membrane.h5` (optional) |
| `dynamic_rotamer_1body` | uint8 | Enable temperature-dependent 1-body rotamer |
| `hb_scale` | float64 | Global scale factor for H-bond potential |
| `env_scale` | float64 | Global scale factor for environment potential |
| `rot_scale` | float64 | Global scale factor for rotamer pair potential |
| `memb_scale` | float64 | Global scale factor for membrane potential |

### 1.3 `/input/potential` group

Contains one sub-group per active force-field node. Each sub-group has at minimum
an `@arguments` attribute listing which upstream node outputs it consumes. A potential may carry
`@integrator_level`, read only under `--integrator mv`: 1 puts it on the slow level (one kick of
inner_step*dt per round), any other value (default 2) on the fast level. Coordinate nodes take their
levels from the potentials that consume them.

#### Node dependency chain (protein run)

```
pos
 |- Distance3D          -> Spring_bond
 |- Angle               -> Spring_angle
 |- Dihedral_phi/psi    -> rama_coord
 |- Dihedral_omega      -> Spring_omega
 |- affine_alignment    -> backbone_pairs
                         -> placement_fixed_point_vector_only (CB)
                         -> placement_fixed_point_vector_scalar (hyd)
                         -> placement_scalar
 |- rama_coord          -> rama_map_pot / rama_map_pot_ref
 |- infer_H_O           -> protein_hbond -> hbond_energy
 |- placement_*         -> hbond_coverage
                         -> hbond_coverage_hydrophobe
                         -> environment_coverage_sc / _hb -> sigmoid/bb coupling
 |- rotamer             (self-consistent; uses placement + pair_interaction table)
```

#### Backbone geometry nodes

| Node | Key datasets | Notes |
|---|---|---|
| `Distance3D` | `id (M,2)` | Atom index pairs |
| `Angle` | `id (M,3)` | Outputs the cosine at the third atom: (a0 - a2).(a1 - a2) normalised |
| `Dihedral_phi` | `id (n,4)` int64; -1 in slot 0 at every chain start | Phi dihedral atom quads |
| `Dihedral_psi` | `id (n,4)` int64; -1 in slot 3 at every chain end | Psi dihedral atom quads |
| `Dihedral_omega` | `id (n,4)` int64 | Omega dihedral atom quads |
| `Spring_bond` | `id (M,)`, `equil_dist (M,)`, `spring_const (M,)` | `@pbc=0`, `@dim1=0` |
| `Spring_angle` | `id (M,)`, `equil_dist (M,)` (cosine), `spring_const (M,)` | Angles stored as cosines |
| `Spring_omega` | `id (M,)`, `equil_dist (M,)` (pi, or 0 before a cis-proline), `spring_const (M,)` | `@pbc=1`, `@box_len=2pi` |
| `Const3D[_name]` | `id (M,)`, `value (M,3)` | Fixed points (walls, springs, nail and group-nail anchors, the curvature centre), written from the config; `id` sizes the node |
| `cb_membrane_potential` | `cb_index (n,)`, `res_type (n,)`, `weights (20,)`, `bl_nodes (n_bl-1,)`, `coeff (20,n_bl,n_node)`; `@z_start`, `@z_scale`, `@bl_range_sharpness`, `@left_right_node`, `@left_right_sharpness`, `@use_curvature`, `@curvature_sign` | Arguments `placement_fixed_point_only_CB`, `environment_coverage_sc`, `Const3D_curvature_center`. One row per residue not in `--membrane-exclude-residues`; the energy of its CB at its depth, blended between burial bins by the environment coverage and between the two leaflets at the bilayer centre. With `@use_curvature` the depth is measured from the sphere around the curvature centre, of radius its \|z\| |
| `hb_membrane_potential` | `hb_index (n,)`, `hb_type (n,)` (0 donor, 1 acceptor), `coeff (2,2,n_node)`; attributes as `cb_membrane_potential` without the burial ones, and `@half_thickness` | Arguments `protein_hbond`, `Const3D_curvature_center`. `@half_thickness` bounds the radius `--curvature-changer-interval` samples (the changer reads it from this node). One row per backbone H and O of `infer_H_O` outside the excluded residues, its energy blended between unpaired and paired by its H-bond |
| `cb_surf_membrane_potential`, `hb_surf_membrane_potential` | as above without `@use_curvature` and `@curvature_sign`, plus `coeff_inner` (`icb_energy`, `ihb_energy`), and for H/O `hb_res (n,)`, `weights (20,)` | The channel potential (`--channel-membrane-potential`), flat membranes only: arguments as above with `surface` in place of `Const3D_curvature_center` (the H/O node also takes `environment_coverage_sc`). Inside the slab \|z\| < half_thickness a residue's energy mixes `coeff` and `coeff_inner` by its surface score; outside it is `coeff` (a hard switch, pass 2 membrane M2) |
| `cat_pos_bb_cb` | `index_pos1 (3 n_res,)`, `index_pos2 (n_res,)`; `@n_dim=3` | The backbone N, CA, C of `pos` followed by the CBs of `placement_fixed_point_only_CB`, for `surface` |
| `surface` | `index_cb (n_res,)` (row of each CB in `cat_pos_bb_cb`), `auxiliary_point (n_aux,2)`, optional `included_list`; `@method_type` (0 or 1), `@exclusion` (0, or 1 with `included_list`), `@n_rotation=20`, `@half_thickness`, `@zbuffer`, `@len_gridx`, `@len_gridy`, `@n_gridz`, `@n_smallgridx/y/z`, `@atom_radius=1.6` | One output per residue: the fraction of `@n_rotation` orientations about the membrane normal in which the residue's CB lies on the protein's lipid-facing surface, found on voxels of the residues in the flat slab \|z\| < half_thickness + zbuffer. Refused with a curved membrane |
| `MovingConst3D` | `id (M,)`, `velocities (M,3)`, `start_pos (M,3)`, `@time_initial`, `@time_step`, `@step_initial` | AFM tips at start_pos + velocities*t, t = time_initial + time_step*steps, the step count starting at step_initial (0 in a fresh config) and advancing one per force evaluation (Verlet, or mv's fast level only); a table without tip positions starts each tip on its CA at time 0; a continuation (`generate_restart_config.py`, example 06) sets step_initial from `/output/tip_step` |
| `affine_alignment` | `atoms (n_res,3)`, `ref_geom (n_res,3,3)` | N/CA/C local frame per residue, its origin the atoms' centroid. Each residue's `ref_geom` must be centred (the engine refuses a centroid beyond 1e-4 A), and every placement point is relative to that centroid |
| `backbone_pairs` | `id (n_res,)` (the `affine_alignment` row), `sequence_id (n_res,)`, `n_atom (n_res,)`, `ref_pos (n_res,4,3)` | Backbone steric wall between the N, CA, C and CB of residues at least 2 apart in `sequence_id` (6, pair exclusion), so the residues across a chain break interact |

#### Rama / secondary structure nodes

| Node | Key datasets/attrs | Notes |
|---|---|---|
| `rama_coord` | `index_pos1 (n_res,)`, `index_pos2 (n_res,)`; `@arguments=[Dihedral_phi, Dihedral_psi]` | Pairs row `index_pos1` of phi with row `index_pos2` of psi into (phi, psi) per residue. `id (n_res,5)`, the atoms of each residue's phi and psi, is read by the `pivot_moves` writer, not by the engine |
| `rama_map_pot` | `residue_id`, `rama_map_id`, `rama_pot (n_res,72,72)`, `rama_map_id_all (n_res,)` | 5-deg grid, one map PER RESIDUE, **not** a library map; `rama_map_id_all` is every residue's layer, read by the `pivot_moves` writer as `pivot_restype`; `--rama-param-deriv` adds `@restype` (the sheet library's residue order), `@sheet_eps` and the `more_sheet_rama_pot_*`, `less_sheet_rama_pot_*` maps for ConDiv. The maps are the node's parameters (`get_param`, `set_param`); their derivative is not implemented. See 2.8a |
| `rama_map_pot_ref` | `residue_id`, `rama_map_id (n_res,)` all 0, `rama_pot (1,72,72)` | Reference-state correction: one map shared by every residue; `@log_pot` |

#### Sidechain rotamer nodes

| Node | Key datasets | Notes |
|---|---|---|
| `placement_fixed_point_vector_only` | `affine_residue (n_sc,)`, `layer_index (n_sc,)`, `placement_data (n_layer,6)` | SC1 bead position (vec3 + unit vec3 direction); `layer_index` selects the library row (86 in `rotamer_center_fixed`) |
| `placement_fixed_point_vector_only_CB` | same | CB position |
| `placement_fixed_point_vector_scalar` | `placement_data (n_layer,7)` | Hydrophobe position (vec3 + dir + scalar), 3 library rows (`hydrophobe_placement`) |
| `placement_scalar` | `affine_residue`, `layer_index`, `rama_residue`, `placement_data (n_layer,36,36,1)` | Rotamer 1-body energy, -log P(rotamer given phi, psi) per library row at the residue's (phi, psi); `placement_fixed_scalar` holds one value per row |
| `rotamer` | `@arguments`, `@damping` (in [0,1); the writer's default 0.4), `@max_iter=1000`, `@tol=1e-3`, `@iteration_chunk_size=2`, all required; sub-group `pair_interaction` | Belief propagation over rotamer states; it stops when the undamped residual max|B(b) - b| is within `tol` |
| `rotamer/pair_interaction` | `index`, `type` (`bead_order` row), `id` (rotamer, n_rot<<4, node<<8: the node holder slot), `sequence_id` (the bead's residue), `interaction_param` (loaded from `sidechain.h5`) | SC-SC interaction graph; the beads of one residue, equal `sequence_id`, do not pair (6, pair exclusion) |

#### H-bond and environment nodes

| Node | Key datasets | Notes |
|---|---|---|
| `infer_H_O` | sub-groups `donors` and `acceptors` with `residue`, `bond_length`, `id` | Backbone H and O position inference |
| `protein_hbond` | `index1`, `id1`, `sequence_id1`, `type1` (donors: `infer_H_O` rows, residues, their sequence ids, parameter rows); `index2`, `id2`, `sequence_id2`, `type2` (acceptors) | Interaction graph; donor-acceptor pairs are found at run time, excluding those within 1 in `sequence_id` (6, pair exclusion) |
| `hbond_energy` | `parameters`, `donor_resid`, `acceptor_resid`, `rama_resid`, optional `residue_class (n_res,)` int32, `@integrator_level=1` | From `hbond.h5` (2.2): 12 values, or 12 + 3 per residue class with `residue_class` |
| `hbond_coverage` | interaction graph: `index1/2`, `id1/2` (residues), `sequence_id1/2`, `type1/2`, `interaction_param` | Side-chain coverage of each H-bond site (`protein_hbond` row); `type1` = 2 x backbone class + (donor 0, acceptor 1). This node and the coverage nodes below exclude the pairs within 2 in `sequence_id` (6, pair exclusion) |
| `hbond_coverage_hydrophobe` | interaction graph: `index1/2`, `id1/2` (residues), `sequence_id1/2`, `type1/2`, `interaction_param` | Side-chain coverage of each residue's three backbone hydrophobe points; `type1` = 3 x backbone class + point |
| `environment_coverage_sc` | interaction graph as above; `@num_aa_types=20`, `aa_types (n_res,)` | Side-chain coverage of each CB, split by the covering residue's type `aa_types[id2]` |
| `environment_coverage_hb` | same | The same for each backbone H and O. Files written before 2026-10-19 hold one type per H/O site, read as master reads them (findings 3.13) |
| `hb_environment_coverage_hn`, `hb_environment_coverage_oc` | interaction graph as above, both groups `infer_H_O` rows (hn: H then O; oc: O then H) | Coverage of each H by the O, and of each O by the H, for `cat_pos_bb_coverage` |
| `hbbb_coverage` | interaction graph as above, group 1 the donor rows of `infer_H_O`, group 2 the N, CA, C of `pos` | Backbone coverage of each H (`--use-heavy-atom-coverage`), read by `get_protection_state.py` only |
| `sigmoid_coupling_environment` | `@integrator_level=1`, `@number_independent_weights` (20: one weight per covering residue type, tied across coupling types; 400: every pair) | SC-environment coupling (ff_2.1) |
| `bb_sigmoid_coupling_environment` | `weights (20,)`, `@num_aa_types=20`, `@scale`, `@center`, `@sharpness`, `@hbond_weight` | Backbone burial (2.4) |

#### MARTINI-specific nodes (pure lipid or hybrid)

| Node | Key datasets/attrs | Notes |
|---|---|---|
| `martini_potential` | `@cutoff`, `@n_points=1000`, `@optimized_format`, `@soft_core_lambda`, `@x_len/y_len/z_len`, `combined_energy_grids (n,1000)`, `combined_energy_end_slopes (n,2)` | Precomputed LJ+Coulomb spline table; per-bead-type-pair energy grid, fitted as a cubic spline clamped to the analytic dV/dr (E_up/A) at `@r_min_ang` and, from inside, at `@r_max_ang`, where the potential-shifted LJ keeps a force. `@soft_core_lambda` 0 is the published potential; 0.01 (stage 6.0 only) is its GROMACS soft core, (1-lambda) V(r_sc) with r_sc = (4 sigma^6 lambda^2 + r^6)^(1/6), as CHARMM-GUI's step 6.0 minimisation uses; the preparation's `inject_particles_table` builds the tables from it, and the engine reads only the tables |
| `martini_sc_table_1body` | `env_atom_index`, `env_target_index`, `row_*`, `grid_ang`, `theta_grid`, `rotamer_count`, `rotamer_full_energy_eup`, `@x_len/y_len/z_len`; optional `@cache_buffer` (default 1 A, not written by the preparation) | Side-chain/environment energy per rotamer row from the `/sc_table` of `martini.h5` (2.6): a clamped cubic spline in r and theta, cut off at the last radial node; each environment atom maps to the target of its bead type and charge |
| `dist_spring` | `id (N,2)`, `equil_dist (N,)`, `spring_const (N,)`, `@x_len/y_len/z_len` | Intra-lipid bond springs, V = k/2 (r - r0)^2 from a canonical spline that spans every minimum-image length |
| `angle_spring` | `id (N,3)` (vertex in the middle), `equil_angle_deg (N,)`, `spring_const (N,)`, `@x_len/y_len/z_len` | Intra-lipid angles, V = k/2 (cos(theta) - cos(theta0))^2 (GROMACS type 2) from a canonical spline in cos(theta) - cos(theta0) over [-2, 2] |

#### Hybrid protein+MARTINI extra nodes

| Node | Notes |
|---|---|
| `affine_alignment` | Protein-only backbone frame (same structure as pure protein) |
| `backbone_pairs` | Protein backbone steric wall |
| `hbond_energy` | Protein H-bond energy from `hbond.h5` |
| `infer_H_O` | Protein H/O inference |
| `protein_hbond` | Protein donor-acceptor pairs |
| `rama_map_pot / _ref` | Protein Ramachandran |
| `rotamer` | Protein rotamer, one row per side-chain bead of every residue exactly as `upside_config` writes it; arguments `placement_fixed_point_vector_only`, `placement_scalar` (Upside's rama-dependent 1-body on `rama_coord`, as `--dynamic-rotamer-1body`), `martini_sc_table_1body` |
| `martini_potential` | MARTINI non-bonded (all bead types) |
| `martini_sc_table_1body` | SC-bilayer 1-body offset, one value per rotamer row; `row_residue_table_index` -1 marks ALA and GLY, which have no MARTINI side chain and no environment term |
| `restraint_position` | Stages 6.2-6.6 only: harmonic PO4 height restraint of CHARMM-GUI's ramp (`BILAYER_LIPIDHEAD_FC` 200, 100, 50, 20, 10 kJ/mol/nm^2, z only), reference the assembled membrane in the stage file's frame |

### 1.4 `/input` groups (MARTINI and hybrid only)

| Group | Key attrs/datasets | Notes |
|---|---|---|
| `barostat` | `@target_p_xy=0`, `@target_p_z=0`, `@interval=10`, `@mc_dmax_xy=0.01`, `@mc_dmax_z=0` (frozen z) | Monte-Carlo semi-isotropic barostat; XY only for membrane. Each trial scales molecule centres (`/input/molecule_ids`, taken whole across the boundary) and accepts on exp(-beta[dU + P dV]) (V'/V)^(N+1). A molecule with a fixed atom is not scaled (the derived BB/O slots do not count); a z-fixed atom blocks only the z move |
| `brownian` | `atom_index (N_dyn,)`, `friction (N_dyn,)`, `protein_contact_atom_index`, `protein_lipid_contact_count`; `@numerical_time_step=0.009`, `@reference_temperature_up=0.8647`, `@lipid_names`, `@lipid_target_diffusion_um2_s`, `@ion_target_diffusion_um2_s`, `@protein_time_ps_per_step=40`, `@martini_time_factor=4`, `@dynamics_phase` (`equilibration` or `production`), `@inner_steps=4`; production adds `@lipid_bead_diffusion_up`, `@lipid_bead_friction_up`, `@ion_bead_friction_up`, `@interface_friction_cutoff_A=12` | g-JF Langevin atom list. Lists **every dynamic particle**. In production each lipid bead has its own molecule's bare-particle friction `kT_ref / D_bead,up`, where `D_bead,up` maps the lipid's measured lateral diffusion (`LIPID_LATERAL_DIFFUSION_UM2_S` in `py/martini_prepare_system_lib.py`, one cited entry per lipid; an untabulated lipid is refused) onto the factor-four, 10 ps-per-step clock; ions take the DOPC entry's clock as an interim rule. A protein N/CA/C carrier has the sum of the frictions of the lipid beads within 12 A (so a carrier with no lipid within 12 A has friction 0 and reduces to velocity Verlet, thermostatted by the OU bath instead), recounted at each stage handoff from `/input/residue_names`. O and sidechain slots are absent because the `placement_*` nodes derive them. The engine scales each friction by `T/T_ref` at runtime so the target diffusion holds at the temperature actually run. In the preparation stages (`equilibration`) every bead has the dry-MARTINI Langevin friction `mass / 4 ps` and the rigid protein's carriers 0. `@inner_steps` splits each outer step into N substeps of `dt/N`; the **outer** dt is the one the friction and the ps-per-step clock are calibrated against and must stay 0.009, which the engine enforces by throwing on mismatch. |
| `hybrid_bb_map` | `atom_indices (n_res,4)` (N, CA, C, O), `atom_mask (n_res,4)`, `weights (n_res,4)`, `bb_atom_index (n_res,)`; metadata the engine does not read: `bb_residue_index`, `bb_type`, `bb_chain_id`, `bb_resseq`, `bb_secondary_structure`, `reference_atom_names` | One row per protein residue. `martini_hybrid_position` places the BB bead at `bb_atom_index` as the weighted centre of N/CA/C and the `infer_H_O` carbonyl O (weights 14/12/12/16 over 54); a residue with no acceptor (a chain's last) uses N/CA/C with the weights renormalised |
| `hybrid_control` | `@activation_stage` (`minimization` or `production`); in the `production` stage also `@po4_z_hold_end_round` | Hybrid stages: see the stage table below |
| `hybrid_env_topology` | `env_atom_indices (N_env,)`, `protein_membership (N,)` | `protein_membership` holds each atom's protein id, or -1 for an environment particle |
| `stage_parameters` | `@current_stage` (`minimization`, `production_handoff` or `production`) | Hybrid stage label |

A hybrid configuration is built for preparation (`activation_stage=minimization`) or for production
(`activation_stage=production`) and runs only in the stages built for it. The engine refuses any other pairing,
including an unknown or misspelled label. The SC-env and BB-env interface is on in every stage.

| `current_stage` | `activation_stage` | Workflow stage | Protein | Lipid PO4 |
|---|---|---|---|---|
| `minimization` | `minimization` | 6.0 soft-core minimisation, 6.1 minimisation, 6.2-6.6 NPT equilibration | N/CA/C one rigid body | free; 6.2-6.6 carry the harmonic `restraint_position` height restraint |
| `production_handoff` | `production` | 7.0 minimization and burn-in | N/CA/C one rigid body | free |
| `production` | `production` | 7.0 production and continuations | flexible | z held in the rounds before `po4_z_hold_end_round`, 150 rounds after the burn-in, by a run that continues the trajectory from its momentum (`--restart-using-momentum`); a run started from positions alone holds nothing |

The workflow follows CHARMM-GUI's Martini protocol: 6.0 minimises on the soft-core tables, 6.1 on the published
ones, and every MD stage runs the published tables. Before each MD stage the workflow refuses a configuration
holding a MARTINI pair whose energy exceeds kT ln(1e9 N_pairs), the energy an equilibrium configuration shows
less than once in 1e9.
`--minimize` minimises under the configuration's own stage. Each stage's seed is a splitmix64 hash of the base
seed and the stage name; production and its continuations share one, and the momentum and its round are carried
through every handoff, so a continued production is the continuous run.

### 1.5 `/input` restart state (hybrid)

A hybrid run restarts from `/input/mom` and its `@round` (1.1), as a protein-only run does, and also from the box
of the frame (in every node that carries one). The PO4 hold ends at a round of that clock, so it is no restart state
of its own. `set_initial_position` in `martini_prepare_system_lib.py` and `generate_restart_config.py` both carry the
momentum, its round and the box.

### 1.6 `/input` Monte Carlo samplers (protein-only)

| Group | Key datasets | Notes |
|---|---|---|
| `pivot_moves` | `pivot_atom (n,5)`, `pivot_range (n,2)`, `pivot_restype (n,)`, `proposal_pot (n_res,72,72)` | MC pivot move table; `proposal_pot` is each residue's own `rama_pot` layer and `pivot_restype` indexes it |
| `jump_moves` | `atom_range (n_chain,2)`, `sigma_trans (n_chain,)`, `sigma_rot (n_chain,)` (radians) | Rigid chain translations and rotations; written instead of `pivot_moves` for multi-chain configs |

### 1.7 `/output` datasets

The engine appends one row per saved frame to each extendable dataset.

| Dataset | Shape | dtype | Notes |
|---|---|---|---|
| `pos` | `(n_frame, 1, N, 3)` | float32 | Coordinates per frame (each system writes its own file): the state at the start of the frame's round, before its Monte Carlo move and re-centering, so a continuation from it applies each once. A hybrid's BB and O slots are written as `martini_hybrid_position` derives them |
| `mom` | `(n_frame, 1, N, 3)` | float32 | Momenta per frame (`--record-momentum`) |
| `round` | `(n_frame,)` | int32 | Integration round of each frame (`--record-momentum`) |
| `replica_index` | `(n_frame, 1)` | int32 | Replica the system holds at the frame (replica exchange) |
| `potential` | `(n_frame, 1)` | float64 | Total potential energy in E_up |
| `kinetic` | `(n_frame, 1)` | float64 | Kinetic energy per three degrees of freedom in E_up (1.5 kT in equilibrium): fixed atoms, the z of z-held atoms and the internal motion of a rigid group are not counted. The final `avg_kinetic_energy/1.5kT` line averages it over the second half of the frames |
| `protein_kinetic`, `lipid_kinetic` | `(n_frame, 1)` | float64 | The same over the PROTEIN and the LIPID particles (hybrid) |
| `temperature` | `(n_frame, 1)` | float64 | Thermostat temperature of the frame in T_up (the annealing schedule's under `--anneal-factor`), not a measured temperature |
| `time` | `(n_frame,)` | float64 | Simulation time |
| `hbond` | `(n_frame, n_donor + n_acceptor)` | float32 | H-bond count of each donor, then of each acceptor (`protein_hbond` column 6), so each bond appears twice; count it from the donors |
| `rama_map_potential` | `(n_frame, n_row)` | float32 | Ramachandran energy per `residue_id` row of `rama_map_pot`: every residue, except that under `--trans-cis` the switching prolines are in `RamaMap2` and absent here |
| `rotamer_potential` | `(n_frame, n_res)` | float32 | Per-residue rotamer pair energy |
| `rotamer_free_energy` | `(n_frame, n_res)` | float32 | Per-residue rotamer free energy |
| `rotamer_entropy` | `(n_frame, n_res)` | float32 | Per-residue rotamer entropy |
| `rotamer_1body_energy0/1/2` | `(n_frame, n_res)` | float32 | Rotamer 1-body components |
| `rotamer_bad_solves_cumulative` | `(n_frame, 1)` | int64 | Rotamer solver failure count |
| `pivot_stats`, `jump_stats` | `(n_frame, 2)` | int32 | Monte Carlo moves [accepted, attempted] since the previous frame (`--monte-carlo-interval`) |
| `box` | `(n_frame, 3)` | float32 | Box dimensions in Angstrom (MARTINI/hybrid) |
| `volume` | `(n_frame, 1)` | float32 | Box volume (MARTINI/hybrid) |
| `curvature_center` | `(n_frame, 1, 3)` | float32 | Centre of curvature as `--curvature-changer-interval` moves it; a continuation copies the frame's value into `Const3D_curvature_center/value` |
| `replica_cumulative_swaps` | `(n_frame, n_partner, 2)` | int32 | Replica exchange with each swap partner [accepted, attempted] since the start of the run |
| `replica_swap_partner` | `(n_partner,)` | int32 | The systems this one exchanges with, written once |
| `tip_pos`, `time_estimate`, `tip_step` | `(n_frame, n_tip, 3)`, `(n_frame, 1)`, `(n_frame, 1)` | float32, float32, int64 | AFM tip positions, the tip's clock and its step count (`MovingConst3D`) |
| `sigmoid_coupling`, `nonlinear_coupling`, `contact_energy`, `centers`, `Const3D_<k>` | `(n_frame, n)` or `(n_frame, n, 3)` | float32 | Per-residue burial energy, contact energy, group centres and constant points (`--log-level detailed`, the default) |
| node group names (`environment_coverage_sc`, `placement_*`, ...), `virtual` | `(n_frame, n_elem[, width])` | float32 | Node outputs written at `--log-level extensive` |

The engine deletes an existing `/output` when it starts (unless `-o` names another file). A continuation script
therefore first renames it to `/output_previous_0` (then `_1`, `_2`, ...), as the examples' `continue_sim` do.
Every run logs the round it starts from, so the first frame of each later segment repeats the last frame of the
one before; `py/mdtraj_upside.py`, `py/extract_vtf.py` and `py/get_info_from_upside_traj.py` read the segments in
order and drop that frame, and `py/martini_extract_vtf.py` reads hybrid trajectories.

---

## 2. Force Field Parameter Files

### 2.1 `sidechain.h5`

Protein sidechain rotamer force field. Root-level datasets:

| Dataset | Shape | dtype | Notes |
|---|---|---|---|
| `bead_order` | `(20,)` | S5 | Residue type ordering, e.g. `[ALA_0, ARG_0, ...]` (20 residue types) |
| `restype_order` | `(20,)` | S3 | Three-letter AA ordering |
| `pair_interaction` | `(20,20,54)` | float64 | SC-SC 2D quadratic spline knots; axes: bead type i, bead type j (`bead_order`; one bead per residue in ff2.1), combined angular+radial knots |
| `coverage_interaction` | `(8,20,50)` | float64 | H-bond coverage interaction; axes: 2 x backbone class + (donor 0, acceptor 1), bead type, knots. Backbone classes: generic, GLY, PRO, and pre-proline (a residue whose chain continues into PRO or CPR; `backbone_group4` in `upside_config.py`) |
| `hydrophobe_interaction` | `(12,20,50)` | float64 | Hydrophobe coverage; axes: 3 x backbone class + hydrophobe point, bead type, knots |
| `hydrophobe_placement` | `(3,7)` | float64 | Fixed-point-vector-scalar placement of each residue's three backbone hydrophobe points, trained sites that started at N, CA and C; columns: xyz + direction_xyz + scalar |
| `rotamer_center_fixed` | `(86,6)` | float64 | Rotamer center positions (vec3 + unit_dir_vec3) per rotamer state |
| `rotamer_start_stop_bead` | `(20,3)` | int64 | [start, stop, n_bead] per residue type into rotamer list |
| `rotamer_prob` | `(36,36,86)` | float32 | Rotamer prior probability (phi_bin, psi_bin, rotamer_idx) |
| `restype_and_chi_and_state` | `(242,4)` | float64 | (restype_idx, chi1, chi2, state_idx) for each rotamer |
| `accuracy_trace` | `(n_iter,4)` | float32 | Training convergence log |

Knot layout for `pair_interaction[i,j,:]`: two angular sections (r1 direction, r2 direction) followed by two radial sections: one isotropic, in ff_2.1 the repulsive core (about +3.5 E_up at 0.7 A, its first coefficient
pinned at 5 by the trainer), and one weighted by the two angular factors, whose sign changes along r. `n_knot_angular=15`, `n_knot_sc=12` => 2×15 + 2×12 = 54 knots total.

`rotamer_center_fixed` and `rotamer_prob` were generated from Dunbrack's backbone-dependent rotamer library by
master's `src/compute_rotamer_centers.cpp` and `src/compute_rotamer_pos.cpp`; Upside 3 does not carry those tools.

**`ff_3.0/sidechain.h5`** is ff_2.1's file with `pair_interaction`, `coverage_interaction` and
`hydrophobe_interaction` retrained by ConDiv (the ff30_gdepth_dt009 release, 4). Every other
dataset is ff_2.1's unchanged, including `accuracy_trace` and the root attribute `program_text`, so
those two describe ff_2.1's rotamer fit, not ff3.0's training.

### 2.2 `hbond.h5`

Root dataset `parameter (12,)` float64, loaded into the `hbond_energy` node as its
`parameters` dataset. **These are not spline knots and not a function of H...O distance.** They
are three energies plus a bias, then four (boundary, sharpness) pairs in **radians** that classify
each residue by its (phi,psi):

| idx | name | ff_2.1 value | meaning |
|---|---|---|---|
| 0 | `E_alpha` | -1.9609 | per-hbond energy when phi < 0 and psi in (-120, 60) deg |
| 1 | `E_beta` | -1.9463 | per-hbond energy when phi < 0 and psi outside that window |
| 2 | `E_other` | -1.7690 | per-hbond energy when **phi in (0, 165) deg** |
| 3 | `E_bias` | -0.4058 | applied via a sigmoid on the per-peptide hbond count |
| 4,5 | `boundary_turn1`, `sharpness_turn1` | 0 deg, 3.8197 | phi window lower edge |
| 6,7 | `boundary_turn2`, `sharpness_turn2` | 165 deg, 3.8197 | phi window upper edge |
| 8,9 | `boundary_helix1`, `sharpness_helix1` | -120 deg, 3.8197 | psi window lower edge |
| 10,11 | `boundary_helix2`, `sharpness_helix2` | 60 deg, 3.8197 | psi window upper edge |

`Ehbond[i] = E_alpha*helix + E_beta*sheet + E_other*turn`, and the three scores sum exactly to 1.
A sharpness of 3.8197 gives a 15 deg ramp at each edge (`compact_sigmoid` saturates at |x| = 1/s).

Two consequences worth knowing. **The "turn" branch is exactly the positive-phi region**, so a
hydrogen bond at phi > 0 is worth **+0.192 E_up less** than the same bond at phi < 0 — in practice
a glycine penalty, since glycine is the residue that populates it. And the total potential is
**exactly linear in `parameters[:4]`** (verified to 7e-7 against the engine), which is what makes
`--hb-scale` and the trainer's `dE/ds = E/s` correct.

**`ff_3.0/hbond.h5`** has the same 12 entries: 0-3 retrained (-1.8150, -1.8294, -1.7335,
-0.5513), 4-11 as ff_2.1. It has no `class_restype`.

**Residue-class offsets.** A file may carry `parameter (12 + 3n,)` and a root dataset
`class_restype (n,)` S3 naming n residue types; entries `12+3(c-1) .. 12+3(c-1)+2` are offsets
added to `E_alpha`, `E_beta`, `E_other` for residues of class c. `upside_config` then writes the
node's `residue_class (n_res,)` int32: 1 + the residue's index in `class_restype`, or 0 for the
shared energies alone. No shipped file uses them. A 12-entry file has no `class_restype`, the
node no `residue_class`, and the engine's value, forces and parameter derivative are bitwise
those of master. The energy stays linear in the offsets, and `--hb-scale` scales them with
`parameters[:4]`.

### 2.3 `environment.h5`

Protein burial/environment potential. Root datasets:

| Dataset | Shape | Notes |
|---|---|---|
| `restype_order` | `(20,)` S3 | AA ordering |
| `center` | `(20,)` | Sigmoid center per residue type |
| `scale` | `(20,)` | Sigmoid scale; which side is favoured depends on its sign together with `sharpness` (PRO has both negative and burial favoured, -0.538 E_up buried in ff_2.1) |
| `sharpness` | `(20,)` | Sigmoid sharpness |
| `coverage_param` | `(20,1,4)` | Coverage function parameters per residue |
| `energies` | `(20,18)` | Energy table indexed by coverage bin; `@inv_dx=2.0`, `@offset=-0.5` |
| `weights` | `(400,)` | 20x20 inter-residue environment weight matrix |

**Trap.** In the `.up` these become `nonlinear_coupling_environment/coeff (20,18)` and
`weights (400,)`, but the engine's parameter vector for that node is **both concatenated,
coeff then weights, 360 + 400 = 760**. Requesting `coeff.shape` alone from `get_param_deriv`
fails with `Wrong number of parameters, expected 760 but got 360`. Request the full vector and
slice.

### 2.4 `bb_env.dat`

Plain-text, one float per line. Four values defining the backbone-environment
sigmoid coupling:

```
line 1: scale        (energy scale factor; ff_2.1 -0.3014, ff_3.0 -0.3300)
line 2: center       (zero-crossing of sigmoid; 2.0)
line 3: sharpness    (sigmoid width; 0.5)
line 4: hbond_weight (H-bond contribution weight; 2.0)
```

This is the order `write_bb_environment` reads (`bb_env_param[0]` is `scale`); an earlier version
of this table had the first two swapped. These parameterize `bb_sigmoid_coupling_environment`, the
FF2 backbone desolvation term. For each backbone NH/CO site it sums the side-chain coverage around
it, weighted by the covering residue's type, plus `hbond_weight` times its H-bond coverage, and
adds `scale * compact_sigmoid(sum - center, sharpness)`. The node holds **20 weights, one per residue
type**: `environment.h5`'s `weights` (2.3), which ff_2.1 and ff_3.0 share across every coupling
type (ff_3.0's 20 rows agree to 4e-6), so the writer stores the first row and refuses a library
whose rows differ by more than 1e-4. Configs
written before 2026-10-19 hold all 400 entries and are read by their first row, bitwise as before.
With `scale` negative (ff_2.1 -0.30) and `compact_sigmoid` equal to 1 at low burial, the term rewards **solvent-exposed** backbone NH/CO and is zero when they are buried: it stabilises the unfolded state. Measured on ubiquitin, -16.8 native against -43.8 expanded.

The node's parameters are `scale`, `center`, `sharpness`, `hbond_weight`, then the 20 weights, all
with exact derivatives (findings.md 3.14 A29). The FF2 trainer (`training/ConDiv.py`
`main_worker`) builds this node, trains `scale`, keeps center, sharpness and hbond weight at
ff_2.1's values, and adds the derivative of the weights to that of the side-chain burial's shared
weights. The retired FF1-form port never built it (findings.md 9t).

### 2.5 `membrane.h5`

Implicit membrane potential (ff_2.1). Root attributes `cb_z_min`, `cb_z_max` (-17, 15) and `hb_z_min`, `hb_z_max`
(-33, 15), in A relative to the membrane surface at +-half_thickness (negative is outside). The writer lays each
z axis on its n nodes with spacing (z_max - z_min)/(n - 2) and `z_start = z_min - half_thickness`; the clamped
spline is flat beyond its end nodes. Root datasets:

| Dataset | Shape | Notes |
|---|---|---|
| `names` | `(20,)` S3 | AA ordering |
| `burial_nodes` | `(1,)` | burial level (2.5) dividing the two burial bins of `cb_energy` |
| `cb_energy` | `(20,2,18)` | CB energy per residue type, burial bin, z node |
| `icb_energy` | `(20,18)` | CB energy of a channel residue that faces no lipid, per residue type, z node |
| `hb_energy` | `(2,2,18)` | backbone H/O energy per donor/acceptor, unpaired/paired, z node |
| `ihb_energy` | `(2,2,18)` | `hb_energy`'s counterpart for a channel residue that faces no lipid |

`icb_energy` and `ihb_energy` meet the clamped spline's end condition (zero slope at both end nodes); `cb_energy` and
`hb_energy` do not, so their force jumps where a CB or H/O crosses an end node's depth (pass 2, membrane M3; a
property of ff_2.1, left for the membrane retraining).

### 2.6 `martini.h5`

MARTINI non-bonded parameter tables. Two sub-groups:

#### `/particles`

Precomputed per-pair LJ + reaction-field Coulomb spline tables (reaction-field
parameters: epsilon_r=15, epsilon_rf=0; potential-shifted LJ; both go to zero
at 1.2 nm). Attrs: `@schema='martini_particles_combined'`, `@r_max_ang=12.0`,
`@r_min_ang=0.3`, `@n_points=1000`.

| Dataset | Shape | Notes |
|---|---|---|
| `unique_charge_product` | `(81,)` float64 | Distinct charge products for all pairs |
| `unique_eps_eup` | `(81,)` float64 | LJ epsilon per unique pair in E_up |
| `unique_sig_ang` | `(81,)` float64 | LJ sigma per unique pair in Angstrom |
| `combined_energy_grids` | `(81,1000)` float64 | Spline knot values; rows index unique pair type |

#### `/sc_table`

Side-chain/environment free-energy tables, one per (residue type, rotamer state, environment target), a target
being a bead type with its topology charge (`Q0+1` is the PC choline bead; the type name alone does not
fix the charge). Attrs document the derivation: `@schema='martini_sc_combined'`,
`@runtime_representation='clamped_cubic_spline_r_theta'`,
`@average='boltzmann_free_energy_over_conformers_and_azimuth'`, `@average_temperature_upside=0.8647`,
`@nonbonded_cutoff_nm=1.2`, `@bead_reach_nm`, `@target_direction_count=2400`,
`@sidechain_bead_frame_count=24`.

A rotamer state's side chain is a set of conformers: ideal-geometry heavy atoms (the bond lengths, angles and
dihedral offsets of master's `src/generate_from_rotamer.cpp`, with which the rotamer library was placed) built on
the reference N, CA, CB at each of the state's (chi1, chi2) rows of `sidechain.h5 restype_and_chi_and_state`, with
chi3 and chi4 at their staggered wells (-60, 60, 180; ARG's chi4 -85, 85, 180; the chi3 of GLU and GLN, which has
no wells, every 30 degrees), syn-pentane combinations of adjacent sp3-sp3 dihedrals left out, all equally weighted
(546 conformers over the 18 residues). Its beads are martinize's martini22 mapping centres, mass-weighted with
martinize's masses. Each cell is the free energy, at the production temperature 0.8647 T_up, of one target over the
state's conformers and the side chain's rotation about the CB vector:
exp(-F/kT) = sum_c w_c <exp(-E_c/kT)>, with E_c the dry-MARTINI energy (potential-shifted LJ and reaction-field
Coulomb) of conformer c's beads with the target and <> the average over 2400 Fibonacci-sphere directions, binned
linearly onto the theta node, and 24 rotations of the beads about the CB vector. F is the target's potential of
mean force in (r, theta): it gives a single target the Boltzmann distribution it has about the side chain with the
conformer and the azimuth integrated out. A rotation about the CB vector is a chi1 rotation, so the table does not
see chi1; states that differ only in chi1 get nearly the same table. r is the distance from the CB, placed as the
runtime places it (relative to the N/CA/C centroid); theta is the angle between the CB vector (CB - CA) and the
CB-to-target direction, 0 on the side-chain side. The radial grid runs from 0 to the largest bead-CB distance plus
1.2 nm and a node, so the last two nodes are exactly zero (the builder asserts it and the engine refuses a
table that is not). The engine fits each table as a clamped cubic spline in r and theta, exact for this
table at both ends: zero value and slope at the radial end, and an azimuthal average is even in theta about
0 and pi. Tables written before 2026-10-20 (cos(theta) grid, 2.5-12 A, one target per bead type) are
refused.

| Dataset | Shape | Notes |
|---|---|---|
| `restype_order` | `(18,)` S4 | Residue types with SC beads (no GLY/ALA) |
| `target_order` | `(13,)` S8 | Target labels, bead type and charge (`C1+0`, `Q0+1`, `Qa-1`, ...) |
| `target_type`, `target_charge` | `(13,)` S8, float32 | The same targets split |
| `grid_ang` | `(N_r,)` float32 | Radial grid in Angstrom, uniform from 0 |
| `theta_grid` | `(N_theta,)` float32 | Angular grid in radians, uniform on [0, pi] |
| `rotamer_count` | `(18,)` float32 | Number of rotamer states per residue type, as in `sidechain.h5` |
| `rotamer_full_energy_eup` | `(18,6,13,N_r,N_theta)` float32 | Free energy in E_up per (residue, rotamer, target, r, theta) |

### 2.7 `parameters/dryMARTINI/` (structure templates)

Holds no parameter table. `DOPC.pdb` is the dry-MARTINI DOPC slab that `example/16.MARTINI/run_sim_hybrid.sh`
passes as its default `--membrane-pdb`; `martini_prepare_system.py` itself requires the option and has no default.
The dry-MARTINI non-bonded parameters are `ff_2.1/martini.h5` (2.6).

### 2.8 `rama.dat`

HDF5 file with two top-level groups:

| Group | Contents |
|---|---|
| `/coil` | `dimer_pot (21,2,22,72,72)` float32, `dimer_weight (21,2,22)` float64. Axes: central_restype (21 including CPR), direction (left/right), neighbor_restype (22 including ALL), phi_bin (72), psi_bin (72). Source: Dunbrack NDRD library. |
| `/sheet` | `dimer_pot (20,2,20,72,72)` float32, `dimer_weight (20,2,20)` float32. Axes: central_restype (20, no CPR), direction, neighbor_restype, phi_bin, psi_bin. Source: PDB beta-sheet subset. |

The 72x72 grid covers phi in [-180, 175] and psi in [-180, 175] in 5-degree
steps. Values are -log(probability), so lower is more probable, **used directly as an energy in
E_up with no scale factor**. Maps are normalised at build time
(`pots -= -log(sum(exp(-pots)))`), an additive constant that cancels in any basin-to-basin
difference. About 4.5% of the coil array is NaN (unpopulated bins); the NaN pattern is itself
mirror-symmetric. The `@restype` attribute lists the residue ordering (the coil neighbour axis has
22 entries, the last being `ALL`, the marginal over neighbours); `@dir=[left, right]` means the
residue to the left/right of the central residue.

**Chirality operation.** `(phi,psi) -> (-phi,-psi)` on this grid is a reversal **with a roll**,
`np.roll(np.roll(m[..., ::-1, ::-1], 1, -2), 1, -1)`, because index `i` maps to `(-i) % 72`. A
plain `[::-1,::-1]` is off by one bin and wrong by ~2.3 in practice.

**Values are `-ln P`, normalised so that `sum(exp(-E)) = 1` over the grid.** That is exact: every
map in `rama.dat` sums to `1.000000`, and `mixture_potential` states the requirement in its
docstring. They are then used **directly as energies in E_up with no scale factor**, because
Upside is defined so `kT = 1 E_up` at `T_up = 1` (350.588 K).

**Writing an externally measured surface into a map therefore means dividing by `kT` at the
temperature it was measured, not by the E_up factor.** A GROMACS PMF in kJ/mol at 300 K enters as
`PMF / 2.494339`; `PMF / 2.914952774272` would be 1.169x too small. The two are easy to swap,
so check the normalisation of the map you are replacing: it is one line, and it settles the
question.

**The 4.5% NaN is one whole neighbour column, not scattered bins.** 4.545% = 1/22 is the `CPR`
column, which is never read: `read_rama_maps_and_weights` maps a cis-proline *neighbour* onto
`PRO`, keeping CPR only as a central residue. A naive `abs(a-b) > tol` comparison across it reports
no difference, because NaN comparisons are False.

**Where each library lives:**

| file | what it is |
|---|---|
| `parameters/common/rama.dat` | ff2.1, as published. Shared, because ff2.1 never trained it |
| `parameters/ff_3.0/rama.dat` | ff3.0's library, trained with the rest of ff3.0 (ff30_gdepth_dt009, 4): `common/rama.dat` with the **37 central-glycine coil maps (`GLY\|left\|X`, `GLY\|right\|X`) shifted by one trained depth pair**, alpha_R -0.1068 and alpha_L +0.5169 (findings 1.26-1.27). The sheet group, every `dimer_weight`, every other map and the pooled `GLY\|ALL` entry are NDRD's |
| `parameters/common/rama31.dat` | The fixed library of the `ff30_bio` training (round 4), not ff3.0's: ff2.1's with the **central-glycine row, coil and sheet, fitted in Upside to BioEmu's plain-MD octapeptides** (amber ff99sb-ildn, 300 K; Zenodo 15641199): iterative Boltzmann inversion at T_up = 1 (the library convention, below) with ff2.1's other terms, on glycines at residues 2-5; pooled `GLY\|X` fitted on X-G-Y (249 glycines), `GLY\|left\|GLY` on G-G-Y (27), `GLY\|right\|GLY` on X-G-G (25), `GLY\|right\|PRO` on X-G-P (10). Fitted where the all-atom data have >= 10 frames per cell; elsewhere NDRD's own coil top with its barrier height above alpha_R kept (Upside never designed the transition region, findings 1.19). Glycine's coil and sheet entries are the same map. Each entry is stored with `rama_map_pot_ref` subtracted. In Upside on the octapeptides it reproduces BioEmu's X-G-Y glycine within 0.01 per basin (ln(aR/aL) -0.55 vs -0.535). Built 2026-10-05 (IBI pass-6 map; scripts and data in `/project/trsosnic/yinhan/checks/gly_bioemu_map/`); root attribute `glycine_row` records its origin. The previous content (AWH capped dipeptides, round 3's library, 2026-10-01) is in git history and `backup/rama31.dat.bak_pre_bioemu_20261005` |

`rama31.dat` differs from `common/rama.dat` only in the central-GLY row of the coil and sheet
groups. **Its GLY entries assume the reference correction**: every config path applies
`--reference-state-rama` (ConDiv, the benchmark, the MARTINI hybrid, the examples), and a config
built without it would give glycine the measured surface minus that correction. The same GLY entry
in coil and sheet makes the coil/sheet mixture leave glycine unchanged; the left/right mixture is
unchanged, so a glycine between a glycine and another residue reads a mixture of two entries.

**A force field's library belongs in its own directory.** A release copies the run's library
into `parameters/<ff>/rama.dat`; `ConDiv.py extract` copies it unchanged, and ff3.0's glycine depth pair was
trained by that run's own `trainer/rama_basin.py` (findings 1.27; not in the repo). A config on ff_3.0 takes
`ff_3.0/rama.dat`; `common/rama.dat` would give it ff2.1's glycine maps.

### 2.8a How a library map becomes the `rama_pot` the engine reads

`rama_pot` in a config is **one map per residue**, not a library map, and it is four
transformations away from the file. Anything differentiating or perturbing the library has to go
through all four, in this order (`read_weighted_maps`, then `write_rama_map_pot`):

1. **Left/right mixture.** `-log(w_l e^{-V_l} + w_r e^{-V_r})` over the two neighbour maps, with
   `w` from `dimer_weight`. Terminal residues use one direction only.
2. **Normalisation.** `pots += log(sum(exp(-pots)))`, applied to the coil and sheet mixtures
   separately, so each is back to `sum(exp(-E)) = 1`.
3. **Coil/sheet mixture.** The same weighted log-sum-exp, with the sheet side weighted by
   `sheet_weight * exp(-sheet_mixing[restype])`.
4. **A Boltzmann-weighted shift:** `rama_pot -= (rama_pot*exp(-rama_pot)).sum()` per map.

**Step 4 is the one that catches people.** It is a constant per map, so it changes no force and
cancels in every basin-to-basin difference, which is why it is invisible in almost all analysis.
But it is a *map-dependent* constant and it enters the total potential, so anything computing
`dE/d(library map)` must include it. Omitting it reproduced the per-residue map with a constant
7.23 offset and got the symmetric-part gradient wrong by 37%.

Two more traps in the same function. `dimer_weight` for glycine is **not** 1.0 (0.908 and 0.948
for the two directions in one measured case); printing it with `precision=0` rounds it to 1 and
invites exactly the wrong simplification. And `read_rama_maps_and_weights` declares `pots` as
`f4`, so the whole chain is float32-limited: a finite difference on a library map cannot resolve
better than about `5e-7` per cell against map values of order 10.

### 2.9 `sheet` (plain text)

One float per line, **20 values: one sheet-mixing energy per residue type**, ordered by the
`@restype` attribute of `rama.dat`'s `/sheet` group. Not spline knots.
`write_rama_map_pot` asserts `size == len(sheet_restype)`, so a file of any other length is
rejected. Each residue's coil and sheet maps are blended as
`mixture_potential([coil_w, sheet_w*exp(-sheet_mixing)], [coil_pot, sheet_pot])`.

Sweeping a residue's value across its whole range barely moves the alpha_R/alpha_L balance,
because the sheet endpoint is empty in both helical basins; it mainly adds beta. Under
`--rama-param-deriv` the config also gains `more_/less_sheet_rama_pot_<TYPE>` pairs at
`eps = 5e-4` for finite-difference training, one pair for each residue type in the sequence
(`PRO` also for a `CPR`), each map set with only that type's value shifted. The engine does not
read them; ConDiv does.

### 2.10 `rama_reference.pkl`

Python 2 pickle (load with `encoding='latin1'`) of one `(72,72)` float64 array: a probability
over the phi/psi grid of 2.8, summing to 1, the same for every residue type (built from a free
alanine chain; Jumper thesis 4.3). `upside_config.py` takes its log, subtracts the mean and writes
the result as the single map of `rama_map_pot_ref`, which every residue reads (1.3).

---

## 3. Hybrid Mapping File (`hybrid_mapping.h5`)

Located in `hybrid_prep/hybrid_mapping.h5` alongside a run directory.
Built by `martini_prepare_system.py`. Contains:

- `/protein_beads`: atom name, residue index, bead type, coordinates for the
  protein pseudo-atoms (N, CA, C, O) that become MARTINI BB beads.
- `/lipid_beads`: MARTINI lipid bead coordinates and types from the initial
  bilayer `.up` or full-atom PDB.
- `/bb_map`: virtual backbone bead map (same content as `input/hybrid_bb_map`
  inside the run `.up`).

---

## 4. Relationship Between Files

```
upside_config.py                   martini_prepare_system.py
   |                                  |
   |-- reads: fasta, initial.npy      |-- reads: martini_prepare_params,
   |-- reads: rama.dat                |           martini.h5, sidechain.h5,
   |-- reads: ff_2.1/{hbond,          |           membrane PDB
   |           sidechain, environment,|
   |           membrane}.h5           |-- writes: hybrid_mapping.h5
   |-- reads: bb_env.dat              |
   |-- reads: sheet, rama_ref.pkl     |-- writes: <protein>.stage_N.up
   |                                       containing all nodes below
   |-- writes: <protein>.up
       /input/args  (CLI metadata)
       /input/pos   (initial coordinates)
       /input/potential/
           backbone geometry nodes
           rama nodes
           H-bond nodes
           sidechain rotamer nodes
           environment nodes
           [MARTINI nodes if hybrid]
       /input/pivot_moves  (MC samplers)
```

```
Force-field files (read-only, never modified at runtime):
  parameters/ff_2.1/
    sidechain.h5   -> rotamer + environment potentials
    hbond.h5       -> H-bond energy
    environment.h5 -> burial/environment potential
    membrane.h5    -> implicit membrane potential
    bb_env.dat     -> backbone-environment coupling (4 floats)
    sheet          -> sheet-mixing energy, one value per residue type (20 floats)
  parameters/common/
    rama.dat       -> Ramachandran library (HDF5)
    rama_reference.pkl -> reference state
  parameters/dryMARTINI/
    DOPC.pdb       -> dry-MARTINI DOPC slab template (2.7)
  parameters/ff_2.1/
    martini.h5     -> MARTINI non-bonded spline tables
  parameters/ff_3.0/
    sidechain.h5   -> trained pair, coverage and hydrophobe interaction tables; the rest ff_2.1's
    environment.h5 -> trained sigmoid burial scale/center/sharpness and weights (type 1)
    hbond.h5       -> trained H-bond energies and bias (entries 0-3); entries 4-11 as ff_2.1;
                      12 entries, no class_restype
    sheet          -> trained per-type sheet mixing energies
    bb_env.dat     -> trained backbone-term scale; center/sharpness/hbond weight as ff_2.1
    rama.dat       -> common/rama.dat with glycine's coil maps shifted by the trained depth pair
                      (2.8); use this, not common/rama.dat
```

**`parameters/ff_3.0` is the release of ff30_gdepth_dt009** (step 75, `epoch_03_minibatch_18`;
converged at its 2026-10-09 gate, findings 1.26), trained from ff2.1 with ff2.1's own FF2
dual-target workflow (`training/ConDiv.py`; the glycine depth pair by the run's `trainer/rama_basin.py`) and
md5-identical to `$P/parameters/ff_3.0_gdepth`. Its `hbond.h5` has ff2.1's 12 entries, no glycine offsets.
The two earlier versions, trained with a port of FF1's trainer, are in `backup/` (gitignored) and
in git history (findings 9t-9v). Build configs from all six files with the default
`--environment-potential-type=1`, and take `membrane.h5` and `martini.h5` from `ff_2.1`. The
examples take their parameter files as arguments, so a protein-only run on ff3.0 passes them explicitly; the hybrid
preparation hard-codes ff_2.1's files and `common/rama.dat` (`martini_prepare_system.py`, `martini_prepare_system_lib.py`),
so it cannot build an ff3.0 hybrid yet.

---

## 5. Atom Indexing and Coordinate Layout

The position array `pos` has shape `(N, 3, n_frame)`. In the input file `n_frame=1`.
The last axis is the frame placeholder; during simulation the engine grows it.

A protein-only configuration holds three atoms per residue:
- index `3*i`: N of residue i
- index `3*i+1`: CA of residue i
- index `3*i+2`: C of residue i

A hybrid configuration holds four backbone atoms per residue, then one MARTINI BB bead per residue, then the
environment (lipids, then ions):
- index `4*i` to `4*i+3`: N, CA, C, O of residue i (`hybrid_bb_map/atom_indices`)
- index `4*n_res + i`: BB bead of residue i (`hybrid_bb_map/bb_atom_index`)

Only N, CA and C are dynamic. `martini_hybrid_position` overwrites the O and BB slots at every evaluation (O
from `infer_H_O`, BB as in `hybrid_bb_map`), and both are held as fixed atoms, so they carry no momentum and
their entries in `pos` keep the input's values; the logged `pos` holds the placed ones. Side chains and CB are not in `pos` in either layout;
the `placement_*` nodes derive them.

---

## 6. Key Conventions and Invariants

**Angle storage**: bond angles in `Spring_angle` are stored as cosines of the
equilibrium angle, not degrees. A tetrahedral angle (109.5 deg) is stored as
cos(109.5 deg) = -0.3338.

**Omega equilibrium**: `Spring_omega` stores the equilibrium at pi (trans),
with `@pbc=1` and `@box_len=2*pi` so the periodic spring wraps correctly.

**GLY Ramachandran**: glycine handedness belongs in the library file, never in
`write_rama_map_pot`. There is no symmetrization step in that function and none should be added.
Mirror symmetry is exact only for a glycine in an achiral environment; a glycine flanked by
L-amino acids sits in a chiral one. ff3.0's library (`parameters/ff_3.0/rama.dat`, 2.8) is NDRD's
with one trained depth pair added to glycine's coil maps. `parameters/common/rama31.dat` (2.8) holds
glycine entries fitted in Upside to unselected all-atom peptides (BioEmu plain MD), fixed and never
trained: the pooled `GLY|X` and, each on its own context, `GLY|left|GLY` (G-G-Y), `GLY|right|GLY`
(X-G-G) and `GLY|right|PRO` (X-G-P). The two `GLY|GLY` entries carry an L neighbour on the other
side, so they are not mirror-symmetric; the data say so (findings 1.19). Glycine's coil and sheet
entries are the same map, so the coil/sheet mixture leaves glycine unchanged. In both libraries
every non-glycine map is the NDRD map unchanged, and the config writer stays ignorant of all of
it. See `GLY_sym.md` (the earlier AWH build) and `training/README.md`.

**Spline tables must reproduce the analytic potential exactly**: verify the
tables against the analytic form (dry-MARTINI: reaction-field Coulomb with
epsilon_r=15, epsilon_rf=0, plus potential-shifted LJ) after unit conversion.
Do not tabulate a bare truncation or smoothed approximation.

**Thermostat decoupling**: `thermostat_timescale` for frozen atoms is set to 1e8
(effectively infinite), meaning those atoms are not thermostated and stay rigid.
Protein atoms use ~5 Upside time units.

**Pair exclusion (`sequence_id`)**: every pair list, `backbone_pairs` and each interaction graph (`protein_hbond`,
`hbond_coverage`, `hbond_coverage_hydrophobe`, `environment_coverage_sc/_hb`, `hb_environment_coverage_hn/_oc`,
`hbbb_coverage`, `rotamer/pair_interaction`), excludes element pairs by a required per-element `sequence_id`
(`sequence_id1`, `sequence_id2` for the two groups of an asymmetric graph); the engine refuses a node without it. `id`
is never the key: it is the node's own payload (the residue that indexes `aa_types`, a packed rotamer id). `sequence_id`
is the element's residue number plus 2 for every chain break before it (`upside_config.sequence_id`), so residues on
the two sides of a break are at least 3 apart, beyond every node's window: pairs within 2 are excluded by the coverage
nodes, within 1 by `protein_hbond` and `backbone_pairs`, and equal ones (the same residue) by the rotamer. In a single
chain `sequence_id` is the residue number.

**NP ensemble barostat**: `@target_p_xy=0` with `@mc_dmax_z=0` gives a tensionless
XY-only MC barostat. The z dimension is frozen. The MC step size (`mc_dmax_xy`) is
tuned to give ~1% box changes per trial.

**NP dt during unfolding**: the NP (protein unfolding) timestep must remain at
0.001 regardless of MARTINI LJ stability. Unfolding drives backbone bonds to
large-amplitude oscillations that require a short dt for numerical accuracy.
