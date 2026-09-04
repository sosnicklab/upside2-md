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
| `mom` | `(N, 3, 1)` | float32 | Initial momenta; zeroed at start |
| `mass` | `(N,)` | float32 | Per-atom mass in m_up |
| `atom_names` | `(N,)` | `S4` | Atom names matching PDB convention (N, CA, C, CGL, PO4, ...) |
| `atom_roles` | `(N,)` | `S4` | Role tag: BB, SC1, SC2, W, CGL, etc. |
| `particle_class` | `(N,)` | `S10` | One of: PROTEIN, LIPID, WATER, ION, OTHER |
| `charges` | `(N,)` | float64 | Partial charge in e |
| `molecule_ids` | `(N,)` | int32 | Per-atom molecule index (0-based, contiguous) |

Additional datasets present in protein or hybrid runs:

| Dataset | Shape | dtype | Notes |
|---|---|---|---|
| `sequence` | `(n_res,)` | `S3` | Three-letter AA codes (ASP, GLU, ...) |
| `type` | `(N,)` | `S8` | MARTINI bead type (Qd, P5, Na, ...) or atom element |
| `vel` | `(N, 3)` | float32 | Velocities for MARTINI integrator |
| `thermostat_timescale` | `(N,)` | float32 | Per-atom Langevin friction timescale; protein ~5, frozen beads ~1e8 |
| `residue_ids` | `(N,)` | int32 | 1-based residue index per atom |

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
an `@arguments` attribute listing which upstream node outputs it consumes, and
`@integrator_level` (0 = evaluated every step, 1 = evaluated at backbone level).

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
| `Angle` | `id (M,3)` | Atom triplets for angle calculation |
| `Dihedral_phi` | `id (n,4)` float64; negative -1 = missing N-terminus | Phi dihedral atom quads |
| `Dihedral_psi` | `id (n,4)` float64 | Psi dihedral atom quads |
| `Dihedral_omega` | `id (n,4)` int64 | Omega dihedral atom quads |
| `Spring_bond` | `id (M,)`, `equil_dist (M,)`, `spring_const (M,)` | `@pbc=0`, `@dim1=0` |
| `Spring_angle` | `id (M,)`, `equil_dist (M,)` (cosine), `spring_const (M,)` | Angles stored as cosines |
| `Spring_omega` | `id (M,)`, `equil_dist (M,)` (~pi), `spring_const (M,)` | `@pbc=1`, `@box_len=2pi` |
| `affine_alignment` | `atoms (n_res,3)`, `ref_geom (n_res,3,3)` | N/CA/C local frame per residue |
| `backbone_pairs` | `id (n_res,)`, `n_atom (n_res,)`, `ref_pos (n_res,4,3)` | Backbone H-bond geometry |

#### Rama / secondary structure nodes

| Node | Key datasets/attrs | Notes |
|---|---|---|
| `rama_coord` | (no datasets; just `@arguments=[pos]`) | Computes phi/psi pair |
| `rama_map_pot` | `residue_id`, `rama_map_id`, `rama_pot (n_res,72,72)` | 5-deg grid, -log(prob); `@restype` string, `@sheet_eps` |
| `rama_map_pot_ref` | `residue_id`, `rama_map_id_all`, `rama_pot (n_res,72,72)` | Reference-state correction; `@log_pot` |

#### Sidechain rotamer nodes

| Node | Key datasets | Notes |
|---|---|---|
| `placement_fixed_point_vector_only` | `affine_residue (n_sc,)`, `layer_index`, `placement_data (n_sc,6)` | SC1 bead position (vec3 + unit vec3 direction) |
| `placement_fixed_point_vector_only_CB` | same | CB position |
| `placement_fixed_point_vector_scalar` | `placement_data (n_sc,7)` | Hydrophobe position (vec3 + dir + scalar) |
| `placement_scalar` | (scalar output per SC) | 1-body rotamer energy coverage |
| `rotamer` | `@arguments`, `@damping=0.4`, `@max_iter`, `@tol`; sub-group `pair_interaction` | Self-consistent mean-field rotamer solver |
| `rotamer/pair_interaction` | (loaded from `sidechain.h5`) | SC-SC interaction table |

#### H-bond and environment nodes

| Node | Key datasets | Notes |
|---|---|---|
| `infer_H_O` | sub-groups `donors` and `acceptors` with `residue`, `bond_length`, `id` | Backbone H and O position inference |
| `protein_hbond` | `id1`, `id2`, `type1`, `type2` | Donor-acceptor index pairs |
| `hbond_energy` | `@integrator_level=1` | Reads `hbond.h5` 12-knot spline |
| `hbond_coverage` | `index_pos`, `index_weight` | HN-atom coverage accumulator |
| `hbond_coverage_hydrophobe` | `index_pos`, `index_weight` | OC-atom hydrophobe coverage |
| `environment_coverage_sc` | `@num_aa_types=20` | SC environment coverage |
| `environment_coverage_hb` | same | HB environment coverage |
| `sigmoid_coupling_environment` | `@integrator_level=1`, `@number_independent_weights` | SC-environment coupling (ff_2.1) |
| `bb_sigmoid_coupling_environment` | `@center`, `@scale`, `@sharpness`, `@hbond_weight` | BB-environment coupling |

#### MARTINI-specific nodes (pure lipid or hybrid)

| Node | Key datasets/attrs | Notes |
|---|---|---|
| `martini_potential` | `@potential_type`, `@cutoff`, `@n_types`, `@n_points=1000`, `@optimized_format` | Precomputed LJ+Coulomb spline table; per-bead-type-pair energy grid |
| `martini_sc_table_1body` | `@energy_conversion_kj_per_eup`, `@length_conversion_angstrom_per_nm` | 1-body SC-bilayer table offset |
| `compose_vector6d` | `direction (N,3)`, `elem_index (N,)`, `orientation_index (N,)`, `display_head_offset_ang`, `display_tail_offset_ang` | Coarse-grained lipid orientational state; attrs document DOPC geometry derivation |
| `cgl_orientation_state` | `@rotational_inertia_up`, `@rotational_thermostat_timescale`, `@conformer_count` | CGL rigid-body orientational integrator |
| `cgl_compaction_state` | `@mass_up`, `@thermostat_timescale`, `@self_coord_min_ang` | CGL compaction coordinate integrator |
| `cg_lipid_pair` | `pair_interaction/id`, `/index`, `/type`, `/interaction_param`, `/reference_energy_eup` | CGL-CGL pair PMF; `@n_radial`, `@n_angular=9`, `@knot_spacing_ang=0.35`, `@schema='cg_lipid_pair_full'` |
| `cg_lipid_rotamer_sc` | `pair_interaction/interaction_param (n_sctypes,1,n_r*n_ang)`, `/delta_compact`, `/delta_compressed` | SC-CGL interaction; 3 compaction states; `@schema='cg_lipid_sc_full_tensor'` |
| `cg_lipid_target` | `pair_interaction` | BB-CGL interaction |
| `cg_lipid_target_base` | `pair_interaction` | BB-CGL base (no compaction overlay) |
| `cg_lipid_compaction_self` | `self_coeff (n_knot,)`, `@self_n_knot`, `@self_coord_min_ang` | CGL self-compaction PMF |
| `dist_spring` | `bonded_atoms (N,)`, `equil_dist (N,)`, `spring_const (N,)` | Intra-lipid bond springs |

#### Hybrid protein+MARTINI extra nodes

| Node | Notes |
|---|---|
| `affine_alignment` | Protein-only backbone frame (same structure as pure protein) |
| `backbone_pairs` | Protein H-bond geometry |
| `hbond_energy` | Protein H-bond energy from `hbond.h5` |
| `infer_H_O` | Protein H/O inference |
| `protein_hbond` | Protein donor-acceptor pairs |
| `rama_map_pot / _ref` | Protein Ramachandran |
| `rotamer` | Protein rotamer |
| `cg_lipid_*` | Lipid potentials (same as pure lipid) |
| `martini_potential` | MARTINI non-bonded (all bead types) |
| `martini_sc_table_1body` | SC-bilayer 1-body offset |

### 1.4 `/input` groups (MARTINI and hybrid only)

| Group | Key attrs/datasets | Notes |
|---|---|---|
| `barostat` | `@type=0`, `@target_p_xy=0`, `@target_p_z=0`, `@compressibility=14.52` (3e-4 bar^-1 in upside units), `@compressibility_z=0` (frozen z), `@semi_isotropic=1`, `@tau_p=4.0`, `@interval=10` | Monte-Carlo semi-isotropic barostat; XY only for membrane |
| `cgl_gle` | `coupling (n_mode,)`, `memory_tau (n_mode,)`, `coupling_scale (7,n_mode)`, `memory_tau_scale (7,n_mode)`, `temperature_grid (7,)`, `atom_index (N_cgl,)`, `aux_momentum (n_mode,N_cgl,3)` | Generalized Langevin (sum-of-exponential) thermostat for CGL beads; `@schema='cgl_exponential_memory_gle'` |
| `hybrid_remap` | `old_to_new (N_allbead,)`, `@n_old`, `@n_new` | Maps full MARTINI bead set to reduced simulation set |
| `hybrid_bb_map` | `atom_indices (n_res,4)`, `bb_atom_index (n_res,)`, `weights (n_res,4)`, `bb_type (n_res,)`, `bb_chain_id`, `bb_resseq`, `bb_secondary_structure`, `bb_comment`, `reference_atom_coords`, `reference_atom_indices` | Virtual backbone bead placement from protein N/CA/C/O; one row per protein residue; `@virtual_backbone_com_mode=1` |
| `hybrid_control` | `@activation_stage`, `@exclude_intra_protein_martini=1`, `@preprod_protein_mode`, `@sc_env_backbone_hold_steps`, `@sc_env_po4_z_hold_steps` | Hybrid staging parameters |
| `hybrid_env_topology` | `env_atom_indices (N_env,)`, `protein_membership (N_prot,)` | Which atoms are in the MARTINI environment layer |
| `stage_parameters` | `@current_stage=b'minimization'` or `b'production'`; sub-groups `minimization_angles`, `minimization_bonds`, `production_angles`, `production_bonds`, each containing `force_constants` | Force constant schedule across stages |

### 1.5 `/input` restart state (hybrid)

| Dataset | Notes |
|---|---|
| `cgl_compaction_mom (N_cgl,)` | CGL compaction momenta from last frame; `@restart_source='output/cgl_compaction_mom[-1]'` |
| `cgl_orientation_mom (N_cgl,3)` | CGL orientational momenta |
| `cgl_gle/aux_momentum (n_mode,N_cgl,3)` | GLE auxiliary momenta |

### 1.6 `/input` Monte Carlo samplers (protein-only)

| Group | Key datasets | Notes |
|---|---|---|
| `pivot_moves` | `pivot_atom (n,5)`, `pivot_range (n,2)`, `pivot_restype (n,)`, `proposal_pot (n_restype,72,72)` | MC pivot move table |

### 1.7 `/output` datasets

The engine appends one row per saved frame to each extendable dataset.

| Dataset | Shape | dtype | Notes |
|---|---|---|---|
| `pos` | `(n_frame, n_replica, N, 3)` | float32 | Coordinates per frame |
| `mom` | `(n_frame, n_replica, N, 3)` | float32 | Momenta per frame |
| `potential` | `(n_frame, n_replica)` | float64 | Total potential energy in E_up |
| `kinetic` | `(n_frame, n_replica)` | float64 | Total kinetic energy in E_up |
| `temperature` | `(n_frame, n_replica)` | float64 | Instantaneous temperature in T_up |
| `time` | `(n_frame,)` | float64 | Simulation time |
| `hbond` | `(n_frame, n_hbond)` | float32 | Per H-bond occupancy |
| `rama_map_potential` | `(n_frame, n_res)` | float32 | Per-residue Ramachandran energy |
| `rotamer_potential` | `(n_frame, n_res)` | float32 | Per-residue rotamer pair energy |
| `rotamer_free_energy` | `(n_frame, n_res)` | float32 | Per-residue rotamer free energy |
| `rotamer_entropy` | `(n_frame, n_res)` | float32 | Per-residue rotamer entropy |
| `rotamer_1body_energy0/1/2` | `(n_frame, n_res)` | float32 | Rotamer 1-body components |
| `rotamer_bad_solves_cumulative` | `(n_frame, 1)` | int64 | Rotamer solver failure count |
| `box` | `(n_frame, 3)` | float32 | Box dimensions in Angstrom (MARTINI/hybrid) |
| `volume` | `(n_frame, 1)` | float32 | Box volume (MARTINI/hybrid) |
| `pressure` | `(n_frame, 2)` | float32 | Lateral and normal pressure (MARTINI/hybrid) |
| `cgl_orientation` | `(n_frame, n_replica, N_cgl, 3)` | float32 | CGL orientational unit vector |
| `cgl_orientation_mom` | `(n_frame, n_replica, N_cgl, 3)` | float32 | CGL orientational momenta |
| `cgl_compaction` | `(n_frame, n_replica, N_cgl)` | float32 | CGL compaction coordinate |
| `cgl_compaction_mom` | `(n_frame, n_replica, N_cgl)` | float32 | CGL compaction momenta |
| `cgl_gle_aux_momentum` | `(n_frame, n_mode, N_cgl, 3)` | float32 | GLE aux momenta |

After a continuation run, the previous `/output` group is renamed to
`/output_previous_0` (then `_1`, `_2`, ...) and a fresh `/output` is created.
`py/run_upside.py:read_output` concatenates all of them.

---

## 2. Force Field Parameter Files

### 2.1 `sidechain.h5`

Protein sidechain rotamer force field. Root-level datasets:

| Dataset | Shape | dtype | Notes |
|---|---|---|---|
| `bead_order` | `(20,)` | S5 | Residue type ordering, e.g. `[ALA_0, ARG_0, ...]` (20 residue types) |
| `restype_order` | `(20,)` | S3 | Three-letter AA ordering |
| `pair_interaction` | `(20,20,54)` | float64 | SC-SC 2D quadratic spline knots; axes: restype_i, restype_j, combined angular+radial knots |
| `coverage_interaction` | `(8,20,50)` | float64 | H-bond coverage interaction; axes: coverage_type, restype, knots |
| `hydrophobe_interaction` | `(12,20,50)` | float64 | Hydrophobe coverage; axes: hyd_type, restype, knots |
| `hydrophobe_placement` | `(3,7)` | float64 | Fixed-point-vector-scalar placement for H/N/C; columns: xyz + direction_xyz + scalar |
| `rotamer_center_fixed` | `(86,6)` | float64 | Rotamer center positions (vec3 + unit_dir_vec3) per rotamer state |
| `rotamer_start_stop_bead` | `(20,3)` | int64 | [start, stop, n_bead] per residue type into rotamer list |
| `rotamer_prob` | `(36,36,86)` | float32 | Rotamer prior probability (phi_bin, psi_bin, rotamer_idx) |
| `restype_and_chi_and_state` | `(242,4)` | float64 | (restype_idx, chi1, chi2, state_idx) for each rotamer |
| `accuracy_trace` | `(n_iter,4)` | float32 | Training convergence log |

Knot layout for `pair_interaction[i,j,:]`: two angular sections (r1 direction, r2 direction) followed by two radial sections (attractive, repulsive). `n_knot_angular=15`, `n_knot_sc=12` => 2×15 + 2×12 = 54 knots total.

### 2.2 `hbond.h5`

Root dataset `parameter (12,)` float64: 12 cubic spline knots parameterizing
the backbone H-bond energy as a function of H...O distance. Loaded directly
into the `hbond_energy` node.

### 2.3 `environment.h5`

Protein burial/environment potential. Root datasets:

| Dataset | Shape | Notes |
|---|---|---|
| `restype_order` | `(20,)` S3 | AA ordering |
| `center` | `(20,)` | Sigmoid center per residue type |
| `scale` | `(20,)` | Sigmoid scale (positive = buried-favored) |
| `sharpness` | `(20,)` | Sigmoid sharpness |
| `coverage_param` | `(20,1,4)` | Coverage function parameters per residue |
| `energies` | `(20,18)` | Energy table indexed by coverage bin; `@inv_dx=2.0`, `@offset=-0.5` |
| `weights` | `(400,)` | 20x20 inter-residue environment weight matrix |

### 2.4 `bb_env.dat`

Plain-text, one float per line. Four values defining the backbone-environment
sigmoid coupling:

```
line 1: center    (zero-crossing of sigmoid)
line 2: scale     (energy scale factor)
line 3: sharpness (sigmoid width)
line 4: hbond_weight (H-bond contribution weight)
```

These parameterize `bb_sigmoid_coupling_environment`. The sign convention is
that large positive coverage (buried) shifts the backbone energy down.

### 2.5 `membrane.h5`

Implicit membrane potential (ff_2.1). Root datasets:

| Dataset | Shape | Notes |
|---|---|---|
| `names` | `(20,)` S3 | AA ordering |
| `burial_nodes` | `(1,)` | z-coordinate node for membrane burial calculation |
| `cb_energy` | `(20,2,18)` | CB burial energy per residue type, leaflet, z-bin |
| `icb_energy` | `(20,18)` | Isotropic CB energy |
| `hb_energy` | `(2,2,18)` | H-bond energy per donor/acceptor type, leaflet, z-bin |
| `ihb_energy` | `(2,2,18)` | Isotropic H-bond energy |

### 2.6 `martini.h5`

MARTINI non-bonded parameter tables. Two sub-groups:

#### `/particles`

Precomputed per-pair LJ + reaction-field Coulomb spline tables (reaction-field
parameters: epsilon_r=15, epsilon_rf=0; potential-shifted LJ; both go to zero
at 1.2 nm). Attrs: `@schema='martini_particles_combined'`, `@r_max_ang=12.0`,
`@r_min_ang=0.3`, `@n_points=1000`.

| Dataset | Shape | Notes |
|---|---|---|
| `type_order` | `(38,)` S8 | MARTINI bead type names (Qda, Qd, Qa, P5, ...) |
| `type_charge` | `(38,)` float32 | Charge per bead type |
| `unique_charge_product` | `(81,)` float64 | Distinct charge products for all pairs |
| `unique_eps_eup` | `(81,)` float64 | LJ epsilon per unique pair in E_up |
| `unique_sig_ang` | `(81,)` float64 | LJ sigma per unique pair in Angstrom |
| `combined_energy_grids` | `(81,1000)` float64 | Spline knot values; rows index unique pair type |

#### `/sc_table`

Azimuthally-averaged SC-bead free energy tables (SC-particle pairs).
Attrs document the derivation: `@schema='martini_sc_combined'`,
`@spline_control_quantity='direct_rotamer_free_energy_kj_mol'`,
`@azimuthal_average='tempered_boltzmann_free_energy'`,
`@azimuthal_average_temperature_upside=25.0`,
`@nonbonded_cutoff_nm=1.2`.

| Dataset | Shape | Notes |
|---|---|---|
| `restype_order` | `(18,)` S4 | Residue types with SC beads (no GLY/ALA) |
| `target_order` | `(38,)` S8 | MARTINI bead types |
| `cos_theta_grid` | `(13,)` float32 | Angular grid points in cos(theta) |
| `grid_ang` | `(96,)` float32 | Radial grid in Angstrom |
| `radial_energy_eup` | `(18,38,96)` float32 | Isotropic radial free energy in E_up |
| `angular_energy_eup` | `(18,38,96)` float32 | Isotropic angular-collapsed energy |
| `angular_profile` | `(18,38,13)` float32 | Angular spline profile per (SC type, bead type) |
| `rotamer_count` | `(18,)` float32 | Number of rotamers per residue type |
| `rotamer_probability_fixed` | `(18,6)` float32 | Prior rotamer probabilities |
| `rotamer_radial_energy_eup` | `(18,6,38,96)` float32 | Per-rotamer radial energy |
| `rotamer_angular_energy_eup` | `(18,6,38,96)` float32 | Per-rotamer angular energy |
| `rotamer_angular_profile` | `(18,6,38,13)` float32 | Per-rotamer angular profile |
| `rotamer_full_energy_eup` | `(18,6,38,96,13)` float32 | Full 2D (radial x angular) per rotamer |

### 2.7 `dopc.h5` (dryMARTINI CGL parameters)

Located in `parameters/dryMARTINI/`. Contains CGL-CGL pair PMFs and
SC-CGL interaction tables built from dry-MARTINI DOPC simulations.
The top-level structure mirrors what is written into `/input/potential/cg_lipid_*`
nodes in a simulation `.up` file. Built by `py/martini_build_tables.py`
and `py/martini_gen_params.py`.

### 2.8 `rama.dat`

HDF5 file with two top-level groups:

| Group | Contents |
|---|---|
| `/coil` | `dimer_pot (21,2,22,72,72)` float32, `dimer_weight (21,2,22)` float64. Axes: central_restype (21 including CPR), direction (left/right), neighbor_restype (22 including ALL), phi_bin (72), psi_bin (72). Source: Dunbrack NDRD library. |
| `/sheet` | `dimer_pot (20,2,20,72,72)` float32, `dimer_weight (20,2,20)` float32. Axes: central_restype (20, no CPR), direction, neighbor_restype, phi_bin, psi_bin. Source: PDB beta-sheet subset. |

The 72x72 grid covers phi in [-180, 175] and psi in [-180, 175] in 5-degree
steps. Values are -log(probability), so lower is more probable. The `@restype`
attribute lists the residue ordering; `@dir=[left, right]` means the residue
to the left/right of the central residue.

### 2.9 `sheet` (plain text)

One float per line. Contains 63 values that are the spline knots of the
sheet-content mixing energy used by `rama_map_pot` to interpolate between
coil and sheet Ramachandran maps. Written as:

```
line 1..63: cubic spline knots for sheet mixing potential
```

### 2.10 `rama_reference.pkl`

Python pickle file containing a dict `{restype_str: log_prob_array}` mapping
each residue type to its 72x72 reference-state log-probability map. Used to
compute the reference-state correction in `rama_map_pot_ref`.

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
   |-- reads: rama.dat                |           dopc.h5, martini.h5
   |-- reads: ff_2.1/{hbond,          |           sidechain.h5
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
    sheet          -> sheet-mixing energy knots (63 floats)
  parameters/common/
    rama.dat       -> Ramachandran library (HDF5)
    rama_reference.pkl -> reference state
  parameters/dryMARTINI/
    dopc.h5        -> CGL pair PMFs and SC-CGL interaction tables
  parameters/ff_2.1/
    martini.h5     -> MARTINI non-bonded spline tables
```

---

## 5. Atom Indexing and Coordinate Layout

The position array `pos` has shape `(N, 3, n_frame)`. In the input file `n_frame=1`.
The last axis is the frame placeholder; during simulation the engine grows it.

Protein atoms are ordered N, CA, C within each residue:
- index `3*i`: N of residue i
- index `3*i+1`: CA of residue i
- index `3*i+2`: C of residue i

Sidechain pseudo-atoms (CB, SC1, ...) are placed at the end of the protein block.
MARTINI lipid beads follow all protein atoms. CGL orientation unit-vectors are
stored in a separate parallel array (`compose_vector6d/direction`) not in `pos`.

The `hybrid_remap` group stores an `old_to_new (N_all,)` index map from the full
all-atom (or all-bead) MARTINI index space to the reduced simulation index space
used by Upside. Atoms mapped to the same index are merged.

---

## 6. Key Conventions and Invariants

**Angle storage**: bond angles in `Spring_angle` are stored as cosines of the
equilibrium angle, not degrees. A tetrahedral angle (109.5 deg) is stored as
cos(109.5 deg) = -0.3338.

**Omega equilibrium**: `Spring_omega` stores the equilibrium at pi (trans),
with `@pbc=1` and `@box_len=2*pi` so the periodic spring wraps correctly.

**CGL target angle convention**: `cg_lipid_target` uses the angular convention
`ang = -n_cgl_dot_n12` (dot product of CGL orientation with the inter-bead
unit vector, negated). The SC-CGL table uses `ang1=-n1_dot_n12; ang2=n2_dot_n12`.
Do NOT negate the `cos_theta_grid` in `_build_cgl_target_table`; the sign is
already embedded in the convention string.

**GLY Ramachandran**: GLY maps must be symmetrized over phi for both the L and D
wings in `write_rama_map_pot`. The symmetrization must apply unconditionally for
all GLY residues; any conditional on phi range or alphaR>alphaL will silently
miss residues and destabilize TM helices.

**Spline tables must reproduce the analytic potential exactly**: verify the
tables against the analytic form (dry-MARTINI: reaction-field Coulomb with
epsilon_r=15, epsilon_rf=0, plus potential-shifted LJ) after unit conversion.
Do not tabulate a bare truncation or smoothed approximation.

**Thermostat decoupling**: `thermostat_timescale` for frozen atoms is set to 1e8
(effectively infinite), meaning those atoms are not thermostated and stay rigid.
Protein atoms use ~5 Upside time units.

**NP ensemble barostat**: `@semi_isotropic=1` with `@compressibility_z=0` gives
a tensionless XY-only MC barostat. The z dimension is frozen. The MC step size
(`mc_dmax_xy`) is tuned to give ~1% box changes per trial.

**NP dt during unfolding**: the NP (protein unfolding) timestep must remain at
0.001 regardless of MARTINI LJ stability. Unfolding drives backbone bonds to
large-amplitude oscillations that require a short dt for numerical accuracy.
