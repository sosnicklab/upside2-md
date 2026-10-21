#!/usr/bin/env python

import argparse
import contextlib
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time
import zlib
from collections import Counter
from pathlib import Path

import numpy as np

from martini_itp_reader import lipid_max_sigma_nm, parse_dry_forcefield, parse_itp_file, parse_lipid_from_itp
from martini_prepare_system_lib import (
    BB_TYPE_CHARGE,
    PROTEIN_RESIDUE_NAMES,
    STAGE_SOFT_CORE_LAMBDA,
    active_lipid_names,
    lipid_diffusion_targets,
    build_backbone_with_virtual_bb,
    center_of_mass,
    convert_stage,
    compute_lipid_residue_indices,
    coords,
    estimate_salt_pairs,
    extract_protein_backbone_atoms_from_aa,
    inject_particles_table,
    inject_stage7_sc_table_nodes,
    lipid_resname,
    measure_belt_hydrophobic_coverage,
    membrane_slab_excluded_volume,
    set_active_lipid_name,
    map_backbone_types_from_structure,
    parse_pdb,
    place_ions,
    reject_inside_slab,
    remove_overlapping_lipids,
    set_initial_position,
    set_box_from_lipid_xy,
    set_coords,
    read_charmm_gui_membrane,
    tile_and_crop_bilayer_lipids,
    validate_hybrid_mapping,
    write_hybrid_mapping_h5,
    write_pdb,
)


PY_DIR = Path(__file__).resolve().parent
REPO_ROOT = PY_DIR.parent

DEFAULT_SC_ENV_PO4_Z_HOLD_STEPS = 150
DEFAULT_NPT_INTERVAL = 10
# Monte-Carlo barostat lateral step: exp(+-mc_dmax_xy) per trial (~1% box change), tuned for a
# reasonable acceptance rate on a dry-MARTINI membrane.
DEFAULT_NPT_MC_DMAX_XY = 0.01
DEFAULT_MARTINI_ENERGY_CONVERSION = 2.914952774272
DEFAULT_MARTINI_LENGTH_CONVERSION = 10.0
MARTINI_MD_INTEGRATOR = "v"
MARTINI_MD_TIME_STEP = 0.009
# Mean amino-acid excluded volume (partial specific volume 0.73 cm^3/g, <MW> ~110 Da) used to
# subtract the protein from the solvent-accessible volume when counting ions.
PROTEIN_RESIDUE_VOLUME_A3 = 130.0
# CHARMM-GUI's Martini bilayer equilibration (step6.2-6.6): the phosphate height restraint
# BILAYER_LIPIDHEAD_FC in kJ/mol/nm^2, released over five stages.
EQUILIBRATION_STAGES = (
    ("6.2", 200.0),
    ("6.3", 100.0),
    ("6.4", 50.0),
    ("6.5", 20.0),
    ("6.6", 10.0),
)
# The stage-7.0 solvation gate (findings 8b). A belt site faces lipid when fewer than BELT_FACING_MAX_BACKBONE
# protein backbone atoms lie within BELT_BURIAL_RADIUS_A of it. Averaged over a burn-in's last quarter, healthy
# productions leave 0.4-3.3% of those sites without an environment bead within BELT_VACUUM_RADIUS_A (1-7.5% without
# an acyl tail), a relaxed glpG prep 4.9% (6.5%), and the insertion void the gate guards 18% (25%).
BELT_BURIAL_RADIUS_A = 10.0
BELT_FACING_MAX_BACKBONE = 60
BELT_VACUUM_RADIUS_A = 8.0
BELT_VACUUM_MAX_FRACTION = 0.10
# An MD start is refused if a pair has an energy the equilibrium ensemble would show in fewer than one
# configuration in this many (assert_no_core_overlap).
CORE_OVERLAP_ODDS = 1.0e9


def find_lipid_itp_path(ff_dir: Path, lipid_name: str) -> Path:
    """Return the ITP under *ff_dir* that defines the *lipid_name* moleculetype."""
    target = lipid_name.strip().upper()
    for itp in sorted(Path(ff_dir).glob("*.itp")):
        if any(mol.upper() == target for mol in parse_itp_file(itp).get("molecules", {})):
            return itp
    raise FileNotFoundError(
        f"Lipid molecule '{lipid_name}' not found in any ITP under '{ff_dir}'."
    )


def derive_lipid_contact_clearance_angstrom(upside_home: Path, lipid_name: str) -> float:
    ff_dir = Path(upside_home) / "example" / "16.MARTINI" / "dryMARTINI_itp"
    dry_ff_path = ff_dir / "dry_martini_v2.1.itp"
    _atomtypes, pair_params = parse_dry_forcefield(dry_ff_path)
    # For a mixed bilayer take the largest per-lipid clearance so the packing cutoff clears
    # the bulkiest bead of any lipid type.
    clearance = 0.0
    for name in (part.strip() for part in str(lipid_name).split(",") if part.strip()):
        lipid = parse_lipid_from_itp(find_lipid_itp_path(ff_dir, name), name)
        max_sigma_nm = lipid_max_sigma_nm(lipid["bead_types"], pair_params)
        clearance = max(clearance, (2.0 ** (1.0 / 6.0)) * max_sigma_nm * DEFAULT_MARTINI_LENGTH_CONVERSION)
    return float(clearance)


def _apolar_bead_indices(lipid: dict) -> list:
    """0-based indices of the truly hydrophobic beads (MARTINI apolar class C1-C5)."""
    return [i for i, t in enumerate(lipid["bead_types"]) if re.fullmatch(r"C[1-5]", str(t).strip())]


def _count_hydrophobic_tails(lipid: dict) -> int:
    """Number of separate acyl chains, from the apolar subgraph of the molecule's bonds."""
    apolar = set(_apolar_bead_indices(lipid))
    adjacency = {i: set() for i in apolar}
    for bond in lipid["bonds"]:
        # parse_lipid_from_itp already returns 0-based bead indices.
        i, j = int(bond[0]), int(bond[1])
        if i in apolar and j in apolar:
            adjacency[i].add(j)
            adjacency[j].add(i)
    seen = set()
    tails = 0
    for start in sorted(apolar):
        if start in seen:
            continue
        tails += 1
        stack = [start]
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            stack.extend(adjacency[node] - seen)
    return tails


def require_bilayer_lipids(upside_home: Path, lipid_name: str):
    """Refuse single-tailed amphiphiles, read from the ITP topology (findings 7.1).

    A one-chain amphiphile with a bulky head has a packing parameter far below 1, so it forms micelles,
    not lamellae: forcing it into a slab pins the area per molecule at the head's cross-section and
    starves the hydrophobic core (DDM: ~14 A against a ~28 A transmembrane belt). Only bilayers are
    built, so a lipid whose acyl-chain count is one is rejected; deciding by tail count rather than a
    name table keeps the rule universal.
    """
    ff_dir = Path(upside_home) / "example" / "16.MARTINI" / "dryMARTINI_itp"
    tails = {}
    for name in (part.strip() for part in str(lipid_name).split(",") if part.strip()):
        lipid = parse_lipid_from_itp(find_lipid_itp_path(ff_dir, name), name)
        tails[name.upper()] = _count_hydrophobic_tails(lipid)
    if not tails:
        raise ValueError("Lipid name must be provided (--lipid-name)")
    single = sorted(n for n, t in tails.items() if t <= 1)
    if single:
        raise ValueError(
            f"{single} are single-tailed amphiphiles, which form micelles, not bilayers; only bilayer-forming "
            f"lipids are supported. Tails per molecule: {tails}."
        )


def derive_lipid_apolar_bead_names(upside_home: Path, lipid_name: str) -> set:
    """Bead NAMES of the hydrophobic beads, for tail-vs-head tests on packed coordinates."""
    ff_dir = Path(upside_home) / "example" / "16.MARTINI" / "dryMARTINI_itp"
    names = set()
    for name in (part.strip() for part in str(lipid_name).split(",") if part.strip()):
        lipid = parse_lipid_from_itp(find_lipid_itp_path(ff_dir, name), name)
        names.update(lipid["atom_names"][i].strip() for i in _apolar_bead_indices(lipid))
    return names


def read_opm_hydrophobic_half_thickness(opm_reference_pdb: Path) -> float:
    """Half of the hydrophobic thickness from the OPM reference's own REMARK.

    OPM publishes this per entry ("1/2 of bilayer thickness: 14.1" for the rhomboid fold, i.e. a 28.2 A
    belt, which matches glpG's measured 28.2 A span). Taking it from the reference keeps the belt
    definition a published measurement rather than a tunable input.
    """
    pattern = re.compile(r"1/2\s+of\s+bilayer\s+thickness:\s*([0-9.]+)", re.IGNORECASE)
    for line in Path(opm_reference_pdb).read_text().splitlines():
        if not line.startswith("REMARK"):
            continue
        match = pattern.search(line)
        if match:
            return float(match.group(1))
    raise ValueError(
        f"No '1/2 of bilayer thickness' REMARK in {opm_reference_pdb}. A micelle build needs the "
        "protein's hydrophobic belt; use an OPM-oriented reference PDB (opm.phar.umich.edu)."
    )


def derive_lipid_net_charges(upside_home: Path, lipid_name: str) -> dict:
    """Per-lipid net charge (electrons) keyed by upper-case residue name, read from the ITPs.

    Anionic lipids such as POPG carry a net charge that counterions must neutralize; neutral
    lipids (DOPC, DDM, POPE) return 0, leaving the single-lipid ion count unchanged.
    """
    ff_dir = Path(upside_home) / "example" / "16.MARTINI" / "dryMARTINI_itp"
    charges = {}
    for name in (part.strip() for part in str(lipid_name).split(",") if part.strip()):
        lipid = parse_lipid_from_itp(find_lipid_itp_path(ff_dir, name), name)
        charges[name.upper()] = int(round(sum(lipid["bead_charges"])))
    return charges


def parse_prepare_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Pack a protein into a CHARMM-GUI bilayer and convert it to an Upside hybrid stage input."
    )
    parser.add_argument("--pdb-id", required=True, help="Runtime PDB id for stage conversion")
    parser.add_argument("--runtime-pdb-output", required=True)
    parser.add_argument("--prepare-structure", type=int, default=1, choices=[0, 1])
    parser.add_argument("--stage", default=None, choices=sorted(STAGE_SOFT_CORE_LAMBDA),
                        help="pair tables of the UPSIDE input conversion")
    parser.add_argument("--run-dir", required=True)

    parser.add_argument("--lipid-name", required=True,
                        help="Membrane lipid moleculetype name(s), comma-separated for a mixed bilayer "
                             "(e.g. DOPC or POPE,POPG); topology read from the input ITP files")
    parser.add_argument("--membrane-pdb", default=None,
                        help="CHARMM-GUI membrane CG PDB (gromacs/step5_charmm2gmx.pdb); must carry a CRYST1 box")
    parser.add_argument("--membrane-top", default=None,
                        help="CHARMM-GUI system.top for the lipid-count cross-check")
    parser.add_argument("--protein-cg-pdb", default=None,
                        help="martinize MARTINI22 CG PDB of the protein (BB+SC envelope; prep-time geometric aid only)")
    parser.add_argument("--opm-reference", default=None,
                        help="OPM membrane-oriented reference PDB for --protein-orientation-mode opm")
    parser.add_argument("--protein-aa-pdb", default=None)
    parser.add_argument("--hybrid-mapping-output", default=None)
    parser.add_argument("--hybrid-bb-map-json-output", default=None)

    parser.add_argument("--xy-scale", type=float, default=1.0)
    parser.add_argument("--box-padding-xy", type=float, default=0.0)
    parser.add_argument("--box-padding-z", type=float, default=0.0)
    parser.add_argument("--membrane-thickness-angstrom", type=float, default=0.0,
                        help="Equilibrated membrane thickness for the ion count; 0 = measure from the "
                             "packed lipids (use the equilibrated dry-MARTINI value so production is 0.15 M)")
    parser.add_argument("--salt-molar", type=float, default=0.15)
    parser.add_argument("--explicit-ions", type=int, choices=[0, 1], default=1)
    parser.add_argument("--ion-cutoff", type=float, default=4.0)
    parser.add_argument("--seed", type=int, default=2026)

    parser.add_argument("--protein-lipid-cutoff", type=float, default=3.0)
    parser.add_argument(
        "--protein-orientation-mode",
        choices=["input", "lay-flat", "opm"],
        default="input",
        help="input: keep the input orientation; lay-flat: PCA thinnest axis -> bilayer normal; "
             "opm: superpose onto an OPM membrane-oriented reference (--opm-reference).",
    )
    parser.add_argument(
        "--protein-pbc-margin",
        type=float,
        default=15.0,
        help="Minimum solvent/lipid belt (Angstrom) added around the protein on each side so it "
             "does not interact with its periodic image (>= the dry-MARTINI nonbonded cutoff).",
    )
    parser.add_argument("--summary-json", default=None)
    return parser.parse_args(argv)


def copy_if_different(src: Path, dst: Path):
    src_resolved = src.expanduser().resolve()
    dst_resolved = dst.expanduser().resolve()
    if src_resolved == dst_resolved:
        return
    shutil.copy2(src_resolved, dst_resolved)


def canonicalize_axis_sign(axis):
    axis = np.asarray(axis, dtype=float)
    dominant = int(np.argmax(np.abs(axis)))
    if axis[dominant] < 0.0:
        axis = -axis
    return axis


def lay_flat_rotation_basis(protein_xyz):
    centered = protein_xyz - center_of_mass(protein_xyz)
    if centered.shape[0] < 3:
        return np.eye(3, dtype=float)

    cov = np.cov(centered.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]

    axis_x = canonicalize_axis_sign(eigvecs[:, order[0]])
    axis_z = canonicalize_axis_sign(eigvecs[:, order[-1]])
    axis_y = np.cross(axis_z, axis_x)
    axis_y_norm = float(np.linalg.norm(axis_y))
    if axis_y_norm <= 1.0e-8:
        axis_y = canonicalize_axis_sign(eigvecs[:, order[1]])
    else:
        axis_y = axis_y / axis_y_norm

    axis_z = np.cross(axis_x, axis_y)
    axis_z = axis_z / np.linalg.norm(axis_z)
    basis = np.column_stack((axis_x, axis_y, axis_z))
    if np.linalg.det(basis) < 0.0:
        basis[:, 1] *= -1.0
    return basis


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H", "CYX": "C", "MSE": "M",
}


def _ordered_ca_seq_xyz(atoms):
    """Return (one-letter sequence, CA xyz array) in residue order for standard-AA CA atoms."""
    seq, xyz, seen = [], [], set()
    for atom in atoms:
        if atom["name"].strip().upper() != "CA":
            continue
        one = THREE_TO_ONE.get(atom["resname"].strip().upper())
        if one is None:
            continue
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key in seen:
            continue
        seen.add(key)
        seq.append(one)
        xyz.append([atom["x"], atom["y"], atom["z"]])
    return "".join(seq), np.asarray(xyz, dtype=float)


def _kabsch(P, Q):
    """Optimal rigid transform (R, t) minimizing ||R.P_i + t - Q_i||; apply as xyz @ R.T + t."""
    Pc, Qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - Pc).T @ (Q - Qc)
    U, _s, Vt = np.linalg.svd(H)
    d = 1.0 if np.linalg.det(Vt.T @ U.T) >= 0.0 else -1.0
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Qc - R @ Pc


def orient_protein_to_opm(protein_backbone_atoms, opm_reference_pdb):
    """Rigid transform aligning the protein onto an OPM membrane-oriented reference.

    Superposes matched Cα atoms (BLOSUM62 global alignment of the shared fold, then Kabsch with
    iterative <=4 A core refinement so flexible loops don't skew the TM-core rotation). The OPM
    reference has its membrane midplane at z=0, so the returned transform places the protein at its
    physical membrane depth/orientation. Returns (R, t); apply as xyz @ R.T + t.
    """
    from Bio import Align
    from Bio.Align import substitution_matrices

    ref_atoms, _ = parse_pdb(Path(opm_reference_pdb))
    ref_seq, ref_xyz = _ordered_ca_seq_xyz(ref_atoms)
    our_seq, our_xyz = _ordered_ca_seq_xyz(protein_backbone_atoms)
    if our_xyz.shape[0] < 3 or ref_xyz.shape[0] < 3:
        raise ValueError("Too few Cα atoms for OPM alignment.")

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"
    alignment = aligner.align(our_seq, ref_seq)[0]

    our_idx, ref_idx = [], []
    for (s1, e1), (s2, e2) in zip(alignment.aligned[0], alignment.aligned[1]):
        for k in range(e1 - s1):
            our_idx.append(s1 + k)
            ref_idx.append(s2 + k)
    if len(our_idx) < 3:
        raise ValueError("OPM sequence alignment produced too few matched residues.")

    P, Q = our_xyz[our_idx], ref_xyz[ref_idx]
    R, t = _kabsch(P, Q)
    for _ in range(5):
        mask = np.linalg.norm(P @ R.T + t - Q, axis=1) < 4.0
        if int(mask.sum()) < 3:
            break
        R, t = _kabsch(P[mask], Q[mask])
    n_core = int((np.linalg.norm(P @ R.T + t - Q, axis=1) < 4.0).sum())
    print(
        f"OPM orientation: {len(our_idx)} aligned Cα, {n_core} within 4 A after Kabsch "
        f"(reference {Path(opm_reference_pdb).name})"
    )
    return R, t


def apply_rigid_transform(atoms, R, t):
    set_coords(atoms, coords(atoms) @ np.asarray(R, dtype=float).T + np.asarray(t, dtype=float))


def compute_protein_orientation(protein_backbone_atoms, args):
    """Rigid transform (R, t) to orient the protein for membrane insertion (apply as xyz @ R.T + t).

    opm: superpose onto the OPM reference (physical depth/orientation).
    lay-flat: PCA thinnest axis -> z, centered at the origin.
    input: keep orientation, centered at the origin.
    """
    mode = args.protein_orientation_mode
    if mode == "opm":
        if not args.opm_reference:
            raise ValueError("--protein-orientation-mode opm requires --opm-reference")
        return orient_protein_to_opm(
            protein_backbone_atoms, Path(args.opm_reference).expanduser().resolve()
        )
    com = center_of_mass(coords(protein_backbone_atoms))
    if mode == "lay-flat":
        R = lay_flat_rotation_basis(coords(protein_backbone_atoms)).T
    elif mode == "input":
        R = np.eye(3, dtype=float)
    else:
        raise ValueError(f"Unsupported protein orientation mode: {mode}")
    return R, -R @ com


def martinize_protein_cg(venv_python, martinize_script, protein_aa_pdb, out_dir):
    """Coarse-grain the all-atom protein to a MARTINI22 CG PDB (BB+SC) via martinize.py.

    The CG beads are used ONLY as a prep-time geometric envelope for insertion + overlap removal;
    the runtime protein stays on the Upside force field. All-coil SS is fine (bead coordinates are
    center-of-mass mappings, independent of the SS assignment)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cg_pdb = out_dir / "protein_cg.pdb"
    cg_top = out_dir / "protein_cg.top"
    run_checked(
        [str(venv_python), str(martinize_script), "-f", str(protein_aa_pdb),
         "-x", str(cg_pdb), "-o", str(cg_top), "-ff", "martini22"],
        cwd=out_dir,
    )
    if not cg_pdb.exists():
        raise FileNotFoundError(f"martinize did not produce a CG PDB: {cg_pdb}")
    return cg_pdb


def pack_bilayer_environment(args, lipid_template, bilayer_box, protein_backbone_atoms, cg_atoms, prot_xyz, cg_xyz):
    """Insert the protein into a periodic lamellar slab by deleting the lipids it overlaps."""
    # Slide the bilayer template in z so its midplane meets the protein's membrane center (opm: z=0
    # by construction; otherwise embed the protein COM in the slab).
    lip_xyz = coords(lipid_template)
    lip_mid_z = 0.5 * float(lip_xyz[:, 2].min() + lip_xyz[:, 2].max())
    membrane_center_z = 0.0 if args.protein_orientation_mode == "opm" else float(prot_xyz[:, 2].mean())
    for atom in lipid_template:
        atom["z"] = float(atom["z"] + (membrane_center_z - lip_mid_z))

    # xy fill window: the protein CG-envelope's larger footprint + a PBC-safe belt on each side, on both
    # axes because the protein turns freely in the membrane plane. Each axis is then rounded up to a whole
    # number of ITS OWN template tiles, so a rectangular tile also repeats seamlessly across both periodic
    # boundaries (findings 7.2).
    fmin = cg_xyz.min(axis=0)
    fmax = cg_xyz.max(axis=0)
    pspan = fmax - fmin
    pcenter_xy = 0.5 * (fmin[:2] + fmax[:2])
    belt = float(args.box_padding_xy) + float(args.protein_pbc_margin)
    tile_xy = np.array([float(bilayer_box[0]), float(bilayer_box[1])])
    fill_side = float(np.max(pspan[:2] + 2.0 * belt)) * max(1.0, float(args.xy_scale))
    target_xy = np.maximum(np.ceil(fill_side / tile_xy), 1.0) * tile_xy
    bilayer_lipids = tile_and_crop_bilayer_lipids(
        bilayer_atoms=lipid_template,
        bilayer_box=bilayer_box,
        target_xy_min=pcenter_xy - 0.5 * target_xy,
        target_xy_max=pcenter_xy + 0.5 * target_xy,
    )

    # Overlap removal against the CG envelope (BB+SC), not the backbone: whole lipids clashing with
    # the protein's true girth (including where the sidechain rotamers are rebuilt) are deleted.
    bilayer_kept, removed_lipids = remove_overlapping_lipids(
        bilayer_atoms=bilayer_lipids,
        protein_atoms=cg_atoms,
        lipid_residues=compute_lipid_residue_indices(bilayer_lipids),
        cutoff=float(args.protein_lipid_cutoff),
    )

    # Measure belt coverage BEFORE the box shift, while z=0 is still the membrane midplane.
    coverage = measure_belt_hydrophobic_coverage(
        protein_backbone_atoms=protein_backbone_atoms,
        env_atoms=bilayer_kept,
        apolar_bead_names=derive_lipid_apolar_bead_names(
            Path(os.environ.get("UPSIDE_HOME", str(REPO_ROOT))), args.lipid_name
        ),
        belt_half_thickness=(
            read_opm_hydrophobic_half_thickness(Path(args.opm_reference).expanduser().resolve())
            if args.opm_reference
            else 0.5 * float(np.ptp(coords(bilayer_kept)[:, 2]))
        ),
    )

    packed_atoms = protein_backbone_atoms + bilayer_kept
    box_lengths = set_box_from_lipid_xy(
        all_atoms=packed_atoms,
        lipid_atoms=bilayer_kept,
        xy_box=target_xy,
        pad_z=float(args.box_padding_z),
        min_box_z=float(pspan[2] + 2.0 * (float(args.box_padding_z) + float(args.protein_pbc_margin))),
    )

    kept_lipid_xyz = coords(bilayer_kept)
    membrane_z_lo = float(kept_lipid_xyz[:, 2].min())
    membrane_z_hi = float(kept_lipid_xyz[:, 2].max())
    packed_thickness_z = membrane_z_hi - membrane_z_lo
    # Count ions against the EQUILIBRATED membrane thickness (the dry-MARTINI membrane relaxes from
    # the compressed CHARMM-GUI start), so the realized concentration is 0.15 M in production. Ions
    # are still placed in the packed-time solvent slabs, outside the current membrane extent.
    count_thickness_z = (
        float(args.membrane_thickness_angstrom)
        if float(args.membrane_thickness_angstrom) > 0.0
        else packed_thickness_z
    )
    # The slab already excludes the protein's transmembrane part; only residues outside it take solvent
    # volume from the ions.
    membrane_mid_z = 0.5 * (membrane_z_lo + membrane_z_hi)
    ca_z = np.array([float(a["z"]) for a in protein_backbone_atoms if a["name"].strip() == "CA"])
    n_residues_outside_slab = int(np.count_nonzero(np.abs(ca_z - membrane_mid_z) > 0.5 * count_thickness_z))
    # The belt is the protein's own hydrophobic band: the OPM reference publishes it, otherwise fall
    # back to the packed slab's extent.
    belt_half_thickness = (
        read_opm_hydrophobic_half_thickness(Path(args.opm_reference).expanduser().resolve())
        if args.opm_reference
        else 0.5 * packed_thickness_z
    )
    geometry = {
        "lipid_residues_removed": int(removed_lipids),
        "belt_half_thickness_a": float(belt_half_thickness),
        "belt_coverage": coverage,
        "base_xy_tile_angstrom": [float(v) for v in tile_xy],
        "target_xy_angstrom": [float(v) for v in target_xy],
        "membrane_thickness_z_angstrom": packed_thickness_z,
        "ion_count_thickness_z_angstrom": count_thickness_z,
        "protein_span_angstrom": [float(v) for v in pspan],
        "protein_residues_outside_ion_slab": n_residues_outside_slab,
        "excluded_volume_a3": membrane_slab_excluded_volume(box_lengths, count_thickness_z)
                              + n_residues_outside_slab * PROTEIN_RESIDUE_VOLUME_A3,
        "membrane_z_range": (membrane_z_lo, membrane_z_hi),
    }
    return bilayer_kept, box_lengths, geometry


def prepare_mixed_structure(args, runtime_pdb):
    if not args.protein_aa_pdb:
        raise ValueError("--protein-aa-pdb is required")
    if not args.membrane_pdb:
        raise ValueError("--membrane-pdb is required")
    if not args.protein_cg_pdb:
        raise ValueError("--protein-cg-pdb (martinize CG envelope) is required")

    protein_aa_pdb = Path(args.protein_aa_pdb).expanduser().resolve()
    membrane_pdb = Path(args.membrane_pdb).expanduser().resolve()
    membrane_top = Path(args.membrane_top).expanduser().resolve() if args.membrane_top else None
    protein_cg_pdb = Path(args.protein_cg_pdb).expanduser().resolve()
    for required in (protein_aa_pdb, membrane_pdb, protein_cg_pdb):
        if not required.exists():
            raise FileNotFoundError(required)

    protein_aa_atoms, _ = parse_pdb(protein_aa_pdb)
    if not protein_aa_atoms:
        raise ValueError(f"No atoms found in protein AA PDB: {protein_aa_pdb}")
    protein_backbone_atoms, _ = extract_protein_backbone_atoms_from_aa(protein_aa_atoms)
    cg_atoms, _ = parse_pdb(protein_cg_pdb)
    if not cg_atoms:
        raise ValueError(f"No atoms found in martinize CG PDB: {protein_cg_pdb}")
    lipid_template, bilayer_box = read_charmm_gui_membrane(membrane_pdb, membrane_top)

    # Orient the protein as a rigid body; apply the SAME transform to the runtime backbone AND the
    # martinize CG envelope (both share the input-PDB frame) so insertion + overlap removal stay
    # geometrically consistent.
    R, t = compute_protein_orientation(protein_backbone_atoms, args)
    apply_rigid_transform(protein_backbone_atoms, R, t)
    apply_rigid_transform(cg_atoms, R, t)

    upside_home = Path(os.environ.get("UPSIDE_HOME", str(REPO_ROOT)))
    require_bilayer_lipids(upside_home, args.lipid_name)
    bilayer_kept, box_lengths, env_geometry = pack_bilayer_environment(
        args, lipid_template, bilayer_box, protein_backbone_atoms, cg_atoms,
        coords(protein_backbone_atoms), coords(cg_atoms)
    )
    packed_atoms = protein_backbone_atoms + bilayer_kept

    # Gate G1: the protein's hydrophobic belt must sit against acyl tails, not headgroups. A lamellar
    # DDM slab fails this (~14 A tail span, half the belt facing maltose), which is what unfolded TM4.
    # Reported, not gated: CHARMM-GUI templates are pre-minimization and laterally compressed, so a
    # packed-state span reads low for every lipid (DOPC 20 A, POPE/POPG 21 A) and cannot be compared with
    # OPM's relaxed hydrophobic thickness. assert_environment_solvation() does the real check on the lipid-
    # facing belt sites over the rigid burn-in's last quarter, which is both local and like-for-like.
    coverage = env_geometry.pop("belt_coverage")
    print(
        f"Hydrophobic coverage (packed state): tail span {coverage['tail_hydrophobic_span_a']:.1f} A, "
        f"belt {coverage['belt_span_a']:.1f} A, {coverage['belt_polar_contact_fraction']:.2f} of "
        f"{coverage['belt_residues']} belt residues nearest a polar bead"
    )

    bb_type_by_residue, bb_secondary_by_residue, bb_type_source_by_residue = \
        map_backbone_types_from_structure(protein_backbone_atoms)
    # Counterions neutralise the charges the topology actually carries. On the protein those are the BB
    # beads' MARTINI charges (the charged chain ends); its side chains carry no explicit charge by design,
    # their environment interaction being the side-chain/lipid table (findings 5.3c).
    protein_charge = int(round(sum(BB_TYPE_CHARGE.get(t, 0.0) for t in bb_type_by_residue.values())))
    lipid_net_charges = derive_lipid_net_charges(upside_home, args.lipid_name)
    membrane_charge = 0
    seen_lipid_molecules = set()
    for atom in bilayer_kept:
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key in seen_lipid_molecules:
            continue
        seen_lipid_molecules.add(key)
        membrane_charge += lipid_net_charges[atom["resname"].upper()]
    system_charge = protein_charge + membrane_charge

    # Ions occupy the solvent-accessible volume only (implicit solvent): the box minus the membrane slab
    # and the protein outside it, and never inside the hydrophobic interior.
    if int(args.explicit_ions):
        salt_pairs = estimate_salt_pairs(
            box_lengths=box_lengths,
            salt_molar=float(args.salt_molar),
            excluded_volume_a3=env_geometry["excluded_volume_a3"],
        )
        n_na = int(salt_pairs + max(0, -system_charge))
        n_cl = int(salt_pairs + max(0, system_charge))
        ion_atoms = place_ions(
            atoms=packed_atoms,
            box_lengths=box_lengths,
            n_na=n_na,
            n_cl=n_cl,
            cutoff=float(args.ion_cutoff),
            rng=np.random.default_rng(int(args.seed)),
            reject=reject_inside_slab(*env_geometry.pop("membrane_z_range")),
        )
    else:
        salt_pairs = 0
        n_na = 0
        n_cl = 0
        ion_atoms = []
        env_geometry.pop("membrane_z_range")

    protein_atoms, bb_entries = build_backbone_with_virtual_bb(
        protein_backbone_atoms, bb_type_by_residue, bb_secondary_by_residue
    )
    all_atoms = protein_atoms + bilayer_kept + ion_atoms
    write_pdb(runtime_pdb, all_atoms, box_lengths)
    protein_xyz_final = coords(protein_backbone_atoms)
    bilayer_xyz_final = coords(bilayer_kept)

    mapping_summary = {}
    if args.hybrid_mapping_output:
        mapping_h5 = Path(args.hybrid_mapping_output).expanduser().resolve()
        mapping_h5.parent.mkdir(parents=True, exist_ok=True)
        write_hybrid_mapping_h5(
            mapping_h5,
            bb_entries=bb_entries,
            total_martini_atoms=len(all_atoms),
            env_atom_indices=list(range(len(protein_atoms), len(all_atoms))),
            n_protein_atoms=len(protein_atoms),
        )
        mapping_summary["mapping_h5"] = str(mapping_h5)
        mapping_summary["bb_secondary_counts"] = dict(sorted(Counter(bb_secondary_by_residue.values()).items()))
        mapping_summary["bb_type_counts"] = dict(sorted(Counter(bb_type_by_residue.values()).items()))
        mapping_summary["bb_type_source_counts"] = dict(sorted(Counter(bb_type_source_by_residue.values()).items()))
        mapping_summary["bb_map_entries"] = int(len(bb_entries))
        if args.hybrid_bb_map_json_output:
            mapping_json = Path(args.hybrid_bb_map_json_output).expanduser().resolve()
            write_summary(mapping_json, {"bb_entries": bb_entries, "count": len(bb_entries)})
            mapping_summary["mapping_json"] = str(mapping_json)

    summary = {
        "input_protein_aa_pdb": str(protein_aa_pdb),
        "input_membrane_pdb": str(membrane_pdb),
        "input_protein_cg_pdb": str(protein_cg_pdb),
        "runtime_pdb": str(runtime_pdb),
        "protein_orientation_mode": str(args.protein_orientation_mode),
        "opm_reference": str(args.opm_reference) if args.opm_reference else None,
        "xy_scale": float(args.xy_scale),
        "box_angstrom": [float(v) for v in box_lengths],
        "protein_charge": int(protein_charge),
        "membrane_charge": int(membrane_charge),
        "explicit_ions": int(args.explicit_ions),
        "salt_pairs_target": int(salt_pairs),
        "na_added": int(n_na),
        "cl_added": int(n_cl),
        "net_charge": int(system_charge + n_na - n_cl),
        "protein_atoms": int(len(protein_atoms)),
        "bilayer_atoms_kept": int(len(bilayer_kept)),
        "ion_atoms_added": int(len(ion_atoms)),
        "total_atoms": int(len(all_atoms)),
        "protein_final_z_min": float(protein_xyz_final[:, 2].min()),
        "protein_final_z_max": float(protein_xyz_final[:, 2].max()),
        "lipid_final_z_min": float(bilayer_xyz_final[:, 2].min()),
        "lipid_final_z_max": float(bilayer_xyz_final[:, 2].max()),
    }
    summary.update(env_geometry)
    summary.update({f"belt_{k}": v for k, v in coverage.items()})
    summary.update(mapping_summary)
    return summary


def write_summary(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")


def run_prepare_command(argv):
    args = parse_prepare_args(argv)
    set_active_lipid_name(args.lipid_name)
    runtime_pdb = Path(args.runtime_pdb_output).expanduser().resolve()
    runtime_pdb.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "pdb_id": args.pdb_id,
        "prepare_structure": bool(args.prepare_structure),
        "stage": args.stage,
        "run_dir": args.run_dir,
    }

    if args.prepare_structure:
        summary.update(prepare_mixed_structure(args, runtime_pdb))
    elif not runtime_pdb.exists():
        raise FileNotFoundError(
            f"Runtime PDB not found for stage conversion: {runtime_pdb}. Run with --prepare-structure 1 first."
        )

    if args.stage:
        with temporary_env(
            {
                "UPSIDE_MARTINI_FF_DIR": os.environ.get(
                    "UPSIDE_MARTINI_FF_DIR", str(REPO_ROOT / "example" / "16.MARTINI" / "dryMARTINI_itp")
                ),
                "UPSIDE_MARTINI_ENERGY_CONVERSION": os.environ.get(
                    "UPSIDE_MARTINI_ENERGY_CONVERSION", str(DEFAULT_MARTINI_ENERGY_CONVERSION)
                ),
                "UPSIDE_MARTINI_LENGTH_CONVERSION": os.environ.get(
                    "UPSIDE_MARTINI_LENGTH_CONVERSION", str(DEFAULT_MARTINI_LENGTH_CONVERSION)
                ),
            }
        ):
            convert_stage(pdb_id=args.pdb_id, stage=args.stage, run_dir=args.run_dir, runtime_pdb=runtime_pdb)
        summary["upside_input"] = str(Path(args.run_dir).expanduser().resolve() / "test.input.up")

    summary_path = (Path(args.summary_json).expanduser().resolve() if args.summary_json
                    else runtime_pdb.with_suffix(".prep_summary.json"))
    write_summary(summary_path, summary)
    print(f"Preparation summary written to: {summary_path}")


SC_LIBRARY_REQUIRED_DATASETS = [
    "target_order",
    "grid_ang",
    "theta_grid",
    "rotamer_count",
    "rotamer_probability_fixed",
    "rotamer_full_energy_eup",
]


@contextlib.contextmanager
def temporary_env(updates):
    previous = {key: os.environ.get(key) for key in updates}
    try:
        for key, value in updates.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = str(value)
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


class Config:
    @staticmethod
    def text(name, default):
        return os.environ.get(name, default)

    @staticmethod
    def integer(name, default):
        return int(os.environ.get(name, str(default)))

    @staticmethod
    def floating(name, default):
        return float(os.environ.get(name, str(default)))


def workflow_path(value, base=REPO_ROOT / "example" / "16.MARTINI"):
    path = Path(value).expanduser()
    return path if path.is_absolute() else base / path


def generate_random_seed():
    seed = random.SystemRandom().randint(1, 2**32 - 1)
    return seed if seed != 0 else 1


def stage_seed(seed: int, stage: str) -> int:
    """Seed of one workflow stage: a splitmix64 hash of (seed, crc32(stage)), the hash main.cpp's
    system_seed applies to (seed, system). Every stage draws its own noise, so no stage replays another's.
    Production and its continuations share the stage name "production", so with the round carried a
    segmented production draws exactly the random numbers of the continuous run."""
    mask = (1 << 64) - 1
    z = (((int(seed) & 0xffffffff) << 32) | zlib.crc32(stage.encode())) + 0x9e3779b97f4a7c15
    z &= mask
    z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & mask
    z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & mask
    return int((z ^ (z >> 31)) >> 32)


def run_checked(cmd, cwd=None, log_file=None, echo=True):
    if echo:
        print(" ".join(str(x) for x in cmd))
    start = time.perf_counter()
    if log_file is None:
        proc = subprocess.Popen(
            [str(x) for x in cmd],
            cwd=str(cwd) if cwd else None,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="")
        rc = proc.wait()
        if rc != 0:
            raise subprocess.CalledProcessError(rc, [str(x) for x in cmd])
        elapsed = time.perf_counter() - start
        if echo:
            print(f"[workflow timing] command finished in {elapsed:.3f} s")
        return

    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("w", encoding="utf-8") as log:
        proc = subprocess.Popen(
            [str(x) for x in cmd],
            cwd=str(cwd) if cwd else None,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            if echo:
                print(line, end="")
            log.write(line)
        rc = proc.wait()
    elapsed = time.perf_counter() - start
    if echo:
        print(f"[workflow timing] command finished in {elapsed:.3f} s")
    if rc != 0:
        raise subprocess.CalledProcessError(rc, [str(x) for x in cmd])


def pdb_min_protein_lipid_distance(pdb_file: Path):
    atoms, _ = parse_pdb(pdb_file)
    protein = coords([a for a in atoms if a["resname"].upper() in PROTEIN_RESIDUE_NAMES])
    lipid = coords([a for a in atoms if lipid_resname(a["resname"])])
    min_d2 = float("inf")
    for i in range(0, lipid.shape[0], 2000):
        d = lipid[i:i + 2000][:, None, :] - protein[None, :, :]
        min_d2 = min(min_d2, float(np.einsum("ijk,ijk->ij", d, d).min()))
    return min_d2 ** 0.5


def prepare_workflow_hybrid_artifacts(args):
    print("=== Stage 0: Hybrid Packing + Mapping Export ===")
    cg_pdb = martinize_protein_cg(
        args.venv_python,
        args.martinize_script,
        workflow_path(args.protein_aa_pdb).resolve(),
        args.hybrid_prep_dir,
    )
    # membrane_top and opm_reference are optional; forwarding str(None) makes the child open a file
    # literally named "None".
    optional_prep_args = []
    if args.membrane_top is not None:
        optional_prep_args += ["--membrane-top", str(args.membrane_top)]
    if args.opm_reference is not None:
        optional_prep_args += ["--opm-reference", str(args.opm_reference)]
    run_prepare_command(
        [
            "--pdb-id", args.runtime_pdb_id,
            "--runtime-pdb-output", str(args.hybrid_packed_pdb),
            "--run-dir", str(args.run_dir),
            "--prepare-structure", "1",
            "--protein-aa-pdb", str(workflow_path(args.protein_aa_pdb).resolve()),
            "--protein-cg-pdb", str(cg_pdb),
            "--hybrid-mapping-output", str(args.hybrid_mapping_file),
            "--hybrid-bb-map-json-output", str(args.hybrid_prep_dir / "hybrid_bb_map.json"),
            "--lipid-name", args.lipid_name,
            "--membrane-pdb", str(args.membrane_pdb),
            *optional_prep_args,
            "--protein-orientation-mode", args.protein_orientation_mode,
            "--protein-pbc-margin", str(args.protein_pbc_margin),
            "--salt-molar", str(args.salt_molar),
            "--explicit-ions", str(args.explicit_ions),
            "--protein-lipid-cutoff", f"{float(args.protein_lipid_cutoff):.6g}",
            "--ion-cutoff", str(args.ion_cutoff),
            "--xy-scale", str(args.xy_scale),
            "--box-padding-xy", str(args.box_padding_xy),
            "--box-padding-z", str(args.box_padding_z),
            "--membrane-thickness-angstrom", str(args.membrane_thickness_angstrom),
            "--seed", str(args.prep_seed),
            "--summary-json", str(args.hybrid_prep_dir / "hybrid_prep_summary.json"),
        ]
    )
    # Overlap removal against the CG envelope is the packing guarantee; the backbone-lipid distance
    # is a sanity log (it can sit near the clearance since a BB bead may be the closest CG bead).
    min_gap = pdb_min_protein_lipid_distance(args.hybrid_packed_pdb)
    print(
        f"Hybrid packing: min protein(backbone)-lipid distance {min_gap:.3f} A "
        f"(CG-envelope clearance {float(args.protein_lipid_cutoff):.3f} A)"
    )
    if min_gap < 2.0:
        raise RuntimeError(f"Hybrid packing left a hard protein-lipid clash ({min_gap:.3f} A).")
    validate_hybrid_mapping(args.hybrid_mapping_file)
    copy_if_different(args.hybrid_packed_pdb, args.runtime_pdb_file)
    print(f"Runtime MARTINI PDB: {args.runtime_pdb_file}")


def h5_as_text(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)


def inject_hybrid_mapping(up_file: Path, mapping_file: Path):
    import h5py

    with h5py.File(mapping_file, "r") as src, h5py.File(up_file, "r+") as dst:
        src_inp = src["/input"]
        dst_inp = dst["/input"]
        n_atom = int(dst_inp["pos"].shape[0])
        if src_inp["hybrid_env_topology"]["protein_membership"].shape[0] != n_atom:
            raise ValueError(
                f"Hybrid mapping n_atom mismatch for {up_file}: up has {n_atom}, mapping has "
                f"{src_inp['hybrid_env_topology']['protein_membership'].shape[0]}"
            )
        groups = ["hybrid_control", "hybrid_bb_map", "hybrid_env_topology"]
        groups += [g for g in ("chain_break",) if g in src_inp]
        for group in groups:
            src.copy(src_inp[group], dst_inp, name=group)


def set_stage_label(up_file: Path, stage_label: str):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        grp = h5.require_group("input").require_group("stage_parameters")
        grp.attrs["current_stage"] = np.bytes_(stage_label)


def set_hybrid_control_mode(up_file: Path, activation_stage: str, preprod_mode="rigid_body"):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        grp = h5.require_group("input").require_group("hybrid_control")
        grp.attrs["activation_stage"] = np.bytes_(activation_stage)
        grp.attrs["preprod_protein_mode"] = np.bytes_(preprod_mode)


def inject_protein_position_restraints(up_file: Path, spring_const: float = 10.0):
    """Hold the protein's N/CA/C carriers during the soft-core minimisation. The O and BB slots are derived
    sites the engine places, so they take no restraint."""
    import h5py

    with h5py.File(up_file, "r+") as h5:
        pm = h5["/input/hybrid_env_topology/protein_membership"][:]
        names = np.char.strip(h5["/input/atom_names"][:].astype(str))
        carriers = np.where((pm >= 0) & np.isin(names, ["N", "CA", "C"]))[0].astype(np.int32)
        grp = h5.create_group("/input/potential/restraint_position_protein")
        grp.attrs["arguments"] = np.array([np.bytes_("pos")])
        grp.create_dataset("restraint_indices", data=carriers)
        grp.create_dataset("ref_pos", data=np.asarray(h5["/input/pos"][:][carriers, :, 0], dtype=np.float32))
        grp.create_dataset("spring_const_xyz", data=np.full((len(carriers), 3), float(spring_const), dtype=np.float32))


def remove_protein_position_restraints(up_file: Path):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        del h5["/input/potential/restraint_position_protein"]


def set_hybrid_production_controls(up_file: Path, args):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        grp = h5.require_group("input").require_group("hybrid_control")
        grp.attrs["sc_env_po4_z_hold_steps"] = np.int32(args.sc_env_po4_z_hold_steps)


def ensure_sc_martini_library(args):
    import h5py

    library = args.martini_h5
    if not library.exists():
        raise FileNotFoundError(library)
    with h5py.File(library, "r") as h5:
        sc_grp = h5["sc_table"] if "sc_table" in h5 else h5
        missing = [name for name in SC_LIBRARY_REQUIRED_DATASETS if name not in sc_grp]
    if missing:
        raise ValueError(f"{library} missing required SC datasets: {','.join(missing)}")


def assert_hybrid_stage_active(
    up_file: Path,
    expected_stage: str,
    expected_activation: str,
    require_interface_nodes: bool = False,
):
    import h5py

    with h5py.File(up_file, "r") as h5:
        inp = h5["/input"]
        pot = inp["potential"]
        stage = h5_as_text(inp["stage_parameters"].attrs.get("current_stage", b"")).strip()
        hy = inp["hybrid_control"].attrs
        activation = h5_as_text(hy.get("activation_stage", b"")).strip()
        if stage != expected_stage or activation != expected_activation:
            raise ValueError(
                f"{up_file}: expected stage={expected_stage}, activation={expected_activation}; "
                f"got stage={stage}, activation={activation}"
            )
        if "particle_class" in inp and "mass" in inp:
            particle_class = np.asarray(
                [h5_as_text(value).strip().upper() for value in inp["particle_class"][:]],
                dtype=object,
            )
            protein_mass = np.asarray(inp["mass"][:], dtype=np.float64)[
                particle_class == "PROTEIN"
            ]
            if protein_mass.size and not np.allclose(
                protein_mass, 1.0, rtol=0.0, atol=1.0e-6
            ):
                raise ValueError(
                    f"{up_file}: protein carriers must retain Upside unit mass; "
                    f"observed range {protein_mass.min():.6g}-{protein_mass.max():.6g}"
                )
        if require_interface_nodes:
            missing = [
                node
                for node in ("martini_potential", "martini_sc_table_1body")
                if node not in pot
            ]
            if missing:
                raise ValueError(
                    f"{up_file}: missing required hybrid interface node(s): {', '.join(missing)}"
                )
        env_membership = inp["hybrid_env_topology"]["protein_membership"][:]
        if not np.any(env_membership < 0):
            raise ValueError(f"{up_file}: no non-protein environment targets found")
        if "hybrid_bb_map" in inp:
            bb_map = inp["hybrid_bb_map"]
            bb_idx = np.asarray(bb_map["bb_atom_index"][:], dtype=np.int32)
            bb_types = np.asarray([h5_as_text(x).strip() for x in bb_map["bb_type"][:]], dtype=object)
            charges = np.asarray(inp["charges"][:], dtype=float) if "charges" in inp else None
            if "chain_break" in inp:
                chain_first = np.asarray(inp["chain_break"]["chain_first_residue"][:], dtype=np.int32)
                first_res = np.concatenate((np.asarray([0], dtype=np.int32), chain_first))
            else:
                first_res = np.asarray([0], dtype=np.int32)
            next_res = np.concatenate((first_res[1:], np.asarray([len(bb_idx)], dtype=np.int32)))
            for start, stop in zip(first_res, next_res):
                if stop <= start:
                    continue
                end = int(stop - 1)
                start = int(start)
                if bb_types[start] != "Qd" or bb_types[end] != "Qa":
                    raise ValueError(
                        f"{up_file}: strand residues {start}-{end} must have Qd/Qa BB endpoints; "
                        f"got {bb_types[start]!r}/{bb_types[end]!r}"
                    )
                if charges is not None:
                    q_start = float(charges[int(bb_idx[start])])
                    q_end = float(charges[int(bb_idx[end])])
                    if abs(q_start - 1.0) > 1.0e-6 or abs(q_end + 1.0) > 1.0e-6:
                        raise ValueError(
                            f"{up_file}: strand residues {start}-{end} BB endpoint charges must be +1/-1; "
                            f"got {q_start:g}/{q_end:g}"
                        )
    print(f"Hybrid activation verified: stage={expected_stage}, activation_stage={expected_activation}")


def stage_conversion_env(args, npt_enable: int, lipidhead_fc: float, dynamics_phase: str):
    return {
        "UPSIDE_HOME": str(args.upside_home),
        "UPSIDE_MARTINI_FF_DIR": str(args.itp_dir),
        "UPSIDE_MARTINI_ENERGY_CONVERSION": str(args.martini_energy_conversion),
        "UPSIDE_MARTINI_LENGTH_CONVERSION": str(args.martini_length_conversion),
        "UPSIDE_MARTINI_DYNAMICS_PHASE": dynamics_phase,
        "UPSIDE_NPT_ENABLE": str(int(npt_enable)),
        # Tensionless (zero lateral pressure) dry-MARTINI membrane barostat.
        "UPSIDE_NPT_TARGET_PXY": "0.0",
        "UPSIDE_NPT_TARGET_PZ": "0.0",
        "UPSIDE_NPT_INTERVAL": str(args.npt_interval),
        # Monte-Carlo barostat steps: couple xy, freeze z (mc_dmax_z=0).
        "UPSIDE_NPT_MC_DMAX_XY": str(args.npt_mc_dmax_xy),
        "UPSIDE_NPT_MC_DMAX_Z": "0.0",
        "UPSIDE_BILAYER_LIPIDHEAD_FC": str(lipidhead_fc),
    }


def inject_hybrid_interface_nodes(args, target_file: Path, current_stage: str, activation_stage: str):
    ensure_sc_martini_library(args)
    set_hybrid_control_mode(target_file, activation_stage)
    set_hybrid_production_controls(target_file, args)
    inject_stage7_sc_table_nodes(
        up_file=target_file,
        martini_h5=args.martini_h5,
        upside_home=args.upside_home,
        rama_library=args.upside_rama_library,
        rama_sheet_mixing=args.upside_rama_sheet_mixing,
        hbond_energy=args.upside_hbond_energy,
        reference_state_rama=args.upside_reference_state_rama,
    )
    assert_hybrid_stage_active(
        target_file,
        current_stage,
        activation_stage,
        require_interface_nodes=True,
    )


def prepare_stage_file(args, target_file: Path, tables: str, npt_enable: int, lipidhead_fc: float,
                       dynamics_phase: str, stage_label: str):
    start = time.perf_counter()
    with temporary_env(stage_conversion_env(args, npt_enable, lipidhead_fc, dynamics_phase)):
        run_prepare_command(
            [
                "--lipid-name", args.lipid_name,
                "--pdb-id", args.runtime_pdb_id,
                "--runtime-pdb-output", str(args.runtime_pdb_file),
                "--prepare-structure", "0",
                "--stage", tables,
                "--run-dir", str(args.run_dir),
                "--summary-json", str(args.hybrid_prep_dir / f"stage_{target_file.stem}.summary.json"),
            ]
        )
    target_file.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(args.run_dir / "test.input.up"), str(target_file))
    inject_hybrid_mapping(target_file, args.hybrid_mapping_file)
    set_stage_label(target_file, stage_label)
    inject_particles_table(up_file=target_file, martini_h5=args.martini_h5,
                           energy_conversion=args.martini_energy_conversion,
                           length_conversion=args.martini_length_conversion)
    activation = "production" if stage_label in {"production", "production_handoff"} else "minimization"
    inject_hybrid_interface_nodes(args, target_file, stage_label, activation)
    print(
        f"[workflow timing] stage {stage_label} preparation/injection finished "
        f"in {time.perf_counter() - start:.3f} s -> {target_file}"
    )


def assert_no_core_overlap(up_file: Path, temperature: float):
    """Refuse to start MD from a configuration with a pair inside the repulsive core.

    The limit is read off the potential. A pair at energy U carries the Boltzmann weight exp(-U/kT), so
    with N pairs inside the cutoff, a pair above U* = kT ln(CORE_OVERLAP_ODDS * N) would appear in fewer
    than one equilibrium configuration in CORE_OVERLAP_ODDS. Each pair is evaluated on the stage file's own
    tables, at the positions MD starts from, for every pair the MARTINI term computes (intra-protein pairs
    are the Upside force field's and are skipped, as the engine skips them).
    """
    import h5py

    kT = float(temperature)
    with h5py.File(up_file, "r") as h5:
        pos = np.asarray(h5["/input/pos"][:, :, 0], dtype=np.float64)
        mp = h5["/input/potential/martini_potential"]
        box = np.array([mp.attrs["x_len"], mp.attrs["y_len"], mp.attrs["z_len"]], dtype=np.float64)
        pairs = mp["pairs"][:]
        coeff_index = mp["coefficient_indices"][:]
        grids = mp["combined_energy_grids"][:]
        r_min, r_max = float(mp.attrs["r_min_ang"]), float(mp.attrs["r_max_ang"])
        membership = h5["/input/hybrid_env_topology/protein_membership"][:]
        names = np.char.strip(h5["/input/atom_names"][:].astype(str))
        types = np.char.strip(h5["/input/type"][:].astype(str))

    keep = ~((membership[pairs[:, 0]] >= 0) & (membership[pairs[:, 0]] == membership[pairs[:, 1]]))
    pairs, coeff_index = pairs[keep], coeff_index[keep]
    n_knot = grids.shape[1]
    energy, dist, atom_i, atom_j = [], [], [], []
    for start in range(0, len(pairs), 2_000_000):
        p = pairs[start:start + 2_000_000]
        delta = pos[p[:, 0]] - pos[p[:, 1]]
        delta -= box * np.rint(delta / box)
        d = np.sqrt(np.sum(delta * delta, axis=1))
        inside = d < r_max
        p, d, rows = p[inside], d[inside], coeff_index[start:start + 2_000_000][inside]
        coord = np.clip((d - r_min) / (r_max - r_min) * (n_knot - 1), 0.0, n_knot - 1 - 1e-9)
        lo = coord.astype(np.int64)
        frac = coord - lo
        energy.append(grids[rows, lo] * (1.0 - frac) + grids[rows, lo + 1] * frac)
        dist.append(d)
        atom_i.append(p[:, 0])
        atom_j.append(p[:, 1])
    energy, dist, atom_i, atom_j = (np.concatenate(x) for x in (energy, dist, atom_i, atom_j))
    limit = kT * np.log(CORE_OVERLAP_ODDS * max(energy.size, 1))
    print(f"Core check {up_file.name}: {energy.size} MARTINI pairs inside {r_max:g} A; largest pair energy "
          f"{energy.max() / kT:.2f} kT (limit {limit / kT:.2f} kT)")
    offenders = np.where(energy > limit)[0]
    if offenders.size:
        worst = offenders[np.argsort(-energy[offenders])][:5]
        detail = "; ".join(f"{atom_i[k]}:{names[atom_i[k]]}({types[atom_i[k]]})-{atom_j[k]}:{names[atom_j[k]]}"
                           f"({types[atom_j[k]]}) at {dist[k]:.3f} A, {energy[k] / kT:.3g} kT" for k in worst)
        raise RuntimeError(
            f"{up_file}: MD would start with {offenders.size} pairs inside the repulsive core, above "
            f"{limit / kT:.2f} kT (an equilibrium configuration shows such a pair less than once in "
            f"{CORE_OVERLAP_ODDS:.0e}); the worst: {detail}. The minimisations that precede this stage did "
            "not remove them."
        )


def refresh_particle_mobility_friction(up_file: Path):
    """Recount each protein carrier's lipid contacts at the stage's positions and set its friction to the sum
    of their frictions, each lipid bead counted at its own molecule's bead friction. Bead frictions themselves
    do not depend on position and are left as built."""
    import h5py

    with h5py.File(up_file, "r+") as h5:
        if "/input/brownian" not in h5:
            return
        brownian = h5["/input/brownian"]
        if h5_as_text(brownian.attrs.get("transport_observable", b"")) != "bare_martini_particle_lateral_diffusion":
            return

        atom_index = brownian["atom_index"][:]
        particle_class = h5["/input/particle_class"][:]
        atom_names = h5["/input/atom_names"][:]
        residue_names = np.char.strip(h5["/input/residue_names"][:].astype(str))
        position = np.asarray(h5["/input/pos"][:], dtype=np.float64).reshape((-1, 3))
        protein_carrier = (
            (particle_class == b"PROTEIN")
            & np.isin(atom_names, np.array([b"N", b"CA", b"C"], dtype=atom_names.dtype))
        )
        potential = h5["/input/potential/martini_potential"]
        box = np.array(
            [potential.attrs["x_len"], potential.attrs["y_len"], potential.attrs["z_len"]],
            dtype=np.float64,
        )
        cutoff = float(brownian.attrs["interface_friction_cutoff_A"])
        carrier_index = np.where(protein_carrier)[0]
        contact_count = np.zeros(particle_class.size, dtype=np.int32)
        carrier_friction = np.zeros(particle_class.size, dtype=np.float64)
        for name, bead_friction in zip(brownian.attrs["lipid_names"], brownian.attrs["lipid_bead_friction_up"]):
            beads = np.where((particle_class == b"LIPID") & (residue_names == h5_as_text(name)))[0]
            count = np.zeros(particle_class.size, dtype=np.int32)
            for start in range(0, beads.size, 512):
                delta = position[carrier_index, None, :] - position[beads[start:start + 512]][None, :, :]
                delta -= box * np.rint(delta / box)
                count[carrier_index] += np.sum(np.sum(delta * delta, axis=2) < cutoff * cutoff, axis=1)
            contact_count += count
            carrier_friction += float(bead_friction) * count

        friction = brownian["friction"][:]
        protein_index = protein_carrier[atom_index]
        friction[protein_index] = carrier_friction[atom_index[protein_index]]
        brownian["friction"][:] = friction

        contact_atom_index = np.where(contact_count > 0)[0].astype(np.int32)
        for name in ("protein_contact_atom_index", "protein_lipid_contact_count"):
            if name in brownian:
                del brownian[name]
        brownian.create_dataset("protein_contact_atom_index", data=contact_atom_index)
        brownian.create_dataset(
            "protein_lipid_contact_count",
            data=contact_count[contact_atom_index],
        )


def handoff_initial_position(input_file: Path, output_file: Path):
    set_initial_position(input_file, output_file)
    refresh_particle_mobility_friction(output_file)


def input_momentum_restart_valid(up_file: Path):
    import h5py

    with h5py.File(up_file, "r") as h5:
        if "/input/mom" not in h5:
            return False
        return int(h5["/input/mom"].attrs.get("restart_valid", 0)) != 0


def output_restart_state_valid(up_file: Path):
    import h5py

    with h5py.File(up_file, "r") as h5:
        if "/output/mom" not in h5 or h5["/output/mom"].shape[0] == 0:
            return False
        return int(h5["/output/mom"].attrs.get("restart_final_state_valid", 0)) != 0


def mark_output_restart_state(up_file: Path, nsteps: int):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        if "/output/mom" not in h5 or h5["/output/mom"].shape[0] == 0:
            raise RuntimeError(f"{up_file} has no output/mom; cannot mark restart state")
        # a continuation runs on from the round its input momentum was recorded at
        start_round = 0
        if int(h5["/input/mom"].attrs.get("restart_valid", 0)) != 0:
            start_round = int(h5["/input/mom"].attrs["round"])
        final_round = int(h5["/output/round"][-1])
        if final_round != start_round + int(nsteps):
            raise RuntimeError(
                f"{up_file} final output round {final_round} does not match expected {start_round + int(nsteps)}"
            )
        h5["/output/mom"].attrs["restart_final_state_valid"] = np.int8(1)


def run_minimization_stage(args, stage_label: str, up_file: Path, max_iter: int):
    print(f"=== Stage {stage_label}: Minimization ===")
    cmd = [
        args.upside_executable,
        up_file,
        "--duration", "0",
        "--frame-interval", "1",
        "--temperature", args.temperature,
        "--thermostat-timescale", args.thermostat_timescale,
        "--thermostat-interval", args.thermostat_interval,
        "--seed", stage_seed(args.seed, stage_label),
        "--disable-recentering",
        "--minimize",
        "--min-max-iter", max_iter,
        "--min-energy-tol", "1e-6",
        "--min-force-tol", "1e-3",
        "--min-step", "0.01",
    ]
    run_checked(cmd, log_file=args.log_dir / f"stage_{stage_label}.min.log")
    promote_minimized_state_to_input(up_file)


def promote_minimized_state_to_input(up_file: Path):
    """The minimised positions become the stage's input. A minimisation does not touch the momentum, so the
    input momentum and its round stay as they were."""
    import h5py

    with h5py.File(up_file, "r+") as h5:
        if "/output/pos" not in h5 or h5["/output/pos"].shape[0] == 0:
            raise RuntimeError(f"{up_file} has no minimized output/pos to promote")
        last_pos = h5["/output/pos"][-1, 0, :, :][:, :, np.newaxis]
        del h5["/input/pos"]
        h5.create_dataset("/input/pos", data=last_pos)
        del h5["/output"]
    print(f"Promoted minimized state to input of {up_file}")
    refresh_particle_mobility_friction(up_file)


def run_md_stage(
    args,
    stage_label: str,
    up_file: Path,
    nsteps: int,
    dt: float,
    frame_steps: int,
    seed_stage: str,
    echo: bool = True,
):
    effective_frame_steps = int(frame_steps)
    if effective_frame_steps >= int(nsteps):
        effective_frame_steps = max(1, int(nsteps) // 10)
        print(f"NOTICE: frame_steps ({frame_steps}) >= nsteps ({nsteps}); using frame_steps={effective_frame_steps}")
    frame_interval = f"{effective_frame_steps * float(dt):.10g}"
    assert_no_core_overlap(up_file, args.temperature)
    print(f"=== Stage {stage_label}: MD ===")
    cmd = [
        args.upside_executable,
        up_file,
        "--duration-steps", nsteps,
        "--frame-interval", frame_interval,
        "--temperature", args.temperature,
        "--time-step", dt,
        "--integrator", MARTINI_MD_INTEGRATOR,
        "--thermostat-timescale", args.thermostat_timescale,
        "--thermostat-interval", args.thermostat_interval,
        "--seed", stage_seed(args.seed, seed_stage),
        "--disable-recentering",
        "--record-momentum",
    ]
    if input_momentum_restart_valid(up_file):
        cmd.append("--restart-using-momentum")
    run_checked(cmd, log_file=args.log_dir / f"stage_{stage_label}.log", echo=echo)
    mark_output_restart_state(up_file, nsteps)


def promote_md_output_state_to_input(up_file: Path):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        if "/output/mom" not in h5 or h5["/output/mom"].shape[0] == 0:
            raise RuntimeError(f"{up_file} has no MD output/mom to promote")
    set_initial_position(up_file, up_file)
    with h5py.File(up_file, "r+") as h5:
        del h5["/output"]
    refresh_particle_mobility_friction(up_file)


def run_stage70_burnin(args, stage_file: Path):
    """The rigid burn-in, then the solvation gate on its last quarter, then its final state becomes production's
    input."""
    nsteps = int(args.prod_70_burnin_nsteps)
    if nsteps <= 0:
        print("=== Stage 7.0: Rigid interface-handoff burn-in skipped (nsteps <= 0) ===")
        assert_environment_solvation(args, stage_file)
        return
    print(f"=== Stage 7.0: Rigid interface-handoff burn-in ({nsteps} steps) ===")
    run_md_stage(args, "7.0.burnin", stage_file, nsteps, args.prod_time_step, args.prod_frame_steps,
                 "7.0.burnin", echo=False)
    assert_environment_solvation(args, stage_file)
    promote_md_output_state_to_input(stage_file)
    print(f"Promoted rigid handoff final state for production: {stage_file}")


def extract_stage_vtf(args, stage_label: str, stage_file: Path, mode: str):
    vtf_file = args.run_dir / f"{args.pdb_id}.stage_{stage_label}.vtf"
    print(f"=== Stage {stage_label}: VTF Extraction (mode {mode}) ===")
    cmd = [
        sys.executable,
        args.extract_vtf_script,
        stage_file,
        vtf_file,
        stage_file,
        args.runtime_pdb_id,
        "--mode",
        mode,
        "--split-segments",
    ]
    run_checked(cmd)


def infer_stage70_label_from_file(stage_file: Path, pdb_id: str):
    match = re.fullmatch(rf"{re.escape(pdb_id)}\.stage_(7\.\d+)\.up", stage_file.name)
    return match.group(1) if match else ""


def infer_next_stage70_label_from_source(source_file: Path, pdb_id: str):
    pattern = re.compile(rf"^{re.escape(pdb_id)}\.stage_7\.(\d+)\.up$")
    match = pattern.fullmatch(source_file.name)
    if not match:
        return "7.1"
    return f"7.{int(match.group(1)) + 1}"


def resolve_continuation_outputs(args):
    if args.continue_stage_70_from is None:
        return None, None, None
    source = args.continue_stage_70_from.resolve()
    label = args.continue_stage_70_label
    if not label and args.continue_stage_70_output:
        label = infer_stage70_label_from_file(args.continue_stage_70_output, args.pdb_id)
        if not label:
            raise ValueError(f"CONTINUE_STAGE_70_OUTPUT must be named {args.pdb_id}.stage_7.N.up")
    if not label:
        label = infer_next_stage70_label_from_source(source, args.pdb_id)
    if not re.fullmatch(r"7\.\d+", label):
        raise ValueError(f"Continuation stage label must be numeric stage_7.N, got {label}")
    output = args.continue_stage_70_output or (args.checkpoint_dir / f"{args.pdb_id}.stage_{label}.up")
    return source, output, label


def run_stage70_continuation(args, source_file: Path, output_file: Path, stage_label: str):
    """Continue production: the source's /input unchanged, with the positions, momentum, round and PO4
    hold clock of its last frame, so the continuation is the continuous run."""
    if not source_file.exists():
        raise FileNotFoundError(source_file)
    assert_hybrid_stage_active(source_file, "production", "production")
    if not output_restart_state_valid(source_file):
        raise ValueError(
            "Production continuation requires a source stage generated by the current workflow with "
            "a validated final restart state. Regenerate the previous production stage with the current workflow."
        )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if source_file.resolve() != output_file.resolve():
        shutil.copy2(source_file, output_file)
    set_initial_position(source_file, output_file, continuation=True)
    run_md_stage(args, stage_label, output_file, args.prod_70_nsteps, args.prod_time_step, args.prod_frame_steps,
                 "production")
    extract_stage_vtf(args, stage_label, output_file, "2")


def normalize_hybrid_workflow_args(args):
    args.upside_home = workflow_path(args.upside_home, REPO_ROOT).resolve()
    args.run_dir = workflow_path(args.run_dir).resolve()
    args.checkpoint_dir = args.run_dir / "checkpoints"
    args.log_dir = args.run_dir / "logs"
    args.hybrid_prep_dir = args.run_dir / "hybrid_prep"
    args.runtime_pdb_file = args.hybrid_prep_dir / f"{args.runtime_pdb_id}.MARTINI.pdb"
    args.hybrid_mapping_file = args.hybrid_prep_dir / "hybrid_mapping.h5"
    args.hybrid_packed_pdb = args.hybrid_prep_dir / "hybrid_packed.MARTINI.pdb"
    args.upside_executable = args.upside_home / "obj" / "upside"
    args.venv_python = args.upside_home / ".venv" / "bin" / "python"
    args.martinize_script = PY_DIR / "martinize.py"
    args.itp_dir = args.upside_home / "example" / "16.MARTINI" / "dryMARTINI_itp"
    args.martini_h5 = args.upside_home / "parameters" / "ff_2.1" / "martini.h5"
    args.upside_rama_library = args.upside_home / "parameters" / "common" / "rama.dat"
    args.upside_rama_sheet_mixing = args.upside_home / "parameters" / "ff_2.1" / "sheet"
    args.upside_hbond_energy = args.upside_home / "parameters" / "ff_2.1" / "hbond.h5"
    args.upside_reference_state_rama = args.upside_home / "parameters" / "common" / "rama_reference.pkl"
    args.sc_env_po4_z_hold_steps = DEFAULT_SC_ENV_PO4_Z_HOLD_STEPS
    args.npt_interval = DEFAULT_NPT_INTERVAL
    args.npt_mc_dmax_xy = DEFAULT_NPT_MC_DMAX_XY
    args.martini_energy_conversion = DEFAULT_MARTINI_ENERGY_CONVERSION
    args.martini_length_conversion = DEFAULT_MARTINI_LENGTH_CONVERSION
    args.extract_vtf_script = workflow_path(args.extract_vtf_script).resolve()
    args.continue_stage_70_from = workflow_path(args.continue_stage_70_from).resolve() if args.continue_stage_70_from else None
    args.continue_stage_70_output = workflow_path(args.continue_stage_70_output).resolve() if args.continue_stage_70_output else None
    if args.prep_seed is None:
        args.prep_seed = generate_random_seed()
    if args.seed is None:
        args.seed = generate_random_seed()
    for path in [args.run_dir, args.checkpoint_dir, args.log_dir, args.hybrid_prep_dir]:
        path.mkdir(parents=True, exist_ok=True)
    return args


def add_hybrid_workflow_arguments(parser):
    parser.add_argument("--pdb-id", default=Config.text("PDB_ID", None))
    parser.add_argument("--runtime-pdb-id", default=Config.text("RUNTIME_PDB_ID", None))
    parser.add_argument("--upside-home", default=Config.text("UPSIDE_HOME", str(REPO_ROOT)))
    parser.add_argument("--run-dir", default=Config.text("RUN_DIR", None))
    parser.add_argument("--protein-aa-pdb", default=Config.text("PROTEIN_AA_PDB", None))
    parser.add_argument("--lipid-name", default=Config.text("LIPID_NAME", None),
                        help="Membrane lipid moleculetype name(s), comma-separated for a mixed "
                             "bilayer (e.g. DOPC or POPE,POPG)")
    parser.add_argument("--charmm-gui-dir", default=Config.text("CHARMM_GUI_DIR", None),
                        help="CHARMM-GUI Martini job dir; membrane read from gromacs/step5_charmm2gmx.pdb + system.top")
    parser.add_argument("--membrane-pdb", default=Config.text("MEMBRANE_PDB", None))
    parser.add_argument("--membrane-top", default=Config.text("MEMBRANE_TOP", None))
    parser.add_argument("--opm-reference", default=Config.text("OPM_REFERENCE", None),
                        help="OPM membrane-oriented reference PDB (for --protein-orientation-mode opm)")
    parser.add_argument("--extract-vtf-script", default=Config.text("EXTRACT_VTF_SCRIPT", str(PY_DIR / "martini_extract_vtf.py")))
    parser.add_argument("--salt-molar", type=float, default=Config.floating("SALT_MOLAR", 0.15))
    parser.add_argument("--explicit-ions", type=int, choices=[0, 1], default=Config.integer("EXPLICIT_IONS", 1))
    parser.add_argument("--protein-lipid-cutoff", type=float, default=Config.floating("PROTEIN_LIPID_CUTOFF", 0.0))
    parser.add_argument("--ion-cutoff", type=float, default=Config.floating("ION_CUTOFF", 10.0))
    parser.add_argument("--xy-scale", type=float, default=Config.floating("XY_SCALE", 1.0))
    parser.add_argument("--box-padding-xy", type=float, default=Config.floating("BOX_PADDING_XY", 0.0))
    parser.add_argument("--box-padding-z", type=float, default=Config.floating("BOX_PADDING_Z", 0.0))
    parser.add_argument("--membrane-thickness-angstrom", type=float, default=Config.floating("MEMBRANE_THICKNESS_ANGSTROM", 0.0),
                        help="Equilibrated dry-MARTINI membrane thickness for the ion count "
                             "(0 = measure from the compressed packed lipids)")
    parser.add_argument("--protein-orientation-mode", choices=["input", "lay-flat", "opm"], default=Config.text("PROTEIN_ORIENTATION_MODE", "opm"))
    parser.add_argument("--protein-pbc-margin", type=float, default=Config.floating("PROTEIN_PBC_MARGIN", 15.0))
    parser.add_argument("--temperature", type=float, default=Config.floating("TEMPERATURE", 0.8647))
    parser.add_argument("--thermostat-timescale", type=float, default=Config.floating("THERMOSTAT_TIMESCALE", 5.0))
    parser.add_argument("--thermostat-interval", type=int, default=Config.integer("THERMOSTAT_INTERVAL", -1))
    # CHARMM-GUI step6.0 and step6.1 each run 5000 steepest-descent steps.
    parser.add_argument("--min-60-max-iter", type=int, default=Config.integer("MIN_60_MAX_ITER", 5000))
    parser.add_argument("--min-61-max-iter", type=int, default=Config.integer("MIN_61_MAX_ITER", 5000))
    parser.add_argument("--min-70-max-iter", type=int, default=Config.integer("MIN_70_MAX_ITER", 500))
    parser.add_argument("--eq-62-nsteps", type=int, default=Config.integer("EQ_62_NSTEPS", 500))
    parser.add_argument("--eq-63-nsteps", type=int, default=Config.integer("EQ_63_NSTEPS", 500))
    parser.add_argument("--eq-64-nsteps", type=int, default=Config.integer("EQ_64_NSTEPS", 500))
    parser.add_argument("--eq-65-nsteps", type=int, default=Config.integer("EQ_65_NSTEPS", 500))
    parser.add_argument("--eq-66-nsteps", type=int, default=Config.integer("EQ_66_NSTEPS", 60000))
    parser.add_argument("--prod-70-burnin-nsteps", type=int, default=Config.integer("PROD_70_BURNIN_NSTEPS", 40000))
    parser.add_argument("--prod-70-nsteps", type=int, default=Config.integer("PROD_70_NSTEPS", 10000),
                        help="production steps of stage 7.0 and of each continuation; 0 ends stage 7.0 after the "
                             "burn-in, leaving its final state as the input of the stage file (a production seed)")
    parser.add_argument("--eq-time-step", type=float, default=Config.floating("EQ_TIME_STEP", MARTINI_MD_TIME_STEP))
    parser.add_argument("--prod-time-step", type=float, default=Config.floating("PROD_TIME_STEP", MARTINI_MD_TIME_STEP))
    parser.add_argument("--eq-frame-steps", type=int, default=Config.integer("EQ_FRAME_STEPS", 1000))
    parser.add_argument("--prod-frame-steps", type=int, default=Config.integer("PROD_FRAME_STEPS", 50))
    parser.add_argument("--prod-70-npt-enable", type=int, default=Config.integer("PROD_70_NPT_ENABLE", 0))
    parser.add_argument("--prep-seed", default=os.environ.get("PREP_SEED"))
    parser.add_argument("--seed", default=os.environ.get("SEED"))
    parser.add_argument("--continue-stage-70-from", default=Config.text("CONTINUE_STAGE_70_FROM", ""))
    parser.add_argument("--continue-stage-70-output", default=Config.text("CONTINUE_STAGE_70_OUTPUT", ""))
    parser.add_argument("--continue-stage-70-label", default=Config.text("CONTINUE_STAGE_70_LABEL", ""))


def parse_hybrid_workflow_args(argv):
    parser = argparse.ArgumentParser(description="Run the hybrid dry-MARTINI workflow.")
    add_hybrid_workflow_arguments(parser)
    args = parser.parse_args(argv)
    if not args.pdb_id or not str(args.pdb_id).strip():
        raise ValueError("A PDB id is required (--pdb-id or PDB_ID, e.g. 1rkl or glpG-RKRK-79HIS)")
    if args.run_dir is None:
        args.run_dir = f"outputs/martini_{args.pdb_id}_hybrid"
    if args.runtime_pdb_id is None:
        args.runtime_pdb_id = f"{args.pdb_id}_hybrid"
    if args.protein_aa_pdb is None:
        args.protein_aa_pdb = f"pdb/{args.pdb_id}.pdb"
    if not args.lipid_name or not str(args.lipid_name).strip():
        raise ValueError("A lipid name is required (--lipid-name, e.g. DOPC or POPE,POPG)")
    set_active_lipid_name(args.lipid_name)
    if args.charmm_gui_dir:
        gmx = Path(args.charmm_gui_dir).expanduser() / "gromacs"
        if args.membrane_pdb is None:
            args.membrane_pdb = str(gmx / "step5_charmm2gmx.pdb")
        if args.membrane_top is None:
            args.membrane_top = str(gmx / "system.top")
    if args.membrane_pdb is None:
        raise ValueError("A CHARMM-GUI membrane is required (--charmm-gui-dir or --membrane-pdb)")
    if args.protein_orientation_mode == "opm" and not args.opm_reference:
        raise ValueError("--protein-orientation-mode opm requires --opm-reference")
    args.prep_seed = int(args.prep_seed) if args.prep_seed not in (None, "") else None
    args.seed = int(args.seed) if args.seed not in (None, "") else None
    args = normalize_hybrid_workflow_args(args)
    if args.protein_lipid_cutoff <= 0.0:
        args.protein_lipid_cutoff = derive_lipid_contact_clearance_angstrom(args.upside_home, args.lipid_name)
    return args


def validate_hybrid_workflow_args(args):
    if not args.upside_executable.exists():
        raise FileNotFoundError(args.upside_executable)
    required = [
        workflow_path(args.protein_aa_pdb).resolve(),
        Path(args.membrane_pdb).expanduser().resolve(),
        args.venv_python,
        args.martinize_script,
        args.martini_h5,
        args.upside_rama_library,
        args.upside_hbond_energy,
        args.upside_reference_state_rama,
        args.extract_vtf_script,
    ]
    if args.opm_reference:
        required.append(Path(args.opm_reference).expanduser().resolve())
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)
    require_bilayer_lipids(args.upside_home, args.lipid_name)
    lipid_diffusion_targets(active_lipid_names())
    # The belt is the protein's published hydrophobic band, so it only exists if an OPM reference was
    # given. Substituting 0.0 makes every |z - midplane| <= 0 test empty, which turned a missing input into
    # a downstream "no backbone sites inside the belt" failure on a correctly inserted protein. None means
    # undefined, and the solvation gate reports that it could not run rather than failing or passing.
    args.belt_half_thickness_a = (
        read_opm_hydrophobic_half_thickness(Path(args.opm_reference).expanduser().resolve())
        if args.opm_reference
        else None
    )


def workflow_stage_files(args):
    files = {label: args.checkpoint_dir / f"{args.pdb_id}.stage_{label}.up" for label in ("6.0", "6.1", "7.0")}
    files.update({label: args.checkpoint_dir / f"{args.pdb_id}.stage_{label}.up" for label, _ in EQUILIBRATION_STAGES})
    return files


def print_hybrid_workflow_summary(args, source=None, output=None):
    print("=== Hybrid Dry MARTINI Workflow ===")
    print(f"Protein ID: {args.pdb_id}")
    print(f"Runtime PDB ID: {args.runtime_pdb_id}")
    print(f"Preparation seed: {args.prep_seed}")
    print(f"Simulation seed: {args.seed}")
    print(f"Run directory: {args.run_dir}")
    if source:
        print("Continuation mode: production only")
        print(f"Continuation source: {source}")
        print(f"Continuation output: {output}")


def run_soft_core_minimization(args, files):
    """CHARMM-GUI step6.0: steepest descent on the soft-core potential, which pushes template overlaps
    apart where the published potential's core is too steep to minimise. The protein carriers are held by a
    position restraint as well as moving as one rigid body."""
    prepare_workflow_hybrid_artifacts(args)
    print("=== Stage 6.0: soft-core minimization ===")
    prepare_stage_file(args, files["6.0"], "soft_core", 0, 0, "equilibration", "minimization")
    inject_protein_position_restraints(files["6.0"])
    run_minimization_stage(args, "6.0", files["6.0"], args.min_60_max_iter)
    remove_protein_position_restraints(files["6.0"])
    extract_stage_vtf(args, "6.0", files["6.0"], "1")


def run_bilayer_equilibration(args, files):
    """CHARMM-GUI step6.1-6.6: steepest descent on the published potential, then NPT dynamics under a
    phosphate height restraint released over five stages."""
    print("=== Stage 6.1: minimization ===")
    prepare_stage_file(args, files["6.1"], "hard", 1, 0, "equilibration", "minimization")
    set_initial_position(files["6.0"], files["6.1"])
    run_minimization_stage(args, "6.1", files["6.1"], args.min_61_max_iter)
    extract_stage_vtf(args, "6.1", files["6.1"], "1")

    previous = files["6.1"]
    for label, lipidhead_fc in EQUILIBRATION_STAGES:
        stage_file = files[label]
        prepare_stage_file(args, stage_file, "hard", 1, lipidhead_fc, "equilibration", "minimization")
        handoff_initial_position(previous, stage_file)
        nsteps = getattr(args, f"eq_{label.replace('.', '')}_nsteps")
        run_md_stage(args, label, stage_file, nsteps, args.eq_time_step, args.eq_frame_steps, label)
        extract_stage_vtf(args, label, stage_file, "1")
        previous = stage_file
    return previous


def assert_environment_solvation(args, up_file: Path):
    """Gate G2: the protein's lipid face must be solvated in the state that enters production.

    Measured on the last quarter of the rigid burn-in (the input state when the burn-in is skipped), on the belt
    sites that face lipid: the inside of a helix bundle is surrounded by protein, not lipid, so it says nothing about
    an insertion void. Two fractions of those sites are taken, averaged over the frames:

      * vacuum -- no environment bead within BELT_VACUUM_RADIUS_A. This is what an unhealed insertion void looks
        like, and it fails the build above BELT_VACUUM_MAX_FRACTION.
      * acyl-tail reach -- no acyl-tail bead within the same radius, the hydrophobic-mismatch signature that
        unfolded glpG TM4 against the DDM slab's headgroups. Reported, with a warning above the same fraction.
    """
    import h5py
    from scipy.spatial import cKDTree

    if args.belt_half_thickness_a is None:
        print(
            "Stage-7.0 solvation gate SKIPPED: the protein's hydrophobic belt is undefined without "
            "--opm-reference, so hydrophobic solvation was NOT verified for this build. Pass "
            "--opm-reference (a membrane-oriented homolog, midplane at z=0) to enable the check.",
            flush=True,
        )
        return

    with h5py.File(up_file, "r") as h5:
        n_frame = h5["/output/pos"].shape[0] if "/output/pos" in h5 else 0
        frames = (np.array(h5["/output/pos"][(3 * n_frame) // 4:, 0]) if n_frame
                  else np.array(h5["/input/pos"])[None, :, :, 0])
        env_index = np.array(h5["/input/hybrid_env_topology/env_atom_indices"])
        bb_index = np.array(h5["/input/hybrid_bb_map/bb_atom_index"])
        atom_names = np.array([n.decode().strip() for n in np.array(h5["/input/atom_names"])])
    apolar_names = derive_lipid_apolar_bead_names(args.upside_home, args.lipid_name)
    is_tail = np.isin(atom_names[env_index], sorted(apolar_names))
    if not is_tail.any():
        raise RuntimeError(
            f"No apolar environment beads {sorted(apolar_names)} found in {up_file.name}; cannot verify "
            "hydrophobic solvation."
        )
    protein = np.setdiff1d(np.arange(atom_names.size), env_index)
    backbone = protein[np.isin(atom_names[protein], ["N", "CA", "C", "O"])]
    half_thickness = float(args.belt_half_thickness_a)

    bare, dry, n_site = [], [], []
    for pos in frames:
        env_xyz, bb_xyz = pos[env_index], pos[bb_index]
        # the bilayer's own midplane, which the equilibration may have moved off the packed one
        midplane = env_xyz[:, 2].mean()
        offset = np.abs(bb_xyz[:, 2] - midplane)
        belt = bb_xyz[offset <= half_thickness]
        if belt.shape[0] == 0:
            raise RuntimeError(
                f"Stage-7.0 solvation gate found no backbone sites inside the belt: no BB site lies within "
                f"{half_thickness:.2f} A of the environment midplane z={midplane:.2f} A. Closest BB site is "
                f"{offset.min():.2f} A away and the BB span is z={bb_xyz[:, 2].min():.2f}..{bb_xyz[:, 2].max():.2f}, "
                "so the protein is not inserted where the OPM reference says its hydrophobic belt is; check the "
                "orientation and insertion depth."
            )
        burial = np.array([len(n) for n in cKDTree(pos[backbone]).query_ball_point(belt, BELT_BURIAL_RADIUS_A)])
        facing = belt[burial < BELT_FACING_MAX_BACKBONE]
        if facing.shape[0] == 0:
            raise RuntimeError(
                f"Stage-7.0 solvation gate found no lipid-facing belt site: each of the {belt.shape[0]} belt sites "
                f"has at least {BELT_FACING_MAX_BACKBONE} protein backbone atoms within {BELT_BURIAL_RADIUS_A:.0f} A."
            )
        bare.append(np.mean(cKDTree(env_xyz).query(facing)[0] > BELT_VACUUM_RADIUS_A))
        dry.append(np.mean(cKDTree(env_xyz[is_tail]).query(facing)[0] > BELT_VACUUM_RADIUS_A))
        n_site.append(facing.shape[0])
    bare, dry = float(np.mean(bare)), float(np.mean(dry))

    print(
        f"Stage-7.0 solvation gate over {len(frames)} frame(s): {np.mean(n_site):.0f} lipid-facing belt sites; "
        f"{100 * bare:.1f}% with no environment bead and {100 * dry:.1f}% with no acyl-tail bead within "
        f"{BELT_VACUUM_RADIUS_A:.1f} A; limit {100 * BELT_VACUUM_MAX_FRACTION:.0f}%"
    )
    if bare > BELT_VACUUM_MAX_FRACTION:
        raise RuntimeError(
            f"{100 * bare:.1f}% of the lipid-facing belt sites enter production with no environment bead within "
            f"{BELT_VACUUM_RADIUS_A:.1f} A (limit {100 * BELT_VACUUM_MAX_FRACTION:.0f}%). The lipid face would "
            "sample against vacuum; lengthen the 6.6 NPT equilibration or check the environment build."
        )
    if dry > BELT_VACUUM_MAX_FRACTION:
        print(
            f"  WARNING: {100 * dry:.1f}% of the lipid-facing belt sites have no acyl-tail bead within "
            f"{BELT_VACUUM_RADIUS_A:.1f} A. This is the hydrophobic-mismatch signature that unfolded glpG TM4 "
            "in the lamellar DDM slab; inspect before trusting the run."
        )


def run_stage70_handoff(args, files, source_stage: Path):
    """Stage 7.0: the production configuration. Minimise and burn in with the protein rigid, then run
    production with it flexible. The lipids are free in both: CHARMM-GUI's production carries no head-group
    restraint, and the 6.6 restraint (10 kJ/mol/nm^2) already leaves the phosphates about 5 A of play."""
    stage_file = files["7.0"]
    prepare_stage_file(args, stage_file, "hard", args.prod_70_npt_enable, 0, "production", "production_handoff")
    handoff_initial_position(source_stage, stage_file)
    assert_hybrid_stage_active(stage_file, "production_handoff", "production")
    run_minimization_stage(args, "7.0", stage_file, args.min_70_max_iter)
    run_stage70_burnin(args, stage_file)
    set_stage_label(stage_file, "production")
    assert_hybrid_stage_active(stage_file, "production", "production")
    if int(args.prod_70_nsteps) <= 0:
        print(f"=== Stage 7.0: production not run (nsteps <= 0); {stage_file} holds the production seed ===")
        return
    run_md_stage(args, "7.0", stage_file, args.prod_70_nsteps, args.prod_time_step, args.prod_frame_steps,
                 "production")
    extract_stage_vtf(args, "7.0", stage_file, "2")


def run_fresh_hybrid_workflow(args):
    files = workflow_stage_files(args)
    run_soft_core_minimization(args, files)
    run_stage70_handoff(args, files, run_bilayer_equilibration(args, files))


def run_hybrid_workflow_command(argv):
    args = parse_hybrid_workflow_args(argv)
    validate_hybrid_workflow_args(args)

    source, output, label = resolve_continuation_outputs(args)
    print_hybrid_workflow_summary(args, source, output)
    if source:
        run_stage70_continuation(args, source, output, label)
    else:
        run_fresh_hybrid_workflow(args)

    print("=== Workflow Complete ===")


if __name__ == "__main__":
    argv = sys.argv[1:]
    if argv and argv[0] == "prepare":
        run_prepare_command(argv[1:])
    elif argv and argv[0] == "run-hybrid-workflow":
        run_hybrid_workflow_command(argv[1:])
    else:
        raise SystemExit(
            "Unsupported command. Use:\n"
            "  martini_prepare_system.py prepare [options]\n"
            "  martini_prepare_system.py run-hybrid-workflow [options]"
        )
