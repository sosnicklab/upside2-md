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
from collections import Counter
from pathlib import Path

import numpy as np

from martini_itp_reader import lipid_max_sigma_nm, parse_dry_forcefield, parse_itp_file, parse_lipid_from_itp
from martini_prepare_system_lib import (
    build_backbone_with_virtual_bb,
    build_detergent_micelle,
    center_of_mass,
    convert_stage,
    compute_lipid_residue_indices,
    coords,
    estimate_salt_pairs,
    extract_protein_backbone_atoms_from_aa,
    infer_protein_charge_from_residues,
    inject_particles_table,
    inject_stage7_sc_table_nodes,
    lipid_resname,
    measure_belt_hydrophobic_coverage,
    membrane_slab_excluded_volume,
    set_active_lipid_name,
    map_backbone_types_from_martinize_fallback,
    map_backbone_types_from_structure,
    parse_pdb,
    place_ions,
    reject_inside_aggregate,
    reject_inside_slab,
    remove_overlapping_lipids,
    set_initial_position,
    set_box_for_micelle,
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
WORKFLOW_DIR = REPO_ROOT / "example" / "16.MARTINI"

DEFAULT_SC_ENV_BACKBONE_HOLD_STEPS = 200
DEFAULT_SC_ENV_PO4_Z_HOLD_STEPS = 150
DEFAULT_NPT_TAU = 4.0
DEFAULT_NPT_INTERVAL = 10
# Monte-Carlo barostat lateral step: exp(+-mc_dmax_xy) per trial (~1% box change), tuned for a
# reasonable acceptance rate on a dry-MARTINI membrane.
DEFAULT_NPT_MC_DMAX_XY = 0.01
DEFAULT_MARTINI_ENERGY_CONVERSION = 2.914952774272
DEFAULT_MARTINI_LENGTH_CONVERSION = 10.0
DEFAULT_COMPRESSIBILITY_3E4_BAR_INV_TO_A3_PER_EUP = 14.521180763676
MARTINI_MD_INTEGRATOR = "v"
MARTINI_MD_TIME_STEP = 0.009
# Mean amino-acid excluded volume (partial specific volume 0.73 cm^3/g, <MW> ~110 Da) used to
# subtract the protein from the solvent-accessible volume when counting ions.
PROTEIN_RESIDUE_VOLUME_A3 = 130.0


def find_lipid_itp_path(ff_dir: Path, lipid_name: str) -> Path:
    """Return the ITP under *ff_dir* that defines the *lipid_name* moleculetype."""
    target = lipid_name.strip().upper()
    for itp in sorted(Path(ff_dir).glob("*.itp")):
        try:
            topology = parse_itp_file(itp)
        except Exception:
            continue
        if any(mol.upper() == target for mol in topology.get("molecules", {})):
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


def derive_environment_morphology(upside_home: Path, lipid_name: str) -> str:
    """"micelle" for a single-tailed amphiphile, "bilayer" otherwise, read from the ITP topology.

    A one-chain amphiphile with a bulky head has a packing parameter far below 1, so it forms micelles,
    not lamellae: forcing it into a slab pins the area per molecule at the head's cross-section and
    starves the hydrophobic core (DDM: ~14 A against a ~28 A transmembrane belt). Deriving this from
    tail count rather than a name table keeps the rule universal and makes the unphysical combination
    impossible to request.
    """
    ff_dir = Path(upside_home) / "example" / "16.MARTINI" / "dryMARTINI_itp"
    tails = {}
    for name in (part.strip() for part in str(lipid_name).split(",") if part.strip()):
        lipid = parse_lipid_from_itp(find_lipid_itp_path(ff_dir, name), name)
        tails[name.upper()] = _count_hydrophobic_tails(lipid)
    if not tails:
        raise ValueError("Lipid name must be provided (--lipid-name)")
    single = {n for n, t in tails.items() if t <= 1}
    if single and len(single) != len(tails):
        raise ValueError(
            f"Cannot mix detergents {sorted(single)} with bilayer-forming lipids "
            f"{sorted(set(tails) - single)}: the aggregate morphology is ambiguous. Tails per "
            f"molecule: {tails}."
        )
    return "micelle" if single else "bilayer"


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
        description=(
            "Unified preparation script for bilayer-only, protein-only, or mixed "
            "protein+bilayer systems. Can optionally convert prepared structure "
            "to UPSIDE input."
        )
    )
    parser.add_argument("--mode", choices=["bilayer", "both"], required=True)
    parser.add_argument("--pdb-id", required=True, help="Runtime PDB id for stage conversion")
    parser.add_argument("--runtime-pdb-output", default=None)
    parser.add_argument("--prepare-structure", type=int, default=1, choices=[0, 1])
    parser.add_argument("--stage", default=None, help="stage name for UPSIDE input conversion")
    parser.add_argument("--run-dir", default="outputs/martini_test")

    parser.add_argument("--lipid-name", required=True,
                        help="Membrane lipid moleculetype name(s), comma-separated for a mixed bilayer "
                             "(e.g. DOPC, DDM, or POPE,POPG); topology read from the input ITP files")
    parser.add_argument("--bilayer-pdb", default=None,
                        help="Bilayer slab PDB for mode=bilayer (default: parameters/dryMARTINI/<lipid-name>.pdb)")
    parser.add_argument("--charmm-gui-dir", default=None,
                        help="CHARMM-GUI Martini job dir (mode=both): membrane read from "
                             "gromacs/step5_charmm2gmx.pdb + gromacs/system.top")
    parser.add_argument("--membrane-pdb", default=None,
                        help="CHARMM-GUI membrane CG PDB (overrides --charmm-gui-dir); must carry a CRYST1 box")
    parser.add_argument("--membrane-top", default=None,
                        help="CHARMM-GUI system.top for the lipid-count cross-check (overrides --charmm-gui-dir)")
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
    parser.add_argument("--protein-net-charge", type=int, default=None)
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


def runtime_paths(args):
    runtime_pdb = (
        Path(args.runtime_pdb_output).expanduser().resolve()
        if args.runtime_pdb_output
        else (WORKFLOW_DIR / "pdb" / f"{args.pdb_id}.MARTINI.pdb")
    )
    return runtime_pdb


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

    # xy fill window: the protein CG-envelope footprint + a PBC-safe belt on each side, rounded up to
    # a whole number of template tiles so the bilayer tiles seamlessly (no lateral gaps/contraction).
    fmin = cg_xyz.min(axis=0)
    fmax = cg_xyz.max(axis=0)
    pspan = fmax - fmin
    pcenter_xy = 0.5 * (fmin[:2] + fmax[:2])
    belt = float(args.box_padding_xy) + float(args.protein_pbc_margin)
    tile_xy = float(bilayer_box[0])
    target_side = max(tile_xy, float(np.max(pspan[:2] + 2.0 * belt)) * max(1.0, float(args.xy_scale)))
    target_side = float(np.ceil(target_side / tile_xy) * tile_xy)
    target_xy_min = np.array([pcenter_xy[0] - 0.5 * target_side, pcenter_xy[1] - 0.5 * target_side], dtype=float)
    target_xy_max = np.array([pcenter_xy[0] + 0.5 * target_side, pcenter_xy[1] + 0.5 * target_side], dtype=float)
    bilayer_lipids = tile_and_crop_bilayer_lipids(
        bilayer_atoms=lipid_template,
        bilayer_box=bilayer_box,
        target_xy_min=target_xy_min,
        target_xy_max=target_xy_max,
    )

    # Overlap removal against the CG envelope (BB+SC), not the backbone: whole lipids clashing with
    # the protein's true girth (including where the sidechain rotamers are rebuilt) are deleted.
    lipid_residues, keep_nonlipid = compute_lipid_residue_indices(bilayer_lipids)
    bilayer_kept, removed_lipids = remove_overlapping_lipids(
        bilayer_atoms=bilayer_lipids,
        protein_atoms=cg_atoms,
        lipid_residues=lipid_residues,
        keep_nonlipid=keep_nonlipid,
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
            else 0.5 * float(np.ptp(coords([a for a in bilayer_kept if lipid_resname(a["resname"])])[:, 2]))
        ),
    )

    packed_atoms = protein_backbone_atoms + bilayer_kept
    min_box_z_target = float(pspan[2] + 2.0 * (float(args.box_padding_z) + float(args.protein_pbc_margin)))
    box_lengths = set_box_from_lipid_xy(
        all_atoms=packed_atoms,
        lipid_atoms=bilayer_kept,
        pad_z=float(args.box_padding_z),
        force_square_xy=True,
        min_box_z=min_box_z_target,
        center_lipid_in_z=True,
        force_xy_box=(target_side, target_side),
    )

    kept_lipid_xyz = coords([a for a in bilayer_kept if lipid_resname(a["resname"])])
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
    # The belt is the protein's own hydrophobic band: the OPM reference publishes it, otherwise fall
    # back to the packed slab's extent.
    belt_half_thickness = (
        read_opm_hydrophobic_half_thickness(Path(args.opm_reference).expanduser().resolve())
        if args.opm_reference
        else 0.5 * packed_thickness_z
    )
    geometry = {
        "morphology": "bilayer",
        "lipid_residues_removed": int(removed_lipids),
        "belt_half_thickness_a": float(belt_half_thickness),
        "belt_coverage": coverage,
        "base_xy_side_angstrom": tile_xy,
        "target_xy_side_angstrom": target_side,
        "membrane_thickness_z_angstrom": packed_thickness_z,
        "ion_count_thickness_z_angstrom": count_thickness_z,
        "protein_span_angstrom": [float(v) for v in pspan],
        "excluded_volume_a3": membrane_slab_excluded_volume(
            box_lengths, count_thickness_z, 0.0
        ),
        "ion_reject": lambda box: reject_inside_slab(membrane_z_lo, membrane_z_hi),
    }
    return bilayer_kept, box_lengths, geometry


def pack_micelle_environment(args, lipid_template, bilayer_box, protein_backbone_atoms, cg_atoms, upside_home):
    """Wrap the protein's hydrophobic belt in a detergent monolayer (tails in, heads out)."""
    if not args.opm_reference:
        raise ValueError(
            "A micelle build needs --opm-reference: the protein's hydrophobic belt comes from the "
            "reference's OPM thickness REMARK."
        )
    if float(args.membrane_thickness_angstrom) > 0.0:
        raise ValueError(
            "--membrane-thickness-angstrom is a lamellar ion-counting input and does not apply to a "
            "micelle; the aggregate's own volume is measured instead."
        )
    belt_half_thickness = read_opm_hydrophobic_half_thickness(
        Path(args.opm_reference).expanduser().resolve()
    )
    ff_dir = Path(upside_home) / "example" / "16.MARTINI" / "dryMARTINI_itp"
    first_lipid = next(part.strip() for part in str(args.lipid_name).split(",") if part.strip())
    lipid = parse_lipid_from_itp(find_lipid_itp_path(ff_dir, first_lipid), first_lipid)

    micelle_atoms, diagnostics = build_detergent_micelle(
        template_atoms=lipid_template,
        template_box=bilayer_box,
        protein_cg_atoms=cg_atoms,
        belt_half_thickness=belt_half_thickness,
        contact_clearance=float(args.protein_lipid_cutoff),
        apolar_bead_indices=_apolar_bead_indices(lipid),
        rng=np.random.default_rng(int(args.seed)),
    )
    tail_covered = diagnostics["exposed_belt_beads"] - diagnostics["uncovered_exposed_belt_beads"]
    print(
        f"Micelle packing: {diagnostics['n_detergent']} {first_lipid} "
        f"({diagnostics['n_refilled']} refilled) around a {2.0 * belt_half_thickness:.1f} A "
        f"hydrophobic belt; {tail_covered}/{diagnostics['exposed_belt_beads']} exposed belt beads "
        f"tail-covered, {diagnostics['vacuum_exposed_belt_beads']} facing vacuum; bead spacing "
        f"min {diagnostics['min_intermolecular_a']:.2f} / median "
        f"{diagnostics['median_intermolecular_a']:.2f} A; min protein-detergent distance "
        f"{diagnostics['min_protein_detergent_distance_a']:.2f} A"
    )
    # Residual bald spots are REPORTED here, not fatal: a rigid conformer cannot always reach into a
    # deeply recessed site that a flexible chain would, and refining the grid only surfaces more such
    # sites rather than filling them. The authoritative solvation check is
    # assert_environment_solvation(), run on the equilibrated coordinates that actually enter
    # production -- by then the ladder has closed single-bead voids, and anything still bare is real.
    if diagnostics["vacuum_exposed_belt_beads"] > 0:
        print(
            f"  note: {diagnostics['vacuum_exposed_belt_beads']} exposed belt bead(s) start beyond "
            f"{diagnostics['coverage_radius_a']:.2f} A of any detergent bead; equilibration must close "
            "this, and the stage-7.0 solvation gate will verify it."
        )

    apolar_names = derive_lipid_apolar_bead_names(upside_home, args.lipid_name)
    # Measure belt coverage BEFORE the box shift, while z=0 is still the membrane midplane.
    coverage = measure_belt_hydrophobic_coverage(
        protein_backbone_atoms=protein_backbone_atoms,
        env_atoms=micelle_atoms,
        apolar_bead_names=apolar_names,
        belt_half_thickness=belt_half_thickness,
    )

    packed_atoms = protein_backbone_atoms + micelle_atoms
    box_lengths = set_box_for_micelle(
        packed_atoms, float(args.box_padding_z) + float(args.protein_pbc_margin)
    )
    tail_xyz = np.array(
        [[a["x"], a["y"], a["z"]] for a in micelle_atoms if a["name"].strip() in apolar_names],
        dtype=float,
    )
    head_xyz = np.array(
        [[a["x"], a["y"], a["z"]] for a in micelle_atoms if a["name"].strip() not in apolar_names],
        dtype=float,
    )
    cg_span = coords(cg_atoms).max(axis=0) - coords(cg_atoms).min(axis=0)
    geometry = {
        "morphology": "micelle",
        "lipid_residues_removed": 0,
        "belt_half_thickness_a": float(belt_half_thickness),
        "belt_coverage": coverage,
        "protein_span_angstrom": [float(v) for v in cg_span],
        "excluded_volume_a3": diagnostics["aggregate_volume_a3"],
        "ion_reject": lambda box: reject_inside_aggregate(tail_xyz, head_xyz, box),
    }
    geometry.update({f"micelle_{k}": v for k, v in diagnostics.items()})
    return micelle_atoms, box_lengths, geometry


def prepare_mixed_structure(args, runtime_pdb):
    if not args.protein_aa_pdb:
        raise ValueError("--protein-aa-pdb is required for mode=both")
    if not args.membrane_pdb:
        raise ValueError("--membrane-pdb (or --charmm-gui-dir) is required for mode=both")
    if not args.protein_cg_pdb:
        raise ValueError("--protein-cg-pdb (martinize CG envelope) is required for mode=both")

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
    prot_xyz = coords(protein_backbone_atoms)
    cg_xyz = coords(cg_atoms)

    upside_home = Path(os.environ.get("UPSIDE_HOME", str(REPO_ROOT)))
    morphology = derive_environment_morphology(upside_home, args.lipid_name)
    print(f"Environment morphology: {morphology} (from acyl-tail count in the lipid ITP)")
    if morphology == "micelle":
        bilayer_kept, box_lengths, env_geometry = pack_micelle_environment(
            args, lipid_template, bilayer_box, protein_backbone_atoms, cg_atoms, upside_home
        )
    else:
        bilayer_kept, box_lengths, env_geometry = pack_bilayer_environment(
            args, lipid_template, bilayer_box, protein_backbone_atoms, cg_atoms, prot_xyz, cg_xyz
        )
    packed_atoms = protein_backbone_atoms + bilayer_kept
    removed_lipids = env_geometry["lipid_residues_removed"]

    # Gate G1: the protein's hydrophobic belt must sit against acyl tails, not headgroups. A lamellar
    # DDM slab fails this (~14 A tail span, half the belt facing maltose), which is what unfolded TM4.
    # Reported, not gated: CHARMM-GUI templates are pre-minimization and laterally compressed, so a
    # packed-state span reads low for every lipid (DOPC 20 A, POPE/POPG 21 A) and cannot be compared with
    # OPM's relaxed hydrophobic thickness. assert_environment_solvation() does the real check per belt
    # residue on equilibrated coordinates, which is both local and like-for-like.
    coverage = env_geometry.pop("belt_coverage")
    print(
        f"Hydrophobic coverage (packed state): tail span {coverage['tail_hydrophobic_span_a']:.1f} A, "
        f"belt {coverage['belt_span_a']:.1f} A, {coverage['belt_polar_contact_fraction']:.2f} of "
        f"{coverage['belt_residues']} belt residues nearest a polar bead"
    )

    protein_charge = (
        int(args.protein_net_charge)
        if args.protein_net_charge is not None
        else int(infer_protein_charge_from_residues(protein_backbone_atoms))
    )
    # Charged lipids (e.g. POPG) contribute a net membrane charge that counterions must
    # neutralize alongside the protein; neutral bilayers add nothing.
    lipid_net_charges = derive_lipid_net_charges(
        Path(os.environ.get("UPSIDE_HOME", str(REPO_ROOT))), args.lipid_name
    )
    membrane_charge = 0
    seen_lipid_molecules = set()
    for atom in bilayer_kept:
        resname = atom["resname"].upper()
        if resname not in lipid_net_charges:
            continue
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key in seen_lipid_molecules:
            continue
        seen_lipid_molecules.add(key)
        membrane_charge += lipid_net_charges[resname]
    system_charge = protein_charge + membrane_charge

    # Ions occupy the solvent-accessible volume only (implicit solvent): the box minus whatever the
    # amphiphile aggregate and the protein exclude, and never inside the hydrophobic interior.
    n_protein_residues = len({(a["chain"], a["resseq"], a["icode"]) for a in protein_backbone_atoms})
    protein_volume_a3 = float(n_protein_residues) * PROTEIN_RESIDUE_VOLUME_A3
    if int(args.explicit_ions):
        salt_pairs = estimate_salt_pairs(
            box_lengths=box_lengths,
            salt_molar=float(args.salt_molar),
            excluded_volume_a3=env_geometry["excluded_volume_a3"] + protein_volume_a3,
        )
        n_na = int(salt_pairs + max(0, -system_charge))
        n_cl = int(salt_pairs + max(0, system_charge))
        rng = np.random.default_rng(int(args.seed))
        ion_atoms = place_ions(
            atoms=packed_atoms,
            box_lengths=box_lengths,
            n_na=n_na,
            n_cl=n_cl,
            cutoff=float(args.ion_cutoff),
            rng=rng,
            reject=env_geometry["ion_reject"](box_lengths),
        )
    else:
        salt_pairs = 0
        n_na = 0
        n_cl = 0
        ion_atoms = []

    (
        bb_type_by_residue,
        bb_secondary_by_residue,
        bb_type_source_by_residue,
    ) = map_backbone_types_from_structure(protein_backbone_atoms)
    protein_atoms, bb_entries = build_backbone_with_virtual_bb(
        protein_backbone_atoms,
        bb_type_by_residue,
        bb_secondary_by_residue=bb_secondary_by_residue,
        bb_type_source_by_residue=bb_type_source_by_residue,
    )
    all_atoms = protein_atoms + bilayer_kept + ion_atoms
    write_pdb(runtime_pdb, all_atoms, box_lengths)
    protein_xyz_final = coords(protein_backbone_atoms)
    bilayer_kept_lipid_atoms = [a for a in bilayer_kept if lipid_resname(a["resname"])]
    bilayer_xyz_final = coords(bilayer_kept_lipid_atoms) if bilayer_kept_lipid_atoms else coords(bilayer_kept)

    mapping_summary = {}
    if args.hybrid_mapping_output:
        mapping_h5 = Path(args.hybrid_mapping_output).expanduser().resolve()
        mapping_h5.parent.mkdir(parents=True, exist_ok=True)
        env_atom_indices = list(range(len(protein_atoms), len(all_atoms)))
        write_hybrid_mapping_h5(
            mapping_h5,
            bb_entries=bb_entries,
            total_martini_atoms=len(all_atoms),
            env_atom_indices=env_atom_indices,
            n_protein_atoms=len(protein_atoms),
        )

        mapping_summary["protein_aa_pdb"] = str(protein_aa_pdb)
        mapping_summary["mapping_h5"] = str(mapping_h5)
        mapping_summary["bb_structure_typed_residues"] = int(
            sum(1 for source in bb_type_source_by_residue.values() if source == "structure_geometry")
        )
        mapping_summary["bb_fallback_typed_residues"] = int(
            sum(1 for source in bb_type_source_by_residue.values() if source == "coil_geometry_fallback")
        )
        mapping_summary["bb_secondary_counts"] = dict(sorted(Counter(bb_secondary_by_residue.values()).items()))
        mapping_summary["bb_type_counts"] = dict(sorted(Counter(bb_type_by_residue.values()).items()))
        mapping_summary["bb_type_source_counts"] = dict(
            sorted(Counter(bb_type_source_by_residue.values()).items())
        )
        mapping_summary["bb_map_entries"] = int(len(bb_entries))

        if args.hybrid_bb_map_json_output:
            mapping_json = Path(args.hybrid_bb_map_json_output).expanduser().resolve()
            write_summary(mapping_json, {"bb_entries": bb_entries, "count": len(bb_entries)})
            mapping_summary["mapping_json"] = str(mapping_json)

    summary = {
        "mode": "both",
        "input_protein_aa_pdb": str(protein_aa_pdb),
        "input_membrane_pdb": str(membrane_pdb),
        "input_protein_cg_pdb": str(protein_cg_pdb),
        "runtime_pdb": str(runtime_pdb),
        "protein_orientation_mode": str(args.protein_orientation_mode),
        "opm_reference": str(args.opm_reference) if args.opm_reference else None,
        "xy_scale": float(args.xy_scale),
        "box_angstrom": [float(v) for v in box_lengths],
        "protein_charge_used": int(protein_charge),
        "membrane_charge": int(membrane_charge),
        "explicit_ions": int(args.explicit_ions),
        "salt_pairs_target": int(salt_pairs),
        "na_added": int(n_na),
        "cl_added": int(n_cl),
        "protein_atoms": int(len(protein_atoms)),
        "bilayer_atoms_kept": int(len(bilayer_kept)),
        "lipid_residues_removed": int(removed_lipids),
        "ion_atoms_added": int(len(ion_atoms)),
        "total_atoms": int(len(all_atoms)),
        "protein_final_z_min": float(protein_xyz_final[:, 2].min()),
        "protein_final_z_max": float(protein_xyz_final[:, 2].max()),
        "lipid_final_z_min": float(bilayer_xyz_final[:, 2].min()),
        "lipid_final_z_max": float(bilayer_xyz_final[:, 2].max()),
    }
    summary.update({k: v for k, v in env_geometry.items() if k != "ion_reject"})
    summary.update({f"belt_{k}": v for k, v in coverage.items()})
    summary.update(mapping_summary)
    return summary


def prepare_bilayer_structure(args, runtime_pdb):
    if not args.bilayer_pdb:
        raise ValueError("--bilayer-pdb is required for mode=bilayer")

    bilayer_pdb = Path(args.bilayer_pdb).expanduser().resolve()
    if not bilayer_pdb.exists():
        raise FileNotFoundError(f"Bilayer PDB not found: {bilayer_pdb}")

    bilayer_atoms, bilayer_box = parse_pdb(bilayer_pdb)
    bilayer_lipid_atoms = [dict(a) for a in bilayer_atoms if lipid_resname(a["resname"])]
    if not bilayer_lipid_atoms:
        raise ValueError("No lipid residues found in bilayer template.")

    # The template is a periodic tile, so its CRYST1 record IS the lateral box. Sizing from the coordinate
    # extent instead is wrong twice over: it loses the sub-bead gap between the outermost bead and the tile
    # edge, and it squares a rectangular tile up to its longer edge (measured 84.41 x 73.10 -> 84.26 x 84.26),
    # which opens a vacuum stripe along the shorter edge and changes the area per lipid by 15%.
    if bilayer_box is None:
        raise ValueError(f"Bilayer template {bilayer_pdb} has no CRYST1 record: the periodic tile's box "
                         "cannot be recovered from the coordinates alone.")
    lipid_xyz = coords(bilayer_lipid_atoms)
    bmin = lipid_xyz.min(axis=0)
    bmax = lipid_xyz.max(axis=0)
    tile_xy = np.array([float(bilayer_box[0]), float(bilayer_box[1])], dtype=float)
    target_xy = tile_xy * float(max(1.0, args.xy_scale))
    if float(args.xy_scale) > 1.000001:
        center_xy = 0.5 * (bmin[:2] + bmax[:2])
        target_xy_min = center_xy - 0.5 * target_xy
        target_xy_max = center_xy + 0.5 * target_xy
        bilayer_lipid_atoms = tile_and_crop_bilayer_lipids(
            bilayer_atoms=bilayer_atoms,
            bilayer_box=bilayer_box,
            target_xy_min=target_xy_min,
            target_xy_max=target_xy_max,
        )

    box_lengths = set_box_from_lipid_xy(
        all_atoms=bilayer_lipid_atoms,
        lipid_atoms=bilayer_lipid_atoms,
        pad_z=float(args.box_padding_z),
        force_square_xy=False,
        center_lipid_in_z=True,
        force_xy_box=target_xy,
    )

    if int(args.explicit_ions):
        lipid_z = coords(bilayer_lipid_atoms)[:, 2]
        membrane_thickness_z = float(lipid_z.max() - lipid_z.min())
        salt_pairs = estimate_salt_pairs(
            box_lengths=box_lengths,
            salt_molar=float(args.salt_molar),
            excluded_volume_a3=membrane_slab_excluded_volume(box_lengths, membrane_thickness_z, 0.0),
        )
        # An anionic lipid (POPG) leaves the membrane charged, so the counterions have to come from the
        # salt budget here too -- otherwise a pure POPG-containing bilayer is built non-neutral.
        lipid_net_charges = derive_lipid_net_charges(
            Path(os.environ.get("UPSIDE_HOME", str(REPO_ROOT))), args.lipid_name
        )
        membrane_charge = 0
        seen = set()
        for atom in bilayer_lipid_atoms:
            key = (atom["chain"], atom["resseq"], atom["icode"])
            if key in seen:
                continue
            seen.add(key)
            membrane_charge += lipid_net_charges.get(atom["resname"].upper(), 0)
        rng = np.random.default_rng(int(args.seed))
        ion_atoms = place_ions(
            atoms=bilayer_lipid_atoms,
            box_lengths=box_lengths,
            n_na=int(salt_pairs + max(0, -membrane_charge)),
            n_cl=int(salt_pairs + max(0, membrane_charge)),
            cutoff=float(args.ion_cutoff),
            rng=rng,
            reject=reject_inside_slab(float(lipid_z.min()), float(lipid_z.max())),
        )
    else:
        salt_pairs = 0
        ion_atoms = []

    all_atoms = bilayer_lipid_atoms + ion_atoms
    write_pdb(runtime_pdb, all_atoms, box_lengths)
    lipid_xyz_final = coords(bilayer_lipid_atoms)
    return {
        "mode": "bilayer",
        "input_bilayer_pdb": str(bilayer_pdb),
        "runtime_pdb": str(runtime_pdb),
        "xy_scale": float(args.xy_scale),
        "template_tile_xy_angstrom": [float(v) for v in tile_xy],
        "target_xy_angstrom": [float(v) for v in target_xy],
        "box_angstrom": [float(v) for v in box_lengths],
        "explicit_ions": int(args.explicit_ions),
        "salt_pairs_target": int(salt_pairs),
        "na_added": int(salt_pairs) if int(args.explicit_ions) else 0,
        "cl_added": int(salt_pairs) if int(args.explicit_ions) else 0,
        "bilayer_atoms_kept": int(len(bilayer_lipid_atoms)),
        "ion_atoms_added": int(len(ion_atoms)),
        "total_atoms": int(len(all_atoms)),
        "lipid_final_z_min": float(lipid_xyz_final[:, 2].min()),
        "lipid_final_z_max": float(lipid_xyz_final[:, 2].max()),
    }


def run_stage_conversion(args, runtime_pdb: Path):
    prev_pdb = os.environ.get("UPSIDE_RUNTIME_PDB_FILE")

    os.environ["UPSIDE_RUNTIME_PDB_FILE"] = str(runtime_pdb)
    os.environ.pop("UPSIDE_RUNTIME_ITP_FILE", None)

    try:
        convert_stage(
            pdb_id=args.pdb_id,
            stage=args.stage,
            run_dir=args.run_dir,
        )
    finally:
        if prev_pdb is None:
            os.environ.pop("UPSIDE_RUNTIME_PDB_FILE", None)
        else:
            os.environ["UPSIDE_RUNTIME_PDB_FILE"] = prev_pdb

        os.environ.pop("UPSIDE_RUNTIME_ITP_FILE", None)


def write_summary(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")


def resolve_charmm_gui_membrane_paths(args):
    """Fill --membrane-pdb/--membrane-top from --charmm-gui-dir (explicit paths win)."""
    if args.charmm_gui_dir:
        gmx = Path(args.charmm_gui_dir).expanduser() / "gromacs"
        if args.membrane_pdb is None:
            args.membrane_pdb = str(gmx / "step5_charmm2gmx.pdb")
        if args.membrane_top is None:
            args.membrane_top = str(gmx / "system.top")
    if args.membrane_pdb is None:
        raise ValueError("mode=both requires --charmm-gui-dir or --membrane-pdb")


def run_prepare_command(argv):
    args = parse_prepare_args(argv)
    set_active_lipid_name(args.lipid_name)
    if args.mode == "bilayer" and args.bilayer_pdb is None:
        args.bilayer_pdb = str(REPO_ROOT / "parameters" / "dryMARTINI" / f"{args.lipid_name}.pdb")
    if args.mode == "both" and args.prepare_structure:
        resolve_charmm_gui_membrane_paths(args)
    runtime_pdb = runtime_paths(args)
    runtime_pdb.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "mode": args.mode,
        "pdb_id": args.pdb_id,
        "prepare_structure": bool(args.prepare_structure),
        "stage": args.stage,
        "run_dir": args.run_dir,
    }

    if args.prepare_structure:
        if args.mode == "bilayer":
            summary.update(prepare_bilayer_structure(args, runtime_pdb))
        else:
            summary.update(prepare_mixed_structure(args, runtime_pdb))
    else:
        if not runtime_pdb.exists():
            raise FileNotFoundError(
                f"Runtime PDB not found for stage conversion: {runtime_pdb}. "
                "Run with --prepare-structure 1 first."
            )

    if args.stage:
        with temporary_env(
            {
                "UPSIDE_HOME": os.environ.get("UPSIDE_HOME", str(REPO_ROOT)),
                "UPSIDE_MARTINI_FF_DIR": os.environ.get(
                    "UPSIDE_MARTINI_FF_DIR",
                    str(REPO_ROOT / "example" / "16.MARTINI" / "dryMARTINI_itp"),
                ),
                "UPSIDE_MARTINI_ENERGY_CONVERSION": os.environ.get(
                    "UPSIDE_MARTINI_ENERGY_CONVERSION",
                    str(DEFAULT_MARTINI_ENERGY_CONVERSION),
                ),
                "UPSIDE_MARTINI_LENGTH_CONVERSION": os.environ.get(
                    "UPSIDE_MARTINI_LENGTH_CONVERSION",
                    str(DEFAULT_MARTINI_LENGTH_CONVERSION),
                ),
            }
        ):
            run_stage_conversion(args, runtime_pdb)
        target_file = Path(args.run_dir).expanduser().resolve() / "test.input.up"
        if args.mode == "bilayer":
            inject_particles_table(
                up_file=target_file,
                martini_h5=REPO_ROOT / "parameters" / "ff_2.1" / "martini.h5",
            )
        summary["upside_input"] = str(target_file)

    if args.summary_json:
        summary_path = Path(args.summary_json).expanduser().resolve()
    else:
        summary_path = runtime_pdb.with_suffix(".prep_summary.json")
    write_summary(summary_path, summary)
    print(f"Preparation summary written to: {summary_path}")


PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX",
}
SC_LIBRARY_REQUIRED_DATASETS = [
    "grid_ang",
    "cos_theta_grid",
    "rotamer_count",
    "rotamer_probability_fixed",
    "rotamer_radial_energy_eup",
    "rotamer_angular_energy_eup",
    "rotamer_angular_profile",
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


def env_default(name, default):
    return Config.text(name, default)


def env_int(name, default):
    return Config.integer(name, default)


def env_float(name, default):
    return Config.floating(name, default)


def workflow_path(value, base=WORKFLOW_DIR):
    path = Path(value).expanduser()
    return path if path.is_absolute() else base / path


def generate_random_seed():
    seed = random.SystemRandom().randint(1, 2**32 - 1)
    return seed if seed != 0 else 1


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
    protein_xyz = []
    lipid_xyz = []
    with pdb_file.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:21].strip().upper()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if resname in PROTEIN_RESIDUES:
                protein_xyz.append(xyz)
            elif lipid_resname(resname):
                lipid_xyz.append(xyz)
    if not protein_xyz or not lipid_xyz:
        return float("nan")
    protein = np.asarray(protein_xyz, dtype=float)
    lipid = np.asarray(lipid_xyz, dtype=float)
    min_d2 = float("inf")
    for i in range(0, lipid.shape[0], 2000):
        block = lipid[i:i + 2000]
        d = block[:, None, :] - protein[None, :, :]
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
            "--mode", "both",
            "--pdb-id", args.runtime_pdb_id,
            "--runtime-pdb-output", str(args.hybrid_packed_pdb),
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
    if not args.hybrid_packed_pdb.exists():
        raise FileNotFoundError(args.hybrid_packed_pdb)
    # Overlap removal against the CG envelope is the packing guarantee; the backbone-lipid distance
    # is a sanity log (it can sit near the clearance since a BB bead may be the closest CG bead).
    min_gap = pdb_min_protein_lipid_distance(args.hybrid_packed_pdb)
    print(
        f"Hybrid packing: min protein(backbone)-lipid distance {min_gap:.3f} A "
        f"(CG-envelope clearance {float(args.protein_lipid_cutoff):.3f} A)"
    )
    if np.isfinite(min_gap) and min_gap < 2.0:
        raise RuntimeError(f"Hybrid packing left a hard protein-lipid clash ({min_gap:.3f} A).")

    if not args.hybrid_mapping_file.exists():
        raise FileNotFoundError(args.hybrid_mapping_file)
    if args.hybrid_validate:
        validate_hybrid_mapping(args.hybrid_mapping_file)
    copy_if_different(args.hybrid_packed_pdb, args.runtime_pdb_file)
    print(f"Runtime MARTINI PDB: {args.runtime_pdb_file}")


def h5_as_text(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)


def inject_hybrid_mapping(up_file: Path, mapping_file: Path):
    import h5py

    groups = ["hybrid_control", "hybrid_bb_map", "hybrid_env_topology"]
    optional_groups = ("chain_break",)

    with h5py.File(mapping_file, "r") as src, h5py.File(up_file, "r+") as dst:
        src_inp = src["/input"]
        dst_inp = dst.require_group("input")
        base_n_atom = int(dst_inp["pos"].shape[0])
        src_mem = src_inp["hybrid_env_topology"]["protein_membership"][:].astype(np.int32)

        if base_n_atom != int(src_mem.shape[0]):
            raise ValueError(
                f"Hybrid mapping n_atom mismatch for {up_file}: up has {base_n_atom}, "
                f"mapping has {src_mem.shape[0]}"
            )

        for group in groups:
            if group not in src_inp:
                raise ValueError(f"Missing mapping group in {mapping_file}: /input/{group}")
            if group in dst_inp:
                del dst_inp[group]
            src.copy(src_inp[group], dst_inp, name=group)
        for group in optional_groups:
            if group not in src_inp:
                continue
            if group in dst_inp:
                del dst_inp[group]
            src.copy(src_inp[group], dst_inp, name=group)

        bb_grp = dst_inp["hybrid_bb_map"]
        env_grp = dst_inp["hybrid_env_topology"]
        if bb_grp["atom_indices"].shape[0] != bb_grp["bb_residue_index"].shape[0]:
            raise ValueError("hybrid_bb_map atom_indices/bb_residue_index size mismatch")
        if bb_grp["atom_indices"].shape[1] != 4:
            raise ValueError("hybrid_bb_map/atom_indices must have shape (n_bb,4)")

        bb_grp.attrs["atom_index_space"] = np.bytes_("stage_runtime")
        bb_grp.attrs["reference_index_space"] = np.bytes_("stage_runtime")
        bb_grp.attrs["reference_index_offset"] = np.int32(0)
        bb_grp.attrs["reference_index_count"] = np.int32(0)
        membership = env_grp["protein_membership"][:].astype(np.int32, copy=False)
        if membership.shape[0] != base_n_atom:
            raise ValueError(f"{up_file}: hybrid_env_topology/protein_membership length mismatch")


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
    import h5py

    with h5py.File(up_file, "r+") as h5:
        pm = h5["/input/hybrid_env_topology/protein_membership"][:]
        protein_atoms = np.where(pm >= 0)[0].astype(np.int32)
        if len(protein_atoms) == 0:
            return
        ref_pos = np.asarray(h5["/input/pos"][:][protein_atoms, :, 0], dtype=np.float32)
        grp = h5.require_group("/input/potential/restraint_position")
        grp.attrs["arguments"] = np.array([np.bytes_("pos")])
        for name in ("restraint_indices", "ref_pos", "spring_const", "spring_const_xyz"):
            if name in grp:
                del grp[name]
        grp.create_dataset("restraint_indices", data=protein_atoms)
        grp.create_dataset("ref_pos", data=ref_pos)
        spring_xyz = np.full((len(protein_atoms), 3), float(spring_const), dtype=np.float32)
        grp.create_dataset("spring_const_xyz", data=spring_xyz)


def remove_protein_position_restraints(up_file: Path):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        path = "/input/potential/restraint_position"
        if path in h5:
            del h5[path]


def set_hybrid_production_controls(up_file: Path, args):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        grp = h5.require_group("input").require_group("hybrid_control")
        grp.attrs["sc_env_backbone_hold_steps"] = np.int32(args.sc_env_backbone_hold_steps)
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


def stage_conversion_env(args, stage_label: str, prepare_stage: str, npt_enable: int, lipidhead_fc: float):
    if stage_label in {"production", "production_handoff"}:
        dynamics_phase = "production"
    elif prepare_stage == "npt_prod":
        dynamics_phase = "overlap_settling"
    else:
        dynamics_phase = "early_equilibration"
    # A finite aggregate has no lateral tension to relax, so box scaling is not a degree of freedom:
    # a tensionless xy barostat would just squeeze the box down onto the micelle. NVT throughout.
    if args.environment_morphology == "micelle":
        npt_enable = 0
    return {
        "UPSIDE_HOME": str(args.upside_home),
        "UPSIDE_MARTINI_FF_DIR": str(args.itp_dir),
        "UPSIDE_MARTINI_ENERGY_CONVERSION": str(args.martini_energy_conversion),
        "UPSIDE_MARTINI_LENGTH_CONVERSION": str(args.martini_length_conversion),
        "UPSIDE_SIMULATION_STAGE": prepare_stage,
        "UPSIDE_MARTINI_DYNAMICS_PHASE": dynamics_phase,
        "UPSIDE_NPT_ENABLE": str(int(npt_enable)),
        # Tensionless (zero lateral pressure) dry-MARTINI membrane barostat.
        "UPSIDE_NPT_TARGET_PXY": "0.0",
        "UPSIDE_NPT_TARGET_PZ": "0.0",
        "UPSIDE_NPT_TAU": str(args.npt_tau),
        "UPSIDE_NPT_COMPRESSIBILITY": str(args.compressibility_3e4_bar_inv_to_a3_per_eup),
        "UPSIDE_NPT_COMPRESSIBILITY_XY": str(args.compressibility_3e4_bar_inv_to_a3_per_eup),
        "UPSIDE_NPT_COMPRESSIBILITY_Z": "0.0",
        "UPSIDE_NPT_INTERVAL": str(args.npt_interval),
        "UPSIDE_NPT_SEMI": "1",
        # Monte-Carlo barostat steps: couple xy, freeze z (mc_dmax_z=0).
        "UPSIDE_NPT_MC_DMAX_XY": str(args.npt_mc_dmax_xy),
        "UPSIDE_NPT_MC_DMAX_Z": "0.0",
        "UPSIDE_NPT_MC_SEED": str(args.prep_seed),
        "UPSIDE_BILAYER_LIPIDHEAD_FC": str(lipidhead_fc),
        "THERMOSTAT_TIMESCALE": str(args.thermostat_timescale),
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


def prepare_stage_file(args, target_file: Path, prepare_stage: str, npt_enable: int, lipidhead_fc: float, stage_label: str):
    start = time.perf_counter()
    with temporary_env(stage_conversion_env(args, stage_label, prepare_stage, npt_enable, lipidhead_fc)):
        run_prepare_command(
            [
                "--mode", args.universal_prep_mode,
                "--lipid-name", args.lipid_name,
                "--pdb-id", args.runtime_pdb_id,
                "--runtime-pdb-output", str(args.runtime_pdb_file),
                "--prepare-structure", "0",
                "--stage", prepare_stage,
                "--run-dir", str(args.run_dir),
                "--summary-json", str(args.hybrid_prep_dir / f"stage_{prepare_stage}.summary.json"),
            ]
        )
    prepared_tmp = args.run_dir / "test.input.up"
    if not prepared_tmp.exists():
        raise FileNotFoundError(prepared_tmp)
    target_file.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(prepared_tmp), str(target_file))
    inject_hybrid_mapping(target_file, args.hybrid_mapping_file)
    set_hybrid_control_mode(target_file, args.hybrid_preprod_activation_stage)
    set_stage_label(target_file, stage_label)
    inject_particles_table(up_file=target_file, martini_h5=args.martini_h5)
    if stage_label in {"production", "production_handoff"}:
        inject_hybrid_interface_nodes(args, target_file, stage_label, "production")
    else:
        inject_hybrid_interface_nodes(
            args,
            target_file,
            stage_label,
            args.hybrid_preprod_activation_stage,
    )
    elapsed = time.perf_counter() - start
    print(
        f"[workflow timing] stage {stage_label} preparation/injection finished "
        f"in {elapsed:.3f} s -> {target_file}"
    )


def refresh_particle_mobility_friction(up_file: Path):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        if "/input/brownian" not in h5:
            return
        brownian = h5["/input/brownian"]
        observable = brownian.attrs.get("transport_observable", b"")
        if isinstance(observable, bytes):
            observable = observable.decode()
        if observable != "bare_martini_particle_lateral_diffusion":
            return

        atom_index = brownian["atom_index"][:]
        particle_class = h5["/input/particle_class"][:]
        atom_names = h5["/input/atom_names"][:]
        position = np.asarray(h5["/input/pos"][:], dtype=np.float64).reshape((-1, 3))
        lipid_beads = np.where(particle_class == b"LIPID")[0]
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
        contact_count = np.zeros(particle_class.size, dtype=np.int32)
        carrier_index = np.where(protein_carrier)[0]
        for start in range(0, lipid_beads.size, 512):
            delta = (
                position[carrier_index, None, :]
                - position[lipid_beads[start:start + 512]][None, :, :]
            )
            delta -= box * np.rint(delta / box)
            contact_count[carrier_index] += np.sum(
                np.sum(delta * delta, axis=2) < cutoff * cutoff,
                axis=1,
            )

        particle_friction = float(brownian.attrs["bare_particle_friction_up"])
        friction = np.full(atom_index.size, particle_friction, dtype=np.float32)
        protein_index = protein_carrier[atom_index]
        friction[protein_index] = (
            particle_friction * contact_count[atom_index[protein_index]]
        )
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


def handoff_initial_position(args, input_file: Path, output_file: Path, mode="default", previous_dt=None):
    preserve_transition = "1" if mode == "production_restart" and previous_dt is not None else "0"
    public_dt = float(previous_dt) if previous_dt is not None else 0.0
    with temporary_env(
        {
            "UPSIDE_SET_INITIAL_STRICT_COPY": str(args.strict_stage_handoff),
            "UPSIDE_SET_INITIAL_RECENTER_PRODUCTION": "0",
            "UPSIDE_SET_INITIAL_PRESERVE_HYBRID_TRANSITION": preserve_transition,
            "UPSIDE_SET_INITIAL_TIME_STEP": str(public_dt),
        }
    ):
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


def mark_output_restart_state(up_file: Path, nsteps: int, dt: float):
    import h5py

    public_dt = float(dt)
    expected_time = float(nsteps) * public_dt
    with h5py.File(up_file, "r+") as h5:
        if "/output/time" not in h5 or h5["/output/time"].shape[0] == 0:
            raise RuntimeError(f"{up_file} has no output/time; cannot mark restart state")
        if "/output/mom" not in h5 or h5["/output/mom"].shape[0] == 0:
            raise RuntimeError(f"{up_file} has no output/mom; cannot mark restart state")
        final_time = float(h5["/output/time"][-1])
        tolerance = max(1e-6, abs(public_dt) * 0.51)
        if abs(final_time - expected_time) > tolerance:
            raise RuntimeError(
                f"{up_file} final output time {final_time:.10g} does not match expected {expected_time:.10g}"
            )
        h5["/output/mom"].attrs["restart_final_state_valid"] = np.int8(1)
        h5["/output/mom"].attrs["restart_duration_steps"] = np.int64(nsteps)
        h5["/output/mom"].attrs["restart_time_step"] = float(dt)
        h5["/output/mom"].attrs["restart_public_time_step"] = public_dt
        h5["/output/mom"].attrs["restart_final_time"] = final_time


def run_minimization_stage(args, stage_label: str, up_file: Path, max_iter: int, preserve_stage: bool = False):
    if int(max_iter) <= 0:
        print(f"=== Stage {stage_label}: Minimization skipped (max_iter <= 0) ===")
        return
    print(f"=== Stage {stage_label}: Minimization ===")
    cmd = [
        args.upside_executable,
        up_file,
        "--duration", "0",
        "--frame-interval", "1",
        "--temperature", args.temperature,
        "--thermostat-timescale", args.thermostat_timescale,
        "--thermostat-interval", args.thermostat_interval,
        "--seed", args.seed,
        "--disable-recentering",
        "--minimize",
        "--min-max-iter", max_iter,
        "--min-energy-tol", "1e-6",
        "--min-force-tol", "1e-3",
        "--min-step", "0.01",
    ]
    if preserve_stage:
        cmd.append("--minimize-preserve-stage")
    run_checked(cmd, log_file=args.log_dir / f"stage_{stage_label}.min.log")
    promote_minimized_state_to_input(up_file)


def promote_minimized_state_to_input(up_file: Path):
    import h5py
    import numpy as np

    with h5py.File(up_file, "r+") as h5:
        if "/output/pos" not in h5 or h5["/output/pos"].shape[0] == 0:
            raise RuntimeError(f"{up_file} has no minimized output/pos to promote")
        last_pos = h5["/output/pos"][-1, 0, :, :][:, :, np.newaxis]
        if "/input/pos" in h5:
            del h5["/input/pos"]
        h5.create_dataset("/input/pos", data=last_pos)

        if "/input/mom" in h5:
            del h5["/input/mom"]
        mom = h5.create_dataset(
            "/input/mom",
            data=np.zeros(last_pos.shape, dtype=last_pos.dtype),
        )
        mom.attrs["restart_valid"] = np.int8(0)

        last_box = None
        if "/output/box" in h5 and h5["/output/box"].shape[0] > 0:
            last_box = np.asarray(h5["/output/box"][-1], dtype=float).reshape(-1)
        if last_box is not None and last_box.size >= 3 and "/input/potential" in h5:
            n_updated = 0
            for pot_grp in h5["/input/potential"].values():
                if not isinstance(pot_grp, h5py.Group):
                    continue
                if all(k in pot_grp.attrs for k in ("x_len", "y_len", "z_len")):
                    pot_grp.attrs["x_len"] = float(last_box[0])
                    pot_grp.attrs["y_len"] = float(last_box[1])
                    pot_grp.attrs["z_len"] = float(last_box[2])
                    n_updated += 1
            print(
                f"Promoted minimized state to input and updated {n_updated} potential node boxes: "
                f"x={last_box[0]:.3f}, y={last_box[1]:.3f}, z={last_box[2]:.3f}"
            )
        else:
            print("Promoted minimized state to input")
    refresh_particle_mobility_friction(up_file)


def run_md_stage(
    args,
    stage_label: str,
    input_file: Path,
    output_file: Path,
    nsteps: int,
    dt: float,
    frame_steps: int,
    echo: bool = True,
):
    effective_frame_steps = int(frame_steps)
    if effective_frame_steps >= int(nsteps):
        effective_frame_steps = max(1, int(nsteps) // 10)
        print(f"NOTICE: frame_steps ({frame_steps}) >= nsteps ({nsteps}); using frame_steps={effective_frame_steps}")
    frame_interval = f"{effective_frame_steps * float(dt):.10g}"
    if input_file.resolve() != output_file.resolve():
        shutil.copy2(input_file, output_file)
        handoff_initial_position(args, input_file, output_file)
    print(f"=== Stage {stage_label}: MD ===")
    cmd = [
        args.upside_executable,
        output_file,
        "--duration-steps", nsteps,
        "--frame-interval", frame_interval,
        "--temperature", args.temperature,
        "--time-step", dt,
        "--integrator", MARTINI_MD_INTEGRATOR,
        "--thermostat-timescale", args.thermostat_timescale,
        "--thermostat-interval", args.thermostat_interval,
        "--seed", args.seed,
        "--disable-recentering",
        "--record-momentum",
    ]
    if input_momentum_restart_valid(output_file):
        cmd.append("--restart-using-momentum")
    run_checked(cmd, log_file=args.log_dir / f"stage_{stage_label}.log", echo=echo)
    mark_output_restart_state(output_file, nsteps, dt)


def promote_md_output_state_to_input(up_file: Path, elapsed_steps: int):
    import h5py
    import numpy as np

    with h5py.File(up_file, "r+") as h5:
        if "/output/pos" not in h5 or h5["/output/pos"].shape[0] == 0:
            raise RuntimeError(f"{up_file} has no MD output/pos to promote")
        if "/output/mom" not in h5 or h5["/output/mom"].shape[0] == 0:
            raise RuntimeError(f"{up_file} has no MD output/mom to promote")

        last_pos = h5["/output/pos"][-1, 0, :, :][:, :, np.newaxis]
        last_mom = h5["/output/mom"][-1, 0, :, :][:, :, np.newaxis]
        if "/input/pos" in h5:
            del h5["/input/pos"]
        if "/input/mom" in h5:
            del h5["/input/mom"]
        h5.create_dataset("/input/pos", data=last_pos)
        mom = h5.create_dataset("/input/mom", data=last_mom)
        mom.attrs["restart_valid"] = np.int8(1)

        if "/input/hybrid_control" in h5:
            ctrl = h5["/input/hybrid_control"].attrs
            stage = h5_as_text(
                h5["/input/stage_parameters"].attrs.get("current_stage", b"")
            ).strip()
            if stage == "production":
                start = int(ctrl.get("sc_env_transition_step_start", 0))
                ctrl["sc_env_transition_step_start"] = np.int32(start + int(elapsed_steps))

        last_box = None
        if "/output/box" in h5 and h5["/output/box"].shape[0] > 0:
            last_box = np.asarray(h5["/output/box"][-1], dtype=float).reshape(-1)
        if last_box is not None and last_box.size >= 3 and "/input/potential" in h5:
            for pot_grp in h5["/input/potential"].values():
                if not isinstance(pot_grp, h5py.Group):
                    continue
                if all(k in pot_grp.attrs for k in ("x_len", "y_len", "z_len")):
                    pot_grp.attrs["x_len"] = float(last_box[0])
                    pot_grp.attrs["y_len"] = float(last_box[1])
                    pot_grp.attrs["z_len"] = float(last_box[2])

        del h5["/output"]
    refresh_particle_mobility_friction(up_file)


def run_stage70_burnin(args, stage_file: Path):
    nsteps = int(args.prod_70_burnin_nsteps)
    if nsteps <= 0:
        print("=== Stage 7.0: Rigid interface-handoff burn-in skipped (nsteps <= 0) ===")
        return
    print(f"=== Stage 7.0: Rigid interface-handoff burn-in ({nsteps} steps) ===")
    run_md_stage(
        args,
        "7.0.burnin",
        stage_file,
        stage_file,
        nsteps,
        args.prod_time_step,
        args.prod_frame_steps,
        echo=False,
    )
    promote_md_output_state_to_input(stage_file, nsteps)
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
    if not source_file.exists():
        raise FileNotFoundError(source_file)
    assert_hybrid_stage_active(source_file, "production", "production")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if source_file.resolve() != output_file.resolve():
        shutil.copy2(source_file, output_file)
    if not output_restart_state_valid(source_file):
        raise ValueError(
            "Production continuation requires a source stage generated by the current workflow with "
            "a validated final restart state. Regenerate the previous production stage with the current workflow."
    )
    handoff_initial_position(args, source_file, output_file, "production_restart", args.prod_time_step)
    inject_hybrid_interface_nodes(args, output_file, "production", "production")
    if not input_momentum_restart_valid(output_file):
        raise ValueError(
            "Production continuation requires restart-valid /output/mom in the source stage. "
            "Regenerate the previous production stage with the current workflow so momentum is recorded."
        )
    assert_hybrid_stage_active(output_file, "production", "production")
    run_md_stage(args, stage_label, output_file, output_file, args.prod_70_nsteps, args.prod_time_step, args.prod_frame_steps)
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
    args.universal_prep_mode = "both"
    args.hybrid_validate = True
    args.hybrid_preprod_activation_stage = "minimization"
    args.sc_env_backbone_hold_steps = DEFAULT_SC_ENV_BACKBONE_HOLD_STEPS
    args.sc_env_po4_z_hold_steps = DEFAULT_SC_ENV_PO4_Z_HOLD_STEPS
    args.npt_tau = DEFAULT_NPT_TAU
    args.npt_interval = DEFAULT_NPT_INTERVAL
    args.npt_mc_dmax_xy = DEFAULT_NPT_MC_DMAX_XY
    args.martini_energy_conversion = DEFAULT_MARTINI_ENERGY_CONVERSION
    args.martini_length_conversion = DEFAULT_MARTINI_LENGTH_CONVERSION
    args.compressibility_3e4_bar_inv_to_a3_per_eup = DEFAULT_COMPRESSIBILITY_3E4_BAR_INV_TO_A3_PER_EUP
    args.extract_vtf_script = workflow_path(args.extract_vtf_script).resolve()
    args.continue_stage_70_from = workflow_path(args.continue_stage_70_from).resolve() if args.continue_stage_70_from else None
    args.continue_stage_70_output = workflow_path(args.continue_stage_70_output).resolve() if args.continue_stage_70_output else None
    if args.prep_seed is None:
        args.prep_seed = generate_random_seed()
    if args.seed is None:
        args.seed = generate_random_seed()
        if args.seed == args.prep_seed:
            args.seed = generate_random_seed()
    for path in [args.run_dir, args.checkpoint_dir, args.log_dir, args.hybrid_prep_dir]:
        path.mkdir(parents=True, exist_ok=True)
    return args


def _decode_stage_strings(values):
    out = []
    for value in np.asarray(values):
        if isinstance(value, (bytes, np.bytes_)):
            out.append(value.decode("utf-8", errors="ignore").strip())
        else:
            out.append(str(value).strip())
    return np.asarray(out, dtype=object)


def _detect_has_bonded_environment_particles(up_file: Path) -> bool:
    """Return True when the stage has real bonded dry-MARTINI environment pairs."""
    import h5py

    with h5py.File(up_file, "r") as h5:
        bond_path = "/input/potential/dist_spring/id"
        if bond_path not in h5:
            print("  No dist_spring bonds found in stage input")
            print("  Extended pre-7.0 equilibrium stages needed: False")
            return False

        bonds = np.asarray(h5[bond_path][:], dtype=np.int64)
        if bonds.ndim != 2 or bonds.shape[1] < 2 or bonds.size == 0:
            print("  Empty dist_spring bond list in stage input")
            print("  Extended pre-7.0 equilibrium stages needed: False")
            return False

        inp = h5["/input"]
        atom_types = (
            _decode_stage_strings(inp["type"][:])
            if "type" in inp
            else np.asarray([""] * int(inp["pos"].shape[0]), dtype=object)
        )
        env_mask = None
        membership_path = "/input/hybrid_env_topology/protein_membership"
        if membership_path in h5:
            membership = np.asarray(h5[membership_path][:], dtype=np.int32)
            env_mask = membership < 0
        elif "particle_class" in inp:
            classes = np.char.upper(_decode_stage_strings(inp["particle_class"][:]).astype(str))
            env_mask = (classes != "PROTEIN") & (classes != "PROTEINAA")
        else:
            env_mask = np.ones(atom_types.shape[0], dtype=bool)

        bonded_env_pairs = []
        for raw_i, raw_j in bonds[:, :2]:
            i = int(raw_i)
            j = int(raw_j)
            if i < 0 or j < 0 or i >= atom_types.shape[0] or j >= atom_types.shape[0]:
                continue
            pair_types = {str(atom_types[i]).strip().upper(), str(atom_types[j]).strip().upper()}
            if bool(env_mask[i]) and bool(env_mask[j]):
                bonded_env_pairs.append((i, j, tuple(sorted(pair_types))))

    print(f"  Real bonded dry-MARTINI environment pairs: {len(bonded_env_pairs)}")
    if bonded_env_pairs:
        preview = ", ".join(
            f"{i}:{j}({'/'.join(pair_types)})" for i, j, pair_types in bonded_env_pairs[:5]
        )
        print(f"  Example bonded environment pairs: {preview}")
    print(f"  Extended pre-7.0 equilibrium stages needed: {bool(bonded_env_pairs)}")
    return bool(bonded_env_pairs)


def ensure_martini_parameter_libraries(args):
    """Check that the consolidated force-field file exists.

    If it is missing, print a clear error telling the user to run the
    generation script.
    """
    if not args.martini_h5.exists():
        raise RuntimeError(
            f"MARTINI force-field file not found: {args.martini_h5}"
            + "\nRun 'python py/martini_gen_params.py --upside-home "
            + str(args.upside_home) + "' to generate it."
        )
    print(f"MARTINI force-field library found: {args.martini_h5}")


def add_hybrid_workflow_arguments(parser):
    parser.add_argument("--pdb-id", default=env_default("PDB_ID", None))
    parser.add_argument("--runtime-pdb-id", default=env_default("RUNTIME_PDB_ID", None))
    parser.add_argument("--upside-home", default=env_default("UPSIDE_HOME", str(REPO_ROOT)))
    parser.add_argument("--run-dir", default=env_default("RUN_DIR", None))
    parser.add_argument("--protein-aa-pdb", default=env_default("PROTEIN_AA_PDB", None))
    parser.add_argument("--lipid-name", default=env_default("LIPID_NAME", None),
                        help="Membrane lipid moleculetype name(s), comma-separated for a mixed "
                             "bilayer (e.g. DOPC, DDM, or POPE,POPG)")
    parser.add_argument("--charmm-gui-dir", default=env_default("CHARMM_GUI_DIR", None),
                        help="CHARMM-GUI Martini job dir; membrane read from gromacs/step5_charmm2gmx.pdb + system.top")
    parser.add_argument("--membrane-pdb", default=env_default("MEMBRANE_PDB", None))
    parser.add_argument("--membrane-top", default=env_default("MEMBRANE_TOP", None))
    parser.add_argument("--opm-reference", default=env_default("OPM_REFERENCE", None),
                        help="OPM membrane-oriented reference PDB (for --protein-orientation-mode opm)")
    parser.add_argument("--extract-vtf-script", default=env_default("EXTRACT_VTF_SCRIPT", str(PY_DIR / "martini_extract_vtf.py")))
    parser.add_argument("--salt-molar", type=float, default=env_float("SALT_MOLAR", 0.15))
    parser.add_argument("--explicit-ions", type=int, choices=[0, 1], default=env_int("EXPLICIT_IONS", 1))
    parser.add_argument("--protein-lipid-cutoff", type=float, default=env_float("PROTEIN_LIPID_CUTOFF", 0.0))
    parser.add_argument("--ion-cutoff", type=float, default=env_float("ION_CUTOFF", 10.0))
    parser.add_argument("--xy-scale", type=float, default=env_float("XY_SCALE", 1.0))
    parser.add_argument("--box-padding-xy", type=float, default=env_float("BOX_PADDING_XY", 0.0))
    parser.add_argument("--box-padding-z", type=float, default=env_float("BOX_PADDING_Z", 0.0))
    parser.add_argument("--membrane-thickness-angstrom", type=float, default=env_float("MEMBRANE_THICKNESS_ANGSTROM", 0.0),
                        help="Equilibrated dry-MARTINI membrane thickness for the ion count "
                             "(0 = measure from the compressed packed lipids)")
    parser.add_argument("--protein-orientation-mode", choices=["input", "lay-flat", "opm"], default=env_default("PROTEIN_ORIENTATION_MODE", "opm"))
    parser.add_argument("--protein-pbc-margin", type=float, default=env_float("PROTEIN_PBC_MARGIN", 15.0))
    parser.add_argument("--temperature", type=float, default=env_float("TEMPERATURE", 0.8647))
    parser.add_argument("--thermostat-timescale", type=float, default=env_float("THERMOSTAT_TIMESCALE", 5.0))
    parser.add_argument("--thermostat-interval", type=int, default=env_int("THERMOSTAT_INTERVAL", -1))
    parser.add_argument("--strict-stage-handoff", type=int, default=env_int("STRICT_STAGE_HANDOFF", 1))
    parser.add_argument("--min-60-max-iter", type=int, default=env_int("MIN_60_MAX_ITER", 500))
    parser.add_argument("--min-61-max-iter", type=int, default=env_int("MIN_61_MAX_ITER", 0))
    parser.add_argument("--min-70-max-iter", type=int, default=env_int("MIN_70_MAX_ITER", 500))
    parser.add_argument("--eq-60-nsteps", type=int, default=env_int("EQ_60_NSTEPS", 500))
    parser.add_argument("--eq-62-nsteps", type=int, default=env_int("EQ_62_NSTEPS", 500))
    parser.add_argument("--eq-63-nsteps", type=int, default=env_int("EQ_63_NSTEPS", 500))
    parser.add_argument("--eq-64-nsteps", type=int, default=env_int("EQ_64_NSTEPS", 500))
    parser.add_argument("--eq-65-nsteps", type=int, default=env_int("EQ_65_NSTEPS", 500))
    parser.add_argument("--eq-66-nsteps", type=int, default=env_int("EQ_66_NSTEPS", 500))
    parser.add_argument("--prod-70-burnin-nsteps", type=int, default=env_int("PROD_70_BURNIN_NSTEPS", 40000))
    parser.add_argument("--prod-70-nsteps", type=int, default=env_int("PROD_70_NSTEPS", 10000))
    parser.add_argument("--eq-time-step", type=float, default=env_float("EQ_TIME_STEP", MARTINI_MD_TIME_STEP))
    parser.add_argument("--prod-time-step", type=float, default=env_float("PROD_TIME_STEP", MARTINI_MD_TIME_STEP))
    parser.add_argument("--eq-frame-steps", type=int, default=env_int("EQ_FRAME_STEPS", 1000))
    parser.add_argument("--prod-frame-steps", type=int, default=env_int("PROD_FRAME_STEPS", 50))
    parser.add_argument("--prod-70-npt-enable", type=int, default=env_int("PROD_70_NPT_ENABLE", 0))
    parser.add_argument("--prep-seed", default=os.environ.get("PREP_SEED"))
    parser.add_argument("--seed", default=os.environ.get("SEED"))
    parser.add_argument("--continue-stage-70-from", default=env_default("CONTINUE_STAGE_70_FROM", ""))
    parser.add_argument("--continue-stage-70-output", default=env_default("CONTINUE_STAGE_70_OUTPUT", ""))
    parser.add_argument("--continue-stage-70-label", default=env_default("CONTINUE_STAGE_70_LABEL", ""))


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
        raise ValueError("A lipid name is required (--lipid-name, e.g. DOPC or DDM)")
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
    # The aggregate morphology follows from the amphiphile's tail count and decides whether the
    # barostat has any meaning, so resolve it once here and carry it on args.
    args.environment_morphology = derive_environment_morphology(args.upside_home, args.lipid_name)
    if args.environment_morphology == "micelle" and not args.opm_reference:
        raise ValueError(
            f"{args.lipid_name} is a single-tailed detergent, so the environment is built as a micelle, "
            "which needs --opm-reference for the protein's hydrophobic belt."
        )
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
    return {
        "prepared_60": args.checkpoint_dir / f"{args.pdb_id}.stage_6.0.prepared.up",
        "stage_60": args.checkpoint_dir / f"{args.pdb_id}.stage_6.0.up",
        "prepared_61": args.checkpoint_dir / f"{args.pdb_id}.stage_6.1.prepared.up",
        "stage_61": args.checkpoint_dir / f"{args.pdb_id}.stage_6.1.up",
        "stage_62": args.checkpoint_dir / f"{args.pdb_id}.stage_6.2.up",
        "prepared_63": args.checkpoint_dir / f"{args.pdb_id}.stage_6.3.prepared.up",
        "stage_63": args.checkpoint_dir / f"{args.pdb_id}.stage_6.3.up",
        "prepared_64": args.checkpoint_dir / f"{args.pdb_id}.stage_6.4.prepared.up",
        "stage_64": args.checkpoint_dir / f"{args.pdb_id}.stage_6.4.up",
        "prepared_65": args.checkpoint_dir / f"{args.pdb_id}.stage_6.5.prepared.up",
        "stage_65": args.checkpoint_dir / f"{args.pdb_id}.stage_6.5.up",
        "prepared_66": args.checkpoint_dir / f"{args.pdb_id}.stage_6.6.prepared.up",
        "stage_66": args.checkpoint_dir / f"{args.pdb_id}.stage_6.6.up",
        "prepared_70": args.checkpoint_dir / f"{args.pdb_id}.stage_7.0.prepared.up",
        "stage_70": args.checkpoint_dir / f"{args.pdb_id}.stage_7.0.up",
    }


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


def run_stage60_relaxation(args, files):
    prepare_workflow_hybrid_artifacts(args)
    print("=== Stage 6.0: rigid-protein NPT box relaxation ===")
    prepare_stage_file(args, files["prepared_60"], "npt_equil", 1, 0, "minimization")
    shutil.copy2(files["prepared_60"], files["stage_60"])
    inject_protein_position_restraints(files["stage_60"])
    run_minimization_stage(args, "6.0", files["stage_60"], args.min_60_max_iter)
    remove_protein_position_restraints(files["stage_60"])
    run_md_stage(
        args,
        "6.0",
        files["stage_60"],
        files["stage_60"],
        args.eq_60_nsteps,
        args.eq_time_step,
        args.eq_frame_steps,
    )
    extract_stage_vtf(args, "6.0", files["stage_60"], "1")
    return files["stage_60"]


def assert_environment_solvation(args, up_file: Path):
    """Gate G2: the protein must be properly solvated in the state that actually enters production.

    Run on the equilibrated stage-7.0 coordinates, because that is the configuration production samples,
    and per belt residue rather than as a global thickness, because that is the only like-for-like
    comparison (a packed-state span reads low for every lipid). Two distinct failures are caught:

      * vacuum -- a belt site with no environment bead at all in reach. This is what an unhealed
        insertion void looks like; the lamellar glpG build entered production with the first shell empty
        (mean nearest environment bead 11.7 A over transmembrane CA).
      * headgroup solvation -- a belt site whose only company is polar. This is the hydrophobic-mismatch
        failure that unfolded TM4: the DDM slab's 14 A tail core left half the belt against maltose, and
        it is invisible to a vacuum test because the headgroups are right there.

    The limit for both is the derived coverage radius the builder uses: twice the force field's contact
    distance, i.e. a gap wide enough for another bead to sit in.
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

    clearance = float(args.protein_lipid_cutoff)
    coverage_radius = 2.0 * clearance
    with h5py.File(up_file, "r") as h5:
        pos = np.array(h5["/input/pos"])[:, :, 0]
        env_index = np.array(h5["/input/hybrid_env_topology/env_atom_indices"])
        bb_index = np.array(h5["/input/hybrid_bb_map/bb_atom_index"])
        atom_names = np.array([n.decode().strip() for n in np.array(h5["/input/atom_names"])])
    env_xyz = pos[env_index]
    apolar_names = derive_lipid_apolar_bead_names(args.upside_home, args.lipid_name)
    is_tail = np.isin(atom_names[env_index], sorted(apolar_names))
    if not is_tail.any():
        raise RuntimeError(
            f"No apolar environment beads {sorted(apolar_names)} found in {up_file.name}; cannot verify "
            "hydrophobic solvation."
        )
    # The aggregate defines its own midplane; a finite micelle is free to drift off the box centre.
    bb_xyz = pos[bb_index]
    half_thickness = float(args.belt_half_thickness_a)
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

    d_any = cKDTree(env_xyz).query(belt)[0]
    d_tail = cKDTree(env_xyz[is_tail]).query(belt)[0]
    bare = int(np.count_nonzero(d_any > coverage_radius))
    dry = int(np.count_nonzero(d_tail > coverage_radius))
    # Local hydrophobic thickness: the z-extent of tail beads in a cylinder around the protein, i.e. the
    # core the belt actually sits in, measured after equilibration so it is not the compressed template.
    lateral = np.linalg.norm(env_xyz[:, :2] - belt[:, :2].mean(axis=0), axis=1)
    local_tails = env_xyz[is_tail & (lateral <= np.linalg.norm(belt[:, :2].std(axis=0)) + 3.0 * clearance)]
    local_span = (
        float(np.percentile(local_tails[:, 2], 95) - np.percentile(local_tails[:, 2], 5))
        if local_tails.shape[0] > 10
        else float("nan")
    )
    belt_span = 2.0 * float(args.belt_half_thickness_a)
    print(
        f"Stage-7.0 solvation gate: {belt.shape[0]} belt sites; nearest environment bead mean "
        f"{d_any.mean():.2f} / max {d_any.max():.2f} A ({bare} bare); nearest acyl tail mean "
        f"{d_tail.mean():.2f} / max {d_tail.max():.2f} A ({dry} beyond reach); local tail core "
        f"{local_span:.1f} A vs {belt_span:.1f} A belt; limit {coverage_radius:.2f} A"
    )
    # Only the vacuum test gates. Tail reach and local core thickness are reported because both are
    # ambiguous on a post-equilibration snapshot: once a helix has already been stripped, detergent flows
    # into the cavity and the statistics recover, so they diagnose rather than decide. Promoting the
    # thickness check to a hard gate needs the DOPC/POPC bilayer regression first (plan.md P7) -- a
    # packed-state span reads ~20 A for every lipid, so an unvalidated threshold would break the
    # lamellar path.
    if bare:
        raise RuntimeError(
            f"{bare} of {belt.shape[0]} hydrophobic-belt sites enter production with no environment bead "
            f"within {coverage_radius:.2f} A. The belt would sample against vacuum; lengthen "
            "equilibration or check the environment build."
        )
    if dry or (local_span == local_span and local_span < belt_span - clearance):
        print(
            f"  WARNING: {dry} belt site(s) beyond acyl-tail reach and a {local_span:.1f} A local tail "
            f"core against a {belt_span:.1f} A belt. This is the hydrophobic-mismatch signature that "
            "unfolded glpG TM4 in the lamellar DDM slab; inspect before trusting the run."
        )


def run_stage70_handoff(args, files, source_stage: Path):
    prepare_stage_file(
        args,
        files["prepared_70"],
        "npt_prod",
        args.prod_70_npt_enable,
        0,
        "production_handoff",
    )
    shutil.copy2(files["prepared_70"], files["stage_70"])
    handoff_initial_position(args, source_stage, files["stage_70"], "production_hybrid")
    assert_hybrid_stage_active(files["stage_70"], "production_handoff", "production")
    run_minimization_stage(args, "7.0", files["stage_70"], args.min_70_max_iter, preserve_stage=True)
    run_stage70_burnin(args, files["stage_70"])
    set_stage_label(files["stage_70"], "production")
    assert_hybrid_stage_active(files["stage_70"], "production", "production")
    assert_environment_solvation(args, files["stage_70"])
    run_md_stage(
        args,
        "7.0",
        files["stage_70"],
        files["stage_70"],
        args.prod_70_nsteps,
        args.prod_time_step,
        args.prod_frame_steps,
    )
    extract_stage_vtf(args, "7.0", files["stage_70"], "2")


def run_pre70_equilibration(args, files):
    print("=== Bonded dry-MARTINI environment detected -> running extended pre-7.0 equilibrium ===")

    prepare_stage_file(args, files["prepared_61"], "npt_prod", 1, 0, "minimization")
    shutil.copy2(files["prepared_61"], files["stage_61"])
    handoff_initial_position(args, files["stage_60"], files["stage_61"])
    run_minimization_stage(args, "6.1", files["stage_61"], args.min_61_max_iter)
    extract_stage_vtf(args, "6.1", files["stage_61"], "1")

    prepare_stage_file(args, files["stage_62"], "npt_equil", 1, 200, "minimization")
    handoff_initial_position(args, files["stage_61"], files["stage_62"])
    run_md_stage(args, "6.2", files["stage_62"], files["stage_62"], args.eq_62_nsteps, args.eq_time_step, args.eq_frame_steps)
    extract_stage_vtf(args, "6.2", files["stage_62"], "1")

    prepare_stage_file(args, files["prepared_63"], "npt_equil_reduced", 1, 100, "minimization")
    shutil.copy2(files["prepared_63"], files["stage_63"])
    handoff_initial_position(args, files["stage_62"], files["stage_63"])
    run_md_stage(args, "6.3", files["stage_63"], files["stage_63"], args.eq_63_nsteps, args.eq_time_step, args.eq_frame_steps)
    extract_stage_vtf(args, "6.3", files["stage_63"], "1")

    prepare_stage_file(args, files["prepared_64"], "npt_prod", 1, 50, "minimization")
    shutil.copy2(files["prepared_64"], files["stage_64"])
    handoff_initial_position(args, files["stage_63"], files["stage_64"])
    run_md_stage(args, "6.4", files["stage_64"], files["stage_64"], args.eq_64_nsteps, args.eq_time_step, args.eq_frame_steps)
    extract_stage_vtf(args, "6.4", files["stage_64"], "1")

    prepare_stage_file(args, files["prepared_65"], "npt_prod", 1, 20, "minimization")
    shutil.copy2(files["prepared_65"], files["stage_65"])
    handoff_initial_position(args, files["stage_64"], files["stage_65"])
    run_md_stage(args, "6.5", files["stage_65"], files["stage_65"], args.eq_65_nsteps, args.eq_time_step, args.eq_frame_steps)
    extract_stage_vtf(args, "6.5", files["stage_65"], "1")

    prepare_stage_file(args, files["prepared_66"], "npt_prod", 1, 10, "minimization")
    shutil.copy2(files["prepared_66"], files["stage_66"])
    handoff_initial_position(args, files["stage_65"], files["stage_66"])
    run_md_stage(args, "6.6", files["stage_66"], files["stage_66"], args.eq_66_nsteps, args.eq_time_step, args.eq_frame_steps)
    extract_stage_vtf(args, "6.6", files["stage_66"], "1")
    return files["stage_66"]


def run_fresh_hybrid_workflow(args):
    ensure_martini_parameter_libraries(args)
    files = workflow_stage_files(args)
    stage60 = run_stage60_relaxation(args, files)

    if _detect_has_bonded_environment_particles(stage60):
        stage70_source = run_pre70_equilibration(args, files)
    else:
        print("=== No bonded dry-MARTINI environment pairs -> handoff from stage 6.0 to stage 7.0 ===")
        stage70_source = stage60

    run_stage70_handoff(args, files, stage70_source)


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
