#!/usr/bin/env python

from __future__ import annotations

import _pickle as cPickle
import importlib.util
import os
import types
import warnings
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path

import math

import h5py
import numpy as np
import tables as tb

from martini_itp_reader import parse_dry_forcefield, parse_itp_atomtype_masses, parse_itp_file, sc_target_label

# NumPy moved VisibleDeprecationWarning under np.exceptions in 1.25 and dropped the
# top-level alias in 2.0. The clusters run 1.23, this machine runs 2.x.
VISIBLE_DEPRECATION_WARNING = getattr(np, "exceptions", np).VisibleDeprecationWarning

PY_DIR = Path(__file__).resolve().parent
REPO_ROOT = PY_DIR.parent
WORKFLOW_DIR = REPO_ROOT / "example" / "16.MARTINI"

NA_AVOGADRO = 6.02214076e23
BB_COMPONENT_NAMES = ("N", "CA", "C", "O")
BB_COMPONENT_MASSES = (14.0, 12.0, 12.0, 16.0)
MARTINIZE_BB_SS_ORDER = ("F", "E", "H", "1", "2", "3", "T", "S", "C")
MARTINIZE_BB_DEFAULT = ("N0", "Nda", "N0", "Nd", "Na", "Nda", "Nda", "P5", "P5")
MARTINIZE_BB_RESIDUE_OVERRIDE = {
    "ALA": ("C5", "N0", "C5", "N0", "N0", "N0", "N0", "P4", "P4"),
    "PRO": ("C5", "N0", "C5", "N0", "Na", "N0", "N0", "Na", "Na"),
    "HYP": ("C5", "N0", "C5", "N0", "N0", "N0", "N0", "Na", "Na"),
}
BB_TYPE_CHARGE = {
    "Qd": 1.0,
    "Qa": -1.0,
    "SQd": 1.0,
    "SQa": -1.0,
    "RQd": 1.0,
    "AQa": -1.0,
}
# Electric conversion factor f = 1/(4 pi eps0) in kJ mol^-1 nm e^-2 and dry MARTINI's relative permittivity.
COULOMB_F_KJ_NM = 138.935458
DRY_MARTINI_EPSILON_R = 15.0
# The soft-core form CHARMM-GUI's step 6.0 minimisation uses (GROMACS free-energy = yes, init-lambda = 0.01,
# sc-alpha = 4, sc-power = 2, sc-coul = yes). Every other stage runs the tables unsoftened.
SOFT_CORE_ALPHA = 4.0
SOFT_CORE_POWER = 2
STAGE_SOFT_CORE_LAMBDA = {
    "soft_core": 0.01,
    "hard": 0.0,
}
# Measured lateral diffusion of each lipid at the friction reference temperature, in um^2/s, keyed by residue
# name. In production every bead of a lipid gets the friction kT_ref / D_bead,up that delivers its own
# molecule's target through the bead clock (findings 2.5). A lipid absent from the table is refused.
LIPID_LATERAL_DIFFUSION_UM2_S = {
    # pfg-NMR at 303 K (Filippov, Oradd and Lindblom, Langmuir 2003, 19:6397), as quoted in Lipid14's Table 9
    # (Dickson et al., JCTC 2014, 10:865). The same group reads 9.8 at 30 wt% water (Biophys J 2003, 84:3079).
    "DOPC": 11.5,
    # No study measures POPE and DOPC by one method. MAS pfg-NMR from one lab puts a fluid PE-rich bilayer
    # (SOPC/POPE 3:7, 8 +- 1 at 306.6 K; Polozov and Gawrisch, Biophys J 2004, 87:1741) at 0.75 +- 0.15 of POPC
    # (19 +- 1 at 322 K, Ea 30.7 kJ/mol; Gaede and Gawrisch, Biophys J 2003, 85:1734), and POPC reads as DOPC
    # within 7% (Filippov 2003): 0.75 x 11.5. Pure POPE at 303 K is 6 K above its gel transition and may be slower.
    "POPE": 8.6,
    # No free-bilayer POPG value at physiological salt: PG diffuses faster than PC at low salt (DOPC/DOPG 3:1 at
    # 1.28x DOPC by GUV FCS; Ramadurai et al., Biophys J 2010, 99:1482) and is slowed by salt (Filippov et al.,
    # Chem Phys Lipids 2009, 159:81), so its ratio to DOPC is taken as 1.0 +- 0.3.
    "POPG": 11.5,
}
# Deliberate interim rule: ions are not lipids and have no measured target here, so they keep the DOPC
# reference clock.
ION_CLOCK_REFERENCE_LIPID = "DOPC"


def lipid_diffusion_targets(lipid_names):
    """Each lipid's measured lateral diffusion (um^2/s); refuse a lipid the table does not hold."""
    missing = sorted(n for n in lipid_names if n.upper() not in LIPID_LATERAL_DIFFUSION_UM2_S)
    if missing:
        raise ValueError(
            f"No measured lateral diffusion for {missing}: add each lipid's value, with its citation, to "
            "LIPID_LATERAL_DIFFUSION_UM2_S in py/martini_prepare_system_lib.py before building it"
        )
    return {n.upper(): LIPID_LATERAL_DIFFUSION_UM2_S[n.upper()] for n in lipid_names}

CANONICAL_AFFINE_REF = np.array(
    [
        [-1.19280531, -0.83127186, 0.0],
        [0.0, 0.0, 0.0],
        [1.25222632, -0.87268266, 0.0],
    ],
    dtype=np.float32,
)
CANONICAL_AFFINE_REF -= CANONICAL_AFFINE_REF.mean(axis=0, keepdims=True)

# CB placement in the affine_alignment frame. That frame's origin is the CENTROID of N/CA/C
# (src/eig.cpp centres the three atoms before the Kabsch fit, and the stored ref_geom is centred to
# match), so a point placed in it must be given relative to the centroid -- which is what
# upside_config.write_environment does. Storing the raw CB instead displaces every sidechain site by the
# centroid, 0.568 A, and CB is the anchor of the whole SC-env term (martini_sc_table_1body takes
# placement_fixed_point_vector_only_CB as its sidechain input). The direction vector is CB-CA, a
# difference, so the centroid cancels there and it is unchanged.
CB_REFERENCE_GEOMETRY = np.array([
    [-1.19280531, -0.83127186, 0.0],        # N
    [ 0.0,         0.0,        0.0],        # CA
    [ 1.25222632, -0.87268266, 0.0],        # C
    [ 0.0,         0.94375626, 1.2068012],  # CB
])
_cb_frame_origin = CB_REFERENCE_GEOMETRY[:3].mean(axis=0)
CB_PLACEMENT = (CB_REFERENCE_GEOMETRY[3] - _cb_frame_origin)[None, :].astype(np.float32)
_cb_direction = CB_REFERENCE_GEOMETRY[3] - CB_REFERENCE_GEOMETRY[1]
CB_VECTOR = (_cb_direction / np.linalg.norm(_cb_direction))[None, :].astype(np.float32)
N_BIT_ROTAMER = 4
# MARTINI 2.2 maps these residues to a backbone bead alone, so they have no side-chain/environment table.
NO_MARTINI_SIDECHAIN = ("ALA", "GLY")
BACKBONE_NODES = [
    "Distance3D",
    "Angle",
    "Dihedral_omega",
    "Dihedral_phi",
    "Dihedral_psi",
    "Spring_bond",
    "Spring_angle",
    "Spring_omega",
    "rama_coord",
    "rama_map_pot",
    "rama_map_pot_ref",
    "infer_H_O",
    "protein_hbond",
    "hbond_energy",
    "backbone_pairs",
]


def parse_pdb(path: Path):
    atoms = []
    cryst1 = None
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            record = line[:6].strip()
            if record == "CRYST1":
                cryst1 = (
                    float(line[6:15]),
                    float(line[15:24]),
                    float(line[24:33]),
                )
                continue
            if record not in {"ATOM", "HETATM"}:
                continue

            atom = {
                "record": record,
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "resname": line[17:21].strip(),
                "chain": line[21].strip(),
                "resseq": int(line[22:26]),
                "icode": line[26].strip(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "occ": float(line[54:60]) if line[54:60].strip() else 1.0,
                "bfac": float(line[60:66]) if line[60:66].strip() else 0.0,
                "segid": line[72:76].strip() if len(line) >= 76 else "",
                "element": line[76:78].strip() if len(line) >= 78 else "",
                "charge": line[78:80].strip() if len(line) >= 80 else "",
            }
            atoms.append(atom)
    return atoms, cryst1


def write_pdb(path: Path, atoms, box_lengths):
    with path.open("w", encoding="utf-8") as f:
        if box_lengths is not None:
            f.write(
                f"CRYST1{box_lengths[0]:9.3f}{box_lengths[1]:9.3f}{box_lengths[2]:9.3f}"
                "  90.00  90.00  90.00 P 1           1\n"
            )
        for idx, atom in enumerate(atoms, start=1):
            chain = atom["chain"] if atom["chain"] else " "
            segid = atom["segid"][:4].ljust(4) if atom["segid"] else "    "
            f.write(
                f"{atom['record']:<6}{idx:5d} {atom['name'][:4]:>4} {atom['resname'][:4]:>4}{chain:1}"
                f"{atom['resseq']:4d}{atom['icode'][:1]:1}   "
                f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
                f"{atom['occ']:6.2f}{atom['bfac']:6.2f}      {segid}"
                f"{atom['element'][:2]:>2}{atom['charge'][:2]:>2}\n"
            )
        f.write("END\n")


def _minimum_image_delta(delta, box_lengths):
    delta = np.asarray(delta, dtype=np.float64).copy()
    if box_lengths is None:
        return delta
    box = np.asarray(box_lengths, dtype=np.float64)
    valid = box > 0.0
    delta[..., valid] -= box[valid] * np.round(delta[..., valid] / box[valid])
    return delta


def coords(atoms):
    return np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)


def set_coords(atoms, xyz):
    for atom, c in zip(atoms, xyz):
        atom["x"], atom["y"], atom["z"] = float(c[0]), float(c[1]), float(c[2])


def center_of_mass(xyz):
    return np.mean(xyz, axis=0)


# Membrane molecule selection. The lipid identity is supplied by the caller
# (--lipid-name, comma-separated for a mixed bilayer); each type's topology is read
# from the input ITP files. No lipid type, ITP filename, or residue name is hardcoded
# here. Names are stored uppercased so residue-name comparisons are case-insensitive.
_ACTIVE_LIPID_NAMES = ()


def set_active_lipid_name(name) -> None:
    global _ACTIVE_LIPID_NAMES
    if name is None:
        raise ValueError("Lipid name must be provided (--lipid-name)")
    parts = str(name).split(",") if isinstance(name, str) else [str(p) for p in name]
    names = tuple(dict.fromkeys(p.strip().upper() for p in parts if p.strip()))
    if not names:
        raise ValueError("Lipid name must be provided (--lipid-name)")
    _ACTIVE_LIPID_NAMES = names


def active_lipid_names() -> tuple:
    if not _ACTIVE_LIPID_NAMES:
        raise ValueError("Active lipid name is not set; call set_active_lipid_name() (--lipid-name)")
    return _ACTIVE_LIPID_NAMES


def active_lipid_name() -> str:
    return "/".join(active_lipid_names())


def lipid_resname(resname: str) -> bool:
    return resname.strip().upper() in active_lipid_names()


def residue_key(atom):
    return (atom["chain"], atom["resseq"], atom["icode"], atom["resname"].upper())


def extract_protein_backbone_atoms_from_aa(aa_atoms):
    protein_res = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "HID",
        "HIE",
        "HIP",
        "HSD",
        "HSE",
        "HSP",
        "CYX",
    }
    backbone_roles = {"N", "CA", "C", "O"}
    by_residue = defaultdict(dict)
    for atom in aa_atoms:
        resname = atom["resname"].upper()
        if resname not in protein_res:
            continue
        role = atom["name"].strip().upper()
        if role not in backbone_roles:
            continue
        key = residue_key(atom)
        if role not in by_residue[key]:
            by_residue[key][role] = deepcopy(atom)

    ordered_keys = sorted(by_residue.keys(), key=lambda x: (x[0], x[1], x[2], x[3]))
    out = []
    residue_index = []
    for rid, key in enumerate(ordered_keys, start=1):
        role_map = by_residue[key]
        if not all(role in role_map for role in BB_COMPONENT_NAMES):
            continue
        for role in BB_COMPONENT_NAMES:
            atom = role_map[role]
            atom["name"] = role
            atom["segid"] = "PROA"
            out.append(atom)
            residue_index.append(rid)
    if not out:
        raise ValueError("No complete protein backbone residues (N/CA/C/O) found in AA PDB.")
    return out, np.asarray(residue_index, dtype=np.int32)


def _complete_backbone_residue_groups(backbone_atoms):
    residue_groups = defaultdict(dict)
    residue_order = []
    for atom in backbone_atoms:
        key = residue_key(atom)
        role = atom["name"].strip().upper()
        if key not in residue_groups:
            residue_order.append(key)
        residue_groups[key][role] = atom

    valid_order = []
    valid_groups = {}
    for key in residue_order:
        group = residue_groups[key]
        if not all(role in group for role in BB_COMPONENT_NAMES):
            continue
        valid_order.append(key)
        valid_groups[key] = group
    if not valid_order:
        raise ValueError("Cannot map BB types: no complete N/CA/C/O residue rows.")
    return valid_order, valid_groups


def _atom_xyz(atom):
    return np.array([float(atom["x"]), float(atom["y"]), float(atom["z"])], dtype=np.float64)


def _martinize_secondary_structure(residue_order, residue_groups):
    """Each residue's martinize secondary-structure class, as `martinize -dssp` assigns it: DSSP (mdtraj's
    implementation) on each chain's backbone, translated and pattern-typed by martinize's own ssClassification
    (helix starts and ends, strands E, turns, bends, coil)."""
    import mdtraj
    from martinize import ssClassification
    labels = {}
    for chain in dict.fromkeys(key[0] for key in residue_order):
        keys = [key for key in residue_order if key[0] == chain]
        top = mdtraj.Topology()
        top_chain = top.add_chain()
        xyz = []
        for key in keys:
            residue = top.add_residue(key[3], top_chain, resSeq=key[1])
            for name in BB_COMPONENT_NAMES:
                top.add_atom(name, mdtraj.element.get_by_symbol(name[0]), residue)
                xyz.append(0.1 * _atom_xyz(residue_groups[key][name]))   # A to nm
        dssp = mdtraj.compute_dssp(mdtraj.Trajectory(np.array(xyz)[None], top), simplified=False)[0]
        labels.update(zip(keys, ssClassification("".join(dssp), "dssp")[1]))
    return labels


def _martinize_bb_type_for_residue(resname, secondary_structure):
    table = MARTINIZE_BB_RESIDUE_OVERRIDE.get(resname.upper(), MARTINIZE_BB_DEFAULT)
    try:
        return table[MARTINIZE_BB_SS_ORDER.index(str(secondary_structure).upper())]
    except ValueError:
        return table[MARTINIZE_BB_SS_ORDER.index("C")]


def map_backbone_types_from_structure(backbone_atoms):
    residue_order, residue_groups = _complete_backbone_residue_groups(backbone_atoms)
    secondary_by_residue = _martinize_secondary_structure(residue_order, residue_groups)
    type_by_residue = {}
    source_by_residue = {}
    for key in residue_order:
        type_by_residue[key] = _martinize_bb_type_for_residue(key[3], secondary_by_residue[key])
        source_by_residue[key] = "dssp"

    if residue_order:
        terminal_overrides = [(residue_order[0], "Qd"), (residue_order[-1], "Qa")]
        for i in range(1, len(residue_order)):
            prev = residue_order[i - 1]
            curr = residue_order[i]
            if curr[0] != prev[0]:
                terminal_overrides.append((curr, "Qd"))
                terminal_overrides.append((prev, "Qa"))
        for key, bb_type in terminal_overrides:
            type_by_residue[key] = bb_type
            source_by_residue[key] = "terminal_charge_override"

    return type_by_residue, secondary_by_residue, source_by_residue


def compute_lipid_residue_indices(bilayer_atoms):
    lipid_residues = defaultdict(list)
    for idx, atom in enumerate(bilayer_atoms):
        key = (atom["chain"], atom["resseq"], atom["icode"])
        lipid_residues[key].append(idx)
    return lipid_residues


def residue_group_atoms(atoms):
    groups = defaultdict(list)
    for atom in atoms:
        key = (atom["chain"], atom["resseq"], atom["icode"], atom["resname"].upper())
        groups[key].append(atom)
    return groups


CHARMM_GUI_WATER_RESNAMES = {"W", "WF"}
CHARMM_GUI_ION_RESNAMES = {"NA", "CL", "SOD", "CLA", "ION"}


def _parse_top_molecule_counts(top_path: Path):
    """Return {RESNAME: count} from the [ molecules ] section of a GROMACS system.top."""
    counts = {}
    in_molecules = False
    with Path(top_path).open("r", encoding="utf-8") as f:
        for line in f:
            text = line.split(";", 1)[0].strip()
            if not text:
                continue
            if text.startswith("["):
                in_molecules = text.replace(" ", "").lower() == "[molecules]"
                continue
            if in_molecules:
                parts = text.split()
                if len(parts) >= 2:
                    counts[parts[0].upper()] = counts.get(parts[0].upper(), 0) + int(parts[1])
    return counts


def read_charmm_gui_membrane(membrane_pdb, membrane_top=None):
    """Read a CHARMM-GUI Martini membrane as the bilayer input for the hybrid system.

    Keeps only the active-lipid residues; water (implicit solvent in dry-MARTINI) is stripped and
    CHARMM-GUI ions are dropped (ions are regenerated from concentration during packing). A CRYST1
    box is required (the gromacs/step5_charmm2gmx.pdb carries it; step5_assembly.pdb does not).
    When a system.top is given, the per-lipid molecule counts are cross-checked against its
    [ molecules ] section. Returns (lipid_atoms, box).
    """
    atoms, box = parse_pdb(Path(membrane_pdb))
    if box is None:
        raise ValueError(
            f"CHARMM-GUI membrane PDB has no CRYST1 box: {membrane_pdb}. "
            "Use the gromacs/step5_charmm2gmx.pdb (it carries the box), not step5_assembly.pdb."
        )
    lipid_atoms = []
    unexpected = set()
    for atom in atoms:
        resname = atom["resname"].strip().upper()
        if lipid_resname(resname):
            lipid_atoms.append(atom)
        elif resname in CHARMM_GUI_WATER_RESNAMES or resname in CHARMM_GUI_ION_RESNAMES:
            continue
        else:
            unexpected.add(resname)
    if unexpected:
        raise ValueError(
            f"Unexpected residues {sorted(unexpected)} in CHARMM-GUI membrane {membrane_pdb}. "
            f"Only active lipids {list(active_lipid_names())}, water, and ions are supported."
        )
    if not lipid_atoms:
        raise ValueError(
            f"No active-lipid residues {list(active_lipid_names())} found in {membrane_pdb}."
        )
    if membrane_top is not None:
        expected = _parse_top_molecule_counts(membrane_top)
        got = {}
        prev_key = None
        for atom in lipid_atoms:
            key = (atom["chain"], atom["resseq"], atom["icode"])
            if key != prev_key:
                resname = atom["resname"].strip().upper()
                got[resname] = got.get(resname, 0) + 1
                prev_key = key
        for name in active_lipid_names():
            if name in expected and expected[name] != got.get(name, 0):
                raise ValueError(
                    f"Lipid count mismatch for {name}: {membrane_pdb} has {got.get(name, 0)} "
                    f"molecules but {membrane_top} [ molecules ] declares {expected[name]}."
                )
    return lipid_atoms, box


def tile_and_crop_bilayer_lipids(
    bilayer_atoms,
    bilayer_box,
    target_xy_min,
    target_xy_max,
):
    # The membrane is read lipids-only (ions are regenerated after packing); its CRYST1 box is the tile.
    lipid_atoms = [a for a in bilayer_atoms if lipid_resname(a["resname"])]
    if not lipid_atoms:
        raise ValueError("No lipid residues found in bilayer input for tiling/cropping.")
    tile_x, tile_y = float(bilayer_box[0]), float(bilayer_box[1])
    if tile_x <= 0.0 or tile_y <= 0.0:
        raise ValueError("Invalid bilayer tile dimensions for XY tiling.")

    target_xy_min = np.array(target_xy_min, dtype=float)
    target_xy_max = np.array(target_xy_max, dtype=float)
    target_center = 0.5 * (target_xy_min + target_xy_max)
    target_span = target_xy_max - target_xy_min

    # Symmetric tiling around the bilayer center.
    # +2 guard ensures crop coverage despite residue COM discretization.
    nx = int(np.ceil(target_span[0] / tile_x)) + 2
    ny = int(np.ceil(target_span[1] / tile_y)) + 2

    base_groups = residue_group_atoms(lipid_atoms)
    tiled_groups = []
    for ix in range(-nx, nx + 1):
        for iy in range(-ny, ny + 1):
            shift = np.array([ix * tile_x, iy * tile_y, 0.0], dtype=float)
            for group in base_groups.values():
                gcopy = []
                for atom in group:
                    a = deepcopy(atom)
                    a["x"] = float(a["x"] + shift[0])
                    a["y"] = float(a["y"] + shift[1])
                    gcopy.append(a)
                tiled_groups.append(gcopy)

    # Recenter tiled lattice so its center matches target center.
    tiled_xyz = np.array(
        [[a["x"], a["y"], a["z"]] for group in tiled_groups for a in group],
        dtype=float,
    )
    tiled_center_xy = center_of_mass(tiled_xyz)[:2]
    recenter_shift_xy = target_center - tiled_center_xy
    for group in tiled_groups:
        for atom in group:
            atom["x"] = float(atom["x"] + recenter_shift_xy[0])
            atom["y"] = float(atom["y"] + recenter_shift_xy[1])

    # Crop by residue center to preserve complete molecules.
    cropped_groups = []
    for group in tiled_groups:
        gxyz = np.array([[a["x"], a["y"]] for a in group], dtype=float)
        gcom = gxyz.mean(axis=0)
        if (
            target_xy_min[0] <= gcom[0] <= target_xy_max[0]
            and target_xy_min[1] <= gcom[1] <= target_xy_max[1]
        ):
            cropped_groups.append(group)

    if not cropped_groups:
        raise RuntimeError("Bilayer tiling/cropping produced no lipid residues.")

    # Renumber residues for unique molecule ids downstream.
    out_atoms = []
    next_resseq = 1
    for group in cropped_groups:
        chain = group[0]["chain"]
        icode = group[0]["icode"]
        for atom in group:
            atom["chain"] = chain
            atom["icode"] = icode
            atom["resseq"] = next_resseq
            atom["segid"] = "MEMB"
            out_atoms.append(atom)
        next_resseq += 1
    return out_atoms


def measure_belt_hydrophobic_coverage(protein_backbone_atoms, env_atoms, apolar_bead_names, belt_half_thickness):
    """Gate G1: is the protein's hydrophobic belt actually against acyl tails?

    Returns the belt-residue fraction whose nearest environment bead is a polar head bead, and the
    tail-bead hydrophobic span. A lamellar DDM slab fails this badly (~14 A span, half the belt facing
    maltose), which is what unfolded glpG TM4.
    """
    env_xyz = coords(env_atoms)
    is_tail = np.array([a["name"].strip() in apolar_bead_names for a in env_atoms], dtype=bool)
    if not is_tail.any():
        raise ValueError("No apolar environment beads found while measuring belt coverage.")

    belt_ca = np.array(
        [
            [a["x"], a["y"], a["z"]]
            for a in protein_backbone_atoms
            if a["name"].strip() == "CA" and abs(float(a["z"])) <= float(belt_half_thickness)
        ],
        dtype=float,
    )
    if belt_ca.shape[0] == 0:
        raise ValueError("No backbone CA inside the hydrophobic belt while measuring coverage.")

    nearest = np.argmin(np.linalg.norm(belt_ca[:, None, :] - env_xyz[None, :, :], axis=2), axis=1)
    tail_z = env_xyz[is_tail][:, 2]
    return {
        "belt_residues": int(belt_ca.shape[0]),
        "belt_polar_contact_fraction": float(np.mean(~is_tail[nearest])),
        "tail_hydrophobic_span_a": float(np.percentile(tail_z, 95) - np.percentile(tail_z, 5)),
        "belt_span_a": 2.0 * float(belt_half_thickness),
    }


def remove_overlapping_lipids(
    bilayer_atoms,
    protein_atoms,
    lipid_residues,
    cutoff,
):
    protein_xyz = coords(protein_atoms)
    cutoff2 = cutoff * cutoff
    keep_atom_idx = set()
    removed_residues = 0
    for res_key, atom_idx_list in lipid_residues.items():
        lipid_xyz = np.array(
            [[bilayer_atoms[i]["x"], bilayer_atoms[i]["y"], bilayer_atoms[i]["z"]] for i in atom_idx_list],
            dtype=float,
        )
        delta = lipid_xyz[:, None, :] - protein_xyz[None, :, :]
        dist2 = np.sum(delta * delta, axis=2)
        if np.min(dist2) < cutoff2:
            removed_residues += 1
            continue
        for i in atom_idx_list:
            keep_atom_idx.add(i)
    kept = [bilayer_atoms[i] for i in sorted(keep_atom_idx)]
    return kept, removed_residues


def set_box_from_lipid_xy(all_atoms, lipid_atoms, xy_box, pad_z, min_box_z):
    """Box for a periodic lamellar fill: lateral edges equal the whole-tile window `xy_box` (x, y), z holds
    every atom with `pad_z` on each side and the bilayer midplane at box_z/2. The lipid COM is centred
    laterally; molecule tails past an edge are the periodic images of the tiling."""
    all_xyz = coords(all_atoms)
    lip_xyz = coords(lipid_atoms)
    box_x, box_y = (float(v) for v in xy_box)
    lip_mid_z = 0.5 * float(lip_xyz[:, 2].min() + lip_xyz[:, 2].max())
    min_z = float(all_xyz[:, 2].min())
    max_z = float(all_xyz[:, 2].max())
    box_z = max(max_z - min_z + 2.0 * float(pad_z), float(min_box_z),
                2.0 * max(lip_mid_z - min_z, max_z - lip_mid_z))

    lip_com = lip_xyz.mean(axis=0)
    shift = np.array([0.5 * box_x - lip_com[0], 0.5 * box_y - lip_com[1], 0.5 * box_z - lip_mid_z])
    set_coords(all_atoms, all_xyz + shift)
    return np.array([box_x, box_y, box_z], dtype=float)


def estimate_salt_pairs(box_lengths, salt_molar, excluded_volume_a3):
    # dry-MARTINI is implicit-solvent, so ions occupy only the solvent-accessible volume: the box
    # minus whatever the amphiphile aggregate and the protein exclude. Using the full box would
    # massively overcount (the box is mostly aggregate + the protein's soluble-domain space).
    box_volume_a3 = float(box_lengths[0]) * float(box_lengths[1]) * float(box_lengths[2])
    volume_a3 = max(0.0, box_volume_a3 - float(excluded_volume_a3))
    pairs = int(round(salt_molar * NA_AVOGADRO * volume_a3 * 1e-27))
    return max(0, pairs)


def membrane_slab_excluded_volume(box_lengths, membrane_thickness_z):
    """Volume ions cannot occupy in a lamellar system: the slab, A_xy * thickness."""
    return float(box_lengths[0]) * float(box_lengths[1]) * float(membrane_thickness_z)


def reject_inside_slab(z_lo, z_hi):
    """Ion veto for a lamellar membrane: no ion between the two headgroup planes."""
    def reject(trial):
        return float(z_lo) <= float(trial[2]) <= float(z_hi)
    return reject


def place_ions(atoms, box_lengths, n_na, n_cl, cutoff, rng, reject=None):
    # reject(trial_xyz) -> bool vetoes a trial that lands in the amphiphile's hydrophobic interior;
    # dry-MARTINI has no explicit water, so nothing else keeps ions out of it.
    existing = coords(atoms)
    placed = []
    cutoff2 = cutoff * cutoff
    box = np.asarray(box_lengths, dtype=float)
    types = [("NA", "NA", 1.0, "IONS")] * n_na + [("CL", "CL", -1.0, "IONS")] * n_cl
    for name, resname, _charge, segid in types:
        accepted = False
        for _ in range(20000):
            trial = rng.uniform([0, 0, 0], box)
            if reject is not None and reject(trial):
                continue
            if existing.size:
                d2 = np.sum(_minimum_image_delta(existing - trial, box) ** 2, axis=1)
                if np.min(d2) < cutoff2:
                    continue
            if placed:
                placed_xyz = np.array([[a["x"], a["y"], a["z"]] for a in placed], dtype=float)
                d2_placed = np.sum(_minimum_image_delta(placed_xyz - trial, box) ** 2, axis=1)
                if np.min(d2_placed) < cutoff2:
                    continue
            placed.append(
                {
                    "record": "HETATM",
                    "serial": 0,
                    "name": name,
                    "resname": resname,
                    "chain": "I",
                    "resseq": len(placed) + 1,
                    "icode": "",
                    "x": float(trial[0]),
                    "y": float(trial[1]),
                    "z": float(trial[2]),
                    "occ": 1.0,
                    "bfac": 0.0,
                    "segid": segid,
                    "element": name[:2],
                    "charge": "",
                }
            )
            accepted = True
            break
        if not accepted:
            raise RuntimeError("Failed to place ions without overlaps; relax cutoff or enlarge box.")
    return placed


def build_backbone_with_virtual_bb(backbone_atoms, bb_type_by_residue, bb_secondary_by_residue):
    """Append one BB proxy per residue after the N/CA/C/O block, placed at the 14/12/12/16 mass centre."""
    residue_groups = defaultdict(dict)
    for atom_idx, atom in enumerate(backbone_atoms):
        residue_groups[residue_key(atom)][atom["name"].strip().upper()] = atom_idx

    masses = np.asarray(BB_COMPONENT_MASSES, dtype=np.float32)
    weights = (masses / masses.sum()).astype(np.float64)
    protein_atoms = [deepcopy(atom) for atom in backbone_atoms]
    bb_entries = []
    for seq_idx, key in enumerate(sorted(residue_groups.keys()), start=1):
        role_map = residue_groups[key]
        if not all(role in role_map for role in BB_COMPONENT_NAMES):
            continue
        idxs = [int(role_map[role]) for role in BB_COMPONENT_NAMES]
        com = weights @ coords([backbone_atoms[i] for i in idxs])
        proxy = deepcopy(backbone_atoms[role_map["CA"]])
        proxy.update(name="BB", x=float(com[0]), y=float(com[1]), z=float(com[2]), segid="PROA", element="C")
        bb_entries.append(
            {
                "bb_residue_index": seq_idx,
                "bb_resseq": int(key[1]),
                "bb_chain": str(key[0]),
                "bb_atom_index": len(protein_atoms),
                "bb_type": str(bb_type_by_residue[key]),
                "bb_secondary_structure": str(bb_secondary_by_residue[key]),
                "atom_indices": idxs,
                "atom_mask": [1, 1, 1, 1],
                "weights": weights.tolist(),
            }
        )
        protein_atoms.append(proxy)
    if not bb_entries:
        raise ValueError("No complete backbone mapping entries could be generated.")
    return protein_atoms, bb_entries


def derive_chain_break_metadata(bb_entries):
    if not bb_entries:
        return np.array([], dtype=np.int32), np.array([], dtype=np.int32)

    # Each chain is one segment: a numbering gap inside a chain is refused, since its two ends are no termini and
    # the missing residues belong in the model
    chain_first_residue = []
    chain_counts = [1]
    prev_chain = str(bb_entries[0].get("bb_chain", "")).strip() or " "
    prev_resseq = int(bb_entries[0].get("bb_resseq", 0))
    for residue_index, entry in enumerate(bb_entries[1:], start=1):
        chain = str(entry.get("bb_chain", "")).strip() or " "
        resseq = int(entry.get("bb_resseq", residue_index))
        if chain != prev_chain:
            chain_first_residue.append(residue_index)
            chain_counts.append(1)
        elif resseq != prev_resseq + 1:
            raise ValueError(f"chain {chain} skips from residue {prev_resseq} to {resseq}; model the missing "
                             "residues before preparing a hybrid")
        prev_chain = chain
        prev_resseq = resseq
    return np.asarray(chain_first_residue, dtype=np.int32), np.asarray(chain_counts, dtype=np.int32)


def write_hybrid_mapping_h5(
    path: Path,
    bb_entries,
    total_martini_atoms,
    env_atom_indices,
    n_protein_atoms,
):
    def strings(values):
        return np.array([str(v).strip() or " " for v in values], dtype=h5py.string_dtype(encoding="utf-8"))

    with h5py.File(path, "w") as h5:
        inp = h5.create_group("input")

        ctrl = inp.create_group("hybrid_control")
        ctrl.attrs["activation_stage"] = b"minimization"

        bb_grp = inp.create_group("hybrid_bb_map")
        bb_grp.create_dataset("bb_residue_index", data=np.array([b["bb_residue_index"] for b in bb_entries], dtype=np.int32))
        bb_grp.create_dataset("bb_resseq", data=np.array([b["bb_resseq"] for b in bb_entries], dtype=np.int32))
        bb_grp.create_dataset("bb_chain_id", data=strings(b["bb_chain"] for b in bb_entries))
        bb_grp.create_dataset("bb_atom_index", data=np.array([b["bb_atom_index"] for b in bb_entries], dtype=np.int32))
        bb_grp.create_dataset("bb_type", data=strings(b["bb_type"] for b in bb_entries))
        bb_grp.create_dataset("bb_secondary_structure", data=strings(b["bb_secondary_structure"] for b in bb_entries))
        bb_grp.create_dataset("atom_indices", data=np.array([b["atom_indices"] for b in bb_entries], dtype=np.int32))
        bb_grp.create_dataset("atom_mask", data=np.array([b["atom_mask"] for b in bb_entries], dtype=np.int8))
        bb_grp.create_dataset("weights", data=np.array([b["weights"] for b in bb_entries], dtype=np.float32))
        bb_grp.create_dataset("reference_atom_names", data=strings(BB_COMPONENT_NAMES))

        env_grp = inp.create_group("hybrid_env_topology")
        env_grp.create_dataset("env_atom_indices", data=np.array(env_atom_indices, dtype=np.int32))
        membership = np.full(total_martini_atoms, -1, dtype=np.int32)
        membership[:n_protein_atoms] = 0
        env_grp.create_dataset("protein_membership", data=membership)

        chain_first_residue, chain_counts = derive_chain_break_metadata(bb_entries)
        if chain_first_residue.size:
            break_grp = inp.create_group("chain_break")
            break_grp.create_dataset("chain_first_residue", data=chain_first_residue, dtype=np.int32)
            break_grp.create_dataset("chain_counts", data=chain_counts, dtype=np.int32)


def read_martini3_nonbond_params(itp_file):
    _, pair_params = parse_dry_forcefield(itp_file)
    martini_table = {
        key: (value["sigma_nm"], value["epsilon_kj_mol"])
        for key, value in pair_params.items()
    }
    print(f"Read {len(martini_table)//2} unique nonbonded parameter pairs")
    return martini_table


def read_martini_masses(ff_file):
    ff_file_path = Path(ff_file).expanduser()
    if not ff_file_path.is_absolute():
        ff_file_path = (WORKFLOW_DIR / ff_file_path).resolve()
    else:
        ff_file_path = ff_file_path.resolve()
    
    if not ff_file_path.exists():
        raise ValueError(
            f"Force-field file '{ff_file}' not found; expected {ff_file_path}"
        )
    
    return parse_itp_atomtype_masses(ff_file_path)


def build_hybrid_mass_array_up(atom_types, particle_class, martini_masses):
    mass = np.empty(len(atom_types), dtype=np.float32)
    for index, (atom_type, class_value) in enumerate(
        zip(atom_types, particle_class)
    ):
        class_name = (
            class_value.decode("ascii")
            if isinstance(class_value, (bytes, np.bytes_))
            else str(class_value)
        ).strip().upper()
        if class_name == "PROTEIN":
            mass[index] = 1.0
            continue
        if atom_type not in martini_masses:
            raise ValueError(
                f"Mass not found for atom type '{atom_type}' (atom index {index}); "
                f"available atom types: {sorted(martini_masses.keys())}"
            )
        mass[index] = martini_masses[atom_type] / 12.0
    return mass


PROTEIN_RESIDUE_NAMES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX',
}


def _load_martini_topology(ff_path, stage_lipidhead_fc):
    """Read the dry MARTINI force field, the lipid and ion topologies and the bead masses."""
    def ff_file(name):
        path = os.path.join(ff_path, name)
        if not os.path.exists(path):
            raise ValueError(f"Required dry MARTINI file '{name}' not found in '{ff_path}'")
        return path

    martini_param_file = ff_file("dry_martini_v2.1.itp")
    martini_table = read_martini3_nonbond_params(martini_param_file)
    if not martini_table:
        raise ValueError(f"Could not read MARTINI parameters from '{martini_param_file}'")

    lipid_preproc_defs = {}
    if stage_lipidhead_fc > 0.0:
        lipid_preproc_defs['BILAYER_LIPIDHEAD_FC'] = stage_lipidhead_fc

    # Resolve each requested lipid moleculetype (one lipid or a comma-separated mix) against every ITP in
    # the force-field directory; bead types, charges and connectivity are kept per residue name.
    itp_paths = [os.path.join(ff_path, f) for f in sorted(os.listdir(ff_path)) if f.lower().endswith(".itp")]
    parsed_itps = [(path, parse_itp_file(path, preprocessor_defines=lipid_preproc_defs)) for path in itp_paths]

    lipids = {}
    for lipid_name in active_lipid_names():
        matches = [(path, mol) for path, topology in parsed_itps
                   for mol in topology.get('molecules', {}) if mol.upper() == lipid_name.upper()]
        if not matches:
            raise ValueError(f"Lipid molecule '{lipid_name}' not found in any ITP under '{ff_path}'")
        lipid_param_file, lipid_molecule = matches[0]
        lipid_topology = parse_itp_file(lipid_param_file, lipid_molecule, preprocessor_defines=lipid_preproc_defs)
        if not lipid_topology['bonds'] or not lipid_topology['angles']:
            raise ValueError(f"{lipid_name} in '{lipid_param_file}' has no bonds or no angles")
        lipids[lipid_molecule.upper()] = {
            'atom_names': [atom['atom'].upper() for atom in lipid_topology['atoms']],
            'atom_to_type': {atom['atom'].upper(): atom['type'] for atom in lipid_topology['atoms']},
            'atom_to_charge': {atom['atom'].upper(): atom['charge'] for atom in lipid_topology['atoms']},
            'bonds': [(b['i'], b['j'], b['r0'], b['k']) for b in lipid_topology['bonds']],
            'angles': [(a['i'], a['j'], a['k'], a['theta0'], a['force_k']) for a in lipid_topology['angles']],
            'position_restraints': lipid_topology.get('position_restraints', []),
        }
        print(f"Read {lipid_molecule} topology from {lipid_param_file}: {len(lipid_topology['atoms'])} beads, "
              f"{len(lipid_topology['bonds'])} bonds, {len(lipid_topology['angles'])} angles")

    ion_topology = parse_itp_file(ff_file("dry_martini_v2.1_ions.itp"))
    ions = {}
    for name in ("NA", "CL"):
        atom = next(a for a in ion_topology['atoms'] if a['atom'].upper() == name)
        ions[name] = (atom['type'], atom['charge'])

    return {
        'martini_table': martini_table,
        'martini_masses': read_martini_masses(martini_param_file),
        'lipids': lipids,
        'ions': ions,
    }


def _parse_runtime_pdb_atoms(input_pdb_file, topo):
    """Parse the runtime PDB, map each atom to its MARTINI bead type and charge, read the box."""
    lipids = topo['lipids']
    atoms, box = parse_pdb(Path(input_pdb_file))
    if box is None:
        raise ValueError(f"No CRYST1 record found in PDB file '{input_pdb_file}'")
    x_len, y_len, z_len = (float(v) for v in box)

    bb_type_by_residue = {}
    if any(atom["resname"].upper() in PROTEIN_RESIDUE_NAMES for atom in atoms):
        protein_backbone_atoms, _ = extract_protein_backbone_atoms_from_aa(atoms)
        bb_type_by_residue = map_backbone_types_from_structure(protein_backbone_atoms)[0]

    atom_types, charges = [], []
    for atom in atoms:
        atom_name = atom["name"].upper()
        residue_name = atom["resname"].upper()
        if residue_name in PROTEIN_RESIDUE_NAMES:
            if atom_name not in {"N", "CA", "C", "O", "BB"}:
                raise ValueError(f"Protein atom '{atom_name}' in residue '{residue_name}' is not one of N/CA/C/O/BB")
            # The residue's MARTINI identity is its BB bead; N/CA/C/O are the atoms BB is built from and
            # carry no MARTINI charge of their own.
            bb_type = bb_type_by_residue[(atom["chain"], atom["resseq"], atom["icode"], residue_name)]
            atom_types.append(bb_type)
            charges.append(BB_TYPE_CHARGE.get(bb_type, 0.0) if atom_name == "BB" else 0.0)
        elif lipid_resname(residue_name):
            lipid_rec = lipids[residue_name]
            if atom_name not in lipid_rec['atom_to_type']:
                raise ValueError(f"Unknown {residue_name} atom '{atom_name}'; the ITP has "
                                 f"{lipid_rec['atom_names']}")
            atom_types.append(lipid_rec['atom_to_type'][atom_name])
            charges.append(lipid_rec['atom_to_charge'][atom_name])
        elif residue_name in topo['ions']:
            bead_type, charge = topo['ions'][residue_name]
            atom_types.append(bead_type)
            charges.append(charge)
        else:
            raise ValueError(f"Unknown residue type '{residue_name}' for atom '{atom_name}'. Supported residue "
                             f"types: PROTEIN, {active_lipid_name()}, NA, CL")

    print(f"Runtime PDB {input_pdb_file}: {len(atoms)} atoms, box {x_len:.3f} x {y_len:.3f} x {z_len:.3f} A")
    return (np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float),
            np.array(atom_types),
            np.array(charges, dtype=float),
            np.array([a["resseq"] for a in atoms], dtype=int),
            np.array([a["name"].upper() for a in atoms]),
            np.array([a["resname"].upper() for a in atoms]),
            np.array([a["chain"] for a in atoms]),
            np.array([a["segid"] for a in atoms]),
            len(atoms), x_len, y_len, z_len)


def _group_atoms_into_molecules(residue_ids, residue_names, atom_names, chain_ids, seg_ids, n_atoms):
    """Group atoms into molecules: one per protein chain, one per residue for everything else."""
    molecules = []
    current = None
    protein_chains = {}
    for i, (resid, resname, atom_name) in enumerate(zip(residue_ids, residue_names, atom_names)):
        if resname in PROTEIN_RESIDUE_NAMES:
            # PDB chain id first: all protein atoms share the segid PROA. The BB proxies follow every chain's
            # N/CA/C/O, so a chain is one molecule wherever its atoms stand.
            key = ('PROTEIN', chain_ids[i] or seg_ids[i])
            if key not in protein_chains:
                protein_chains[key] = (key, [], [])
                molecules.append(protein_chains[key])
            current = protein_chains[key]
        else:
            key = (resname, resid)
            if current is None or key != current[0]:
                current = (key, [], [])
                molecules.append(current)
        current[1].append(atom_name)
        current[2].append(i)
    molecules = [(key[0], names, indices) for key, names, indices in molecules]

    molecule_ids = np.zeros(n_atoms, dtype=np.int32)
    for mol_idx, (_, _, atom_indices) in enumerate(molecules):
        molecule_ids[atom_indices] = mol_idx

    protein_sequence = []
    seen = set()
    for idx, (resid, resname) in enumerate(zip(residue_ids, residue_names)):
        key = (chain_ids[idx] or seg_ids[idx], int(resid))
        if resname in PROTEIN_RESIDUE_NAMES and key not in seen:
            seen.add(key)
            protein_sequence.append(normalize_resname(str(resname)))

    mol_counts = Counter(mol_type for mol_type, _, _ in molecules)
    print("Molecule summary: " + ", ".join(f"{t} {n}" for t, n in mol_counts.items()))
    lipid_count = sum(mol_counts.get(name, 0) for name in active_lipid_names())
    return molecules, molecule_ids, protein_sequence, lipid_count


def _build_bonded_terms(molecules, topo, reference_positions, bond_conversion, angle_conversion):
    """Lipid bonds, angles and head-group position restraints.

    Each ITP term names its beads by ITP index; the index is resolved to the bead's NAME and the name to the
    molecule's PDB atom, so the terms do not depend on the order the coordinates list the beads in. A
    molecule whose bead names are not exactly the ITP's is refused. Restraint references are taken from
    `reference_positions`, the coordinates the stage file stores, so reference and position share a frame.
    """
    lipids = topo['lipids']
    bonds, bond_lengths, bond_k = [], [], []
    angles, angle_equil_deg, angle_k = [], [], []
    restraint_indices, restraint_ref_pos, restraint_spring_xyz = [], [], []

    active_lipids = set(active_lipid_names())
    lipid_molecules = [mol for mol in molecules if mol[0] in active_lipids]
    for mol_type, atom_names_mol, atom_indices in lipid_molecules:
        rec = lipids[mol_type]
        itp_names = rec['atom_names']
        if sorted(atom_names_mol) != sorted(itp_names):
            raise ValueError(f"{mol_type} molecule at atoms {atom_indices[0]}-{atom_indices[-1]} has beads "
                             f"{list(atom_names_mol)}; its ITP defines {itp_names}")
        name_to_idx = dict(zip(atom_names_mol, atom_indices))

        def atom(itp_index):
            if not 0 <= itp_index < len(itp_names):
                raise ValueError(f"{mol_type} ITP term refers to bead {itp_index + 1} of {len(itp_names)}")
            return name_to_idx[itp_names[itp_index]]

        for i, j, r0_nm, k_kj in rec['bonds']:
            bonds.append([atom(i), atom(j)])
            bond_lengths.append(r0_nm * 10.0)
            bond_k.append(k_kj * bond_conversion)
        for i, j, k, theta0_deg, k_kj in rec['angles']:
            angles.append([atom(i), atom(j), atom(k)])
            angle_equil_deg.append(theta0_deg)
            angle_k.append(k_kj * angle_conversion)
        for restraint in rec['position_restraints']:
            atom_idx = atom(restraint['i'])
            restraint_indices.append(atom_idx)
            restraint_ref_pos.append(reference_positions[atom_idx].tolist())
            restraint_spring_xyz.append([restraint['fx'] * bond_conversion,
                                         restraint['fy'] * bond_conversion,
                                         restraint['fz'] * bond_conversion])

    print(f"Created {len(bonds)} {active_lipid_name()} bonds, {len(angles)} angles and "
          f"{len(restraint_indices)} head-group restraints across {len(lipid_molecules)} lipids")
    return (bonds, bond_lengths, bond_k, angles, angle_equil_deg, angle_k,
            restraint_indices, restraint_ref_pos, restraint_spring_xyz)


def _write_nonbonded_pairs(t, martini_potential, n_atoms, atom_types, charges, bonds_list, martini_table,
                           energy_conversion, length_conversion, martini_site_mask):
    """Write the optimized non-bonded pair table (pairs + coefficient indices).

    Only the sites in martini_site_mask carry MARTINI non-bonded interactions. On the protein that is the
    BB bead alone: dry-MARTINI represents a residue backbone as one bead, BB is built as the mass-weighted
    centre of N/CA/C/O and carries the secondary-structure-dependent BB type, so giving N, CA, C and O
    their own full-strength bead interactions as well would count the backbone against the environment
    five times over. Intra-protein physics belongs to the Upside force field and is excluded separately.
    """
    # Create atom indices and charges arrays for the potential
    t.create_array(martini_potential, 'atom_indices', obj=np.arange(n_atoms))
    t.create_array(martini_potential, 'charges', obj=charges)

    # Create pairs and optimized coefficient indices for non-bonded interactions.
    # The full coefficient table is O(N^2) and can OOM for large membrane boxes;
    # store unique coefficient rows plus one int index per pair instead.
    # MARTINI uses nrexcl=1: only 1-2 bonded pairs are excluded. The dry-MARTINI lipid ITPs carry no
    # [ exclusions ] section.
    bonded_pairs_12 = set((min(b[0], b[1]), max(b[0], b[1])) for b in bonds_list)

    # Generate all unique pairs (i < j) over the interacting sites, with nrexcl=1 exclusions
    active_sites = np.where(np.asarray(martini_site_mask, dtype=bool))[0].astype(np.int32)
    excluded_12_count = 0
    n_active = active_sites.size
    total_candidate_pairs = n_active * (n_active - 1) // 2
    total_pairs_written = 0
    unique_coeffs = []
    coeff_to_index = {}
    atom_type_strings = [
        atom_type.decode('utf-8') if isinstance(atom_type, bytes) else str(atom_type)
        for atom_type in atom_types
    ]
    atom_classes = [(atom_type_strings[i], float(charges[i])) for i in range(n_atoms)]
    class_to_id = {}
    atom_class_ids = np.empty(n_atoms, dtype=np.int32)
    class_values = []
    for i, cls in enumerate(atom_classes):
        class_id = class_to_id.get(cls)
        if class_id is None:
            class_id = len(class_values)
            class_to_id[cls] = class_id
            class_values.append(cls)
        atom_class_ids[i] = class_id

    def coefficient_index_for_classes(class_i, class_j, atom_i):
        type1, q1_raw = class_values[class_i]
        type2, q2_raw = class_values[class_j]
        if (type1, type2) in martini_table:
            sigma_nm, epsilon_kj = martini_table[(type1, type2)]
        elif (type2, type1) in martini_table:
            sigma_nm, epsilon_kj = martini_table[(type2, type1)]
        else:
            available_types = sorted(set([t[0] for t in martini_table.keys()] + [t[1] for t in martini_table.keys()]))
            raise ValueError(
                f"Missing interaction parameters for bead type pair ({type1}, {type2}) "
                f"at atom index {atom_i}; available bead types: {available_types}"
            )

        epsilon = np.float32(epsilon_kj / energy_conversion)
        sigma = np.float32(sigma_nm * length_conversion)
        q1 = np.float32(q1_raw)
        q2 = np.float32(q2_raw)
        coeff_key = (float(epsilon), float(sigma), float(q1), float(q2))
        coeff_index = coeff_to_index.get(coeff_key)
        if coeff_index is None:
            coeff_index = len(unique_coeffs)
            coeff_to_index[coeff_key] = coeff_index
            unique_coeffs.append((float(epsilon), float(sigma), float(q1), float(q2)))
        return coeff_index

    pairs_array = t.create_earray(
        martini_potential,
        'pairs',
        atom=tb.Int32Atom(),
        shape=(0, 2),
        expectedrows=total_candidate_pairs,
    )
    coeff_index_array = t.create_earray(
        martini_potential,
        'coefficient_indices',
        atom=tb.Int64Atom(),
        shape=(0,),
        expectedrows=total_candidate_pairs,
    )

    bonded_by_i = defaultdict(set)
    for i, j in bonded_pairs_12:
        bonded_by_i[i].add(j)

    for position in range(n_active):
        i = int(active_sites[position])
        js = active_sites[position + 1:]
        if js.size == 0:
            continue
        bonded_js = bonded_by_i.get(i)
        if bonded_js:
            bonded_mask = np.isin(js, np.fromiter(bonded_js, dtype=np.int32), assume_unique=True)
            excluded_12_count += int(np.count_nonzero(bonded_mask))
            js = js[~bonded_mask]
        if js.size == 0:
            continue

        pairs_chunk = np.empty((js.size, 2), dtype=np.int32)
        pairs_chunk[:, 0] = i
        pairs_chunk[:, 1] = js

        coeff_indices = np.empty(js.size, dtype=np.int64)
        class_i = int(atom_class_ids[i])
        js_class_ids = atom_class_ids[js]
        for class_j in np.unique(js_class_ids):
            class_mask = js_class_ids == class_j
            coeff_index = coefficient_index_for_classes(class_i, int(class_j), i)
            coeff_indices[class_mask] = coeff_index

        pairs_array.append(pairs_chunk)
        coeff_index_array.append(coeff_indices)
        total_pairs_written += int(js.size)

    print(f"MARTINI non-bonded sites: {n_active} of {n_atoms} "
          f"({n_atoms - n_active} excluded; protein carries BB only)")
    print(f"Excluded {excluded_12_count} 1-2 bonded pairs from non-bonded interactions (nrexcl=1)")
    print(f"Wrote {total_pairs_written} non-bonded pairs with {len(unique_coeffs)} unique coefficient rows")
    martini_potential._v_attrs.optimized_format = 1
    coeff_array = np.array(unique_coeffs, dtype='f4')
    if coeff_array.size == 0:
        coeff_array = np.zeros((0, 4), dtype='f4')
    t.create_array(martini_potential, 'coefficients', obj=coeff_array)

    return excluded_12_count


def convert_stage(pdb_id, stage, run_dir, runtime_pdb):
    """Build the HDF5 input for one MARTINI workflow stage from the runtime PDB.

    `stage` selects the pair tables: "soft_core" for CHARMM-GUI's soft-core minimisation, "hard" for the
    published potential. The force-field directory, unit conversions, barostat, head-group restraint
    constant and Brownian clock come from the UPSIDE_* environment the workflow sets for the stage.
    """
    if stage not in STAGE_SOFT_CORE_LAMBDA:
        valid = ", ".join(sorted(STAGE_SOFT_CORE_LAMBDA))
        raise ValueError(f"Unknown MARTINI preparation stage {stage!r}; expected one of: {valid}")
    soft_core_lambda = STAGE_SOFT_CORE_LAMBDA[stage]
    stage_lipidhead_fc = float(os.environ.get('UPSIDE_BILAYER_LIPIDHEAD_FC', '0'))
    os.makedirs(run_dir, exist_ok=True)
    print(f"Preparing MARTINI stage {stage} for {pdb_id} in {run_dir}")

    ff_path = str(Path(os.environ['UPSIDE_MARTINI_FF_DIR']).expanduser().resolve())
    topo = _load_martini_topology(ff_path, stage_lipidhead_fc)
    martini_masses = topo['martini_masses']

    energy_conversion = float(os.environ['UPSIDE_MARTINI_ENERGY_CONVERSION'])
    length_conversion = float(os.environ['UPSIDE_MARTINI_LENGTH_CONVERSION'])
    if energy_conversion <= 0.0 or length_conversion <= 0.0:
        raise ValueError("UPSIDE_MARTINI_ENERGY_CONVERSION and UPSIDE_MARTINI_LENGTH_CONVERSION must be positive")
    bond_conversion = 1.0 / (energy_conversion * length_conversion ** 2)
    angle_conversion = 1.0 / energy_conversion
    print(f"Unit conversions: energy={energy_conversion:g} kJ/mol per E_up, length={length_conversion:g} A/nm")

    (initial_positions, atom_types, charges, residue_ids, atom_names,
     residue_names, chain_ids, seg_ids, n_atoms, x_len, y_len, z_len) = \
        _parse_runtime_pdb_atoms(runtime_pdb, topo)
    net_charge = float(charges.sum())
    if abs(net_charge) > 1e-6:
        raise ValueError(f"{runtime_pdb} carries a net charge of {net_charge:+.6g} e; the ions must neutralise "
                         "the charges the topology assigns")

    molecules, molecule_ids, protein_sequence, lipid_count = _group_atoms_into_molecules(
        residue_ids, residue_names, atom_names, chain_ids, seg_ids, n_atoms
    )

    # Centre on the mean of all atoms and wrap each molecule whole, by its centre, into [-L/2, L/2): the protein's
    # Upside terms are not periodic, so a protein atom moved by a box length would be a bond one box long. Every
    # array that stores coordinates (positions and restraint references alike) uses this frame.
    box = np.array([x_len, y_len, z_len])
    centered_positions = initial_positions - np.mean(initial_positions, axis=0)
    final_positions = centered_positions.copy()
    for _, _, indices in molecules:
        centre = centered_positions[indices].mean(axis=0)
        final_positions[indices] -= box * np.floor((centre + 0.5 * box) / box)

    (bonds_list, bond_lengths_list, bond_force_constants_list,
     angles_list, angle_equil_deg_list, angle_force_constants_list,
     lipid_restraint_indices, lipid_restraint_ref_pos, lipid_restraint_spring_xyz) = _build_bonded_terms(
        molecules, topo, final_positions, bond_conversion, angle_conversion,
    )

    particle_class = np.empty(n_atoms, dtype='S10')
    for i, resname in enumerate(residue_names):
        if resname in PROTEIN_RESIDUE_NAMES:
            particle_class[i] = b"PROTEIN"
        elif lipid_resname(resname):
            particle_class[i] = b"LIPID"
        else:
            particle_class[i] = b"ION"

    input_file = f"{run_dir}/test.input.up"
    with tb.open_file(input_file, 'w') as t:
        input_grp = t.create_group(t.root, 'input')

        pos = np.zeros((n_atoms, 3, 1), dtype='f4')
        pos[:,:,0] = final_positions
        pos_array = t.create_array(input_grp, 'pos', obj=pos)
        pos_array._v_attrs.arguments = np.array([b'pos'])
        pos_array._v_attrs.shape = pos.shape
        pos_array._v_attrs.n_atoms = n_atoms
        pos_array._v_attrs.n_frames = 1
        pos_array._v_attrs.dim = 3
        pos_array._v_attrs.initialized = True

        momentum = np.zeros((n_atoms, 3, 1), dtype='f4')
        mom_array = t.create_array(input_grp, 'mom', obj=momentum)
        mom_array._v_attrs.arguments = np.array([b'mom'])
        mom_array._v_attrs.shape = momentum.shape
        mom_array._v_attrs.n_atoms = n_atoms
        mom_array._v_attrs.dim = 3
        mom_array._v_attrs.initialized = True

        mass = build_hybrid_mass_array_up(
            atom_types,
            particle_class,
            martini_masses,
        )
        mass_array = t.create_array(input_grp, 'mass', obj=mass)
        mass_array._v_attrs.arguments = np.array([b'mass'])
        mass_array._v_attrs.shape = mass.shape
        mass_array._v_attrs.n_atoms = n_atoms
        mass_array._v_attrs.initialized = True
        mass_array._v_attrs.protein_mass_up = np.float32(1.0)
        mass_array._v_attrs.protein_mass_source = b"upside_core_unit_mass"
        mass_array._v_attrs.environment_mass_source = b"native_dry_martini_mass_div_12"

        stage_grp = t.create_group(input_grp, 'stage_parameters')
        stage_grp._v_attrs.current_stage = b'minimization'

        barostat_enable = int(os.environ.get('UPSIDE_NPT_ENABLE', '0'))
        if barostat_enable:
            # Semi-isotropic Monte-Carlo barostat for the waterless (implicit-solvent) membrane:
            # xy is coupled to zero lateral tension; z is frozen (mc_dmax_z=0) because there is no
            # solvent to set the box normal. COM-based moves preserve the rigid protein.
            barostat_grp = t.create_group(input_grp, 'barostat')
            barostat_grp._v_attrs.target_p_xy = float(os.environ['UPSIDE_NPT_TARGET_PXY'])
            barostat_grp._v_attrs.target_p_z = float(os.environ['UPSIDE_NPT_TARGET_PZ'])
            barostat_grp._v_attrs.interval = int(os.environ['UPSIDE_NPT_INTERVAL'])
            barostat_grp._v_attrs.mc_dmax_xy = float(os.environ['UPSIDE_NPT_MC_DMAX_XY'])
            barostat_grp._v_attrs.mc_dmax_z = float(os.environ['UPSIDE_NPT_MC_DMAX_Z'])
            print(f"Semi-isotropic Monte-Carlo barostat: Pxy {barostat_grp._v_attrs.target_p_xy}, "
                  f"dmax xy/z {barostat_grp._v_attrs.mc_dmax_xy}/{barostat_grp._v_attrs.mc_dmax_z}, "
                  f"interval {barostat_grp._v_attrs.interval}")
        else:
            print("NPT barostat disabled (NVT mode)")

        type_array = t.create_array(input_grp, 'type', obj=atom_types.astype('S4'))
        type_array._v_attrs.arguments = np.array([b'type'])
        type_array._v_attrs.shape = atom_types.shape
        type_array._v_attrs.n_atoms = n_atoms
        type_array._v_attrs.initialized = True
        type_array._v_attrs.description = b"Interaction matrix type names (e.g., C1, C2, Qa, Qd)"

        particle_class_array = t.create_array(input_grp, 'particle_class', obj=particle_class)
        particle_class_array._v_attrs.arguments = np.array([b'particle_class'])
        particle_class_array._v_attrs.shape = particle_class.shape
        particle_class_array._v_attrs.n_atoms = n_atoms
        particle_class_array._v_attrs.initialized = True
        particle_class_array._v_attrs.description = b"Per-atom class: PROTEIN, LIPID, ION"

        # Every physical degree of freedom uses one g-JF step. Regenerated O and BB sites are excluded:
        # their positions and BB-env derivatives are projected through N/CA/C. Production gives each lipid
        # bead the bare-particle friction of its own molecule's diffusion target on the requested clock, and a
        # protein carrier the sum of the frictions of the lipid beads inside the protein--environment spline
        # range. The preparation stages use the dry-MARTINI Langevin relaxation, friction = mass / 4 ps.
        lipid_beads = np.where(particle_class == b"LIPID")[0].astype('i4')
        if lipid_beads.size:
            protein_carrier = (
                (particle_class == b"PROTEIN")
                & np.isin(atom_names, np.array(["N", "CA", "C"], dtype=atom_names.dtype))
            )
            dynamic_atom = (particle_class != b"PROTEIN") | protein_carrier
            brownian_atom_index = np.where(dynamic_atom)[0].astype('i4')
            numerical_time_step = float(os.environ.get("UPSIDE_MARTINI_TIME_STEP_UP", "0.009"))
            protein_time_ps = float(os.environ.get("UPSIDE_PROTEIN_TIME_PS_PER_STEP", "40.0"))
            martini_time_factor = float(os.environ.get("UPSIDE_MARTINI_TIME_FACTOR", "4.0"))
            reference_temperature = float(os.environ.get("UPSIDE_DOPC_REFERENCE_TEMPERATURE_UP", "0.8647"))
            raw_relaxation_ps = float(os.environ.get("UPSIDE_DRY_MARTINI_RELAXATION_PS", "4.0"))
            dynamics_phase = os.environ.get("UPSIDE_MARTINI_DYNAMICS_PHASE", "production")
            lipid_names = list(active_lipid_names())
            lipid_targets = lipid_diffusion_targets(lipid_names)
            ion_target = LIPID_LATERAL_DIFFUSION_UM2_S[ION_CLOCK_REFERENCE_LIPID]
            # Substeps of numerical_time_step/inner_steps inside each outer step. The outer step is
            # what the friction and the ps-per-step clock are calibrated against and stays 0.009;
            # this only integrates the same Langevin equation more accurately. Required because
            # mass-1 protein carriers at dt = 0.009 pick up a kinetic-temperature excess against the
            # steep MARTINI core: measured on glpG-RKRK-79ALA at T = 0.7215, the backbone runs
            # +10.1% above its set point at inner_steps = 1, +3.5% at 2 and +1.1% at 4, while the
            # lipids go 1.024 -> 1.009, for 1.27x and 2.05x the wall time (findings 3.10).
            inner_steps = int(os.environ.get("UPSIDE_MARTINI_INNER_STEPS", "4"))
            if inner_steps < 1:
                raise ValueError("UPSIDE_MARTINI_INNER_STEPS must be >= 1")
            if numerical_time_step <= 0.0 or protein_time_ps <= 0.0 or martini_time_factor <= 0.0:
                raise ValueError("MARTINI clock parameters must be positive")
            if reference_temperature <= 0.0 or raw_relaxation_ps <= 0.0:
                raise ValueError("The friction reference temperature and relaxation time must be positive")
            if dynamics_phase not in {"equilibration", "production"}:
                raise ValueError(f"Unknown MARTINI dynamics phase '{dynamics_phase}'")
            raw_martini_time_ps = protein_time_ps / martini_time_factor
            time_unit_ps = 1.0e12 * math.sqrt(
                (0.012 * 1.0e-20) / (energy_conversion * 1000.0))

            def bead_diffusion_up(target_um2_s):
                """D_bead,up in A^2 per time unit: the factor-four raw MARTINI diffusion on the raw step clock."""
                return target_um2_s * martini_time_factor * 1.0e-4 * raw_martini_time_ps / numerical_time_step

            lipid_bead_diffusion = {n: bead_diffusion_up(d) for n, d in lipid_targets.items()}
            lipid_bead_friction = {n: reference_temperature / d for n, d in lipid_bead_diffusion.items()}
            ion_friction = reference_temperature / bead_diffusion_up(ion_target)

            interface_cutoff = 12.0
            contact_count = np.zeros(n_atoms, dtype=np.int32)
            protein_friction = protein_carrier[brownian_atom_index]
            if dynamics_phase == "production":
                friction = np.full(brownian_atom_index.size, ion_friction, dtype=np.float64)
                bead_resname = residue_names[brownian_atom_index]
                carrier_index = np.where(protein_carrier)[0]
                carrier_friction = np.zeros(n_atoms, dtype=np.float64)
                box = np.array([x_len, y_len, z_len], dtype=np.float64)
                for name in lipid_names:
                    friction[bead_resname == name] = lipid_bead_friction[name]
                    beads = lipid_beads[residue_names[lipid_beads] == name]
                    count = np.zeros(n_atoms, dtype=np.int32)
                    for start in range(0, beads.size, 512):
                        delta = (
                            final_positions[carrier_index, None, :]
                            - final_positions[beads[start:start + 512]][None, :, :]
                        )
                        delta -= box * np.rint(delta / box)
                        count[carrier_index] += np.sum(
                            np.sum(delta * delta, axis=2) < interface_cutoff * interface_cutoff,
                            axis=1,
                        )
                    contact_count += count
                    carrier_friction += lipid_bead_friction[name] * count
                friction[protein_friction] = carrier_friction[brownian_atom_index[protein_friction]]
            else:
                relaxation_up = raw_relaxation_ps / time_unit_ps
                friction = mass[brownian_atom_index] / relaxation_up
                friction[protein_friction] = 0.0
            brownian_grp = t.create_group(input_grp, 'brownian')
            t.create_array(brownian_grp, 'atom_index', obj=brownian_atom_index)
            t.create_array(brownian_grp, 'friction', obj=friction.astype(np.float32))
            brownian_grp._v_attrs.numerical_time_step = np.float64(numerical_time_step)
            brownian_grp._v_attrs.inner_steps = np.int32(inner_steps)
            brownian_grp._v_attrs.protein_time_ps_per_step = np.float64(protein_time_ps)
            brownian_grp._v_attrs.martini_time_factor = np.float64(martini_time_factor)
            brownian_grp._v_attrs.raw_martini_time_ps_per_step = np.float64(raw_martini_time_ps)
            brownian_grp._v_attrs.lipid_names = np.array([n.encode() for n in lipid_names])
            brownian_grp._v_attrs.lipid_target_diffusion_um2_s = np.array(
                [lipid_targets[n] for n in lipid_names], dtype=np.float64)
            brownian_grp._v_attrs.ion_target_diffusion_um2_s = np.float64(ion_target)
            brownian_grp._v_attrs.reference_temperature_up = np.float64(reference_temperature)
            brownian_grp._v_attrs.dynamics_phase = dynamics_phase.encode()
            if dynamics_phase == "production":
                contact_atom_index = np.where(contact_count > 0)[0].astype('i4')
                t.create_array(brownian_grp, 'protein_contact_atom_index', obj=contact_atom_index)
                t.create_array(
                    brownian_grp,
                    'protein_lipid_contact_count',
                    obj=contact_count[contact_atom_index],
                )
                brownian_grp._v_attrs.lipid_bead_diffusion_up = np.array(
                    [lipid_bead_diffusion[n] for n in lipid_names], dtype=np.float64)
                brownian_grp._v_attrs.lipid_bead_friction_up = np.array(
                    [lipid_bead_friction[n] for n in lipid_names], dtype=np.float64)
                brownian_grp._v_attrs.ion_bead_friction_up = np.float64(ion_friction)
                brownian_grp._v_attrs.interface_friction_cutoff_A = np.float64(interface_cutoff)
                brownian_grp._v_attrs.transport_observable = b"bare_martini_particle_lateral_diffusion"
            else:
                brownian_grp._v_attrs.raw_langevin_relaxation_ps = np.float64(raw_relaxation_ps)
                brownian_grp._v_attrs.relaxation_time_up = np.float64(relaxation_up)
                brownian_grp._v_attrs.time_unit_ps = np.float64(time_unit_ps)
                brownian_grp._v_attrs.transport_observable = b"native_dry_martini_relaxation"

        charge_array = t.create_array(input_grp, 'charges', obj=charges)
        charge_array._v_attrs.arguments = np.array([b'charges'])
        charge_array._v_attrs.shape = charges.shape
        charge_array._v_attrs.n_atoms = n_atoms
        charge_array._v_attrs.initialized = True

        residue_array = t.create_array(input_grp, 'residue_ids', obj=residue_ids)
        residue_array._v_attrs.arguments = np.array([b'residue_ids'])
        residue_array._v_attrs.shape = residue_ids.shape
        residue_array._v_attrs.n_atoms = n_atoms
        residue_array._v_attrs.initialized = True
        residue_array._v_attrs.description = b"PDB residue indices per atom"

        molecule_array = t.create_array(input_grp, 'molecule_ids', obj=molecule_ids)
        molecule_array._v_attrs.arguments = np.array([b'molecule_ids'])
        molecule_array._v_attrs.shape = molecule_ids.shape
        molecule_array._v_attrs.n_atoms = n_atoms
        molecule_array._v_attrs.initialized = True
        molecule_array._v_attrs.description = b"Per-atom molecule index (0-based, contiguous)"

        atom_name_array = t.create_array(input_grp, 'atom_names', obj=atom_names.astype('S4'))
        atom_name_array._v_attrs.arguments = np.array([b'atom_names'])
        atom_name_array._v_attrs.shape = atom_names.shape
        atom_name_array._v_attrs.n_atoms = n_atoms
        atom_name_array._v_attrs.initialized = True

        atom_roles = atom_names.astype('S4')
        atom_role_array = t.create_array(input_grp, 'atom_roles', obj=atom_roles)
        atom_role_array._v_attrs.arguments = np.array([b'atom_roles'])
        atom_role_array._v_attrs.shape = atom_roles.shape
        atom_role_array._v_attrs.n_atoms = n_atoms
        atom_role_array._v_attrs.initialized = True
        atom_role_array._v_attrs.description = b"PDB atom role names (BB, SC1, SC2, W, etc.)"

        residue_name_array = t.create_array(input_grp, 'residue_names', obj=residue_names.astype('S4'))
        residue_name_array._v_attrs.description = b"PDB residue name per atom; names each lipid bead's molecule"

        if protein_sequence:
            t.create_array(input_grp, 'sequence', obj=np.asarray([x.encode('ascii') for x in protein_sequence], dtype='S3'))

        potential_grp = t.create_group(input_grp, 'potential')

        # The BB proxy is the mass-weighted centre of N/CA/C/O, and O is Upside's own derived site, so the
        # node consumes infer_H_O rather than rebuilding O from a local frame. That makes the proxy a plain
        # linear combination of nodes the engine differentiates: the O share of any environment force on BB
        # is handed back to infer_H_O, which carries it to CA, C and the next residue's N exactly.
        hybrid_position = t.create_group(potential_grp, 'martini_hybrid_position')
        hybrid_position._v_attrs.arguments = np.array([b'pos', b'infer_H_O'])

        martini_potential = t.create_group(potential_grp, 'martini_potential')
        martini_potential._v_attrs.arguments = np.array([b'martini_hybrid_position'])
        martini_potential._v_attrs.cutoff = 12.0
        martini_potential._v_attrs.cache_buffer = 2.0   # Verlet-list skin; rebuild is O(N) cell list, so keep the active list lean
        martini_potential._v_attrs.x_len = x_len
        martini_potential._v_attrs.y_len = y_len
        martini_potential._v_attrs.z_len = z_len
        # inject_particles_table evaluates the published pair potential at the soft-core radius for this
        # lambda; 0 is the potential itself.
        martini_potential._v_attrs.soft_core_lambda = soft_core_lambda

        # The protein presents one MARTINI site per residue, the BB bead. N, CA, C and O are the atoms BB
        # is built from, so they carry no bead interaction of their own with the environment.
        martini_site_mask = (particle_class != b"PROTEIN") | np.isin(
            atom_names, np.array(["BB"], dtype=atom_names.dtype))
        excluded_12_count = _write_nonbonded_pairs(
            t, martini_potential, n_atoms, atom_types, charges,
            bonds_list, topo['martini_table'],
            energy_conversion, length_conversion, martini_site_mask,
        )

        bond_group = t.create_group(potential_grp, 'dist_spring')
        bond_group._v_attrs.arguments = np.array([b'pos'])
        bond_group._v_attrs.x_len = x_len
        bond_group._v_attrs.y_len = y_len
        bond_group._v_attrs.z_len = z_len
        t.create_array(bond_group, 'id', obj=np.array(bonds_list, dtype=int))
        t.create_array(bond_group, 'equil_dist', obj=np.array(bond_lengths_list, dtype='f4'))
        t.create_array(bond_group, 'spring_const', obj=np.array(bond_force_constants_list, dtype='f4'))

        angle_group = t.create_group(potential_grp, 'angle_spring')
        angle_group._v_attrs.arguments = np.array([b'pos'])
        angle_group._v_attrs.x_len = x_len
        angle_group._v_attrs.y_len = y_len
        angle_group._v_attrs.z_len = z_len
        t.create_array(angle_group, 'id', obj=np.array(angles_list, dtype=int))
        t.create_array(angle_group, 'equil_angle_deg', obj=np.array(angle_equil_deg_list, dtype='f4'))
        t.create_array(angle_group, 'spring_const', obj=np.array(angle_force_constants_list, dtype='f4'))

        # Head-group height restraint of CHARMM-GUI's 6.2-6.6 ramp (BILAYER_LIPIDHEAD_FC in the lipid ITP).
        # Its reference is the assembled membrane in the stage file's own frame, as gmx grompp -r
        # step5_charmm2gmx.pdb uses the assembly.
        if lipid_restraint_indices:
            restraint_group = t.create_group(potential_grp, 'restraint_position')
            restraint_group._v_attrs.arguments = np.array([b'pos'])
            t.create_array(
                restraint_group, 'restraint_indices', obj=np.array(lipid_restraint_indices, dtype='i4')
            )
            t.create_array(
                restraint_group, 'ref_pos', obj=np.array(lipid_restraint_ref_pos, dtype='f4')
            )
            t.create_array(restraint_group, 'spring_const_xyz', obj=np.array(lipid_restraint_spring_xyz, dtype='f4'))

    print(f"Created UPSIDE input file: {input_file}")
    with open(f"{run_dir}/preparation_summary.txt", 'w') as f:
        f.write("MARTINI preparation summary\n")
        f.write(f"PDB ID: {pdb_id}\n")
        f.write(f"Stage tables: {stage} (soft-core lambda {soft_core_lambda:g})\n")
        f.write(f"Total atoms: {n_atoms}\n")
        f.write(f"Box dimensions: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} A\n")
        f.write(f"Total bonds: {len(bonds_list)}\n")
        f.write(f"Total angles: {len(angles_list)}\n")
        f.write(f"{active_lipid_name()} molecules: {lipid_count}\n")
        f.write(f"1-2 exclusions (nrexcl=1): {excluded_12_count}\n")


def require_group(root, path):
    if path not in root:
        raise ValueError(f"Missing group: {path}")
    return root[path]


def require_dataset(group, name):
    if name not in group:
        raise ValueError(f"Missing dataset: {group.name}/{name}")
    return group[name]


def validate_hybrid_mapping(mapping_h5: Path, n_atom: int | None = None):
    path = Path(mapping_h5).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(path)

    with h5py.File(path, "r") as h5:
        inp = require_group(h5, "/input")
        ctrl = require_group(inp, "hybrid_control")
        bb = require_group(inp, "hybrid_bb_map")
        env = require_group(inp, "hybrid_env_topology")

        if "activation_stage" not in ctrl.attrs:
            raise ValueError("Missing hybrid_control attr: activation_stage")

        bb_atom_idx = require_dataset(bb, "bb_atom_index")[:]
        bb_atom_map = require_dataset(bb, "atom_indices")[:]
        bb_mask = require_dataset(bb, "atom_mask")[:]
        bb_w = require_dataset(bb, "weights")[:]

        if bb_atom_map.ndim != 2 or bb_atom_map.shape[1] != 4:
            raise ValueError("hybrid_bb_map/atom_indices must have shape (n_bb,4)")
        if bb_mask.shape != bb_atom_map.shape or bb_w.shape != bb_atom_map.shape:
            raise ValueError("hybrid_bb_map mask/weights shapes must match atom_indices")

        n_bb = bb_atom_map.shape[0]
        if bb_atom_idx.shape != (n_bb,):
            raise ValueError("hybrid_bb_map/bb_atom_index must have shape (n_bb,)")
        for i in range(n_bb):
            mask = bb_mask[i].astype(bool)
            if mask.any():
                wsum = float(bb_w[i][mask].sum())
                if abs(wsum - 1.0) > 1e-4:
                    raise ValueError(f"BB weights do not sum to 1 for row {i}: {wsum}")

        membership = require_dataset(env, "protein_membership")[:]
        if membership.ndim != 1:
            raise ValueError("hybrid_env_topology/protein_membership must be 1D")
        if n_atom is not None and membership.shape[0] != n_atom:
            raise ValueError(
                f"protein_membership length ({membership.shape[0]}) != expected n_atom ({n_atom})"
            )
        n_atom_runtime = int(membership.shape[0])

        if "chain_break" in inp:
            chain_break = require_group(inp, "chain_break")
            chain_first_residue = require_dataset(chain_break, "chain_first_residue")[:]
            chain_counts = require_dataset(chain_break, "chain_counts")[:]
            if chain_first_residue.ndim != 1:
                raise ValueError("chain_break/chain_first_residue must be 1D")
            if chain_counts.ndim != 1:
                raise ValueError("chain_break/chain_counts must be 1D")
            if chain_counts.shape[0] != chain_first_residue.shape[0] + 1:
                raise ValueError("chain_break/chain_counts length must equal number of chains")
            if chain_first_residue.size and np.any(chain_first_residue[:-1] >= chain_first_residue[1:]):
                raise ValueError("chain_break/chain_first_residue must be strictly increasing")
            if chain_first_residue.size and (
                np.any(chain_first_residue <= 0) or np.any(chain_first_residue >= n_bb)
            ):
                raise ValueError("chain_break/chain_first_residue entries must be within (0, n_bb)")
            if np.any(chain_counts <= 0):
                raise ValueError("chain_break/chain_counts entries must be positive")

        for i in range(n_bb):
            bb_i = int(bb_atom_idx[i])
            if bb_i >= 0:
                if bb_i >= n_atom_runtime:
                    raise ValueError(f"BB proxy index out of bounds at row {i}: {bb_i}")
                if membership[bb_i] < 0:
                    raise ValueError(f"BB proxy index is not protein atom at row {i}: {bb_i}")
            elif bb_i != -1:
                raise ValueError(f"BB proxy index must be >=0 or -1 sentinel at row {i}: {bb_i}")

            if int(bb_mask[i, 1]) == 0:
                raise ValueError(f"BB map row {i} must include a CA carrier (atom_mask[:,1]==1)")
            for j in range(4):
                if int(bb_mask[i, j]) == 0:
                    continue
                ai = int(bb_atom_map[i, j])
                if ai < 0 or ai >= n_atom_runtime:
                    raise ValueError(f"BB target index out of bounds at row {i}, col {j}: {ai}")
                if membership[ai] < 0:
                    raise ValueError(f"BB target index is not protein atom at row {i}, col {j}: {ai}")

        n_protein = int(np.sum(membership >= 0))
        n_env = int(np.sum(membership < 0))
        print(f"OK: {path}")
        print(f"  n_atom={membership.shape[0]} n_protein={n_protein} n_env={n_env}")
        print(f"  n_bb={n_bb}")


def set_initial_position(input_file, output_file):
    """Hand the final state of `input_file` to `output_file`.

    Positions, the box and the momentum, with the round it was recorded at, come from the source's last
    output frame. A source that ran no MD (a minimisation) has no output momentum, so its own input momentum
    and round pass through unchanged. That is the whole restart state, so a production run in segments is the
    continuous run.
    """
    with h5py.File(input_file, "r") as f:
        has_output = "/output/pos" in f and f["/output/pos"].shape[0] > 0
        last_pos = (f["/output/pos"][-1, 0, :, :] if has_output else f["/input/pos"][:, :, 0])[:, :, np.newaxis]
        # the box the source ended in: its barostat's last box, or the fixed box it ran at
        if "/output/box" in f and f["/output/box"].shape[0] > 0:
            last_box = np.asarray(f["/output/box"][-1], dtype=float).reshape(-1)[:3]
        else:
            attrs = f["/input/potential/martini_potential"].attrs
            last_box = np.array([attrs["x_len"], attrs["y_len"], attrs["z_len"]], dtype=float)
        if "/output/mom" in f and f["/output/mom"].shape[0] > 0:
            last_mom = f["/output/mom"][-1, 0, :, :][:, :, np.newaxis]
            last_round = int(f["/output/round"][-1])
            restart_valid = 1
        else:
            last_mom = f["/input/mom"][:]
            last_round = int(f["/input/mom"].attrs.get("round", 0))
            restart_valid = int(f["/input/mom"].attrs.get("restart_valid", 0))

    with h5py.File(output_file, "r+") as f:
        if f["/input/pos"].shape[0] != last_pos.shape[0]:
            raise ValueError(f"{input_file} has {last_pos.shape[0]} atoms but {output_file} has "
                             f"{f['/input/pos'].shape[0]}; stages must share one topology")
        del f["/input/pos"]
        f.create_dataset("/input/pos", data=last_pos)
        del f["/input/mom"]
        mom_ds = f.create_dataset("/input/mom", data=last_mom)
        mom_ds.attrs["restart_valid"] = np.int8(restart_valid)
        mom_ds.attrs["round"] = np.int64(last_round)
        set_potential_box(f, last_box)


def set_potential_box(h5, box):
    """Write the box into every potential node that carries one."""
    n_updated = 0
    for pot_grp in h5["/input/potential"].values():
        if isinstance(pot_grp, h5py.Group) and all(k in pot_grp.attrs for k in ("x_len", "y_len", "z_len")):
            pot_grp.attrs["x_len"], pot_grp.attrs["y_len"], pot_grp.attrs["z_len"] = (float(v) for v in box)
            n_updated += 1
    print(f"Updated box dimensions for {n_updated} potential nodes: "
          f"x={box[0]:.3f}, y={box[1]:.3f}, z={box[2]:.3f}")


def normalize_resname(name: str) -> str:
    name = name.upper()
    aliases = {
        "HSD": "HIS",
        "HSE": "HIS",
        "HSP": "HIS",
        "HID": "HIS",
        "HIE": "HIS",
        "HIP": "HIS",
        "CYX": "CYS",
    }
    return aliases.get(name, name)


def decode_string_array(dataset):
    return [x.decode("ascii") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in dataset[:]]


def unique_preserving_order(values):
    out = []
    seen = set()
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def resolve_sequence(inp, residue_count: int):
    sequence = [normalize_resname(x) for x in decode_string_array(inp["sequence"])]
    if len(sequence) != residue_count:
        raise ValueError(f"/input/sequence has {len(sequence)} residues but hybrid_bb_map has {residue_count}")
    return sequence


def build_affine_atoms(inp):
    if "hybrid_bb_map" not in inp:
        raise ValueError("Missing /input/hybrid_bb_map in stage file")
    bb_grp = inp["hybrid_bb_map"]
    if "atom_indices" not in bb_grp:
        raise ValueError("hybrid_bb_map/atom_indices is required for stage-7 CB placement")

    bb_residue_raw = bb_grp["bb_residue_index"][:].astype(np.int32)
    bb_atom_idx = bb_grp["atom_indices"][:].astype(np.int32)
    residue_ids = unique_preserving_order(int(x) for x in bb_residue_raw.tolist())

    residue_to_ncac = {}
    n_atom = int(inp["pos"].shape[0])
    for resid, atom_row in zip(bb_residue_raw.tolist(), bb_atom_idx.tolist()):
        rid = int(resid)
        n_idx, ca_idx, c_idx = [int(atom_row[0]), int(atom_row[1]), int(atom_row[2])]
        if n_idx < 0 or ca_idx < 0 or c_idx < 0:
            continue
        for idx in (n_idx, ca_idx, c_idx):
            if idx < 0 or idx >= n_atom:
                raise ValueError(
                    f"Backbone carrier index out of bounds for residue {rid}: idx={idx}, n_atom={n_atom}"
                )
        residue_to_ncac.setdefault(rid, (n_idx, ca_idx, c_idx))

    missing = [rid for rid in residue_ids if rid not in residue_to_ncac]
    if missing:
        raise ValueError(f"Missing N/CA/C runtime mapping for residues: {missing[:8]}")

    affine_atoms = np.zeros((len(residue_ids), 3), dtype=np.int32)
    for seq_idx, resid in enumerate(residue_ids):
        affine_atoms[seq_idx, :] = residue_to_ncac[resid]
    return residue_ids, affine_atoms


def load_sidechain_rotamer_payload(sidechain_lib: Path, sequence, restype_to_index):
    if not sidechain_lib.exists():
        raise ValueError(f"Sidechain library not found: {sidechain_lib}")

    with h5py.File(sidechain_lib, "r") as sclib:
        sc_restype_order = decode_string_array(sclib["restype_order"])
        sc_restype_num = {name: i for i, name in enumerate(sc_restype_order)}
        start_stop_bead = sclib["rotamer_start_stop_bead"][:].astype(np.int32)
        rot_center_fixed = sclib["rotamer_center_fixed"][:, :6].astype(np.float32)
        # Upside's rotamer 1-body, -log P(rotamer | phi, psi) on the library's Ramachandran grid, exactly as
        # upside_config writes it under --dynamic-rotamer-1body: one (phi, psi) table per library row.
        rot_energy_rama = -np.log(sclib["rotamer_prob"][:].transpose((2, 0, 1)))[..., None]
        bead_order = decode_string_array(sclib["bead_order"])
        bead_num = {name: i for i, name in enumerate(bead_order)}
        pair_interaction = sclib["pair_interaction"][:].astype(np.float32)

    count_by_n_rot = {}
    affine_residue = []
    layer_index = []
    bead_type_index = []
    id_seq = []
    row_rotamer_index = []
    row_residue_table_index = []

    # Every residue has its Upside side chain, as upside_config writes it. A residue without a MARTINI side
    # chain carries table index -1, which martini_sc_table_1body reads as no environment term.
    for seq_idx, raw_resname in enumerate(sequence):
        resname = normalize_resname(raw_resname)
        if resname not in sc_restype_num:
            raise ValueError(f"Residue {seq_idx} ({resname}) is not in the sidechain library {sidechain_lib}")
        if resname in restype_to_index:
            residue_table_idx = restype_to_index[resname]
        elif resname in NO_MARTINI_SIDECHAIN:
            residue_table_idx = -1
        else:
            raise ValueError(f"Residue {seq_idx} ({resname}) has a MARTINI side chain but no side-chain/environment table")

        start, stop, n_bead = [int(x) for x in start_stop_bead[sc_restype_num[resname]]]
        n_rot = (stop - start) // n_bead

        if n_rot not in count_by_n_rot:
            count_by_n_rot[n_rot] = 0
        base_id = (count_by_n_rot[n_rot] << N_BIT_ROTAMER) + n_rot
        count_by_n_rot[n_rot] += 1

        for rel in range(stop - start):
            lid = start + rel
            rot_idx = rel // n_bead
            bead_name = f"{resname}_{rel % n_bead}"
            if bead_name not in bead_num:
                raise ValueError(f"Missing bead '{bead_name}' in sidechain library bead_order")
            affine_residue.append(seq_idx)
            layer_index.append(lid)
            bead_type_index.append(bead_num[bead_name])
            id_seq.append(rot_idx + (base_id << N_BIT_ROTAMER))
            row_rotamer_index.append(rot_idx)
            row_residue_table_index.append(residue_table_idx)

    if not layer_index:
        raise ValueError("No production rotamer rows could be generated from the sidechain library")

    # layer_index holds absolute library rows, so the placement tables are the whole library, as
    # upside_config writes them; the engine reads placement_data[layer_index[row]].
    return {
        "pair_interaction": pair_interaction,
        "affine_residue": np.asarray(affine_residue, dtype=np.int32),
        "layer_index": np.asarray(layer_index, dtype=np.int32),
        "bead_type_index": np.asarray(bead_type_index, dtype=np.int32),
        "id_seq": np.asarray(id_seq, dtype=np.int32),
        "placement_data": rot_center_fixed,
        "placement_scalar_data": rot_energy_rama,
        "row_rotamer_index": np.asarray(row_rotamer_index, dtype=np.int32),
        "row_residue_table_index": np.asarray(row_residue_table_index, dtype=np.int32),
    }


def build_env_rows(inp, target_to_index):
    if "hybrid_env_topology" not in inp:
        raise ValueError("Missing /input/hybrid_env_topology in stage file")
    env_grp = inp["hybrid_env_topology"]
    if "protein_membership" not in env_grp:
        raise ValueError("Missing hybrid_env_topology/protein_membership in stage file")

    protein_membership = env_grp["protein_membership"][:].astype(np.int32)
    atom_types = decode_string_array(inp["type"])
    charges = inp["charges"][:]

    # A table target is a bead type with its charge: the PC choline bead is Q0 with charge +1.
    env_atom_index = []
    env_target_index = []
    for atom_idx, (membership, atom_type, charge) in enumerate(zip(protein_membership.tolist(), atom_types, charges)):
        if membership >= 0:
            continue
        label = sc_target_label(atom_type, charge)
        if label not in target_to_index:
            raise ValueError(f"environment atom {atom_idx} ({label}) has no side-chain table target; "
                             "rebuild martini.h5 with its topology beside the force field")
        env_atom_index.append(atom_idx)
        env_target_index.append(target_to_index[label])

    return np.asarray(env_atom_index, dtype=np.int32), np.asarray(env_target_index, dtype=np.int32)


def recreate_group(parent, name):
    if name in parent:
        del parent[name]
    return parent.create_group(name)


def read_chain_first_residue(inp):
    """The first residue of every chain after the first: /input/chain_break, else the chain ids of hybrid_bb_map."""
    if "chain_break" in inp and "chain_first_residue" in inp["chain_break"]:
        return inp["chain_break"]["chain_first_residue"][:].astype(np.int32)
    if "hybrid_bb_map" in inp and "bb_chain_id" in inp["hybrid_bb_map"]:
        bb_chain_ids = decode_string_array(inp["hybrid_bb_map"]["bb_chain_id"])
        return derive_chain_break_metadata([{"bb_chain": chain_id} for chain_id in bb_chain_ids])[0]
    return np.array([], dtype=np.int32)


def load_upside_config(upside_home: Path):
    upside_config_py = upside_home / "py" / "upside_config.py"
    if not upside_config_py.exists():
        raise ValueError(f"upside_config.py not found: {upside_config_py}")
    spec = importlib.util.spec_from_file_location("upside_config_runtime", str(upside_config_py))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def remap_atom_index_array(arr, atom_map):
    """Map the 3-atoms-per-residue indices upside_config writes onto the hybrid's atoms; -1 stays -1."""
    out = arr.copy()
    valid = out >= 0
    mapped = atom_map[out[valid]]
    if np.any(mapped < 0):
        raise ValueError(f"Backbone atom remap has no target for reference index {int(out[valid][mapped < 0][0])}")
    out[valid] = mapped
    return out


def inject_backbone_nodes(
    up_file: Path,
    sequence,
    affine_atoms,
    chain_first_residue,
    rama_library: Path,
    rama_sheet_mixing: Path,
    hbond_energy: Path,
    reference_state_rama: Path,
    uc,
):
    fasta_seq = np.asarray(sequence)
    ref_n_atom = 3 * len(sequence)
    spring_args = types.SimpleNamespace(bond_stiffness=48.0, angle_stiffness=175.0, omega_stiffness=30.0)
    n_chains = int(chain_first_residue.size) + 1
    chain_starts = np.append(np.array([0], dtype=np.int32), chain_first_residue) * 3

    with tb.open_file(str(up_file), mode="a") as tf:
        uc.t = tf
        uc.potential = tf.root.input.potential
        uc.n_atom = ref_n_atom
        uc.n_chains = n_chains
        uc.chain_first_residue = chain_first_residue.astype(np.int32, copy=False)
        uc.chain_starts = chain_starts.astype(np.int32, copy=False)
        uc.rl_chains = None
        uc.use_intensive_memory = False

        uc.write_dist_spring(spring_args)
        uc.write_angle_spring(spring_args)
        uc.write_omega_spring1(spring_args, fasta_seq)
        uc.write_rama_map_pot(
            fasta_seq,
            str(rama_library),
            str(rama_sheet_mixing),
            secstr_bias="",
            mode="mixture",
            param_deriv=False,
        )

        with open(reference_state_rama, "rb") as fh:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=r"dtype\(\): align should be passed as Python or NumPy boolean.*",
                    category=VISIBLE_DEPRECATION_WARNING,
                )
                ref_state_raw = cPickle.load(fh, encoding="latin1")
        ref_state_cor = np.log(np.asarray(ref_state_raw, dtype=np.float64))
        ref_state_cor -= ref_state_cor.mean()
        grp = tf.create_group(tf.root.input.potential, "rama_map_pot_ref")
        grp._v_attrs.arguments = np.array([b"rama_coord"])
        grp._v_attrs.log_pot = 0
        uc.create_array(grp, "residue_id", obj=np.arange(len(fasta_seq), dtype=np.int32))
        uc.create_array(grp, "rama_map_id", obj=np.zeros(len(fasta_seq), dtype=np.int32))
        uc.create_array(grp, "rama_pot", obj=ref_state_cor[None])

        uc.write_infer_H_O(fasta_seq, [])
        uc.write_count_hbond(fasta_seq, False)
        uc.write_short_hbond(fasta_seq, str(hbond_energy))
        uc.write_rama_coord2()
        uc.write_backbone_pair(fasta_seq)

    atom_map = np.full((ref_n_atom,), -1, dtype=np.int64)
    for seq_idx, (n_idx, ca_idx, c_idx) in enumerate(affine_atoms.tolist()):
        atom_map[3 * seq_idx + 0] = n_idx
        atom_map[3 * seq_idx + 1] = ca_idx
        atom_map[3 * seq_idx + 2] = c_idx

    remap_datasets = [
        ("Distance3D", "id"),
        ("Angle", "id"),
        ("Dihedral_omega", "id"),
        ("Dihedral_phi", "id"),
        ("Dihedral_psi", "id"),
        ("rama_coord", "id"),
        ("infer_H_O/donors", "id"),
        ("infer_H_O/acceptors", "id"),
    ]
    with h5py.File(up_file, "r+") as up:
        pot = up["/input/potential"]
        for node_name in BACKBONE_NODES:
            if node_name not in pot:
                raise ValueError(f"Missing generated backbone node in stage file {up_file}: {node_name}")
        for grp_name, dset_name in remap_datasets:
            ds = pot[grp_name][dset_name]
            ds[...] = remap_atom_index_array(ds[:], atom_map).astype(ds.dtype, copy=False)


def published_pair_potential(r, eps, sig, qq, coulomb_k, r_cut):
    """The dry-MARTINI pair potential the particle tables hold: potential-shifted LJ plus reaction-field
    Coulomb with a conducting boundary (eps_rf -> infinity), both zero at and beyond r_cut."""
    r = np.asarray(r, dtype=np.float64)
    energy = np.zeros_like(r)
    inside = r < r_cut
    ri = r[inside]
    if eps != 0.0 and sig != 0.0:
        energy[inside] += 4.0 * eps * ((sig / ri) ** 12 - (sig / ri) ** 6) \
            - 4.0 * eps * ((sig / r_cut) ** 12 - (sig / r_cut) ** 6)
    if qq != 0.0:
        energy[inside] += coulomb_k * qq * (1.0 / ri + ri * ri / (2.0 * r_cut ** 3) - 1.5 / r_cut)
    return energy


def soft_core_radius(r, sig, soft_core_lambda):
    # (sc_alpha sigma^6 lambda^2 + r^6)^(1/6), finite at r = 0; lambda = 0 is r itself
    r = np.asarray(r, dtype=np.float64)
    if soft_core_lambda == 0.0:
        return r
    return (SOFT_CORE_ALPHA * soft_core_lambda ** SOFT_CORE_POWER * sig ** 6 + r ** 6) ** (1.0 / 6.0)


def soft_core_pair_potential(r, eps, sig, qq, coulomb_k, r_cut, soft_core_lambda):
    """GROMACS soft core (sc-power 2, sc-coul yes) of the published potential at coupling lambda:
    (1 - lambda) V(r_sc), r_sc = (sc_alpha sigma^6 lambda^2 + r^6)^(1/6). lambda = 0 is V itself."""
    r_sc = soft_core_radius(np.asarray(r, dtype=np.float64), sig, soft_core_lambda)
    return (1.0 - soft_core_lambda) * published_pair_potential(r_sc, eps, sig, qq, coulomb_k, r_cut)


def pair_table_end_slopes(r_ends, eps, sig, qq, coulomb_k, r_cut, soft_core_lambda):
    """dV/dr at the table ends, where its spline is clamped: (1 - lambda) V'(r_sc) (r / r_sc)^5, V' the derivative
    of the published expression without its cutoff. At r_max that is the limit from inside, where the
    potential-shifted LJ keeps a force; the soft core's r_sc passes the cutoff a few 1e-6 A short of r_max, so its
    table ends on the slope of its last smooth piece."""
    r = np.asarray(r_ends, dtype=np.float64)
    r_sc = soft_core_radius(r, sig, soft_core_lambda)
    deriv = np.zeros_like(r)
    if eps != 0.0 and sig != 0.0:
        deriv += 4.0 * eps * (-12.0 * (sig / r_sc) ** 12 + 6.0 * (sig / r_sc) ** 6) / r_sc
    if qq != 0.0:
        deriv += coulomb_k * qq * (-1.0 / (r_sc * r_sc) + r_sc / r_cut ** 3)
    return (1.0 - soft_core_lambda) * deriv * (r / r_sc) ** 5


def inject_particles_table(up_file: Path, martini_h5: Path, energy_conversion: float, length_conversion: float):
    """Copy the particle pair tables of `martini_h5` that the stage's coefficient rows use into the stage.

    A soft-core stage replaces each row by the soft core of the published potential, tabulated from r = 0, where
    it is finite with zero slope. Before that, every row of the library is checked against the published form, so
    the soft core is built from the same potential the library tabulates and lambda = 0 reproduces the library row.
    """
    up_file = Path(up_file).expanduser().resolve()
    martini_h5 = Path(martini_h5).expanduser().resolve()

    with h5py.File(martini_h5, "r") as src:
        pg = src["particles"]
        eps_arr = pg["unique_eps_eup"][:].astype(np.float64)
        sig_arr = pg["unique_sig_ang"][:].astype(np.float64)
        qq_arr = pg["unique_charge_product"][:].astype(np.float64)
        library_grids = pg["combined_energy_grids"][:].astype(np.float64)
        # converted once here, with the factors the stage was built with (the library stores it as float32)
        coulomb_k = COULOMB_F_KJ_NM / DRY_MARTINI_EPSILON_R * length_conversion / energy_conversion
        # The radial domain is the library's own; its outer edge is the cutoff. The attributes are float32,
        # so the knots are rebuilt from their shortest decimal spelling, the value the builder tabulated at.
        grid_n = int(pg.attrs["n_points"])
        r_min = float(str(np.float32(pg.attrs["r_min_ang"])))
        r_max = float(str(np.float32(pg.attrs["r_max_ang"])))
    if library_grids.shape[1] != grid_n:
        raise ValueError(f"{martini_h5} declares n_points={grid_n} but its grids have {library_grids.shape[1]} columns")
    r = np.linspace(r_min, r_max, grid_n)

    with h5py.File(up_file, "r+") as up:
        g = up["/input/potential/martini_potential"]
        soft_core_lambda = float(g.attrs["soft_core_lambda"])
        # the published potential diverges at r = 0, so its table starts at the library's r_min; the soft core does not
        r_lo = 0.0 if soft_core_lambda > 0.0 else r_min
        coeffs = g["coefficients"][:].astype(np.float64)
        rows = []
        for coeff_i, (eps, sig, q1, q2) in enumerate(coeffs[:, :4]):
            match = np.where(
                np.isclose(eps_arr, eps, rtol=1e-5, atol=1e-7)
                & np.isclose(sig_arr, sig, rtol=1e-5, atol=1e-6)
                & np.isclose(qq_arr, q1 * q2, rtol=1e-5, atol=1e-7)
            )[0]
            if match.size == 0:
                raise ValueError(f"{up_file}: no particle spline row for local coefficient {coeff_i}: "
                                 f"eps={eps:g}, sig={sig:g}, q1*q2={q1 * q2:g}")
            rows.append(int(match[0]))
        eps_arr, sig_arr, qq_arr, grids = eps_arr[rows], sig_arr[rows], qq_arr[rows], library_grids[rows]

        for i in range(len(rows)):
            reference = published_pair_potential(r, eps_arr[i], sig_arr[i], qq_arr[i], coulomb_k, r_max)
            scale = np.maximum(np.abs(reference), max(abs(eps_arr[i]), 1.0))
            deviation = float(np.max(np.abs(grids[i] - reference) / scale))
            if deviation > 1e-9:
                raise ValueError(f"{martini_h5} row (eps={eps_arr[i]:g}, sig={sig_arr[i]:g}, qq={qq_arr[i]:g}) "
                                 f"departs from the published dry-MARTINI form by {deviation:.3e}")
            if soft_core_lambda > 0.0:
                grids[i] = soft_core_pair_potential(np.linspace(r_lo, r_max, grid_n), eps_arr[i], sig_arr[i], qq_arr[i],
                                                    coulomb_k, r_max, soft_core_lambda)
            if grids[i, -1] != 0.0 and abs(grids[i, -1]) > 1e-9 * max(abs(eps_arr[i]), 1.0):
                raise ValueError(f"{up_file}: pair table row {i} is {grids[i, -1]:.3e} E_up at the {r_max:g} A cutoff")
        end_slopes = np.array([pair_table_end_slopes([r_lo, r_max], eps_arr[i], sig_arr[i], qq_arr[i], coulomb_k,
                                                     r_max, soft_core_lambda) for i in range(len(rows))])

        g.create_dataset("unique_eps_eup", data=eps_arr, dtype=np.float64)
        g.create_dataset("unique_sig_ang", data=sig_arr, dtype=np.float64)
        g.create_dataset("unique_charge_product", data=qq_arr, dtype=np.float64)
        g.create_dataset("combined_energy_grids", data=grids, dtype=np.float64)
        g.create_dataset("combined_energy_end_slopes", data=end_slopes, dtype=np.float64)
        g.attrs["r_min_ang"] = np.float32(r_lo)
        g.attrs["r_max_ang"] = np.float32(r_max)
        g.attrs["n_points"] = np.int32(grid_n)


def inject_stage7_sc_table_nodes(
    up_file: Path,
    martini_h5: Path,
    upside_home: Path,
    rama_library: Path,
    rama_sheet_mixing: Path,
    hbond_energy: Path,
    reference_state_rama: Path,
):
    up_file = Path(up_file).expanduser().resolve()
    martini_h5 = Path(martini_h5).expanduser().resolve()
    upside_home = Path(upside_home).expanduser().resolve()
    rama_library = Path(rama_library).expanduser().resolve()
    rama_sheet_mixing = Path(rama_sheet_mixing).expanduser().resolve()
    hbond_energy = Path(hbond_energy).expanduser().resolve()
    reference_state_rama = Path(reference_state_rama).expanduser().resolve()
    sidechain_lib = (upside_home / "parameters" / "ff_2.1" / "sidechain.h5").resolve()

    if not up_file.exists():
        raise SystemExit(f"ERROR: stage file not found: {up_file}")
    if not martini_h5.exists():
        raise SystemExit(f"ERROR: martini.h5 not found: {martini_h5}")
    for path in [rama_library, rama_sheet_mixing, hbond_energy, reference_state_rama, sidechain_lib]:
        if not path.exists():
            raise SystemExit(f"ERROR: required Upside input not found: {path}")

    with h5py.File(martini_h5, "r") as sc_lib:
        sc_grp = sc_lib["sc_table"] if "sc_table" in sc_lib else sc_lib
        required_sc_datasets = [
            "restype_order",
            "target_order",
            "grid_ang",
            "theta_grid",
            "rotamer_count",
            "rotamer_full_energy_eup",
        ]
        missing_sc_datasets = [name for name in required_sc_datasets if name not in sc_grp]
        if missing_sc_datasets:
            missing_text = ", ".join(missing_sc_datasets)
            raise SystemExit(
                f"ERROR: {martini_h5} is missing required rotamer-resolved SC datasets: {missing_text}. "
                "Regenerate martini.h5 via build_martini_h5()."
            )
        restype_order = decode_string_array(sc_grp["restype_order"])
        target_order = decode_string_array(sc_grp["target_order"])
        grid_ang = sc_grp["grid_ang"][:].astype(np.float32)
        theta_grid = sc_grp["theta_grid"][:].astype(np.float32)
        rotamer_count = sc_grp["rotamer_count"][:].astype(np.int32)
        rotamer_full_energy_eup = sc_grp["rotamer_full_energy_eup"][:].astype(np.float32)

    restype_to_index = {name: i for i, name in enumerate(restype_order)}
    target_to_index = {name: i for i, name in enumerate(target_order)}
    uc = load_upside_config(upside_home)

    with h5py.File(up_file, "r+") as up:
        inp = up["input"]
        pot = inp["potential"]
        martini_potential = pot["martini_potential"]

        residue_ids, affine_atoms = build_affine_atoms(inp)
        sequence = resolve_sequence(inp, len(residue_ids))
        chain_first_residue = read_chain_first_residue(inp)
        env_atom_index, env_target_index = build_env_rows(inp, target_to_index)
        rotamer_payload = load_sidechain_rotamer_payload(sidechain_lib, sequence, restype_to_index)

        if "hybrid_sc_map" in inp:
            del inp["hybrid_sc_map"]
        if "sequence" in inp:
            del inp["sequence"]
        inp.create_dataset("sequence", data=np.asarray([np.bytes_(x) for x in sequence], dtype="S3"))

        g_aff = recreate_group(pot, "affine_alignment")
        g_aff.attrs["arguments"] = np.asarray([np.bytes_("pos")])
        g_aff.create_dataset("atoms", data=affine_atoms, dtype=np.int32)
        g_aff.create_dataset(
            "ref_geom",
            data=np.repeat(CANONICAL_AFFINE_REF[None, :, :], len(sequence), axis=0),
            dtype=np.float32,
        )

        g_sc_place = recreate_group(pot, "placement_fixed_point_vector_only")
        g_sc_place.attrs["arguments"] = np.asarray([np.bytes_("affine_alignment")])
        g_sc_place.create_dataset("affine_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_sc_place.create_dataset("layer_index", data=rotamer_payload["layer_index"], dtype=np.int32)
        g_sc_place.create_dataset("placement_data", data=rotamer_payload["placement_data"], dtype=np.float32)

        g_pl = recreate_group(pot, "placement_scalar")
        g_pl.attrs["arguments"] = np.asarray([np.bytes_("affine_alignment"), np.bytes_("rama_coord")])
        g_pl.create_dataset("rama_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_pl.create_dataset("affine_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_pl.create_dataset("layer_index", data=rotamer_payload["layer_index"], dtype=np.int32)
        g_pl.create_dataset("placement_data", data=rotamer_payload["placement_scalar_data"], dtype=np.float32)

        g_rot = recreate_group(pot, "rotamer")
        g_rot.attrs["arguments"] = np.asarray(
            [
                np.bytes_("placement_fixed_point_vector_only"),
                np.bytes_("placement_scalar"),
                np.bytes_("martini_sc_table_1body"),
            ]
        )
        g_rot.attrs["integrator_level"] = np.int32(1)
        g_rot.attrs["max_iter"] = np.int32(1000)
        g_rot.attrs["tol"] = np.float32(1.0e-3)
        g_rot.attrs["damping"] = np.float32(0.4)
        g_rot.attrs["iteration_chunk_size"] = np.int32(2)

        g_pair = g_rot.create_group("pair_interaction")
        g_pair.create_dataset("interaction_param", data=rotamer_payload["pair_interaction"], dtype=np.float32)
        g_pair.create_dataset(
            "index", data=np.arange(len(rotamer_payload["id_seq"]), dtype=np.int32), dtype=np.int32
        )
        g_pair.create_dataset("type", data=rotamer_payload["bead_type_index"], dtype=np.int32)
        g_pair.create_dataset("id", data=rotamer_payload["id_seq"], dtype=np.int32)
        g_pair.create_dataset(
            "sequence_id", data=uc.sequence_id(rotamer_payload["affine_residue"], chain_first_residue), dtype=np.int32
        )

        g_cb = recreate_group(pot, "placement_fixed_point_vector_only_CB")
        g_cb.attrs["arguments"] = np.asarray([np.bytes_("affine_alignment")])
        g_cb.create_dataset("affine_residue", data=np.arange(len(sequence), dtype=np.int32), dtype=np.int32)
        g_cb.create_dataset("layer_index", data=np.zeros(len(sequence), dtype=np.int32), dtype=np.int32)
        g_cb.create_dataset(
            "placement_data", data=np.concatenate([CB_PLACEMENT, CB_VECTOR], axis=1), dtype=np.float32
        )

    inject_backbone_nodes(
        up_file=up_file,
        sequence=sequence,
        affine_atoms=affine_atoms,
        chain_first_residue=chain_first_residue,
        rama_library=rama_library,
        rama_sheet_mixing=rama_sheet_mixing,
        hbond_energy=hbond_energy,
        reference_state_rama=reference_state_rama,
        uc=uc,
    )

    with h5py.File(up_file, "r+") as up:
        inp = up["input"]
        pot = inp["potential"]
        martini_potential = pot["martini_potential"]
        g_sc = recreate_group(pot, "martini_sc_table_1body")
        g_sc.attrs["arguments"] = np.asarray([np.bytes_("pos"), np.bytes_("placement_fixed_point_vector_only_CB")])
        g_sc.attrs["x_len"] = np.float32(martini_potential.attrs["x_len"])
        g_sc.attrs["y_len"] = np.float32(martini_potential.attrs["y_len"])
        g_sc.attrs["z_len"] = np.float32(martini_potential.attrs["z_len"])
        g_sc.create_dataset("row_residue_index", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_sc.create_dataset("row_rotamer_index", data=rotamer_payload["row_rotamer_index"], dtype=np.int32)
        g_sc.create_dataset(
            "row_residue_table_index", data=rotamer_payload["row_residue_table_index"], dtype=np.int32
        )
        g_sc.create_dataset("env_atom_index", data=env_atom_index, dtype=np.int32)
        g_sc.create_dataset("env_target_index", data=env_target_index, dtype=np.int32)
        g_sc.create_dataset("grid_ang", data=grid_ang, dtype=np.float32)
        g_sc.create_dataset("theta_grid", data=theta_grid, dtype=np.float32)
        g_sc.create_dataset("rotamer_count", data=rotamer_count.astype(np.int32), dtype=np.int32)
        g_sc.create_dataset("rotamer_full_energy_eup", data=rotamer_full_energy_eup, dtype=np.float32)

    print(
        f"Injected rotamer-weighted martini_sc_table_1body into {up_file}: "
        f"n_rows={len(rotamer_payload['id_seq'])} "
        f"without_environment_table={int(np.sum(rotamer_payload['row_residue_table_index'] < 0))} "
        f"n_env={len(env_atom_index)} n_restypes={len(restype_order)} n_targets={len(target_order)}"
    )


if __name__ == "__main__":
    raise SystemExit("Use martini_prepare_system.py as the workflow preparation entrypoint.")
