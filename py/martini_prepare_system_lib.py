#!/usr/bin/env python3
from __future__ import annotations

import _pickle as cPickle
import importlib.util
import json
import os
import sys
import types
import warnings
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path

import math

import h5py
import numpy as np
import tables as tb

PY_DIR = Path(__file__).resolve().parent
REPO_ROOT = PY_DIR.parent
WORKFLOW_DIR = REPO_ROOT / "example" / "16.MARTINI"

NA_AVOGADRO = 6.02214076e23
BB_COMPONENT_NAMES = ("N", "CA", "C", "O")
BB_COMPONENT_MASSES = (14.0, 12.0, 12.0, 16.0)
BB_TYPE_CHARGE = {
    "Qd": 1.0,
    "Qa": -1.0,
    "SQd": 1.0,
    "SQa": -1.0,
    "RQd": 1.0,
    "AQa": -1.0,
}
TWOPI = 2.0 * np.pi
CANONICAL_AFFINE_REF = np.array(
    [
        [-1.19280531, -0.83127186, 0.0],
        [0.0, 0.0, 0.0],
        [1.25222632, -0.87268266, 0.0],
    ],
    dtype=np.float32,
)
CANONICAL_AFFINE_REF -= CANONICAL_AFFINE_REF.mean(axis=0, keepdims=True)


def _cgl_effective_lj_sigma_cap_nm() -> float:
    raw = os.environ.get("UPSIDE_CG_LIPID_MAX_EFFECTIVE_SIGMA_NM", "0.9")
    try:
        cap = float(raw)
    except ValueError as exc:
        raise ValueError(f"Invalid UPSIDE_CG_LIPID_MAX_EFFECTIVE_SIGMA_NM={raw!r}") from exc
    if cap <= 0.0:
        return float("inf")
    return cap
CB_PLACEMENT = np.array([[0.0, 0.94375626, 1.2068012]], dtype=np.float32)
CB_VECTOR = CB_PLACEMENT / np.linalg.norm(CB_PLACEMENT, axis=1, keepdims=True)
N_BIT_ROTAMER = 4
LEGACY_STAGE7_NODES = [
    "rotamer",
    "placement_fixed_scalar",
    "placement_fixed_point_vector_only",
    "martini_sc_table_potential",
    "martini_sc_table_1body",
]
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


def _decode_h5_text_array(values):
    arr = np.asarray(values)
    out = []
    for value in arr:
        if isinstance(value, (bytes, np.bytes_)):
            out.append(value.decode("utf-8", errors="ignore").strip())
        else:
            out.append(str(value).strip())
    return np.asarray(out, dtype=object)


def _stage_pos_array(raw):
    pos = np.asarray(raw)
    if pos.ndim == 3 and pos.shape[1] == 3 and pos.shape[2] >= 1:
        return pos[:, :, 0].astype(np.float64)
    if pos.ndim == 4 and pos.shape[1] == 1 and pos.shape[2] == 3:
        return pos[:, 0, :, 0].astype(np.float64)
    if pos.ndim == 2 and pos.shape[1] == 3:
        return pos.astype(np.float64)
    raise ValueError(f"Unsupported /input/pos shape for debug PDB export: {pos.shape}")


def _debug_box_lengths(input_grp):
    pot = input_grp.get("potential")
    if pot is not None:
        for node_name in ("martini_potential", "dist_spring", "compose_vector6d"):
            node = pot.get(node_name)
            if node is None:
                continue
            attrs = node.attrs
            if all(name in attrs for name in ("x_len", "y_len", "z_len")):
                box = np.array(
                    [attrs["x_len"], attrs["y_len"], attrs["z_len"]],
                    dtype=np.float64,
                )
                if np.all(np.isfinite(box)) and np.all(box > 0.0):
                    return box
    return None


def _minimum_image_delta(delta, box_lengths):
    delta = np.asarray(delta, dtype=np.float64).copy()
    if box_lengths is None:
        return delta
    box = np.asarray(box_lengths, dtype=np.float64)
    valid = box > 0.0
    delta[..., valid] -= box[valid] * np.round(delta[..., valid] / box[valid])
    return delta


def _debug_element(atom_name, atom_type):
    name = str(atom_name).strip().upper()
    atype = str(atom_type).strip().upper()
    for value in (name, atype):
        if value in {"NA", "CL"}:
            return value.title() if value == "NA" else "Cl"
    letters = "".join(ch for ch in name if ch.isalpha())
    if not letters:
        letters = "".join(ch for ch in atype if ch.isalpha())
    if not letters:
        return ""
    if letters.startswith("CL"):
        return "Cl"
    if letters.startswith("NA"):
        return "Na"
    return letters[0].upper()


def _debug_residue_name(atom_type, atom_class):
    atype = str(atom_type).strip().upper()
    cls = str(atom_class).strip().upper()
    if atype == "CGLD":
        return "CGLD"
    if atype == "CGL" or cls == "LIPID":
        return "DOPC"
    if cls == "PROTEIN":
        return "PRO"
    if cls == "ION":
        return atype if atype in {"NA", "CL"} else "ION"
    if cls == "WATER":
        return "W"
    return atype[:4] if atype else "UNK"


def _debug_chain_id(atom_type, atom_class):
    atype = str(atom_type).strip().upper()
    cls = str(atom_class).strip().upper()
    if atype == "CGLD":
        return "D"
    if cls == "PROTEIN":
        return "P"
    if cls == "LIPID":
        return "L"
    if cls == "ION":
        return "I"
    if cls == "WATER":
        return "W"
    return "X"


def _debug_records(
    pos,
    atom_types,
    atom_names,
    particle_classes,
    residue_ids,
    charges,
    include_hidden,
    include_indices=None,
):
    include_set = None if include_indices is None else set(int(i) for i in include_indices)
    atoms = []
    for idx, xyz in enumerate(pos):
        if include_set is not None and idx not in include_set:
            continue
        atom_type = str(atom_types[idx]).strip()
        atom_name = str(atom_names[idx]).strip() if atom_names is not None else atom_type
        atom_class = str(particle_classes[idx]).strip() if particle_classes is not None else "OTHER"
        if not include_hidden and atom_type.upper() == "CGLD":
            continue
        pdb_name = "CGLD" if atom_type.upper() == "CGLD" else (atom_name or atom_type or "X")
        atoms.append(
            {
                "record": "ATOM" if atom_class.upper() == "PROTEIN" else "HETATM",
                "name": pdb_name[:4],
                "resname": _debug_residue_name(atom_type, atom_class),
                "chain": _debug_chain_id(atom_type, atom_class),
                "resseq": int(residue_ids[idx]) if residue_ids is not None else idx + 1,
                "icode": "",
                "x": float(xyz[0]),
                "y": float(xyz[1]),
                "z": float(xyz[2]),
                "occ": 1.0,
                "bfac": float(charges[idx]) if charges is not None else 0.0,
                "segid": atom_class[:4].upper(),
                "element": _debug_element(pdb_name, atom_type),
                "charge": "",
            }
        )
    return atoms


def _nearest_xy_stats(xy, box_lengths):
    if xy.shape[0] < 2:
        return {"min": None, "p05": None, "median": None}
    values = []
    box_xy = np.asarray(box_lengths[:2], dtype=np.float64) if box_lengths is not None else None
    for i in range(xy.shape[0]):
        delta = xy - xy[i]
        if box_xy is not None and np.all(box_xy > 0.0):
            delta -= box_xy * np.round(delta / box_xy)
        dist = np.sqrt(np.sum(delta * delta, axis=1))
        dist[i] = np.inf
        values.append(float(np.min(dist)))
    arr = np.asarray(values, dtype=np.float64)
    return {
        "min": float(np.min(arr)),
        "p05": float(np.percentile(arr, 5.0)),
        "median": float(np.median(arr)),
    }


def _debug_leaflet_stats(pos, atom_types, box_lengths):
    cgl_idx = np.where(np.asarray([str(x).strip().upper() == "CGL" for x in atom_types]))[0]
    if cgl_idx.size == 0:
        return {}
    cgl_pos = pos[cgl_idx]
    split = float(np.median(cgl_pos[:, 2]))
    out = {
        "cgl_count": int(cgl_idx.size),
        "z_min": float(np.min(cgl_pos[:, 2])),
        "z_max": float(np.max(cgl_pos[:, 2])),
        "z_mean": float(np.mean(cgl_pos[:, 2])),
        "z_std": float(np.std(cgl_pos[:, 2])),
        "median_split_z": split,
    }
    for name, mask in (("lower", cgl_pos[:, 2] <= split), ("upper", cgl_pos[:, 2] > split)):
        if not np.any(mask):
            out[name] = {"count": 0}
            continue
        leaflet = cgl_pos[mask]
        out[name] = {
            "count": int(leaflet.shape[0]),
            "z_mean": float(np.mean(leaflet[:, 2])),
            "z_std": float(np.std(leaflet[:, 2])),
            "nearest_xy": _nearest_xy_stats(leaflet[:, :2], box_lengths),
        }
    return out


def _debug_cgl_partner_stats(
    pos,
    atom_types,
    particle_classes,
    box_lengths,
    martini_pairs=None,
    martini_coefficient_indices=None,
    martini_coefficients=None,
):
    cgl_mask = np.asarray([str(x).strip().upper() == "CGL" for x in atom_types], dtype=bool)
    if not np.any(cgl_mask):
        return {}

    class_upper = np.asarray([str(x).strip().upper() for x in particle_classes], dtype=object)
    partner_masks = {
        "ion": class_upper == "ION",
        "protein": class_upper == "PROTEIN",
    }
    out = {}
    cgl_idx = np.where(cgl_mask)[0]
    cgl_pos = pos[cgl_idx]
    for name, mask in partner_masks.items():
        partner_idx = np.where(mask)[0]
        if partner_idx.size == 0:
            out[name] = {
                "partner_count": 0,
                "min_distance_ang": None,
                "lj_pair_count": 0,
                "lj_sum_eup": None,
                "lj_max_eup": None,
            }
            continue

        partner_pos = pos[partner_idx]
        delta = partner_pos[:, None, :] - cgl_pos[None, :, :]
        delta = _minimum_image_delta(delta.reshape(-1, 3), box_lengths).reshape(
            partner_idx.size, cgl_idx.size, 3
        )
        dist = np.sqrt(np.sum(delta * delta, axis=2))
        dist_flat = dist.reshape(-1)
        stats = {
            "partner_count": int(partner_idx.size),
            "min_distance_ang": float(np.min(dist)),
            "cgl_count_within_ang": {
                f"{radius:g}": int(np.count_nonzero(np.any(dist <= radius, axis=0)))
                for radius in (6.0, 8.0, 10.0, 12.0, 14.0, 16.0)
            },
            "pair_annulus_counts_ang": {
                f"{lo:g}-{hi:g}": int(np.count_nonzero((dist_flat >= lo) & (dist_flat < hi)))
                for lo, hi in ((0.0, 10.0), (10.0, 15.0), (15.0, 20.0),
                               (20.0, 25.0), (25.0, 30.0), (30.0, 40.0),
                               (40.0, 55.0))
            },
            "lj_pair_count": 0,
            "lj_sum_eup": None,
            "lj_max_eup": None,
        }

        if (
            martini_pairs is not None
            and martini_coefficient_indices is not None
            and martini_coefficients is not None
        ):
            pair_cgl = cgl_mask[martini_pairs[:, 0]] | cgl_mask[martini_pairs[:, 1]]
            pair_partner = mask[martini_pairs[:, 0]] | mask[martini_pairs[:, 1]]
            pair_mask = pair_cgl & pair_partner
            if np.any(pair_mask):
                pairs = martini_pairs[pair_mask]
                coeff_idx = martini_coefficient_indices[pair_mask]
                coeff = martini_coefficients[coeff_idx]
                d = pos[pairs[:, 1]] - pos[pairs[:, 0]]
                d = _minimum_image_delta(d, box_lengths)
                r = np.maximum(np.sqrt(np.sum(d * d, axis=1)), 1.0e-6)
                eps = coeff[:, 0].astype(np.float64)
                sigma = coeff[:, 1].astype(np.float64)
                active = (eps != 0.0) & (sigma != 0.0)
                lj = np.zeros(r.shape[0], dtype=np.float64)
                sr = sigma[active] / r[active]
                sr2 = sr * sr
                sr6 = sr2 * sr2 * sr2
                lj[active] = 4.0 * eps[active] * (sr6 * sr6 - sr6)
                stats.update({
                    "lj_pair_count": int(r.shape[0]),
                    "lj_sum_eup": float(np.sum(lj)),
                    "lj_max_eup": float(np.max(lj)),
                })
        out[name] = stats
    return out


def _debug_scalar_attr_value(value):
    if isinstance(value, bytes):
        return value.decode("ascii", errors="replace")
    arr = np.asarray(value)
    if arr.shape == ():
        item = arr.item()
        if isinstance(item, bytes):
            return item.decode("ascii", errors="replace")
        if isinstance(item, np.generic):
            return item.item()
        return item
    return [
        x.decode("ascii", errors="replace") if isinstance(x, bytes) else x.item() if isinstance(x, np.generic) else x
        for x in arr.tolist()
    ]


def _debug_selected_attrs(group, attr_names):
    out = {}
    for attr_name in attr_names:
        if attr_name in group.attrs:
            out[attr_name] = _debug_scalar_attr_value(group.attrs[attr_name])
    return out


def _cg_lipid_sc_runtime_cutoff_ang(sc_attrs):
    if "cutoff_ang" not in sc_attrs:
        return None
    cutoff_ang = float(sc_attrs["cutoff_ang"])
    taper_width_ang = float(sc_attrs.get("taper_width_ang", sc_attrs.get("knot_spacing_ang", 0.0)))
    if taper_width_ang <= 0.0:
        return cutoff_ang
    if "fit_r_max_nm" not in sc_attrs:
        return cutoff_ang
    fitted_support_ang = float(sc_attrs["fit_r_max_nm"]) * 10.0 + taper_width_ang
    return min(cutoff_ang, fitted_support_ang)


def write_stage_debug_outputs(up_file: Path, debug_dir: Path | None = None, prefix: str | None = None):
    """Write PDB/topology diagnostics for a generated MARTINI stage input."""
    enabled = os.environ.get("UPSIDE_WRITE_DEBUG_PDB", "1").strip().lower()
    if enabled in {"0", "false", "no", "off"}:
        return {}

    up_file = Path(up_file)
    if debug_dir is None:
        debug_dir = up_file.parent / "debug"
    debug_dir = Path(debug_dir)
    debug_dir.mkdir(parents=True, exist_ok=True)
    prefix = prefix or up_file.stem

    with h5py.File(up_file, "r") as h5:
        input_grp = h5["/input"]
        pos = _stage_pos_array(input_grp["pos"][:])
        n_atom = int(pos.shape[0])
        atom_types = _decode_h5_text_array(input_grp["type"][:])
        atom_names = (
            _decode_h5_text_array(input_grp["atom_names"][:])
            if "atom_names" in input_grp else atom_types.copy()
        )
        particle_classes = (
            _decode_h5_text_array(input_grp["particle_class"][:])
            if "particle_class" in input_grp else np.asarray(["OTHER"] * n_atom, dtype=object)
        )
        residue_ids = (
            np.asarray(input_grp["residue_ids"][:], dtype=np.int32)
            if "residue_ids" in input_grp else np.arange(1, n_atom + 1, dtype=np.int32)
        )
        molecule_ids = (
            np.asarray(input_grp["molecule_ids"][:], dtype=np.int32)
            if "molecule_ids" in input_grp else np.full(n_atom, -1, dtype=np.int32)
        )
        charges = (
            np.asarray(input_grp["charges"][:], dtype=np.float64)
            if "charges" in input_grp else np.zeros(n_atom, dtype=np.float64)
        )
        box_lengths = _debug_box_lengths(input_grp)

        pot = input_grp.get("potential")
        bond_ids = np.zeros((0, 2), dtype=np.int32)
        bond_r0 = np.zeros(0, dtype=np.float64)
        bond_k = np.zeros(0, dtype=np.float64)
        if pot is not None and "dist_spring" in pot:
            spring = pot["dist_spring"]
            bond_ids = np.asarray(spring["id"][:], dtype=np.int32)
            if "equil_dist" in spring:
                bond_r0 = np.asarray(spring["equil_dist"][:], dtype=np.float64)
            if "spring_const" in spring:
                bond_k = np.asarray(spring["spring_const"][:], dtype=np.float64)

        martini_pairs = None
        martini_coefficient_indices = None
        martini_coefficients = None
        if pot is not None and "martini_potential" in pot:
            mpot = pot["martini_potential"]
            if (
                "pairs" in mpot
                and "coefficient_indices" in mpot
                and "coefficients" in mpot
            ):
                martini_pairs = np.asarray(mpot["pairs"][:], dtype=np.int32)
                martini_coefficient_indices = np.asarray(
                    mpot["coefficient_indices"][:], dtype=np.int64
                )
                martini_coefficients = np.asarray(mpot["coefficients"][:], dtype=np.float64)

        cg_lipid_sc_attrs = {}
        if pot is not None and "cg_lipid_sc" in pot:
            cg_lipid_sc_attrs = _debug_selected_attrs(
                pot["cg_lipid_sc"],
                (
                    "schema",
                    "radial_mode",
                    "angle_convention",
                    "fit_r_max_nm",
                    "knot_spacing_ang",
                    "cutoff_ang",
                    "taper_width_ang",
                    "n_modes",
                    "n_radial",
                    "n_angular",
                ),
            )

        orientation_pairs = np.zeros((0, 2), dtype=np.int32)
        orientation_lengths = np.zeros(0, dtype=np.float64)
        if pot is not None and "compose_vector6d" in pot:
            compose = pot["compose_vector6d"]
            if "elem_index" in compose and "orientation_index" in compose:
                elem = np.asarray(compose["elem_index"][:], dtype=np.int32)
                orient = np.asarray(compose["orientation_index"][:], dtype=np.int32)
                if elem.shape == orient.shape:
                    orientation_pairs = np.column_stack([elem, orient]).astype(np.int32)
                    d = _minimum_image_delta(pos[orient] - pos[elem], box_lengths)
                    orientation_lengths = np.sqrt(np.sum(d * d, axis=1))

    all_pdb = debug_dir / f"{prefix}.input_debug.all.pdb"
    visible_pdb = debug_dir / f"{prefix}.input_debug.visible.pdb"
    bilayer_pdb = debug_dir / f"{prefix}.input_debug.bilayer.pdb"
    bonds_tsv = debug_dir / f"{prefix}.input_debug.bonds.tsv"
    summary_json = debug_dir / f"{prefix}.input_debug.summary.json"
    cgl_indices = np.where(np.asarray([str(x).strip().upper() == "CGL" for x in atom_types]))[0]

    write_pdb(
        all_pdb,
        _debug_records(pos, atom_types, atom_names, particle_classes, residue_ids, charges, True),
        box_lengths,
    )
    write_pdb(
        visible_pdb,
        _debug_records(pos, atom_types, atom_names, particle_classes, residue_ids, charges, False),
        box_lengths,
    )
    write_pdb(
        bilayer_pdb,
        _debug_records(
            pos,
            atom_types,
            atom_names,
            particle_classes,
            residue_ids,
            charges,
            False,
            include_indices=cgl_indices,
        ),
        box_lengths,
    )

    warnings_list = []
    ion_bond_count = 0
    cgl_cgld_bond_count = 0
    unexpected_bond_count = 0
    with bonds_tsv.open("w", encoding="utf-8") as out:
        out.write(
            "bond_index\ti\tj\tserial_i\tserial_j\tname_i\tname_j\ttype_i\ttype_j\t"
            "class_i\tclass_j\tresidue_i\tresidue_j\tmolecule_i\tmolecule_j\t"
            "distance_ang\tequil_dist_ang\tspring_const\n"
        )
        for b_idx, (i_raw, j_raw) in enumerate(bond_ids):
            i = int(i_raw)
            j = int(j_raw)
            if i < 0 or j < 0 or i >= n_atom or j >= n_atom:
                warnings_list.append(f"Bond {b_idx} has out-of-range atom ids: {i},{j}")
                continue
            distance = float(np.linalg.norm(_minimum_image_delta(pos[j] - pos[i], box_lengths)))
            r0 = float(bond_r0[b_idx]) if b_idx < bond_r0.shape[0] else float("nan")
            spring_k = float(bond_k[b_idx]) if b_idx < bond_k.shape[0] else float("nan")
            class_i = str(particle_classes[i])
            class_j = str(particle_classes[j])
            type_i = str(atom_types[i])
            type_j = str(atom_types[j])
            if class_i.upper() == "ION" or class_j.upper() == "ION":
                ion_bond_count += 1
            pair_types = {type_i.upper(), type_j.upper()}
            if pair_types == {"CGL", "CGLD"}:
                cgl_cgld_bond_count += 1
            else:
                unexpected_bond_count += 1
            out.write(
                f"{b_idx}\t{i}\t{j}\t{i + 1}\t{j + 1}\t{atom_names[i]}\t{atom_names[j]}\t"
                f"{type_i}\t{type_j}\t{class_i}\t{class_j}\t{int(residue_ids[i])}\t"
                f"{int(residue_ids[j])}\t{int(molecule_ids[i])}\t{int(molecule_ids[j])}\t"
                f"{distance:.6f}\t{r0:.6f}\t{spring_k:.6f}\n"
            )

    type_counts = {str(k): int(v) for k, v in Counter(str(x) for x in atom_types).items()}
    class_counts = {str(k): int(v) for k, v in Counter(str(x) for x in particle_classes).items()}
    atom_name_counts = {str(k): int(v) for k, v in Counter(str(x) for x in atom_names).items()}
    cgl_count = int(type_counts.get("CGL", 0))
    cgld_count = int(type_counts.get("CGLD", 0))
    if ion_bond_count:
        warnings_list.append(f"ERROR: {ion_bond_count} true dist_spring bonds involve ions")
    if cgl_count != cgld_count:
        warnings_list.append(f"CGL/CGLD count mismatch: CGL={cgl_count}, CGLD={cgld_count}")
    if cgld_count and cgl_cgld_bond_count != cgld_count:
        warnings_list.append(
            f"CGL-CGLD bond count {cgl_cgld_bond_count} does not match CGLD count {cgld_count}"
        )
    if unexpected_bond_count:
        warnings_list.append(f"{unexpected_bond_count} dist_spring bonds are not CGL-CGLD orientation bonds")
    if orientation_lengths.size:
        target = None
        with h5py.File(up_file, "r") as h5:
            compose = h5["/input/potential/compose_vector6d"]
            if "orientation_length_ang" in compose.attrs:
                target = float(compose.attrs["orientation_length_ang"])
        if target is not None:
            max_dev = float(np.max(np.abs(orientation_lengths - target)))
            if max_dev > 0.25:
                warnings_list.append(
                    f"CGLD orientation length max deviation {max_dev:.6f} A from target {target:.6f} A"
                )
    cgl_cgld_same_molecule_count = 0
    cgl_cgld_molecule_mismatch_count = 0
    if orientation_pairs.size:
        same_mol = molecule_ids[orientation_pairs[:, 0]] == molecule_ids[orientation_pairs[:, 1]]
        cgl_cgld_same_molecule_count = int(np.count_nonzero(same_mol))
        cgl_cgld_molecule_mismatch_count = int(same_mol.size - cgl_cgld_same_molecule_count)
        if cgl_cgld_molecule_mismatch_count:
            warnings_list.append(
                f"{cgl_cgld_molecule_mismatch_count} CGLD sites do not share the parent CGL molecule_id; "
                "regenerate this input with the patched converter"
            )
    cgl_partner_stats = _debug_cgl_partner_stats(
        pos,
        atom_types,
        particle_classes,
        box_lengths,
        martini_pairs=martini_pairs,
        martini_coefficient_indices=martini_coefficient_indices,
        martini_coefficients=martini_coefficients,
    )
    ion_lj_max = cgl_partner_stats.get("ion", {}).get("lj_max_eup")
    if ion_lj_max is not None and ion_lj_max > 100.0:
        warnings_list.append(
            f"CGL-ion LJ max {ion_lj_max:.3f} E_up exceeds 100 E_up; "
            "check ion placement cutoff and CGL effective LJ sigma"
        )

    summary = {
        "up_file": str(up_file),
        "outputs": {
            "all_pdb": str(all_pdb),
            "visible_pdb": str(visible_pdb),
            "bilayer_pdb": str(bilayer_pdb),
            "bonds_tsv": str(bonds_tsv),
            "summary_json": str(summary_json),
        },
        "n_atom": n_atom,
        "box_lengths": None if box_lengths is None else [float(x) for x in box_lengths],
        "type_counts": type_counts,
        "atom_name_counts": atom_name_counts,
        "particle_class_counts": class_counts,
        "molecule_id_count": int(len(set(int(x) for x in molecule_ids))),
        "bond_count": int(bond_ids.shape[0]),
        "ion_bond_count": int(ion_bond_count),
        "cgl_count": cgl_count,
        "cgld_count": cgld_count,
        "cgl_cgld_bond_count": int(cgl_cgld_bond_count),
        "cgl_cgld_same_molecule_count": int(cgl_cgld_same_molecule_count),
        "cgl_cgld_molecule_mismatch_count": int(cgl_cgld_molecule_mismatch_count),
        "unexpected_bond_count": int(unexpected_bond_count),
        "orientation_pair_count": int(orientation_pairs.shape[0]),
        "orientation_length": {
            "min": float(np.min(orientation_lengths)) if orientation_lengths.size else None,
            "max": float(np.max(orientation_lengths)) if orientation_lengths.size else None,
            "mean": float(np.mean(orientation_lengths)) if orientation_lengths.size else None,
        },
        "leaflet_stats": _debug_leaflet_stats(pos, atom_types, box_lengths),
        "cgl_partner_stats": cgl_partner_stats,
        "cg_lipid_sc_attrs": cg_lipid_sc_attrs,
        "hidden_particle_types": {
            "CGLD": int(cgld_count),
        },
        "warnings": warnings_list,
    }
    with summary_json.open("w", encoding="utf-8") as out:
        json.dump(summary, out, indent=2, sort_keys=True)
        out.write("\n")

    print(f"Debug initial-structure PDBs written: {visible_pdb} and {all_pdb}")
    print(f"Debug single-particle lipid bilayer PDB written: {bilayer_pdb}")
    print(f"Debug bond report written: {bonds_tsv}")
    return summary


def write_initial_structure_debug_outputs(
    debug_dir: Path,
    prefix: str,
    pos,
    atom_types,
    atom_names,
    particle_classes,
    residue_ids,
    molecule_ids,
    charges,
    bond_ids=None,
    bond_r0=None,
    bond_k=None,
    box_lengths=None,
):
    """Write PDB/topology diagnostics from in-memory preparation arrays."""
    debug_dir = Path(debug_dir)
    debug_dir.mkdir(parents=True, exist_ok=True)
    pos = np.asarray(pos, dtype=np.float64)
    atom_types = np.asarray(atom_types, dtype=object)
    atom_names = np.asarray(atom_names, dtype=object)
    particle_classes = np.asarray(particle_classes, dtype=object)
    residue_ids = np.asarray(residue_ids, dtype=np.int32)
    molecule_ids = np.asarray(molecule_ids, dtype=np.int32)
    charges = np.asarray(charges, dtype=np.float64)
    bond_ids = np.zeros((0, 2), dtype=np.int32) if bond_ids is None else np.asarray(bond_ids, dtype=np.int32)
    bond_r0 = np.zeros(0, dtype=np.float64) if bond_r0 is None else np.asarray(bond_r0, dtype=np.float64)
    bond_k = np.zeros(0, dtype=np.float64) if bond_k is None else np.asarray(bond_k, dtype=np.float64)
    box_lengths = None if box_lengths is None else np.asarray(box_lengths, dtype=np.float64)

    all_pdb = debug_dir / f"{prefix}.input_debug.all.pdb"
    visible_pdb = debug_dir / f"{prefix}.input_debug.visible.pdb"
    bilayer_pdb = debug_dir / f"{prefix}.input_debug.bilayer.pdb"
    bonds_tsv = debug_dir / f"{prefix}.input_debug.bonds.tsv"
    summary_json = debug_dir / f"{prefix}.input_debug.summary.json"
    cgl_indices = np.where(np.asarray([str(x).strip().upper() == "CGL" for x in atom_types]))[0]

    write_pdb(
        all_pdb,
        _debug_records(pos, atom_types, atom_names, particle_classes, residue_ids, charges, True),
        box_lengths,
    )
    write_pdb(
        visible_pdb,
        _debug_records(pos, atom_types, atom_names, particle_classes, residue_ids, charges, False),
        box_lengths,
    )
    write_pdb(
        bilayer_pdb,
        _debug_records(
            pos,
            atom_types,
            atom_names,
            particle_classes,
            residue_ids,
            charges,
            False,
            include_indices=cgl_indices,
        ),
        box_lengths,
    )

    warnings_list = []
    ion_bond_count = 0
    cgl_cgld_bond_count = 0
    unexpected_bond_count = 0
    with bonds_tsv.open("w", encoding="utf-8") as out:
        out.write(
            "bond_index\ti\tj\tserial_i\tserial_j\tname_i\tname_j\ttype_i\ttype_j\t"
            "class_i\tclass_j\tresidue_i\tresidue_j\tmolecule_i\tmolecule_j\t"
            "distance_ang\tequil_dist_ang\tspring_const\n"
        )
        for b_idx, (i_raw, j_raw) in enumerate(bond_ids):
            i = int(i_raw)
            j = int(j_raw)
            if i < 0 or j < 0 or i >= pos.shape[0] or j >= pos.shape[0]:
                warnings_list.append(f"Bond {b_idx} has out-of-range atom ids: {i},{j}")
                continue
            distance = float(np.linalg.norm(_minimum_image_delta(pos[j] - pos[i], box_lengths)))
            r0 = float(bond_r0[b_idx]) if b_idx < bond_r0.shape[0] else float("nan")
            spring_k = float(bond_k[b_idx]) if b_idx < bond_k.shape[0] else float("nan")
            class_i = str(particle_classes[i])
            class_j = str(particle_classes[j])
            type_i = str(atom_types[i])
            type_j = str(atom_types[j])
            if class_i.upper() == "ION" or class_j.upper() == "ION":
                ion_bond_count += 1
            pair_types = {type_i.upper(), type_j.upper()}
            if pair_types == {"CGL", "CGLD"}:
                cgl_cgld_bond_count += 1
            else:
                unexpected_bond_count += 1
            out.write(
                f"{b_idx}\t{i}\t{j}\t{i + 1}\t{j + 1}\t{atom_names[i]}\t{atom_names[j]}\t"
                f"{type_i}\t{type_j}\t{class_i}\t{class_j}\t{int(residue_ids[i])}\t"
                f"{int(residue_ids[j])}\t{int(molecule_ids[i])}\t{int(molecule_ids[j])}\t"
                f"{distance:.6f}\t{r0:.6f}\t{spring_k:.6f}\n"
            )

    type_counts = {str(k): int(v) for k, v in Counter(str(x) for x in atom_types).items()}
    class_counts = {str(k): int(v) for k, v in Counter(str(x) for x in particle_classes).items()}
    cgl_count = int(type_counts.get("CGL", 0))
    cgld_count = int(type_counts.get("CGLD", 0))
    if ion_bond_count:
        warnings_list.append(f"ERROR: {ion_bond_count} true dist_spring bonds involve ions")
    if cgl_count != cgld_count:
        warnings_list.append(f"CGL/CGLD count mismatch: CGL={cgl_count}, CGLD={cgld_count}")
    if unexpected_bond_count:
        warnings_list.append(f"{unexpected_bond_count} dist_spring bonds are not CGL-CGLD orientation bonds")

    summary = {
        "outputs": {
            "all_pdb": str(all_pdb),
            "visible_pdb": str(visible_pdb),
            "bilayer_pdb": str(bilayer_pdb),
            "bonds_tsv": str(bonds_tsv),
            "summary_json": str(summary_json),
        },
        "n_atom": int(pos.shape[0]),
        "box_lengths": None if box_lengths is None else [float(x) for x in box_lengths],
        "type_counts": type_counts,
        "particle_class_counts": class_counts,
        "molecule_id_count": int(len(set(int(x) for x in molecule_ids))),
        "bond_count": int(bond_ids.shape[0]),
        "ion_bond_count": int(ion_bond_count),
        "cgl_count": cgl_count,
        "cgld_count": cgld_count,
        "cgl_cgld_bond_count": int(cgl_cgld_bond_count),
        "unexpected_bond_count": int(unexpected_bond_count),
        "leaflet_stats": _debug_leaflet_stats(pos, atom_types, box_lengths),
        "cgl_partner_stats": _debug_cgl_partner_stats(pos, atom_types, particle_classes, box_lengths),
        "warnings": warnings_list,
        "source": "in_memory_initial_structure",
    }
    with summary_json.open("w", encoding="utf-8") as out:
        json.dump(summary, out, indent=2, sort_keys=True)
        out.write("\n")
    print(f"Initial-structure debug PDBs written: {visible_pdb} and {all_pdb}")
    print(f"Initial single-particle lipid bilayer PDB written: {bilayer_pdb}")
    print(f"Initial-structure bond report written: {bonds_tsv}")
    return summary


def coords(atoms):
    return np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)


def set_coords(atoms, xyz):
    for atom, c in zip(atoms, xyz):
        atom["x"], atom["y"], atom["z"] = float(c[0]), float(c[1]), float(c[2])


def center_of_mass(xyz):
    return np.mean(xyz, axis=0)


def lipid_resname(resname: str) -> bool:
    return resname.upper() in {"DOP", "DOPC"}


def canonical_lipid_resname(resname: str) -> str:
    if lipid_resname(resname):
        return "DOPC"
    return resname


# DOPC 14-bead atom names in ITP order (indices 0-13)
_DOPC_ATOM_NAMES = [
    "NC3", "PO4", "GL1", "GL2",
    "C1A", "C2A", "D3A", "C4A", "C5A",
    "C1B", "C2B", "D3B", "C4B", "C5B",
]


def build_cg_lipid_array(initial_positions, atom_types, charges, residue_ids,
                          atom_names, residue_names, chain_ids, seg_ids):
    """Coarse-grain 14-bead DOPC molecules into single 6D vector particles.

    For each DOPC residue:
      - COM: geometric center of all 14 beads
      - Direction: normalized vector from head (NC3) to tail midpoint (C5A, C5B)

    Returns
    -------
    new_positions : (n_new, 3) float
    new_atom_types : (n_new,) S
    new_charges : (n_new,) float
    new_residue_ids : (n_new,) int
    new_atom_names : (n_new,) S
    new_residue_names : (n_new,) S
    new_chain_ids : (n_new,) S
    new_seg_ids : (n_new,) S
    lipid_directions : (n_lipid, 3) float  -- unit direction vectors
    cg_lipid_indices : (n_lipid,) int     -- indices of CG lipids in the new arrays
    lipid_to_atom_map : list[list[int]]    -- original atom index per lipid
    ref_bead_positions : (n_lipid, 14, 3)  -- bead positions relative to COM (Ångstrom)
    """
    n_atoms = len(initial_positions)
    is_dopc = np.array([
        residue_names[i].upper() in ("DOPC", "DOP")
        for i in range(n_atoms)
    ])

    if not np.any(is_dopc):
        n_a = len(initial_positions)
        return (initial_positions, atom_types, charges, residue_ids,
                atom_names, residue_names, chain_ids, seg_ids,
                np.empty((0, 3), dtype=float),
                np.empty(0, dtype=int),
                [],
                np.empty((0, 14, 3), dtype=float),
                np.arange(n_a, dtype=int))

    # Group DOPC atoms by residue
    dopc_indices = np.where(is_dopc)[0]
    lipid_groups = []
    current_group = [dopc_indices[0]]
    for idx in dopc_indices[1:]:
        if residue_ids[idx] == residue_ids[current_group[-1]]:
            current_group.append(idx)
        else:
            lipid_groups.append(current_group)
            current_group = [idx]
    lipid_groups.append(current_group)

    n_lipid = len(lipid_groups)
    print(f"Found {n_lipid} DOPC lipids ({sum(len(g) for g in lipid_groups)} beads)")

    # Build atom-name → position map for each lipid to get beads in ITP order
    lipid_directions = np.zeros((n_lipid, 3), dtype=float)
    lipid_to_atom_map = []
    ref_bead_positions = np.zeros((n_lipid, 14, 3), dtype=float)
    cg_positions = np.zeros((n_lipid, 3), dtype=float)

    for li, group in enumerate(lipid_groups):
        name_to_pos = {}
        for ai in group:
            name_to_pos[atom_names[ai].upper()] = initial_positions[ai]
        lipid_to_atom_map.append(group)

        # Get bead positions in ITP order
        bead_pos = np.zeros((14, 3), dtype=float)
        for bi, aname in enumerate(_DOPC_ATOM_NAMES):
            if aname in name_to_pos:
                bead_pos[bi] = name_to_pos[aname]
            else:
                raise ValueError(
                    f"DOPC lipid residue {residue_ids[group[0]]} missing atom '{aname}'"
                )

        # COM: geometric center
        com = np.mean(bead_pos, axis=0)
        cg_positions[li] = com

        # Reference bead positions relative to COM
        ref_bead_positions[li] = bead_pos - com

        # Direction: head (NC3) → tail midpoint (C5A, C5B)
        head_pos = bead_pos[0]  # NC3
        tail_mid = (bead_pos[8] + bead_pos[13]) / 2.0  # C5A + C5B
        direction = tail_mid - head_pos
        norm = np.linalg.norm(direction)
        if norm < 1e-8:
            # Fallback: use first principal component of bead positions
            centered = bead_pos - com
            _, _, vh = np.linalg.svd(centered, full_matrices=False)
            direction = vh[0]
            norm = np.linalg.norm(direction)
        lipid_directions[li] = direction / norm

    # Build new arrays: keep non-DOPC atoms, replace each DOPC group with 1 CG particle
    n_new = n_atoms - sum(len(g) for g in lipid_groups) + n_lipid
    old_to_new = np.full(n_atoms, -1, dtype=int)

    new_positions = np.zeros((n_new, 3), dtype=float)
    new_atom_types = np.empty(n_new, dtype='<U8')
    new_charges = np.zeros(n_new, dtype=float)
    new_residue_ids = np.zeros(n_new, dtype=int)
    new_atom_names = np.empty(n_new, dtype='<U8')
    new_residue_names = np.empty(n_new, dtype='<U8')
    new_chain_ids = np.empty(n_new, dtype='<U8')
    new_seg_ids = np.empty(n_new, dtype='<U8')

    cg_lipid_indices = np.zeros(n_lipid, dtype=int)
    write_pos = 0

    for li, group in enumerate(lipid_groups):
        # Copy non-DOPC atoms before this group
        prev_end = lipid_groups[li - 1][-1] + 1 if li > 0 else 0
        group_start = group[0]
        if group_start > prev_end:
            count = group_start - prev_end
            old_to_new[prev_end:group_start] = np.arange(write_pos, write_pos + count)
            new_positions[write_pos:write_pos + count] = initial_positions[prev_end:group_start]
            new_atom_types[write_pos:write_pos + count] = atom_types[prev_end:group_start]
            new_charges[write_pos:write_pos + count] = charges[prev_end:group_start]
            new_residue_ids[write_pos:write_pos + count] = residue_ids[prev_end:group_start]
            new_atom_names[write_pos:write_pos + count] = atom_names[prev_end:group_start]
            new_residue_names[write_pos:write_pos + count] = residue_names[prev_end:group_start]
            new_chain_ids[write_pos:write_pos + count] = chain_ids[prev_end:group_start]
            new_seg_ids[write_pos:write_pos + count] = seg_ids[prev_end:group_start]
            write_pos += count

        # Insert CG lipid particle
        cg_lipid_indices[li] = write_pos
        for ai in group:
            old_to_new[ai] = write_pos
        new_positions[write_pos] = cg_positions[li]
        new_atom_types[write_pos] = "CGL"
        new_charges[write_pos] = 0.0
        new_residue_ids[write_pos] = residue_ids[group[0]]
        new_atom_names[write_pos] = "CGL"
        new_residue_names[write_pos] = "DOPC"
        new_chain_ids[write_pos] = chain_ids[group[0]]
        new_seg_ids[write_pos] = seg_ids[group[0]]
        write_pos += 1

    # Copy remaining non-DOPC atoms after the last lipid group
    last_end = lipid_groups[-1][-1] + 1
    if last_end < n_atoms:
        count = n_atoms - last_end
        old_to_new[last_end:n_atoms] = np.arange(write_pos, write_pos + count)
        new_positions[write_pos:write_pos + count] = initial_positions[last_end:n_atoms]
        new_atom_types[write_pos:write_pos + count] = atom_types[last_end:n_atoms]
        new_charges[write_pos:write_pos + count] = charges[last_end:n_atoms]
        new_residue_ids[write_pos:write_pos + count] = residue_ids[last_end:n_atoms]
        new_atom_names[write_pos:write_pos + count] = atom_names[last_end:n_atoms]
        new_residue_names[write_pos:write_pos + count] = residue_names[last_end:n_atoms]
        new_chain_ids[write_pos:write_pos + count] = chain_ids[last_end:n_atoms]
        new_seg_ids[write_pos:write_pos + count] = seg_ids[last_end:n_atoms]
        write_pos += count

    print(f"CG lipid coarse-graining: {n_atoms} → {n_new} atoms "
          f"({n_lipid} CG lipids, {n_atoms - sum(len(g) for g in lipid_groups)} non-lipid)")

    return (new_positions, new_atom_types, new_charges, new_residue_ids,
            new_atom_names, new_residue_names, new_chain_ids, new_seg_ids,
            lipid_directions, cg_lipid_indices, lipid_to_atom_map, ref_bead_positions,
            old_to_new)


def infer_protein_charge_from_residues(protein_atoms):
    charged_res = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}
    seen = set()
    total = 0
    for atom in protein_atoms:
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key in seen:
            continue
        seen.add(key)
        total += charged_res.get(atom["resname"].upper(), 0)
    return total


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


def _collect_complete_backbone_residue_order(backbone_atoms):
    residue_groups = defaultdict(dict)
    residue_order = []
    for atom in backbone_atoms:
        key = residue_key(atom)
        role = atom["name"].strip().upper()
        if key not in residue_groups:
            residue_order.append(key)
        residue_groups[key][role] = 1

    valid_order = []
    for key in residue_order:
        group = residue_groups[key]
        if not all(role in group for role in BB_COMPONENT_NAMES):
            continue
        valid_order.append(key)
    if not valid_order:
        raise ValueError("Cannot map BB types: no complete N/CA/C/O residue rows.")
    return valid_order


def map_backbone_types_from_martinize_fallback(backbone_atoms):
    # Use martinize fallback behavior with secondary structure fixed to coil ("C"):
    # default/backbone override table + charged termini/chain-break endpoints.
    residue_order = _collect_complete_backbone_residue_order(backbone_atoms)
    coil_column = 8
    bb_default = ("N0", "Nda", "N0", "Nd", "Na", "Nda", "Nda", "P5", "P5")
    bb_residue_override = {
        "ALA": ("C5", "N0", "C5", "N0", "N0", "N0", "N0", "P4", "P4"),
        "PRO": ("C5", "N0", "C5", "N0", "Na", "N0", "N0", "Na", "Na"),
        "HYP": ("C5", "N0", "C5", "N0", "N0", "N0", "N0", "Na", "Na"),
    }
    out = {}
    for key in residue_order:
        resname = key[3].upper()
        table = bb_residue_override.get(resname, bb_default)
        out[key] = table[coil_column]

    # Match martinize charged-termini behavior: each fragment endpoint is
    # assigned charged backbone bead types (Qd/Qa) in residue order.
    if residue_order:
        out[residue_order[0]] = "Qd"
        out[residue_order[-1]] = "Qa"
        for i in range(1, len(residue_order)):
            prev = residue_order[i - 1]
            curr = residue_order[i]
            chain_break = curr[0] != prev[0]
            seq_break = curr[1] != (prev[1] + 1)
            if chain_break or seq_break:
                out[curr] = "Qd"
                out[prev] = "Qa"
    return out


def compute_lipid_residue_indices(bilayer_atoms):
    lipid_residues = defaultdict(list)
    keep_nonlipid = []
    for idx, atom in enumerate(bilayer_atoms):
        if lipid_resname(atom["resname"]):
            key = (atom["chain"], atom["resseq"], atom["icode"])
            lipid_residues[key].append(idx)
        else:
            keep_nonlipid.append(idx)
    return lipid_residues, keep_nonlipid


def residue_group_atoms(atoms):
    groups = defaultdict(list)
    for atom in atoms:
        key = (atom["chain"], atom["resseq"], atom["icode"], atom["resname"].upper())
        groups[key].append(atom)
    return groups


def tile_and_crop_bilayer_lipids(
    bilayer_atoms,
    bilayer_box,
    target_xy_min,
    target_xy_max,
):
    # Keep only lipid residues for tiling; ions are re-generated after packing.
    lipid_atoms = [a for a in bilayer_atoms if lipid_resname(a["resname"])]
    if not lipid_atoms:
        raise ValueError("No lipid residues found in bilayer input for tiling/cropping.")

    lip_xyz = coords(lipid_atoms)
    lip_center = center_of_mass(lip_xyz)

    if bilayer_box is not None:
        tile_x, tile_y = float(bilayer_box[0]), float(bilayer_box[1])
    else:
        lip_min = lip_xyz.min(axis=0)
        lip_max = lip_xyz.max(axis=0)
        span = lip_max - lip_min
        tile_x, tile_y = float(span[0]), float(span[1])
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
                    a["resname"] = canonical_lipid_resname(a["resname"])
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
            atom["resname"] = canonical_lipid_resname(atom["resname"])
            atom["segid"] = "MEMB"
            out_atoms.append(atom)
        next_resseq += 1
    return out_atoms


def remove_overlapping_lipids(
    bilayer_atoms,
    protein_atoms,
    lipid_residues,
    keep_nonlipid,
    cutoff,
):
    protein_xyz = coords(protein_atoms)
    cutoff2 = cutoff * cutoff
    keep_atom_idx = set(keep_nonlipid)
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


def resize_and_shift(all_atoms, pad_xy, pad_z):
    xyz = coords(all_atoms)
    min_c = xyz.min(axis=0)
    max_c = xyz.max(axis=0)
    spans = max_c - min_c
    box = np.array(
        [spans[0] + 2.0 * pad_xy, spans[1] + 2.0 * pad_xy, spans[2] + 2.0 * pad_z],
        dtype=float,
    )
    shift = np.array([pad_xy, pad_xy, pad_z], dtype=float) - min_c
    shifted = xyz + shift
    set_coords(all_atoms, shifted)
    return box


def set_box_from_lipid_xy(
    all_atoms,
    lipid_atoms,
    pad_z,
    force_square_xy=True,
    min_box_z=None,
    center_lipid_in_z=True,
):
    all_xyz = coords(all_atoms)
    lip_xyz = coords(lipid_atoms)
    lip_min = lip_xyz.min(axis=0)
    lip_max = lip_xyz.max(axis=0)

    # Force XY box edges to coincide with bilayer edges.
    span_x = float(lip_max[0] - lip_min[0])
    span_y = float(lip_max[1] - lip_min[1])
    box_x = span_x
    box_y = span_y
    if force_square_xy:
        square_side = max(span_x, span_y)
        box_x = square_side
        box_y = square_side
    if span_x <= 0.0 or span_y <= 0.0:
        raise ValueError("Invalid lipid XY span while defining box edges.")

    min_z = float(all_xyz[:, 2].min())
    max_z = float(all_xyz[:, 2].max())
    lip_mid_z = 0.5 * float(lip_min[2] + lip_max[2])
    span_z = max_z - min_z
    box_z = float(span_z + 2.0 * pad_z)
    if min_box_z is not None:
        box_z = float(max(box_z, float(min_box_z)))
    if center_lipid_in_z:
        # Keep lipid bilayer centered around box_z/2 while ensuring all atoms remain inside [0, box_z].
        needed_z = 2.0 * max(lip_mid_z - min_z, max_z - lip_mid_z)
        box_z = float(max(box_z, needed_z))

    shift = np.array([0.0, 0.0, 0.0], dtype=float)
    shift[0] = -lip_min[0]
    shift[1] = -lip_min[1]
    if center_lipid_in_z:
        shift[2] = 0.5 * box_z - lip_mid_z
    else:
        shift[2] = pad_z - min_z

    shifted = all_xyz + shift
    set_coords(all_atoms, shifted)
    return np.array([box_x, box_y, box_z], dtype=float)


def estimate_salt_pairs(box_lengths, salt_molar, effective_volume_fraction=1.0):
    volume_a3 = float(box_lengths[0] * box_lengths[1] * box_lengths[2])
    volume_a3 *= float(max(0.0, min(1.0, effective_volume_fraction)))
    volume_l = volume_a3 * 1e-27
    pairs = int(round(salt_molar * NA_AVOGADRO * volume_l))
    return max(0, pairs)


def infer_effective_ion_volume_fraction_from_template(
    bilayer_atoms,
    bilayer_box,
    salt_molar,
):
    # Calibrate effective ion-accessible volume from the template system ions.
    # This avoids overestimating ion count by using full geometric box volume.
    if bilayer_box is None or salt_molar <= 0.0:
        return 1.0

    n_na = sum(1 for a in bilayer_atoms if a["resname"].upper() == "NA")
    n_cl = sum(1 for a in bilayer_atoms if a["resname"].upper() == "CL")
    base_pairs = min(n_na, n_cl)
    if base_pairs <= 0:
        return 1.0

    base_volume_l = float(bilayer_box[0] * bilayer_box[1] * bilayer_box[2]) * 1e-27
    expected_pairs = salt_molar * NA_AVOGADRO * base_volume_l
    if expected_pairs <= 0.0:
        return 1.0

    frac = float(base_pairs) / float(expected_pairs)
    return float(max(0.0, min(1.0, frac)))


def place_ions(atoms, box_lengths, n_na, n_cl, cutoff, rng):
    existing = coords(atoms)
    placed = []
    cutoff2 = cutoff * cutoff
    types = [("NA", "NA", 1.0, "IONS")] * n_na + [("CL", "CL", -1.0, "IONS")] * n_cl
    for name, resname, _charge, segid in types:
        accepted = False
        for _ in range(20000):
            trial = rng.uniform([0, 0, 0], box_lengths)
            if existing.size:
                d2 = np.sum((existing - trial) ** 2, axis=1)
                if np.min(d2) < cutoff2:
                    continue
            if placed:
                placed_xyz = np.array([[a["x"], a["y"], a["z"]] for a in placed], dtype=float)
                d2_placed = np.sum((placed_xyz - trial) ** 2, axis=1)
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


def collect_backbone_only_bb_map(backbone_atoms, bb_type_by_residue):
    residue_groups = defaultdict(dict)
    for atom_idx, atom in enumerate(backbone_atoms):
        key = residue_key(atom)
        role = atom["name"].strip().upper()
        residue_groups[key][role] = atom_idx

    ordered_keys = sorted(residue_groups.keys(), key=lambda x: (x[0], x[1], x[2], x[3]))
    bb_entries = []
    for seq_idx, key in enumerate(ordered_keys, start=1):
        role_map = residue_groups[key]
        if not all(role in role_map for role in BB_COMPONENT_NAMES):
            continue
        idxs = [int(role_map[role]) for role in BB_COMPONENT_NAMES]
        masses = np.asarray(BB_COMPONENT_MASSES, dtype=np.float32)
        weights = (masses / masses.sum()).tolist()
        coords_row = []
        for role in BB_COMPONENT_NAMES:
            atom = backbone_atoms[role_map[role]]
            coords_row.append([float(atom["x"]), float(atom["y"]), float(atom["z"])])
        bb_type = str(bb_type_by_residue.get(key, "P5")).strip()
        bb_entries.append(
            {
                "bb_residue_index": seq_idx,
                "bb_resseq": int(key[1]),
                "bb_chain": str(key[0]),
                "bb_icode": str(key[2]),
                "bb_atom_index": -1,
                "bb_type": bb_type,
                "atom_indices": idxs,
                "atom_mask": [1, 1, 1, 1],
                "weights": weights,
                "reference_atom_indices": idxs,
                "reference_atom_coords": coords_row,
                "bb_comment": (
                    f"Backbone-only BB map resid={seq_idx} resseq={key[1]} chain={key[0]} "
                    f"type={bb_type} atom_indices={idxs}"
                ),
            }
        )
    if not bb_entries:
        raise ValueError("No complete backbone mapping entries could be generated.")
    return bb_entries


def build_backbone_with_virtual_bb(backbone_atoms, bb_type_by_residue):
    residue_groups = defaultdict(dict)
    for atom_idx, atom in enumerate(backbone_atoms):
        key = residue_key(atom)
        role = atom["name"].strip().upper()
        residue_groups[key][role] = atom_idx

    ordered_keys = sorted(residue_groups.keys(), key=lambda x: (x[0], x[1], x[2], x[3]))
    protein_atoms = [deepcopy(atom) for atom in backbone_atoms]
    bb_entries = []
    for seq_idx, key in enumerate(ordered_keys, start=1):
        role_map = residue_groups[key]
        if not all(role in role_map for role in BB_COMPONENT_NAMES):
            continue
        idxs = [int(role_map[role]) for role in BB_COMPONENT_NAMES]
        masses = np.asarray(BB_COMPONENT_MASSES, dtype=np.float32)
        weights = (masses / masses.sum()).tolist()
        coords_row = []
        com = np.zeros(3, dtype=np.float64)
        wsum = 0.0
        for role, weight in zip(BB_COMPONENT_NAMES, weights):
            atom = backbone_atoms[role_map[role]]
            xyz = np.array([float(atom["x"]), float(atom["y"]), float(atom["z"])], dtype=np.float64)
            coords_row.append([float(xyz[0]), float(xyz[1]), float(xyz[2])])
            com += float(weight) * xyz
            wsum += float(weight)
        if wsum <= 0.0:
            raise ValueError("Backbone BB COM weights sum to zero.")
        if abs(wsum - 1.0) > 1.0e-8:
            com /= wsum

        bb_type = str(bb_type_by_residue.get(key, "P5")).strip()
        proxy = deepcopy(backbone_atoms[role_map["CA"]])
        proxy["name"] = "BB"
        proxy["x"] = float(com[0])
        proxy["y"] = float(com[1])
        proxy["z"] = float(com[2])
        proxy["segid"] = "PROA"
        proxy["element"] = "C"
        bb_atom_index = len(protein_atoms)
        protein_atoms.append(proxy)

        bb_entries.append(
            {
                "bb_residue_index": seq_idx,
                "bb_resseq": int(key[1]),
                "bb_chain": str(key[0]),
                "bb_icode": str(key[2]),
                "bb_atom_index": bb_atom_index,
                "bb_type": bb_type,
                "atom_indices": idxs,
                "atom_mask": [1, 1, 1, 1],
                "weights": weights,
                "reference_atom_indices": idxs,
                "reference_atom_coords": coords_row,
                "bb_comment": (
                    f"Backbone virtual BB map resid={seq_idx} resseq={key[1]} chain={key[0]} "
                    f"type={bb_type} bb_atom_index={bb_atom_index} atom_indices={idxs}"
                ),
            }
        )
    if not bb_entries:
        raise ValueError("No complete backbone mapping entries could be generated.")
    return protein_atoms, bb_entries


def derive_chain_break_metadata(bb_entries):
    if not bb_entries:
        return np.array([], dtype=np.int32), np.array([], dtype=np.int32)

    chain_ids = [str(entry.get("bb_chain", "")).strip() or " " for entry in bb_entries]
    chain_first_residue = [
        residue_index
        for residue_index in range(1, len(chain_ids))
        if chain_ids[residue_index] != chain_ids[residue_index - 1]
    ]
    n_chains = len(chain_first_residue) + 1
    return np.asarray(chain_first_residue, dtype=np.int32), np.ones(n_chains, dtype=np.int32)


def write_hybrid_mapping_h5(
    path: Path,
    bb_entries,
    total_martini_atoms,
    env_atom_indices,
    n_protein_atoms,
):
    with h5py.File(path, "w") as h5:
        inp = h5.create_group("input")

        ctrl = inp.create_group("hybrid_control")
        ctrl.attrs["enable"] = np.int8(1)
        ctrl.attrs["activation_stage"] = b"minimization"
        ctrl.attrs["preprod_protein_mode"] = b"rigid_body"
        ctrl.attrs["preprod_lipid_headgroup_roles"] = b"PO4"
        ctrl.attrs["exclude_intra_protein_martini"] = np.int8(1)
        ctrl.attrs["production_nonprotein_hard_sphere"] = np.int8(0)
        ctrl.attrs["protein_env_interface_scale"] = np.float32(1.0)
        ctrl.attrs["virtual_backbone_com_mode"] = np.int8(1)
        ctrl.attrs["schema_version"] = np.int32(1)

        bb_grp = inp.create_group("hybrid_bb_map")
        bb_grp.attrs["atom_index_space"] = b"protein_aa_pdb_0based"
        bb_grp.attrs["runtime_index_space"] = b"stage_runtime_after_injection"
        bb_grp.create_dataset(
            "bb_residue_index",
            data=np.array([b.get("bb_residue_index", b["bb_resseq"]) for b in bb_entries], dtype=np.int32),
        )
        bb_grp.create_dataset(
            "bb_resseq",
            data=np.array([b["bb_resseq"] for b in bb_entries], dtype=np.int32),
        )
        bb_grp.create_dataset(
            "bb_chain_id",
            data=np.array(
                [str(b.get("bb_chain", "")).strip() or " " for b in bb_entries],
                dtype=h5py.string_dtype(encoding="utf-8"),
            ),
        )
        bb_grp.create_dataset(
            "bb_icode",
            data=np.array(
                [str(b.get("bb_icode", "")).strip() or " " for b in bb_entries],
                dtype=h5py.string_dtype(encoding="utf-8"),
            ),
        )
        bb_grp.create_dataset(
            "bb_atom_index",
            data=np.array([b["bb_atom_index"] for b in bb_entries], dtype=np.int32),
        )
        bb_grp.create_dataset(
            "bb_type",
            data=np.array(
                [str(b.get("bb_type", "P5")).strip() for b in bb_entries],
                dtype=h5py.string_dtype(encoding="utf-8"),
            ),
        )
        bb_grp.create_dataset(
            "atom_indices",
            data=np.array([b["atom_indices"] for b in bb_entries], dtype=np.int32),
        )
        bb_grp.create_dataset(
            "atom_mask",
            data=np.array([b["atom_mask"] for b in bb_entries], dtype=np.int8),
        )
        bb_grp.create_dataset(
            "weights",
            data=np.array([b["weights"] for b in bb_entries], dtype=np.float32),
        )
        bb_grp.create_dataset(
            "protein_id",
            data=np.zeros(len(bb_entries), dtype=np.int32),
        )
        bb_grp.attrs["virtual_backbone_com_mode"] = np.int8(1)
        # Optional reference backbone metadata for auditability in .up/.h5 artifacts.
        bb_grp.create_dataset(
            "reference_atom_names",
            data=np.array(["N", "CA", "C", "O"], dtype=h5py.string_dtype(encoding="utf-8")),
        )
        bb_grp.create_dataset(
            "reference_atom_indices",
            data=np.array([b["reference_atom_indices"] for b in bb_entries], dtype=np.int32).reshape(
                len(bb_entries), 4
            ),
        )
        bb_grp.create_dataset(
            "reference_atom_coords",
            data=np.array([b["reference_atom_coords"] for b in bb_entries], dtype=np.float32).reshape(
                len(bb_entries), 4, 3
            ),
        )
        bb_grp.create_dataset(
            "bb_comment",
            data=np.array([b["bb_comment"] for b in bb_entries], dtype=h5py.string_dtype(encoding="utf-8")),
        )

        env_grp = inp.create_group("hybrid_env_topology")
        env_grp.create_dataset(
            "env_atom_indices",
            data=np.array(env_atom_indices, dtype=np.int32),
        )
        membership = np.full(total_martini_atoms, -1, dtype=np.int32)
        membership[:n_protein_atoms] = 0
        env_grp.create_dataset(
            "protein_membership",
            data=membership,
        )

        chain_first_residue, chain_counts = derive_chain_break_metadata(bb_entries)
        if chain_first_residue.size:
            break_grp = inp.create_group("chain_break")
            break_grp.create_dataset("chain_first_residue", data=chain_first_residue, dtype=np.int32)
            break_grp.create_dataset("chain_counts", data=chain_counts, dtype=np.int32)


# -----------------------------------------------------------------------------
# Stage Conversion Helpers
# -----------------------------------------------------------------------------
def runtime_input_pdb_path(script_dir, pdb_id):
    """Resolve runtime protein+bilayer PDB path with optional override."""
    override = os.environ.get("UPSIDE_RUNTIME_PDB_FILE", "").strip()
    if override:
        return os.path.abspath(os.path.expanduser(override))
    return os.path.join(script_dir, f"pdb/{pdb_id}.pdb")


def read_martini3_nonbond_params(itp_file):
    """
    Read MARTINI nonbonded parameters from the main .itp file
    Returns a dictionary mapping (type1, type2) tuples to (sigma, epsilon) values
    """
    martini_table = {}
    
    if not os.path.exists(itp_file):
        print(f"Error: MARTINI parameter file '{itp_file}' not found!")
        return martini_table
    
    with open(itp_file, 'r') as f:
        lines = f.readlines()
    
    # Parse #define macros for dry MARTINI style nonbond_params entries
    macro_params = {}
    for raw in lines:
        line = raw.split(';', 1)[0].strip()
        if not line.startswith('#define'):
            continue
        parts = line.split()
        # Format: #define m_VI 0.47 2.700
        if len(parts) >= 4:
            macro_name = parts[1]
            try:
                sigma = float(parts[2])
                epsilon = float(parts[3])
                macro_params[macro_name] = (sigma, epsilon)
            except ValueError:
                continue

    # Find the nonbond_params section
    in_nonbond_params = False
    param_count = 0
    
    for i, line in enumerate(lines):
        line = line.split(';', 1)[0].strip()
        
        # Check for nonbond_params section start
        if line == '[ nonbond_params ]' or line == '[nonbond_params]':
            in_nonbond_params = True
            continue
        elif line.startswith('[') and line.endswith(']') and line != '[ nonbond_params ]' and line != '[nonbond_params]':
            if in_nonbond_params:
                break
            in_nonbond_params = False
            continue
        
        # Parse nonbond_params lines
        if in_nonbond_params and line and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 5:
                type1 = parts[0]
                type2 = parts[1]
                try:
                    _func = int(parts[2])
                except ValueError:
                    continue
                try:
                    sigma = float(parts[3])  # nm
                    epsilon = float(parts[4])  # kJ/mol
                except ValueError:
                    # Dry MARTINI often uses a macro token in column 4
                    macro = parts[3]
                    if macro not in macro_params:
                        continue
                    sigma, epsilon = macro_params[macro]

                # Store both orientations for easy lookup
                martini_table[(type1, type2)] = (sigma, epsilon)
                martini_table[(type2, type1)] = (sigma, epsilon)
                param_count += 1
            elif len(parts) >= 4:
                # Alternative dry MARTINI format: type1 type2 func macro
                type1 = parts[0]
                type2 = parts[1]
                try:
                    _func = int(parts[2])
                except ValueError:
                    continue
                macro = parts[3]
                if macro not in macro_params:
                    continue
                sigma, epsilon = macro_params[macro]
                martini_table[(type1, type2)] = (sigma, epsilon)
                martini_table[(type2, type1)] = (sigma, epsilon)
                param_count += 1
    
    print(f"Read {len(martini_table)//2} unique nonbonded parameter pairs")
    return martini_table

# --- Helpers: read MARTINI protein topology (ITP) for protein bead typing and connectivity ---
def parse_itp_file(itp_file, target_molecule=None, preprocessor_defines=None):
    """
    Universal ITP parser to read MARTINI topology files.
    Returns a dictionary with parsed sections: atoms, bonds, angles, dihedrals, etc.
    If target_molecule is specified, only returns data for that specific molecule type.
    """
    topology = {
        'atoms': [],
        'bonds': [], 
        'angles': [],
        'dihedrals': [],
        'position_restraints': [],
        'exclusions': [],
        'moleculetype': None,
        'molecules': {}  # Store multiple molecule types
    }
    
    if not os.path.exists(itp_file):
        print(f"Warning: ITP file {itp_file} not found")
        return topology
    
    current_section = None
    current_molecule = None
    current_mol_data = None
    macro_defs = {}
    if preprocessor_defines:
        for macro_name, macro_value in preprocessor_defines.items():
            if isinstance(macro_value, (list, tuple)):
                macro_defs[macro_name] = [float(x) for x in macro_value]
            else:
                macro_defs[macro_name] = [float(macro_value)]
    pp_stack = []
    current_active = True

    def parse_macro_value(tok, macro_index=None):
        # Try direct float
        try:
            return float(tok)
        except ValueError:
            pass
        # Try macro lookup
        values = macro_defs.get(tok)
        if values is None:
            raise ValueError(f"Unknown macro '{tok}'")
        if macro_index is None:
            if not values:
                raise ValueError(f"Macro '{tok}' has no values")
            return values[0]
        if macro_index < len(values):
            return values[macro_index]
        if values:
            return values[-1]
        raise ValueError(f"Macro '{tok}' has no values")
    
    with open(itp_file, 'r') as f:
        for line in f:
            line = line.split(';', 1)[0].strip()
            
            # Skip empty lines and comments
            if not line:
                continue

            # Handle a minimal preprocessor subset so dry MARTINI #ifndef/#else blocks
            # resolve to a single active topology branch.
            if line.startswith('#ifdef') or line.startswith('#ifndef'):
                parts = line.split()
                macro_name = parts[1] if len(parts) >= 2 else ""
                is_defined = macro_name in macro_defs
                cond = is_defined if line.startswith('#ifdef') else (not is_defined)
                pp_stack.append((current_active, cond))
                current_active = current_active and cond
                continue
            if line.startswith('#else'):
                if pp_stack:
                    parent_active, cond = pp_stack[-1]
                    cond = not cond
                    pp_stack[-1] = (parent_active, cond)
                    current_active = parent_active and cond
                continue
            if line.startswith('#endif'):
                if pp_stack:
                    parent_active, _cond = pp_stack.pop()
                    current_active = parent_active
                continue
            if not current_active:
                continue

            # Parse and store macro definitions used in dry MARTINI lipids
            if line.startswith('#define'):
                parts = line.split()
                if len(parts) >= 3:
                    macro_name = parts[1]
                    macro_vals = []
                    for tok in parts[2:]:
                        try:
                            macro_vals.append(float(tok))
                        except ValueError:
                            break
                    if macro_vals:
                        macro_defs[macro_name] = macro_vals
                continue

            # Ignore preprocessor directives like #ifdef/#ifndef/#endif/#include
            if line.startswith('#'):
                continue
            
            # Check for section headers
            if line.startswith('[') and line.endswith(']'):
                section_name = line[1:-1].strip().lower()
                current_section = section_name
                continue
            
            # Parse based on current section
            if current_section == 'moleculetype':
                parts = line.split()
                if len(parts) >= 1:
                    current_molecule = parts[0]
                    if topology['moleculetype'] is None:
                        topology['moleculetype'] = current_molecule
                    
                    # Initialize data structure for this molecule
                    current_mol_data = {
                        'atoms': [],
                        'bonds': [],
                        'angles': [],
                        'dihedrals': [],
                        'position_restraints': [],
                        'exclusions': []
                    }
                    topology['molecules'][current_molecule] = current_mol_data
            
            elif current_section == 'atoms' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        atom_data = {
                            'id': int(parts[0]),
                            'type': parts[1],
                            'resnr': int(parts[2]),
                            'residue': parts[3],
                            'atom': parts[4],
                            'cgnr': int(parts[5]),
                            'charge': 0.0,
                            'mass': 0.0
                        }
                        # Parse charge and mass more carefully
                        # Charge is typically in column 6 (0-indexed), mass in column 7
                        if len(parts) >= 7:
                            try:
                                atom_data['charge'] = float(parts[6])
                            except ValueError:
                                pass
                        if len(parts) >= 8:
                            try:
                                atom_data['mass'] = float(parts[7])
                            except ValueError:
                                pass
                        current_mol_data['atoms'].append(atom_data)
                        topology['atoms'].append(atom_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'bonds' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        bond_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'j': int(parts[1]) - 1,
                            'func': int(parts[2]),
                            'r0': 0.0,
                            'k': 0.0
                        }
                        if len(parts) >= 5:
                            bond_data['r0'] = parse_macro_value(parts[3], macro_index=0)  # nm
                            bond_data['k'] = parse_macro_value(parts[4], macro_index=1)   # kJ/mol/nm²
                        elif len(parts) >= 4:
                            # Dry MARTINI alias form: i j func macro
                            bond_data['r0'] = parse_macro_value(parts[3], macro_index=0)
                            bond_data['k'] = parse_macro_value(parts[3], macro_index=1)
                        current_mol_data['bonds'].append(bond_data)
                        topology['bonds'].append(bond_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'angles' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        angle_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'j': int(parts[1]) - 1,
                            'k': int(parts[2]) - 1,
                            'func': int(parts[3]),
                            'theta0': 0.0,  # degrees
                            'force_k': 0.0  # kJ/mol/rad²
                        }
                        if len(parts) >= 6:
                            angle_data['theta0'] = parse_macro_value(parts[4], macro_index=0)
                            angle_data['force_k'] = parse_macro_value(parts[5], macro_index=1)
                        else:
                            # Dry MARTINI alias form: i j k func macro
                            angle_data['theta0'] = parse_macro_value(parts[4], macro_index=0)
                            angle_data['force_k'] = parse_macro_value(parts[4], macro_index=1)
                        current_mol_data['angles'].append(angle_data)
                        topology['angles'].append(angle_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'dihedrals' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        dihedral_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'j': int(parts[1]) - 1,
                            'k': int(parts[2]) - 1,
                            'l': int(parts[3]) - 1,
                            'func': int(parts[4]),
                            'phi0': 0.0,
                            'k': 0.0,
                            'mult': 1
                        }
                        if len(parts) >= 7:
                            dihedral_data['phi0'] = float(parts[5])  # degrees
                            dihedral_data['k'] = float(parts[6])     # kJ/mol
                        if len(parts) >= 8:
                            dihedral_data['mult'] = int(parts[7])
                        current_mol_data['dihedrals'].append(dihedral_data)
                        topology['dihedrals'].append(dihedral_data)
                    except (ValueError, IndexError):
                        continue

            elif current_section == 'position_restraints' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        restraint_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'func': int(parts[1]),
                            'fx': parse_macro_value(parts[2]),
                            'fy': parse_macro_value(parts[3]),
                            'fz': parse_macro_value(parts[4]),
                        }
                        current_mol_data['position_restraints'].append(restraint_data)
                        topology['position_restraints'].append(restraint_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'exclusions' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        # Convert to 0-indexed and add all pairs
                        atoms = [int(part)-1 for part in parts]
                        for i in range(len(atoms)):
                            for j in range(i+1, len(atoms)):
                                exclusion = (atoms[i], atoms[j])
                                current_mol_data['exclusions'].append(exclusion)
                                topology['exclusions'].append(exclusion)
                    except ValueError:
                        continue
    
    # If target_molecule is specified, return only that molecule's data
    if target_molecule and target_molecule in topology['molecules']:
        mol_data = topology['molecules'][target_molecule]
        return {
            'atoms': mol_data['atoms'],
            'bonds': mol_data['bonds'],
            'angles': mol_data['angles'],
            'dihedrals': mol_data['dihedrals'],
            'position_restraints': mol_data['position_restraints'],
            'exclusions': mol_data['exclusions'],
            'moleculetype': target_molecule,
            'molecules': {target_molecule: mol_data}
        }
    
    return topology

def read_martini_masses(ff_file):
    """Read atom type masses from MARTINI force field file"""
    masses = {}
    ff_file_path = Path(ff_file).expanduser()
    if not ff_file_path.is_absolute():
        ff_file_path = (WORKFLOW_DIR / ff_file_path).resolve()
    else:
        ff_file_path = ff_file_path.resolve()
    
    if not ff_file_path.exists():
        raise ValueError(f"FATAL ERROR: Force field file '{ff_file}' not found.\n"
                        f"  Full path: {ff_file_path}\n"
                        f"  This file is required for atom type masses.\n"
                        f"  Please ensure the MARTINI force field file exists and is readable.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    in_atomtypes = False
    with ff_file_path.open('r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('[ atomtypes ]'):
                in_atomtypes = True
                continue
            elif line.startswith('['):
                in_atomtypes = False
                continue
            
            if in_atomtypes and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 2:
                    atom_type = parts[0]
                    try:
                        mass = float(parts[1])
                        masses[atom_type] = mass
                    except ValueError:
                        continue
    
    return masses

def convert_stage(pdb_id=None, stage='minimization', run_dir=None):
    """
    Main preparation function with stage-specific parameterization
    
    Args:
        stage: Simulation stage ('minimization', 'npt_equil', 'npt_equil_reduced', 'npt_prod')
        run_dir: Optional run directory (default: outputs/martini_test)
    """

    # Get UPSIDE home directory
    upside_path = os.environ['UPSIDE_HOME']

    # Resolve PDB ID from explicit argument first, then legacy CLI fallback.
    if pdb_id is None:
        if len(sys.argv) > 1:
            pdb_id = sys.argv[1]
        else:
            raise ValueError("FATAL ERROR: No PDB ID provided.\n"
                            f"  Usage: python {sys.argv[0]} <pdb_id>\n"
                            f"  Example: python {sys.argv[0]} 1rkl\n"
                            f"  Aborting to prevent incorrect simulation results.")

    # Get run directory from function parameter or use default
    if run_dir is None:
        if len(sys.argv) > 2 and sys.argv[2] != '--stage':
            run_dir = sys.argv[2]
        else:
            run_dir = "outputs/martini_test"
    os.makedirs(run_dir, exist_ok=True)

    # Get stage from environment variable or use default
    stage = os.environ.get('UPSIDE_SIMULATION_STAGE', stage)
    print(f"Preparing for stage: {stage}")

    # Stage-specific parameterization
    stage_params = {
        'minimization': {
            'lj_soften': 1,
            'lj_alpha': 0.2,
            'coulomb_soften': 1,
            'slater_alpha': 2.0,
            'barostat_type': 0  # Berendsen
        },
        'npt_equil': {
            'lj_soften': 1,
            'lj_alpha': 0.2,
            'coulomb_soften': 1,
            'slater_alpha': 2.0,
            'barostat_type': 0  # Berendsen
        },
        'npt_equil_reduced': {
            'lj_soften': 1,
            'lj_alpha': 0.05,
            'coulomb_soften': 1,
            'slater_alpha': 0.5,
            'barostat_type': 0  # Berendsen
        },
        'npt_prod': {
            'lj_soften': 0,
            'lj_alpha': 0.0,
            'coulomb_soften': 0,
            'slater_alpha': 0.0,
            'barostat_type': 1  # Parrinello-Rahman
        }
    }

    # Get parameters for current stage
    params = stage_params.get(stage, stage_params['npt_prod'])
    stage_lipidhead_fc = float(os.environ.get('UPSIDE_BILAYER_LIPIDHEAD_FC', '0'))


    # Configuration
    strict_from_martini_pdb = True
    include_protein = True
    
    print("=== Dry MARTINI Protein-Lipid System Preparation ===")
    print(f"PDB ID: {pdb_id}")
    print(f"Output directory: {run_dir}")
    
    workflow_dir = str(WORKFLOW_DIR)
    # Read dry MARTINI parameter files
    print("\n=== Reading Dry MARTINI Parameters ===")
    ff_dir = Path(
        os.environ.get('UPSIDE_MARTINI_FF_DIR', str(REPO_ROOT / "parameters" / "dryMARTINI"))
    ).expanduser()
    if not ff_dir.is_absolute():
        ff_dir = (REPO_ROOT / ff_dir).resolve()
    else:
        ff_dir = ff_dir.resolve()
    ff_path = str(ff_dir)

    if not os.path.isdir(ff_path):
        raise ValueError(f"FATAL ERROR: Force-field directory '{ff_path}' not found.\n"
                        "  Set UPSIDE_MARTINI_FF_DIR to a valid dry-MARTINI force-field directory.")

    ff_files = sorted(os.listdir(ff_path))
    print(f"Using force-field directory: {ff_path}")

    def pick_ff_file(name, required=True):
        path = os.path.join(ff_path, name)
        if os.path.exists(path):
            return path
        if required:
            raise ValueError(f"FATAL ERROR: Required dry MARTINI file '{name}' not found in '{ff_path}'.\n"
                            f"  Available: {ff_files}")
        return None

    # Read nonbonded parameters
    martini_param_file = pick_ff_file("dry_martini_v2.1.itp")
    martini_table = read_martini3_nonbond_params(martini_param_file)
    
    if not martini_table:
        raise ValueError(f"FATAL ERROR: Could not read MARTINI parameters from '{martini_param_file}'\n"
                        f"  This file is required for proper force field parameterization.\n"
                        f"  Please ensure the dry MARTINI parameter file exists and is readable.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    print("=== Mixed Protein-Lipid System Detected ===")
    print("Protein runtime representation: AA backbone carriers only (N/CA/C/O)")
    
    # For both protein and lipid systems, we need DOPC parameters
    # Parse DOPC topology from ITP file
    dopc_param_file = pick_ff_file("dry_martini_v2.1_lipids.itp")
    lipid_preproc_defs = {}
    if stage_lipidhead_fc > 0.0:
        lipid_preproc_defs['BILAYER_LIPIDHEAD_FC'] = stage_lipidhead_fc
    full_topology = parse_itp_file(dopc_param_file, preprocessor_defines=lipid_preproc_defs)
    
    # Try to find DOPC or similar molecule
    dopc_molecule = None
    for mol_name in full_topology['molecules'].keys():
        if 'DOPC' in mol_name.upper() or 'DOP' in mol_name.upper():
            dopc_molecule = mol_name
            break
    
    if dopc_molecule:
        dopc_topology = parse_itp_file(
            dopc_param_file, dopc_molecule, preprocessor_defines=lipid_preproc_defs
        )
    else:
        available_molecules = list(full_topology['molecules'].keys())
        raise ValueError(f"FATAL ERROR: DOPC molecule not found in '{dopc_param_file}'.\n"
                        f"  Available molecules: {available_molecules}\n"
                        f"  Please ensure DOPC is defined in the phospholipid parameter file.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    dopc_bead_types = [atom['type'] for atom in dopc_topology['atoms']]
    dopc_charges = [atom['charge'] for atom in dopc_topology['atoms']]
    # Create mapping from atom names to bead types and charges
    dopc_atom_to_type = {atom['atom']: atom['type'] for atom in dopc_topology['atoms']}
    dopc_atom_to_charge = {atom['atom']: atom['charge'] for atom in dopc_topology['atoms']}
    print(f"Read DOPC topology: {len(dopc_bead_types)} bead types from {dopc_param_file}")
    
    # Parse ion topologies from ITP file
    ion_param_file = pick_ff_file("dry_martini_v2.1_ions.itp", required=False)
    if ion_param_file:
        ion_topology = parse_itp_file(ion_param_file)
        # Extract NA and CL atoms specifically
        # In MARTINI ion ITP files, residue name is "ION" and atom name is "NA" or "CL"
        na_atoms = [atom for atom in ion_topology['atoms'] if atom['atom'].upper() == 'NA']
        cl_atoms = [atom for atom in ion_topology['atoms'] if atom['atom'].upper() == 'CL']
        
        # For standard ions, use the first occurrence (standard chloride is TQ5, not SQ5n which is for acetate)
        na_bead_types = [na_atoms[0]['type']] if na_atoms else []
        na_charges = [na_atoms[0]['charge']] if na_atoms else []
        cl_bead_types = [cl_atoms[0]['type']] if cl_atoms else []  # Use first CL (TQ5), not acetate CL (SQ5n)
        cl_charges = [cl_atoms[0]['charge']] if cl_atoms else []
        print(f"Ion topology loaded: NA={len(na_bead_types)} type(s), CL={len(cl_bead_types)} type(s)")
    else:
        na_bead_types, na_charges = [], []
        cl_bead_types, cl_charges = [], []
        print("Ion topology file not found in selected FF")
    
    # Parse water topology from ITP file
    water_param_file = pick_ff_file("dry_martini_v2.1_solvents.itp", required=False)
    if water_param_file:
        water_topology = parse_itp_file(water_param_file)
        water_atoms = [atom for atom in water_topology['atoms'] if atom['residue'].upper() == 'W']
        water_bead_types = [atom['type'] for atom in water_atoms]
        water_charges = [atom['charge'] for atom in water_atoms]
        print(f"Read water topology: {len(water_bead_types)} bead types from {water_param_file}")
    else:
        water_bead_types, water_charges = [], []
        print("Water topology file not found in selected FF")
    
    # Read bead masses from force field file
    mass_file = martini_param_file
    martini_masses = read_martini_masses(mass_file)
    print(f"Read {len(martini_masses)} atom type masses from force field file")
    
    # Read DOPC bonds and angles from parsed topology (for both lipid and mixed systems)
    dopc_bonds = [(bond['i'], bond['j']) for bond in dopc_topology['bonds']]
    dopc_bond_lengths = [bond['r0'] for bond in dopc_topology['bonds']]  # nm
    dopc_bond_force_constants = [bond['k'] for bond in dopc_topology['bonds']]  # kJ/mol/nm²
    
    dopc_angles = [(angle['i'], angle['j'], angle['k']) for angle in dopc_topology['angles']]
    dopc_angle_equil_deg = [angle['theta0'] for angle in dopc_topology['angles']]  # degrees
    dopc_angle_force_constants = [angle['force_k'] for angle in dopc_topology['angles']]  # kJ/mol/rad²
    dopc_position_restraints = dopc_topology.get('position_restraints', [])
    
    print(f"Read DOPC connectivity: {len(dopc_bonds)} bonds, {len(dopc_angles)} angles")
    
    # Validate that required topology data was found
    if not dopc_bonds:
        raise ValueError(f"FATAL ERROR: No DOPC bonds found in topology from '{dopc_param_file}'.\n"
                        f"  This indicates incomplete molecule definition.\n"
                        f"  Please ensure DOPC bonds are properly defined in the phospholipid parameter file.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    if not dopc_angles:
        raise ValueError(f"FATAL ERROR: No DOPC angles found in topology from '{dopc_param_file}'.\n"
                        f"  This indicates incomplete molecule definition.\n"
                        f"  Please ensure DOPC angles are properly defined in the phospholipid parameter file.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    # Unit conversions for mapping native MARTINI units into the active
    # simulation unit system. These must be provided explicitly rather than
    # baked into the generator.
    energy_conversion_raw = os.environ.get('UPSIDE_MARTINI_ENERGY_CONVERSION', '').strip()
    length_conversion_raw = os.environ.get('UPSIDE_MARTINI_LENGTH_CONVERSION', '').strip()
    if not energy_conversion_raw:
        raise ValueError("Missing required environment variable UPSIDE_MARTINI_ENERGY_CONVERSION")
    if not length_conversion_raw:
        raise ValueError("Missing required environment variable UPSIDE_MARTINI_LENGTH_CONVERSION")
    energy_conversion = float(energy_conversion_raw)
    length_conversion = float(length_conversion_raw)
    coulomb_constant_native = float(os.environ.get('UPSIDE_MARTINI_COULOMB_CONSTANT_NATIVE', str(138.935458 / 15.0)))
    if energy_conversion <= 0.0:
        raise ValueError("UPSIDE_MARTINI_ENERGY_CONVERSION must be positive")
    if length_conversion <= 0.0:
        raise ValueError("UPSIDE_MARTINI_LENGTH_CONVERSION must be positive")

    # Pressure conversion (for NPT simulations):
    # 1 atm = 101325 Pa = 101.325 kJ/m³
    # 1 atm = 1.01325e-28 kJ/Å³
    # 1 atm = (1.01325e-28 kJ/Å³) × (6.02214e23 mol⁻¹) / (2.914952774272 kJ/mol)
    # 1 atm = 0.000020933 E_up/Å³
    # 1 bar = 0.986923 atm = 0.000020659 E_up/Å³
    # For LAMMPS real units compatibility, use 1 bar = 0.000020659 E_up/Å³
    pressure_conversion_bar_to_eup = 0.000020659  # bar → E_up/Å³

    # Bonds: kJ/mol/nm² → E_up/Å²
    # = (kJ/mol → E_up) / (nm² → Å²)
    # = (1/energy_conversion) / (length_conversion²)
    bond_conversion = 1.0 / (energy_conversion * length_conversion ** 2)  # ≈ 0.003406

    # Angles: kJ/mol/deg² → E_up/deg²
    # = (kJ/mol → E_up) (degrees stay the same)
    angle_conversion = 1.0 / energy_conversion  # ≈ 0.343

    # Dihedrals: kJ/mol → E_up
    dihedral_conversion = 1.0 / energy_conversion  # ≈ 0.343

    print("\n=== MARTINI Unit Conversions ===")
    print(f"Bond lengths (nm -> Å): 0.40 nm -> 4.0 Å")
    print(f"Bond force constants (kJ/mol/nm² -> E_up/Å²): 7000.0 -> {7000.0 * bond_conversion:.3f}")
    print(f"Angle equilibrium (degrees): 108.0°")
    print(f"Angle force constants (kJ/mol/deg² -> E_up/deg²): 21.5 -> {21.5 * angle_conversion:.6f}")
    print(f"Dihedral force constants (kJ/mol -> E_up): 400.0 -> {400.0 * dihedral_conversion:.6f}")
    print(f"Pressure (bar -> E_up/Å³): 1 bar -> {pressure_conversion_bar_to_eup:.9f}")
    print(f"Energy conversion factor: {energy_conversion} (kJ/mol -> E_up)")
    print(f"Length conversion: 1 nm = {length_conversion} Å, so 1 nm² = {length_conversion**2} Å²")
    print(f"Bond conversion factor: {bond_conversion:.6f} (divide by energy_conv × length_conv²)")
    print(f"Angle conversion factor: {angle_conversion:.6f} (divide by energy_conv)")
    print(f"Dihedral conversion factor: {dihedral_conversion:.6f} (divide by energy_conv)")
    
    # Read PDB file
    if strict_from_martini_pdb or include_protein:
        input_pdb_file = runtime_input_pdb_path(workflow_dir, pdb_id)
        print(f"\nUsing MARTINI PDB as base structure: {input_pdb_file}")
    
    # Read PDB and populate arrays
    print(f"\n=== Reading PDB Structure ===")
    initial_positions = []
    atom_types = []
    charges = []
    residue_ids = []
    atom_names = []
    residue_names = []
    chain_ids = []
    seg_ids = []
    
    protein_bonds, protein_angles, protein_dihedrals, protein_constraints = [], [], [], []
    protein_exclusions = []
    protein_position_restraints = []

    protein_atoms_for_mapping, _ = parse_pdb(Path(input_pdb_file))
    protein_like_atoms = [
        atom for atom in protein_atoms_for_mapping
        if atom["resname"].upper() in {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
            'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX',
        }
    ]
    if protein_like_atoms:
        protein_backbone_atoms_for_mapping, _ = extract_protein_backbone_atoms_from_aa(
            protein_atoms_for_mapping
        )
        bb_type_by_residue = map_backbone_types_from_martinize_fallback(
            protein_backbone_atoms_for_mapping
        )
    else:
        protein_backbone_atoms_for_mapping = []
        bb_type_by_residue = {}
    
    
    protein_residue_names = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
        'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
        'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX'
    }
    with open(input_pdb_file, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
                
            # Parse PDB line
            atom_name_raw = line[12:16].strip()
            atom_name = atom_name_raw.upper()
            atom_names.append(atom_name)
            residue_id = int(line[22:26])
            residue_ids.append(residue_id)
            # Read 4-character residue field to support DOPC while preserving
            # standard 3-letter protein residues (e.g., " ASN" -> "ASN").
            residue_name = line[17:21].strip().upper()
            residue_names.append(residue_name)
            chain_id = line[21:22].strip()
            chain_ids.append(chain_id)
            seg_id = line[72:76].strip()
            seg_ids.append(seg_id)
            
            
            
            # Extract coordinates
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            initial_positions.append([x, y, z])
            
            is_protein = (residue_name in protein_residue_names)
            
            # Map to MARTINI type based on context
            if is_protein:
                if atom_name not in {"N", "CA", "C", "O", "BB"}:
                    raise ValueError(
                        f"FATAL ERROR: Protein atom '{atom_name}' in residue '{residue_name}' is not one of N/CA/C/O/BB.\n"
                        "Runtime AA protein representation must be backbone-only."
                    )
                icode = line[26:27].strip()
                res_key = (chain_id, residue_id, icode, residue_name)
                martini_type = bb_type_by_residue.get(res_key, "P5")
                charge = float(BB_TYPE_CHARGE.get(martini_type, 0.0))
            elif residue_name == 'DOPC' or residue_name == 'DOP':
                # For DOPC, use the topology from parameter file (for both lipid and mixed systems)
                if atom_name in dopc_atom_to_type:
                    martini_type = dopc_atom_to_type[atom_name]
                    charge = dopc_atom_to_charge[atom_name]
                else:
                    available_atom_names = sorted(dopc_atom_to_type.keys())
                    raise ValueError(f"FATAL ERROR: Unknown DOPC atom '{atom_name}' in residue '{residue_name}'.\n"
                                   f"  Available DOPC atom names: {available_atom_names}\n"
                                   f"  This indicates incomplete DOPC topology mapping.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
            elif residue_name == 'W':
                if not water_bead_types:
                    raise ValueError(f"FATAL ERROR: No water bead types found in topology.\n"
                                   f"  This indicates incomplete water parameter file.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                martini_type = water_bead_types[0]
                charge = water_charges[0]
            elif residue_name == 'NA':
                if not na_bead_types:
                    raise ValueError(f"FATAL ERROR: No sodium bead types found in topology.\n"
                                   f"  This indicates incomplete ion parameter file.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                martini_type = na_bead_types[0]
                charge = na_charges[0]
            elif residue_name == 'CL':
                if not cl_bead_types:
                    raise ValueError(f"FATAL ERROR: No chloride bead types found in topology.\n"
                                   f"  This indicates incomplete ion parameter file.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                martini_type = cl_bead_types[0]
                charge = cl_charges[0]
            else:
                raise ValueError(f"FATAL ERROR: Unknown residue type '{residue_name}' for atom '{atom_name}'.\n"
                               f"  Supported residue types: PROTEIN, DOPC, W, NA, CL\n"
                               f"  This indicates incomplete system definition.\n"
                               f"  Aborting to prevent incorrect simulation results.")
            
            atom_types.append(martini_type)
            charges.append(charge)
    
    # Convert to numpy arrays
    initial_positions = np.array(initial_positions, dtype=float)
    atom_types = np.array(atom_types)
    charges = np.array(charges, dtype=float)
    residue_ids = np.array(residue_ids, dtype=int)
    atom_names = np.array(atom_names)
    residue_names = np.array(residue_names)
    chain_ids = np.array(chain_ids)
    seg_ids = np.array(seg_ids)
    n_atoms = len(initial_positions)

    # --- CG lipid coarse-graining: collapse 14-bead DOPC → single 6D vector particles ---
    (initial_positions, atom_types, charges, residue_ids,
     atom_names, residue_names, chain_ids, seg_ids,
     lipid_directions, cg_lipid_indices, lipid_to_atom_map,
     ref_bead_positions, old_to_new_map) = build_cg_lipid_array(
        initial_positions, atom_types, charges, residue_ids,
        atom_names, residue_names, chain_ids, seg_ids)
    n_atoms = len(initial_positions)
    n_cg_lipids = len(cg_lipid_indices)

    cg_lipid_orientation_indices = np.zeros(0, dtype=np.int32)
    cg_lipid_orientation_length_ang = float(
        os.environ.get('UPSIDE_CG_LIPID_ORIENTATION_LENGTH', '5.0')
    )
    cg_lipid_orientation_bond_fc = float(
        os.environ.get('UPSIDE_CG_LIPID_ORIENTATION_BOND_FC', '20.0')
    )
    cg_lipid_orientation_mass_g = 0.0
    if n_cg_lipids > 0:
        bead_masses = np.array(
            [martini_masses.get(bt, 0.0) for bt in dopc_bead_types],
            dtype=np.float64,
        )
        inertia_values = []
        for rel_pos, direction in zip(ref_bead_positions, lipid_directions):
            direction = np.asarray(direction, dtype=np.float64)
            direction /= max(float(np.linalg.norm(direction)), 1e-12)
            rel_pos = np.asarray(rel_pos, dtype=np.float64)
            parallel = rel_pos @ direction
            r2_perp = np.sum(rel_pos * rel_pos, axis=1) - parallel * parallel
            inertia_values.append(float(np.sum(bead_masses * np.maximum(r2_perp, 0.0))))
        mean_inertia = float(np.mean(inertia_values)) if inertia_values else 0.0
        cg_lipid_orientation_mass_g = max(
            mean_inertia / max(cg_lipid_orientation_length_ang ** 2, 1e-12),
            12.0,
        )
        cg_lipid_orientation_indices = np.arange(
            n_atoms, n_atoms + n_cg_lipids, dtype=np.int32
        )
        orientation_positions = (
            initial_positions[cg_lipid_indices]
            + cg_lipid_orientation_length_ang * lipid_directions
        )
        initial_positions = np.vstack([initial_positions, orientation_positions])
        atom_types = np.concatenate([
            atom_types.astype('<U4'),
            np.array(['CGLD'] * n_cg_lipids, dtype='<U4'),
        ])
        charges = np.concatenate([charges, np.zeros(n_cg_lipids, dtype=charges.dtype)])
        residue_ids = np.concatenate([residue_ids, residue_ids[cg_lipid_indices]])
        atom_names = np.concatenate([
            atom_names.astype('<U4'),
            np.array(['CGLD'] * n_cg_lipids, dtype='<U4'),
        ])
        residue_names = np.concatenate([residue_names, residue_names[cg_lipid_indices]])
        chain_ids = np.concatenate([chain_ids, chain_ids[cg_lipid_indices]])
        seg_ids = np.concatenate([seg_ids, seg_ids[cg_lipid_indices]])
        n_atoms = len(initial_positions)
        print(
            f"Added {n_cg_lipids} hidden CGLD orientation sites "
            f"(length={cg_lipid_orientation_length_ang:.3f} Å, "
            f"mass={cg_lipid_orientation_mass_g:.3f} g/mol)"
        )

    # Read box dimensions from CRYST1 record
    print(f"Reading box dimensions from {input_pdb_file}...")
    with open(input_pdb_file, 'r') as f:
        for line in f:
            if line.startswith('CRYST1'):
                fields = line.split()
                if len(fields) >= 4:
                    x_len = float(fields[1])
                    y_len = float(fields[2])
                    z_len = float(fields[3])
                    print(f"Found CRYST1 record: X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f} Angstroms")
                    break
        else:
            raise ValueError(f"FATAL ERROR: No CRYST1 record found in PDB file '{input_pdb_file}'.\n"
                           f"  Box dimensions are required for proper simulation setup.\n"
                           f"  Please ensure the PDB file contains a CRYST1 record with box dimensions.\n"
                           f"  Aborting to prevent incorrect simulation results.")
    
    # Print system parameters
    print(f"Box dimensions: X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f} Angstroms")
    print(f"Box volume: {x_len * y_len * z_len:.1f} Å³")
    if n_cg_lipids > 0 and int(os.environ.get("UPSIDE_CG_LIPID_CONDITION_INITIAL", "1")):
        target_nn = float(os.environ.get("UPSIDE_CG_LIPID_MIN_LEAFLET_NN", "7.0"))
        max_step = float(os.environ.get("UPSIDE_CG_LIPID_CONDITION_MAX_STEP", "0.25"))
        n_iter = int(os.environ.get("UPSIDE_CG_LIPID_CONDITION_STEPS", "100"))
        cgl_pos = initial_positions[cg_lipid_indices]
        leaflet_split = float(np.median(cgl_pos[:, 2]))
        leaflet_masks = (cgl_pos[:, 2] <= leaflet_split, cgl_pos[:, 2] > leaflet_split)

        def _same_leaflet_nn_stats() -> tuple[float, float]:
            values = []
            for leaflet_mask in leaflet_masks:
                ids = np.where(leaflet_mask)[0]
                if ids.size < 2:
                    continue
                xy = initial_positions[cg_lipid_indices[ids], :2]
                for local_i in range(ids.size):
                    delta = xy - xy[local_i]
                    delta[:, 0] -= x_len * np.round(delta[:, 0] / x_len)
                    delta[:, 1] -= y_len * np.round(delta[:, 1] / y_len)
                    dist = np.sqrt(np.sum(delta * delta, axis=1))
                    dist[local_i] = np.inf
                    values.append(float(np.min(dist)))
            if not values:
                return float("nan"), float("nan")
            arr = np.asarray(values, dtype=np.float64)
            return float(np.min(arr)), float(np.percentile(arr, 5.0))

        start_min, start_p05 = _same_leaflet_nn_stats()
        for _ in range(max(n_iter, 0)):
            delta_xy = np.zeros((n_cg_lipids, 2), dtype=np.float64)
            for leaflet_mask in leaflet_masks:
                ids = np.where(leaflet_mask)[0]
                for local_a in range(ids.size):
                    ia = int(ids[local_a])
                    for local_b in range(local_a + 1, ids.size):
                        ib = int(ids[local_b])
                        dxy = initial_positions[cg_lipid_indices[ib], :2] - initial_positions[cg_lipid_indices[ia], :2]
                        dxy[0] -= x_len * np.round(dxy[0] / x_len)
                        dxy[1] -= y_len * np.round(dxy[1] / y_len)
                        dist = float(np.sqrt(np.dot(dxy, dxy)))
                        if dist <= 1e-8 or dist >= target_nn:
                            continue
                        push = 0.5 * (target_nn - dist) * dxy / dist
                        delta_xy[ia] -= push
                        delta_xy[ib] += push
            norms = np.sqrt(np.sum(delta_xy * delta_xy, axis=1))
            max_norm = float(np.max(norms)) if norms.size else 0.0
            if max_norm <= 1e-8:
                break
            if max_norm > max_step:
                delta_xy *= max_step / max_norm
            initial_positions[cg_lipid_indices, :2] += delta_xy
            if cg_lipid_orientation_indices.size == n_cg_lipids:
                initial_positions[cg_lipid_orientation_indices, :2] += delta_xy
            initial_positions[cg_lipid_indices, 0] %= x_len
            initial_positions[cg_lipid_indices, 1] %= y_len
            if cg_lipid_orientation_indices.size == n_cg_lipids:
                initial_positions[cg_lipid_orientation_indices, 0] %= x_len
                initial_positions[cg_lipid_orientation_indices, 1] %= y_len
        end_min, end_p05 = _same_leaflet_nn_stats()
        print(
            "Conditioned initial CGL same-leaflet XY spacing: "
            f"min/p05 {start_min:.3f}/{start_p05:.3f} -> {end_min:.3f}/{end_p05:.3f} Å"
        )
    print(f"Total atoms: {n_atoms}")
    
    # Group atoms into molecules by chain (for proteins) or residue ID (for other molecules)
    molecules = []
    current_mol_atoms = []
    current_mol_names = []
    current_mol_indices = []
    current_resid = None
    current_resname = None
    current_chain = None
    
    # Define protein residue names
    protein_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX',
    }
    
    def pick_chain_id(idx):
        # Prefer PDB chain ID to preserve distinct protein chains in multi-chain
        # systems where segid may be shared (e.g., "PROA" for all protein atoms).
        if chain_ids[idx]:
            return chain_ids[idx]
        return seg_ids[idx]

    for i, (resid, resname, atom_name) in enumerate(zip(residue_ids, residue_names, atom_names)):
        # Determine molecule type
        if resname in protein_residues:
            mol_type = 'PROTEIN'
        else:
            # Normalize DOP (resname in PDB) to DOPC for reporting and selection
            mol_type = 'DOPC' if resname == 'DOP' else resname
            
        # Start new molecule if protein chain changes or residue ID/name changes for non-proteins
        chain_id = pick_chain_id(i) if mol_type == 'PROTEIN' else ""
        start_new = False
        if mol_type == 'PROTEIN':
            if chain_id != current_chain:
                start_new = True
        else:
            if resid != current_resid or resname != current_resname:
                start_new = True

        if start_new:
            if current_mol_atoms:
                molecules.append((current_mol_type, current_mol_atoms, current_mol_indices))
            current_mol_atoms = [atom_name]
            current_mol_names = [atom_name]
            current_mol_indices = [i]
            current_resid = resid
            current_resname = resname
            current_mol_type = mol_type
            current_chain = chain_id
        else:
            current_mol_atoms.append(atom_name)
            current_mol_indices.append(i)
    
    if current_mol_atoms:
        molecules.append((current_mol_type, current_mol_atoms, current_mol_indices))

    # Build per-atom molecule indices
    molecule_ids = np.zeros(n_atoms, dtype=np.int32)
    for mol_idx, (_, _, atom_indices) in enumerate(molecules):
        for atom_idx in atom_indices:
            molecule_ids[atom_idx] = mol_idx
    if cg_lipid_orientation_indices.size == n_cg_lipids:
        molecule_ids[cg_lipid_orientation_indices] = molecule_ids[cg_lipid_indices]

    # Validate that each non-protein molecule contains only atoms of the same residue type
    for i, (mol_type, atoms, indices) in enumerate(molecules):
        residue_names_in_mol = [residue_names[idx] for idx in indices]
        unique_resnames = set(residue_names_in_mol)
        if mol_type != 'PROTEIN':
            if len(unique_resnames) > 1:
                raise ValueError(f"FATAL ERROR: Molecule {i} contains mixed residue types: {unique_resnames}\n"
                               f"  This indicates incorrect molecule grouping.\n"
                               f"  Molecule atoms: {atoms}\n"
                               f"  Residue names: {residue_names_in_mol}\n"
                               f"  Aborting to prevent incorrect simulation results.")
        else:
            pass
    
    # Count molecules by type, but group all protein residues together
    mol_counts = Counter()
    protein_residue_count = 0
    protein_sequence = []
    protein_sequence_seen = set()
    for idx, (resid, resname) in enumerate(zip(residue_ids, residue_names)):
        if resname not in protein_residues:
            continue
        key = (pick_chain_id(idx), int(resid))
        if key in protein_sequence_seen:
            continue
        protein_sequence_seen.add(key)
        protein_sequence.append(normalize_resname(str(resname)))

    protein_residue_keys = set()
    for resid, resname, idx in zip(residue_ids, residue_names, range(n_atoms)):
        if resname in protein_residues:
            protein_residue_keys.add((pick_chain_id(idx), resid))

    for mol_type, _, _ in molecules:
        if mol_type == 'PROTEIN':
            mol_counts[mol_type] += 1
        else:
            mol_counts[mol_type] += 1
    
    protein_residue_count = len(protein_residue_keys)
    if protein_residue_count > 0:
        protein_chains = mol_counts.get('PROTEIN', 0)
        mol_counts['PROTEIN'] = f"{protein_chains} chain(s) ({protein_residue_count} residues)"
    
    dopc_count = mol_counts.get('DOPC', 0)
    water_count = mol_counts.get('W', 0)
    if n_cg_lipids > 0:
        mol_counts['DOPC'] = n_cg_lipids
        dopc_count = n_cg_lipids
    
    print(f"\n=== Molecule Summary ===")
    for moltype, count in mol_counts.items():
        print(f"{moltype}: {count} molecules")
    
    # Create bonds and angles
    print(f"\n=== Creating Connectivity ===")
    
    # Initialize lists for bonds and angles
    bonds_list = []
    bond_lengths_list = []
    bond_force_constants_list = []
    angles_list = []
    angle_equil_deg_list = []
    angle_force_constants_list = []
    dihedrals_list = []
    dihedral_equil_deg_list = []
    dihedral_force_constants_list = []
    dihedral_type_list = []
    lipid_restraint_indices = []
    lipid_restraint_ref_pos = []
    lipid_restraint_spring_xyz = []
    
    # CG lipids have no internal bonds/angles — skip per-bead DOPC topology.
    if n_cg_lipids > 0:
        print(f"CG lipid mode: {n_cg_lipids} DOPC lipids coarse-grained to single 6D vector particles "
              f"with dynamic orientation sites")
        for cgl_i, orient_i in zip(cg_lipid_indices, cg_lipid_orientation_indices):
            bonds_list.append([int(cgl_i), int(orient_i)])
            bond_lengths_list.append(cg_lipid_orientation_length_ang)
            bond_force_constants_list.append(cg_lipid_orientation_bond_fc)
        print(
            f"Added {n_cg_lipids} CGL-CGLD orientation bonds "
            f"(r0={cg_lipid_orientation_length_ang:.3f} Å, "
            f"k={cg_lipid_orientation_bond_fc:.3f} E_up/Å²)"
        )
    
    # Create protein connectivity if available
    protein_bond_count = 0
    protein_angle_count = 0
    protein_dihedral_count = 0
    protein_constraint_count = 0
    
    if protein_bonds or protein_constraints:
        print("\n=== Protein Connectivity ===")
        print(f"Found {len(protein_bonds)} bonds, {len(protein_angles)} angles, {len(protein_dihedrals)} dihedrals")
        print(f"Found {len(protein_constraints)} constraints, {len(protein_position_restraints)} position restraints")
        
        # Add protein bonds to the bond list
        for i, j, r0_nm, k_kj in protein_bonds:
            # Convert MARTINI units to UPSIDE units
            r0_angstrom = r0_nm * 10.0  # nm to Å
            k_upside = k_kj * bond_conversion  # kJ/mol/nm² to E_up/Å²
            
            # Add to bond list (assuming protein atoms come first in the system)
            bonds_list.append([i, j])
            bond_lengths_list.append(r0_angstrom)
            bond_force_constants_list.append(k_upside)
            protein_bond_count += 1
        
        # Add protein constraints as bonds with large spring constants
        for i, j, r0_nm, k_kj in protein_constraints:
            # Convert MARTINI units to UPSIDE units
            r0_angstrom = r0_nm * 10.0  # nm to Å
            k_upside = k_kj * bond_conversion  # kJ/mol/nm² to E_up/Å²
            
            # Add to bond list
            bonds_list.append([i, j])
            bond_lengths_list.append(r0_angstrom)
            bond_force_constants_list.append(k_upside)
            protein_constraint_count += 1
        
        # Add protein angles to the angle list
        for i, j, k, theta0_deg, k_kj in protein_angles:
            # Convert MARTINI units to UPSIDE units
            theta0_upside = theta0_deg  # degrees (same unit)
            k_upside = k_kj * angle_conversion  # kJ/mol/deg² to E_up/deg²
            
            # Add to angle list
            angles_list.append([i, j, k])
            angle_equil_deg_list.append(theta0_upside)
            angle_force_constants_list.append(k_upside)
            protein_angle_count += 1
        
        # Add protein dihedrals to the dihedral list
        for i, j, k, l, phi0_deg, k_kj, func_type in protein_dihedrals:
            # Convert MARTINI units to UPSIDE units
            phi0_upside = phi0_deg  # degrees (same unit)
            k_upside = k_kj * dihedral_conversion  # kJ/mol to E_up
            
            # Add to dihedral list
            dihedrals_list.append([i, j, k, l])
            dihedral_equil_deg_list.append(phi0_upside)
            dihedral_force_constants_list.append(k_upside)
            dihedral_type_list.append(func_type)
            protein_dihedral_count += 1
        
        print(f"Added {protein_bond_count} protein bonds")
        print(f"Added {protein_constraint_count} protein constraints (as bonds with large spring constants)")
        print(f"Added {protein_angle_count} protein angles")
        print(f"Added {protein_dihedral_count} protein dihedrals")
    else:
        print("\nNo explicit protein MARTINI connectivity is used for AA backbone runtime mode")
    
    print(f"Total system bonds: {len(bonds_list)}")
    print(f"Total system angles: {len(angles_list)}")
    print(f"Total system dihedrals: {len(dihedrals_list)}")
    
    # Center and wrap positions
    print(f"\n=== Preparing Final Structure ===")
    center_mask = atom_types != "CGLD"
    center = np.mean(initial_positions[center_mask], axis=0)
    centered_positions = initial_positions - center
    half_box = np.array([x_len/2, y_len/2, z_len/2])
    centered_positions = (centered_positions + half_box) % (2*half_box) - half_box
    final_positions = centered_positions

    particle_class = np.empty(n_atoms, dtype='S10')
    for i, resname in enumerate(residue_names):
        if resname in protein_residues:
            particle_class[i] = b"PROTEIN"
        elif resname == 'DOP' or resname == 'DOPC':
            particle_class[i] = b"LIPID"
        elif resname == 'W':
            particle_class[i] = b"WATER"
        elif resname in ('NA', 'CL'):
            particle_class[i] = b"ION"
        else:
            particle_class[i] = b"OTHER"

    if int(os.environ.get("UPSIDE_INITIAL_DEBUG_ONLY", "0")):
        debug_prefix = os.environ.get("UPSIDE_INITIAL_DEBUG_PREFIX", f"{pdb_id}.{stage}")
        write_initial_structure_debug_outputs(
            debug_dir=Path(run_dir) / "debug",
            prefix=debug_prefix,
            pos=final_positions,
            atom_types=atom_types,
            atom_names=atom_names,
            particle_classes=_decode_h5_text_array(particle_class),
            residue_ids=residue_ids,
            molecule_ids=molecule_ids,
            charges=charges,
            bond_ids=np.array(bonds_list, dtype=np.int32),
            bond_r0=np.array(bond_lengths_list, dtype=np.float64),
            bond_k=np.array(bond_force_constants_list, dtype=np.float64),
            box_lengths=np.array([x_len, y_len, z_len], dtype=np.float64),
        )
        print("Initial debug-only mode complete; skipping UPSIDE HDF5 pair-table generation.")
        return
    
    # Create UPSIDE input file
    print(f"\n=== Creating UPSIDE Input File ===")
    input_file = f"{run_dir}/test.input.up"
    
    with tb.open_file(input_file, 'w') as t:
        # Create input group (required by UPSIDE)
        input_grp = t.create_group(t.root, 'input')
        
        # Create position array with correct format
        pos = np.zeros((n_atoms, 3, 1), dtype='f4')
        pos[:,:,0] = final_positions
        pos_array = t.create_array(input_grp, 'pos', obj=pos)
        pos_array._v_attrs.arguments = np.array([b'pos'])
        pos_array._v_attrs.shape = pos.shape
        pos_array._v_attrs.n_atoms = n_atoms
        pos_array._v_attrs.n_frames = 1
        pos_array._v_attrs.dim = 3
        pos_array._v_attrs.initialized = True
        
        # Create velocity array (required by UPSIDE)
        velocity = np.zeros((n_atoms, 3), dtype='f4')
        vel_array = t.create_array(input_grp, 'vel', obj=velocity)
        vel_array._v_attrs.arguments = np.array([b'vel'])
        vel_array._v_attrs.shape = velocity.shape
        vel_array._v_attrs.n_atoms = n_atoms
        vel_array._v_attrs.dim = 3
        vel_array._v_attrs.initialized = True
        
        # Create momentum array (required by UPSIDE)
        momentum = np.zeros((n_atoms, 3, 1), dtype='f4')
        mom_array = t.create_array(input_grp, 'mom', obj=momentum)
        mom_array._v_attrs.arguments = np.array([b'mom'])
        mom_array._v_attrs.shape = momentum.shape
        mom_array._v_attrs.n_atoms = n_atoms
        mom_array._v_attrs.dim = 3
        mom_array._v_attrs.initialized = True
        
        # Create mass array (required by UPSIDE)
        # Register CG lipid mass (sum of 14 DOPC bead masses) if not already present
        if n_cg_lipids > 0 and "CGL" not in martini_masses:
            cgl_mass = sum(martini_masses.get(bt, 0.0) for bt in dopc_bead_types)
            martini_masses["CGL"] = cgl_mass
            print(f"CG lipid mass: {cgl_mass:.1f} g/mol (sum of {len(dopc_bead_types)} DOPC beads)")
        if n_cg_lipids > 0:
            martini_masses["CGLD"] = cg_lipid_orientation_mass_g
            print(f"CG lipid orientation-site mass: {cg_lipid_orientation_mass_g:.3f} g/mol")

        # Extend martini_table with effective CGL LJ parameters for isotropic
        # CG lipid ↔ dryMARTINI particle interactions (ions, water, BB, env).
        if n_cg_lipids > 0 and ("CGL", "CGL") not in martini_table:
            cgl_sigma_cap_nm = _cgl_effective_lj_sigma_cap_nm()
            martini_h5_path = Path(run_dir) / "martini.h5"
            if martini_h5_path.exists():
                # Read effective LJ from pre-built martini.h5 to ensure
                # consistency with the combined energy grids already stored there.
                with h5py.File(martini_h5_path, "r") as _mh5:
                    eff_grp = _mh5["cg_lipid_table/effective_lj"]
                    target_types = [t.decode("ascii") if isinstance(t, bytes) else str(t)
                                  for t in eff_grp["target_types"][:]]
                    sigmas = eff_grp["sigma_nm"][:]
                    epsilons = eff_grp["epsilon_kj_mol"][:]
                for tgt_type, sigma_nm, epsilon_kj in zip(target_types, sigmas, epsilons):
                    capped_sigma = min(float(sigma_nm), cgl_sigma_cap_nm)
                    martini_table[("CGL", tgt_type)] = (capped_sigma, float(epsilon_kj))
                    martini_table[(tgt_type, "CGL")] = (capped_sigma, float(epsilon_kj))
            else:
                from martini_build_tables import _compute_cgl_effective_lj_params
                ref_nm = ref_bead_positions[0].astype(np.float64) * 0.1  # first lipid, Å → nm
                ff_pair_params = {}
                for (t1, t2), (sigma, eps) in martini_table.items():
                    ff_pair_params[(t1, t2)] = {"sigma_nm": sigma, "epsilon_kj_mol": eps}
                effective_lj = _compute_cgl_effective_lj_params(
                    ref_bead_positions_nm=ref_nm,
                    bead_types=list(dopc_bead_types),
                    pair_params=ff_pair_params,
                )
                for tgt_type, eff in effective_lj.items():
                    s = min(float(eff["sigma_nm"]), cgl_sigma_cap_nm)
                    e = eff["epsilon_kj_mol"]
                    martini_table[("CGL", tgt_type)] = (s, e)
                    martini_table[(tgt_type, "CGL")] = (s, e)
            # CGL↔CGL is represented by the full directional cg_lipid_pair node.
            # CGL↔SC is represented by the full directional cg_lipid_sc node.
            # Zero both in the particles table — the quadsplines are the sole
            # interaction; an isotropic MartiniPotential on top would double-count.
            n_sc_zeroed = 0
            if martini_h5_path.exists():
                with h5py.File(martini_h5_path, "r") as _mh5:
                    cglt = _mh5.get("cg_lipid_table")
                    if cglt is not None:
                        sc_grp = cglt.get("cg_lipid_sc")
                        if sc_grp is not None and "sc_bead_types" in sc_grp:
                            for bt_bytes in sc_grp["sc_bead_types"][:]:
                                bt = bt_bytes.decode("ascii") if isinstance(bt_bytes, bytes) else str(bt_bytes)
                                if ("CGL", bt) in martini_table:
                                    martini_table[("CGL", bt)] = (0.0, 0.0)
                                    martini_table[(bt, "CGL")] = (0.0, 0.0)
                                    n_sc_zeroed += 1
            martini_table[("CGL", "CGL")] = (0.0, 0.0)
            n_types = len([k for k in martini_table if k[0] == "CGL"])
            print(f"Extended martini_table with CGL entries for {n_types} target types "
                  f"(CGL↔CGL and CGL↔SC({n_sc_zeroed} types) supplied by quadsplines)")
        if n_cg_lipids > 0:
            all_martini_types = set()
            for t1, t2 in martini_table.keys():
                all_martini_types.add(t1)
                all_martini_types.add(t2)
            all_martini_types.update(str(x) for x in np.unique(atom_types))
            all_martini_types.add("CGLD")
            for tgt_type in all_martini_types:
                martini_table[("CGLD", tgt_type)] = (0.0, 0.0)
                martini_table[(tgt_type, "CGLD")] = (0.0, 0.0)

        mass = np.zeros(n_atoms, dtype='f4')
        for i, atom_type in enumerate(atom_types):
            # Get mass from force field file, raise error if not found
            if atom_type not in martini_masses:
                raise ValueError(f"FATAL ERROR: Mass not found for atom type '{atom_type}' (atom index {i}).\n"
                                f"  Available atom types with masses: {sorted(martini_masses.keys())}\n"
                                f"  This indicates incomplete force field parameters.\n"
                                f"  Aborting to prevent incorrect simulation results.")
            # Divide by 12.0 for reduced mass units (1 unit = 12 g/mol)
            mass[i] = martini_masses[atom_type] / 12.0
        
        mass_array = t.create_array(input_grp, 'mass', obj=mass)
        mass_array._v_attrs.arguments = np.array([b'mass'])
        mass_array._v_attrs.shape = mass.shape
        mass_array._v_attrs.n_atoms = n_atoms
        mass_array._v_attrs.initialized = True
        
        print("Hybrid stage files use AA backbone carriers (N/CA/C/O) for protein runtime representation")
        
        # Create stage-specific parameters group (always create this)
        stage_grp = t.create_group(input_grp, 'stage_parameters')
        stage_grp._v_attrs.enable = 1
        stage_grp._v_attrs.current_stage = b'minimization'
        
        # Store minimization stage bond parameters (large spring constants)
        min_bonds_grp = t.create_group(stage_grp, 'minimization_bonds')
        min_bond_fc = np.array([1000000.0] * len(protein_bonds), dtype='f4')
        t.create_array(min_bonds_grp, 'force_constants', obj=min_bond_fc)
        
        # Store production stage bond parameters (regular spring constants)
        prod_bonds_grp = t.create_group(stage_grp, 'production_bonds')
        prod_bond_fc = np.array([bond[3] for bond in protein_bonds], dtype='f4')  # Regular k values
        t.create_array(prod_bonds_grp, 'force_constants', obj=prod_bond_fc)
        
        # Store minimization stage angle parameters (NORMANG angles)
        min_angles_grp = t.create_group(stage_grp, 'minimization_angles')
        min_angle_fc = np.array([angle[4] for angle in protein_angles], dtype='f4')
        t.create_array(min_angles_grp, 'force_constants', obj=min_angle_fc)
        
        # Store production stage angle parameters (regular angles)
        prod_angles_grp = t.create_group(stage_grp, 'production_angles')
        # For production, we'd use different angle parameters (not NORMANG)
        # For now, use the same as minimization but this could be different
        prod_angle_fc = np.array([angle[4] for angle in protein_angles], dtype='f4')
        t.create_array(prod_angles_grp, 'force_constants', obj=prod_angle_fc)
        
        print(f"Stage-specific parameters: minimization bonds={len(min_bond_fc)}, production bonds={len(prod_bond_fc)}")
        print(f"Stage-specific parameters: minimization angles={len(min_angle_fc)}, production angles={len(prod_angle_fc)}")
        
        # ===================== NPT BAROSTAT CONFIGURATION =====================
        # Create barostat configuration group for NPT simulations
        # Settings are read from environment variables (set by run_sim_bilayer.sh)
        barostat_enable = int(os.environ.get('UPSIDE_NPT_ENABLE', '0'))
        if barostat_enable:
            print(f"\n=== Creating NPT Barostat Configuration ===")
            barostat_grp = t.create_group(input_grp, 'barostat')
            barostat_grp._v_attrs.enable = barostat_enable
            # Default: 1 bar = 0.000020659 E_up/Angstrom^3 (from 1 atm = 0.000020933215)
            barostat_grp._v_attrs.target_p_xy = float(os.environ.get('UPSIDE_NPT_TARGET_PXY', '0.000020659'))
            barostat_grp._v_attrs.target_p_z = float(os.environ.get('UPSIDE_NPT_TARGET_PZ', '0.000020659'))
            barostat_grp._v_attrs.tau_p = float(os.environ.get('UPSIDE_NPT_TAU', '1.0'))
            # Legacy isotropic compressibility (kept for compatibility)
            legacy_compressibility = float(os.environ.get('UPSIDE_NPT_COMPRESSIBILITY', '14.521180763676'))
            barostat_grp._v_attrs.compressibility = legacy_compressibility
            # Axis-specific compressibility for semi-isotropic membrane coupling
            barostat_grp._v_attrs.compressibility_xy = float(
                os.environ.get('UPSIDE_NPT_COMPRESSIBILITY_XY', str(legacy_compressibility))
            )
            barostat_grp._v_attrs.compressibility_z = float(
                os.environ.get('UPSIDE_NPT_COMPRESSIBILITY_Z', str(legacy_compressibility))
            )
            barostat_grp._v_attrs.interval = int(os.environ.get('UPSIDE_NPT_INTERVAL', '10'))
            barostat_grp._v_attrs.semi_isotropic = int(os.environ.get('UPSIDE_NPT_SEMI', '1'))
            barostat_grp._v_attrs.debug = int(os.environ.get('UPSIDE_NPT_DEBUG', '1'))
            # Barostat type: 0 = Berendsen (default), 1 = Parrinello-Rahman
            barostat_grp._v_attrs.type = params['barostat_type']
            print(f"  Enabled: {barostat_enable}")
            print(f"  Type: {'Parrinello-Rahman' if barostat_grp._v_attrs.type == 1 else 'Berendsen'}")
            print(f"  Target Pxy: {barostat_grp._v_attrs.target_p_xy} E_up/Angstrom^3 (~1 bar)")
            print(f"  Target Pz: {barostat_grp._v_attrs.target_p_z} E_up/Angstrom^3 (~1 bar)")
            print(f"  Compressibility XY: {barostat_grp._v_attrs.compressibility_xy} Angstrom^3/E_up")
            print(f"  Compressibility Z: {barostat_grp._v_attrs.compressibility_z} Angstrom^3/E_up")
            print(f"  Tau_p: {barostat_grp._v_attrs.tau_p}")
            print(f"  Interval: {barostat_grp._v_attrs.interval} steps")
        else:
            print(f"\n=== NPT Barostat Disabled (NVT mode) ===")
        
        # Create type array
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
        particle_class_array._v_attrs.description = b"Per-atom class: PROTEIN, LIPID, WATER, ION, OTHER"
        
        # Create charges array
        charge_array = t.create_array(input_grp, 'charges', obj=charges)
        charge_array._v_attrs.arguments = np.array([b'charges'])
        charge_array._v_attrs.shape = charges.shape
        charge_array._v_attrs.n_atoms = n_atoms
        charge_array._v_attrs.initialized = True
        
        # Create residue IDs array
        residue_array = t.create_array(input_grp, 'residue_ids', obj=residue_ids)
        residue_array._v_attrs.arguments = np.array([b'residue_ids'])
        residue_array._v_attrs.shape = residue_ids.shape
        residue_array._v_attrs.n_atoms = n_atoms
        residue_array._v_attrs.initialized = True
        residue_array._v_attrs.description = b"PDB residue indices per atom"

        # Create molecule IDs array (per-atom molecule index)
        molecule_array = t.create_array(input_grp, 'molecule_ids', obj=molecule_ids)
        molecule_array._v_attrs.arguments = np.array([b'molecule_ids'])
        molecule_array._v_attrs.shape = molecule_ids.shape
        molecule_array._v_attrs.n_atoms = n_atoms
        molecule_array._v_attrs.initialized = True
        molecule_array._v_attrs.description = b"Per-atom molecule index (0-based, contiguous)"

        # Store atom names for MARTINI-specific selections (e.g., BB backbone)
        atom_name_array = t.create_array(input_grp, 'atom_names', obj=atom_names.astype('S4'))
        atom_name_array._v_attrs.arguments = np.array([b'atom_names'])
        atom_name_array._v_attrs.shape = atom_names.shape
        atom_name_array._v_attrs.n_atoms = n_atoms
        atom_name_array._v_attrs.initialized = True

        # Store atom role names from PDB (BB, SC1, SC2, W, etc.)
        atom_roles = atom_names.astype('S4')
        atom_role_array = t.create_array(input_grp, 'atom_roles', obj=atom_roles)
        atom_role_array._v_attrs.arguments = np.array([b'atom_roles'])
        atom_role_array._v_attrs.shape = atom_roles.shape
        atom_role_array._v_attrs.n_atoms = n_atoms
        atom_role_array._v_attrs.initialized = True
        atom_role_array._v_attrs.description = b"PDB atom role names (BB, SC1, SC2, W, etc.)"

        if protein_sequence:
            t.create_array(input_grp, 'sequence', obj=np.asarray([x.encode('ascii') for x in protein_sequence], dtype='S3'))
        
        # Create potential group (required by UPSIDE)
        potential_grp = t.create_group(input_grp, 'potential')
        
        # Create MARTINI potential with proper parameters
        martini_potential = t.create_group(potential_grp, 'martini_potential')
        martini_potential._v_attrs.arguments = np.array([b'pos'])
        martini_potential._v_attrs.potential_type = b'lj_coulomb'
        # Note: epsilon and sigma attributes are required by C++ interface but not used in computation
        # They will be calculated from the coefficients array after it's created
        martini_potential._v_attrs.lj_cutoff = 12.0
        martini_potential._v_attrs.coul_cutoff = 12.0
        martini_potential._v_attrs.energy_conversion_kj_per_eup = energy_conversion
        martini_potential._v_attrs.length_conversion_angstrom_per_nm = length_conversion
        martini_potential._v_attrs.coulomb_constant_native_kj_mol_nm_e2 = coulomb_constant_native
        martini_potential._v_attrs.n_types = 1
        martini_potential._v_attrs.n_params = 4
        martini_potential._v_attrs.cutoff = 12.0
        martini_potential._v_attrs.cache_buffer = 1.0
        martini_potential._v_attrs.initialized = True
        force_cap = float(os.environ.get('UPSIDE_FORCE_CAP', '0'))
        martini_potential._v_attrs.force_cap = force_cap
        
        martini_potential._v_attrs.x_len = x_len
        martini_potential._v_attrs.y_len = y_len
        martini_potential._v_attrs.z_len = z_len
        # Optional softening controls via environment variables
        # UPSIDE_SOFTEN_COULOMB: 1 to enable Slater softening for Coulomb
        # UPSIDE_SLATER_ALPHA: float value for Slater alpha (1/Angstrom)
        # UPSIDE_SOFTEN_LJ: 1 to enable soft-core LJ
        # UPSIDE_LJ_ALPHA: float value for LJ softening alpha (dimensionless)
        soften_coul = int(os.environ.get('UPSIDE_SOFTEN_COULOMB', '0'))
        slater_alpha = float(os.environ.get('UPSIDE_SLATER_ALPHA', '1.0'))
        soften_lj = int(os.environ.get('UPSIDE_SOFTEN_LJ', '0'))
        lj_alpha = float(os.environ.get('UPSIDE_LJ_ALPHA', '1.0'))
        
        # Set stage-specific softening parameters
        martini_potential._v_attrs.coulomb_soften = params['coulomb_soften']
        if params['coulomb_soften']:
            martini_potential._v_attrs.slater_alpha = params['slater_alpha']
        martini_potential._v_attrs.lj_soften = params['lj_soften']
        if params['lj_soften']:
            martini_potential._v_attrs.lj_soften_alpha = params['lj_alpha']
        
        # Periodic boundary potential removed - using NVT ensemble without boundaries


        # Create atom indices and charges arrays for the potential
        t.create_array(martini_potential, 'atom_indices', obj=np.arange(n_atoms))
        t.create_array(martini_potential, 'charges', obj=charges)
        
        # Create pairs and optimized coefficient indices for non-bonded interactions.
        # The full coefficient table is O(N^2) and can OOM for large membrane boxes;
        # store unique coefficient rows plus one int index per pair instead.
        # Create sets for exclusions (MARTINI uses nrexcl=1, so only 1-2 exclusions)
        bonded_pairs_12 = set()  # Directly bonded (1-2) - full exclusion
        additional_exclusions = set()  # Additional exclusions from ITP file
        
        # Add 1-2 exclusions from bond list
        for bond in bonds_list:
            sorted_bond = (min(bond[0], bond[1]), max(bond[0], bond[1]))
            bonded_pairs_12.add(sorted_bond)
        
        # Add additional exclusions from protein ITP file
        if protein_exclusions:
            for exclusion in protein_exclusions:
                sorted_exclusion = (min(exclusion[0], exclusion[1]), max(exclusion[0], exclusion[1]))
                additional_exclusions.add(sorted_exclusion)
        
        # Generate all unique pairs (i < j) with proper exclusions (nrexcl=1)
        excluded_12_count = 0
        excluded_additional_count = 0
        total_candidate_pairs = n_atoms * (n_atoms - 1) // 2
        total_pairs_written = 0
        unique_coeffs = []
        coeff_to_index = {}
        coeff_counts = []
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
                raise ValueError(f"FATAL ERROR: Missing interaction parameters for bead type pair ({type1}, {type2})\n"
                               f"  Atom index: {atom_i} ({type1})\n"
                               f"  This indicates incomplete MARTINI force field parameters.\n"
                               f"  Available bead types in parameter table: {available_types}\n"
                               f"  Aborting to prevent incorrect simulation results.")

            epsilon = np.float32(epsilon_kj / energy_conversion)
            sigma = np.float32(sigma_nm * length_conversion)
            q1 = np.float32(q1_raw)
            q2 = np.float32(q2_raw)
            coeff_key = (float(epsilon), float(sigma), float(q1), float(q2))
            coeff_index = coeff_to_index.get(coeff_key)
            if coeff_index is None:
                coeff_index = len(unique_coeffs)
                coeff_to_index[coeff_key] = coeff_index
                unique_coeffs.append(coeff_key)
                coeff_counts.append(0)
            return coeff_index

        def weighted_median(values_and_counts, total_count):
            midpoint = total_count // 2
            running = 0
            for value, count in sorted(values_and_counts):
                running += count
                if running > midpoint:
                    return value
            return values_and_counts[-1][0]

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
        additional_by_i = defaultdict(set)
        for i, j in additional_exclusions:
            additional_by_i[i].add(j)

        for i in range(n_atoms):
            if i + 1 >= n_atoms:
                continue
            js = np.arange(i + 1, n_atoms, dtype=np.int32)
            bonded_js = bonded_by_i.get(i)
            additional_js = additional_by_i.get(i)
            if bonded_js:
                bonded_mask = np.isin(js, np.fromiter(bonded_js, dtype=np.int32), assume_unique=True)
                excluded_12_count += int(np.count_nonzero(bonded_mask))
                js = js[~bonded_mask]
            if additional_js and js.size:
                additional_candidates = np.fromiter(
                    (j for j in additional_js if not bonded_js or j not in bonded_js),
                    dtype=np.int32,
                )
                if additional_candidates.size:
                    additional_mask = np.isin(js, additional_candidates, assume_unique=True)
                    excluded_additional_count += int(np.count_nonzero(additional_mask))
                    js = js[~additional_mask]
            if js.size == 0:
                continue

            pairs_chunk = np.empty((js.size, 2), dtype=np.int32)
            pairs_chunk[:, 0] = i
            pairs_chunk[:, 1] = js

            coeff_indices = np.empty(js.size, dtype=np.int64)
            class_i = int(atom_class_ids[i])
            js_class_ids = atom_class_ids[js]
            for class_j in np.unique(js_class_ids):
                mask = js_class_ids == class_j
                coeff_index = coefficient_index_for_classes(class_i, int(class_j), i)
                count = int(np.count_nonzero(mask))
                coeff_indices[mask] = coeff_index
                coeff_counts[coeff_index] += count

            pairs_array.append(pairs_chunk)
            coeff_index_array.append(coeff_indices)
            total_pairs_written += int(js.size)
        
        print(f"Excluded {excluded_12_count} 1-2 bonded pairs from non-bonded interactions (nrexcl=1)")
        print(f"Excluded {excluded_additional_count} additional pairs from ITP exclusions")
        print(f"Wrote {total_pairs_written} non-bonded pairs with {len(unique_coeffs)} unique coefficient rows")
        martini_potential._v_attrs.optimized_format = 1
        t.create_array(martini_potential, 'coefficients', obj=np.array(unique_coeffs, dtype='f4'))
        
        # Calculate representative epsilon and sigma from the coefficients array
        # These are required by the C++ interface but not used in computation
        if total_pairs_written > 0:
            epsilon_counts = [(coeff[0], coeff_counts[idx]) for idx, coeff in enumerate(unique_coeffs)]
            sigma_counts = [(coeff[1], coeff_counts[idx]) for idx, coeff in enumerate(unique_coeffs)]
            median_epsilon = weighted_median(epsilon_counts, total_pairs_written)
            median_sigma = weighted_median(sigma_counts, total_pairs_written)
            
            martini_potential._v_attrs.epsilon = median_epsilon
            martini_potential._v_attrs.sigma = median_sigma
        else:
            raise ValueError("FATAL ERROR: No interaction coefficients found.\n"
                           f"  This indicates no non-bonded interactions were defined.\n"
                           f"  Aborting to prevent incorrect simulation results.")

        # Add bonded potentials mirroring original run_martini.py
        # Bonds: dist_spring
        if bonds_list:
            bond_group = t.create_group(potential_grp, 'dist_spring')
            bond_group._v_attrs.arguments = np.array([b'pos'])
            bond_group._v_attrs.initialized = True

            bond_group._v_attrs.x_len = x_len
            bond_group._v_attrs.y_len = y_len
            bond_group._v_attrs.z_len = z_len
            t.create_array(bond_group, 'id', obj=np.array(bonds_list, dtype=int))
            t.create_array(bond_group, 'equil_dist', obj=np.array(bond_lengths_list, dtype='f4'))
            t.create_array(bond_group, 'spring_const', obj=np.array(bond_force_constants_list, dtype='f4'))
            # Compatibility dataset expected by some builds
            t.create_array(bond_group, 'bonded_atoms', obj=np.ones(len(bonds_list), dtype='i4'))

        # Angles: angle_spring
        if angles_list:
            angle_group = t.create_group(potential_grp, 'angle_spring')
            angle_group._v_attrs.arguments = np.array([b'pos'])
            angle_group._v_attrs.initialized = True

            angle_group._v_attrs.x_len = x_len
            angle_group._v_attrs.y_len = y_len
            angle_group._v_attrs.z_len = z_len
            t.create_array(angle_group, 'id', obj=np.array(angles_list, dtype=int))
            t.create_array(angle_group, 'equil_angle_deg', obj=np.array(angle_equil_deg_list, dtype='f4'))
            t.create_array(angle_group, 'spring_const', obj=np.array(angle_force_constants_list, dtype='f4'))

        # Dihedrals: dihedral_spring (if present)
        if dihedrals_list:
            dihedral_group = t.create_group(potential_grp, 'dihedral_spring')
            dihedral_group._v_attrs.arguments = np.array([b'pos'])
            dihedral_group._v_attrs.initialized = True

            dihedral_group._v_attrs.x_len = x_len
            dihedral_group._v_attrs.y_len = y_len
            dihedral_group._v_attrs.z_len = z_len
            t.create_array(dihedral_group, 'id', obj=np.array(dihedrals_list, dtype=int))
            # Some builds of UPSIDE expect 'equil_dist' for this potential; write both names for compatibility
            eq_deg = np.array(dihedral_equil_deg_list, dtype='f4')
            t.create_array(dihedral_group, 'equil_angle_deg', obj=eq_deg)
            t.create_array(dihedral_group, 'equil_dist', obj=eq_deg)
            t.create_array(dihedral_group, 'spring_const', obj=np.array(dihedral_force_constants_list, dtype='f4'))
            # Store dihedral type information (1=periodic, 2=harmonic)
            t.create_array(dihedral_group, 'dihedral_type', obj=np.array(dihedral_type_list, dtype=int))

        # Position restraints (dry MARTINI lipid-head ramp before production)
        if lipid_restraint_indices:
            restraint_group = t.create_group(potential_grp, 'restraint_position')
            restraint_group._v_attrs.arguments = np.array([b'pos'])
            restraint_group._v_attrs.initialized = True
            t.create_array(
                restraint_group, 'restraint_indices', obj=np.array(lipid_restraint_indices, dtype='i4')
            )
            t.create_array(
                restraint_group, 'ref_pos', obj=np.array(lipid_restraint_ref_pos, dtype='f4')
            )
            spring_xyz = np.array(lipid_restraint_spring_xyz, dtype='f4')
            t.create_array(restraint_group, 'spring_const_xyz', obj=spring_xyz)
            # Backward-compatible scalar spring constant for older readers.
            t.create_array(restraint_group, 'spring_const', obj=np.max(spring_xyz, axis=1).astype('f4'))

        # --- CG Lipid 6D vector particle nodes ---
        if n_cg_lipids > 0:
            # compose_vector6d: combines pos (3D) with static direction (3D) → 6D
            compose_grp = t.create_group(potential_grp, 'compose_vector6d')
            compose_grp._v_attrs.arguments = np.array([b'pos'])
            compose_grp._v_attrs.initialized = True
            compose_grp._v_attrs.x_len = x_len
            compose_grp._v_attrs.y_len = y_len
            compose_grp._v_attrs.z_len = z_len
            compose_grp._v_attrs.orientation_length_ang = cg_lipid_orientation_length_ang
            t.create_array(compose_grp, 'elem_index',
                           obj=cg_lipid_indices.astype('i4'))
            t.create_array(compose_grp, 'direction',
                           obj=lipid_directions.astype('f4'))
            t.create_array(compose_grp, 'orientation_index',
                           obj=cg_lipid_orientation_indices.astype('i4'))
            print(f"Injected compose_vector6d node: {n_cg_lipids} CG lipid 6D particles")

        # Store old→new atom index remapping for hybrid mapping injection
        if n_cg_lipids > 0:
            remap_grp = t.create_group(input_grp, 'hybrid_remap')
            t.create_array(remap_grp, 'old_to_new', obj=old_to_new_map.astype('i4'))
            remap_grp._v_attrs.n_old = int(len(old_to_new_map))
            remap_grp._v_attrs.n_new = int(n_atoms)

    print(f"Created UPSIDE input file: {input_file}")
    print(f"Preparation complete!")
    
    # Save preparation summary
    summary_file = f"{run_dir}/preparation_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=== MARTINI 3.0 Preparation Summary ===\n")
        f.write(f"PDB ID: {pdb_id}\n")
        f.write(f"Total atoms: {n_atoms}\n")
        f.write(f"Box dimensions: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å\n")
        f.write(f"Total bonds: {len(bonds_list)}\n")
        f.write(f"Total angles: {len(angles_list)}\n")
        f.write(f"Total dihedrals: {len(dihedrals_list)}\n")
        f.write(f"Protein bonds: {protein_bond_count}\n")
        f.write(f"Protein constraints: {protein_constraint_count}\n")
        f.write(f"Protein angles: {protein_angle_count}\n")
        f.write(f"Protein dihedrals: {protein_dihedral_count}\n")
        f.write(f"DOPC lipids: {dopc_count}\n")
        f.write(f"Water molecules: {water_count}\n")
        f.write(f"1-2 exclusions (nrexcl=1): {excluded_12_count}\n")
        f.write(f"Additional exclusions: {excluded_additional_count}\n")

    write_stage_debug_outputs(Path(input_file), debug_dir=Path(run_dir) / "debug")
    
    print(f"Preparation summary saved to: {summary_file}")



def require_group(root, path):
    if path not in root:
        raise ValueError(f"Missing group: {path}")
    return root[path]


def require_dataset(group, name):
    if name not in group:
        raise ValueError(f"Missing dataset: {group.name}/{name}")
    return group[name]


def decode_attr_string(value, default=""):
    if value is None:
        return default
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)


def validate_hybrid_mapping(mapping_h5: Path, n_atom: int | None = None):
    path = Path(mapping_h5).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(path)

    with h5py.File(path, "r") as h5:
        inp = require_group(h5, "/input")
        ctrl = require_group(inp, "hybrid_control")
        bb = require_group(inp, "hybrid_bb_map")
        env = require_group(inp, "hybrid_env_topology")

        for attr in [
            "enable",
            "activation_stage",
            "preprod_protein_mode",
            "exclude_intra_protein_martini",
            "schema_version",
        ]:
            if attr not in ctrl.attrs:
                raise ValueError(f"Missing hybrid_control attr: {attr}")

        if "protein_env_interface_scale" in ctrl.attrs:
            interface_scale = float(ctrl.attrs["protein_env_interface_scale"])
            if not np.isfinite(interface_scale) or interface_scale <= 0.0:
                raise ValueError(
                    "hybrid_control/protein_env_interface_scale must be finite and > 0"
                )

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
        bb_index_space = decode_attr_string(bb.attrs.get("atom_index_space", "runtime_n_atom"))
        bb_reference_space = bb_index_space == "protein_aa_pdb_0based"

        if "reference_atom_names" in bb and bb["reference_atom_names"][:].shape != (4,):
            raise ValueError("hybrid_bb_map/reference_atom_names must have shape (4,)")
        if "reference_atom_indices" in bb and bb["reference_atom_indices"][:].shape != (n_bb, 4):
            raise ValueError("hybrid_bb_map/reference_atom_indices must have shape (n_bb,4)")
        if "reference_atom_coords" in bb and bb["reference_atom_coords"][:].shape != (n_bb, 4, 3):
            raise ValueError("hybrid_bb_map/reference_atom_coords must have shape (n_bb,4,3)")
        if "bb_comment" in bb and bb["bb_comment"][:].shape != (n_bb,):
            raise ValueError("hybrid_bb_map/bb_comment must have shape (n_bb,)")

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
                if bb_reference_space:
                    if ai < 0:
                        raise ValueError(
                            f"BB reference target index invalid at row {i}, col {j}: {ai}"
                        )
                else:
                    if ai < 0 or ai >= n_atom_runtime:
                        raise ValueError(f"BB target index out of bounds at row {i}, col {j}: {ai}")
                    if membership[ai] < 0:
                        raise ValueError(f"BB target index is not protein atom at row {i}, col {j}: {ai}")

        n_protein = int(np.sum(membership >= 0))
        n_env = int(np.sum(membership < 0))
        print(f"OK: {path}")
        print(f"  n_atom={membership.shape[0]} n_protein={n_protein} n_env={n_env}")
        print(f"  n_bb={n_bb} bb_index_space={bb_index_space}")


def decode_stage_label(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore").strip().lower()
    return str(value).strip().lower()


def env_bool(name, default=False):
    raw = os.getenv(name)
    if raw is None:
        return default
    value = str(raw).strip().lower()
    return value not in ("0", "false", "no", "off", "")


def recenter_protein_for_production(h5f, pos, box_lengths):
    if box_lengths is None:
        return pos
    if "/input/stage_parameters" not in h5f:
        return pos
    stage = decode_stage_label(h5f["/input/stage_parameters"].attrs.get("current_stage", b""))
    if stage != "production":
        return pos
    if "/input/hybrid_bb_map/bb_atom_index" not in h5f:
        return pos

    box = np.asarray(box_lengths, dtype=np.float64).reshape(3)
    if not np.all(np.isfinite(box)) or np.any(box <= 0.0):
        return pos

    bb_idx = h5f["/input/hybrid_bb_map/bb_atom_index"][:].astype(np.int32)
    bb_idx = bb_idx[(bb_idx >= 0) & (bb_idx < pos.shape[0])]
    if bb_idx.size == 0:
        return pos

    bb_xyz = pos[bb_idx, :, 0].astype(np.float64)
    bb_center = np.zeros(3, dtype=np.float64)
    for ax in range(3):
        coords_ax = np.mod(bb_xyz[:, ax], box[ax])
        angles = TWOPI * coords_ax / box[ax]
        s = np.mean(np.sin(angles))
        c = np.mean(np.cos(angles))
        if abs(s) < 1e-12 and abs(c) < 1e-12:
            bb_center[ax] = float(np.mean(coords_ax))
        else:
            ang = np.arctan2(s, c)
            if ang < 0.0:
                ang += TWOPI
            bb_center[ax] = box[ax] * ang / TWOPI

    target_center = 0.5 * box
    delta = target_center - bb_center
    pos[:, :, 0] = pos[:, :, 0] + delta.astype(pos.dtype, copy=False)[None, :]
    print(
        "Recentered production protein BB to box center: "
        f"before=({bb_center[0]:.3f}, {bb_center[1]:.3f}, {bb_center[2]:.3f}) "
        f"target=({target_center[0]:.3f}, {target_center[1]:.3f}, {target_center[2]:.3f})"
    )
    return pos


def set_initial_position(input_file, output_file):
    strict_copy = env_bool("UPSIDE_SET_INITIAL_STRICT_COPY", False)
    apply_recenter_production = env_bool(
        "UPSIDE_SET_INITIAL_RECENTER_PRODUCTION", default=(not strict_copy)
    )
    preserve_hybrid_transition = env_bool("UPSIDE_SET_INITIAL_PRESERVE_HYBRID_TRANSITION", False)
    time_step = float(os.environ.get("UPSIDE_SET_INITIAL_TIME_STEP", "0") or "0")

    with h5py.File(input_file, "r") as f:
        if "/output/pos" in f and f["/output/pos"].shape[0] > 0:
            last_pos = f["/output/pos"][-1, 0, :, :]
            last_pos = last_pos[:, :, np.newaxis]
        else:
            last_pos = f["/input/pos"][:, :, 0]
            last_pos = last_pos[:, :, np.newaxis]

        last_box = None
        if "/output/box" in f:
            box_data = f["/output/box"][:]
            if box_data.size > 0:
                last_box = box_data[-1]
                if len(last_box.shape) == 2 and last_box.shape[1] == 3:
                    last_box = last_box[0]
        if last_box is None and "/input/potential/martini_potential" in f:
            pot_grp = f["/input/potential/martini_potential"]
            if all(k in pot_grp.attrs for k in ("x_len", "y_len", "z_len")):
                last_box = np.array(
                    [pot_grp.attrs["x_len"], pot_grp.attrs["y_len"], pot_grp.attrs["z_len"]]
                )
        last_mom = None
        if "/output/mom" in f and f["/output/mom"].shape[0] > 0:
            last_mom = f["/output/mom"][-1, 0, :, :]
            last_mom = last_mom[:, :, np.newaxis]

        source_stage = ""
        if "/input/stage_parameters" in f:
            source_stage = f["/input/stage_parameters"].attrs.get("current_stage", b"")
            if isinstance(source_stage, (bytes, np.bytes_)):
                source_stage = source_stage.decode("utf-8", errors="ignore")
            else:
                source_stage = str(source_stage)
            source_stage = source_stage.strip()
        source_transition_start = 0
        if "/input/hybrid_control" in f:
            source_transition_start = int(
                f["/input/hybrid_control"].attrs.get("sc_env_transition_step_start", 0)
            )
        source_elapsed_steps = 0
        if preserve_hybrid_transition and source_stage == "production" and time_step > 0.0:
            if "/output/time" in f and f["/output/time"].shape[0] > 0:
                last_time = float(f["/output/time"][-1])
                source_elapsed_steps = max(0, int(round(last_time / time_step)))
        exact_production_restart = preserve_hybrid_transition and source_stage == "production" and last_mom is not None

    with h5py.File(output_file, "r+") as f:
        target_stage = ""
        if "/input/stage_parameters" in f:
            target_stage = f["/input/stage_parameters"].attrs.get("current_stage", b"")
            if isinstance(target_stage, (bytes, np.bytes_)):
                target_stage = target_stage.decode("utf-8", errors="ignore")
            else:
                target_stage = str(target_stage)
            target_stage = target_stage.strip()

        target_pos = f["/input/pos"][:]
        target_n = target_pos.shape[0]
        source_n = last_pos.shape[0]
        if source_n != target_n:
            merged = target_pos.copy()
            merged[: min(source_n, target_n), :, :] = last_pos[: min(source_n, target_n), :, :]
            last_pos = merged

        if strict_copy and not apply_recenter_production:
            print("Strict handoff mode: preserving exact coordinates from previous stage output.")
        if exact_production_restart and apply_recenter_production:
            raise ValueError(
                "Production restart must preserve saved coordinates exactly when saved momentum is reused."
            )

        release_from_preprod = (source_stage != "production" and target_stage == "production")
        if apply_recenter_production:
            last_pos = recenter_protein_for_production(f, last_pos, last_box)

        if "/input/pos" in f:
            del f["/input/pos"]
        f.create_dataset("/input/pos", data=last_pos)

        if last_mom is not None:
            if last_mom.shape[0] != target_n:
                merged_mom = f["/input/mom"][:] if "/input/mom" in f else np.zeros((target_n, 3, 1), dtype=last_mom.dtype)
                merged_mom[: min(last_mom.shape[0], target_n), :, :] = last_mom[: min(last_mom.shape[0], target_n), :, :]
                last_mom = merged_mom
            if release_from_preprod and "/input/hybrid_env_topology/protein_membership" in f:
                protein_membership = np.asarray(
                    f["/input/hybrid_env_topology/protein_membership"][:], dtype=np.int32
                )
                if protein_membership.shape[0] == last_mom.shape[0]:
                    protein_mask = protein_membership >= 0
                    if np.any(protein_mask):
                        last_mom[protein_mask, :, :] = 0.0
            if "/input/mom" in f:
                del f["/input/mom"]
            mom_ds = f.create_dataset("/input/mom", data=last_mom)
            mom_ds.attrs["restart_valid"] = np.int8(1)
        elif "/input/mom" in f:
            f["/input/mom"].attrs["restart_valid"] = np.int8(0)

        if preserve_hybrid_transition and source_stage == "production" and "/input/hybrid_control" in f:
            transition_step = source_transition_start + source_elapsed_steps
            f["/input/hybrid_control"].attrs["sc_env_transition_step_start"] = np.int32(transition_step)

        if last_box is not None and "/input/potential/martini_potential" in f:
            pot_grp = f["/input/potential/martini_potential"]
            pot_grp.attrs["x_len"] = float(last_box[0])
            pot_grp.attrs["y_len"] = float(last_box[1])
            pot_grp.attrs["z_len"] = float(last_box[2])
            print(
                f"Updated box dimensions: x={last_box[0]:.3f}, y={last_box[1]:.3f}, z={last_box[2]:.3f}"
            )


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


def parse_itp_residue_names(path: Path):
    resnames = []
    seen = set()
    in_atoms = False
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            low = line.lower()
            if low in {"[ atoms ]", "[atoms]"}:
                in_atoms = True
                continue
            if line.startswith("[") and line.endswith("]"):
                in_atoms = False
                continue
            if not in_atoms:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                resnr = int(parts[2])
            except ValueError:
                continue
            if parts[4].strip().upper() != "BB" or resnr in seen:
                continue
            seen.add(resnr)
            resnames.append(parts[3].strip().upper())
    return resnames


def resolve_sequence(inp, residue_count: int, protein_itp: Path | None):
    mismatch_notes = []
    if "sequence" in inp:
        sequence = [normalize_resname(x) for x in decode_string_array(inp["sequence"])]
        if len(sequence) == residue_count:
            return sequence
        mismatch_notes.append(f"/input/sequence has {len(sequence)} residues")

    if protein_itp is not None:
        if not protein_itp.exists():
            raise ValueError(f"Protein ITP not found for sequence fallback: {protein_itp}")
        sequence = [normalize_resname(x) for x in parse_itp_residue_names(protein_itp)]
        if len(sequence) == residue_count:
            return sequence
        mismatch_notes.append(f"{protein_itp} has {len(sequence)} residues")

    detail = ", ".join(mismatch_notes) if mismatch_notes else "no sequence source was available"
    raise ValueError(
        f"Could not resolve a residue sequence matching hybrid_bb_map residue count {residue_count}: {detail}"
    )


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
        if "rotamer_prob_fixed" in sclib:
            rot_energy_fixed = sclib["rotamer_prob_fixed"][:].astype(np.float32).reshape(-1)
        elif "rotamer_prob" in sclib:
            rot_prob = sclib["rotamer_prob"][:].astype(np.float64)
            if rot_prob.ndim != 3:
                raise ValueError(f"Unsupported rotamer_prob shape in {sidechain_lib}: {rot_prob.shape}")
            rot_prob_mean = np.clip(rot_prob.mean(axis=(0, 1)), 1.0e-12, None)
            rot_energy_fixed = (-np.log(rot_prob_mean)).astype(np.float32)
        else:
            raise ValueError(
                f"Missing rotamer probability tables in {sidechain_lib}: need rotamer_prob_fixed or rotamer_prob"
            )
        bead_order = decode_string_array(sclib["bead_order"])
        bead_num = {name: i for i, name in enumerate(bead_order)}
        pair_interaction = sclib["pair_interaction"][:].astype(np.float32)

    count_by_n_rot = {}
    affine_residue = []
    layer_index = []
    beadtype_seq = []
    bead_type_index = []
    id_seq = []
    row_rotamer_index = []
    row_residue_table_index = []
    skipped = []

    for seq_idx, raw_resname in enumerate(sequence):
        resname = normalize_resname(raw_resname)
        restype_idx = sc_restype_num.get(resname)
        residue_table_idx = restype_to_index.get(resname)
        if restype_idx is None or residue_table_idx is None:
            skipped.append((seq_idx, resname))
            continue

        start, stop, n_bead = [int(x) for x in start_stop_bead[restype_idx]]
        if n_bead <= 0 or stop <= start:
            skipped.append((seq_idx, resname))
            continue
        n_rot = (stop - start) // n_bead
        if n_rot <= 0:
            skipped.append((seq_idx, resname))
            continue

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
            beadtype_seq.append(bead_name)
            bead_type_index.append(bead_num[bead_name])
            id_seq.append(rot_idx + (base_id << N_BIT_ROTAMER))
            row_rotamer_index.append(rot_idx)
            row_residue_table_index.append(residue_table_idx)

    if not layer_index:
        raise ValueError("No production rotamer rows could be generated from the sidechain library")

    layer_index_arr = np.asarray(layer_index, dtype=np.int32)
    return {
        "pair_interaction": pair_interaction,
        "affine_residue": np.asarray(affine_residue, dtype=np.int32),
        "layer_index": layer_index_arr,
        "beadtype_seq": np.asarray([np.bytes_(x) for x in beadtype_seq], dtype="S16"),
        "bead_type_index": np.asarray(bead_type_index, dtype=np.int32),
        "id_seq": np.asarray(id_seq, dtype=np.int32),
        "placement_data": rot_center_fixed[layer_index_arr].astype(np.float32),
        "placement_scalar_data": rot_energy_fixed[layer_index_arr][:, None].astype(np.float32),
        "row_rotamer_index": np.asarray(row_rotamer_index, dtype=np.int32),
        "row_residue_table_index": np.asarray(row_residue_table_index, dtype=np.int32),
        "skipped": skipped,
    }


def build_env_rows(inp, target_to_index):
    if "hybrid_env_topology" not in inp:
        raise ValueError("Missing /input/hybrid_env_topology in stage file")
    env_grp = inp["hybrid_env_topology"]
    if "protein_membership" not in env_grp:
        raise ValueError("Missing hybrid_env_topology/protein_membership in stage file")

    protein_membership = env_grp["protein_membership"][:].astype(np.int32)
    atom_types = decode_string_array(inp["type"])

    env_atom_index = []
    env_target_index = []
    for atom_idx, (membership, atom_type) in enumerate(zip(protein_membership.tolist(), atom_types)):
        if membership >= 0:
            continue
        target_idx = target_to_index.get(atom_type)
        if target_idx is None:
            continue
        env_atom_index.append(atom_idx)
        env_target_index.append(target_idx)

    if not env_atom_index:
        raise ValueError("No non-protein dry-MARTINI particles matched martini.h5 sc_table target_order")

    return np.asarray(env_atom_index, dtype=np.int32), np.asarray(env_target_index, dtype=np.int32)


def recreate_group(parent, name):
    if name in parent:
        del parent[name]
    return parent.create_group(name)


def load_upside_config(upside_home: Path):
    upside_config_py = upside_home / "py" / "upside_config.py"
    if not upside_config_py.exists():
        raise ValueError(f"upside_config.py not found: {upside_config_py}")
    spec = importlib.util.spec_from_file_location("upside_config_runtime", str(upside_config_py))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def remap_atom_index_array(arr, atom_map):
    out = arr.copy()
    if out.size == 0:
        return out
    finite_mask = np.isfinite(out)
    if not np.any(finite_mask):
        return out
    nonneg_mask = finite_mask & (out >= 0)
    if not np.any(nonneg_mask):
        return out
    idx = out[nonneg_mask].astype(np.int64)
    if np.any(idx >= atom_map.shape[0]):
        raise ValueError(
            f"Backbone atom remap index exceeds reference map size: max={int(idx.max())}, map_size={atom_map.shape[0]}"
        )
    mapped = atom_map[idx]
    if np.any(mapped < 0):
        bad = np.where(mapped < 0)[0][0]
        raise ValueError(
            f"Backbone atom remap produced negative target for reference index {int(idx[bad])}"
        )
    out[nonneg_mask] = mapped.astype(out.dtype, copy=False)
    return out


def inject_backbone_nodes(
    up_file: Path,
    sequence,
    affine_atoms,
    rama_library: Path,
    rama_sheet_mixing: Path,
    hbond_energy: Path,
    reference_state_rama: Path,
    upside_home: Path,
):
    uc = load_upside_config(upside_home)
    fasta_seq = np.asarray(sequence)
    ref_n_atom = 3 * len(sequence)
    spring_args = types.SimpleNamespace(bond_stiffness=48.0, angle_stiffness=175.0, omega_stiffness=30.0)
    chain_first_residue = np.array([], dtype=np.int32)

    with h5py.File(up_file, "r") as up:
        inp = up["input"]
        if "chain_break" in inp and "chain_first_residue" in inp["chain_break"]:
            chain_first_residue = inp["chain_break"]["chain_first_residue"][:].astype(np.int32)
        elif "hybrid_bb_map" in inp and "bb_chain_id" in inp["hybrid_bb_map"]:
            bb_chain_ids = decode_string_array(inp["hybrid_bb_map"]["bb_chain_id"])
            chain_first_residue, _ = derive_chain_break_metadata(
                [{"bb_chain": chain_id} for chain_id in bb_chain_ids]
            )
    n_chains = int(chain_first_residue.size) + 1
    chain_starts = np.append(np.array([0], dtype=np.int32), chain_first_residue) * 3
    # Mirror upside_config chain-break handling: exclude residues adjacent to
    # each chain boundary from inferred H-bond donor/acceptor construction.
    # For each chain first residue i (except chain 0), exclude residues i-1 and i.
    hbond_excluded_residues = np.array(
        sorted({int(i + j) for i in chain_first_residue.tolist() for j in (-1, 0)}),
        dtype=np.int32,
    )

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
                    category=np.exceptions.VisibleDeprecationWarning,
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

        uc.write_infer_H_O(fasta_seq, hbond_excluded_residues)
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


def inject_particles_table(up_file: Path, martini_h5: Path):
    up_file = Path(up_file).expanduser().resolve()
    martini_h5 = Path(martini_h5).expanduser().resolve()

    if not martini_h5.exists():
        raise SystemExit(f"ERROR: martini.h5 not found: {martini_h5}")

    with h5py.File(martini_h5, "r") as src:
        pg = src["particles"]
        eps_arr = pg["unique_eps_eup"][:].astype(np.float64)
        sig_arr = pg["unique_sig_ang"][:].astype(np.float64)
        qq_arr = pg["unique_charge_product"][:].astype(np.float64)
        combined_grids = pg["combined_energy_grids"][:].astype(np.float64)
        n_triples = len(eps_arr)

    has_zero_triple = np.any(
        (np.abs(eps_arr) < 1e-12)
        & (np.abs(sig_arr) < 1e-12)
        & (np.abs(qq_arr) < 1e-12)
    )
    if not has_zero_triple:
        eps_arr = np.concatenate([eps_arr, np.array([0.0], dtype=np.float64)])
        sig_arr = np.concatenate([sig_arr, np.array([0.0], dtype=np.float64)])
        qq_arr = np.concatenate([qq_arr, np.array([0.0], dtype=np.float64)])
        combined_grids = np.vstack([
            combined_grids,
            np.zeros((1, combined_grids.shape[1]), dtype=np.float64),
        ])
        n_triples = len(eps_arr)

    GRID_N = 1000
    R_MIN = 0.0
    R_MAX = 12.0

    with h5py.File(up_file, "r+") as up:
        g = up["/input/potential/martini_potential"]

        lj_soften = bool(int(g.attrs.get("lj_soften", 0)))
        lj_soften_alpha = float(g.attrs.get("lj_soften_alpha", 0.0))
        coulomb_soften = bool(int(g.attrs.get("coulomb_soften", 0)))
        slater_alpha = float(g.attrs.get("slater_alpha", 0.0))
        has_softening = (lj_soften and lj_soften_alpha > 0.0) or coulomb_soften

        if has_softening:
            coulomb_k_native = float(g.attrs["coulomb_constant_native_kj_mol_nm_e2"])
            energy_conv = float(g.attrs["energy_conversion_kj_per_eup"])
            length_conv = float(g.attrs["length_conversion_angstrom_per_nm"])
            coulomb_k = coulomb_k_native * (length_conv / energy_conv)
            lj_alpha = float(lj_soften_alpha) if lj_soften else 0.0
            slater = float(slater_alpha) if coulomb_soften else 0.0

            for idx in range(n_triples):
                eps = float(eps_arr[idx])
                sig = float(sig_arr[idx])
                qq = float(qq_arr[idx])
                for i in range(GRID_N):
                    r = R_MIN + i * (R_MAX - R_MIN) / (GRID_N - 1)
                    if r == 0.0:
                        r = 1.0e-6

                    # LJ component with optional soft-core
                    if eps != 0.0 and sig != 0.0:
                        if lj_soften and lj_alpha > 0.0:
                            x = r / sig
                            t = (x * x) ** 3 + lj_alpha
                            inv_t = 1.0 / t
                            lj = 4.0 * eps * (inv_t * inv_t - inv_t)
                        else:
                            x = sig / r
                            x6 = (x * x) ** 3
                            lj = 4.0 * eps * (x6 * x6 - x6)
                    else:
                        lj = 0.0

                    # Coulomb component with optional softening
                    if abs(qq) > 1e-10:
                        coul = coulomb_k * qq / r
                        if coulomb_soften:
                            alpha_r = slater * r
                            coul *= 1.0 - (1.0 + alpha_r * 0.5) * math.exp(-alpha_r)
                    else:
                        coul = 0.0

                    combined_grids[idx, i] = lj + coul

        g.create_dataset("unique_eps_eup", data=eps_arr, dtype=np.float64)
        g.create_dataset("unique_sig_ang", data=sig_arr, dtype=np.float64)
        g.create_dataset("unique_charge_product", data=qq_arr, dtype=np.float64)
        g.create_dataset("combined_energy_grids", data=combined_grids, dtype=np.float64)
        g.attrs["n_combined_triples"] = np.int32(n_triples)
        g.attrs["r_min_ang"] = np.float32(R_MIN)
        g.attrs["r_max_ang"] = np.float32(R_MAX)
        g.attrs["n_points"] = np.int32(GRID_N)


def inject_cg_lipid_nodes(
    up_file: Path,
    martini_h5: Path,
):
    """Inject CG lipid quadspline interaction nodes into the .up file.

    Reads interaction_param from martini.h5 (cg_lipid_table group) and
    writes cg_lipid_pair and cg_lipid_sc potential nodes.
    """
    up_file = Path(up_file).expanduser().resolve()
    martini_h5 = Path(martini_h5).expanduser().resolve()

    if not up_file.exists():
        raise SystemExit(f"ERROR: stage file not found: {up_file}")
    if not martini_h5.exists():
        raise SystemExit(f"ERROR: martini.h5 not found: {martini_h5}")

    with h5py.File(martini_h5, "r") as mh5:
        if "cg_lipid_table" not in mh5:
            print("  CG lipid tables not found in martini.h5, skipping CG lipid node injection")
            return
        cg_tab = mh5["cg_lipid_table"]

        has_pair = "cg_lipid_pair" in cg_tab
        has_sc = "cg_lipid_sc" in cg_tab

    with h5py.File(up_file, "r+") as up:
        inp = up["input"]
        pot = inp["potential"]
        box_attrs = {}
        if "martini_potential" in pot:
            martini_pot = pot["martini_potential"]
            for attr_name in ("x_len", "y_len", "z_len"):
                if attr_name in martini_pot.attrs:
                    box_attrs[attr_name] = np.float32(martini_pot.attrs[attr_name])

        if "compose_vector6d" not in pot:
            raise SystemExit(
                "ERROR: compose_vector6d node not found in .up file. "
                "Run convert_stage() with CG lipids first."
            )

        with h5py.File(martini_h5, "r") as mh5:
            if has_pair:
                pair_grp = mh5["cg_lipid_table/cg_lipid_pair"]

                # cg_lipid_pair: symmetric directional residual over CG lipid vectors.
                cg_pair = pot.create_group("cg_lipid_pair")
                cg_pair.attrs["initialized"] = True
                cg_pair.attrs["arguments"] = np.array([b"compose_vector6d"])
                for attr_name, attr_value in box_attrs.items():
                    cg_pair.attrs[attr_name] = attr_value
                cg_pair.attrs["angle_convention"] = "ang1=-n1_dot_n12;ang2=n2_dot_n12"
                cg_pair.attrs["radial_mode"] = "residual_zero_isotropic"
                for attr_name in ("knot_spacing_ang", "cutoff_ang", "taper_width_ang"):
                    if attr_name in pair_grp.attrs:
                        cg_pair.attrs[attr_name] = np.float32(pair_grp.attrs[attr_name])
                for attr_name in ("n_modes", "n_radial", "n_angular"):
                    if attr_name in pair_grp.attrs:
                        cg_pair.attrs[attr_name] = np.int32(pair_grp.attrs[attr_name])

                pi = cg_pair.create_group("pair_interaction")
                pi.create_dataset(
                    "interaction_param",
                    data=pair_grp["interaction_param"][:].astype(np.float32),
                )

                # Element mapping: CG lipids use identity mapping (1 CG type)
                # Shift ids by 4 bits so all shifted ids are unique —
                # PosQuadSplineInteraction::acceptable_id_pair rejects pairs
                # where (id1>>4) == (id2>>4).
                with h5py.File(up_file, "r") as up_r:
                    n_cg = up_r["input/potential/compose_vector6d/elem_index"].shape[0]
                pi.create_dataset("index", data=np.arange(n_cg, dtype=np.int32))
                pi.create_dataset("type", data=np.zeros(n_cg, dtype=np.int32))
                pi.create_dataset("id", data=(np.arange(n_cg, dtype=np.int32) << 4))

                n_param = int(pair_grp["interaction_param"].shape[-1])
                print(f"  Injected cg_lipid_pair: {n_cg} CG lipids, 1×1×{n_param} params")

            # Orientation spring: penalize angular deviation of CGL direction
            # from its initial reference. The single harmonic bond between CGL
            # and CGLD provides only quartic angular stiffness (∝ dθ⁴), allowing
            # free wobbling at small angles. This spring adds genuine quadratic
            # stiffness: E = k * (1 − n·n_ref) ≈ ½k·θ² for small θ.
            orient_k = float(os.environ.get("UPSIDE_CG_LIPID_ORIENT_SPRING_K", "50.0"))
            if orient_k > 0.0:
                with h5py.File(up_file, "r") as up_r:
                    compose_grp = up_r["input/potential/compose_vector6d"]
                    ref_dir = compose_grp["direction"][:].astype(np.float32)
                orient_grp = pot.create_group("cg_lipid_orientation_spring")
                orient_grp.attrs["initialized"] = True
                orient_grp.attrs["arguments"] = np.array([b"compose_vector6d"])
                orient_grp.attrs["k_orient"] = np.float32(orient_k)
                orient_grp.create_dataset("ref_dir", data=ref_dir)
                print(f"  Injected cg_lipid_orientation_spring: k={orient_k:.1f} E_up, "
                      f"{ref_dir.shape[0]} directions")

            if has_sc:
                # Check for SC placement node before creating cg_lipid_sc node
                with h5py.File(up_file, "r") as up_r:
                    has_sc_node = "placement_fixed_point_vector_only_CB" in up_r["input/potential"]
                    if has_sc_node:
                        sc_place = up_r["input/potential/placement_fixed_point_vector_only_CB"]
                        n_sc = int(sc_place["affine_residue"].shape[0])
                        cb_residue_idx = sc_place["affine_residue"][:].astype(np.int32)
                        sequence = [
                            s.decode("ascii") if isinstance(s, bytes) else str(s)
                            for s in up_r["input/sequence"][:]
                        ]
                    else:
                        n_sc = 0
                        cb_residue_idx = np.empty(0, dtype=np.int32)
                        sequence = []
                    n_cg = up_r["input/potential/compose_vector6d/elem_index"].shape[0]

                if n_sc > 0:
                    sc_grp = mh5["cg_lipid_table/cg_lipid_sc"]
                    sc_restype_order = [
                        r.decode("ascii") if isinstance(r, bytes) else str(r)
                        for r in sc_grp["restype_order"][:]
                    ]
                    # Build restype → type_index map for CG↔SC table
                    restype_to_sc_type = {r: i for i, r in enumerate(sc_restype_order)}

                    # Filter SC particles to only residues with CG↔SC table entries
                    sc_indices = []
                    sc_types = []
                    sc_ids = []
                    for cb_idx in range(n_sc):
                        res_idx = int(cb_residue_idx[cb_idx])
                        resname = sequence[res_idx] if res_idx < len(sequence) else ""
                        type_idx = restype_to_sc_type.get(resname)
                        if type_idx is not None:
                            sc_indices.append(cb_idx)
                            sc_types.append(type_idx)
                            sc_ids.append(res_idx)

                    injected_cutoff_ang = None
                    if not sc_indices:
                        print("  cg_lipid_sc: no residues with CG↔SC table entries, skipping SC↔CG interaction")
                    else:
                        sc_indices = np.array(sc_indices, dtype=np.int32)
                        sc_types = np.array(sc_types, dtype=np.int32)
                        sc_ids = np.array(sc_ids, dtype=np.int32)

                        cg_sc = pot.create_group("cg_lipid_sc")
                        cg_sc.attrs["initialized"] = True
                        cg_sc.attrs["arguments"] = np.array([
                            b"placement_fixed_point_vector_only_CB",
                            b"compose_vector6d",
                        ])
                        for attr_name, attr_value in box_attrs.items():
                            cg_sc.attrs[attr_name] = attr_value
                        if "schema" in sc_grp.attrs:
                            cg_sc.attrs["schema"] = sc_grp.attrs["schema"]
                        cg_sc.attrs["angle_convention"] = "ang1=-n1_dot_n12;ang2=n2_dot_n12"
                        cg_sc.attrs["radial_mode"] = sc_grp.attrs.get("radial_mode", "full_multimode")
                        for attr_name in ("knot_spacing_ang", "taper_width_ang", "fit_r_min_nm", "fit_r_max_nm"):
                            if attr_name in sc_grp.attrs:
                                cg_sc.attrs[attr_name] = np.float32(sc_grp.attrs[attr_name])
                        cutoff_ang = _cg_lipid_sc_runtime_cutoff_ang(sc_grp.attrs)
                        if cutoff_ang is not None:
                            cg_sc.attrs["cutoff_ang"] = np.float32(cutoff_ang)
                            injected_cutoff_ang = float(cutoff_ang)
                        for attr_name in ("n_modes", "n_radial", "n_angular"):
                            if attr_name in sc_grp.attrs:
                                cg_sc.attrs[attr_name] = np.int32(sc_grp.attrs[attr_name])

                        psi = cg_sc.create_group("pair_interaction")
                        psi.create_dataset(
                            "interaction_param",
                            data=sc_grp["interaction_param"][:].astype(np.float32),
                        )

                        psi.create_dataset("index1", data=sc_indices)
                        psi.create_dataset("type1", data=sc_types)
                        # Shift ids by 4 bits: PosQuadSplineInteraction::acceptable_id_pair
                        # accepts only pairs where (id1>>4) != (id2>>4).
                        # SC ids use residue_idx << 4 (range 0–480).
                        psi.create_dataset("id1", data=(sc_ids << 4).astype(np.int32))

                        psi.create_dataset("index2", data=np.arange(n_cg, dtype=np.int32))
                        psi.create_dataset("type2", data=np.zeros(n_cg, dtype=np.int32))
                        # CG lipid ids shifted into disjoint range (100000+cg_idx)<<4
                        # so no SC↔CG pair is rejected by the id filter.
                        psi.create_dataset("id2", data=((np.arange(n_cg, dtype=np.int32) + 100000) << 4))

                    if injected_cutoff_ang is not None:
                        print(
                            f"  Injected cg_lipid_sc: {len(sc_indices)} SC × {n_cg} CG lipids, "
                            f"cutoff={injected_cutoff_ang:.3f} Å"
                        )
                    else:
                        print(f"  Injected cg_lipid_sc: {len(sc_indices)} SC × {n_cg} CG lipids")
                else:
                    print("  cg_lipid_sc: no SC placement node found, skipping SC↔CG interaction")


def inject_stage7_sc_table_nodes(
    up_file: Path,
    martini_h5: Path,
    upside_home: Path,
    rama_library: Path,
    rama_sheet_mixing: Path,
    hbond_energy: Path,
    reference_state_rama: Path,
    protein_itp: Path | None = None,
):
    up_file = Path(up_file).expanduser().resolve()
    martini_h5 = Path(martini_h5).expanduser().resolve()
    upside_home = Path(upside_home).expanduser().resolve()
    rama_library = Path(rama_library).expanduser().resolve()
    rama_sheet_mixing = Path(rama_sheet_mixing).expanduser().resolve()
    hbond_energy = Path(hbond_energy).expanduser().resolve()
    reference_state_rama = Path(reference_state_rama).expanduser().resolve()
    protein_itp = Path(protein_itp).expanduser().resolve() if protein_itp else None
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
            "grid_nm",
            "cos_theta_grid",
            "rotamer_count",
            "rotamer_probability_fixed",
            "rotamer_radial_energy_kj_mol",
            "rotamer_angular_energy_kj_mol",
            "rotamer_angular_profile",
        ]
        missing_sc_datasets = [name for name in required_sc_datasets if name not in sc_grp]
        if missing_sc_datasets:
            missing_text = ", ".join(missing_sc_datasets)
            raise SystemExit(
                f"ERROR: {martini_h5} is missing required rotamer-resolved SC datasets: {missing_text}. "
                "Regenerate martini.h5 via build_martini_tables()."
            )
        restype_order = decode_string_array(sc_grp["restype_order"])
        target_order = decode_string_array(sc_grp["target_order"])
        grid_nm = sc_grp["grid_nm"][:].astype(np.float32)
        cos_theta_grid = sc_grp["cos_theta_grid"][:].astype(np.float32)
        rotamer_count = sc_grp["rotamer_count"][:].astype(np.int32)
        rotamer_probability_fixed = sc_grp["rotamer_probability_fixed"][:].astype(np.float32)
        rotamer_radial_energy_kj_mol = sc_grp["rotamer_radial_energy_kj_mol"][:].astype(np.float32)
        rotamer_angular_energy_kj_mol = sc_grp["rotamer_angular_energy_kj_mol"][:].astype(np.float32)
        rotamer_angular_profile = sc_grp["rotamer_angular_profile"][:].astype(np.float32)

    restype_to_index = {name: i for i, name in enumerate(restype_order)}
    target_to_index = {name: i for i, name in enumerate(target_order)}

    with h5py.File(up_file, "r+") as up:
        inp = up["input"]
        pot = inp["potential"]
        martini_potential = pot["martini_potential"]

        residue_ids, affine_atoms = build_affine_atoms(inp)
        sequence = resolve_sequence(inp, len(residue_ids), protein_itp)
        env_atom_index, env_target_index = build_env_rows(inp, target_to_index)
        rotamer_payload = load_sidechain_rotamer_payload(sidechain_lib, sequence, restype_to_index)

        if "hybrid_sc_map" in inp:
            del inp["hybrid_sc_map"]
        if "sequence" in inp:
            del inp["sequence"]
        inp.create_dataset("sequence", data=np.asarray([np.bytes_(x) for x in sequence], dtype="S3"))

        for node_name in [*LEGACY_STAGE7_NODES, *BACKBONE_NODES]:
            if node_name in pot:
                del pot[node_name]

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
        g_sc_place.create_dataset("rama_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_sc_place.create_dataset("affine_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_sc_place.create_dataset("layer_index", data=rotamer_payload["layer_index"], dtype=np.int32)
        g_sc_place.create_dataset("placement_data", data=rotamer_payload["placement_data"], dtype=np.float32)
        g_sc_place.create_dataset("beadtype_seq", data=rotamer_payload["beadtype_seq"])
        g_sc_place.create_dataset("id_seq", data=rotamer_payload["id_seq"], dtype=np.int32)
        g_sc_place.create_dataset("fix_rotamer", data=np.zeros((0, 2), dtype=np.int32), dtype=np.int32)

        g_pl = recreate_group(pot, "placement_fixed_scalar")
        g_pl.attrs["arguments"] = np.asarray([np.bytes_("affine_alignment")])
        g_pl.create_dataset("rama_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_pl.create_dataset("affine_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_pl.create_dataset("layer_index", data=rotamer_payload["layer_index"], dtype=np.int32)
        g_pl.create_dataset("placement_data", data=rotamer_payload["placement_scalar_data"], dtype=np.float32)

        g_rot = recreate_group(pot, "rotamer")
        g_rot.attrs["arguments"] = np.asarray(
            [
                np.bytes_("placement_fixed_point_vector_only"),
                np.bytes_("placement_fixed_scalar"),
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
        rama_library=rama_library,
        rama_sheet_mixing=rama_sheet_mixing,
        hbond_energy=hbond_energy,
        reference_state_rama=reference_state_rama,
        upside_home=upside_home,
    )

    with h5py.File(up_file, "r+") as up:
        inp = up["input"]
        pot = inp["potential"]
        martini_potential = pot["martini_potential"]
        g_sc = recreate_group(pot, "martini_sc_table_1body")
        g_sc.attrs["arguments"] = np.asarray([np.bytes_("pos"), np.bytes_("placement_fixed_point_vector_only_CB")])
        g_sc.attrs["energy_conversion_kj_per_eup"] = np.float32(martini_potential.attrs["energy_conversion_kj_per_eup"])
        g_sc.attrs["length_conversion_angstrom_per_nm"] = np.float32(
            martini_potential.attrs["length_conversion_angstrom_per_nm"]
        )
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
        g_sc.create_dataset("grid_nm", data=grid_nm, dtype=np.float32)
        g_sc.create_dataset("cos_theta_grid", data=cos_theta_grid, dtype=np.float32)
        g_sc.create_dataset("rotamer_count", data=rotamer_count.astype(np.int32), dtype=np.int32)
        g_sc.create_dataset("rotamer_probability_fixed", data=rotamer_probability_fixed, dtype=np.float32)
        g_sc.create_dataset("rotamer_radial_energy_kj_mol", data=rotamer_radial_energy_kj_mol, dtype=np.float32)
        g_sc.create_dataset("rotamer_angular_energy_kj_mol", data=rotamer_angular_energy_kj_mol, dtype=np.float32)
        g_sc.create_dataset("rotamer_angular_profile", data=rotamer_angular_profile, dtype=np.float32)
        g_sc.create_dataset("restype_order", data=np.asarray([np.bytes_(x) for x in restype_order], dtype="S4"))
        g_sc.create_dataset("target_order", data=np.asarray([np.bytes_(x) for x in target_order], dtype="S8"))

    print(
        f"Injected rotamer-weighted martini_sc_table_1body into {up_file}: "
        f"n_rows={len(rotamer_payload['id_seq'])} skipped={len(rotamer_payload['skipped'])} "
        f"n_env={len(env_atom_index)} n_restypes={len(restype_order)} n_targets={len(target_order)}"
    )


if __name__ == "__main__":
    raise SystemExit("Use martini_prepare_system.py as the workflow preparation entrypoint.")
