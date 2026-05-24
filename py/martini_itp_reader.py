#!/usr/bin/env python3
"""ITP readers used by the dry-MARTINI workflow."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any, Dict, List, Tuple

CANONICAL_RESIDUES = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
)


def infer_charge_from_atomtype(bead_type: str) -> float:
    """Derive nominal charge from a MARTINI atom-type name.

    Convention:
      Qda, Qd, SQda, SQd  ->  +1
      Qa, SQa              ->  -1
      everything else      ->   0
    """
    bt = bead_type.strip()
    if bt in ("Qda", "Qd", "SQda", "SQd"):
        return 1.0
    if bt in ("Qa", "SQa"):
        return -1.0
    return 0.0


def parse_itp_defines(itp_path: str | Path) -> Dict[str, Tuple[float, ...]]:
    """Extract ``#define`` macros from an ITP file.

    Returns ``{NAME: (val1, val2, ...)}``.  Bond macros have two floats
    (r0_nm, k_kj_mol_nm2); angle macros have two floats (theta0_deg, k_kj_mol).
    """
    itp_path = Path(itp_path).expanduser().resolve()
    macros: Dict[str, Tuple[float, ...]] = {}
    with itp_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            stripped = raw.split(";", 1)[0].strip()
            if not stripped.startswith("#define"):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                continue
            name = parts[1]
            vals = []
            for tok in parts[2:]:
                try:
                    vals.append(float(tok))
                except ValueError:
                    break
            if vals:
                macros[name] = tuple(vals)
    return macros


def parse_dry_forcefield(
    ff_path: str | Path,
) -> Tuple[List[str], Dict[Tuple[str, str], Dict[str, float]]]:
    """Parse ``[atomtypes]`` and ``[nonbond_params]`` from a dry-MARTINI ITP.

    Returns ``(atomtype_names, pair_params)`` where *pair_params* maps
    ``(type_i, type_j)`` to ``{"sigma_nm": ..., "epsilon_kj_mol": ...}``.
    Both orderings are stored for symmetric lookup.
    """
    ff_path = Path(ff_path).expanduser().resolve()
    macros = parse_itp_defines(ff_path)
    atomtypes: List[str] = []
    pair_params: Dict[Tuple[str, str], Dict[str, float]] = {}
    section = ""

    with ff_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            stripped = raw.split(";", 1)[0].strip()
            if not stripped:
                continue
            if stripped.startswith("#define"):
                continue
            if stripped.startswith("[") and stripped.endswith("]"):
                section = stripped[1:-1].strip().lower()
                continue
            parts = stripped.split()
            if section == "atomtypes":
                atomtypes.append(parts[0])
                continue
            if section != "nonbond_params":
                continue
            if len(parts) < 4 or parts[2] != "1":
                continue
            type_i, type_j = parts[0], parts[1]
            if len(parts) == 4:
                macro = parts[3]
                if macro not in macros:
                    raise RuntimeError(
                        f"Unknown dry-MARTINI macro '{macro}' in {ff_path}"
                    )
                sigma_nm, epsilon_kj = macros[macro]
            else:
                sigma_nm = float(parts[3])
                epsilon_kj = float(parts[4])
            payload = {"sigma_nm": sigma_nm, "epsilon_kj_mol": epsilon_kj}
            pair_params[(type_i, type_j)] = payload
            pair_params[(type_j, type_i)] = payload
    return atomtypes, pair_params


def parse_itp_atomtype_masses(ff_path: str | Path) -> Dict[str, float]:
    """Read ``[atomtypes]`` masses from an ITP force-field file.

    Returns ``{type_name: mass_g_mol}``.
    """
    ff_path = Path(ff_path).expanduser().resolve()
    masses: Dict[str, float] = {}
    in_atomtypes = False
    with ff_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("["):
                in_atomtypes = line.replace(" ", "") == "[atomtypes]"
                continue
            if not in_atomtypes:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                masses[parts[0]] = float(parts[1])
            except ValueError:
                continue
    return masses


def _empty_topology() -> Dict[str, Any]:
    return {
        "atoms": [],
        "bonds": [],
        "angles": [],
        "dihedrals": [],
        "position_restraints": [],
        "exclusions": [],
        "moleculetype": None,
        "molecules": {},
    }


def _macro_map(preprocessor_defines: Dict[str, Any] | None) -> Dict[str, List[float]]:
    macros: Dict[str, List[float]] = {}
    if preprocessor_defines is None:
        return macros
    for name, value in preprocessor_defines.items():
        if isinstance(value, (list, tuple)):
            macros[name] = [float(x) for x in value]
        else:
            macros[name] = [float(value)]
    return macros


def _macro_value(token: str, macros: Dict[str, List[float]], index: int | None = None) -> float:
    try:
        return float(token)
    except ValueError:
        pass
    values = macros.get(token)
    if not values:
        raise ValueError(f"Unknown macro '{token}'")
    if index is None:
        return values[0]
    return values[min(index, len(values) - 1)]


def parse_itp_file(
    itp_file: str | Path,
    target_molecule: str | None = None,
    preprocessor_defines: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """Read atoms, bonded terms, restraints, and exclusions from a MARTINI ITP."""
    itp_path = Path(itp_file).expanduser()
    if not itp_path.exists():
        print(f"Warning: ITP file {itp_file} not found")
        return _empty_topology()

    topology = _empty_topology()
    macros = _macro_map(preprocessor_defines)
    pp_stack: List[Tuple[bool, bool]] = []
    current_active = True
    current_section = ""
    current_molecule = None
    current_mol_data = None

    with itp_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.split(";", 1)[0].strip()
            if not line:
                continue

            if line.startswith("#ifdef") or line.startswith("#ifndef"):
                parts = line.split()
                name = parts[1] if len(parts) >= 2 else ""
                is_defined = name in macros
                cond = is_defined if line.startswith("#ifdef") else not is_defined
                pp_stack.append((current_active, cond))
                current_active = current_active and cond
                continue
            if line.startswith("#else"):
                if pp_stack:
                    parent_active, cond = pp_stack[-1]
                    cond = not cond
                    pp_stack[-1] = (parent_active, cond)
                    current_active = parent_active and cond
                continue
            if line.startswith("#endif"):
                if pp_stack:
                    parent_active, _ = pp_stack.pop()
                    current_active = parent_active
                continue
            if not current_active:
                continue

            if line.startswith("#define"):
                parts = line.split()
                if len(parts) >= 3:
                    values = []
                    for token in parts[2:]:
                        try:
                            values.append(float(token))
                        except ValueError:
                            break
                    if values:
                        macros[parts[1]] = values
                continue
            if line.startswith("#"):
                continue

            if line.startswith("[") and line.endswith("]"):
                current_section = line[1:-1].strip().lower()
                continue

            parts = line.split()
            if current_section == "moleculetype":
                if not parts:
                    continue
                current_molecule = parts[0]
                if topology["moleculetype"] is None:
                    topology["moleculetype"] = current_molecule
                current_mol_data = {
                    "atoms": [],
                    "bonds": [],
                    "angles": [],
                    "dihedrals": [],
                    "position_restraints": [],
                    "exclusions": [],
                }
                topology["molecules"][current_molecule] = current_mol_data
                continue

            if current_mol_data is None:
                continue

            if current_section == "atoms" and len(parts) >= 6:
                try:
                    atom = {
                        "id": int(parts[0]),
                        "type": parts[1],
                        "resnr": int(parts[2]),
                        "residue": parts[3],
                        "atom": parts[4],
                        "cgnr": int(parts[5]),
                        "charge": 0.0,
                        "mass": 0.0,
                    }
                    if len(parts) >= 7:
                        atom["charge"] = float(parts[6])
                    if len(parts) >= 8:
                        atom["mass"] = float(parts[7])
                    current_mol_data["atoms"].append(atom)
                    topology["atoms"].append(atom)
                except (ValueError, IndexError):
                    continue

            elif current_section == "bonds" and len(parts) >= 3:
                try:
                    bond = {
                        "i": int(parts[0]) - 1,
                        "j": int(parts[1]) - 1,
                        "func": int(parts[2]),
                        "r0": 0.0,
                        "k": 0.0,
                    }
                    if len(parts) >= 5:
                        bond["r0"] = _macro_value(parts[3], macros, 0)
                        bond["k"] = _macro_value(parts[4], macros, 1)
                    elif len(parts) >= 4:
                        bond["r0"] = _macro_value(parts[3], macros, 0)
                        bond["k"] = _macro_value(parts[3], macros, 1)
                    current_mol_data["bonds"].append(bond)
                    topology["bonds"].append(bond)
                except (ValueError, IndexError):
                    continue

            elif current_section == "angles" and len(parts) >= 5:
                try:
                    angle = {
                        "i": int(parts[0]) - 1,
                        "j": int(parts[1]) - 1,
                        "k": int(parts[2]) - 1,
                        "func": int(parts[3]),
                        "theta0": 0.0,
                        "force_k": 0.0,
                    }
                    if len(parts) >= 6:
                        angle["theta0"] = _macro_value(parts[4], macros, 0)
                        angle["force_k"] = _macro_value(parts[5], macros, 1)
                    else:
                        angle["theta0"] = _macro_value(parts[4], macros, 0)
                        angle["force_k"] = _macro_value(parts[4], macros, 1)
                    current_mol_data["angles"].append(angle)
                    topology["angles"].append(angle)
                except (ValueError, IndexError):
                    continue

            elif current_section == "dihedrals" and len(parts) >= 5:
                try:
                    dihedral = {
                        "i": int(parts[0]) - 1,
                        "j": int(parts[1]) - 1,
                        "k": int(parts[2]) - 1,
                        "l": int(parts[3]) - 1,
                        "func": int(parts[4]),
                        "phi0": float(parts[5]) if len(parts) >= 7 else 0.0,
                        "k": float(parts[6]) if len(parts) >= 7 else 0.0,
                        "mult": int(parts[7]) if len(parts) >= 8 else 1,
                    }
                    current_mol_data["dihedrals"].append(dihedral)
                    topology["dihedrals"].append(dihedral)
                except (ValueError, IndexError):
                    continue

            elif current_section == "position_restraints" and len(parts) >= 5:
                try:
                    restraint = {
                        "i": int(parts[0]) - 1,
                        "func": int(parts[1]),
                        "fx": _macro_value(parts[2], macros),
                        "fy": _macro_value(parts[3], macros),
                        "fz": _macro_value(parts[4], macros),
                    }
                    current_mol_data["position_restraints"].append(restraint)
                    topology["position_restraints"].append(restraint)
                except (ValueError, IndexError):
                    continue

            elif current_section == "exclusions" and len(parts) >= 2:
                try:
                    atoms = [int(part) - 1 for part in parts]
                    for i, atom_i in enumerate(atoms):
                        for atom_j in atoms[i + 1:]:
                            exclusion = (atom_i, atom_j)
                            current_mol_data["exclusions"].append(exclusion)
                            topology["exclusions"].append(exclusion)
                except ValueError:
                    continue

    if target_molecule and target_molecule in topology["molecules"]:
        mol_data = topology["molecules"][target_molecule]
        return {
            "atoms": mol_data["atoms"],
            "bonds": mol_data["bonds"],
            "angles": mol_data["angles"],
            "dihedrals": mol_data["dihedrals"],
            "position_restraints": mol_data["position_restraints"],
            "exclusions": mol_data["exclusions"],
            "moleculetype": target_molecule,
            "molecules": {target_molecule: mol_data},
        }

    return topology


def parse_dopc_from_itp(
    lipids_itp_path: str | Path,
) -> dict:
    """Parse the DOPC molecule definition from the lipids ITP.

    Returns a dict with keys:
      ``atom_names`` (list[str]), ``bead_types`` (list[str]),
      ``bead_charges`` (list[float]), ``bonds`` (list[tuple]),
      ``angles`` (list[tuple]).
    """
    topology = parse_itp_file(lipids_itp_path, "DOPC")
    dopc_atoms = topology["atoms"]
    dopc_bonds = [
        (bond["i"], bond["j"], bond["r0"], bond["k"])
        for bond in topology["bonds"]
    ]
    dopc_angles = [
        (angle["i"], angle["j"], angle["k"], angle["theta0"], angle["force_k"])
        for angle in topology["angles"]
    ]

    return {
        "atom_names": [a["atom"] for a in dopc_atoms],
        "bead_types": [a["type"] for a in dopc_atoms],
        "bead_charges": [a["charge"] for a in dopc_atoms],
        "bonds": dopc_bonds,
        "angles": dopc_angles,
    }


def load_martini_forcefield(
    martinize_path: str | Path,
    forcefield_name: str = "martini22",
) -> Dict[str, List[str]]:
    """Import *martinize_path* and return the residue-to-bead-type mapping.

    Returns ``{RESNAME: [bead_type, ...]}`` for each canonical residue.
    """
    martinize_path = Path(martinize_path).expanduser().resolve()
    spec = importlib.util.spec_from_file_location(
        "sc_training_martinize_runtime", martinize_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load martinize module from {martinize_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, forcefield_name):
        raise RuntimeError(
            f"Forcefield '{forcefield_name}' not found in {martinize_path}"
        )
    ff = getattr(module, forcefield_name)()
    residue_map: Dict[str, List[str]] = {}
    for residue in CANONICAL_RESIDUES:
        raw = ff.sidechains.get(residue, [])
        if not raw:
            residue_map[residue] = []
            continue
        bead_tokens = [str(tok).strip() for tok in raw[0] if str(tok).strip()]
        residue_map[residue] = [tok for tok in bead_tokens if tok != "D"]
    return residue_map
