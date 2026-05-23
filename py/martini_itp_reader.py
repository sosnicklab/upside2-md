#!/usr/bin/env python3
"""Pure ITP file parser for dry-MARTINI force field parameters.

No physics constants or hardcoded bead types live here.  Every value is read
from the .itp files on disk.  This is the single source of truth for all
ITP-derived data consumed by the table builders and workflow scripts.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Dict, List, Tuple

CANONICAL_RESIDUES = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
)


# ---------------------------------------------------------------------------
# Charge inference (follows MARTINI naming convention — no hardcoded list)
# ---------------------------------------------------------------------------

def infer_charge_from_atomtype(bead_type: str) -> float:
    """Derive nominal charge from a MARTINI atom-type name.

    Convention:
      Qda, Qd, SQda, SQd  →  +1
      Qa, SQa              →  -1
      everything else      →   0
    """
    bt = bead_type.strip()
    if bt in ("Qda", "Qd", "SQda", "SQd"):
        return 1.0
    if bt in ("Qa", "SQa"):
        return -1.0
    return 0.0


# ---------------------------------------------------------------------------
# #define macro parsing
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Force-field parsing (atomtypes + nonbond_params)
# ---------------------------------------------------------------------------

def parse_dry_forcefield(
    ff_path: str | Path,
) -> Tuple[List[str], Dict[Tuple[str, str], Dict[str, float]]]:
    """Parse ``[atomtypes]`` and ``[nonbond_params]`` from a dry-MARTINI ITP.

    Returns ``(atomtype_names, pair_params)`` where *pair_params* maps
    ``(type_i, type_j)`` → ``{"sigma_nm": ..., "epsilon_kj_mol": ...}``.
    Both orderings are stored for symmetric lookup.
    """
    ff_path = Path(ff_path).expanduser().resolve()
    macros: Dict[str, Tuple[float, float]] = {}
    atomtypes: List[str] = []
    pair_params: Dict[Tuple[str, str], Dict[str, float]] = {}
    section = ""

    with ff_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            stripped = raw.split(";", 1)[0].strip()
            if not stripped:
                continue
            if stripped.startswith("#define"):
                parts = stripped.split()
                if len(parts) == 4:
                    try:
                        macros[parts[1]] = (float(parts[2]), float(parts[3]))
                    except ValueError:
                        pass
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


# ---------------------------------------------------------------------------
# Atom-type masses
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Molecule-type parsing (reads [atoms], [bonds], [angles] for one molecule)
# ---------------------------------------------------------------------------

def _skip_preprocessor_block(lines: List[str], cursor: int, directive: str) -> int:
    """Advance *cursor* past a ``#ifdef`` / ``#ifndef`` ... ``#endif`` block.

    Only skips when the directive evaluates to False (i.e. the block is inactive).
    Nested conditionals are tracked via a depth counter.
    """
    depth = 1
    cursor += 1
    while cursor < len(lines) and depth > 0:
        stripped = lines[cursor].split(";", 1)[0].strip()
        if stripped.startswith("#ifdef") or stripped.startswith("#ifndef"):
            depth += 1
        elif stripped.startswith("#endif"):
            depth -= 1
        elif stripped.startswith("#else") and depth == 1:
            # When we hit #else at the level we're skipping, switch to active.
            return cursor + 1  # resume reading from after #else
        cursor += 1
    return cursor


def _read_itp_lines(itp_path: Path) -> List[str]:
    """Read an ITP file, returning non-empty, non-comment lines."""
    with itp_path.open("r", encoding="utf-8", errors="ignore") as fh:
        return [raw.split(";", 1)[0].strip() for raw in fh
                if raw.split(";", 1)[0].strip()]


def parse_molecule_atoms(
    itp_path: str | Path,
    molecule_name: str,
    *,
    defined_symbols: frozenset = frozenset(),
) -> List[dict]:
    """Parse the ``[atoms]`` block of *molecule_name* from a lipids ITP.

    Handles ``#ifndef`` / ``#ifdef`` / ``#else`` / ``#endif`` conditional
    blocks.  When *defined_symbols* is non-empty, ``#ifdef X`` is True iff
    ``X in defined_symbols`` and ``#ifndef X`` is True iff
    ``X not in defined_symbols``.

    Returns a list of dicts with keys:
      ``atom_idx`` (int), ``type`` (str), ``residue`` (str),
      ``atom_name`` (str), ``charge_group`` (int), ``charge`` (float).
    """
    itp_path = Path(itp_path).expanduser().resolve()
    lines = _read_itp_lines(itp_path)

    # Locate the [ moleculetype ] block
    mol_start = None
    for idx, line in enumerate(lines):
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip()
            if section == "moleculetype":
                # Next non-comment line is the molecule name + nrexcl
                for j in range(idx + 1, len(lines)):
                    if lines[j]:
                        mol_name_candidate = lines[j].split()[0]
                        if mol_name_candidate == molecule_name:
                            mol_start = j + 1
                            break
                        else:
                            break  # wrong molecule, keep scanning
                if mol_start is not None:
                    break

    if mol_start is None:
        raise RuntimeError(
            f"Molecule type '{molecule_name}' not found in {itp_path}"
        )

    # Walk through the molecule's sections, respecting preprocessor conditionals
    atoms: List[dict] = []
    current_section = ""
    i = mol_start

    while i < len(lines):
        line = lines[i]

        # Handle preprocessor directives
        if line.startswith("#ifdef "):
            symbol = line.split()[1]
            if symbol not in defined_symbols:
                i = _skip_preprocessor_block(lines, i, "#ifdef")
                continue
            i += 1
            continue
        if line.startswith("#ifndef "):
            symbol = line.split()[1]
            if symbol in defined_symbols:
                i = _skip_preprocessor_block(lines, i, "#ifndef")
                continue
            i += 1
            continue
        if line.startswith("#else"):
            # We're in the active branch; skip to #endif
            depth = 1
            i += 1
            while i < len(lines) and depth > 0:
                stripped = lines[i].split(";", 1)[0].strip()
                if stripped.startswith("#ifdef") or stripped.startswith("#ifndef"):
                    depth += 1
                elif stripped.startswith("#endif"):
                    depth -= 1
                i += 1
            continue
        if line.startswith("#endif"):
            i += 1
            continue

        # Section headers
        if line.startswith("[") and line.endswith("]"):
            current_section = line[1:-1].strip()
            i += 1
            continue

        if current_section == "atoms":
            parts = line.split()
            if len(parts) >= 7:
                atoms.append({
                    "atom_idx": int(parts[0]),
                    "type": parts[1],
                    "residue": parts[3],
                    "atom_name": parts[4],
                    "charge_group": int(parts[5]),
                    "charge": float(parts[6]),
                })
        elif current_section != "" and current_section != "atoms":
            # We've moved past [atoms] to another section of this molecule.
            # Stop when we reach the next [ moleculetype ] or end of relevant content
            if current_section in ("bonds", "angles", "dihedrals", "constraints",
                                   "pairs", "exclusions", "position_restraints",
                                   "virtual_sitesn", "virtual_sites2",
                                   "virtual_sites3", "virtual_sites4"):
                # Still within the same molecule
                pass
            else:
                break

        i += 1

    return atoms


def parse_molecule_bonds(
    itp_path: str | Path,
    molecule_name: str,
    macros: Dict[str, Tuple[float, ...]] | None = None,
    *,
    defined_symbols: frozenset = frozenset(),
) -> List[Tuple[int, int, float, float]]:
    """Parse the ``[bonds]`` block of *molecule_name*.

    Returns ``[(i_0based, j_0based, r0_nm, k_kj_mol_nm2), ...]``.
    Bond macro aliases (e.g. ``mb_np``) are resolved via *macros*.
    """
    if macros is None:
        macros = parse_itp_defines(itp_path)

    itp_path = Path(itp_path).expanduser().resolve()
    lines = _read_itp_lines(itp_path)

    # Locate the molecule
    mol_start = None
    for idx, line in enumerate(lines):
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip()
            if section == "moleculetype":
                for j in range(idx + 1, len(lines)):
                    if lines[j]:
                        mol_name_candidate = lines[j].split()[0]
                        if mol_name_candidate == molecule_name:
                            mol_start = j + 1
                            break
                        else:
                            break
                if mol_start is not None:
                    break

    if mol_start is None:
        raise RuntimeError(
            f"Molecule type '{molecule_name}' not found in {itp_path}"
        )

    bonds: List[Tuple[int, int, float, float]] = []
    current_section = ""
    i = mol_start

    while i < len(lines):
        line = lines[i]

        if line.startswith("#ifdef "):
            symbol = line.split()[1]
            if symbol not in defined_symbols:
                i = _skip_preprocessor_block(lines, i, "#ifdef")
                continue
            i += 1
            continue
        if line.startswith("#ifndef "):
            symbol = line.split()[1]
            if symbol in defined_symbols:
                i = _skip_preprocessor_block(lines, i, "#ifndef")
                continue
            i += 1
            continue
        if line.startswith("#else"):
            depth = 1
            i += 1
            while i < len(lines) and depth > 0:
                stripped = lines[i].split(";", 1)[0].strip()
                if stripped.startswith("#ifdef") or stripped.startswith("#ifndef"):
                    depth += 1
                elif stripped.startswith("#endif"):
                    depth -= 1
                i += 1
            continue
        if line.startswith("#endif"):
            i += 1
            continue

        if line.startswith("[") and line.endswith("]"):
            current_section = line[1:-1].strip()
            if current_section not in ("atoms", "bonds", "angles", "dihedrals",
                                        "constraints", "pairs", "exclusions",
                                        "position_restraints", "virtual_sitesn",
                                        "virtual_sites2", "virtual_sites3",
                                        "virtual_sites4"):
                break
            i += 1
            continue

        if current_section == "bonds":
            parts = line.split()
            if len(parts) >= 4:
                i_idx = int(parts[0]) - 1
                j_idx = int(parts[1]) - 1
                alias = parts[3]
                if alias in macros:
                    r0, k = macros[alias]
                else:
                    r0 = float(parts[3])
                    k = float(parts[4]) if len(parts) >= 5 else 0.0
                bonds.append((i_idx, j_idx, r0, k))

        i += 1

    return bonds


def parse_molecule_angles(
    itp_path: str | Path,
    molecule_name: str,
    macros: Dict[str, Tuple[float, ...]] | None = None,
    *,
    defined_symbols: frozenset = frozenset(),
) -> List[Tuple[int, int, int, float, float]]:
    """Parse the ``[angles]`` block of *molecule_name*.

    Returns ``[(i_0based, j_0based, k_0based, theta0_deg, k_kj_mol), ...]``.
    Angle macro aliases (e.g. ``ma_pgg``) are resolved via *macros*.
    """
    if macros is None:
        macros = parse_itp_defines(itp_path)

    itp_path = Path(itp_path).expanduser().resolve()
    lines = _read_itp_lines(itp_path)

    mol_start = None
    for idx, line in enumerate(lines):
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip()
            if section == "moleculetype":
                for j in range(idx + 1, len(lines)):
                    if lines[j]:
                        mol_name_candidate = lines[j].split()[0]
                        if mol_name_candidate == molecule_name:
                            mol_start = j + 1
                            break
                        else:
                            break
                if mol_start is not None:
                    break

    if mol_start is None:
        raise RuntimeError(
            f"Molecule type '{molecule_name}' not found in {itp_path}"
        )

    angles: List[Tuple[int, int, int, float, float]] = []
    current_section = ""
    i = mol_start

    while i < len(lines):
        line = lines[i]

        if line.startswith("#ifdef "):
            symbol = line.split()[1]
            if symbol not in defined_symbols:
                i = _skip_preprocessor_block(lines, i, "#ifdef")
                continue
            i += 1
            continue
        if line.startswith("#ifndef "):
            symbol = line.split()[1]
            if symbol in defined_symbols:
                i = _skip_preprocessor_block(lines, i, "#ifndef")
                continue
            i += 1
            continue
        if line.startswith("#else"):
            depth = 1
            i += 1
            while i < len(lines) and depth > 0:
                stripped = lines[i].split(";", 1)[0].strip()
                if stripped.startswith("#ifdef") or stripped.startswith("#ifndef"):
                    depth += 1
                elif stripped.startswith("#endif"):
                    depth -= 1
                i += 1
            continue
        if line.startswith("#endif"):
            i += 1
            continue

        if line.startswith("[") and line.endswith("]"):
            current_section = line[1:-1].strip()
            if current_section not in ("atoms", "bonds", "angles", "dihedrals",
                                        "constraints", "pairs", "exclusions",
                                        "position_restraints", "virtual_sitesn",
                                        "virtual_sites2", "virtual_sites3",
                                        "virtual_sites4"):
                break
            i += 1
            continue

        if current_section == "angles":
            parts = line.split()
            if len(parts) >= 5:
                i_idx = int(parts[0]) - 1
                j_idx = int(parts[1]) - 1
                k_idx = int(parts[2]) - 1
                alias = parts[4]
                if alias in macros:
                    theta0, k = macros[alias]
                else:
                    theta0 = float(parts[3])
                    k = float(parts[4])
                angles.append((i_idx, j_idx, k_idx, theta0, k))

        i += 1

    return angles


# ---------------------------------------------------------------------------
# Convenience: extract full DOPC parameter set from ITP files
# ---------------------------------------------------------------------------

def parse_dopc_from_itp(
    lipids_itp_path: str | Path,
) -> dict:
    """Parse the DOPC molecule definition from the lipids ITP.

    Returns a dict with keys:
      ``atom_names`` (list[str]), ``bead_types`` (list[str]),
      ``bead_charges`` (list[float]), ``bonds`` (list[tuple]),
      ``angles`` (list[tuple]).
    """
    lipids_itp_path = Path(lipids_itp_path).expanduser().resolve()
    macros = parse_itp_defines(lipids_itp_path)

    dopc_atoms = parse_molecule_atoms(lipids_itp_path, "DOPC")
    dopc_bonds = parse_molecule_bonds(lipids_itp_path, "DOPC", macros)
    dopc_angles = parse_molecule_angles(lipids_itp_path, "DOPC", macros)

    return {
        "atom_names": [a["atom_name"] for a in dopc_atoms],
        "bead_types": [a["type"] for a in dopc_atoms],
        "bead_charges": [a["charge"] for a in dopc_atoms],
        "bonds": dopc_bonds,
        "angles": dopc_angles,
    }


# ---------------------------------------------------------------------------
# Martinize forcefield loader (sidechain bead type mappings)
# ---------------------------------------------------------------------------

def load_martini_forcefield(
    martinize_path: str | Path,
    forcefield_name: str = "martini22",
) -> Dict[str, List[str]]:
    """Import *martinize_path* and return the residue→bead_type mapping.

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
