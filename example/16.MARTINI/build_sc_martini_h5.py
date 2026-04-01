#!/usr/bin/env python3

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


DEFAULT_SC_TABLE_JSON = Path(
    "/Users/yinhan/Documents/SC-training/runs/default/results/assembled/sc_table.json"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a native-unit martini.h5 SC table from assembled SC-training JSON."
    )
    parser.add_argument(
        "--sc-table-json",
        default=str(DEFAULT_SC_TABLE_JSON),
        help="Path to assembled sc_table.json from SC-training.",
    )
    parser.add_argument(
        "--output-h5",
        default="parameters/ff_2.1/martini.h5",
        help="Output HDF5 file path.",
    )
    return parser.parse_args()


def require_h5import():
    exe = shutil.which("h5import")
    if not exe:
        raise SystemExit("ERROR: h5import was not found in PATH")
    return exe


def read_sc_table(path: Path):
    if not path.exists():
        raise SystemExit(f"ERROR: SC table JSON not found: {path}")
    return json.loads(path.read_text())


def normalize_string(value):
    return "" if value is None else str(value)


def build_dense_table(sc_table):
    tables_by_residue = sc_table.get("tables_by_residue")
    if not isinstance(tables_by_residue, dict) or not tables_by_residue:
        raise SystemExit("ERROR: sc_table.json is missing non-empty tables_by_residue")

    residues = list(tables_by_residue)
    first_residue = residues[0]
    targets = list(tables_by_residue[first_residue])
    if not targets:
        raise SystemExit("ERROR: first residue has no target tables")

    reference_grid = np.asarray(
        tables_by_residue[first_residue][targets[0]]["grid_nm"], dtype=np.float32
    )
    if reference_grid.ndim != 1 or reference_grid.size < 2:
        raise SystemExit("ERROR: invalid radial grid in sc_table.json")

    energy = np.zeros((len(residues), len(targets), reference_grid.size), dtype=np.float32)

    for ri, residue in enumerate(residues):
        residue_tables = tables_by_residue[residue]
        missing_targets = sorted(set(targets) - set(residue_tables))
        extra_targets = sorted(set(residue_tables) - set(targets))
        if missing_targets or extra_targets:
            raise SystemExit(
                f"ERROR: target mismatch for residue {residue}: "
                f"missing={missing_targets} extra={extra_targets}"
            )

        for ti, target in enumerate(targets):
            entry = residue_tables[target]
            grid = np.asarray(entry["grid_nm"], dtype=np.float32)
            if grid.shape != reference_grid.shape or not np.allclose(grid, reference_grid):
                raise SystemExit(
                    f"ERROR: non-uniform radial grid for residue {residue} target {target}"
                )
            energy_row = np.asarray(entry["energy_kj_mol"], dtype=np.float32)
            if energy_row.shape != reference_grid.shape:
                raise SystemExit(
                    f"ERROR: energy grid shape mismatch for residue {residue} target {target}: "
                    f"{energy_row.shape} vs {reference_grid.shape}"
                )
            energy[ri, ti, :] = energy_row

    return residues, targets, reference_grid, energy


def write_text_dataset(
    h5import_exe: str,
    tmpdir: Path,
    output_h5: Path,
    dataset_path: str,
    values,
):
    txt_path = tmpdir / f"{dataset_path.strip('/').replace('/', '__')}.txt"
    cfg_path = tmpdir / f"{dataset_path.strip('/').replace('/', '__')}.cfg"
    txt_path.write_text("".join(f"{normalize_string(v)}\n" for v in values))
    cfg_path.write_text(f"PATH {dataset_path}\nINPUT-CLASS STR\n")
    subprocess.run(
        [h5import_exe, str(txt_path), "-c", str(cfg_path), "-o", str(output_h5)],
        check=True,
    )


def write_float_dataset(
    h5import_exe: str,
    tmpdir: Path,
    output_h5: Path,
    dataset_path: str,
    array: np.ndarray,
):
    arr = np.asarray(array, dtype=np.float32)
    bin_path = tmpdir / f"{dataset_path.strip('/').replace('/', '__')}.bin"
    cfg_path = tmpdir / f"{dataset_path.strip('/').replace('/', '__')}.cfg"
    arr.tofile(bin_path)
    dims = " ".join(str(x) for x in arr.shape)
    cfg_path.write_text(
        "\n".join(
            [
                f"PATH {dataset_path}",
                "INPUT-CLASS FP",
                "INPUT-SIZE 32",
                "INPUT-BYTE-ORDER LE",
                f"RANK {arr.ndim}",
                f"DIMENSION-SIZES {dims}",
                "OUTPUT-CLASS FP",
                "OUTPUT-SIZE 32",
                "OUTPUT-ARCHITECTURE NATIVE",
                "",
            ]
        )
    )
    subprocess.run(
        [h5import_exe, str(bin_path), "-c", str(cfg_path), "-o", str(output_h5)],
        check=True,
    )


def main():
    args = parse_args()
    h5import_exe = require_h5import()
    sc_table_path = Path(args.sc_table_json).expanduser().resolve()
    output_h5 = Path(args.output_h5).expanduser().resolve()
    output_h5.parent.mkdir(parents=True, exist_ok=True)
    if output_h5.exists():
        output_h5.unlink()

    sc_table = read_sc_table(sc_table_path)
    residues, targets, grid_nm, energy_kj_mol = build_dense_table(sc_table)

    with tempfile.TemporaryDirectory(prefix="build_sc_martini_h5.") as tmp:
        tmpdir = Path(tmp)
        write_text_dataset(h5import_exe, tmpdir, output_h5, "/schema", [sc_table.get("schema", "")])
        write_text_dataset(
            h5import_exe,
            tmpdir,
            output_h5,
            "/forcefield_name",
            [sc_table.get("forcefield_name", "")],
        )
        write_text_dataset(
            h5import_exe,
            tmpdir,
            output_h5,
            "/created_at_utc",
            [sc_table.get("created_at_utc", "")],
        )
        write_text_dataset(
            h5import_exe,
            tmpdir,
            output_h5,
            "/manifest_path",
            [sc_table.get("manifest_path", "")],
        )
        write_text_dataset(h5import_exe, tmpdir, output_h5, "/restype_order", residues)
        write_text_dataset(h5import_exe, tmpdir, output_h5, "/target_order", targets)
        write_float_dataset(h5import_exe, tmpdir, output_h5, "/grid_nm", grid_nm)
        write_float_dataset(
            h5import_exe,
            tmpdir,
            output_h5,
            "/energy_kj_mol",
            energy_kj_mol,
        )

    print(
        f"Built {output_h5} from {sc_table_path} "
        f"with {len(residues)} residues, {len(targets)} targets, {grid_nm.size} radial points"
    )


if __name__ == "__main__":
    main()
