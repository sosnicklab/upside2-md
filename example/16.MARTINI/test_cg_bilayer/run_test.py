#!/usr/bin/env python3
"""Small bilayer-only test for single-particle DOPC lipids."""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "py"))

TEST_DIR = Path(__file__).resolve().parent
RUN_DIR = TEST_DIR / "run"
SOURCE_PDB = REPO_ROOT / "parameters" / "dryMARTINI" / "DOPC.pdb"
RUNTIME_PDB = RUN_DIR / "dopc_only.pdb"
PDB_ID = "test_cg_bilayer"
MARTINI_H5 = RUN_DIR / "martini.h5"
UP_FILE = RUN_DIR / "test.input.up"
SIM_FILE = RUN_DIR / f"{PDB_ID}.dyn.up"
UPSIDE_BIN = REPO_ROOT / "obj" / "upside"

ENERGY_CONVERSION_KJ_PER_EUP = "2.914952774272"
LENGTH_CONVERSION_ANG_PER_NM = "10.0"


def _is_atom_line(line: str) -> bool:
    return line.startswith(("ATOM  ", "HETATM"))


def _resname(line: str) -> str:
    return line[17:21].strip().upper()


def _atom_name(line: str) -> str:
    return line[12:16].strip().upper()


def _coords_ang(line: str) -> np.ndarray:
    return np.array(
        [float(line[30:38]), float(line[38:46]), float(line[46:54])],
        dtype=np.float64,
    )


def write_dopc_only_pdb() -> Path:
    """Strip ions/water from the dry-MARTINI DOPC reference PDB."""
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    kept = []
    cryst1 = None
    with SOURCE_PDB.open() as f:
        for line in f:
            if line.startswith("CRYST1"):
                cryst1 = line
            if _is_atom_line(line) and _resname(line) in {"DOP", "DOPC"}:
                kept.append(line)

    if not kept:
        raise RuntimeError(f"No DOPC atoms found in {SOURCE_PDB}")

    with RUNTIME_PDB.open("w") as f:
        if cryst1 is not None:
            f.write(cryst1)
        for line in kept:
            f.write(line)
        f.write("END\n")

    print(f"Wrote {RUNTIME_PDB} with {len(kept)} DOPC beads")
    return RUNTIME_PDB


def _read_first_dopc_reference() -> tuple[np.ndarray, list[str]]:
    atoms = []
    with RUNTIME_PDB.open() as f:
        for line in f:
            if not (_is_atom_line(line) and _resname(line) in {"DOP", "DOPC"}):
                continue
            atoms.append((_atom_name(line), _coords_ang(line)))

    if len(atoms) < 14:
        raise RuntimeError("DOPC-only PDB does not contain a complete lipid")

    first = atoms[:14]
    coords = np.array([coord for _, coord in first], dtype=np.float64)
    com = coords.mean(axis=0)
    ref_bead_positions_nm = (coords - com) * 0.1

    from martini_build_tables import _DOPC_BEAD_TYPES

    return ref_bead_positions_nm, list(_DOPC_BEAD_TYPES)


def build_tables(resolution: str) -> None:
    from martini_build_tables import build_martini_tables

    ref_bead_positions_nm, bead_types = _read_first_dopc_reference()
    dry_ff = REPO_ROOT / "parameters" / "dryMARTINI" / "dry_martini_v2.1.itp"
    martinize = REPO_ROOT / "py" / "martinize.py"
    sc_lib = REPO_ROOT / "parameters" / "ff_2.1" / "sidechain.h5"

    if resolution == "coarse":
        r_count = 40
        direction_count = 8
        cos_theta_count = 7
    else:
        r_count = 96
        direction_count = 24
        cos_theta_count = 13

    print(f"Building MARTINI tables ({resolution})")
    build_martini_tables(
        output_path=MARTINI_H5,
        dry_ff_path=dry_ff,
        martinize_path=martinize,
        sidechain_lib_path=sc_lib,
        forcefield_name="martini22",
        active_residue_names=[],
        active_atom_types=set(bead_types),
        r_count=r_count,
        direction_count=direction_count,
        cos_theta_count=cos_theta_count,
        cg_lipid_config={
            "ref_bead_positions_nm": ref_bead_positions_nm,
            "bead_types": bead_types,
        },
    )


def convert_and_inject() -> None:
    from martini_prepare_system_lib import (
        convert_stage,
        inject_cg_lipid_nodes,
        inject_particles_table,
        write_stage_debug_outputs,
    )

    os.environ["UPSIDE_RUNTIME_PDB_FILE"] = str(RUNTIME_PDB.resolve())
    os.environ["UPSIDE_SIMULATION_STAGE"] = "minimization"
    os.environ["UPSIDE_MARTINI_ENERGY_CONVERSION"] = ENERGY_CONVERSION_KJ_PER_EUP
    os.environ["UPSIDE_MARTINI_LENGTH_CONVERSION"] = LENGTH_CONVERSION_ANG_PER_NM

    convert_stage(pdb_id=PDB_ID, stage="minimization", run_dir=str(RUN_DIR))
    inject_particles_table(UP_FILE, MARTINI_H5)
    inject_cg_lipid_nodes(UP_FILE, MARTINI_H5)
    write_stage_debug_outputs(UP_FILE, debug_dir=RUN_DIR / "debug")


def _normalize_pos(pos: np.ndarray) -> np.ndarray:
    arr = np.asarray(pos)
    if arr.ndim == 4 and arr.shape[1] == 1:
        arr = arr[:, 0, :, :]
    if arr.ndim == 3 and arr.shape[-1] == 1:
        arr = arr[:, :, 0]
    if arr.ndim == 3 and arr.shape[1] == 3:
        arr = np.transpose(arr, (0, 2, 1))
    return arr


def _nearest_xy(xy: np.ndarray, box_xy: tuple[float, float]) -> np.ndarray:
    n = xy.shape[0]
    out = np.full(n, np.nan, dtype=np.float64)
    box = np.array(box_xy, dtype=np.float64)
    for i in range(n):
        delta = xy - xy[i]
        delta -= box * np.round(delta / box)
        d = np.sqrt(np.sum(delta * delta, axis=1))
        d[i] = np.inf
        out[i] = d.min()
    return out


def _directions_from_sites(
    pos: np.ndarray,
    cgl_index: np.ndarray,
    orientation_index: np.ndarray | None,
    stored_dirs: np.ndarray,
    box: tuple[float, float, float],
) -> np.ndarray:
    if orientation_index is None:
        return stored_dirs.astype(np.float64)
    dr = pos[orientation_index] - pos[cgl_index]
    box_arr = np.array(box, dtype=np.float64)
    dr -= box_arr * np.round(dr / box_arr)
    norm = np.linalg.norm(dr, axis=1)
    dirs = stored_dirs.astype(np.float64).copy()
    mask = norm > 1e-12
    dirs[mask] = dr[mask] / norm[mask, None]
    return dirs


def _load_cgl_geometry(up_file: Path) -> tuple[np.ndarray, np.ndarray, tuple[float, float, float]]:
    with h5py.File(up_file, "r") as h5:
        pos = _normalize_pos(h5["input/pos"][:])
        cv = h5["input/potential/compose_vector6d"]
        cgl_index = cv["elem_index"][:].astype(np.int64)
        stored_dirs = cv["direction"][:].astype(np.float64)
        orientation_index = (
            cv["orientation_index"][:].astype(np.int64)
            if "orientation_index" in cv
            else None
        )
        pot = h5["input/potential/martini_potential"]
        box = (
            float(pot.attrs["x_len"]),
            float(pot.attrs["y_len"]),
            float(pot.attrs["z_len"]),
        )
    dirs = _directions_from_sites(pos, cgl_index, orientation_index, stored_dirs, box)
    return pos[cgl_index], dirs, box


def report_initial_geometry(up_file: Path) -> None:
    cgl_pos, dirs, box = _load_cgl_geometry(up_file)
    z = cgl_pos[:, 2]
    median_z = float(np.median(z))
    print("Initial CGL placement")
    print(f"  n={len(cgl_pos)} box=({box[0]:.3f}, {box[1]:.3f}, {box[2]:.3f})")
    print(
        "  z min/mean/max/std="
        f"{z.min():.3f}/{z.mean():.3f}/{z.max():.3f}/{z.std():.3f}"
    )
    for name, mask in (("lower", z <= median_z), ("upper", z > median_z)):
        leaflet = cgl_pos[mask]
        leaflet_dirs = dirs[mask]
        nn = _nearest_xy(leaflet[:, :2], box[:2])
        print(
            f"  {name}: n={leaflet.shape[0]} "
            f"z_mean={leaflet[:, 2].mean():.3f} "
            f"z_min={leaflet[:, 2].min():.3f} "
            f"z_max={leaflet[:, 2].max():.3f} "
            f"dir_z_mean={leaflet_dirs[:, 2].mean():.3f} "
            f"nn_min/p05/mean={nn.min():.3f}/{np.percentile(nn, 5):.3f}/{nn.mean():.3f}"
        )


def report_table() -> None:
    with h5py.File(MARTINI_H5, "r") as h5:
        param = h5["cg_lipid_table/cg_lipid_pair/interaction_param"][:].reshape(-1)
        attrs = dict(h5["cg_lipid_table/cg_lipid_pair"].attrs)
        n_radial = int(attrs.get("n_radial", 12))
        n_modes = int(attrs.get("n_modes", 1))
        n_angular = int(attrs.get("n_angular", 15))
    print("CG lipid table")
    if attrs.get("schema") == "cg_lipid_quadspline_v3":
        v0 = param[:n_radial]
        mode_stride = 2 * n_angular + n_radial
        mode_radial = np.array([
            param[n_radial + m * mode_stride + 2 * n_angular:
                  n_radial + m * mode_stride + 2 * n_angular + n_radial]
            for m in range(n_modes)
        ])
        print(f"  v0 maxabs={float(np.max(np.abs(v0))):.6g}")
        print(f"  mode radial maxabs={float(np.max(np.abs(mode_radial))):.6g}")
    else:
        radial = param[30:42]
        angular = param[42:54]
        print(f"  radial maxabs={float(np.max(np.abs(radial))):.6g}")
        print(f"  angular maxabs={float(np.max(np.abs(angular))):.6g}")
    for key in ("schema", "radial_mode", "angle_convention", "n_modes", "n_radial", "cutoff_ang", "knot_spacing_ang"):
        if key in attrs:
            print(f"  {key}={attrs[key]}")


def run_short_dynamics(steps: int, dt: float, frame_steps: int) -> bool:
    shutil.copy2(UP_FILE, SIM_FILE)
    frame_interval = dt * frame_steps
    print(f"Running {steps} steps, dt={dt}, frame_interval={frame_interval}")
    result = subprocess.run(
        [
            str(UPSIDE_BIN),
            str(SIM_FILE),
            "--duration-steps",
            str(steps),
            "--frame-interval",
            str(frame_interval),
            "--temperature",
            "0.8647",
            "--time-step",
            str(dt),
            "--thermostat-timescale",
            "5.0",
            "--thermostat-interval",
            "-1",
            "--seed",
            "42",
            "--integrator",
            "v",
            "--disable-recentering",
            "--record-momentum",
        ],
        cwd=str(RUN_DIR),
        capture_output=True,
        text=True,
        timeout=600,
    )
    if result.stdout:
        print(result.stdout[-4000:])
    if result.returncode != 0:
        print(result.stderr[-4000:])
        print(f"upside exited with {result.returncode}")
        return False
    return True


def report_trajectory() -> None:
    with h5py.File(SIM_FILE, "r") as h5:
        out = _normalize_pos(h5["output/pos"][:]).astype(np.float64)
        cv = h5["input/potential/compose_vector6d"]
        cgl_index = cv["elem_index"][:].astype(np.int64)
        stored_dirs = cv["direction"][:].astype(np.float64)
        orientation_index = (
            cv["orientation_index"][:].astype(np.int64)
            if "orientation_index" in cv
            else None
        )
        pot = h5["input/potential/martini_potential"]
        box = (
            float(pot.attrs["x_len"]),
            float(pot.attrs["y_len"]),
            float(pot.attrs["z_len"]),
        )
        cgl = out[:, cgl_index, :]
        dir0 = _directions_from_sites(out[0], cgl_index, orientation_index, stored_dirs, box)
        start = cgl[0]
        times = h5["output/time"][:]

    print("CGL trajectory")
    for frame in sorted({0, 1, min(5, len(cgl) - 1), len(cgl) - 1}):
        disp = np.sqrt(np.sum((cgl[frame] - start) ** 2, axis=1))
        z = cgl[frame, :, 2]
        dirs = _directions_from_sites(out[frame], cgl_index, orientation_index, stored_dirs, box)
        dir_norm = np.linalg.norm(dirs, axis=1)
        dir_dot = np.sum(dirs * dir0, axis=1)
        dir_angle = np.degrees(np.arccos(np.clip(dir_dot, -1.0, 1.0)))
        print(
            f"  frame={frame} time={float(times[frame]):.6g} "
            f"finite={bool(np.isfinite(cgl[frame]).all())} "
            f"z_min/z_max/z_std={z.min():.3f}/{z.max():.3f}/{z.std():.3f} "
            f"max_disp={disp.max():.3f} mean_disp={disp.mean():.3f} "
            f"dir_norm_min/max={dir_norm.min():.6f}/{dir_norm.max():.6f} "
            f"dir_angle_max/mean={dir_angle.max():.3f}/{dir_angle.mean():.3f}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--skip-stage", action="store_true")
    parser.add_argument("--skip-run", action="store_true")
    parser.add_argument("--resolution", choices=("coarse", "full"), default="coarse")
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--dt", type=float, default=0.002)
    parser.add_argument("--frame-steps", type=int, default=10)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    os.environ.setdefault("UPSIDE_HOME", str(REPO_ROOT))
    os.environ["PYTHONPATH"] = str(REPO_ROOT / "py") + os.pathsep + os.environ.get("PYTHONPATH", "")

    if not args.skip_build:
        write_dopc_only_pdb()
        build_tables(args.resolution)
    elif not RUNTIME_PDB.exists():
        write_dopc_only_pdb()

    if not args.skip_stage:
        convert_and_inject()

    report_initial_geometry(UP_FILE)
    report_table()

    if not args.skip_run:
        ok = run_short_dynamics(args.steps, args.dt, args.frame_steps)
        if ok:
            report_trajectory()
        return 0 if ok else 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
