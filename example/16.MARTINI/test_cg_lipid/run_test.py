#!/usr/bin/env python3
"""Minimal CG DOPC lipid test — builds tables, creates .up file, runs short minimization."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "py"))

TEST_DIR = Path(__file__).resolve().parent
RUN_DIR = TEST_DIR / "run"
PDB_ID = "test_cg_lipid"
RUNTIME_PDB = TEST_DIR / "test_dopc_gly.pdb"
MARTINI_H5 = RUN_DIR / "martini.h5"
UP_FILE = RUN_DIR / "test.input.up"
UPSIDE_BIN = REPO_ROOT / "obj" / "upside"


def step1_build_tables():
    """Build martini.h5 with CG lipid quadspline + particles table extension."""
    print("=== Step 1: Build martini tables ===")
    RUN_DIR.mkdir(parents=True, exist_ok=True)

    from martini_build_tables import build_martini_tables, _compute_lipid_bonded_energy

    # Parse DOPC reference bead positions from the test PDB
    atoms = []
    with open(RUNTIME_PDB) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:21].strip().upper()
            if resname not in ("DOP", "DOPC"):
                continue
            aname = line[12:16].strip().upper()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append({"name": aname, "x": x, "y": y, "z": z})

    n_per_lipid = 14
    n_lipids = len(atoms) // n_per_lipid
    assert n_lipids > 0, "No DOPC lipids found in PDB"
    print(f"  Found {n_lipids} DOPC lipids ({len(atoms)} atoms)")

    # Use first lipid as reference
    first = atoms[:n_per_lipid]
    com = np.mean([[a["x"], a["y"], a["z"]] for a in first], axis=0)
    ref_bead_positions = np.array(
        [[a["x"] - com[0], a["y"] - com[1], a["z"] - com[2]] for a in first]
    )
    ref_bead_positions_nm = ref_bead_positions * 0.1  # Å → nm

    # Verify bond energy of reference positions is reasonable
    ref_bonded = _compute_lipid_bonded_energy(ref_bead_positions_nm)
    print(f"  Reference bonded energy: {ref_bonded:.2f} kJ/mol (should be small)")
    if ref_bonded > 500.0:
        print(f"  WARNING: Reference positions have high bonded energy; PDB may need relaxation")

    bead_types = [
        "Q0", "Qa", "Na", "Na",
        "C1", "C1", "C3", "C1", "C1",
        "C1", "C1", "C3", "C1", "C1",
    ]

    ff_dir = Path(os.environ.get("UPSIDE_MARTINI_FF_DIR", REPO_ROOT / "parameters" / "dryMARTINI"))
    dry_ff = ff_dir / "dry_martini_v2.1.itp"
    sc_lib = REPO_ROOT / "parameters" / "ff_2.1" / "sidechain.h5"
    martinize = REPO_ROOT / "py" / "martinize.py"

    cg_lipid_config = {
        "ref_bead_positions_nm": ref_bead_positions_nm,
        "bead_types": bead_types,
    }

    build_martini_tables(
        output_path=MARTINI_H5,
        dry_ff_path=dry_ff,
        martinize_path=martinize,
        sidechain_lib_path=sc_lib,
        forcefield_name="martini22",
        cg_lipid_config=cg_lipid_config,
    )
    print(f"  Built: {MARTINI_H5}")


def step2_convert_stage():
    """Convert PDB to .up file via convert_stage()."""
    print("=== Step 2: Convert stage ===")
    os.environ["UPSIDE_RUNTIME_PDB_FILE"] = str(RUNTIME_PDB.resolve())
    os.environ["UPSIDE_SIMULATION_STAGE"] = "minimization"
    os.environ["UPSIDE_MARTINI_ENERGY_CONVERSION"] = "2.914952774272"
    os.environ["UPSIDE_MARTINI_LENGTH_CONVERSION"] = "10.0"

    from martini_prepare_system_lib import convert_stage

    convert_stage(
        pdb_id=PDB_ID,
        stage="minimization",
        run_dir=str(RUN_DIR),
    )
    print(f"  Created: {UP_FILE}")


def step3_inject_tables():
    """Inject particles table and CG lipid nodes into .up file."""
    print("=== Step 3: Inject tables ===")
    from martini_prepare_system_lib import inject_particles_table, inject_cg_lipid_nodes

    inject_particles_table(UP_FILE, MARTINI_H5)
    print("  Injected particles table (with CGL type)")

    inject_cg_lipid_nodes(UP_FILE, MARTINI_H5)
    print("  Injected CG lipid nodes (quadspline + compose_vector6d)")


def step4_run_minimization():
    """Run a short minimization with the upside engine."""
    print("=== Step 4: Run minimization ===")
    checkpoints = RUN_DIR / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        [
            str(UPSIDE_BIN),
            str(UP_FILE),
            "--minimize",
            "--min-max-iter", "100",
            "--temperature", "0.8647",
            "--time-step", "0.01",
            "--seed", "42",
            "--frame-interval", "1.0",
            "--duration-steps", "0",
            "--disable-recentering",
            "-o", str(RUN_DIR / "minimize"),
        ],
        cwd=str(RUN_DIR),
        capture_output=True,
        text=True,
        timeout=300,
    )
    print(result.stdout[-5000:] if len(result.stdout) > 5000 else result.stdout)
    if result.returncode != 0:
        print(f"STDERR (last 3000 chars):\n{result.stderr[-3000:]}")
        print(f"Exit code: {result.returncode}")
        return False
    print("  Minimization completed successfully!")
    return True


def main():
    os.environ.setdefault("UPSIDE_MARTINI_FF_DIR", str(REPO_ROOT / "parameters" / "dryMARTINI"))
    os.environ["PYTHONPATH"] = str(REPO_ROOT / "py") + ":" + os.environ.get("PYTHONPATH", "")

    import argparse
    parser = argparse.ArgumentParser(description="CG lipid minimal test")
    parser.add_argument("--skip-build", action="store_true", help="Skip table building")
    parser.add_argument("--skip-stage", action="store_true", help="Skip stage conversion")
    parser.add_argument("--skip-inject", action="store_true", help="Skip table injection")
    parser.add_argument("--min-only", action="store_true", help="Only run minimization")
    args = parser.parse_args()

    if not args.min_only:
        if not args.skip_build:
            step1_build_tables()
        if not args.skip_stage:
            step2_convert_stage()
        if not args.skip_inject:
            step3_inject_tables()

    step4_run_minimization()


if __name__ == "__main__":
    main()
