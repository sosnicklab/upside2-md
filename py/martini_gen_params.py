#!/usr/bin/env python3
"""One-shot pre-generation of MARTINI parameter .h5 files.

Generates the four parameter files under ``parameters/dryMARTINI/``:

    particle.h5      - all 38 particle-type pair LJ+Coulomb energy grids
    sidechain.h5     - all 20 residues x 38 target SC orientation tables
    dopc.h5          - DOPC CG lipid directional splines
    interlipid.h5    - cross-lipid placeholder (empty)

Run this once after cloning or when the .itp files change:

    python py/martini_gen_params.py --upside-home /path/to/repo

Pass ``--force`` to regenerate files that already exist.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Pre-generate MARTINI parameter .h5 files under parameters/dryMARTINI/"
    )
    parser.add_argument(
        "--upside-home",
        default=None,
        help="Path to the Upside repository root (default: auto-detect from this script)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing .h5 files instead of skipping them",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help=(
            "Parallel table-build workers. Defaults to UPSIDE_MARTINI_TABLE_WORKERS, "
            "then Slurm CPU allocation, then local CPU count."
        ),
    )
    parser.add_argument(
        "--bead-frame-count",
        type=int,
        default=None,
        help="Set optional bead-frame samples around each sampled direction vector.",
    )
    parser.add_argument(
        "--sc-bead-frame-count",
        type=int,
        default=None,
        help="Set optional SC bead-frame samples for SC-particle and SC-CGL tables.",
    )
    parser.add_argument(
        "--cgl-bead-frame-count",
        type=int,
        default=None,
        help="Set optional CGL bead-frame samples for CGL-particle, SC-CGL, and CGL-CGL tables.",
    )
    args = parser.parse_args(argv)

    if args.workers is not None:
        os.environ["UPSIDE_MARTINI_TABLE_WORKERS"] = str(max(1, int(args.workers)))
    if args.bead_frame_count is not None:
        os.environ["UPSIDE_MARTINI_BEAD_FRAME_COUNT"] = str(max(1, int(args.bead_frame_count)))
    if args.sc_bead_frame_count is not None:
        os.environ["UPSIDE_MARTINI_SC_BEAD_FRAME_COUNT"] = str(max(1, int(args.sc_bead_frame_count)))
    if args.cgl_bead_frame_count is not None:
        os.environ["UPSIDE_MARTINI_CGL_BEAD_FRAME_COUNT"] = str(max(1, int(args.cgl_bead_frame_count)))

    if args.upside_home:
        repo_root = Path(args.upside_home).expanduser().resolve()
    else:
        repo_root = Path(__file__).resolve().parent.parent

    if not repo_root.exists():
        print(f"ERROR: UPSIDE_HOME does not exist: {repo_root}", file=sys.stderr)
        return 1

    # Resolve required paths
    dry_ff_path = repo_root / "parameters" / "dryMARTINI" / "dry_martini_v2.1.itp"
    lipids_itp_path = repo_root / "parameters" / "dryMARTINI" / "dry_martini_v2.1_lipids.itp"
    dopc_pdb_path = repo_root / "parameters" / "dryMARTINI" / "DOPC.pdb"
    martinize_path = repo_root / "py" / "martinize.py"
    sidechain_lib_path = repo_root / "parameters" / "ff_2.1" / "sidechain.h5"
    output_dir = repo_root / "parameters" / "dryMARTINI"

    for path in [dry_ff_path, lipids_itp_path, dopc_pdb_path, martinize_path, sidechain_lib_path]:
        if not path.exists():
            print(f"ERROR: required file not found: {path}", file=sys.stderr)
            return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    files_to_generate = {
        output_dir / "particle.h5": "build_particle_h5",
        output_dir / "sidechain.h5": "build_sidechain_h5",
        output_dir / "dopc.h5": "build_dopc_h5",
        output_dir / "interlipid.h5": "build_interlipid_h5",
    }

    from martini_build_tables import (
        build_particle_h5,
        build_sidechain_h5,
        build_dopc_h5,
        build_interlipid_h5,
    )

    builders: dict = {
        "build_particle_h5": lambda: build_particle_h5(
            output_path=output_dir / "particle.h5",
            dry_ff_path=dry_ff_path,
        ),
        "build_sidechain_h5": lambda: build_sidechain_h5(
            output_path=output_dir / "sidechain.h5",
            dry_ff_path=dry_ff_path,
            martinize_path=martinize_path,
            sidechain_lib_path=sidechain_lib_path,
        ),
        "build_dopc_h5": lambda: build_dopc_h5(
            output_path=output_dir / "dopc.h5",
            dry_ff_path=dry_ff_path,
            lipids_itp_path=lipids_itp_path,
            martinize_path=martinize_path,
            sidechain_lib_path=sidechain_lib_path,
            dopc_pdb_path=dopc_pdb_path,
        ),
        "build_interlipid_h5": lambda: build_interlipid_h5(
            output_path=output_dir / "interlipid.h5",
        ),
    }

    print(f"Upside home: {repo_root}")
    print(f"Output directory: {output_dir}")
    print(
        "Parallel workers: "
        f"{os.environ.get('UPSIDE_MARTINI_TABLE_WORKERS', 'auto')}; "
        "optional bead-frame samples: "
        f"SC={os.environ.get('UPSIDE_MARTINI_SC_BEAD_FRAME_COUNT', os.environ.get('UPSIDE_MARTINI_BEAD_FRAME_COUNT', '1'))}, "
        f"CGL={os.environ.get('UPSIDE_MARTINI_CGL_BEAD_FRAME_COUNT', os.environ.get('UPSIDE_MARTINI_BEAD_FRAME_COUNT', '1'))}"
    )
    print()

    for output_path, builder_name in files_to_generate.items():
        if output_path.exists() and not args.force:
            print(f"Skipping {output_path.name} (already exists, use --force to overwrite)")
            continue
        print(f"Generating {output_path.name} ...")
        builders[builder_name]()
        print()

    print("All MARTINI parameter files generated successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
