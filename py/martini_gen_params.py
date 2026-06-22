#!/usr/bin/env python3
"""Build dry-MARTINI parameter HDF5 files.

The generated files live under ``parameters/dryMARTINI``:

    particle.h5
    sidechain.h5
    dopc.h5

Typical use:

    python py/martini_gen_params.py --upside-home /path/to/repo

Use ``--force`` to replace existing files.
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
    parser.add_argument(
        "--isolated-dopc-conformer-count",
        type=int,
        default=2,
        help="Number of isolated bonded DOPC conformers for the default CGL projection.",
    )
    parser.add_argument(
        "--isolated-dopc-conformer-seed",
        type=int,
        default=1729,
        help="Random seed for isolated bonded DOPC conformer sampling.",
    )
    parser.add_argument(
        "--isolated-dopc-conformer-mc-steps",
        type=int,
        default=2000,
        help="Metropolis steps between isolated bonded DOPC conformer samples.",
    )
    parser.add_argument(
        "--isolated-dopc-conformer-sigma-nm",
        type=float,
        default=0.025,
        help="Cartesian proposal sigma in nm for isolated bonded DOPC conformer sampling.",
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

    from martini_build_tables import (
        build_particle_h5,
        build_sidechain_h5,
        build_dopc_h5,
    )

    builders = [
        (output_dir / "particle.h5", build_particle_h5, {
            "output_path": output_dir / "particle.h5",
            "dry_ff_path": dry_ff_path,
        }),
        (output_dir / "sidechain.h5", build_sidechain_h5, {
            "output_path": output_dir / "sidechain.h5",
            "dry_ff_path": dry_ff_path,
            "martinize_path": martinize_path,
            "sidechain_lib_path": sidechain_lib_path,
        }),
        (output_dir / "dopc.h5", build_dopc_h5, {
            "output_path": output_dir / "dopc.h5",
            "dry_ff_path": dry_ff_path,
            "lipids_itp_path": lipids_itp_path,
            "martinize_path": martinize_path,
            "sidechain_lib_path": sidechain_lib_path,
            "dopc_pdb_path": dopc_pdb_path,
            "isolated_conformer_count": max(1, int(args.isolated_dopc_conformer_count)),
            "isolated_conformer_seed": int(args.isolated_dopc_conformer_seed),
            "isolated_conformer_mc_steps": max(1, int(args.isolated_dopc_conformer_mc_steps)),
            "isolated_conformer_proposal_sigma_nm": float(args.isolated_dopc_conformer_sigma_nm),
        }),
    ]

    print(f"Upside home: {repo_root}")
    print(f"Output directory: {output_dir}")
    print(
        "Parallel workers: "
        f"{os.environ.get('UPSIDE_MARTINI_TABLE_WORKERS', 'auto')}; "
        "optional bead-frame samples: "
        f"SC={os.environ.get('UPSIDE_MARTINI_SC_BEAD_FRAME_COUNT', os.environ.get('UPSIDE_MARTINI_BEAD_FRAME_COUNT', '1'))}, "
        f"CGL={os.environ.get('UPSIDE_MARTINI_CGL_BEAD_FRAME_COUNT', os.environ.get('UPSIDE_MARTINI_BEAD_FRAME_COUNT', '4'))}"
    )
    print()

    for output_path, builder, kwargs in builders:
        if output_path.exists() and not args.force:
            print(f"Skipping {output_path.name} (already exists, use --force to overwrite)")
            continue
        print(f"Generating {output_path.name} ...")
        builder(**kwargs)
        print()

    print("All MARTINI parameter files generated successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
