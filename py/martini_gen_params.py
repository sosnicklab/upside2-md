#!/usr/bin/env python3
"""Build dry-MARTINI parameter HDF5 files under ``parameters/dryMARTINI``."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
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
        help="Optional bead-frame samples around each sampled direction vector.",
    )
    parser.add_argument(
        "--sc-bead-frame-count",
        type=int,
        default=None,
        help="Optional SC bead-frame samples for SC-particle and SC-CGL tables.",
    )
    parser.add_argument(
        "--cgl-bead-frame-count",
        type=int,
        default=None,
        help="Optional CGL bead-frame samples for CGL-particle, SC-CGL, and CGL-CGL tables.",
    )
    parser.add_argument(
        "--isolated-dopc-conformer-count",
        type=int,
        default=8,
        help="Representative isolated DOPC conformers used for the default CGL projection.",
    )
    parser.add_argument(
        "--isolated-dopc-conformer-pool-count",
        type=int,
        default=32,
        help="Isolated DOPC conformers sampled before representative selection.",
    )
    parser.add_argument(
        "--isolated-dopc-conformer-seed",
        type=int,
        default=1729,
        help="Random seed for isolated bonded DOPC conformer sampling.",
    )
    parser.add_argument(
        "--isolated-dopc-conformer-burnin-steps",
        type=int,
        default=5000,
        help="Metropolis burn-in steps before recording isolated DOPC conformers.",
    )
    parser.add_argument(
        "--isolated-dopc-conformer-mc-steps",
        type=int,
        default=2000,
        help="Metropolis steps between isolated DOPC conformer samples.",
    )
    parser.add_argument(
        "--isolated-dopc-conformer-sigma-nm",
        type=float,
        default=0.025,
        help="Cartesian proposal sigma in nm for isolated DOPC conformer sampling.",
    )
    return parser.parse_args(argv)


def configure_environment(args: argparse.Namespace) -> None:
    env_updates = {
        "UPSIDE_MARTINI_TABLE_WORKERS": args.workers,
        "UPSIDE_MARTINI_BEAD_FRAME_COUNT": args.bead_frame_count,
        "UPSIDE_MARTINI_SC_BEAD_FRAME_COUNT": args.sc_bead_frame_count,
        "UPSIDE_MARTINI_CGL_BEAD_FRAME_COUNT": args.cgl_bead_frame_count,
    }
    for name, value in env_updates.items():
        if value is not None:
            os.environ[name] = str(max(1, int(value)))


def resolve_repo_root(args: argparse.Namespace) -> Path:
    repo_root = (
        Path(args.upside_home).expanduser().resolve()
        if args.upside_home
        else Path(__file__).resolve().parent.parent
    )
    if not repo_root.exists():
        raise FileNotFoundError(f"UPSIDE_HOME does not exist: {repo_root}")
    return repo_root


def required_inputs(repo_root: Path) -> dict[str, Path]:
    paths = {
        "dry_ff_path": repo_root / "parameters" / "dryMARTINI" / "dry_martini_v2.1.itp",
        "lipids_itp_path": repo_root / "parameters" / "dryMARTINI" / "dry_martini_v2.1_lipids.itp",
        "dopc_pdb_path": repo_root / "parameters" / "dryMARTINI" / "DOPC.pdb",
        "martinize_path": repo_root / "py" / "martinize.py",
        "sidechain_lib_path": repo_root / "parameters" / "ff_2.1" / "sidechain.h5",
    }
    for path in paths.values():
        if not path.exists():
            raise FileNotFoundError(path)
    return paths


def build_targets(repo_root: Path, args: argparse.Namespace) -> list[tuple[Path, object, dict]]:
    output_dir = repo_root / "parameters" / "dryMARTINI"
    output_dir.mkdir(parents=True, exist_ok=True)

    inputs = required_inputs(repo_root)

    from martini_build_tables import build_dopc_h5, build_particle_h5, build_sidechain_h5

    return [
        (
            output_dir / "particle.h5",
            build_particle_h5,
            {
                "output_path": output_dir / "particle.h5",
                "dry_ff_path": inputs["dry_ff_path"],
            },
        ),
        (
            output_dir / "sidechain.h5",
            build_sidechain_h5,
            {
                "output_path": output_dir / "sidechain.h5",
                "dry_ff_path": inputs["dry_ff_path"],
                "martinize_path": inputs["martinize_path"],
                "sidechain_lib_path": inputs["sidechain_lib_path"],
            },
        ),
        (
            output_dir / "dopc.h5",
            build_dopc_h5,
            {
                "output_path": output_dir / "dopc.h5",
                "dry_ff_path": inputs["dry_ff_path"],
                "lipids_itp_path": inputs["lipids_itp_path"],
                "martinize_path": inputs["martinize_path"],
                "sidechain_lib_path": inputs["sidechain_lib_path"],
                "dopc_pdb_path": inputs["dopc_pdb_path"],
                "isolated_conformer_count": max(1, int(args.isolated_dopc_conformer_count)),
                "isolated_conformer_pool_count": max(
                    int(args.isolated_dopc_conformer_count),
                    int(args.isolated_dopc_conformer_pool_count),
                ),
                "isolated_conformer_seed": int(args.isolated_dopc_conformer_seed),
                "isolated_conformer_burnin_steps": max(
                    0, int(args.isolated_dopc_conformer_burnin_steps)
                ),
                "isolated_conformer_mc_steps": max(
                    1, int(args.isolated_dopc_conformer_mc_steps)
                ),
                "isolated_conformer_proposal_sigma_nm": float(
                    args.isolated_dopc_conformer_sigma_nm
                ),
            },
        ),
    ]


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    configure_environment(args)

    try:
        repo_root = resolve_repo_root(args)
        targets = build_targets(repo_root, args)
    except FileNotFoundError as exc:
        print(f"ERROR: required file not found: {exc}", file=sys.stderr)
        return 1

    output_dir = repo_root / "parameters" / "dryMARTINI"
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

    for output_path, builder, kwargs in targets:
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
