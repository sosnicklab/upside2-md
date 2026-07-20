#!/usr/bin/env python
"""Build the consolidated dry-MARTINI force-field file ``parameters/ff_2.1/martini.h5``."""

import argparse
import os
import sys
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pre-generate the consolidated MARTINI force-field file parameters/ff_2.1/martini.h5"
    )
    parser.add_argument(
        "--upside-home",
        default=None,
        help="Path to the Upside repository root (default: auto-detect from this script)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing martini.h5 instead of skipping it",
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
        help="Optional SC bead-frame samples for the SC-particle table.",
    )
    return parser.parse_args(argv)


def configure_environment(args: argparse.Namespace) -> None:
    env_updates = {
        "UPSIDE_MARTINI_TABLE_WORKERS": args.workers,
        "UPSIDE_MARTINI_BEAD_FRAME_COUNT": args.bead_frame_count,
        "UPSIDE_MARTINI_SC_BEAD_FRAME_COUNT": args.sc_bead_frame_count,
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
    itp_dir = repo_root / "example" / "16.MARTINI" / "dryMARTINI_itp"
    paths = {
        "dry_ff_path": itp_dir / "dry_martini_v2.1.itp",
        "martinize_path": repo_root / "py" / "martinize.py",
        "sidechain_lib_path": repo_root / "parameters" / "ff_2.1" / "sidechain.h5",
    }
    for path in paths.values():
        if not path.exists():
            raise FileNotFoundError(path)
    return paths


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    configure_environment(args)

    try:
        repo_root = resolve_repo_root(args)
        inputs = required_inputs(repo_root)
    except FileNotFoundError as exc:
        print(f"ERROR: required file not found: {exc}", file=sys.stderr)
        return 1

    output_dir = repo_root / "parameters" / "ff_2.1"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "martini.h5"

    print(f"Upside home: {repo_root}")
    print(f"Output file: {output_path}")
    print(
        "Parallel workers: "
        f"{os.environ.get('UPSIDE_MARTINI_TABLE_WORKERS', 'auto')}; "
        "optional bead-frame samples: "
        f"SC={os.environ.get('UPSIDE_MARTINI_SC_BEAD_FRAME_COUNT', os.environ.get('UPSIDE_MARTINI_BEAD_FRAME_COUNT', '1'))}"
    )
    print()

    if output_path.exists() and not args.force:
        print(f"Skipping {output_path.name} (already exists, use --force to overwrite)")
        return 0

    from martini_build_tables import build_martini_h5

    print(f"Generating {output_path.name} (/particles + /sc_table) ...")
    build_martini_h5(
        output_path=output_path,
        dry_ff_path=inputs["dry_ff_path"],
        martinize_path=inputs["martinize_path"],
        sidechain_lib_path=inputs["sidechain_lib_path"],
    )
    print()
    print("MARTINI force-field file generated successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
