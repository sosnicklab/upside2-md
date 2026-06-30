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
        default=8,
        help="Number of representative isolated DOPC conformers used for the default CGL projection.",
    )
    parser.add_argument(
        "--isolated-dopc-conformer-pool-count",
        type=int,
        default=32,
        help="Number of isolated bonded DOPC conformers sampled before selecting representative CGL reference conformers.",
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
        help="Metropolis burn-in steps before recording isolated bonded DOPC conformers.",
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
    parser.add_argument(
        "--dopc-target-overlay-base-h5",
        default=None,
        help=(
            "Optional base dopc.h5 to preserve pair/SC/compaction tables from while "
            "rebuilding only effective_lj and cg_lipid_target from the representative "
            "isolated-DOPC ensemble."
        ),
    )
    parser.add_argument(
        "--dopc-target-overlay-output",
        default=None,
        help="Optional output path for the target-overlay dopc.h5.",
    )
    parser.add_argument(
        "--dopc-target-overlay-rebuild-sc",
        action="store_true",
        help=(
            "Also rebuild cg_lipid_sc from the representative isolated-DOPC ensemble "
            "instead of preserving the base SC table."
        ),
    )
    parser.add_argument(
        "--dopc-target-overlay-preserve-target",
        action="store_true",
        help=(
            "Preserve the base effective_lj and cg_lipid_target tables instead of "
            "rebuilding them in the overlay."
        ),
    )
    parser.add_argument(
        "--dopc-target-overlay-reference-up-file",
        action="append",
        default=None,
        help=(
            "Optional full-resolution Upside trajectory to use as the interface-reference "
            "DOPC conformer source for the overlay. Provide multiple times to pool "
            "bulk-bilayer conformers across trajectories."
        ),
    )
    parser.add_argument(
        "--dopc-target-overlay-reference-conformer-count",
        type=int,
        default=12,
        help=(
            "Number of conformers to retain from the pooled full-resolution overlay "
            "reference trajectories when --dopc-target-overlay-reference-up-file is used."
        ),
    )
    parser.add_argument(
        "--dopc-target-overlay-reference-frame-start-fraction",
        type=float,
        default=0.5,
        help=(
            "Fraction of the full-resolution reference trajectory to discard before "
            "sampling overlay DOPC conformers."
        ),
    )
    parser.add_argument(
        "--dopc-target-overlay-reference-min-protein-xy-distance-ang",
        type=float,
        default=0.0,
        help=(
            "If positive, exclude full-resolution overlay reference lipids whose center "
            "lies within this XY distance of any non-lipid particle."
        ),
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
        build_dopc_target_overlay_h5,
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
            "isolated_conformer_pool_count": max(
                int(args.isolated_dopc_conformer_count),
                int(args.isolated_dopc_conformer_pool_count),
            ),
            "isolated_conformer_seed": int(args.isolated_dopc_conformer_seed),
            "isolated_conformer_burnin_steps": max(0, int(args.isolated_dopc_conformer_burnin_steps)),
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

    if bool(args.dopc_target_overlay_base_h5) != bool(args.dopc_target_overlay_output):
        print(
            "ERROR: --dopc-target-overlay-base-h5 and --dopc-target-overlay-output "
            "must be provided together",
            file=sys.stderr,
        )
        return 1
    if args.dopc_target_overlay_base_h5 and args.dopc_target_overlay_output:
        overlay_output = Path(args.dopc_target_overlay_output).expanduser().resolve()
        if overlay_output.exists() and not args.force:
            print(f"Skipping {overlay_output.name} (already exists, use --force to overwrite)")
        else:
            print(f"Generating {overlay_output.name} ...")
            build_dopc_target_overlay_h5(
                base_h5_path=Path(args.dopc_target_overlay_base_h5).expanduser().resolve(),
                output_path=overlay_output,
                dry_ff_path=dry_ff_path,
                lipids_itp_path=lipids_itp_path,
                martinize_path=martinize_path,
                sidechain_lib_path=sidechain_lib_path,
                rebuild_sc=bool(args.dopc_target_overlay_rebuild_sc),
                rebuild_target=not bool(args.dopc_target_overlay_preserve_target),
                conformer_upside_h5_paths=(
                    [Path(path).expanduser().resolve() for path in args.dopc_target_overlay_reference_up_file]
                    if args.dopc_target_overlay_reference_up_file
                    else None
                ),
                conformer_count=max(1, int(args.dopc_target_overlay_reference_conformer_count)),
                conformer_frame_start_fraction=float(
                    args.dopc_target_overlay_reference_frame_start_fraction
                ),
                conformer_min_nonlipid_xy_distance_ang=float(
                    args.dopc_target_overlay_reference_min_protein_xy_distance_ang
                ),
                isolated_conformer_count=max(1, int(args.isolated_dopc_conformer_count)),
                isolated_conformer_pool_count=max(
                    int(args.isolated_dopc_conformer_count),
                    int(args.isolated_dopc_conformer_pool_count),
                ),
                isolated_conformer_seed=int(args.isolated_dopc_conformer_seed),
                isolated_conformer_burnin_steps=max(0, int(args.isolated_dopc_conformer_burnin_steps)),
                isolated_conformer_mc_steps=max(1, int(args.isolated_dopc_conformer_mc_steps)),
                isolated_conformer_proposal_sigma_nm=float(args.isolated_dopc_conformer_sigma_nm),
            )
            print()

    print("All MARTINI parameter files generated successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
