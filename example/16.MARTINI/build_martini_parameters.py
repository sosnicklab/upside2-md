#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path

import prepare_system as prep


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Build refined dry-MARTINI/Upside parameter tables and publish the "
            "runtime cross artifact to parameters/ff_2.1/martini.h5."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR / "outputs" / "parameter_build",
        help="Directory for intermediate CSV/JSON/HDF5 build outputs.",
    )
    parser.add_argument(
        "--publish-output",
        type=Path,
        default=REPO_ROOT / "parameters" / "ff_2.1" / "martini.h5",
        help="Published runtime artifact path.",
    )
    parser.add_argument(
        "--bilayer-pdb",
        type=Path,
        default=SCRIPT_DIR / "pdb" / "bilayer.MARTINI.pdb",
    )
    parser.add_argument(
        "--lipid-itp",
        type=Path,
        default=SCRIPT_DIR / "ff_dry" / "dry_martini_v2.1_lipids.itp",
    )
    parser.add_argument(
        "--ff-itp",
        type=Path,
        default=SCRIPT_DIR / "ff_dry" / "dry_martini_v2.1.itp",
    )
    parser.add_argument(
        "--membrane-h5",
        type=Path,
        default=REPO_ROOT / "parameters" / "ff_2.1" / "membrane.h5",
    )
    parser.add_argument("--radial-min", type=float, default=2.0)
    parser.add_argument("--radial-max", type=float, default=12.0)
    parser.add_argument("--sample-step", type=float, default=0.1)
    parser.add_argument("--n-radial", type=int, default=12)
    parser.add_argument("--ridge", type=float, default=1e-4)
    parser.add_argument("--validate", type=int, choices=[0, 1], default=1)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    output_dir = args.output_dir.expanduser().resolve()
    publish_output = args.publish_output.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    publish_output.parent.mkdir(parents=True, exist_ok=True)

    depth_csv = output_dir / "depth_interaction_table.csv"
    depth_meta = output_dir / "depth_interaction_table.meta.json"
    cross_csv = output_dir / "backbone_cross_interaction_table.csv"
    cross_meta = output_dir / "backbone_cross_interaction_table.meta.json"
    cross_h5 = output_dir / "backbone_cross_interaction_table.h5"

    prep.run_build_depth_table_command(
        [
            "--bilayer-pdb",
            str(args.bilayer_pdb.expanduser().resolve()),
            "--lipid-itp",
            str(args.lipid_itp.expanduser().resolve()),
            "--membrane-h5",
            str(args.membrane_h5.expanduser().resolve()),
            "--output-csv",
            str(depth_csv),
            "--output-json",
            str(depth_meta),
        ]
    )

    prep.run_build_backbone_cross_table_command(
        [
            "--ff-itp",
            str(args.ff_itp.expanduser().resolve()),
            "--depth-table-csv",
            str(depth_csv),
            "--depth-table-meta",
            str(depth_meta),
            "--output-csv",
            str(cross_csv),
            "--output-json",
            str(cross_meta),
            "--output-h5",
            str(cross_h5),
            "--radial-min",
            str(args.radial_min),
            "--radial-max",
            str(args.radial_max),
            "--sample-step",
            str(args.sample_step),
            "--n-radial",
            str(args.n_radial),
            "--ridge",
            str(args.ridge),
        ]
    )

    if args.validate:
        prep.run_validate_backbone_only_command(
            [
                "--table-csv",
                str(depth_csv),
                "--table-meta",
                str(depth_meta),
                "--cross-table-csv",
                str(cross_csv),
                "--cross-table-meta",
                str(cross_meta),
                "--cross-artifact",
                str(cross_h5),
            ]
        )

    shutil.copy2(cross_h5, publish_output)

    if args.validate:
        prep.run_validate_backbone_only_command(
            [
                "--cross-artifact",
                str(publish_output),
            ]
        )

    print(f"Published martini cross artifact to: {publish_output}")
    print(f"Build outputs written to: {output_dir}")


if __name__ == "__main__":
    main()
