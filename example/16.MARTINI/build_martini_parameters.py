#!/usr/bin/env python3

import argparse
import importlib.util
import json
import os
import pickle as cp
import shutil
import sys
from pathlib import Path

import prepare_system as prep


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]


def _guess_condiv_worker_script(checkpoint: Path) -> Path:
    candidates = [
        checkpoint.parent / "ConDiv_mem.py",
        checkpoint.parent.parent / "ConDiv_mem.py",
        checkpoint.parent.parent.parent / "ConDiv_mem.py",
        REPO_ROOT / "ConDiv_symlay" / "ConDiv_mem.py",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise RuntimeError(f"Could not locate ConDiv_mem.py for checkpoint {checkpoint}")


def _load_condiv_module(worker_script: Path):
    os.environ.setdefault("CONDIV_PROJECT_ROOT", str(REPO_ROOT))
    script_dir = str(worker_script.parent)
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    spec = importlib.util.spec_from_file_location("ConDiv_mem", str(worker_script))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load ConDiv_mem module from {worker_script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["ConDiv_mem"] = module
    spec.loader.exec_module(module)
    return module


def materialize_membrane_from_checkpoint(checkpoint: Path, output_path: Path) -> dict:
    checkpoint = checkpoint.expanduser().resolve()
    worker_script = _guess_condiv_worker_script(checkpoint)
    module = _load_condiv_module(worker_script)
    with checkpoint.open("rb") as fh:
        state = cp.load(fh)

    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    module.expand_param(state["param"], state["init_param_files"], {"memb": str(output_path)})
    return {
        "condiv_checkpoint": str(checkpoint),
        "worker_script": str(worker_script),
        "materialized_membrane_h5": str(output_path),
        "layer_manifest_path": state.get("layer_manifest_path", ""),
        "pair_ff_itp": state.get("pair_ff_itp", ""),
        "seed_source_init_dir": state.get("seed_source_init_dir", ""),
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Build the published MARTINI force-field artifact with both the "
            "dry-dry nonbond table and the refined runtime backbone cross table."
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
    parser.add_argument(
        "--condiv-checkpoint",
        type=Path,
        default=None,
        help="Optional ConDiv_symlay checkpoint to materialize and use as the membrane teacher source.",
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
    trained_source = None
    membrane_source = args.membrane_h5.expanduser().resolve()
    if args.condiv_checkpoint is not None:
        materialized_membrane = output_dir / "trained_membrane" / "membrane.h5"
        trained_source = materialize_membrane_from_checkpoint(args.condiv_checkpoint, materialized_membrane)
        membrane_source = materialized_membrane

    prep.run_build_depth_table_command(
        [
            "--bilayer-pdb",
            str(args.bilayer_pdb.expanduser().resolve()),
            "--lipid-itp",
            str(args.lipid_itp.expanduser().resolve()),
            "--membrane-h5",
            str(membrane_source),
            "--output-csv",
            str(depth_csv),
            "--output-json",
            str(depth_meta),
        ]
    )
    if trained_source is not None:
        depth_meta_payload = json.loads(depth_meta.read_text(encoding="utf-8"))
        depth_meta_payload["trained_source"] = trained_source
        depth_meta.write_text(json.dumps(depth_meta_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

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

    prep.write_dry_nonbond_forcefield_artifact(
        cross_h5,
        args.ff_itp.expanduser().resolve(),
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

    print(f"Published MARTINI force-field artifact to: {publish_output}")
    print(f"Build outputs written to: {output_dir}")
    if trained_source is not None:
        print(f"Trained ConDiv checkpoint source: {trained_source['condiv_checkpoint']}")
        print(f"Materialized membrane teacher: {trained_source['materialized_membrane_h5']}")


if __name__ == "__main__":
    main()
