#!/usr/bin/env python3
"""Build the DOPC topology-slot manifest for the ConDiv_symlay workflow."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from symlay_utils import (
    build_layer_manifest,
    load_layer_template,
    repo_relative,
    write_layer_manifest_csv,
    write_layer_manifest_plot,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the ConDiv_symlay DOPC layer manifest")
    parser.add_argument(
        "--template",
        type=Path,
        default=Path(__file__).resolve().parent / "layer_template.json",
        help="Path to the layer template JSON",
    )
    parser.add_argument(
        "--bilayer-pdb",
        type=Path,
        default=Path("example/16.MARTINI/pdb/bilayer.MARTINI.pdb"),
        help="Bilayer template PDB used to measure bead depths",
    )
    parser.add_argument(
        "--lipid-itp",
        type=Path,
        default=Path("example/16.MARTINI/ff_dry/dry_martini_v2.1_lipids.itp"),
        help="Dry-MARTINI lipid ITP containing the DOPC topology",
    )
    parser.add_argument("--output-json", type=Path, required=True, help="Manifest JSON output path")
    parser.add_argument("--output-csv", type=Path, required=True, help="Manifest CSV output path")
    parser.add_argument("--output-png", type=Path, required=True, help="Manifest PNG output path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    template = load_layer_template(args.template)
    template["_source_path"] = str(args.template.expanduser().resolve())
    manifest = build_layer_manifest(template, args.lipid_itp, args.bilayer_pdb)

    output_json = args.output_json.expanduser().resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_layer_manifest_csv(args.output_csv, manifest)
    write_layer_manifest_plot(args.output_png, manifest)

    print("Built ConDiv_symlay layer manifest")
    print(f"  template: {repo_relative(args.template)}")
    print(f"  lipid itp: {repo_relative(args.lipid_itp)}")
    print(f"  bilayer pdb: {repo_relative(args.bilayer_pdb)}")
    print(f"  output json: {repo_relative(output_json)}")
    print(f"  output csv: {repo_relative(args.output_csv)}")
    print(f"  output png: {repo_relative(args.output_png)}")
    print("  full type sequence: " + "-".join(manifest["full_type_sequence"]))
    print(
        "  positive projection depths (A): "
        + ", ".join(f"{depth:.3f}" for depth in manifest["positive_projection_depths"])
    )


if __name__ == "__main__":
    main()
