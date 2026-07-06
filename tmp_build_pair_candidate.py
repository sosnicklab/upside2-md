#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import h5py

import py.martini_build_tables as martini_build_tables


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--skip-pair-relax", action="store_true")
    parser.add_argument("--preserve-existing-pair-states", action="store_true")
    args = parser.parse_args()

    base = Path(args.base).expanduser().resolve()
    out = Path(args.out).expanduser().resolve()
    shutil.copy2(base, out)
    if args.skip_pair_relax:
        with h5py.File(out, "r+") as h5:
            comp = h5["cg_lipid_table"]["cg_lipid_compaction"]
            source = comp.attrs.get("source", "")
            if isinstance(source, bytes):
                source = source.decode("utf-8", errors="ignore")
            source = str(source)
            comp.attrs["source"] = source.replace("pair_relaxation_", "").replace("pair_relax", "")
    corrections = martini_build_tables._build_single_cgl_compaction_corrections_from_base_h5(
        base_h5_path=out,
        include_sc=False,
        include_target=False,
        include_pair_compaction=True,
    )
    if args.preserve_existing_pair_states:
        with h5py.File(base, "r") as h5:
            comp = h5["cg_lipid_table"]["cg_lipid_compaction"]
            preserved_datasets = (
                "delta_extended_extended",
                "delta_extended_compact",
                "delta_compact_compact",
                "grid_extended_extended_kj_mol",
                "grid_extended_compact_kj_mol",
                "grid_compact_compact_kj_mol",
                "grid_average_kj_mol",
                "face_mask",
            )
            for name in preserved_datasets:
                if name in comp:
                    corrections["cg_lipid_compaction"][name] = comp[name][:]
            for attr_name in (
                "pair_reference_extended_center_ang",
                "pair_reference_compact_center_ang",
                "pair_compaction_state_source",
                "correction_center_mode",
                "pair_state_model",
                "face_cos_min",
                "radial_cutoff_nm",
                "mask_mode",
            ):
                if attr_name in comp.attrs:
                    corrections["cg_lipid_compaction"][attr_name] = comp.attrs[attr_name]
    with h5py.File(out, "r+") as h5:
        martini_build_tables._apply_single_cgl_compaction_corrections_to_h5(
            h5["cg_lipid_table"],
            corrections,
        )
    print(out)


if __name__ == "__main__":
    main()
