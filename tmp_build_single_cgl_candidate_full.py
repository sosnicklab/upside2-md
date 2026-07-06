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
    args = parser.parse_args()

    base = Path(args.base).expanduser().resolve()
    out = Path(args.out).expanduser().resolve()
    shutil.copy2(base, out)

    corrections = martini_build_tables._build_single_cgl_compaction_corrections_from_base_h5(
        base_h5_path=out,
        include_sc=True,
        include_target=True,
        include_pair_compaction=True,
    )
    with h5py.File(out, "r+") as h5:
        martini_build_tables._apply_single_cgl_compaction_corrections_to_h5(
            h5["cg_lipid_table"],
            corrections,
        )
    print(out)


if __name__ == "__main__":
    main()
