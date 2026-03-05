#!/usr/bin/env python3

import argparse
from pathlib import Path

import h5py
import numpy as np


def require_group(root, path):
    if path not in root:
        raise ValueError(f"Missing group: {path}")
    return root[path]


def require_dataset(group, name):
    if name not in group:
        raise ValueError(f"Missing dataset: {group.name}/{name}")
    return group[name]


def decode_attr_string(value, default=""):
    if value is None:
        return default
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)


def main():
    parser = argparse.ArgumentParser(
        description="Validate backbone-only mapping HDF5 schema and index consistency."
    )
    parser.add_argument("mapping_h5", type=Path, help="Path to mapping HDF5 file.")
    parser.add_argument("--n-atom", type=int, default=None, help="Optional expected n_atom value.")
    args = parser.parse_args()

    p = args.mapping_h5.expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(p)

    with h5py.File(p, "r") as h5:
        inp = require_group(h5, "/input")
        bb = require_group(inp, "hybrid_bb_map")
        env = require_group(inp, "hybrid_env_topology")

        bb_atom_idx = require_dataset(bb, "bb_atom_index")[:]
        bb_atom_map = require_dataset(bb, "atom_indices")[:]
        bb_mask = require_dataset(bb, "atom_mask")[:]
        bb_w = require_dataset(bb, "weights")[:]

        if bb_atom_map.ndim != 2 or bb_atom_map.shape[1] != 4:
            raise ValueError("hybrid_bb_map/atom_indices must have shape (n_bb,4)")
        if bb_mask.shape != bb_atom_map.shape or bb_w.shape != bb_atom_map.shape:
            raise ValueError("hybrid_bb_map mask/weights shapes must match atom_indices")

        n_bb = bb_atom_map.shape[0]
        if bb_atom_idx.shape != (n_bb,):
            raise ValueError("hybrid_bb_map/bb_atom_index must have shape (n_bb,)")
        bb_index_space = decode_attr_string(bb.attrs.get("atom_index_space", "runtime_n_atom"))
        bb_reference_space = bb_index_space == "protein_aa_pdb_0based"

        if "reference_atom_names" in bb:
            ref_names = bb["reference_atom_names"][:]
            if ref_names.shape != (4,):
                raise ValueError("hybrid_bb_map/reference_atom_names must have shape (4,)")
        if "reference_atom_indices" in bb:
            ref_idx = bb["reference_atom_indices"][:]
            if ref_idx.shape != (n_bb, 4):
                raise ValueError("hybrid_bb_map/reference_atom_indices must have shape (n_bb,4)")
        if "reference_atom_coords" in bb:
            ref_xyz = bb["reference_atom_coords"][:]
            if ref_xyz.shape != (n_bb, 4, 3):
                raise ValueError("hybrid_bb_map/reference_atom_coords must have shape (n_bb,4,3)")
        if "bb_comment" in bb:
            bb_comment = bb["bb_comment"][:]
            if bb_comment.shape != (n_bb,):
                raise ValueError("hybrid_bb_map/bb_comment must have shape (n_bb,)")

        for i in range(n_bb):
            mask = bb_mask[i].astype(bool)
            if mask.any():
                wsum = float(bb_w[i][mask].sum())
                if abs(wsum - 1.0) > 1e-4:
                    raise ValueError(f"BB weights do not sum to 1 for row {i}: {wsum}")

        membership = require_dataset(env, "protein_membership")[:]
        if membership.ndim != 1:
            raise ValueError("hybrid_env_topology/protein_membership must be 1D")
        if args.n_atom is not None and membership.shape[0] != args.n_atom:
            raise ValueError(
                f"protein_membership length ({membership.shape[0]}) != expected n_atom ({args.n_atom})"
            )
        n_atom = int(membership.shape[0])

        for i in range(n_bb):
            bb_i = int(bb_atom_idx[i])
            if bb_i < 0 or bb_i >= n_atom:
                raise ValueError(f"BB proxy index out of bounds at row {i}: {bb_i}")
            if membership[bb_i] < 0:
                raise ValueError(f"BB proxy index is not protein atom at row {i}: {bb_i}")
            for j in range(4):
                if int(bb_mask[i, j]) == 0:
                    continue
                ai = int(bb_atom_map[i, j])
                if bb_reference_space:
                    if ai < 0:
                        raise ValueError(f"BB reference target index invalid at row {i}, col {j}: {ai}")
                else:
                    if ai < 0 or ai >= n_atom:
                        raise ValueError(f"BB target index out of bounds at row {i}, col {j}: {ai}")
                    if membership[ai] < 0:
                        raise ValueError(f"BB target index is not protein atom at row {i}, col {j}: {ai}")

        n_protein = int(np.sum(membership >= 0))
        n_env = int(np.sum(membership < 0))
        print(f"OK: {p}")
        print(f"  n_atom={membership.shape[0]} n_protein={n_protein} n_env={n_env}")
        print(f"  n_bb={n_bb} bb_index_space={bb_index_space}")


if __name__ == "__main__":
    main()
