#!/usr/bin/env python3
"""Upgrade a hybrid .up written before the BB proxy moved onto Upside's derived-site path.

The BB proxy used to be placed from a stored reference frame, so `martini_hybrid_position` consumed only
`pos`. It is now the mass-weighted centre of N/CA/C/O with O taken from `infer_H_O`, which makes it a
linear combination of nodes the engine differentiates and leaves the proxy with no inertia of its own.
The node therefore takes a second argument, and a file written by the old prep declares only the first --
the engine rejects it with "expected 1 arguments but got 2".

Nothing else in those files needs changing: `infer_H_O` is already present in every hybrid config, and the
BB map still supplies the carrier indices and mass weights. So this rewrites one attribute.

Idempotent, and refuses a file that has no infer_H_O node rather than producing a config that cannot load.
"""
import argparse
import sys

import h5py
import numpy as np

NODE = "/input/potential/martini_hybrid_position"
WANT = [b"pos", b"infer_H_O"]


def upgrade(path, dry_run=False):
    with h5py.File(path, "r" if dry_run else "r+") as h5:
        if NODE not in h5:
            return "not a hybrid config"
        node = h5[NODE]
        current = [x if isinstance(x, bytes) else str(x).encode() for x in node.attrs.get("arguments", [])]
        if current == WANT:
            return "already upgraded"
        if current != [b"pos"]:
            raise SystemExit("%s: unexpected arguments %s" % (path, current))
        if "/input/potential/infer_H_O" not in h5:
            raise SystemExit("%s: has no infer_H_O node to source the backbone O from" % path)
        if dry_run:
            return "would upgrade"
        node.attrs["arguments"] = np.array(WANT)
        return "upgraded"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+", help=".up files to upgrade in place")
    parser.add_argument("--dry-run", action="store_true", help="report what would change and write nothing")
    args = parser.parse_args()

    for path in args.files:
        print("%-58s %s" % (path.split("/")[-1], upgrade(path, args.dry_run)), flush=True)


if __name__ == "__main__":
    main()
