#!/usr/bin/env python3
"""Join the chunked output groups of a chained REMD run into one analysis trajectory.

A chained run does not leave its trajectory where an analysis can find it. `run_remd.py` reseeds every
chunk from the last frame and renames `output` to the lowest free `output_previous_<n>`, so after a
36 h block `/output` holds only the final ~300 frames and everything else is scattered across history
groups. The HDX tools (`martini_hdx_project.py`, `get_info_from_upside_traj.py`) read `/output`, so the
chunks have to be concatenated first -- this is the step the DDM campaign did by hand as "trimming".

Chunks are contiguous: each starts from the previous one's last frame, so concatenating them in restart
order reproduces the trajectory. Velocities are re-thermalised at each boundary, which is why `time`
restarts at zero per chunk and is renumbered here rather than copied.

A chunk whose potential is not finite anywhere is dropped whole, and which ones were dropped is printed.
That is not a threshold: a rolled-back chunk is a chunk the driver already declared destroyed and
re-ran, so its frames are not samples of anything and averaging them in would be wrong. Nothing is
clipped, rescaled or repaired -- a chunk is either sound or absent.
"""
import argparse
import re

import h5py
import numpy as np


def chunk_groups(h5):
    """Every completed chunk plus the live one, in restart order."""
    prev = sorted((int(m.group(1)), name) for name in h5["/"].keys()
                  if (m := re.fullmatch(r"output_previous_(\d+)", name)))
    return [name for _, name in prev] + (["output"] if "output" in h5 else [])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", help="chained REMD .up")
    parser.add_argument("dest", help="analysis .up to write (/input copied, /output concatenated)")
    parser.add_argument("--discard-frames", type=int, default=0,
                        help="(default 0) drop this many frames from the front of the joined trajectory")
    parser.add_argument("--stride", type=int, default=1, help="(default 1) keep every Nth frame")
    args = parser.parse_args()

    with h5py.File(args.source, "r") as src, h5py.File(args.dest, "w") as dst:
        src.copy("/input", dst, name="input")
        groups = chunk_groups(src)
        if not groups:
            raise SystemExit("%s has no output groups" % args.source)

        keep, dropped = [], []
        for g in groups:
            pot = np.asarray(src[g]["potential"]).reshape(-1)
            (keep if pot.size and np.isfinite(pot).all() else dropped).append(g)
        if not keep:
            raise SystemExit("%s: every chunk is non-finite" % args.source)

        names = [k for k in src[keep[0]].keys() if all(k in src[g] for g in keep)]
        out = dst.create_group("output")
        total = 0
        # Appended chunk by chunk rather than concatenated in one call: a full 12-block chain is several
        # GB per replica for `pos` alone, and 48 of these run in parallel.
        for name in names:
            written = 0
            seen = 0                      # frames passed over, so the stride stays continuous across chunks
            dset = None
            for g in keep:
                part = np.asarray(src[g][name])
                # The stride is over the joined trajectory, so its phase has to carry across a chunk
                # boundary; restarting it at each chunk would resample the boundary frames.
                if seen >= args.discard_frames:
                    start = (args.discard_frames - seen) % args.stride
                else:
                    start = args.discard_frames - seen
                take = part[start::args.stride] if start < part.shape[0] else part[:0]
                seen += part.shape[0]
                if take.shape[0] == 0:
                    continue
                if dset is None:
                    dset = out.create_dataset(name, data=take,
                                              maxshape=(None,) + take.shape[1:], chunks=True)
                    written = take.shape[0]
                else:
                    dset.resize(written + take.shape[0], axis=0)
                    dset[written:] = take
                    written += take.shape[0]
            total = max(total, written)

        if "time" in out:
            del out["time"]
            out.create_dataset("time", data=np.arange(total, dtype=np.float32))

        print("%s -> %s: %d chunks kept, %d frames%s"
              % (args.source.split("/")[-1], args.dest.split("/")[-1], len(keep), total,
                 "" if not dropped else "; DROPPED non-finite %s" % ",".join(dropped)), flush=True)


if __name__ == "__main__":
    main()
