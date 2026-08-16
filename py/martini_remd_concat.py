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

Destroyed frames are dropped individually, not by the chunk. A blow-up usually happens partway through a
chunk -- the local ladder lost its top rung 2333 frames into a 2448-frame segment -- so discarding the
whole chunk would throw away thousands of perfectly good samples to remove a hundred bad ones. The test is
the physical one used to pick restart points: a frame is kept if its potential is finite and negative. A
condensed bilayer plus protein sits near -2.2e4 E_up, so a positive total potential means overlapping
cores, and such a frame is not a sample of the ensemble. Note that finiteness alone is not enough -- the
same replica stayed finite for 96 frames at +1.9e6 before reaching NaN. Nothing is clipped, rescaled or
repaired; a frame is either sound or absent, and the count removed is reported.

Chunks at a different temperature are dropped for the same reason. The production seeds carry an
`/output` from their equilibration stage, which the first reseed rotates to `output_previous_0`, so the
oldest chunk of every replica is an equilibration run at a single shared temperature rather than the
replica's rung. Joining it on would put samples from two thermodynamic states into one trajectory, and --
because `get_info_from_upside_traj.py` reads the temperature from the first frame -- would label all 48
replicas with the equilibration temperature, collapsing the MBAR ladder to one state and returning uniform
weights. The production temperature is taken from the newest chunk and every chunk must match it.
"""
import argparse
import re

import h5py
import numpy as np


def chunk_groups(h5, completed_only=False):
    """Chunks in restart order; `completed_only` drops the live `output` group.

    On a job that is still running, the driver's next reseed renames `output` to the lowest free
    `output_previous_<n>`, so that group can disappear between being listed and being read. Analysing a
    live run therefore has to skip it -- the completed chunks are never renamed again.
    """
    prev = sorted((int(m.group(1)), name) for name in h5["/"].keys()
                  if (m := re.fullmatch(r"output_previous_(\d+)", name)))
    groups = [name for _, name in prev]
    if not completed_only and "output" in h5:
        groups.append("output")
    return groups


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", help="chained REMD .up")
    parser.add_argument("dest", help="analysis .up to write (/input copied, /output concatenated)")
    parser.add_argument("--discard-frames", type=int, default=0,
                        help="(default 0) drop this many frames from the front of the joined trajectory")
    parser.add_argument("--stride", type=int, default=1, help="(default 1) keep every Nth frame")
    parser.add_argument("--completed-only", action="store_true",
                        help="skip the live `output` group; required when the run is still going")
    args = parser.parse_args()

    with h5py.File(args.source, "r") as src, h5py.File(args.dest, "w") as dst:
        src.copy("/input", dst, name="input")
        groups = chunk_groups(src, args.completed_only)
        if not groups:
            raise SystemExit("%s has no output groups" % args.source)

        def chunk_temperature(g):
            if "temperature" not in src[g]:
                return None
            t = np.asarray(src[g]["temperature"]).reshape(-1)
            return None if t.size == 0 else float(t[-1])

        production_t = chunk_temperature(groups[-1])

        keep, healthy, off_temp, n_dropped = [], {}, [], 0
        for g in groups:
            t = chunk_temperature(g)
            if production_t is not None and t is not None and abs(t - production_t) > 1e-6:
                off_temp.append("%s@%.4f" % (g, t))
                continue
            pot = np.asarray(src[g]["potential"]).reshape(-1)
            ok = np.isfinite(pot) & (pot < 0.0)
            if not ok.any():
                n_dropped += pot.size
                continue
            n_dropped += int((~ok).sum())
            keep.append(g)
            healthy[g] = np.where(ok)[0]
        if not keep:
            raise SystemExit("%s: no healthy on-temperature frames" % args.source)

        # Only per-frame datasets can be joined. A chunk also carries a few whole-run records whose first
        # axis is not the frame count (e.g. a single row of replica bookkeeping); indexing those with frame
        # numbers is meaningless, so they are left out rather than mangled.
        n_frames = {g: np.asarray(src[g]["potential"]).reshape(-1).size for g in keep}
        names = [k for k in src[keep[0]].keys()
                 if all(k in src[g] and src[g][k].shape[0] == n_frames[g] for g in keep)]
        skipped = sorted(set(src[keep[0]].keys()) - set(names))
        out = dst.create_group("output")
        total = 0
        # Appended chunk by chunk rather than concatenated in one call: a full 12-block chain is several
        # GB per replica for `pos` alone, and 48 of these run in parallel.
        for name in names:
            written = 0
            seen = 0                      # frames passed over, so the stride stays continuous across chunks
            dset = None
            for g in keep:
                part = np.asarray(src[g][name])[healthy[g]]     # destroyed frames removed first
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

        print("%s -> %s: %d chunks kept at T=%s, %d frames%s%s"
              % (args.source.split("/")[-1], args.dest.split("/")[-1], len(keep),
                 "?" if production_t is None else "%.4f" % production_t, total,
                 "" if not n_dropped else "; DROPPED %d destroyed frames" % n_dropped,
                 "" if not off_temp else "; DROPPED off-temperature %s" % ",".join(off_temp)), flush=True)
        if skipped:
            print("    not per-frame, not joined: %s" % ",".join(skipped), flush=True)


if __name__ == "__main__":
    main()
