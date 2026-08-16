#!/usr/bin/env python3
"""Render residue-vs-dG in the reference layout (Downloads/79HIS_0.90.png).

Reads the per-temperature profiles that calc_hdx_ht.py saves and draws them the way the reference figure
does: one series per temperature with its own marker, and amides that never exchange (the +1000 sentinel,
p_f = 1) drawn as lines running off the top of the axis rather than as markers parked on the edge, so a
non-exchanging residue reads as off-scale instead of as a data point at y = 30.
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SENTINEL = 100.0          # calc_hdx_ht.py uses +1000 for p_f == 1 and -100 for p_f == 0
Y_LIMITS = (-20.0, 30.0)
STYLE = [("blue", "o"), ("red", "s"), ("green", "^"), ("darkorange", "D"), ("purple", "v")]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("npz", help="<pdb>_dG_profiles.npz from calc_hdx_ht.py")
    parser.add_argument("--out", required=True, help="figure to write")
    parser.add_argument("--temperatures", default="0.75,0.80,0.85",
                        help="(default 0.75,0.80,0.85) which tabulated temperatures to draw")
    parser.add_argument("--title", default="Residue ID vs Delta G at Different Temperatures")
    args = parser.parse_args()

    d = np.load(args.npz)
    residue, temps, dG = d["res"], np.asarray(d["temps"], dtype=float), d["dG"]
    wanted = [float(x) for x in args.temperatures.split(",")]

    figure, axis = plt.subplots(figsize=(14, 9))
    for position, target in enumerate(wanted):
        j = int(np.argmin(np.abs(temps - target)))
        colour, marker = STYLE[position % len(STYLE)]
        y = np.array(dG[j], dtype=float)
        # Draw through the sentinels so they leave the axis as vertical excursions, which is what the
        # reference does; the y-limit then clips them instead of them being plotted at the boundary.
        drawn = np.where(y >= SENTINEL, Y_LIMITS[1] * 3.0,
                         np.where(y <= -SENTINEL, Y_LIMITS[0] * 3.0, y))
        finite = np.abs(y) < SENTINEL
        axis.plot(residue, drawn, "-", color=colour, linewidth=1.4, zorder=2)
        axis.plot(residue[finite], y[finite], marker, color=colour, markersize=5,
                  linestyle="none", label="T = %g" % temps[j], zorder=3)

    axis.set_xlabel("Residue ID", fontsize=16)
    axis.set_ylabel("DG (kcal/mol)", fontsize=16)
    axis.set_title(args.title, fontsize=17)
    axis.set_ylim(*Y_LIMITS)
    axis.set_xlim(residue.min() - 6, residue.max() + 6)
    axis.grid(True, alpha=0.25)
    axis.tick_params(labelsize=13)
    axis.legend(fontsize=14, loc="upper left")
    figure.tight_layout()
    figure.savefig(args.out, dpi=200)
    print("wrote %s" % args.out)
    for position, target in enumerate(wanted):
        j = int(np.argmin(np.abs(temps - target)))
        y = np.array(dG[j], dtype=float)
        fin = y[np.abs(y) < SENTINEL]
        if fin.size == 0:
            # Every amide off-scale means the reweighting resolved nothing at this target temperature,
            # which is a result worth stating rather than a crash in the summary line.
            print("  T=%.2f: 0 of %d exchange | ALL residues off-scale (%d at p_f = 1, %d at p_f = 0)"
                  % (temps[j], y.size, int((y >= SENTINEL).sum()), int((y <= -SENTINEL).sum())))
            continue
        print("  T=%.2f: %d of %d exchange | dG min %.2f max %.2f median %.2f | never-exchanging %d"
              % (temps[j], fin.size, y.size, fin.min(), fin.max(), np.median(fin),
                 int((y >= SENTINEL).sum())))


if __name__ == "__main__":
    main()
