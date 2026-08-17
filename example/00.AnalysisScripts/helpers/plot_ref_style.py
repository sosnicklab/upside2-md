#!/usr/bin/env python3
"""Render residue-vs-dG in the reference layout (Downloads/79HIS_0.85.png).

Reads the per-temperature profiles that calc_hdx_ht.py saves and draws the unclipped profile: an amide that
never exchanges in the reweighted ensemble runs off the top of the axis, so a non-exchanging transmembrane
helix reads as one continuous excursion rather than a capped plateau. That is how these profiles are
conventionally read, and it is only correct once the membrane term is in the protection state -- without it
(martini_hdx_membrane_accessibility.py, findings 113) a lipid-facing amide drops out of the protected state
every time its backbone H-bond flickers, and the same rendering breaks each helix into isolated needles.

Everything below the excursions is the value the estimator actually returns, up to ~21 kcal/mol. Values above
0.001987 * temp_scale * ln(ESS) (4.7-6.0 kcal/mol on the glpG POPE/POPG ladder, printed per temperature) rest
on a reweighted (1 - p_f) smaller than one effective frame, so they order amides correctly but must be read as
lower bounds, not quoted as dG.

The figure carries nothing but the profiles: bare axis, one small legend of temperatures, no shading and no
marker variants. Annotation was tried on it and blocks the region the eye needs. The two things worth knowing
per region are printed instead -- the DSSP helix bounds, and the medians split into helix interior, helix
N-terminal cap and loop. Both matter for reading the figure. The secondary structure bundled in a hybrid
trajectory (`hybrid_bb_map/bb_secondary_structure`) is an idealised pattern, not DSSP, and merges glpG's 20-23
residue transmembrane helices into 41- and 50-residue blocks; and an amide in the first four positions of a
helix has no i-4 carbonyl to bond to, so it exchanges fast for a reason that has nothing to do with fold
quality. Without those two, both read as broken helix.
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md

Y_LIMITS = (-20.0, 30.0)  # the axis the reference figures use, so the two can be read side by side
STYLE = [("blue", "o"), ("red", "s"), ("green", "^"), ("darkorange", "D"), ("purple", "v")]
HELIX_CODES = ("H", "G", "I")
MIN_HELIX = 5             # shorter runs are single turns, not helices worth reporting
N_CAP = 4                 # an amide this far into a helix is the first one with an i-4 acceptor
SENTINEL = 100.0          # calc_hdx_ht.py marks a non-exchanging amide with +1000


def helix_segments(structure_path, n_residue):
    """DSSP helices of the reference coordinates, as inclusive (start, end) residue-index pairs."""
    structure = md.load(structure_path)
    if structure.n_residues < n_residue:
        raise ValueError('structure %s has %d residues, fewer than the %d the profile covers'
                         % (structure_path, structure.n_residues, n_residue))
    code = md.compute_dssp(structure, simplified=False)[0]
    helical = np.isin(code, HELIX_CODES)
    segments = []
    start = None
    for index, flag in enumerate(helical):
        if flag and start is None:
            start = index
        elif not flag and start is not None:
            if index - start >= MIN_HELIX:
                segments.append((start, index - 1))
            start = None
    if start is not None and len(helical) - start >= MIN_HELIX:
        segments.append((start, len(helical) - 1))
    return segments


def cap_mask(residue, segments):
    """True for donors in the first N_CAP positions of a helix, which have no i-4 acceptor."""
    mask = np.zeros(residue.size, dtype=bool)
    for start, end in segments:
        mask |= (residue >= start) & (residue < min(start + N_CAP, end + 1))
    return mask


def interior_mask(residue, segments):
    """True for donors inside a helix past its N-terminal cap -- where protection is expected."""
    mask = np.zeros(residue.size, dtype=bool)
    for start, end in segments:
        mask |= (residue >= start + N_CAP) & (residue <= end)
    return mask


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("npz", help="<pdb>_dG_profiles.npz from calc_hdx_ht.py")
    parser.add_argument("--out", required=True, help="figure to write")
    parser.add_argument("--temperatures", default="0.75,0.80,0.85",
                        help="(default 0.75,0.80,0.85) which tabulated temperatures to draw")
    parser.add_argument("--structure", default=None,
                        help="reference PDB; when given, the DSSP helix bounds and the helix "
                             "interior/cap/loop medians are printed")
    parser.add_argument("--title", default="Residue ID vs Delta G at Different Temperatures")
    args = parser.parse_args()

    d = np.load(args.npz)
    for required in ("dG_unclipped", "dg_limit"):
        if required not in d.files:
            raise ValueError('%s has no %r; re-run 4.calc_D_uptake.py to regenerate it'
                             % (args.npz, required))
    residue = d["res"]
    temps = np.asarray(d["temps"], dtype=float)
    dG = d["dG_unclipped"]
    dg_limit = np.asarray(d["dg_limit"], dtype=float)
    wanted = [float(x) for x in args.temperatures.split(",")]

    figure, axis = plt.subplots(figsize=(15, 10))
    for position, target in enumerate(wanted):
        j = int(np.argmin(np.abs(temps - target)))
        colour, marker = STYLE[position % len(STYLE)]
        y = np.array(dG[j], dtype=float)
        # A non-exchanging amide is joined to a point far off the axis, which draws the vertical excursion
        # the figure is read by; no marker is placed there, since there is no value to mark.
        high = y >= SENTINEL
        low = y <= -SENTINEL
        drawn = np.where(high, Y_LIMITS[1] * 3.0, np.where(low, Y_LIMITS[0] * 3.0, y))
        resolved = ~high & ~low
        axis.plot(residue, drawn, "-", color=colour, linewidth=2.0, zorder=2)
        axis.plot(residue[resolved], y[resolved], marker, color=colour, markersize=8,
                  linestyle="none", label="T = %g" % temps[j], zorder=3)

    axis.set_xlabel("Residue ID", fontsize=26)
    axis.set_ylabel("DG (kcal/mol)", fontsize=26)
    axis.set_title(args.title, fontsize=27)
    axis.set_ylim(*Y_LIMITS)
    axis.set_xlim(residue.min() - 6, residue.max() + 6)
    axis.grid(True, color="0.85", linewidth=0.8)
    axis.set_axisbelow(True)
    axis.tick_params(labelsize=20)
    handles, labels = axis.get_legend_handles_labels()
    # Lower left, not the reference figure's upper left. On this axis nothing falls below -1 kcal/mol, so the
    # bottom corner is the one region guaranteed to be empty; the reference's own corner sits under the
    # excursions off helix 1 here.
    axis.legend(handles, labels, fontsize=22, loc="lower left", framealpha=1.0)
    figure.tight_layout()
    figure.savefig(args.out, dpi=200)
    print("wrote %s" % args.out)

    segments = helix_segments(args.structure, int(residue.max()) + 1) if args.structure else []
    if segments:
        cap = cap_mask(residue, segments)
        interior = interior_mask(residue, segments)
        print("  DSSP helices: %s" % ", ".join("%d-%d" % s for s in segments))
        print("  donors: %d helix interior, %d helix N-terminal cap, %d loop"
              % (int(interior.sum()), int(cap.sum()), int((~interior & ~cap).sum())))
    for target in wanted:
        j = int(np.argmin(np.abs(temps - target)))
        y = np.array(dG[j], dtype=float)
        resolved = np.abs(y) < SENTINEL
        print("  T=%.2f: %d of %d non-exchanging (off scale), %d resolved, resolution limit %.2f"
              % (temps[j], int((~resolved).sum()), y.size, int(resolved.sum()), dg_limit[j]))
        if resolved.any():
            print("          resolved dG min %.2f max %.2f median %.2f"
                  % (y[resolved].min(), y[resolved].max(), np.median(y[resolved])))
        if segments:
            deep = resolved & interior
            loop = resolved & ~interior & ~cap
            print("          helix interior median %.2f | cap median %.2f | loop median %.2f"
                  % (np.median(y[deep]) if deep.any() else np.nan,
                     np.median(y[resolved & cap]) if (resolved & cap).any() else np.nan,
                     np.median(y[loop]) if loop.any() else np.nan))
            outlier = residue[deep & (y < 0.0)]
            print("          helix-interior amides below 0 kcal/mol: %s"
                  % (outlier.tolist() if outlier.size else "none"))


if __name__ == "__main__":
    main()
