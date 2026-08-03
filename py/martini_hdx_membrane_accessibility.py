#!/usr/bin/env python
"""Per-amide membrane (water) accessibility for hybrid dry-MARTINI HDX.

The stock protein-only protection (get_protection_state.py) scores an amide as protected only if it is
backbone/side-chain H-bonded or buried by other *protein* atoms. For a membrane protein that is wrong:
a transmembrane amide facing the lipid is shielded from water by the lipid tails, yet the protein-only
calculation sees no burial there, so it drops out of the protected state whenever its backbone H-bond
flickers -- breaking the continuous +inf that a non-exchanging TM helix should show.

This script recovers the missing term from the *full* hybrid trajectory (which still has the lipid):
for each amide backbone N, per frame, it counts lipid hydrophobic-tail beads within `cutoff`; if at
least `min_contacts` are present the amide is treated as water-INACCESSIBLE (buried in the tail core).
The output is a (n_frame, n_donor) accessibility array (1 = water-accessible, 0 = inaccessible) aligned
with get_protection_state's donor ordering, to be fed to combine_hdx_protection.py --water-accessibility
(which sets protection = 1 - (1 - protein_protection) * accessibility).

Donor ordering and frame set are taken to match get_protection_state.py: donors come from the HDX
topology's infer_H_O.donors.residue, and frames are the concatenated output groups (output_previous_0..N
then output) with the same --start/--stride.
"""
import argparse
import os

import numpy as np
import tables as tb
from scipy.spatial import cKDTree


def output_group_names(h5):
    names = []
    i = 0
    while hasattr(h5.root, "output_previous_{}".format(i)):
        names.append("output_previous_{}".format(i))
        i += 1
    if hasattr(h5.root, "output"):
        names.append("output")
    return names


def read_all_positions(traj_h5, start, stride):
    # Mirror mdtraj_upside.load_upside_traj (what get_protection_state reads): concatenate the output
    # groups, dropping the first frame of every group after the first (it duplicates the pre-restart
    # frame), then apply start/stride. This keeps the accessibility frame set aligned with the
    # protein-protection frame set so combine_hdx_protection can multiply them element-wise.
    names = output_group_names(traj_h5)
    if not names:
        raise ValueError("hybrid trajectory has no output groups")
    frames = []
    for g_no, name in enumerate(names):
        pos = traj_h5.get_node("/" + name).pos[:, 0, :, :]  # (n_frame, n_atom, 3)
        frames.append(pos[1:] if g_no else pos)
    return np.concatenate(frames, axis=0)[start::stride]


def box_lengths(struct_h5):
    g = struct_h5.root.input.potential.martini_potential
    return np.array([float(g._v_attrs.x_len), float(g._v_attrs.y_len), float(g._v_attrs.z_len)])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("hdx_topology", help="protein-only HDX topology .up (for donor residues)")
    parser.add_argument("hybrid_trajectory", help="full hybrid trajectory .up (protein + lipid)")
    parser.add_argument("output_npy", help="output accessibility array (n_frame, n_donor)")
    parser.add_argument("--structure-up", default=None,
                        help="source of the periodic box (default: hybrid_trajectory)")
    parser.add_argument("--lipid-class", default="LIPID", help="(default LIPID) particle_class of lipid beads")
    parser.add_argument("--tail-bead-names", default="C1,C2,C3",
                        help="(default C1,C2,C3) comma-separated hydrophobic-tail bead names")
    parser.add_argument("--cutoff", type=float, default=8.0,
                        help="(default 8.0 A) tail-bead contact radius around each amide N")
    parser.add_argument("--min-contacts", type=int, default=5,
                        help="(default 5) tail beads within cutoff for the amide to be water-inaccessible")
    parser.add_argument("--start", type=int, default=0, help="(default 0) initial frame (match get_protection_state)")
    parser.add_argument("--stride", type=int, default=1, help="(default 1) frame stride (match get_protection_state)")
    args = parser.parse_args()

    tail_names = [s.strip() for s in args.tail_bead_names.split(",") if s.strip()]
    structure_path = args.structure_up or args.hybrid_trajectory

    with tb.open_file(args.hdx_topology, "r") as top:
        donor = top.root.input.potential.infer_H_O.donors.residue[:]

    with tb.open_file(structure_path, "r") as s:
        box = box_lengths(s)

    with tb.open_file(args.hybrid_trajectory, "r") as t:
        particle_class = np.array([x.decode() for x in t.root.input.particle_class[:]])
        atom_names = np.array([x.decode() for x in t.root.input.atom_names[:]])
        residue_ids = t.root.input.residue_ids[:]
        pos = read_all_positions(t, args.start, args.stride)

    is_amide_N = (particle_class == "PROTEIN") & (atom_names == "N")
    res_to_N = {int(residue_ids[i]): int(i) for i in np.where(is_amide_N)[0]}
    missing = [int(d) for d in donor if int(d) not in res_to_N]
    if missing:
        raise ValueError("no backbone N for donor residues {}".format(missing[:10]))
    amide_idx = np.array([res_to_N[int(d)] for d in donor])

    is_tail = (particle_class == args.lipid_class) & np.isin(atom_names, tail_names)
    tail_idx = np.where(is_tail)[0]
    if tail_idx.size == 0:
        raise ValueError("no lipid tail beads found (class {} names {})".format(args.lipid_class, tail_names))

    n_frame = pos.shape[0]
    accessibility = np.ones((n_frame, donor.size), dtype=np.float32)
    for f in range(n_frame):
        wrapped = np.mod(pos[f], box)
        tree = cKDTree(wrapped[tail_idx], boxsize=box)
        counts = tree.query_ball_point(wrapped[amide_idx], r=args.cutoff, return_length=True)
        accessibility[f, np.asarray(counts) >= args.min_contacts] = 0.0

    output_base = args.output_npy[:-4] if args.output_npy.endswith(".npy") else args.output_npy
    np.save(output_base, accessibility)
    buried = (accessibility == 0.0).any(axis=0).sum()
    print("wrote {} : {} frames x {} donors; {} donors lipid-buried in >=1 frame".format(
        output_base + ".npy", n_frame, donor.size, int(buried)))


if __name__ == "__main__":
    main()
