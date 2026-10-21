#!/usr/bin/env python
"""Per-amide membrane (water) accessibility for hybrid dry-MARTINI HDX.

The stock protein-only protection (get_protection_state.py) scores an amide as protected only if it is
backbone/side-chain H-bonded or buried by other *protein* atoms. For a membrane protein that is wrong:
a transmembrane amide facing the lipid is shielded from water by the lipid tails, yet the protein-only
calculation sees no burial there, so it drops out of the protected state whenever its backbone H-bond
flickers -- breaking the continuous +inf that a non-exchanging TM helix should show. Master's equivalent
for an implicit membrane is get_protection_state.py --use-TM-region, which reads the `surface` node; that
node is implicit-membrane machinery and does not exist in a hybrid HDX topology.

This script recovers the missing term from the *full* hybrid trajectory, which still has the lipid: an
amide backbone N with at least one lipid hydrophobic-tail bead inside the tails' own first coordination
shell is treated as water-INACCESSIBLE. The output is a (n_frame, n_donor) accessibility array
(1 = water-accessible, 0 = inaccessible) aligned with get_protection_state's donor ordering, to be fed to
combine_hdx_protection.py --water-accessibility (which sets
protection = 1 - (1 - protein_protection) * accessibility).

The criterion is the bilayer boundary, and it is measured from the bilayer rather than chosen. The radius
is the flat first minimum of the intermolecular tail-tail g(r) of this system's own lipids, so "inside the
radius" means "inside the tail region's first neighbour shell"; on the glpG POPE/POPG patch that comes out
at 7.00 A on every replica of the ladder. One contact is then the boundary itself: measured at that radius,
a phosphate bead -- the conventional edge of the bilayer -- has a median of zero intermolecular tail
neighbours (ester beads 2, first tail bead 6, terminal tail bead 11), so the first tail contact is the
first step past the head groups.

The test is deliberately local rather than a slab between the two leaflet planes. A slab also protects
every amide lining the protein's own polar interior, which for glpG merges the short interfacial loops
into the helices and gives contiguous +inf runs of 40-50 residues. Master ANDs its slab with a
lipid-facing surface test for the same reason.

Donor ordering and frame set are taken to match get_protection_state.py: donors come from the HDX
topology's infer_H_O.donors.residue, and frames are the concatenated output groups (output_previous_0..N
then output) with the same --start/--stride.
"""
import argparse

import numpy as np
import re
import tables as tb
from scipy.spatial import cKDTree

RDF_EDGES = np.arange(2.0, 12.01, 0.25)   # tail-tail g(r) bins used to locate the first coordination shell
RDF_FRAMES = 60                           # frames subsampled for the g(r); it converges well before this
RDF_FLAT = 1.05                           # bins within this factor of the minimum count as part of it


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


def tail_mask(is_lipid, atom_names, tail_names):
    """Hydrophobic-tail beads, by explicit name or by the MARTINI apolar naming pattern."""
    if tail_names:
        return is_lipid & np.isin(atom_names, tail_names), "names " + ",".join(tail_names)
    # MARTINI names an apolar tail bead C<n> or D<n> (D marks the unsaturated one), optionally suffixed
    # with a chain letter: C1/C2/C3 for a single-tailed detergent, C1A..C4A/C1B..C5B and D3B for a
    # two-tailed lipid. Detecting the pattern rather than defaulting to one system's bead list keeps this
    # usable for any lipid; a lipid that breaks the pattern is named with --tail-bead-names.
    pattern = re.compile(r"^[CD]\d[A-Z]?$")
    is_tail = is_lipid & np.array([bool(pattern.match(n)) for n in atom_names])
    return is_tail, "auto-detected " + ",".join(sorted(set(atom_names[is_tail])))


def contact_radius(pos, box, tail_idx, tail_molecule):
    """The flat first minimum of the intermolecular tail-tail g(r) -- the tail region's first shell."""
    frames = pos[:: max(1, pos.shape[0] // RDF_FRAMES)]
    histogram = np.zeros(RDF_EDGES.size - 1)
    for f in range(frames.shape[0]):
        wrapped = np.mod(frames[f], box)
        tree = cKDTree(wrapped[tail_idx], boxsize=box)
        pairs = tree.query_pairs(RDF_EDGES[-1], output_type="ndarray")
        pairs = pairs[tail_molecule[pairs[:, 0]] != tail_molecule[pairs[:, 1]]]
        delta = wrapped[tail_idx[pairs[:, 0]]] - wrapped[tail_idx[pairs[:, 1]]]
        delta -= box * np.round(delta / box)
        histogram += np.histogram(np.linalg.norm(delta, axis=1), RDF_EDGES)[0]

    centre = 0.5 * (RDF_EDGES[1:] + RDF_EDGES[:-1])
    shell = 4.0 * np.pi / 3.0 * (RDF_EDGES[1:] ** 3 - RDF_EDGES[:-1] ** 3)
    g = histogram / shell
    peak = int(np.argmax(g))
    beyond = g[peak:]
    # The minimum spans two bins and which of them is lowest flips with sampling noise, which would give
    # different replicas of one ladder different boundaries. Averaging over the bins within RDF_FLAT of
    # the minimum makes the radius reproducible: 7.00 A on all 16 glpG POPE/POPG replicas.
    flat = np.where(beyond <= beyond.min() * RDF_FLAT)[0] + peak
    return float(0.5 * (centre[flat[0]] + centre[flat[-1]])), float(centre[peak])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("hdx_topology", help="protein-only HDX topology .up (for donor residues)")
    parser.add_argument("hybrid_trajectory", help="full hybrid trajectory .up (protein + lipid)")
    parser.add_argument("output_npy", help="output accessibility array (n_frame, n_donor)")
    parser.add_argument("--structure-up", default=None,
                        help="source of the periodic box (default: hybrid_trajectory)")
    parser.add_argument("--lipid-class", default="LIPID", help="(default LIPID) particle_class of lipid beads")
    parser.add_argument("--tail-bead-names", default="",
                        help="comma-separated hydrophobic-tail bead names; by default they are detected "
                             "from the trajectory as the MARTINI apolar naming pattern C<n>/D<n> with an "
                             "optional chain letter (C1, C3 for a detergent; C1A, D3B, C5B for a "
                             "two-tailed lipid)")
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

    is_lipid = particle_class == args.lipid_class
    is_tail, selected = tail_mask(is_lipid, atom_names, tail_names)
    tail_idx = np.where(is_tail)[0]
    if tail_idx.size == 0:
        raise ValueError(
            "no lipid tail beads found (class {}, {}); lipid bead names present: {}".format(
                args.lipid_class, selected, ",".join(sorted(set(atom_names[is_lipid])))))
    print("tail beads: {} ({} beads)".format(selected, tail_idx.size))

    radius, peak = contact_radius(pos, box, tail_idx, residue_ids[tail_idx])
    print("tail-tail g(r): first peak {:.2f} A, flat first minimum {:.2f} A".format(peak, radius))

    accessibility = np.ones((pos.shape[0], donor.size), dtype=np.float32)
    for f in range(pos.shape[0]):
        wrapped = np.mod(pos[f], box)
        tree = cKDTree(wrapped[tail_idx], boxsize=box)
        counts = tree.query_ball_point(wrapped[amide_idx], r=radius, return_length=True)
        accessibility[f, np.asarray(counts) > 0] = 0.0

    output_base = args.output_npy[:-4] if args.output_npy.endswith(".npy") else args.output_npy
    np.save(output_base, accessibility)
    buried = (accessibility == 0.0).any(axis=0).sum()
    print("wrote {} : {} frames x {} donors; {} donors lipid-buried in >=1 frame, {:.1%} of amide-frames"
          .format(output_base + ".npy", pos.shape[0], donor.size, int(buried), 1.0 - accessibility.mean()))


if __name__ == "__main__":
    main()
