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
shell is treated as water-INACCESSIBLE. The output is one (n_frame, n_donor) accessibility array per
replica (1 = water-accessible, 0 = inaccessible) aligned with get_protection_state's donor ordering, to be
fed to combine_hdx_protection.py --water-accessibility (which sets
protection = 1 - (1 - protein_protection) * accessibility).

The criterion is the bilayer boundary, and it is measured from the bilayer rather than chosen. The radius
is the flat first minimum of the intermolecular tail-tail g(r) of this system's own lipids, so "inside the
radius" means "inside the tail region's first neighbour shell"; on the glpG POPE/POPG patch that comes out
at 7.00 A. One contact is then the boundary itself: measured at that radius, a phosphate bead -- the
conventional edge of the bilayer -- has a median of zero intermolecular tail neighbours (ester beads 2,
first tail bead 6, terminal tail bead 11), so the first tail contact is the first step past the head groups.

The radius is measured once, from the g(r) pooled over every replica of the analysis, and that one radius is
applied to all of them: MBAR averages the protection state across the ladder, so it must be one function of
the coordinates on every rung. Measured per replica it moved between 6.6 and 7.25 A on 1rkl runs, enough to
shift the shielded fraction from 0.65 to 0.79. The radius is written beside the arrays.

The test is deliberately local rather than a slab between the two leaflet planes. A slab also protects
every amide lining the protein's own polar interior, which for glpG merges the short interfacial loops
into the helices and gives contiguous +inf runs of 40-50 residues. Master ANDs its slab with a
lipid-facing surface test for the same reason.

Each replica is named by its HDX view (martini_hdx_project.py), which records the hybrid segments it was
projected from. Frames are those segments' output groups in order, the first frame of every group after the
first dropped as the view's readers drop it, with the same --start/--stride as get_protection_state.py; each
frame is evaluated in its own box (/output/box under a barostat). Donors come from the HDX topology's
infer_H_O.donors.residue, and a donor residue's amide N is found through /input/hybrid_bb_map, whose row i is
Upside residue i.
"""
import argparse
import os

import numpy as np
import re
import tables as tb
from scipy.spatial import cKDTree

RDF_EDGES = np.arange(2.0, 12.01, 0.25)   # tail-tail g(r) bins used to locate the first coordination shell
RDF_FRAMES = 60                           # frames subsampled per replica for the g(r)
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


def view_frame_count(view_h5):
    groups = [view_h5.get_node("/" + name) for name in output_group_names(view_h5)]
    return sum(g.pos.shape[0] for g in groups) - (len(groups) - 1)


def hybrid_frames(view_path, start, stride):
    """Positions and per-frame boxes of the full hybrid trajectory behind one HDX view, frame for frame."""
    with tb.open_file(view_path, "r") as view:
        sources = [x.decode() for x in view.root._v_attrs.hdx_hybrid_sources]
        n_view = view_frame_count(view)
    positions = []
    boxes = []
    for source in sources:
        with tb.open_file(source, "r") as t:
            g = t.root.input.potential.martini_potential
            input_box = np.array([float(g._v_attrs.x_len), float(g._v_attrs.y_len), float(g._v_attrs.z_len)])
            for name in output_group_names(t):
                group = t.get_node("/" + name)
                pos = group.pos[:, 0, :, :]
                box = group.box[:] if hasattr(group, "box") else np.tile(input_box, (pos.shape[0], 1))
                skip = 1 if positions else 0
                positions.append(pos[skip:])
                boxes.append(np.asarray(box, dtype=np.float64)[skip:])
    pos = np.concatenate(positions, axis=0)
    box = np.concatenate(boxes, axis=0)
    if pos.shape[0] != n_view:
        raise ValueError("{} has {} frames but its hybrid sources give {}".format(view_path, n_view, pos.shape[0]))
    return sources[0], pos[start::stride], box[start::stride]


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


def tail_tail_histogram(pos, boxes, tail_idx, tail_molecule):
    """Intermolecular tail-tail pair distances, binned on RDF_EDGES, over a subsample of the frames."""
    histogram = np.zeros(RDF_EDGES.size - 1)
    for f in range(0, pos.shape[0], max(1, pos.shape[0] // RDF_FRAMES)):
        box = boxes[f]
        wrapped = np.mod(pos[f], box)
        tree = cKDTree(wrapped[tail_idx], boxsize=box)
        pairs = tree.query_pairs(RDF_EDGES[-1], output_type="ndarray")
        pairs = pairs[tail_molecule[pairs[:, 0]] != tail_molecule[pairs[:, 1]]]
        delta = wrapped[tail_idx[pairs[:, 0]]] - wrapped[tail_idx[pairs[:, 1]]]
        delta -= box * np.round(delta / box)
        histogram += np.histogram(np.linalg.norm(delta, axis=1), RDF_EDGES)[0]
    return histogram


def first_minimum(histogram):
    """The flat first minimum of the tail-tail g(r) -- the tail region's first shell -- and its first peak."""
    centre = 0.5 * (RDF_EDGES[1:] + RDF_EDGES[:-1])
    shell = 4.0 * np.pi / 3.0 * (RDF_EDGES[1:] ** 3 - RDF_EDGES[:-1] ** 3)
    g = histogram / shell
    peak = int(np.argmax(g))
    beyond = g[peak:]
    # The minimum spans two bins and which of them is lowest flips with sampling noise. Averaging over the
    # bins within RDF_FLAT of the minimum makes the radius reproducible.
    flat = np.where(beyond <= beyond.min() * RDF_FLAT)[0] + peak
    return float(0.5 * (centre[flat[0]] + centre[flat[-1]])), float(centre[peak])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("hdx_topology", help="protein-only HDX topology .up (for donor residues)")
    parser.add_argument("output_pattern",
                        help="accessibility array path containing {replica}, the view's position in the list")
    parser.add_argument("hdx_views", nargs="+",
                        help="HDX views of every replica of the analysis (martini_hdx_project.py output)")
    parser.add_argument("--lipid-class", default="LIPID", help="(default LIPID) particle_class of lipid beads")
    parser.add_argument("--tail-bead-names", default="",
                        help="comma-separated hydrophobic-tail bead names; by default they are detected "
                             "from the trajectory as the MARTINI apolar naming pattern C<n>/D<n> with an "
                             "optional chain letter (C1, C3 for a detergent; C1A, D3B, C5B for a "
                             "two-tailed lipid)")
    parser.add_argument("--start", type=int, default=0, help="(default 0) initial frame (match get_protection_state)")
    parser.add_argument("--stride", type=int, default=1, help="(default 1) frame stride (match get_protection_state)")
    args = parser.parse_args()

    if "{replica}" not in args.output_pattern or not args.output_pattern.endswith(".npy"):
        raise ValueError("output_pattern must contain {replica} and end in .npy")
    tail_names = [s.strip() for s in args.tail_bead_names.split(",") if s.strip()]

    with tb.open_file(args.hdx_topology, "r") as top:
        donor = top.root.input.potential.infer_H_O.donors.residue[:]
        n_residue = top.root.input.sequence.shape[0]

    def replica_system(view_path):
        source, pos, boxes = hybrid_frames(view_path, args.start, args.stride)
        with tb.open_file(source, "r") as t:
            particle_class = np.array([x.decode() for x in t.root.input.particle_class[:]])
            atom_names = np.array([x.decode() for x in t.root.input.atom_names[:]])
            residue_ids = t.root.input.residue_ids[:]
            atom_map = t.root.input.hybrid_bb_map.atom_indices[:]
        if atom_map.shape != (n_residue, 4):
            raise ValueError("hybrid_bb_map has shape {} but the HDX topology has {} residues".format(
                atom_map.shape, n_residue))
        is_lipid = particle_class == args.lipid_class
        is_tail, selected = tail_mask(is_lipid, atom_names, tail_names)
        tail_idx = np.where(is_tail)[0]
        if tail_idx.size == 0:
            raise ValueError(
                "no lipid tail beads found (class {}, {}); lipid bead names present: {}".format(
                    args.lipid_class, selected, ",".join(sorted(set(atom_names[is_lipid])))))
        return pos, boxes, tail_idx, residue_ids[tail_idx], atom_map[donor, 0], selected

    histogram = np.zeros(RDF_EDGES.size - 1)
    for view_path in args.hdx_views:
        pos, boxes, tail_idx, tail_molecule, _, selected = replica_system(view_path)
        histogram += tail_tail_histogram(pos, boxes, tail_idx, tail_molecule)
    print("tail beads: {} ({} beads)".format(selected, tail_idx.size))
    radius, peak = first_minimum(histogram)
    print("tail-tail g(r) pooled over {} replicas: first peak {:.2f} A, flat first minimum {:.2f} A".format(
        len(args.hdx_views), peak, radius))

    radius_path = args.output_pattern.replace("{replica}", "all")[:-4] + "_radius.txt"
    os.makedirs(os.path.dirname(os.path.abspath(radius_path)), exist_ok=True)
    with open(radius_path, "w") as handle:
        handle.write("# tail-tail g(r) pooled over {} replicas: first peak (A), contact radius (A)\n".format(
            len(args.hdx_views)))
        handle.write("{:.3f} {:.3f}\n".format(peak, radius))

    for replica, view_path in enumerate(args.hdx_views):
        pos, boxes, tail_idx, _, amide_idx, _ = replica_system(view_path)
        accessibility = np.ones((pos.shape[0], donor.size), dtype=np.float32)
        for f in range(pos.shape[0]):
            wrapped = np.mod(pos[f], boxes[f])
            tree = cKDTree(wrapped[tail_idx], boxsize=boxes[f])
            counts = tree.query_ball_point(wrapped[amide_idx], r=radius, return_length=True)
            accessibility[f, np.asarray(counts) > 0] = 0.0
        output_path = args.output_pattern.replace("{replica}", str(replica))
        np.save(output_path, accessibility)
        buried = (accessibility == 0.0).any(axis=0).sum()
        print("wrote {} : {} frames x {} donors; {} donors lipid-buried in >=1 frame, {:.1%} of amide-frames"
              .format(output_path, pos.shape[0], donor.size, int(buried), 1.0 - accessibility.mean()))


if __name__ == "__main__":
    main()
