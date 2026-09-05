#!/usr/bin/env python3
"""Add Upside's intra-protein sidechain coverage nodes to a hybrid MARTINI .up.

A hybrid config gives the rotamer solver one 1-body field, the MARTINI SC-env table, and so carries
neither of the two coverage terms a standard Upside config supplies: `hbond_coverage`, which lets a
sidechain compete with a backbone hydrogen bond, and `hbond_coverage_hydrophobe`, which covers the
oriented backbone atoms themselves. Both are protein-internal and MARTINI has no hydrogen bond, so
restoring them neither replaces nor double-counts a protein-environment term. `RotamerSidechain` sums a
variable-length list of 1-body nodes, so no C++ change is needed; the two nodes are appended to the
rotamer argument list, which keeps `martini_sc_table_1body` in place.

The node schemas and the force-field tables are those of upside_config.write_rotamer_backbone. The
sidechain bead list, the hydrogen-bond donor and acceptor lists, the residue sequence and the residue
count all come out of the file being edited, so nothing here knows which protein it is looking at.

Idempotent: the nodes are rewritten from scratch and the rotamer argument list is rebuilt rather than
appended to. No MARTINI node is read or written.
"""
import argparse

import h5py
import numpy as np

# upside_config.backbone_group4, the grouping the coverage tables were trained against. Copied rather
# than imported because upside_config needs pytables, which the cluster runtime does not have.
backbone_group4 = { 'ALA': 0,  'ARG': 0,  'ASN': 0,  'ASP': 0,
                    'CYS': 0,  'GLN': 0,  'GLU': 0,  'GLY': 1,
                    'HIS': 0,  'ILE': 0,  'LEU': 0,  'LYS': 0,
                    'MET': 0,  'PHE': 0,  'PRO': 2,  'SER': 0,
                    'THR': 0,  'TRP': 0,  'TYR': 0,  'VAL': 0  }

POTENTIAL = "/input/potential"
BB_NODE = "placement_fixed_point_vector_scalar"
COVERAGE = "hbond_coverage"
HYDROPHOBE = "hbond_coverage_hydrophobe"


def backbone_types(sequence):
    n_group = len(set(backbone_group4.values()))
    types = np.zeros(len(sequence), dtype=np.int32)
    for resid, res in enumerate(sequence):
        if res not in backbone_group4:
            raise SystemExit("residue %i is %s, which has no backbone group" % (resid, res))
        if resid < len(sequence)-1 and sequence[resid+1] in ("PRO", "CPR"):
            types[resid] = n_group
        else:
            types[resid] = backbone_group4[res]
    return types


def recreate_group(parent, name):
    if name in parent:
        del parent[name]
    return parent.create_group(name)


def bstrings(names):
    return np.asarray([np.bytes_(x) for x in names])


def inject(path, sidechain_lib):
    with h5py.File(sidechain_lib, "r") as lib:
        bead_num = dict((k.decode(), i) for i, k in enumerate(lib["bead_order"][:]))
        coverage_interaction = lib["coverage_interaction"][:]
        hydrophobe_placement = lib["hydrophobe_placement"][:]
        hydrophobe_interaction = lib["hydrophobe_interaction"][:]

    with h5py.File(path, "r+") as h5:
        pot = h5[POTENTIAL]
        for name in ("affine_alignment", "infer_H_O", "protein_hbond", "rotamer"):
            if name not in pot:
                raise SystemExit("%s: has no %s node, so it is not a hybrid config" % (path, name))
        if "sequence" not in h5["/input"]:
            raise SystemExit("%s: has no /input/sequence to take backbone groups from" % path)

        # the rotamer solver's first argument is the sidechain placement node, and the coverage nodes
        # must score the same beads it does
        rotamer = pot["rotamer"]
        rotamer_args = [x.decode() for x in rotamer.attrs["arguments"]]
        sc_name = rotamer_args[0]

        sequence = [x.decode() for x in h5["/input/sequence"][:]]
        n_res = pot["affine_alignment/ref_geom"].shape[0]
        if len(sequence) != n_res:
            raise SystemExit("%s: /input/sequence has %i residues but affine_alignment has %i"
                             % (path, len(sequence), n_res))
        bb_type = backbone_types(sequence)

        d_residues = pot["infer_H_O/donors/residue"][:]
        a_residues = pot["infer_H_O/acceptors/residue"][:]
        hb_index = np.arange(len(d_residues) + len(a_residues))

        sc_node = pot[sc_name]
        sc_resnum = sc_node["affine_residue"][:]
        sc_type = np.array([bead_num[s.decode()] for s in sc_node["beadtype_seq"][:]])
        sc_index = np.arange(len(sc_type))

        # sc-hbond interaction
        cgrp = recreate_group(pot, COVERAGE)
        cgrp.attrs["arguments"] = bstrings(["protein_hbond", sc_name])
        cgrp.create_dataset("interaction_param", data=coverage_interaction, dtype=np.float32)
        # group1 is the HBond partners, donor is 0 and acceptor is 1
        cgrp.create_dataset("index1", data=hb_index, dtype=np.int32)
        cgrp.create_dataset("type1", data=np.concatenate([2*bb_type[d_residues],
                                                          2*bb_type[a_residues] + 1]), dtype=np.int32)
        cgrp.create_dataset("id1", data=np.concatenate([d_residues, a_residues]), dtype=np.int32)
        # group2 is the sc
        cgrp.create_dataset("index2", data=sc_index, dtype=np.int32)
        cgrp.create_dataset("type2", data=sc_type, dtype=np.int32)
        cgrp.create_dataset("id2", data=sc_resnum, dtype=np.int32)

        # the oriented backbone atoms
        bb_index = np.arange(3*n_res)
        bb_resnum = bb_index//3
        grp = recreate_group(pot, BB_NODE)
        grp.attrs["arguments"] = bstrings(["affine_alignment"])
        grp.create_dataset("affine_residue", data=bb_resnum, dtype=np.int32)
        grp.create_dataset("layer_index", data=bb_index % 3, dtype=np.int32)
        grp.create_dataset("placement_data", data=hydrophobe_placement, dtype=np.float32)

        # sc-backbone interaction
        cgrp = recreate_group(pot, HYDROPHOBE)
        cgrp.attrs["arguments"] = bstrings([BB_NODE, sc_name])
        cgrp.create_dataset("interaction_param", data=hydrophobe_interaction, dtype=np.float32)
        # group1 is the backbone partners
        cgrp.create_dataset("index1", data=bb_index, dtype=np.int32)
        cgrp.create_dataset("type1", data=3*bb_type[bb_resnum] + bb_index % 3, dtype=np.int32)
        cgrp.create_dataset("id1", data=bb_resnum, dtype=np.int32)
        # group2 is the sidechains
        cgrp.create_dataset("index2", data=sc_index, dtype=np.int32)
        cgrp.create_dataset("type2", data=sc_type, dtype=np.int32)
        cgrp.create_dataset("id2", data=sc_resnum, dtype=np.int32)

        kept = [x for x in rotamer_args if x not in (COVERAGE, HYDROPHOBE)]
        rotamer.attrs["arguments"] = bstrings(kept + [COVERAGE, HYDROPHOBE])

        return "%i hbond partners, %i backbone atoms, %i sidechain beads" % (
                len(hb_index), len(bb_index), len(sc_index))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sidechain", required=True, help="sidechain.h5 to take the coverage tables from")
    parser.add_argument("files", nargs="+", help=".up files to edit in place")
    args = parser.parse_args()

    for path in args.files:
        print("%-58s %s" % (path.split("/")[-1], inject(path, args.sidechain)), flush=True)


if __name__ == "__main__":
    main()
