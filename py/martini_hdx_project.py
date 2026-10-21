#!/usr/bin/env python

import argparse
import os
import shutil

import h5py
import numpy as np
import tables as tb


def output_group_names(h5):
    names = []
    index = 0
    while hasattr(h5.root, 'output_previous_{}'.format(index)):
        names.append('output_previous_{}'.format(index))
        index += 1
    if hasattr(h5.root, 'output'):
        names.append('output')
    return names


def copy_attributes(source, destination):
    for name in source._v_attrs._f_list():
        destination._v_attrs[name] = source._v_attrs[name]


def copy_array(source, destination_group, name, filters):
    atom = tb.Atom.from_dtype(source.dtype)
    output = destination_group._v_file.create_carray(
        destination_group,
        name,
        atom=atom,
        shape=source.shape,
        filters=filters,
    )
    chunk = max(1, min(256, source.shape[0]))
    for start in range(0, source.shape[0], chunk):
        stop = min(start + chunk, source.shape[0])
        output[start:stop] = source[start:stop]
    copy_attributes(source, output)


def copy_projected_positions(source, destination_group, atom_indices, filters):
    n_frame = source.shape[0]
    n_atom = atom_indices.size
    output = destination_group._v_file.create_carray(
        destination_group,
        'pos',
        atom=tb.Float32Atom(),
        shape=(n_frame, 1, n_atom, 3),
        filters=filters,
    )
    chunk = max(1, min(64, n_frame))
    for start in range(0, n_frame, chunk):
        stop = min(start + chunk, n_frame)
        positions = source[start:stop, 0]
        output[start:stop, 0] = positions[:, atom_indices]
    copy_attributes(source, output)


def read_pdb_numbering(path):
    """(chain, residue number) of every residue, from the .chi file PDB_to_initial_structure.py writes."""
    chains = []
    numbers = []
    with open(path) as handle:
        next(handle)
        for line in handle:
            # residue restype chain resnum chi1 chi2, the chain column empty for a PDB without chain ids
            fields = line.split()
            chains.append(fields[2] if len(fields) == 6 else '')
            numbers.append(int(fields[-3]))
    return np.array(chains), np.array(numbers)


def validate_mapping(source_path, source, destination, pdb_numbering):
    """N/CA/C indices of the hybrid residues, refused unless the hybrid protein is the HDX topology's protein."""
    if not hasattr(source.root.input, 'hybrid_bb_map'):
        raise ValueError('hybrid trajectory is missing /input/hybrid_bb_map')
    mapping = source.root.input.hybrid_bb_map
    atom_map = np.asarray(mapping.atom_indices[:], dtype=np.int64)
    n_residue = destination.root.input.sequence.shape[0]
    if atom_map.shape != (n_residue, 4):
        raise ValueError('hybrid_bb_map/atom_indices has shape {} but the HDX topology has {} residues'.format(
            atom_map.shape, n_residue))
    if not np.array_equal(source.root.input.sequence[:], destination.root.input.sequence[:]):
        raise ValueError('hybrid sequence differs from the HDX topology sequence')

    n_source_atom = source.root.input.pos.shape[0]
    backbone = atom_map[:, :3]
    if np.any(backbone < 0) or np.any(backbone >= n_source_atom):
        raise ValueError('hybrid N/CA/C mapping contains an out-of-range atom index')
    roles = np.array([x.decode() for x in source.root.input.atom_roles[:]])
    if not np.all(roles[backbone] == np.array(['N', 'CA', 'C'])):
        raise ValueError('hybrid_bb_map/atom_indices do not point at N, CA, C for every residue')

    # The hybrid's own H-bond log is copied into the view, so its donors and acceptors must be the HDX topology's.
    for group in ('donors', 'acceptors'):
        hybrid = getattr(source.root.input.potential.infer_H_O, group).residue[:]
        topology = getattr(destination.root.input.potential.infer_H_O, group).residue[:]
        if not np.array_equal(hybrid, topology):
            raise ValueError('hybrid infer_H_O {} differ from the HDX topology ({} against {})'.format(
                group, hybrid.size, topology.size))

    chains, numbers = pdb_numbering
    with h5py.File(source_path, 'r') as h5:
        bb_resseq = h5['input/hybrid_bb_map/bb_resseq'][:]
        bb_chain = np.array([x.decode() if isinstance(x, bytes) else str(x)
                             for x in h5['input/hybrid_bb_map/bb_chain_id'][:]])
    named = chains != ''
    if not (np.array_equal(bb_resseq, numbers) and np.array_equal(bb_chain[named], chains[named])):
        raise ValueError('hybrid_bb_map residue numbering differs from the PDB numbering of the HDX topology')
    return backbone.reshape(-1)


def project(hdx_topology, hybrid_segments, output_path, pdb_numbering_path, overwrite=False):
    """Write the protein-only view of one replica: its hybrid segments, in order, as consecutive output groups.

    Every group after the first starts at the frame that ended the group before it, the convention the HDX
    readers rely on when they drop that first frame, and a segment that does not is refused."""
    for path, description in [(hdx_topology, 'HDX topology'), (pdb_numbering_path, 'PDB numbering')] + [
            (path, 'hybrid trajectory') for path in hybrid_segments]:
        if not os.path.isfile(path):
            raise FileNotFoundError('cannot read {}: {}'.format(description, path))
    if os.path.abspath(output_path) in [os.path.abspath(p) for p in [hdx_topology] + list(hybrid_segments)]:
        raise ValueError('output path must differ from every input path')
    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError('output already exists: {}'.format(output_path))

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    temporary_path = output_path + '.tmp'
    if os.path.exists(temporary_path):
        os.remove(temporary_path)

    shutil.copyfile(hdx_topology, temporary_path)
    filters = tb.Filters(complevel=5, complib='zlib', shuffle=True)
    # The potential is the coupled protein + environment energy, copied unchanged: temperature
    # reweighting across the ladder needs the whole system's energy (findings 5.1, 3.5).
    copied_arrays = ('hbond', 'potential', 'temperature', 'time')
    pdb_numbering = read_pdb_numbering(pdb_numbering_path)

    try:
        with tb.open_file(temporary_path, 'a') as destination:
            for name in output_group_names(destination):
                destination.remove_node('/' + name, recursive=True)

            n_group = 0
            previous_last_frame = None
            for segment in hybrid_segments:
                with tb.open_file(segment, 'r') as source:
                    atom_indices = validate_mapping(segment, source, destination, pdb_numbering)
                    names = output_group_names(source)
                    if not names:
                        raise ValueError('{} has no output groups'.format(segment))
                    for name in names:
                        source_group = source.get_node('/' + name)
                        missing = [array for array in copied_arrays if not hasattr(source_group, array)]
                        if missing:
                            raise ValueError('{}:{} is missing {}'.format(segment, name, ', '.join(missing)))
                        if previous_last_frame is not None and not np.array_equal(
                                source_group.pos[0, 0], previous_last_frame):
                            raise ValueError('{}:{} does not start at the last frame of the group before it'.format(
                                segment, name))
                        previous_last_frame = source_group.pos[-1, 0]

                        destination_group = destination.create_group('/', 'output_previous_{}'.format(n_group))
                        n_group += 1
                        copy_attributes(source_group, destination_group)
                        copy_projected_positions(source_group.pos, destination_group, atom_indices, filters)
                        for array_name in copied_arrays:
                            copy_array(getattr(source_group, array_name), destination_group, array_name, filters)
                        if hasattr(source_group, 'replica_index'):
                            copy_array(source_group.replica_index, destination_group, 'replica_index', filters)
            destination.rename_node('/output_previous_{}'.format(n_group - 1), 'output')

            destination.root._v_attrs.hdx_trajectory_kind = 'martini_hybrid_projection'
            destination.root._v_attrs.hdx_hybrid_sources = np.array(
                [os.path.abspath(path).encode() for path in hybrid_segments])
            destination.root._v_attrs.hdx_topology_source = os.path.abspath(hdx_topology)
            destination.root._v_attrs.hdx_projected_atom_count = int(atom_indices.size)
            destination.flush()
        os.replace(temporary_path, output_path)
    except Exception:
        if os.path.exists(temporary_path):
            os.remove(temporary_path)
        raise


def main():
    parser = argparse.ArgumentParser(
        description='Project a full dry-MARTINI hybrid trajectory onto the standard Upside HDX trajectory contract.')
    parser.add_argument('hdx_topology', help='protein-only *-HDX.up configuration')
    parser.add_argument('output_h5', help='projected protein-only trajectory-view .up file')
    parser.add_argument('hybrid_segments', nargs='+',
                        help='the hybrid trajectory of one replica: one .up file, or the consecutive segments of one '
                             'run in order (a production stage_7.0, stage_7.1, ...)')
    parser.add_argument('--pdb-numbering', required=True,
                        help='.chi file of the PDB the HDX topology was built from; the hybrid residue numbering '
                             '(hybrid_bb_map/bb_resseq) must equal it')
    parser.add_argument('--overwrite', action='store_true', help='replace an existing output file')
    args = parser.parse_args()

    project(args.hdx_topology, args.hybrid_segments, args.output_h5, args.pdb_numbering, overwrite=args.overwrite)
    print('HDX trajectory view written to {}'.format(args.output_h5))


if __name__ == '__main__':
    main()
