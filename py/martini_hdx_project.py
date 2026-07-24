#!/usr/bin/env python

import argparse
import os
import shutil

import numpy as np
import tables as tb

import upside_engine as ue


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


def write_projected_potential(engine, source_pos, atom_indices, destination_group, filters):
    """Re-score the protein-only Upside potential for each projected frame.

    The hybrid /output/potential is the total (protein + lipid + ion) energy, which is
    dominated by the lipid bath and nearly temperature-independent. HDX temperature
    reweighting needs the protein subsystem energy, so it is recomputed here with the
    protein-only HDX engine on the projected backbone atoms.
    """
    n_frame = source_pos.shape[0]
    output = destination_group._v_file.create_carray(
        destination_group,
        'potential',
        atom=tb.Float32Atom(),
        shape=(n_frame, 1),
        filters=filters,
    )
    chunk = max(1, min(64, n_frame))
    for start in range(0, n_frame, chunk):
        stop = min(start + chunk, n_frame)
        frames = source_pos[start:stop, 0][:, atom_indices]
        for offset in range(stop - start):
            output[start + offset, 0] = engine.energy(frames[offset])


def validate_mapping(source, destination):
    if not hasattr(source.root.input, 'hybrid_bb_map'):
        raise ValueError('hybrid trajectory is missing /input/hybrid_bb_map')
    mapping = source.root.input.hybrid_bb_map
    if not hasattr(mapping, 'atom_indices'):
        raise ValueError('hybrid trajectory is missing /input/hybrid_bb_map/atom_indices')

    atom_map = np.asarray(mapping.atom_indices[:], dtype=np.int64)
    if atom_map.ndim != 2 or atom_map.shape[1] != 4:
        raise ValueError('hybrid backbone atom_indices must have shape (n_residue, 4)')

    n_residue = destination.root.input.sequence.shape[0]
    if atom_map.shape[0] != n_residue:
        raise ValueError(
            'hybrid mapping has {} residues but HDX topology has {}'.format(
                atom_map.shape[0], n_residue))

    n_source_atom = source.root.input.pos.shape[0]
    backbone = atom_map[:, :3].reshape(-1)
    if np.any(backbone < 0) or np.any(backbone >= n_source_atom):
        raise ValueError('hybrid N/CA/C mapping contains an out-of-range atom index')
    if np.unique(backbone).size != backbone.size:
        raise ValueError('hybrid N/CA/C mapping contains duplicate atom indices')
    return backbone


def project(hdx_topology, hybrid_trajectory, output_path, overwrite=False):
    for path, description in ((hdx_topology, 'HDX topology'),
                              (hybrid_trajectory, 'hybrid trajectory')):
        if not os.path.isfile(path):
            raise FileNotFoundError('cannot read {}: {}'.format(description, path))
    if os.path.abspath(output_path) in (os.path.abspath(hdx_topology),
                                        os.path.abspath(hybrid_trajectory)):
        raise ValueError('output path must differ from both input paths')
    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError('output already exists: {}'.format(output_path))

    output_directory = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(output_directory, exist_ok=True)
    temporary_path = output_path + '.tmp'
    if os.path.exists(temporary_path):
        os.remove(temporary_path)

    engine = ue.Upside(hdx_topology)

    shutil.copyfile(hdx_topology, temporary_path)
    filters = tb.Filters(complevel=5, complib='zlib', shuffle=True)
    copied_arrays = ('hbond', 'temperature', 'time')

    try:
        with tb.open_file(hybrid_trajectory, 'r') as source, tb.open_file(temporary_path, 'a') as destination:
            atom_indices = validate_mapping(source, destination)
            group_names = output_group_names(source)
            if not group_names:
                raise ValueError('hybrid trajectory has no output groups')

            for name in output_group_names(destination):
                destination.remove_node('/' + name, recursive=True)

            for name in group_names:
                source_group = source.get_node('/' + name)
                missing = [array for array in copied_arrays if not hasattr(source_group, array)]
                if missing:
                    raise ValueError('{} is missing {}'.format(source_group._v_pathname, ', '.join(missing)))
                if source_group.pos.ndim != 4 or source_group.pos.shape[1] != 1:
                    raise ValueError('{} pos must have shape (n_frame, 1, n_atom, 3)'.format(
                        source_group._v_pathname))

                destination_group = destination.create_group('/', name)
                copy_attributes(source_group, destination_group)
                copy_projected_positions(source_group.pos, destination_group, atom_indices, filters)
                write_projected_potential(engine, source_group.pos, atom_indices, destination_group, filters)
                for array_name in copied_arrays:
                    copy_array(getattr(source_group, array_name), destination_group, array_name, filters)
                if hasattr(source_group, 'replica_index'):
                    copy_array(source_group.replica_index, destination_group, 'replica_index', filters)

            destination.root._v_attrs.hdx_trajectory_kind = 'martini_hybrid_projection'
            destination.root._v_attrs.hdx_hybrid_source = os.path.abspath(hybrid_trajectory)
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
    parser.add_argument('hybrid_trajectory', help='full hybrid trajectory .up file')
    parser.add_argument('output_h5', help='projected protein-only trajectory-view .up file')
    parser.add_argument('--overwrite', action='store_true', help='replace an existing output file')
    args = parser.parse_args()

    project(args.hdx_topology, args.hybrid_trajectory, args.output_h5, overwrite=args.overwrite)
    print('HDX trajectory view written to {}'.format(args.output_h5))


if __name__ == '__main__':
    main()
