import os
import sys

import mdtraj as md
import mdtraj_upside as mu
import numpy as np
import tables as tb

from helpers.advanced_analysis_utils import bool_env, build_paths


pdb_id = os.environ.get('pdb_id', 'Pab1_RRM1')  # CHECKME
sim_id = os.environ.get('sim_id', 'REMD')  # CHECKME
n_rep = int(os.environ.get('n_rep', '48'))  # CHECKME
rmsd_limit = float(os.environ.get('rmsd_limit', '10'))  # CHECKME
stop_once_exceed = bool_env('stop_once_exceed', True)  # CHECKME

paths = build_paths(pdb_id, sim_id)


def _output_groups(tables_file):
    i = 0
    while 'output_previous_{}'.format(i) in tables_file.root:
        yield tables_file.get_node('/output_previous_{}'.format(i))
        i += 1
    if 'output' in tables_file.root:
        yield tables_file.get_node('/output')


def write_h5(source_h5_path, destination_h5_path, selected_indices):
    with tb.open_file(source_h5_path, 'r') as src_file:
        with tb.open_file(destination_h5_path, 'w') as dest_file:
            src_file.copy_node('/input', newparent=dest_file.root, recursive=True)
            if '/output' not in dest_file:
                dest_file.create_group('/', 'output')

            lengths = []
            data = []
            node_names = []
            for output_group in _output_groups(src_file):
                lengths.append(len(output_group.pos))
                for node in src_file.list_nodes(output_group):
                    data.append(node.read())
                    node_names.append(node.name)

            total_length = sum(lengths)
            node_count = int(len(node_names) / (src_file.root._v_nchildren - 1))
            if len(node_names[:node_count]) != len(set(node_names[:node_count])):
                raise SystemExit('Output nodes do not match across groups in {}'.format(source_h5_path))
            node_names = node_names[:node_count]

            combined_data = []
            for i in range(src_file.root._v_nchildren - 1):
                if i == 0:
                    combined_data = data[:node_count]
                    continue

                node_data = data[node_count * i:node_count * (i + 1)]
                previous_data = combined_data
                combined_data = []
                for old_data, new_data in zip(previous_data, node_data):
                    combined_data.append(np.concatenate([old_data, new_data]))

            for combined, node_name in zip(combined_data, node_names):
                if combined.shape[0] == total_length:
                    filtered = combined[selected_indices]
                else:
                    filtered = combined
                dest_file.create_array('/output', node_name, obj=filtered)


target = md.load('{}/{}.pdb'.format(paths['pdb_dir'], pdb_id))

for i in range(n_rep):
    traj_path = '{}/{}.run.{}.up'.format(paths['run_dir'], pdb_id, i)
    selected_path = '{}/{}.selected.{}.up'.format(paths['run_dir'], pdb_id, i)

    traj = mu.load_upside_traj(traj_path)
    sele_ca = traj.topology.select('name CA')
    sele_ca_target = target.topology.select('name CA')
    rmsd_values = md.rmsd(traj.atom_slice(sele_ca), target.atom_slice(sele_ca_target)) * 10

    below_limit_indices = np.where(rmsd_values <= rmsd_limit)[0]
    if stop_once_exceed and np.any(rmsd_values > rmsd_limit):
        first_exceed_index = np.argmax(rmsd_values > rmsd_limit)
        below_limit_indices = below_limit_indices[below_limit_indices < first_exceed_index]

    write_h5(traj_path, selected_path, below_limit_indices)
