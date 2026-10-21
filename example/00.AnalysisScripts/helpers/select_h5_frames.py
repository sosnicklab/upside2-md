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
            dest_file.create_group('/', 'output')

            # Join each dataset by name across the restart segments. The first frame of every segment after the
            # first repeats the last frame before the restart; mdtraj_upside.load_upside_traj drops it, and
            # selected_indices count frames that way. Datasets that are not per-frame are taken once.
            per_frame = {}
            once = {}
            groups = list(_output_groups(src_file))
            for g_no, output_group in enumerate(groups):
                n_frame = len(output_group.pos)
                for node in src_file.list_nodes(output_group):
                    if node.shape[0] == n_frame:
                        per_frame.setdefault(node.name, []).append(node.read()[1 if g_no else 0:])
                    else:
                        once.setdefault(node.name, node.read())

            for node_name, parts in per_frame.items():
                if len(parts) != len(groups):
                    raise SystemExit('{} is not logged in every output group of {}'.format(node_name, source_h5_path))
                dest_file.create_array('/output', node_name, obj=np.concatenate(parts)[selected_indices])
            for node_name, data in once.items():
                dest_file.create_array('/output', node_name, obj=data)


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
