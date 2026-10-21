import csv
import os

import numpy as np


def csv_list_env(name, default=None):
    raw = os.environ.get(name)
    if raw is None:
        return list(default) if default is not None else []
    return [item.strip() for item in raw.split(',') if item.strip()]


def float_list_env(name, default):
    raw = os.environ.get(name)
    if raw is None:
        return np.asarray(default, dtype=float)
    return np.asarray([float(item.strip()) for item in raw.split(',') if item.strip()], dtype=float)


def bool_env(name, default=False):
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {'1', 'true', 'yes', 'y', 'on'}


def build_paths(pdb_id, sim_id, work_dir='./'):
    output_dir = '{}/outputs'.format(work_dir)
    return {
        'work_dir': work_dir,
        'input_dir': '{}/inputs'.format(work_dir),
        'output_dir': output_dir,
        'result_dir': '{}/results'.format(work_dir),
        'run_dir': '{}/{}'.format(output_dir, sim_id),
        'pdb_dir': '{}/pdb'.format(work_dir),
        'pdb_id': pdb_id,
        'sim_id': sim_id,
    }


def _load_energy(path):
    energy = np.load(path)
    if getattr(energy, 'ndim', 1) > 1:
        return energy[:, 0]
    return energy


def load_replica_results(pdb_id, sim_id, result_dir, n_rep):
    Pot = []
    Rg = []
    Rmsd = []
    Hb = []
    Ts = []
    T = []
    PS = []

    for i in range(n_rep):
        prefix = '{}/{}_{}_{}'.format(result_dir, pdb_id, sim_id, i)
        energy = _load_energy('{}_Energy.npy'.format(prefix))
        rg = np.load('{}_Rg.npy'.format(prefix))
        hb = np.load('{}_Hbond.npy'.format(prefix))
        rmsd = np.load('{}_Rmsd.npy'.format(prefix))
        ps = np.load('{}_PS.npy'.format(prefix))
        t = np.load('{}_T.npy'.format(prefix))

        Pot.append(energy)
        Rg.append(rg)
        Hb.append(hb)
        Rmsd.append(rmsd)
        PS.append(ps)

        scalar_t = float(np.asarray(t).reshape(-1)[0])
        Ts.append(np.zeros(energy.size) + scalar_t)
        T.append(scalar_t)

    min_idx = min(len(arr) for arr in Rmsd)

    return {
        'Pot': np.asarray([arr[:min_idx] for arr in Pot]),
        'Rg': np.asarray([arr[:min_idx] for arr in Rg]),
        'Hb': np.asarray([arr[:min_idx] for arr in Hb]),
        'Rmsd': np.asarray([arr[:min_idx] for arr in Rmsd]),
        'Ts': np.asarray([arr[:min_idx] for arr in Ts]),
        'PS': np.asarray([arr[:min_idx] for arr in PS]),
        'T': np.asarray(T, dtype=float),
        'res': np.loadtxt('{}/{}.resid'.format(result_dir, pdb_id), dtype=int),
        'min_idx': min_idx,
    }


def resolve_start_frame(num_frames, default_fraction=0.20):
    raw = os.environ.get('start_frame')
    if raw is not None:
        return max(0, min(int(raw), max(num_frames - 1, 0)))
    return max(0, min(int(round(num_frames * default_fraction)), max(num_frames - 1, 0)))


def require_on_ladder(T, temps, name):
    """MBAR reweights within the sampled temperatures only; the ladder is logged in float32."""
    for value in np.atleast_1d(np.asarray(temps, dtype=float)):
        if not np.float32(np.min(T)) <= np.float32(value) <= np.float32(np.max(T)):
            raise SystemExit('{} = {:.4f} lies outside the replica ladder {:.4f}-{:.4f}'.format(
                name, value, np.min(T), np.max(T)))


def load_pdb_numbering(input_dir, pdb_id):
    """(chain, residue number) in the PDB of every Upside residue, from the .chi file PDB_to_initial_structure.py
    writes. For a hybrid view martini_hdx_project.py has checked it equals the hybrid's hybrid_bb_map/bb_resseq."""
    path = '{}/{}.chi'.format(input_dir, pdb_id)
    if not os.path.exists(path):
        raise SystemExit('missing {}, the PDB numbering PDB_to_initial_structure.py writes'.format(path))
    with open(path) as handle:
        rows = [line.split() for line in handle.read().splitlines()[1:] if line.strip()]
    # residue restype chain resnum chi1 chi2, the chain column empty for a PDB without chain ids
    return (np.array([row[2] if len(row) == 6 else '' for row in rows]),
            np.array([int(row[-3]) for row in rows]))


def residue_labels(chains, numbers):
    """PDB residue labels, prefixed with the chain when there is more than one."""
    if np.unique(chains).size > 1:
        return np.array(['{}{}'.format(chain, number) for chain, number in zip(chains, numbers)])
    return np.array([str(number) for number in numbers])


def label_residue_axis(ax, labels):
    """Label an axis drawn against the Upside residue index with the PDB numbering of each residue."""
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    def label(x, position):
        index = int(round(x))
        return labels[index] if abs(x - index) < 1e-6 and 0 <= index < len(labels) else ''

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_major_formatter(FuncFormatter(label))
    ax.set_xlabel('Residue (PDB numbering)')


def save_csv_rows(path, rows):
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle)
        for row in rows:
            writer.writerow(row)


def load_csv_rows(path):
    with open(path, 'r', encoding='utf-8-sig') as handle:
        return list(csv.reader(handle))


def load_optional_numeric_csv(path, skip_header=0, delimiter=','):
    if not os.path.exists(path):
        return None
    return np.genfromtxt(path, skip_header=skip_header, delimiter=delimiter)
