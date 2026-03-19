import os
from math import ceil

import matplotlib
if 'MPLBACKEND' not in os.environ:
    matplotlib.use('Agg')
import matplotlib.backends.backend_pdf as pdf
import matplotlib.pyplot as plt
import mdtraj_upside as mu
import numpy as np
from scipy.stats import gaussian_kde

try:
    from sklearn.decomposition import PCA as SklearnPCA
except ImportError:
    SklearnPCA = None

from helpers.advanced_analysis_utils import build_paths, load_replica_results, resolve_start_frame


def reduce_cartesian_to_two_components(xyz):
    flattened = xyz.reshape(xyz.shape[0], xyz.shape[1] * 3)
    if SklearnPCA is not None:
        return SklearnPCA(n_components=2).fit_transform(flattened)

    centered = flattened - np.mean(flattened, axis=0, keepdims=True)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    return centered @ vt[:2].T


pdb_id = os.environ.get('pdb_id', 'Pab1_RRM1')  # CHECKME
sim_id = os.environ.get('sim_id', 'REMD')  # CHECKME
n_rep = int(os.environ.get('n_rep', '48'))  # CHECKME

paths = build_paths(pdb_id, sim_id)
results = load_replica_results(pdb_id, sim_id, paths['result_dir'], n_rep)

Rmsd = results['Rmsd']
Hb = results['Hb']
min_idx = results['min_idx']
start_frame = resolve_start_frame(min_idx)

ncols = 4
nrows = -(-n_rep // ncols)
pdf_traj = pdf.PdfPages('{}/{}_traj_RMSD_Hbond_PCA.pdf'.format(paths['result_dir'], pdb_id))

fig = plt.figure(figsize=(4 * ncols, 3 * nrows), facecolor='w')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
for i in range(n_rep):
    plt.subplot(int(ceil(n_rep / ncols)), ncols, i + 1)
    plt.plot(Rmsd[i, start_frame:])
    plt.title('Trajectory {}'.format(i))
    plt.xlim([0, max(min_idx - start_frame, 1)])
plt.tight_layout()
pdf_traj.savefig(fig, dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(4 * ncols, 3 * nrows), facecolor='w')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
for i in range(n_rep):
    x = Hb[i, start_frame:]
    y = Rmsd[i, start_frame:]
    plt.subplot(int(ceil(n_rep / ncols)), ncols, i + 1)
    try:
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        plt.scatter(x, y, c=z, cmap='jet', s=8)
    except Exception:
        plt.scatter(x, y, s=8)
    plt.title('Trajectory {}'.format(i))
    plt.xlim([0, max(100, np.nanmax(x) + 5)])
    plt.ylim([0, max(10, np.nanmax(y) + 5)])
plt.tight_layout()
pdf_traj.savefig(fig, dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(4 * ncols, 3 * nrows), facecolor='w')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
for i in range(n_rep):
    traj = mu.load_upside_traj('{}/{}.run.{}.up'.format(paths['run_dir'], pdb_id, i))
    traj.superpose(traj, 0)
    reduced_cartesian = reduce_cartesian_to_two_components(traj.xyz)

    plt.subplot(int(ceil(n_rep / ncols)), ncols, i + 1)
    plt.scatter(reduced_cartesian[:, 0], reduced_cartesian[:, 1], c=traj.time, cmap='jet', s=8)
    plt.title('Trajectory {}'.format(i))
plt.tight_layout()
pdf_traj.savefig(fig, dpi=300, bbox_inches='tight')
plt.close()

pdf_traj.close()
