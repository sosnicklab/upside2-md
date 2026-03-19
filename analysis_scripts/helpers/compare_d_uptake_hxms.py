import csv
import os
import sys

import matplotlib
if 'MPLBACKEND' not in os.environ:
    matplotlib.use('Agg')
import matplotlib.backends.backend_pdf as pdf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp

from helpers.advanced_analysis_utils import build_paths, csv_list_env


sim_id = os.environ.get('sim_id', 'REMD')  # CHECKME
HXMS_method = os.environ.get('HXMS_method', 'stretch_exp')  # CHECKME
pdb_ids = csv_list_env('pdb_ids', ['Scer_Pab1', 'Skud_Pab1'])  # CHECKME
protein_states = csv_list_env('protein_states', ['monomer', 'monomer'])  # CHECKME

if len(pdb_ids) != 2:
    raise SystemExit('Set pdb_ids to exactly two comma-separated protein IDs for 5.compare_D_uptake_HXMS.py.')

if len(protein_states) == 1:
    protein_states = protein_states * 2
if len(protein_states) != 2:
    raise SystemExit('Set protein_states to one or two comma-separated entries.')

paths = build_paths(pdb_ids[0], sim_id)
output_dir = paths['output_dir']
fig_dir = paths['result_dir']

peps_HXMS = []
time_peps_HXMS = []
d_norm_peps_HXMS = []
pep_starts = []
pep_ends = []
pep_ids = []

for pdb_id, protein_state in zip(pdb_ids, protein_states):
    peps_HXMS.append(np.load('{}/{}/{}_{}_peps_{}.npy'.format(output_dir, sim_id, pdb_id, HXMS_method, protein_state)))
    time_peps_HXMS.append(np.load('{}/{}/{}_{}_time_peps_{}.npy'.format(output_dir, sim_id, pdb_id, HXMS_method, protein_state)))
    d_norm_peps_HXMS.append(np.load('{}/{}/{}_{}_d_norm_peps_{}.npy'.format(output_dir, sim_id, pdb_id, HXMS_method, protein_state)))

    with open('{}/{}/{}_pep_ids.csv'.format(output_dir, sim_id, pdb_id), 'r', encoding='utf-8-sig') as handle:
        rows = list(csv.reader(handle))
    pep_starts.append(np.array(rows[0]))
    pep_ends.append(np.array(rows[1]))
    pep_ids.append(np.array(rows[2]))

peps_all = np.concatenate(peps_HXMS)
seen = set()
peps_uniq = [pep for pep in peps_all if pep not in seen and not seen.add(pep)]

hxms_dfs = pd.DataFrame({'peps': peps_uniq})
for pdb_id, peps, starts, ends, times, d_norm in zip(pdb_ids, peps_HXMS, pep_starts, pep_ends, time_peps_HXMS, d_norm_peps_HXMS):
    frame = pd.DataFrame(
        {
            'peps': peps,
            'peps_start_{}'.format(pdb_id): starts,
            'peps_end_{}'.format(pdb_id): ends,
            'time_peps_{}'.format(pdb_id): [times for _ in range(len(peps))],
            'd_norms_HXMS_{}'.format(pdb_id): d_norm.tolist(),
        }
    )
    hxms_dfs = pd.merge(hxms_dfs, frame, on='peps')

pdb_1, pdb_2 = pdb_ids
hxms_pep_start = hxms_dfs['peps_start_{}'.format(pdb_1)]
hxms_pep_end = hxms_dfs['peps_end_{}'.format(pdb_1)]
d_norms_hxms_pdb_1 = hxms_dfs['d_norms_HXMS_{}'.format(pdb_1)]
d_norms_hxms_pdb_2 = hxms_dfs['d_norms_HXMS_{}'.format(pdb_2)]

pdf_d = pdf.PdfPages('{}/D_uptake_{}_{}.pdf'.format(fig_dir, pdb_1, pdb_2))
for pep, p_start, p_end, d_norm_hxms_1, d_norm_hxms_2 in zip(
    hxms_dfs.peps,
    hxms_pep_start,
    hxms_pep_end,
    d_norms_hxms_pdb_1,
    d_norms_hxms_pdb_2,
):
    d_norm_hxms_1 = np.asarray(d_norm_hxms_1, dtype=float)
    d_norm_hxms_2 = np.asarray(d_norm_hxms_2, dtype=float)
    if np.all(np.isnan(d_norm_hxms_1)) or np.all(np.isnan(d_norm_hxms_2)):
        continue

    d_norm_1 = (d_norm_hxms_1 - d_norm_hxms_1[0]) / (np.nanmax(d_norm_hxms_1) - d_norm_hxms_1[0]) * 100
    d_norm_2 = (d_norm_hxms_2 - d_norm_hxms_2[0]) / (np.nanmax(d_norm_hxms_2) - d_norm_hxms_2[0]) * 100

    interpolator = sp.interpolate.interp1d(np.arange(d_norm_2.size), d_norm_2)
    d_norm_2_plot = interpolator(np.linspace(0, d_norm_2.size - 1, d_norm_1.size))
    d_norm_stats = sp.stats.linregress(d_norm_1, d_norm_2_plot)
    r_square = d_norm_stats.rvalue ** 2
    if np.isnan(r_square):
        continue

    fig = plt.figure()
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.scatter(d_norm_1, d_norm_2_plot)
    plt.xlabel('{} %D Uptake'.format(pdb_1))
    plt.ylabel('{} %D Uptake'.format(pdb_2))
    plt.title('{}-{} {}'.format(p_start, p_end, pep))
    diagonal = np.linspace(0, 100, 100)
    plt.plot(diagonal, diagonal, color='k', lw=2.5)
    plt.text(5, 95, '$R^2$ = {}'.format(round(r_square, 8)))
    pdf_d.savefig(fig)
    plt.close()

pdf_d.close()
