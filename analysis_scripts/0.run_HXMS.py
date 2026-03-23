import os

import numpy as np
import scipy as sp
import matplotlib
if 'MPLBACKEND' not in os.environ:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib.font_manager import FontProperties
import pymbar  # for MBAR analysis
from pymbar import timeseries  # for timeseries analysis
import csv
import pandas as pd
from collections import OrderedDict
import math
import matplotlib.backends.backend_pdf as pdf
from scipy.optimize import curve_fit
from math import ceil

SANS_SERIF_FONTS = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'Helvetica']

try:
    from helpers.function import str_exp, plot_uptake
except ImportError:
    # Fallback to the standard stretched-exponential form when the helper module is absent.
    def str_exp(t, k, b):
        return 100 * (1 - np.exp(-np.power(k * t, b)))

    def plot_uptake(k, b, t):
        return str_exp(t, k, b)


pdb_id = os.environ.get('pdb_id', 'glpG-RKRK-79HIS')  # CHECKME
sim_id = os.environ.get('sim_id', 'memb_test')  # CHECKME
hxms_subworkflow = os.environ.get('HXMS_SUBWORKFLOW', 'all').strip().lower()  # CHECKME

work_dir = './'
input_dir = '{}/inputs'.format(work_dir)
output_dir = '{}/outputs'.format(work_dir)
result_dir = '{}/results'.format(work_dir)
run_dir = '{}/{}'.format(output_dir, sim_id)
fig_dir = '{}/{}'.format(output_dir, sim_id)
pdb_dir = '{}/pdb'.format(work_dir)

os.makedirs(fig_dir, exist_ok=True)


def load_uptake_dataframe():
    uptake_path = '{}/{}_rawuptake_HXMS.csv'.format(pdb_dir, pdb_id)
    return pd.read_csv(uptake_path)


def build_state_frames(uptake_df):
    protein_states = np.array(np.unique(uptake_df['Protein State'].astype('str')))
    state_frames = []
    for protein_state in protein_states:
        state_df = uptake_df.loc[uptake_df['Protein State'] == protein_state].copy()
        for column in ['Deut Time (sec)', '%D', 'Theor Uptake #D', 'maxD', 'Start', 'End']:
            state_df[column] = pd.to_numeric(state_df[column], errors='coerce')
        state_frames.append(state_df)
    return protein_states, state_frames


def build_peptide_metadata(reference_df):
    peptides = list(OrderedDict.fromkeys(reference_df['Sequence'].tolist()))
    peptides_plot = np.array(peptides)

    pep_starts = []
    pep_ends = []
    pep_nums = []
    for peptide in peptides:
        pep_start = reference_df.loc[reference_df['Sequence'] == peptide]['Start'].tolist()[0]
        pep_end = reference_df.loc[reference_df['Sequence'] == peptide]['End'].tolist()[0]
        pep_starts.append(pep_start)
        pep_ends.append(pep_end)
        pep_nums.append(str(pep_start) + '--' + str(pep_end))

    return peptides, peptides_plot, np.array(pep_starts), np.array(pep_ends), pep_nums


def export_peptide_ids(pep_starts, pep_ends, peptides_plot):
    with open('{}/{}/{}_pep_ids.csv'.format(output_dir, sim_id, pdb_id), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(pep_starts)
        writer.writerow(pep_ends)
        writer.writerow(peptides_plot)


def build_timepoints(state_frames):
    raw_times = []
    for df in state_frames:
        state_times = pd.to_numeric(df['Deut Time (sec)'], errors='coerce').dropna().tolist()
        raw_times.extend(state_times)

    if not raw_times:
        return np.array([])

    return np.array(sorted(np.unique(raw_times)))


def build_uptake_arrays(state_frames, peptides, timepoints):
    df_arrays = []
    df_kchems = []
    df_maxs = []

    for df in state_frames:
        df_array = np.full((len(peptides), len(timepoints)), np.nan)
        df_kchem = np.full((len(peptides), len(timepoints)), np.nan)
        df_max = np.full((len(peptides), len(timepoints)), np.nan)

        for p, peptide in enumerate(peptides):
            curr_df = df.loc[df['Sequence'] == peptide]
            for t, time in enumerate(timepoints):
                time_data = curr_df.loc[curr_df['Deut Time (sec)'] == time]
                if len(time_data) == 0 or pd.isna(time_data['%D'].iloc[0]):
                    continue

                df_array[p][t] = float(time_data['%D'].iloc[0])
                df_kchem[p][t] = float(time_data['Theor Uptake #D'].iloc[0]) / float(time_data['maxD'].iloc[0]) * 100
                df_max[p][t] = float(time_data['maxD'].iloc[0])

        df_arrays.append(df_array)
        df_kchems.append(df_kchem)
        df_maxs.append(df_max)

    return df_arrays, df_kchems, df_maxs


def build_peptide_indexes(pep_nums):
    pep_list = pep_nums
    pep_ind = np.zeros(len(pep_list))
    for p, peptide in enumerate(pep_list):
        pep_ind[p] = pep_list.index(peptide)
    return pep_list, pep_ind.astype(int)


def safe_normalize(perD, ref_value, maxD):
    with np.errstate(divide='ignore', invalid='ignore'):
        perDnorm = (perD - ref_value) / (maxD - ref_value) * 100
    perDnorm[~np.isfinite(perDnorm)] = np.nan
    return perDnorm


def run_normalized_workflow(state_names, peptides_plot, pep_ind, pep_list, timepoints, df_arrays, df_kchems):
    time_s = 1
    plot_time_s = 2
    plot_time_e = -1
    tp = np.array(timepoints)

    allnorms = []
    df_backexs = []
    for df_array in df_arrays:
        state_backex = np.nanmax(df_array[pep_ind], axis=1)
        for d in range(len(state_backex)):
            state_max = np.nanmax(df_array[pep_ind], axis=1)[d]
            if state_backex[d] < state_max:
                state_backex[d] = state_max
        df_backexs.append(state_backex)

        allnorm = []
        for p in pep_ind:
            perD = df_array[p]
            maxD = state_backex[p]

            if math.isnan(perD[time_s]):
                new = np.nanmin(perD[time_s:])
                perDnorm = safe_normalize(perD, new, maxD)
            else:
                perDnorm = safe_normalize(perD, perD[time_s], maxD)
            allnorm.append(perDnorm)

        allnorms.append(np.array(allnorm))

    n = 4
    m = -(-len(pep_ind) // n)
    maxv = 110
    minv = -5
    reference_kchem = df_kchems[0]

    pdf_d_norm = pdf.PdfPages(fig_dir + '/{}_D_norm_uptake_HXMS.pdf'.format(pdb_id))
    for state_name, df_array, allnorm in zip(state_names, df_arrays, allnorms):
        fig = plt.figure(figsize=(4 * n, 3 * m), facecolor='w')
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = SANS_SERIF_FONTS

        for p in pep_ind:
            masks = ~np.isnan(allnorm[p][plot_time_s:plot_time_e])
            kchem = reference_kchem[p]

            ax = plt.subplot(int(ceil(len(pep_ind) / n)), n, p + 1)
            ax.plot(np.log10(tp[plot_time_s:plot_time_e]), kchem[plot_time_s:plot_time_e], linestyle='dashed', color='grey', label='Theoretical')
            ax.plot(np.log10(tp[plot_time_s:plot_time_e][masks]), allnorm[p][plot_time_s:plot_time_e][masks], '-o', color='black', alpha=0.8, label='{}'.format(state_name))
            ax.plot(6, allnorm[p][plot_time_e], 'o', color='black', alpha=0.8)

            ax.set_xscale('log')
            ax.set_ylim([minv, maxv])
            plt.title(str(pep_list[p] + ' ' + peptides_plot[pep_ind[p]]))

        fig.text(0.5, 1, '%D Uptake for {}'.format(state_name), fontsize=16, ha='center', va='center')
        fig.text(0.0, 0.5, 'Log Time (s)', fontsize=16, ha='center', va='center', rotation='vertical')
        fig.text(0.5, 0, 'Normalized %D', fontsize=16, ha='center', va='center')
        plt.tight_layout()
        pdf_d_norm.savefig(fig, dpi=300, bbox_inches='tight')
        plt.close()

    pdf_d_norm.close()

    for state_name, allnorm in zip(state_names, allnorms):
        np.save('{}/{}/{}_D_norm_uptake_peps_{}.npy'.format(output_dir, sim_id, pdb_id, state_name), [peptides_plot[pep_ind[p]] for p in pep_ind])
        np.save('{}/{}/{}_D_norm_uptake_time_peps_{}.npy'.format(output_dir, sim_id, pdb_id, state_name), [np.log10(tp[plot_time_s:plot_time_e]) for p in pep_ind])
        np.save('{}/{}/{}_D_norm_uptake_d_norm_peps_{}.npy'.format(output_dir, sim_id, pdb_id, state_name), [allnorm[p][plot_time_s:plot_time_e] for p in pep_ind])
        np.save('{}/{}/{}_D_norm_uptake_d_norm_theor_{}.npy'.format(output_dir, sim_id, pdb_id, state_name), [reference_kchem[p][plot_time_s:plot_time_e] for p in pep_ind])


def run_stretched_exp_workflow(state_names, peptides_plot, pep_ind, pep_list, timepoints, df_arrays, df_kchems):
    time_s = 1
    time_e = -1
    tp = np.array(timepoints)
    tp_fit = np.array(timepoints[time_s:time_e])
    t_int = np.geomspace(tp_fit[0], tp_fit[-1].item(), np.size(pep_list))

    df_plots = [df_array[pep_ind] for df_array in df_arrays]
    df_fits = [df_array[pep_ind, time_s:time_e] for df_array in df_arrays]
    df_theors = [df_kchem[pep_ind, time_s:time_e] for df_kchem in df_kchems]

    df_backexs = []
    for df_array in df_arrays:
        state_backex = df_array[pep_ind, -1].copy()
        for d in range(len(state_backex)):
            if math.isnan(state_backex[d]):
                state_backex[d] = np.nanmax(df_array[d])
        df_backexs.append(state_backex)

    smasks = []
    df_fit_fits = []
    df_theor_fits = []
    for df_fit, df_theor in zip(df_fits, df_theors):
        state_masks = []
        state_fit = []
        state_theor = []
        for p in range(len(pep_list)):
            mask = np.isfinite(df_fit[p])
            state_masks.append(mask)
            state_fit.append(df_fit[p][mask])
            state_theor.append(df_theor[p][mask])
        smasks.append(state_masks)
        df_fit_fits.append(state_fit)
        df_theor_fits.append(state_theor)

    fit_Ds = []
    fit_Ths = []
    for df_fit_fit, df_theor_fit, smask, df_backex in zip(df_fit_fits, df_theor_fits, smasks, df_backexs):
        fit_D = []
        fit_Th = []
        for p in range(len(pep_list)):
            if len(df_fit_fit[p]) >= 3 and df_fit_fit[p][-1] > 20:
                try:
                    popt, pcov = curve_fit(
                        str_exp,
                        tp_fit[smask[p]],
                        df_fit_fit[p] / (df_backex[p] / 100),
                        p0=[0.01, 1],
                        bounds=([0, 0], [np.inf, np.inf]),
                        maxfev=5000000,
                    )
                    fit_D.append(plot_uptake(popt[0], popt[1], t_int) * (df_backex[p] / 100))

                    thopt, thcov = curve_fit(
                        str_exp,
                        tp_fit[smask[p]],
                        df_theor_fit[p],
                        p0=[0.01, 1],
                        maxfev=5000000,
                    )
                    fit_Th.append(plot_uptake(thopt[0], thopt[1], t_int))
                except Exception as exc:
                    print('--> WARNING: Could not fit peptide {} ({}). Error: {}'.format(p, pep_list[p], exc))
                    fit_D.append([np.nan for _ in range(len(peptides_plot))])
                    fit_Th.append([np.nan for _ in range(len(peptides_plot))])
            else:
                fit_D.append([np.nan for _ in range(len(peptides_plot))])
                fit_Th.append([np.nan for _ in range(len(peptides_plot))])

        fit_Ds.append(np.array(fit_D))
        fit_Ths.append(np.array(fit_Th))

    n = 4
    m = -(-len(pep_ind) // n)
    maxv = 110
    minv = -5

    pdf_stretch_exp = pdf.PdfPages(fig_dir + '/{}_stretch_exp_HXMS.pdf'.format(pdb_id))
    for state_name, df_plot, fit_d, fit_th in zip(state_names, df_plots, fit_Ds, fit_Ths):
        fig = plt.figure(figsize=(4 * n, 3 * m), facecolor='w')
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = SANS_SERIF_FONTS

        for p in pep_ind:
            series = df_plot[p][1:].astype(np.double)
            masks = np.isfinite(series)

            ax = plt.subplot(int(ceil(len(pep_ind) / n)), n, p + 1)
            ax.plot(t_int, fit_d[p], color='black')
            ax.plot(t_int, fit_th[p], linestyle='dashed', color='grey')
            ax.plot(tp[1:][masks], df_plot[p][1:][masks], 'o', color='black')

            ax.set_xscale('log')
            ax.set_ylim([minv, maxv])
            plt.title(str(pep_list[p] + ' ' + peptides_plot[pep_ind[p]]))

        fig.text(0.5, 1, '%D Uptake for {}'.format(state_name), fontsize=16, ha='center', va='center')
        fig.text(0.0, 0.5, 'Log Time (s)', fontsize=16, ha='center', va='center', rotation='vertical')
        fig.text(0.5, 0, 'Normalized %D', fontsize=16, ha='center', va='center')
        plt.tight_layout()
        pdf_stretch_exp.savefig(fig, dpi=300, bbox_inches='tight')
        plt.close()

    pdf_stretch_exp.close()

    for state_name, fit_d, fit_th in zip(state_names, fit_Ds, fit_Ths):
        np.save('{}/{}/{}_stretch_exp_peps_{}.npy'.format(output_dir, sim_id, pdb_id, state_name), [peptides_plot[p] for p in range(len(pep_list))])
        np.save('{}/{}/{}_stretch_exp_time_peps_{}.npy'.format(output_dir, sim_id, pdb_id, state_name), np.log10(t_int))
        np.save('{}/{}/{}_stretch_exp_d_norm_peps_{}.npy'.format(output_dir, sim_id, pdb_id, state_name), [fit_d[p] for p in range(len(pep_list))])
        np.save('{}/{}/{}_stretch_exp_d_norm_theor_{}.npy'.format(output_dir, sim_id, pdb_id, state_name), [fit_th[p] for p in range(len(pep_list))])


uptake_df = load_uptake_dataframe()
state_names, state_frames = build_state_frames(uptake_df)
peptides, peptides_plot, pep_starts, pep_ends, pep_nums = build_peptide_metadata(state_frames[0])
timepoints = build_timepoints(state_frames)

export_peptide_ids(pep_starts, pep_ends, peptides_plot)

df_arrays, df_kchems, df_maxs = build_uptake_arrays(state_frames, peptides, timepoints)
pep_list, pep_ind = build_peptide_indexes(pep_nums)

if hxms_subworkflow not in {'all', 'normalized', 'stretched'}:
    raise SystemExit('Unsupported HXMS_SUBWORKFLOW: {}'.format(hxms_subworkflow))

if hxms_subworkflow in {'all', 'normalized'}:
    run_normalized_workflow(state_names, peptides_plot, pep_ind, pep_list, timepoints, df_arrays, df_kchems)

if hxms_subworkflow in {'all', 'stretched'}:
    run_stretched_exp_workflow(state_names, peptides_plot, pep_ind, pep_list, timepoints, df_arrays, df_kchems)
