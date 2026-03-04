import mdtraj as md
import mdtraj_upside as mu
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib.font_manager import FontProperties
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis
import os # for high-throughput analysis
import csv # for high-throughput analysis
from scipy.integrate import trapezoid # for COF integration
import pandas as pd
import sys
import math
import matplotlib.backends.backend_pdf as pdf # for D uptake and COF
import matplotlib.patches as patches # for peptide plots
from scipy.stats import percentileofscore # for COF comparison
from matplotlib import colors # for COF colormap
from pyhdx.plot import apply_cmap # for COF colormap
from function import *
import matplotlib.gridspec as gridspec # for COF subplots

pdb_id      = 'glpG-RKRK-79HIS'
#pdb_id      = os.environ['pdb_id']
sim_id      = 'memb_test'
#sim_id      = os.environ['sim_id'] #CHECKME

main_pdb = 'glpG-RKRK-79HIS' #CHECKME
HXMS_method = 'stretch_exp' #CHECKME
#HXMS_method = 'D_norm_uptake'
protein_state = 'pd9' #CHECKME

work_dir         = './'
n_rep            = 48     #CHECKME # replica number
#n_rep            = os.environ['n_rep'] #CHECKME

start_frame = 100 #CHECKME

pdb_input=[
"glpG-RKRK-79HIS",
]

pdb_ids = np.array(pdb_input)

output_dir = "{}/outputs".format(work_dir)
pdb_dir    = "{}/pdb".format(work_dir)
fig_dir    = "results"

input_dirs = []
result_dirs = []
run_dirs = []
for pdb_id in pdb_ids:
    input_dir  = "{}/inputs/{}".format(work_dir, pdb_id)
    result_dir = "{}/results".format(work_dir)
    run_dir    = "{}/{}/{}".format(output_dir, sim_id, pdb_id)

    input_dirs.append(input_dir)
    result_dirs.append(result_dir)
    run_dirs.append(run_dir)

#============================================
# load data
#============================================

# D uptake Upside predictions across proteins/protein domains
times = []
d_norms = []
for p, (pdb_id, result_dir) in enumerate(zip(pdb_ids, result_dirs)):
    with open('{}/{}_percentD.csv'.format(result_dir, pdb_id), 'r', encoding='utf-8-sig') as f:
        rows = csv.reader(f)

        time_row = []
        d_norm_row = []
        for i, row in enumerate(rows):
            if i == 0:
                time_row = row

            elif i == 1:
                d_norm_row = row

            else:
                break

    times.append(time_row)
    d_norms.append(d_norm_row)

# D uptake Upside predictions for experimental peptides
peps = []
peps_seq = np.array([])
time_peps = []
d_norm_peps = []
d_norm_peps_theors = []
for p, (pdb_id, result_dir) in enumerate(zip(pdb_ids, result_dirs)):
    pep = np.load('{}/{}_percentD_peps.npy'.format(result_dir, pdb_id))
    peps.append(pep)
    peps_seq = np.concatenate([peps_seq, pep])
    time_peps.append(np.load('{}/{}_percentD_time_peps.npy'.format(result_dir, pdb_id)))
    d_norm_peps.append(np.load('{}/{}_percentD_d_norm_peps.npy'.format(result_dir, pdb_id)))
    d_norm_peps_theors.append(np.load('{}/{}_percentD_d_norm_theor.npy'.format(result_dir, pdb_id)))

# identifiers for Upside peptides
#pep_starts = np.array([])
#pep_ends = np.array([])
#pep_ids = np.array([])
#for p, (pdb_id, result_dir) in enumerate(zip(pdb_ids, result_dirs)):
    #with open('{}/{}_pep_ids.csv'.format(result_dir, pdb_id), 'r', encoding='utf-8-sig') as f:
        #rows = csv.reader(f)

        #start_row = []
        #end_row = []
        #pep_row = []
        #for i, row in enumerate(rows):
            #if i == 0:
                #start_row = row

            #elif i == 1:
                #end_row = row

            #elif i == 2:
                #pep_row = row

            #else:
                #break

    #pep_starts = np.concatenate([pep_starts, np.array(start_row)])
    #pep_ends = np.concatenate([pep_ends, np.array(end_row)])
    #pep_ids = np.concatenate([pep_ids, np.array(pep_row)])

# normalized D uptake / stretched exponential for experimental peptides
peps_HXMS = np.load('{}/{}/{}_{}_peps_{}.npy'.format(output_dir, sim_id, main_pdb, HXMS_method, protein_state))

time_peps_HXMS = np.load('{}/{}/{}_{}_time_peps_{}.npy'.format(output_dir, sim_id, main_pdb, HXMS_method, protein_state))

d_norm_peps_HXMS = np.load('{}/{}/{}_{}_d_norm_peps_{}.npy'.format(output_dir, sim_id, main_pdb, HXMS_method, protein_state))

d_norm_theor_HXMS = np.load('{}/{}/{}_{}_d_norm_theor_{}.npy'.format(output_dir, sim_id, main_pdb, HXMS_method, protein_state))

# identifiers for experimental peptides
with open('{}/{}/{}_pep_ids.csv'.format(output_dir, sim_id, main_pdb), 'r', encoding='utf-8-sig') as f:
    rows = csv.reader(f)

    pep_start = []
    pep_end = []
    pep_id = []
    for i, row in enumerate(rows):
        if i == 0:
            pep_start = row

        elif i == 1:
            pep_end = row

        elif i == 2:
            pep_id = row

        else:
            break

pep_start = np.array(pep_start)
pep_end = np.array(pep_end)
pep_id = np.array(pep_id)

fasta = np.loadtxt('{}/inputs/{}.fasta'.format(work_dir, main_pdb, main_pdb), skiprows = 1, dtype=str)

if np.size(fasta) == 1:
    sequence = fasta.tolist()

else:
    sequence = ''.join(fasta)

#============================================
# calculate dG_HX at different T 
#============================================

dGhx_Ts = np.array([])
for p, (pdb_id, result_dir) in enumerate(zip(pdb_ids, result_dirs)):
    PS = []
    T  = []
    for i in range(n_rep):
        PS.append( np.load('{}/{}_{}_{}_PS.npy'  .format(result_dir, pdb_id, sim_id, i)) ) 
        T.append(np.load( '{}/{}_{}_{}_T.npy'.format(result_dir, pdb_id, sim_id, i) ))

    min_idx = np.min([len(PS[i][:,0]) for i in range(n_rep)])
    PS_idx = [PS[i][:min_idx] for i in range(n_rep)]
    PS   = np.array(PS_idx)
    T    = np.array(T)

    res = np.loadtxt('{}/{}.resid'.format(result_dir, pdb_id), dtype=int)
    n_res = res.size

    dGhx_T = []
    for k in range(T.size):
        t = T[k]
        tt = t/0.85*298.

        dG = np.zeros(n_res)
        for r in range(n_res):
            pf_i = PS[:,start_frame:,r].flatten()
            mean_pf = np.average(pf_i)
            #mean_pf = np.average(pf_i, weights=w1)
            if mean_pf == 1:
            #print(k, r)
                dG[r] = 1000.
            else:
                dG[r] = 0.001987*tt*np.log((mean_pf/(1.-mean_pf)))
        dGhx_T.append(dG)
    dGhx_T = np.array(dGhx_T)

    dGhx_Ts = np.concatenate([dGhx_Ts, dGhx_T[0]], axis=0)

dGhx_T_plots = []
for r in range(len(sequence)):
    if r in res:
        dGhx_T_plot = dGhx_Ts[np.argwhere(res==r)][0][0]
        dGhx_T_plots.append(dGhx_T_plot)

    else:
        dGhx_T_plot = np.nan
        dGhx_T_plots.append(dGhx_T_plot)

dGhx_T_plots = np.array(dGhx_T_plots)

#============================================
# compare D uptake at peptide level
#============================================

# build dataframe of experimental start and stop indices for peptides
pep_exp = {"peps": pep_id, "pep_starts": pep_start, "pep_ends": pep_end}
pep_exp_df = pd.DataFrame(pep_exp)

# normalize the stretched exponential data
d_norms_HXMS = normalize_se(d_norm_peps_HXMS)

# build dataframe of D uptakes for experimental peptides
HXMS = {"peps": peps_HXMS, "d_norms_HXMS": d_norms_HXMS, "d_theor_HXMS": d_norm_theor_HXMS.tolist()}
HXMS_df = pd.DataFrame(HXMS)

# build dataframe of computational start and stop indices for peptides
#pep_ind = {"peps": pep_ids, "pep_starts": pep_starts, "pep_ends": pep_ends}
#pep_ind_df = pd.DataFrame(pep_ind)

# build dataframe of D uptakes and COFs for computational peptides
Up_dfs = []
for p, (p_id_Up, time_pep, d_norm_Up, d_norm_theor) in enumerate(zip(peps, time_peps, d_norm_peps, d_norm_peps_theors)):

    # normalize the stretched exponential data
    d_norms_theor = normalize_se(d_norm_theor)
    d_norms_Up = normalize_se(d_norm_Up)

    # get dataframes of D uptakes for each protein/domain
    Up = {"peps": p_id_Up, "time_peps": [time_pep for i in enumerate(p_id_Up)], "d_norms_Up": d_norms_Up, "d_norms_theor": d_norms_theor}
    #Up = {"peps": p_id_Up, "time_peps": [time_pep for i in enumerate(p_id_Up)], "d_norms_Up": d_norm_Up.tolist(), "d_norms_theor": d_norm_theor.tolist()}
    Up_df = pd.DataFrame(Up)
    Up_dfs.append(Up_df)

# append dataframes for all protein domains
Up_dfs = pd.concat(Up_dfs)

# merge dataframes of experimental and computational peptides
merge_df = pd.merge(pep_exp_df, HXMS_df, on="peps")
merge_dfs = pd.merge(merge_df, Up_dfs, on="peps")

#============================================
# calculate cooperativity factor (COF)
#============================================

# calculate the squared derivative and its integral for experimental peptides
integrals_pep_HXMS = []
for i, (p_id, d_norm_pep) in enumerate(zip(merge_dfs.peps, merge_dfs.d_norms_HXMS)):

    if np.size(d_norm_pep) == np.size(time_peps_HXMS): # ensure sizes match
        derivative = numerical_derivative(time_peps_HXMS, d_norm_pep)
        squared_derivative = np.square(derivative) # square the derivative

        # calculate the integral of the squared derivative using trap rule
        integral = trapezoid(squared_derivative, time_peps_HXMS)

    else:
        break

    integrals_pep_HXMS.append(integral)

merge_dfs['COF_HXMS'] = integrals_pep_HXMS

# calculate the squared derivative and its integral for computational peptides
integrals_pep_Upside = []
for p, (p_id, d_norm_pep) in enumerate(zip(merge_dfs.peps, merge_dfs.d_norms_Up)):

    if np.size(d_norm_pep) == np.size(time_peps[0]): # ensure sizes match
        derivative = numerical_derivative(time_peps[0], d_norm_pep)
        squared_derivative = np.square(derivative) # square the derivative

        # calculate the integral of the squared derivative using trap rule
        integral = trapezoid(squared_derivative, time_peps[0])

    else:
        break

    integrals_pep_Upside.append(integral)

merge_dfs['COF_Up'] = integrals_pep_Upside

#============================================
# plot D uptake curves
#============================================

maxv = 100
minv = 0

# plot Upside predicted uptake curves
#fig = plt.figure()
#ax = plt.subplot(111)
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = ['Arial']

#for i, (pdb_id, time, d_norm) in enumerate(zip(pdb_ids, times, d_norms)):
    #ax.plot(np.power(10, np.array(time, dtype=float)), np.array(d_norm, dtype=float), label='{}'.format(pdb_id))
#ax.set_xscale('log')
#x_ticks = np.arange(min(times), max(times)+np.average(np.diff(times)), np.average(np.diff(times)))
#ax.set_xticks(x_ticks, [f'{x:.2e}' for x in x_ticks])
#ax.set_ylim(minv, maxv)
#ax.set_xlabel('Log Time (s)')
#ax.set_ylabel('Normalized %D')
#ax.set_title('Upside %D Uptake Across Protein Domains')

#fontP = FontProperties()
#fontP.set_size('x-small')

#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop=fontP, columnspacing=0.8)
#plt.savefig(fig_dir + '/D_uptake_Upside_{}.png'.format(main_pdb))
#plt.close()

# plot D uptake for HX-MS vs. Upside at peptide-level
# pdf_d = pdf.PdfPages(fig_dir + '/D_uptake_{}.pdf'.format(main_pdb))

# df_nan_bool = np.logical_and(np.logical_and(~(np.array([np.isnan(x[0]) for x in merge_dfs.d_norms_HXMS.values])), ~(np.array([np.isnan(x[0]) for x in merge_dfs.d_theor_HXMS.values]))), np.logical_and(~(np.array([np.isnan(x[0]) for x in merge_dfs.d_norms_Up.values])), ~(np.array([np.isnan(x[0]) for x in merge_dfs.d_norms_theor.values]))))

# for p, (pep, time_pep, d_norm_HXMS, d_theor_HXMS, d_norm_Up, d_norm_theor, p_start, p_end) in enumerate(zip(merge_dfs.peps[df_nan_bool], merge_dfs.time_peps[df_nan_bool], merge_dfs.d_norms_HXMS.values[df_nan_bool], merge_dfs.d_theor_HXMS[df_nan_bool], merge_dfs.d_norms_Up[df_nan_bool], merge_dfs.d_norms_theor[df_nan_bool], merge_dfs.pep_starts[df_nan_bool], merge_dfs.pep_ends[df_nan_bool])):

#     d_norm_stats = sp.stats.linregress(d_norm_HXMS, d_norm_Up)
#     r_square = d_norm_stats.rvalue**2

#     fig = plt.figure()
#     ax = plt.subplot(111)
#     plt.rcParams['font.family'] = 'sans-serif'
#     plt.rcParams['font.sans-serif'] = ['Arial']

#     #plt.plot(np.log10(time_peps_HXMS), d_norm_HXMS, label='HX-MS Experiment')
#     #plt.plot(np.log10(time_peps, d_norm_Up), label='Upside Simulation')
#     ax.plot(np.power(10, time_peps_HXMS), d_norm_HXMS, label='HX-MS Experiment')
#     ax.plot(np.power(10, time_peps_HXMS), d_theor_HXMS, label='HX-MS Theoretical', linestyle='dashed', color='black')
#     ax.plot(np.power(10, time_pep), d_norm_Up, label='Upside Simulation')
#     ax.plot(np.power(10, time_pep), d_norm_theor, label='Upside Theoretical', linestyle='dashed', color='grey')
#     ax.set_xscale('log')
#     x_ticks = np.arange(min(time_pep), max(time_pep)+np.average(np.diff(time_pep)), np.average(np.diff(time_pep)))
#     ax.set_xticks(x_ticks, [f'{x:.2e}' for x in x_ticks])
#     ax.set_ylim(minv, maxv)
#     ax.legend()
#     ax.set_xlabel('Log Time (s)')
#     ax.set_ylabel('%D Uptake')
#     ax.set_title('{}—{} {}'.format(p_start, p_end, pep))

#     plt.tight_layout()
#     pdf_d.savefig(fig, dpi=300, bbox_inches = 'tight')
#     plt.close()

# pdf_d.close()

#============================================
# plot bar charts and histograms of COFs
#============================================

# plot COFs for HX-MS vs. Upside
pdf_cof = pdf.PdfPages(fig_dir + '/COF_{}.pdf'.format(main_pdb))

# plot bar chart of COF integrals
#fig, axes = plt.subplots(1, 2, figsize=(15, 5))
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = ['Arial']

#axes[0].bar(range(1, np.size(integrals_pep_HXMS) + 1), integrals_pep_HXMS, edgecolor='black')
#axes[0].set_xlabel('Fit Number')
#axes[0].set_ylabel('Integral of Squared Derivative')
#axes[0].set_title('COFs of HX-MS Peptides')
#axes[0].grid(axis='y')

#axes[1].bar(range(1, np.size(integrals_seq_Upside) + 1), integrals_seq_Upside, edgecolor='black')
#axes[1].set_xlabel('Fit Number')
#axes[1].set_ylabel('Integral of Squared Derivative')
#axes[1].set_title('COFs of Upside Peptides')
#axes[1].grid(axis='y')

#plt.tight_layout()
#pdf_cof.savefig(fig, dpi=300, bbox_inches = 'tight')
#plt.close()

# plot histogram of COF integrals
#fig, axes = plt.subplots(1, 2, figsize=(15, 5))
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = ['Arial']

#axes[0].hist(integrals_pep_HXMS, bins=40, edgecolor='black', density=True)
#axes[0].set_xlim(0, np.nanmax(integrals_pep_HXMS))
#axes[0].set_xlabel('Cooperativity Factor')
#axes[0].set_ylabel('Normalized Frequency')
#axes[0].set_title('COFs of HX-MS Peptides')
#axes[0].grid(True)

#axes[1].hist(integrals_seq_Upside, bins=40, edgecolor='black', density=True)
#axes[1].set_xlim(0, np.nanmax(integrals_seq_Upside))
#axes[1].set_xlabel('Cooperativity Factor')
#axes[1].set_ylabel('Normalized Frequency')
#axes[1].set_title('COFs of Upside Peptides')
#axes[1].grid(True)

#plt.tight_layout()
#pdf_cof.savefig(fig, bbox_inches = 'tight')
#plt.close()

#============================================
# plot peptide plots across protein
#============================================

def peptide_colors(merged_df, rows, df_row, quant):
    """
    set colors based on distributions of COFs from multiple proteins
    """
    colors = []
    for index, row in merged_df.iterrows():
        if quant == True: # use quantile coloring
            if pd.isna(row['{}'.format(rows)]):
                colors.append('darkblue')
            elif row['{}'.format(rows)] < np.nanquantile(df_row, 0.10):
                colors.append('mediumblue')
            elif np.nanquantile(df_row, 0.10) <= row['{}'.format(rows)] < np.nanquantile(df_row, 0.20):
                colors.append('blue')
            elif np.nanquantile(df_row, 0.20) <= row['{}'.format(rows)] < np.nanquantile(df_row, 0.30):
                colors.append('slateblue')
            elif np.nanquantile(df_row, 0.30) <= row['{}'.format(rows)] < np.nanquantile(df_row, 0.40):
                colors.append('mediumslateblue')
            elif np.nanquantile(df_row, 0.40) <= row['{}'.format(rows)] < np.nanquantile(df_row, 0.50):
                colors.append('blueviolet')
            elif np.nanquantile(df_row, 0.50) <= row['{}'.format(rows)] < np.nanquantile(df_row, 0.60):
                colors.append('darkorchid')
            elif np.nanquantile(df_row, 0.60) <= row['{}'.format(rows)] < np.nanquantile(df_row, 0.70):
                colors.append('mediumorchid')
            elif np.nanquantile(df_row, 0.70) <= row['{}'.format(rows)] < np.nanquantile(df_row, 0.80):
                colors.append('orchid')
            elif np.nanquantile(df_row, 0.80) <= row['{}'.format(rows)] < np.nanquantile(df_row, 0.90):
                colors.append('violet')
            else:
                colors.append('magenta')

        elif quant == False: # use absolute coloring
            if pd.isna(row['{}'.format(rows)]):
                colors.append('blue') # PD6 b factor = zero
            elif row['{}'.format(rows)] < 500:
                colors.append('darkturquoise') # PD6 b factor within 20% of kchem
            elif 500 <= row['{}'.format(rows)] < 1000:
                colors.append('darkorchid') # PD6 b factor ≥ kchem
            else:
                colors.append('magenta') # PD6 b factor < kchem

        else:
            break

    return colors

#colors_HXMS_quant = peptide_colors(merge_df, 'COF_HXMS', merge_df.COF_HXMS, True)
#colors_HXMS = peptide_colors(merge_dfs, 'COF_HXMS', merge_dfs.COF_HXMS, False)
#colors_Upside_quant = peptide_colors(merge_dfs, 'COF_Up', merge_dfs.COF_Up, True)
#colors_Upside = peptide_colors(merge_dfs, 'COF_Up', merge_dfs.COF_Up, False)

# define colormap and normalization
cmap = matplotlib.colormaps.get_cmap('bwr').reversed()
#norm_HXMS_res = colors.TwoSlopeNorm(vmin=0, vcenter=np.nanmedian(merge_dfs['COF_HXMS']), vmax=np.nanmax(merge_dfs['COF_HXMS']))
norm_HXMS_res = colors.TwoSlopeNorm(vmin=-np.nanmax(merge_dfs['COF_HXMS']), vcenter=-np.nanmedian(merge_dfs['COF_HXMS']), vmax=-np.nanmin(merge_dfs['COF_HXMS']))
#norm_Up_res = colors.TwoSlopeNorm(vmin=0, vcenter=np.nanmedian(merge_dfs['COF_Up']), vmax=np.nanmax(merge_dfs['COF_Up']))
norm_Up_res = colors.TwoSlopeNorm(vmin=-np.nanmax(merge_dfs['COF_Up']), vcenter=-np.nanmedian(merge_dfs['COF_Up']), vmax=-np.nanmin(merge_dfs['COF_Up']))

colors_HXMS = apply_cmap(-merge_dfs['COF_HXMS'], cmap, norm_HXMS_res)
colors_Upside = apply_cmap(-merge_dfs['COF_Up'], cmap, norm_Up_res)

# set segment length, # of residues in a single plot for the map
segment_length = 250 #CHECKME

# break sequence into segments based on segment length
sequence_segments = [sequence[i:i+segment_length] for i in range(0, len(sequence), segment_length)]
num_segments = len(sequence_segments)

# set initial height and increment for staggered bars
initial_height = 0.02
increment = 0.02

def bars_overlap(start1, end1, start2, end2):
    """
    check whether two bars overlap to stack both correctly
    """
    return not (end1 < start2 or end2 < start1)

# calculate adjusted heights
heights_HXMS = [initial_height] * len(merge_dfs)
heights_Upside = [initial_height] * len(merge_dfs)

def height_overlap(merged_df, heights):
    """
    determine whether there are overlaps and adjust heights
    """
    for i in range(len(merged_df)):
        for j in range(i + 1, len(merged_df)):
            if bars_overlap(merged_df.iloc[i]['pep_starts'], merged_df.iloc[i]['pep_ends'],
                        merged_df.iloc[j]['pep_starts'], merged_df.iloc[j]['pep_ends']):
                heights[j] = max(heights[j], heights[i] + increment)

height_overlap(merge_dfs, heights_HXMS)
height_overlap(merge_dfs, heights_Upside)

def peptide_plot(merged_df, heights, colors):
    """
    plot peptide plots across protein sequence
    """
    for segment_index, segment in enumerate(sequence_segments):
        fig, ax = plt.subplots(figsize=(60, 10))
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial']

        start_pos_segment = segment_index * segment_length
        end_pos_segment = start_pos_segment + segment_length - 1

        for index, row in merged_df.iterrows():
            if start_pos_segment <= int(row['pep_starts']) <= end_pos_segment:
                start_pos = max(int(row['pep_starts']), start_pos_segment) - start_pos_segment
                end_pos = min(int(row['pep_ends']), end_pos_segment) - start_pos_segment
                height = heights[index] # use the calculated height
                bar = patches.Rectangle((start_pos, height), end_pos - start_pos + 1, 0.01, color=colors[index])
                ax.add_patch(bar)

        # set axis limits and labels for each segment
        ax.set_xlim(0, segment_length + 1)
        ax.set_ylim(0, max(heights) + 0.1) # increase y-space to avoid clipping

        # set x-ticks and x-tick labels to match the sequence positions
        x_ticks = range(1, segment_length + 1)
        sequence_segment = segment
        x_tick_labels = [sequence_segment[i] if i < len(sequence_segment) else '' for i in range(len(sequence_segment))]

        # test whether x_tick_labels have the same length as x_ticks
        if len(x_tick_labels) < len(x_ticks):
            x_tick_labels.extend([''] * (len(x_ticks) - len(x_tick_labels)))

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels[:len(x_ticks)], fontsize=20)

        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        plt.tight_layout()
        pdf_cof.savefig(fig, dpi=300, bbox_inches = 'tight')
        plt.close()

#peptide_plot(merge_df, heights_HXMS, colors_HXMS_quant)
peptide_plot(merge_dfs, heights_HXMS, colors_HXMS.values)
#peptide_plot(merge_dfs, heights_Upside, colors_Upside_quant)
peptide_plot(merge_dfs, heights_Upside, colors_Upside.values)

def reg_plot(HXMS_value, Up_value, reg_method):
    """
    plot regression plots across protein sequence
    """
    # determine r square value
    COF_stats = sp.stats.linregress(HXMS_value.values[np.logical_and(~(np.isnan(HXMS_value.values)), ~(np.isnan(Up_value.values)))], Up_value.values[np.logical_and(~(np.isnan(HXMS_value.values)), ~(np.isnan(Up_value.values)))])
    r_square = COF_stats.rvalue**2

    # plot regression for COF HXMS vs. Upside
    fig = plt.figure()
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.scatter(HXMS_value.values[np.logical_and(~(np.isnan(HXMS_value.values)), ~(np.isnan(Up_value.values)))], Up_value.values[np.logical_and(~(np.isnan(HXMS_value.values)), ~(np.isnan(Up_value.values)))])
    plt.xlabel('COF HXMS')
    plt.ylabel('COF Upside')
    plt.title('COF for HXMS vs. Upside')

    a, _, _, _ = np.linalg.lstsq(HXMS_value.values[:,np.newaxis][np.logical_and(~(np.isnan(HXMS_value.values)), ~(np.isnan(Up_value.values)))], Up_value.values[np.logical_and(~(np.isnan(HXMS_value.values)), ~(np.isnan(Up_value.values)))], rcond=None)
    line = np.linspace(0, 8000, 100)
    plt.plot(line, a*line, color='k', lw=2.5)

    plt.text(0.5, 7500, '$R^2$ = ' + str(round(r_square, 8)))

    plt.xlim(0, 8000)
    plt.ylim(0, 8000)
    plt.savefig(fig_dir + '/COF_HXMS_Up_R2_{}_{}.png'.format(reg_method, main_pdb))
    plt.close()

# plot regression for COF HXMS vs. Upside at peptide-level
reg_plot(merge_dfs.COF_HXMS, merge_dfs.COF_Up, 'pep')

comb_res_HXMS, comb_res_idx_HXMS, comb_res_COF_HXMS = average_int(merge_dfs, 'COF_HXMS', 'peps', 'pep_starts')

comb_res_Up, comb_res_idx_Up, comb_res_COF_Up = average_int(merge_dfs, 'COF_Up', 'peps', 'pep_starts')

# build dataframe of integrals for experimental residues
res_int_HXMS = {"res": comb_res_HXMS, "res_idx": comb_res_idx_HXMS, "res_COF_HXMS": comb_res_COF_HXMS}
res_int_df_HXMS = pd.DataFrame(res_int_HXMS)

res_avg_int_df_HXMS = res_int_df_HXMS.groupby(['res_idx', 'res'])[['res_COF_HXMS']].mean()

# build dataframe of integrals for computational residues
res_int_Up = {"res": comb_res_Up, "res_idx": comb_res_idx_Up, "res_COF_Up": comb_res_COF_Up}
res_int_df_Up = pd.DataFrame(res_int_Up)

res_avg_int_df_Up = res_int_df_Up.groupby(['res_idx', 'res'])[['res_COF_Up']].mean()

# merge dataframes of integrals
res_avg_int_df = pd.merge(res_avg_int_df_Up, res_avg_int_df_HXMS, on=['res_idx', 'res'])

# plot regression for COF HXMS vs. Upside at residue-level
reg_plot(res_avg_int_df.res_COF_HXMS, res_avg_int_df.res_COF_Up, 'res')

#import pdb
#pdb.set_trace()

# define colormap and normalization
cmap = matplotlib.colormaps.get_cmap('bwr').reversed()
#norm_HXMS = colors.TwoSlopeNorm(vmin=0, vcenter=np.nanmedian(res_avg_int_df_HXMS['res_COF_HXMS']), vmax=np.nanmax(res_avg_int_df_HXMS['res_COF_HXMS']))
norm_HXMS = colors.TwoSlopeNorm(vmin=-np.nanmax(res_avg_int_df_HXMS['res_COF_HXMS']), vcenter=-np.nanmedian(res_avg_int_df_HXMS['res_COF_HXMS']), vmax=-np.nanmin(res_avg_int_df_HXMS['res_COF_HXMS']))
#norm_Up = colors.TwoSlopeNorm(vmin=0, vcenter=np.nanmedian(res_avg_int_df_Up['res_COF_Up']), vmax=np.nanmax(res_avg_int_df_Up['res_COF_Up']))
norm_Up = colors.TwoSlopeNorm(vmin=-np.nanmax(res_avg_int_df_Up['res_COF_Up']), vcenter=-np.nanmedian(res_avg_int_df_Up['res_COF_Up']), vmax=-np.nanmin(res_avg_int_df_Up['res_COF_Up']))

diff_color_HXMS = apply_cmap(-res_avg_int_df_HXMS, cmap, norm_HXMS)
diff_color_Up = apply_cmap(-res_avg_int_df_Up, cmap, norm_Up)

def peptide_plot_res(merged_df, colors, subplot):
    """
    plot peptide plots across protein sequence
    """
    fig, ax = plt.subplots(figsize=(60, 10))
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']

    for index, row in merged_df.reset_index().iterrows():
        for segment_index, segment in enumerate(sequence):
            if (int(row['res_idx']) == segment_index+1) and (row['res'] == segment):
                height = (index + 1) * 0.02 # use the calculated height
                bar = patches.Rectangle((index, height), 0.01, 0.01, color=colors[index])
                ax.add_patch(bar)

    # set axis limits and labels for each segment
    ax.set_xlim(0, len(sequence) + 1)
    ax.set_ylim(0, ((len(colors) + 1) * 0.02) + 0.1) # increase y-space to avoid clipping

    # set x-ticks and x-tick labels to match the sequence positions
    x_ticks = range(1, len(sequence) + 1)
    x_tick_labels = [c for i, c in enumerate(sequence)]

    # test whether x_tick_labels have the same length as x_ticks
    if len(x_tick_labels) < len(x_ticks):
        x_tick_labels.extend([''] * (len(x_ticks) - len(x_tick_labels)))

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels[:len(x_ticks)], fontsize=20)

    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    if subplot == 'stability':
        twin = ax.twinx()
        for i, c in enumerate(sequence):
            twin.plot(dGhx_T_plots[i])

        twin.set_ylim(0, np.nanmax(dGhx_T_plots))
        twin.set_ylabel('∆G (kcal/mol)')

        plt.tight_layout()
        pdf_cof.savefig(fig, dpi=300, bbox_inches = 'tight')
        plt.close()

    elif subplot == 'helicity':
        twin = ax.twinx()

        helicity_avg = []
        for p, (pdb_id, run_dir) in enumerate(zip(pdb_ids, run_dirs)):
            helicity = []
            for i in range(n_rep):
                traj = mu.load_upside_traj('{}/{}.run.{}.up'.format(run_dir, pdb_id, i))
                dssp = md.compute_dssp(traj)
                time = traj._time * 1e-6
                #dssp.shape, time.shape
                helicity_matrix = (dssp=='H').astype('int').T
                helicity.append(np.average(helicity_matrix, axis=1)*100)

            helicity_avg.append(np.average(helicity, axis=0))

        helicity_avg = np.array(helicity_avg)

        for i, c in enumerate(sequence):
            twin.plot(helicity_avg[i])

        twin.set_ylim(0, 100)
        twin.set_ylabel("Helicity (%)")

        plt.tight_layout()
        pdf_cof.savefig(fig, dpi=300, bbox_inches = 'tight')
        plt.close()

    else:
        plt.tight_layout()
        pdf_cof.savefig(fig, dpi=300, bbox_inches = 'tight')
        plt.close()


pdf_cof.close()

#import pdb
#pdb.set_trace()

#peptide_plot_res(res_int_df_HXMS[['res_idx', 'res']].drop_duplicates(), diff_color_HXMS['res_COF_HXMS'].values, 'stability')
#peptide_plot_res(res_int_df_HXMS[['res_idx', 'res']].drop_duplicates(), diff_color_HXMS['res_COF_HXMS'].values, 'helicity')
#peptide_plot_res(res_int_df_Up[['res_idx', 'res']].drop_duplicates(), diff_color_Up['res_COF_Up'].values, 'stability')
#peptide_plot_res(res_int_df_Up[['res_idx', 'res']].drop_duplicates(), diff_color_Up['res_COF_Up'].values, 'helicity')

# export colors to .npy file
np.save('{}/{}/{}_colors_res_HXMS.npy'.format(work_dir, fig_dir, main_pdb), res_int_df_HXMS['res_idx'].unique())
np.save('{}/{}/{}_colors_HXMS.npy'.format(work_dir, fig_dir, main_pdb), np.array([str(i) for i in diff_color_HXMS['res_COF_HXMS'].values]))

np.save('{}/{}/{}_colors_res_Up.npy'.format(work_dir, fig_dir, main_pdb), res_int_df_Up['res_idx'].unique())
np.save('{}/{}/{}_colors_Upside.npy'.format(work_dir, fig_dir, main_pdb), np.array([str(i) for i in diff_color_Up['res_COF_Up'].values]))
