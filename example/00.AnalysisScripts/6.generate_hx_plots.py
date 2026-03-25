import os
import pickle
import sqlite3
from pathlib import Path

import matplotlib

if 'MPLBACKEND' not in os.environ:
    matplotlib.use('Agg')

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler

from helpers import hxfunctions_clean as hxf


DEFAULT_PREFIX = 'glpg'
DEFAULT_STATE = 'pd9'
DEFAULT_TEMPERATURE_K = 296.0
DEFAULT_PD_CORR = 9.0
DEFAULT_RESIDUE_OFFSET = 66
DEFAULT_DSEQ = -1
DEFAULT_SEQUENCE = (
    "GSSHHHHHHSSGLVPRGSHMAALRERAGPVTWVMMIACVVVFIAMQILGDQEVMLWLAWPFDPTLKFEFWRYFTHALMHFSLMHILFNLLWWWYLGGAVEKRLGSGKLIVITLISALLSGYVQQKFSGPWFGGLTGVVYALMGYVWLRGERDPQSGIYLQRGLIIFALIWIVAGWFDLFGMSMANGAHIAGLAVGLAMAFVDSLNARKRK"
)
DEFAULT_SS_RANGES = [[29, 47], [82, 102], [106, 126], [135, 151], [162, 175], [183, 207]]
DEFAULT_UPSIDE_TEMPS = [0.70, 0.75, 0.80, 0.85, 0.90]


mycolors = mpl.color_sequences['tab10']
plt.rc('axes', prop_cycle=(cycler(color=mycolors)))

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Tahoma', 'DejaVu Sans', 'Lucida Grande', 'Verdana']
plt.rcParams['font.family'] = 'monospace'
plt.rcParams['font.monospace'] = ['Source Code Pro', 'Consolas', 'Courier New']
plt.rcParams['font.size'] = 18


def parse_float_list(raw_value, default_values):
    if raw_value is None or raw_value.strip() == '':
        return list(default_values)
    return [float(item.strip()) for item in raw_value.split(',') if item.strip()]


def parse_ss_ranges(raw_value):
    if raw_value is None or raw_value.strip() == '':
        return [list(item) for item in DEFAULT_SS_RANGES]

    ranges = []
    for item in raw_value.split(','):
        item = item.strip()
        if not item:
            continue
        start_text, end_text = item.split('-', 1)
        ranges.append([int(start_text), int(end_text)])
    return ranges


def resolve_path(path_value, default_path):
    if path_value is None or path_value.strip() == '':
        return default_path
    return Path(path_value).expanduser()


def require_file(path_obj, description):
    if not path_obj.is_file():
        raise FileNotFoundError(f"Missing {description}: {path_obj}")


def load_fit_table(db_path):
    connection = sqlite3.connect(str(db_path))
    try:
        cursor = connection.cursor()
        cursor.execute('SELECT key, value FROM Dict')
        rows = cursor.fetchall()
    finally:
        connection.close()

    loaded = {}
    for key_blob, value_blob in rows:
        loaded[key_blob] = pickle.loads(value_blob)

    if b'fits' in loaded:
        print("Successfully loaded fit data using key b'fits'")
        return loaded[b'fits']
    if 'fits' in loaded:
        print("Successfully loaded fit data using key 'fits'")
        return loaded['fits']

    available_keys = list(loaded.keys())
    raise KeyError(f"'fits' key not found in {db_path}. Available keys: {available_keys}")


def save_figure(fig, output_path):
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    print(f"Saved plot: {output_path}")


def main():
    prefix = os.environ.get('hx_plot_prefix', DEFAULT_PREFIX).strip()
    state = os.environ.get('hx_plot_state', os.environ.get('protein_state', DEFAULT_STATE)).strip()
    pdb_id = os.environ.get('pdb_id', 'glpG-RKRK-79HIS').strip()

    work_dir = Path(os.environ.get('hx_plot_work_dir', '.')).expanduser().resolve()
    output_dir = resolve_path(os.environ.get('hx_plot_output_dir'), work_dir).resolve()
    results_dir = resolve_path(os.environ.get('hx_plot_results_dir'), work_dir / 'results').resolve()

    dfout_path = resolve_path(os.environ.get('hx_plot_dfout_path'), work_dir / f'{prefix}-dfout-{state}.csv').resolve()
    fitdata_path = resolve_path(os.environ.get('hx_plot_fitdata_path'), work_dir / f'{prefix}-fitdata-{state}').resolve()
    dg_path = resolve_path(os.environ.get('hx_plot_dg_path'), work_dir / f'{prefix}-dg.csv').resolve()
    resid_path = resolve_path(os.environ.get('hx_plot_resid_path'), results_dir / f'{pdb_id}.resid').resolve()

    sequence = os.environ.get('hx_plot_sequence', DEFAULT_SEQUENCE).strip()
    ss_ranges = parse_ss_ranges(os.environ.get('hx_plot_ss_ranges'))
    pD_corr = float(os.environ.get('hx_plot_pd_corr', str(DEFAULT_PD_CORR)))
    temperature_k = float(os.environ.get('hx_plot_temperature_k', str(DEFAULT_TEMPERATURE_K)))
    residue_offset = int(os.environ.get('hx_plot_residue_offset', str(DEFAULT_RESIDUE_OFFSET)))
    dseq = int(os.environ.get('hx_plot_dseq', str(DEFAULT_DSEQ)))
    upside_temps = parse_float_list(os.environ.get('hx_plot_upside_temps'), DEFAULT_UPSIDE_TEMPS)

    require_file(dfout_path, 'HX fit summary CSV')
    require_file(fitdata_path, 'HX fit database')
    require_file(dg_path, 'simulation dG CSV')
    require_file(resid_path, 'residue ID file')
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f'Loading HX plot inputs from {work_dir}')
    print(f'Using dfout file: {dfout_path}')
    print(f'Using fit database: {fitdata_path}')
    print(f'Using dG file: {dg_path}')
    print(f'Using resid file: {resid_path}')

    kchem = hxf.getkchem(sequence, pD_corr, temperature_k)
    dfout = pd.read_csv(dfout_path)
    fit_table = load_fit_table(fitdata_path)

    mybkc = []
    mybkcerr = []
    for _, row in dfout.iterrows():
        label = row['lab']
        if label not in fit_table:
            print(f"Warning: Label {label} not found in fit data. Skipping.")
            mybkc.append(np.nan)
            mybkcerr.append(np.nan)
            continue

        mybkc.append(fit_table[label]['kchem'][2][0][1])
        mybkcerr.append(np.sqrt(np.diag(fit_table[label]['kchem'][2][1]))[1])

    dfout['bkc'] = mybkc
    dfout['bkcerr'] = mybkcerr
    dfout['dgs'] = 0.001987 * temperature_k * np.log(dfout['pf'])
    dfout['dgerrs'] = abs(0.001987 * temperature_k * dfout['pferr'] / dfout['pf'])

    dcut2 = 1.1
    deltaD = 0.4
    df1 = dfout[(dfout['maxD'] >= deltaD) & (dfout['deltaD'] >= deltaD) & (dfout['lastD'] <= dcut2)]

    dGups = np.loadtxt(dg_path, delimiter=',')
    if dGups.ndim == 1:
        dGups = dGups[np.newaxis, :]
    print(f"Simulation data loaded. Shape: {dGups.shape}")

    resis = np.loadtxt(resid_path)
    print(f"Loaded residues from {resid_path}")

    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = cm.plasma_r
    norm = mcolors.Normalize(vmin=0, vmax=5800)

    for ss_start, ss_end in ss_ranges:
        bar = patches.Rectangle(
            (ss_start - 0.5 + residue_offset, 0),
            ss_end - ss_start,
            20,
            facecolor='#f3f3f3',
            edgecolor='#a6a6a6',
            linewidth=0.5,
            linestyle='-',
            zorder=0,
        )
        ax.add_patch(bar)

    for index, (temp, dgs) in enumerate(zip(upside_temps, dGups)):
        ax.plot(
            resis + residue_offset + 1,
            dgs,
            'o-',
            label='T = {:.2f}'.format(temp),
            ms=6,
            markerfacecolor=mycolors[index % len(mycolors)],
            markeredgecolor='k',
            mew=0.3,
            lw=0.5,
        )

    for _, row in df1.iterrows():
        dg_value = row['dgs']
        if pd.isna(dg_value):
            continue
        try:
            x1, x2 = list(map(int, row['lab'].split('-')))
        except ValueError:
            continue
        y1 = x1 - 0.5 + residue_offset + 2
        y2 = x2 + dseq + 0.5 + residue_offset
        ax.hlines(dg_value, y1, y2, colors='k', linewidth=3)

    ax.set_ylim(0, 17)
    ax.set_xlim(residue_offset - 5, len(sequence) + residue_offset + 5)
    exp_patch = patches.Patch(color='k', lw=1, label=f'{state} exp')
    handles, labels = ax.get_legend_handles_labels()
    handles.append(exp_patch)
    labels.append(exp_patch.get_label())
    ax.set_xlabel('residue')
    ax.set_ylabel('$\\Delta$G$_{HX,ave}$ (kcal/mol)')
    ax.legend(handles=handles, labels=labels, handleheight=0.1, handlelength=1, bbox_to_anchor=(1.05, 1), loc='upper left')
    save_figure(fig, output_dir / 'hx_overview_plot.png')

    valid_indices = resis.astype(int)
    valid_indices = valid_indices[(valid_indices >= 0) & (valid_indices < len(kchem))]
    kchem0 = kchem[valid_indices]

    if len(kchem0) < len(resis):
        resis = resis[:len(kchem0)]
        dGups = dGups[:, :len(kchem0)]

    mythxs = 10.0 ** np.arange(-2, 7, 0.2)
    DGpred = []

    print('Calculating predictions...')
    for _, row in df1.iterrows():
        try:
            x1, x2 = list(map(int, row['lab'].split('-')))
            mask = (resis >= x1 + 2 + dseq) & (resis < x2 + dseq + 1)
            nonzero_ks = kchem0[mask]
            if len(nonzero_ks) == 0:
                continue

            smin = int(np.floor(np.log10(1 / nonzero_ks.max()))) - 0.5
            smax = int(np.ceil(np.log10(1 / nonzero_ks.min()))) + 0.5
            tkcs = 10.0 ** np.arange(smin, smax, 0.2)

            kchem_curve = np.mean(1 - np.exp(-np.outer(nonzero_ks, tkcs)), axis=0)
            skc, bkc = hxf.fitstretchsig(kchem_curve, tkcs, True)[0]

            dgpred = []
            for temp, dgs in zip(upside_temps, dGups):
                rt = 0.001987 * temp / 0.85 * 298.0
                dgnonzero = dgs[mask]
                with np.errstate(over='ignore'):
                    kobs1 = nonzero_ks / (np.exp(dgnonzero / rt) + 1)

                pred_curve = np.mean(1 - np.exp(-np.outer(kobs1, mythxs)), axis=0)
                sdg, bdg = hxf.fitstretchsig(pred_curve, mythxs, True)[0]
                mypf = hxf.pfave(skc, bkc, sdg, bdg)
                dgpred.append([x1, x2, rt * np.log(mypf)])

            DGpred.append(dgpred)
        except Exception as exc:
            print(f"Error processing segment {row['lab']}: {exc}")

    if len(DGpred) == 0:
        print('Warning: No predictions generated. Check segment alignments or kchem data.')
        return

    DGpred = np.transpose(np.asarray(DGpred), axes=(1, 0, 2))

    fig, ax = plt.subplots(figsize=(6, 6))
    for index, (temp, dg_values) in enumerate(zip(upside_temps, DGpred[:, :, 2])):
        min_len = min(len(df1['dgs']), len(dg_values))
        ax.plot(
            df1['dgs'][:min_len],
            dg_values[:min_len],
            'o',
            ms=6,
            markerfacecolor='none',
            mew=1.5,
            color=mycolors[index % len(mycolors)],
            label='T = {:.2f}'.format(temp),
        )

    ax.plot([-2.5, 17], [-2.5, 17], '--', color='k', zorder=0)
    ax.hlines(0, -2.5, 17, color='gray', ls='-')
    ax.vlines(0, -2.5, 17, color='gray', ls='-')
    ax.set_aspect(1)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_ylim(-2.5, 17)
    ax.set_xlim(-2.5, 17)
    ax.set_yticks([0, 5, 10, 15])
    ax.set_xlabel('$\\Delta$G$_{HX,ave}^{exp}$ (kcal/mol)')
    ax.set_ylabel('$\\Delta$G$_{HX,ave}^{pred}$ (kcal/mol)')
    save_figure(fig, output_dir / 'hx_correlation_plot.png')

    fig, ax = plt.subplots(figsize=(8, 5))
    for ss_start, ss_end in ss_ranges:
        bar = patches.Rectangle(
            (ss_start - 0.5 + residue_offset, 0),
            ss_end - ss_start,
            20,
            facecolor='#f3f3f3',
            edgecolor='#a6a6a6',
            linewidth=0.5,
            linestyle='-',
            zorder=0,
        )
        ax.add_patch(bar)

    legend_patches = []
    for index, (temp, dgpred_for_temp) in enumerate(zip(upside_temps, DGpred)):
        for segment in dgpred_for_temp:
            x1, x2, dgu = segment
            y1 = x1 - 0.5 + residue_offset + 2
            y2 = x2 + dseq + 0.5 + residue_offset
            ax.hlines(dgu, y1, y2, colors=mycolors[index % len(mycolors)], linewidth=3)
        legend_patches.append(patches.Patch(color=mycolors[index % len(mycolors)], lw=1, label='T = {:.2f}'.format(temp)))

    for _, row in df1.iterrows():
        dg_value = row['dgs']
        if pd.isna(dg_value):
            continue
        try:
            x1, x2 = list(map(int, row['lab'].split('-')))
        except ValueError:
            continue
        y1 = x1 - 0.5 + residue_offset + 2
        y2 = x2 + dseq + 0.5 + residue_offset
        ax.hlines(dg_value, y1, y2, colors='k', linewidth=3)

    ax.set_ylim(0, 17)
    ax.set_xlim(residue_offset - 5, len(sequence) + residue_offset + 5)
    ax.set_xlabel('residue')
    ax.set_ylabel('$\\Delta$G$_{HX,ave}$ (kcal/mol)')
    legend_patches.append(patches.Patch(color='k', lw=1, label=f'{state} exp'))
    ax.legend(
        handles=legend_patches,
        labels=[patch.get_label() for patch in legend_patches],
        handleheight=0.1,
        handlelength=1,
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
    )
    save_figure(fig, output_dir / 'hx_prediction_segments.png')


if __name__ == '__main__':
    main()
