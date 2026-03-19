import os

import matplotlib
if 'MPLBACKEND' not in os.environ:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

from helpers.advanced_analysis_utils import build_paths, csv_list_env, load_csv_rows, load_optional_numeric_csv


default_pdb = os.environ.get('pdb_id', 'Pab1_RRM1')  # CHECKME
pdb_ids = csv_list_env('pdb_ids', [default_pdb])  # CHECKME
sim_id = os.environ.get('sim_id', 'REMD')  # CHECKME

paths = build_paths(default_pdb, sim_id)
fig_dir = paths['result_dir']
pdb_dir = paths['pdb_dir']


def load_profile_rows(pdb_id):
    profile_rows = load_csv_rows('{}/{}.csv'.format(fig_dir, pdb_id))
    ranked_rows = load_csv_rows('{}/{}_ranked.csv'.format(fig_dir, pdb_id))

    residue_ids = np.asarray(profile_rows[0], dtype=float)
    sim_profile = np.asarray(profile_rows[1], dtype=float)
    sim_ranked = np.asarray(ranked_rows[0], dtype=float)
    exp_ranked = np.asarray(ranked_rows[1], dtype=float) if len(ranked_rows) > 1 else None
    return residue_ids, sim_profile, sim_ranked, exp_ranked


loaded = []
for pdb_id in pdb_ids:
    try:
        loaded.append((pdb_id, *load_profile_rows(pdb_id)))
    except FileNotFoundError:
        print('WARNING: Missing exported dG summary for {}. Skipping.'.format(pdb_id))

if len(loaded) == 0:
    raise SystemExit('No exported dG summary CSVs were found in results/.')

fig = plt.figure()
ax = plt.subplot(111)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
for pdb_id, residue_ids, sim_profile, sim_ranked, exp_ranked in loaded:
    ax.plot(sim_ranked, label='{}'.format(pdb_id))
ax.set_xlabel('Residue (ranked by ∆G)')
ax.set_ylabel('∆G (kcal/mol)')
ax.set_title('Ranked ∆G Across Analyses')
fontP = FontProperties()
fontP.set_size('x-small')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop=fontP, columnspacing=0.8)
plt.savefig('{}/DG_ranked_summary.png'.format(fig_dir), dpi=300, bbox_inches='tight')
plt.close()

if any(exp_ranked is not None for _, _, _, _, exp_ranked in loaded):
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    for pdb_id, residue_ids, sim_profile, sim_ranked, exp_ranked in loaded:
        ax.plot(sim_ranked, label='{} Sim'.format(pdb_id))
        if exp_ranked is not None:
            ax.plot(exp_ranked, linestyle='dashed', label='{} Exp'.format(pdb_id))
    ax.set_xlabel('Residue (ranked by ∆G)')
    ax.set_ylabel('∆G (kcal/mol)')
    ax.set_title('Ranked ∆G: Simulation vs Experimental Reference')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop=fontP, columnspacing=0.8)
    plt.savefig('{}/DG_ranked_summary_with_reference.png'.format(fig_dir), dpi=300, bbox_inches='tight')
    plt.close()

fig = plt.figure()
ax = plt.subplot(111)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
for pdb_id, residue_ids, sim_profile, sim_ranked, exp_ranked in loaded:
    ax.plot(residue_ids, sim_profile, '-o', label='{}'.format(pdb_id))
ax.set_xlabel('Residue Number')
ax.set_ylabel('∆G (kcal/mol)')
ax.set_title('Residue-Level ∆G Across Analyses')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop=fontP, columnspacing=0.8)
plt.savefig('{}/DG_residue_summary.png'.format(fig_dir), dpi=300, bbox_inches='tight')
plt.close()


PASTA_IDS = {
    "A0A1G7LVP6.1_499-560",
    "A0A1G4FHL9.1_670-724",
    "A0A1H1G7N6.1_452-512",
    "A0A200HJR9.1_427-487",
    "A0A1I3SSH5.1_420-48",
    "A0A2A9E9G9.1_632-690",
    "A0A1A9GI13.1_689-74",
    "J4UDX3.1_675-727",
}
COLD_SHOCK_IDS = {
    "A0A0N9HSA5.1_3-64",
    "A0A1F3JSH7.1_3-62",
}


def classify_protein_family(pdb_id):
    if "HHH" in pdb_id:
        return "HHH"
    if "HEEH" in pdb_id:
        return "HEEH"
    if "EEHEE" in pdb_id:
        return "EEHEE"
    if "EHEE" in pdb_id:
        return "EHEE"
    if "PDB" in pdb_id:
        return "PDB"
    if pdb_id in COLD_SHOCK_IDS:
        return "Cold-Shock"
    if "R6" in pdb_id or pdb_id in PASTA_IDS:
        return "PASTA"
    if "LysM" in pdb_id:
        return "LysM"
    return "LysM"


def family_label(family_key):
    labels = {
        "HHH": "ααα",
        "HEEH": "αßßα",
        "EEHEE": "ßßαßß",
        "EHEE": "ßαßß",
        "PDB": "PDB",
        "Cold-Shock": "Cold-Shock",
        "PASTA": "ßαßßß",
        "LysM": "ßααß",
    }
    return labels.get(family_key, family_key)


def plot_family_rankorder(family_key, family_entries):
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']

    for pdb_id, residue_ids, sim_profile, sim_ranked, exp_ranked in family_entries:
        axes[0].plot(sim_ranked, label='{}'.format(pdb_id))
        if exp_ranked is not None:
            axes[1].plot(exp_ranked, label='{}'.format(pdb_id))

    axes[0].set_xlabel('Residue (ranked by ∆G)')
    axes[0].set_ylabel('∆G (kcal/mol)')
    axes[0].set_title('∆G Ranked for Upside Simulations ({})'.format(family_label(family_key)))

    axes[1].set_xlabel('Residue (ranked by ∆G)')
    axes[1].set_ylabel('∆G (kcal/mol)')
    axes[1].set_title('∆G Ranked for Experimental Reference ({})'.format(family_label(family_key)))

    fontP = FontProperties()
    fontP.set_size('x-small')

    for ax in axes:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height * 0.4])
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=max(1, min(5, len(family_entries))), prop=fontP, columnspacing=0.8)

    plt.tight_layout()
    plt.savefig('{}/DG_res_rankorder_{}.png'.format(fig_dir, family_key), dpi=300, bbox_inches='tight')
    plt.close()


def plot_family_residues(family_key, family_entries):
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']

    for pdb_id, residue_ids, sim_profile, sim_ranked, exp_ranked in family_entries:
        ax.plot(residue_ids, sim_profile, '-o', label='{}'.format(pdb_id))

    ax.set_xlabel('Residue Number (from N-term)')
    ax.set_ylabel('∆G (kcal/mol)')
    ax.set_title('∆G for Upside Simulations ({})'.format(family_label(family_key)))

    fontP = FontProperties()
    fontP.set_size('x-small')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop=fontP, columnspacing=0.8)

    plt.savefig('{}/DG_res_{}.png'.format(fig_dir, family_key), dpi=300, bbox_inches='tight')
    plt.close()


families = {}
for entry in loaded:
    families.setdefault(classify_protein_family(entry[0]), []).append(entry)

for family_key, family_entries in families.items():
    if len(family_entries) == 0:
        continue
    plot_family_rankorder(family_key, family_entries)
    plot_family_residues(family_key, family_entries)

temp_summary = load_optional_numeric_csv('{}/temp_random100.csv'.format(pdb_dir), skip_header=1, delimiter=',')
if temp_summary is not None and getattr(temp_summary, 'ndim', 0) == 2 and temp_summary.shape[1] >= 2:
    family_temps = {family_key: [] for family_key in families}
    for row in temp_summary:
        pdb_name = str(row[0]) if not isinstance(row[0], str) else row[0]
        temp_value = float(row[1])
        family_key = classify_protein_family(pdb_name)
        if family_key in family_temps:
            family_temps[family_key].append(temp_value)

    fig = plt.figure()
    ax = plt.subplot(111)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    for family_key, temperatures in family_temps.items():
        if len(temperatures) == 0:
            continue
        ax.plot(np.asarray(temperatures, dtype=float), '-o', label=family_key)

    if any(len(temperatures) > 0 for temperatures in family_temps.values()):
        ax.set_xlabel('Protein Family')
        ax.set_ylabel('T (K)')
        ax.set_title('Upside Temperatures Across Protein Families')

        fontP = FontProperties()
        fontP.set_size('x-small')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop=fontP, columnspacing=0.8)
        plt.savefig('{}/temp_Ups_100.png'.format(fig_dir), dpi=300, bbox_inches='tight')
    plt.close()
