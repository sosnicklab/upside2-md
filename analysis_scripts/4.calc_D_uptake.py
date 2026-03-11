import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib.font_manager import FontProperties
import matplotlib.cm as cm
from matplotlib.ticker import LogLocator
import pymbar  # for MBAR analysis
from pymbar import timeseries  # for timeseries analysis
import os
import csv
import pandas as pd
import sys
import matplotlib.backends.backend_pdf as pdf
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

pdb_id = 'glpG-RKRK-79HIS'  # CHECKME
# pdb_id = os.environ['pdb_id']
sim_id = 'memb_test'
# sim_id = os.environ['sim_id']  # CHECKME
start_frame = 100

main_pdb = 'glpG-RKRK-79HIS'  # CHECKME
HXMS_pdb = 'glpG-RKRK-79HIS'  # CHECKME
# HXMS_method = 'stretch_exp'  # CHECKME
HXMS_method = 'D_norm_uptake'
protein_state = 'pd9'
exp_data_file = 'GlpG psWT Sub final peptides up sum 11192024.csv'

work_dir = './'
n_rep = 48  # CHECKME
# n_rep = os.environ['n_rep']  # CHECKME

input_dir = '{}/inputs'.format(work_dir)
output_dir = '{}/outputs'.format(work_dir)
result_dir = '{}/results'.format(work_dir)
run_dir = '{}/{}'.format(output_dir, sim_id)
fig_dir = 'results'
pdb_dir = '{}/pdb'.format(work_dir)

os.makedirs(fig_dir, exist_ok=True)
os.makedirs(result_dir, exist_ok=True)

# ============================================
# set kchem params
# ============================================

mol_type = 'poly'  # CHECKME

aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

Ea_a = 14
Ea_b = 17
Ea_w = 19

Ea_asp = 1
Ea_glu = 1.083
Ea_his = 7.5

single_T = True  # CHECKME
legacy_T_range = np.array([0.85], dtype=float)
# legacy_T_range = np.array([0.84, 0.85, 0.92, 0.99, 1.1], dtype=float)

T_search_range = (0.7, 1.2)

sim_T_data = np.array([0.7, 0.71, 0.73, 0.74, 0.76, 0.77, 0.79, 0.8, 0.82, 0.84, 0.85, 0.87, 0.89, 0.9, 0.92, 0.94, 0.96, 0.98, 1, 1.02, 1.04, 1.06, 1.08, 1.2])
real_K_data = np.array([245.411765, 248.917647, 255.929412, 259.435294, 266.447059, 269.952941, 276.964706, 280.470588, 287.482353, 294.494118, 298, 305.011765, 312.023529, 315.529412, 322.541176, 329.552941, 336.564706, 343.576471, 350.588235, 357.6, 364.611765, 371.623529, 378.635294, 420.705882])
sim_T_to_real_K = interp1d(sim_T_data, real_K_data, kind='linear', fill_value='extrapolate')

pKD2O_temps = np.array([293, 298, 323, 348, 373], dtype=float)
pKD2O_values = np.array([15.05, 14.96, 14.18, 13.57, 13.11], dtype=float)
real_K_to_pKD2O = interp1d(pKD2O_temps, pKD2O_values, kind='linear', fill_value='extrapolate')

T_ref_K = 293.0
T_acid_K = 278.0

# for calculation of protein HX in D2O
pD_corr = 6.70  # CHECKME pD_read+0.4
legacy_Kws = np.array([14.96], dtype=float)
D_plus = np.power(10, -pD_corr)
legacy_od_minus_map = {
    float(T_target): np.power(10, pD_corr - Kw)
    for T_target, Kw in zip(legacy_T_range, legacy_Kws)
}

ka_poly = (np.power(10, 1.62)) / 60
kb_poly = (np.power(10, 10.18)) / 60
kw_poly = (np.power(10, -1.5)) / 60

ka_oligo = ka_poly * 2.34
kb_oligo = kb_poly * 1.35
kw_oligo = kw_poly * 1.585

if mol_type == 'poly':
    ka = ka_poly
    kb = kb_poly
    kw = kw_poly
elif mol_type == 'oligo':
    ka = ka_oligo
    kb = kb_oligo
    kw = kw_oligo
else:
    sys.exit()

# ============================================
# load data
# ============================================

Pot = []
Rg = []
Rmsd = []
Hb = []
Ts = []
T = []
PS = []
for i in range(n_rep):
    Pot.append(np.load('{}/{}_{}_{}_Energy.npy'.format(result_dir, pdb_id, sim_id, i))[:, 0])
    Rg.append(np.load('{}/{}_{}_{}_Rg.npy'.format(result_dir, pdb_id, sim_id, i)))
    Hb.append(np.load('{}/{}_{}_{}_Hbond.npy'.format(result_dir, pdb_id, sim_id, i)))
    Rmsd.append(np.load('{}/{}_{}_{}_Rmsd.npy'.format(result_dir, pdb_id, sim_id, i)))
    PS.append(np.load('{}/{}_{}_{}_PS.npy'.format(result_dir, pdb_id, sim_id, i)))

    t = np.load('{}/{}_{}_{}_T.npy'.format(result_dir, pdb_id, sim_id, i))
    nsize = Pot[-1].size
    Ts.append(np.zeros(nsize) + t)
    T.append(t)

min_idx = np.min([len(Rmsd[i]) for i in range(n_rep)])

Rmsd_idx = [Rmsd[i][:min_idx] for i in range(n_rep)]
Pot_idx = [Pot[i][:min_idx] for i in range(n_rep)]
Rg_idx = [Rg[i][:min_idx] for i in range(n_rep)]
Hb_idx = [Hb[i][:min_idx] for i in range(n_rep)]
Ts_idx = [Ts[i][:min_idx] for i in range(n_rep)]
PS_idx = [PS[i][:min_idx] for i in range(n_rep)]

Rmsd = np.array(Rmsd_idx)
Pot = np.array(Pot_idx)
Rg = np.array(Rg_idx)
Hb = np.array(Hb_idx)
Ts = np.array(Ts_idx)
T = np.array(T)
PS = np.array(PS_idx)

res = np.loadtxt('{}/{}.resid'.format(result_dir, pdb_id), dtype=int)
n_res = res.size
res_to_idx = {int(residue): idx for idx, residue in enumerate(res)}

fasta = np.loadtxt('{}/{}.fasta'.format(input_dir, pdb_id), skiprows=1, dtype=str)
fasta_main = fasta

if np.size(fasta) == 1:
    sequence = fasta.tolist()
    seq = fasta.reshape(1)[0]
else:
    sequence = ''.join(fasta)
    seq = np.array(sequence).reshape(1)[0]

if np.size(fasta_main) == 1:
    main_sequence = fasta_main.tolist()
    main_seq = fasta_main.reshape(1)[0]
else:
    main_sequence = ''.join(fasta_main)
    main_seq = np.array(main_sequence).reshape(1)[0]


def build_time_grid_from_log_times(log_times):
    log_times = np.asarray(log_times, dtype=float)
    if log_times.ndim > 1:
        log_times = log_times[0]

    log_times = log_times[~np.isnan(log_times)]
    log_times = np.unique(log_times)
    log_times.sort()

    if log_times.size == 0:
        return np.array([])
    if log_times.size == 1:
        return np.power(10, log_times)

    step = float(np.mean(np.diff(log_times)))
    return np.asarray([10 ** s for s in np.arange(np.min(log_times), np.max(log_times) + step, step)])


def build_default_time_grid(log_start=-1.0, log_stop=5.0, log_step=0.02):
    return np.asarray([10 ** s for s in np.arange(log_start, log_stop, log_step)])


def sanitize_peptides(peptide_starts, peptide_ends, peptide_labels, sequence_str, source_name):
    starts = np.asarray(peptide_starts)
    ends = np.asarray(peptide_ends)
    labels = np.asarray(peptide_labels, dtype=str)

    valid_starts = []
    valid_ends = []
    valid_labels = []

    for start, end, label in zip(starts, ends, labels):
        try:
            start = int(start)
            end = int(end)
        except (TypeError, ValueError):
            print('WARNING: Skipping {} peptide with non-integer bounds: {} {}'.format(source_name, start, end))
            continue

        if start < 1 or end < start or end > len(sequence_str):
            print('WARNING: Skipping {} peptide outside sequence bounds: {}-{}'.format(source_name, start, end))
            continue

        expected_label = sequence_str[start - 1:end]
        label = label.strip()
        if label and label != expected_label:
            print(
                'WARNING: {} peptide {}-{} sequence mismatch ({} != {}). Using residue bounds.'.format(
                    source_name,
                    start,
                    end,
                    label,
                    expected_label,
                )
            )

        valid_starts.append(start)
        valid_ends.append(end)
        valid_labels.append(label if label else expected_label)

    return (
        np.asarray(valid_starts, dtype=int),
        np.asarray(valid_ends, dtype=int),
        np.asarray(valid_labels, dtype=str),
    )


def load_legacy_inputs():
    try:
        legacy_time_log = np.load(
            '{}/{}/{}_{}_time_peps_{}.npy'.format(output_dir, sim_id, HXMS_pdb, HXMS_method, protein_state)
        )
    except FileNotFoundError:
        print('WARNING: Legacy HXMS time grid not found. Skipping legacy workflow.')
        return np.array([]), np.array([]), np.array([]), np.array([])

    try:
        with open('{}/{}/{}_pep_ids.csv'.format(output_dir, sim_id, HXMS_pdb), 'r', encoding='utf-8-sig') as f:
            rows = list(csv.reader(f))
    except FileNotFoundError:
        print('WARNING: Legacy peptide ID file not found. Skipping legacy workflow.')
        return np.array([]), np.array([]), np.array([]), np.array([])

    if len(rows) < 3:
        print('WARNING: Legacy peptide ID file is incomplete. Skipping legacy workflow.')
        return np.array([]), np.array([]), np.array([]), np.array([])

    legacy_times = build_time_grid_from_log_times(legacy_time_log)
    legacy_starts, legacy_ends, legacy_labels = sanitize_peptides(rows[0], rows[1], rows[2], seq, 'legacy')
    return legacy_times, legacy_starts, legacy_ends, legacy_labels


def load_experimental_data():
    col_names = [
        'Protein State',
        'Protein',
        'Start',
        'End',
        'Sequence',
        'Peptide Mass',
        'RT (min)',
        'Deut Time (sec)',
        'maxD',
        'Theor Uptake #D',
        '#D',
        '%D',
        'Conf Interval (#D)',
        '#Rep',
        'Confidence',
        'Stddev',
    ]

    try:
        exp_df = pd.read_csv(
            exp_data_file,
            skiprows=3,
            header=None,
            names=col_names,
            usecols=range(16),
            encoding='utf-8-sig',
        )
    except FileNotFoundError:
        print('WARNING: Experimental data file not found: {}'.format(exp_data_file))
        return pd.DataFrame()
    except Exception as exc:
        print('WARNING: Failed to load experimental data: {}'.format(exc))
        return pd.DataFrame()

    exp_df.columns = exp_df.columns.str.strip()
    exp_df['Protein State'] = exp_df['Protein State'].astype(str).str.strip()
    exp_df = exp_df[exp_df['Protein State'] == protein_state].copy()
    exp_df['Deut Time (sec)'] = pd.to_numeric(exp_df['Deut Time (sec)'], errors='coerce')
    exp_df['%D'] = pd.to_numeric(exp_df['%D'], errors='coerce')
    exp_df['Start'] = pd.to_numeric(exp_df['Start'], errors='coerce')
    exp_df['End'] = pd.to_numeric(exp_df['End'], errors='coerce')
    exp_df = exp_df.dropna(subset=['Deut Time (sec)', '%D', 'Start', 'End'])

    if exp_df.empty:
        print("WARNING: Experimental data loaded, but no valid rows matched protein state '{}'.".format(protein_state))
        return pd.DataFrame()

    print('Successfully loaded experimental data from {}.'.format(exp_data_file))
    print("Found {} unique peptides for '{}'.".format(len(exp_df[['Start', 'End']].drop_duplicates()), protein_state))
    return exp_df


def build_experimental_peptides(exp_df):
    if exp_df.empty:
        return np.array([]), np.array([]), np.array([])

    exp_peptides = exp_df[['Start', 'End', 'Sequence']].drop_duplicates().astype({'Start': int, 'End': int})
    return sanitize_peptides(
        exp_peptides['Start'].values,
        exp_peptides['End'].values,
        exp_peptides['Sequence'].values,
        seq,
        'experimental',
    )


def export_peptide_ids(output_path, peptide_starts, peptide_ends, peptide_labels):
    with open(output_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(peptide_starts)
        writer.writerow(peptide_ends)
        writer.writerow(peptide_labels)


def calculate_temperature_parameters(T_target, mode):
    if mode == 'legacy':
        real_K = float((T_target / 0.85) * 298.0)
        OD_minus = legacy_od_minus_map.get(float(T_target))
        if OD_minus is None:
            OD_minus = np.power(10, pD_corr - float(real_K_to_pKD2O(real_K)))
    elif mode == 'exp':
        real_K = float(sim_T_to_real_K(T_target))
        OD_minus = np.power(10, pD_corr - float(real_K_to_pKD2O(real_K)))
    else:
        raise ValueError('Unknown temperature mode: {}'.format(mode))

    dT_R = ((1 / real_K) - (1 / T_ref_K)) / 0.001987

    Ft_a = np.exp(-Ea_a * dT_R)
    Ft_b = np.exp(-Ea_b * dT_R)
    Ft_w = np.exp(-Ea_w * dT_R)

    pKc_asp = -1 * np.log10(np.power(10, -4.48) * np.exp(-Ea_asp * (((1 / real_K) - (1 / T_acid_K)) / 0.001987)))
    pKc_glu = -1 * np.log10(np.power(10, -4.93) * np.exp(-Ea_glu * (((1 / real_K) - (1 / T_acid_K)) / 0.001987)))
    pKc_his = -1 * np.log10(np.power(10, -7.42) * np.exp(-Ea_his * (((1 / real_K) - (1 / T_acid_K)) / 0.001987)))

    return Ft_a, Ft_b, Ft_w, OD_minus, pKc_asp, pKc_glu, pKc_his


def normalize_uptake(actual, theoretical):
    denom_theor = np.subtract(theoretical[-1], actual[0])
    if denom_theor == 0:
        norm_theoretical = np.zeros_like(theoretical)
    else:
        norm_theoretical = np.subtract(theoretical, actual[0]) / denom_theor * 100

    with np.errstate(divide='ignore', invalid='ignore'):
        norm_actual = norm_theoretical * (actual / theoretical)
        norm_actual[~np.isfinite(norm_actual)] = 0

    return norm_actual, norm_theoretical


def is_non_decreasing(values, atol=1e-10):
    if values.size < 2:
        return True
    return bool(np.all(np.diff(values) >= -atol))


def calculate_all_hdx_data(T_target, times, peptide_starts, peptide_ends, temperature_mode):
    if len(times) == 0:
        print('WARNING: T={} skipped because no time grid is available.'.format(T_target))
        return None

    try:
        Ft_a, Ft_b, Ft_w, OD_minus, pKc_asp, pKc_glu, pKc_his = calculate_temperature_parameters(T_target, temperature_mode)
    except Exception as exc:
        print('WARNING: T={} failed during temperature parameter calculation: {}'.format(T_target, exc))
        return None

    lambda_ma = [0, -0.54, np.log10(np.power(10, -0.90 - pD_corr) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr)) + np.power(10, 0.90 - pKc_asp) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr))), np.log10(np.power(10, -0.60 - pD_corr) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr)) + np.power(10, -0.90 - pKc_glu) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr))), -0.52, -0.22, np.log10(np.power(10, -0.80 - pD_corr) / (np.power(10, -pKc_his) + np.power(10, -pD_corr)) + np.power(10, -pKc_his) / (np.power(10, -pKc_his) + np.power(10, -pD_corr))), -0.91, -0.56, -0.57, -0.64, -0.58, 0, -0.47, -0.59, -0.437992278, -0.79, -0.739022273, -0.4, -0.41]
    rho_ma = [0, -0.46, np.log10(np.power(10, -0.12 - pD_corr) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr)) + np.power(10, 0.58 - pKc_asp) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr))), np.log10(np.power(10, -0.27 - pD_corr) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr)) + np.power(10, 0.31 - pKc_glu) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr))), -0.43, 0.218176047, np.log10(np.power(10, -0.51 - pD_corr) / (np.power(10, -pKc_his) + np.power(10, -pD_corr)) + np.power(10, -pKc_his) / (np.power(10, -pKc_his) + np.power(10, -pD_corr))), -0.59, -0.29, -0.13, -0.28, -0.13, -0.194773472, -0.27, -0.32, -0.388518935, -0.468073126, -0.3, -0.44, -0.37]
    lambda_mb = [0, 0.62, np.log10(np.power(10, 0.69 - pD_corr) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr)) + np.power(10, 0.10 - pKc_asp) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr))), np.log10(np.power(10, 0.24 - pD_corr) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr)) + np.power(10, -0.11 - pKc_glu) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr))), -0.235859464, -0.03, np.log10(np.power(10, 0.80 - pD_corr) / (np.power(10, -pKc_his) + np.power(10, -pD_corr)) + np.power(10, -0.10 - pKc_his) / (np.power(10, -pKc_his) + np.power(10, -pD_corr))), -0.73, -0.04, -0.576252728, -0.008954843, 0.49, 0, 0.06, 0.076712254, 0.37, -0.06625798, -0.701934483, -0.41, -0.27]
    rho_mb = [0, 0.55, np.log10(np.power(10, 0.60 - pD_corr) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr)) + np.power(10, -0.18 - pKc_asp) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr))), np.log10(np.power(10, 0.39 - pD_corr) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr)) + np.power(10, -0.15 - pKc_glu) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr))), 0.063131587, 0.17, np.log10(np.power(10, 0.83 - pD_corr) / (np.power(10, -pKc_his) + np.power(10, -pD_corr)) + np.power(10, 0.14 - pKc_his) / (np.power(10, -pKc_his) + np.power(10, -pD_corr))), -0.23, 0.12, -0.21, 0.11, 0.32, -0.24, 0.2, 0.22, 0.299550286, 0.2, -0.14, -0.11, 0.05]
    p_df = pd.DataFrame(
        {
            'aa': aa,
            'lambda_ma': lambda_ma,
            'rho_ma': rho_ma,
            'lambda_mb': lambda_mb,
            'rho_mb': rho_mb,
        }
    )

    Fa_a = []
    Fb_b = []
    Fb_w = []
    for i, s in enumerate(sequence):
        if s == 'P':
            Fa = 0
            Fb = 0
        elif i == 0:
            Fa = 0
            Fb = 0
        elif i == 1:
            Fa = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_ma'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_ma'].iloc[0] - 1.32)
            Fb = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_mb'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_mb'].iloc[0] + 1.62)
        elif i == len(sequence) - 1:
            Fa = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_ma'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_ma'].iloc[0] + 0.95)
            Fb = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_mb'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_mb'].iloc[0] - 1.80)
        else:
            Fa = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_ma'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_ma'].iloc[0])
            Fb = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_mb'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_mb'].iloc[0])

        Fa_a.append(Fa * D_plus * ka * Ft_a)
        Fb_b.append(Fb * OD_minus * kb * Ft_b)
        Fb_w.append(Fb * kw * Ft_w)

    Fa_a = np.asarray(Fa_a)
    Fb_b = np.asarray(Fb_b)
    Fb_w = np.asarray(Fb_w)
    k_chem = Fa_a + Fb_b + Fb_w

    k_chem_Up = k_chem[res]
    u_n = (cE0 / (T_target * kB)).flatten()
    log_w1 = mbar0._computeUnnormalizedLogWeights(u_n)
    w1 = np.exp(log_w1)
    w1 /= np.sum(w1)

    k_obs = np.zeros(n_res)
    mean_pf = np.zeros(n_res)
    for r in range(n_res):
        pf_i = PS[:, start_frame:, r].flatten()
        mean_pf[r] = np.average(pf_i, weights=w1)
        if mean_pf[r] == 1:
            k_obs[r] = k_chem_Up[r] / 1000.0
        else:
            k_obs[r] = k_chem_Up[r] * (1 - mean_pf[r])

    if np.min(np.subtract(k_chem_Up, k_obs)) < 0:
        print('WARNING: T={} skipped because k_obs exceeded k_chem.'.format(T_target))
        return None

    D = []
    D_theor = []
    for residue_index in range(len(sequence)):
        if residue_index in res_to_idx:
            idx = res_to_idx[residue_index]
            D.append(1 - np.exp(-k_obs[idx] * times))
            D_theor.append(1 - np.exp(-k_chem_Up[idx] * times))
        else:
            zero_curve = 1 - np.exp(0 * times)
            D.append(zero_curve)
            D_theor.append(zero_curve)

    D = np.asarray(D)
    D_theor = np.asarray(D_theor)

    if len(peptide_starts) > 0:
        D_peps = []
        D_peps_theors = []
        D_norm_peps = []
        D_norm_peps_theors = []
        for start, end in zip(peptide_starts, peptide_ends):
            start_idx = start - 1
            end_idx = end

            D_pep = np.sum(D[start_idx:end_idx], axis=0)
            D_pep_theor = np.sum(D_theor[start_idx:end_idx], axis=0)
            D_norm_pep, D_norm_pep_theor = normalize_uptake(D_pep, D_pep_theor)

            if not is_non_decreasing(D_pep) or not is_non_decreasing(D_pep_theor):
                print('WARNING: T={} skipped because peptide uptake was not monotonic.'.format(T_target))
                return None
            if not is_non_decreasing(D_norm_pep) or not is_non_decreasing(D_norm_pep_theor):
                print('WARNING: T={} skipped because normalized peptide uptake was not monotonic.'.format(T_target))
                return None

            D_peps.append(D_pep)
            D_peps_theors.append(D_pep_theor)
            D_norm_peps.append(D_norm_pep)
            D_norm_peps_theors.append(D_norm_pep_theor)

        D_peps = np.asarray(D_peps)
        D_peps_theors = np.asarray(D_peps_theors)
        D_norm_peps = np.asarray(D_norm_peps)
        D_norm_peps_theors = np.asarray(D_norm_peps_theors)

        if np.nanmin(np.subtract(D_peps_theors, D_peps)) < 0:
            print('WARNING: T={} skipped because peptide uptake exceeded theoretical uptake.'.format(T_target))
            return None
        if np.nanmin(np.subtract(D_norm_peps_theors, D_norm_peps)) < 0:
            print('WARNING: T={} skipped because normalized peptide uptake exceeded theoretical uptake.'.format(T_target))
            return None
    else:
        empty_shape = (0, len(times))
        D_peps = np.empty(empty_shape)
        D_peps_theors = np.empty(empty_shape)
        D_norm_peps = np.empty(empty_shape)
        D_norm_peps_theors = np.empty(empty_shape)

    D_res = np.sum(D, axis=0)
    D_res_theor = np.sum(D_theor, axis=0)
    D_norm_res, D_norm_res_theor = normalize_uptake(D_res, D_res_theor)

    if not is_non_decreasing(D_res) or not is_non_decreasing(D_res_theor):
        print('WARNING: T={} skipped because full-protein uptake was not monotonic.'.format(T_target))
        return None
    if not is_non_decreasing(D_norm_res) or not is_non_decreasing(D_norm_res_theor):
        print('WARNING: T={} skipped because normalized full-protein uptake was not monotonic.'.format(T_target))
        return None
    if np.nanmin(np.subtract(D_res_theor, D_res)) < 0:
        print('WARNING: T={} skipped because full-protein uptake exceeded theoretical uptake.'.format(T_target))
        return None
    if np.nanmin(np.subtract(D_norm_res_theor, D_norm_res)) < 0:
        print('WARNING: T={} skipped because normalized full-protein uptake exceeded theoretical uptake.'.format(T_target))
        return None

    return {
        'D_peps': D_peps,
        'D_peps_theors': D_peps_theors,
        'D_norm_peps': D_norm_peps,
        'D_norm_peps_theors': D_norm_peps_theors,
        'D_res': D_res,
        'D_res_theor': D_res_theor,
        'D_norm_res': D_norm_res,
        'D_norm_res_theor': D_norm_res_theor,
    }


def collect_temperature_slices(T_targets, times, peptide_starts, peptide_ends, temperature_mode):
    collected = {
        'D_peps': [],
        'D_peps_theors': [],
        'D_norm_peps': [],
        'D_norm_peps_theors': [],
        'D_res': [],
        'D_res_theor': [],
        'D_norm_res': [],
        'D_norm_res_theor': [],
    }
    valid_T_targets = []

    for T_target in np.asarray(T_targets, dtype=float):
        hdx_data = calculate_all_hdx_data(T_target, times, peptide_starts, peptide_ends, temperature_mode)
        if hdx_data is None:
            continue

        valid_T_targets.append(float(T_target))
        for key in collected:
            collected[key].append(hdx_data[key])

    for key in collected:
        collected[key] = np.asarray(collected[key])

    return np.asarray(valid_T_targets, dtype=float), collected


def build_colors(num_colors):
    if num_colors <= 0:
        return []
    if num_colors == 1:
        return cm.viridis(np.array([0.6]))
    return cm.viridis(np.linspace(0, 0.9, num_colors))


def has_peptide_plot_data(peptide_labels, uptake_plot):
    return len(peptide_labels) > 0 and getattr(uptake_plot, 'ndim', 0) >= 3 and uptake_plot.shape[0] > 0 and uptake_plot.shape[1] > 0


def has_full_protein_plot_data(uptake_plot):
    return getattr(uptake_plot, 'ndim', 0) >= 2 and uptake_plot.shape[0] > 0


def plot_peptide_pdf(output_path, times, T_targets, peptide_labels, peptide_starts, peptide_ends, uptake_plot, theor_plot, ylabel, title_mode, exp_df=None):
    pdf_pages = pdf.PdfPages(output_path)

    if not has_peptide_plot_data(peptide_labels, uptake_plot):
        print('Skipping {} because no peptide data is available.'.format(output_path))
        pdf_pages.close()
        return

    colors = build_colors(len(T_targets))
    for p, peptide_label in enumerate(peptide_labels):
        fig = plt.figure()
        ax = plt.subplot(111)
        plt.rcParams['font.family'] = 'sans-serif'

        if np.isnan(uptake_plot[:, p]).any() or np.isnan(theor_plot[:, p]).any():
            plt.close()
            continue

        for tt, (T_target, uptake, theoretical) in enumerate(zip(T_targets, uptake_plot[:, p], theor_plot[:, p])):
            color = colors[tt]
            ax.plot(times, uptake, label='T = {:.3f} Uptake'.format(T_target), color=color)
            ax.plot(times, theoretical, label='T = {:.3f} Theoretical'.format(T_target), linestyle='dashed', color=color)

        if exp_df is not None and not exp_df.empty:
            start = peptide_starts[p]
            end = peptide_ends[p]
            peptide_exp_data = exp_df[(exp_df['Start'] == start) & (exp_df['End'] == end)]
            if not peptide_exp_data.empty:
                marker = 'o' if ylabel == 'Normalized %D' else 'x'
                color = 'red' if ylabel == 'Normalized %D' else 'green'
                label = 'Experimental ({})'.format(protein_state)
                if ylabel == 'D':
                    label = 'Experimental ({}, %D)'.format(protein_state)
                ax.scatter(peptide_exp_data['Deut Time (sec)'], peptide_exp_data['%D'], label=label, color=color, zorder=10, marker=marker)

        ax.set_xscale('log')
        ax.xaxis.set_major_locator(LogLocator(base=10, numticks=10))
        ax.set_xlabel('Log Time (s)')
        ax.set_ylabel(ylabel)
        if ylabel == 'Normalized %D':
            ax.set_ylim(0, 100)

        if ylabel == 'D':
            title_prefix = 'D Uptake Prediction'
        else:
            title_prefix = '%D Uptake Prediction'

        if title_mode == 'legacy':
            ax.set_title('{} for {} at Upside T Slices'.format(title_prefix, peptide_label))
        else:
            ax.set_title('{} for res {}-{} (Optimal T)'.format(title_prefix, peptide_starts[p], peptide_ends[p]))

        plt.legend()
        plt.tight_layout()
        pdf_pages.savefig(fig, dpi=300, bbox_inches='tight')
        plt.close()

    pdf_pages.close()


def plot_full_protein_png(output_path, times, T_targets, uptake_plot, theor_plot, ylabel, title):
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.rcParams['font.family'] = 'sans-serif'

    if uptake_plot.shape[0] == 0:
        print('Skipping {} because no full-protein data is available.'.format(output_path))
        plt.close()
        return

    colors = build_colors(len(T_targets))
    for tt, (T_target, uptake, theoretical) in enumerate(zip(T_targets, uptake_plot, theor_plot)):
        color = colors[tt]
        ax.plot(times, uptake, label='T = {:.3f}'.format(T_target), color=color)
        ax.plot(times, theoretical, linestyle='dashed', color=color if len(T_targets) == 1 else 'grey')

    ax.set_xscale('log')
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=10))
    ax.set_xlabel('Log Time (s)')
    ax.set_ylabel(ylabel)
    if ylabel == 'Normalized %D':
        ax.set_ylim(0, 100)
    ax.set_title(title)

    fontP = FontProperties()
    fontP.set_size('x-small')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop=fontP, columnspacing=0.8)

    plt.savefig(output_path)
    plt.close()


def export_legacy_outputs(times, peptide_labels, peptide_starts, peptide_ends, T_targets, plot_data):
    has_peptide_data = has_peptide_plot_data(peptide_labels, plot_data['D_peps'])
    has_full_protein_data = has_full_protein_plot_data(plot_data['D_res'])

    if not has_peptide_data and not has_full_protein_data:
        print('Legacy workflow produced no output files.')
        return

    if has_peptide_data:
        export_peptide_ids('{}/{}_pep_ids.csv'.format(result_dir, pdb_id), peptide_starts, peptide_ends, peptide_labels)

        plot_peptide_pdf(
            fig_dir + '/D_uptake_{}.pdf'.format(pdb_id),
            times,
            T_targets,
            peptide_labels,
            peptide_starts,
            peptide_ends,
            plot_data['D_peps'],
            plot_data['D_peps_theors'],
            'D',
            'legacy',
        )
        plot_peptide_pdf(
            fig_dir + '/D_uptake_norm_{}.pdf'.format(pdb_id),
            times,
            T_targets,
            peptide_labels,
            peptide_starts,
            peptide_ends,
            plot_data['D_norm_peps'],
            plot_data['D_norm_peps_theors'],
            'Normalized %D',
            'legacy',
        )
    else:
        print('WARNING: No peptide definitions available. Skipping peptide-level legacy exports.')

    if has_full_protein_data:
        plot_full_protein_png(
            fig_dir + '/D_uptake.png',
            times,
            T_targets,
            plot_data['D_res'],
            plot_data['D_res_theor'],
            'D',
            'D Uptake Prediction at Upside T Slices',
        )
        plot_full_protein_png(
            fig_dir + '/D_uptake_norm.png',
            times,
            T_targets,
            plot_data['D_norm_res'],
            plot_data['D_norm_res_theor'],
            'Normalized %D',
            '%D Uptake Prediction at Upside T Slices',
        )
    else:
        print('WARNING: No full-protein legacy uptake data available.')

    if single_T and len(T_targets) == 1 and has_full_protein_data:
        if has_peptide_data:
            np.save('{}/{}_percentD_peps.npy'.format(result_dir, pdb_id), peptide_labels)
            np.save('{}/{}_percentD_time_peps.npy'.format(result_dir, pdb_id), np.log10(times))
            np.save('{}/{}_percentD_d_norm_peps.npy'.format(result_dir, pdb_id), plot_data['D_norm_peps'][0])
            np.save('{}/{}_percentD_d_norm_theor.npy'.format(result_dir, pdb_id), plot_data['D_norm_peps_theors'][0])

        with open('{}/{}_percentD.csv'.format(result_dir, pdb_id), 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(np.log10(times))
            writer.writerow(plot_data['D_norm_res'][0])
            writer.writerow(plot_data['D_norm_res_theor'][0])


def export_experimental_outputs(times, peptide_labels, peptide_starts, peptide_ends, T_targets, plot_data, exp_df):
    has_peptide_data = has_peptide_plot_data(peptide_labels, plot_data['D_peps'])
    has_full_protein_data = has_full_protein_plot_data(plot_data['D_res'])

    if not has_peptide_data and not has_full_protein_data:
        print('Experimental workflow produced no output files.')
        return

    if has_peptide_data:
        plot_peptide_pdf(
            fig_dir + '/D_uptake_{}_EXP.pdf'.format(pdb_id),
            times,
            T_targets,
            peptide_labels,
            peptide_starts,
            peptide_ends,
            plot_data['D_peps'],
            plot_data['D_peps_theors'],
            'D',
            'exp',
            exp_df=exp_df,
        )
        plot_peptide_pdf(
            fig_dir + '/D_uptake_norm_{}_EXP.pdf'.format(pdb_id),
            times,
            T_targets,
            peptide_labels,
            peptide_starts,
            peptide_ends,
            plot_data['D_norm_peps'],
            plot_data['D_norm_peps_theors'],
            'Normalized %D',
            'exp',
            exp_df=exp_df,
        )
    else:
        print('WARNING: No peptide definitions available. Skipping peptide-level experimental exports.')

    if has_full_protein_data:
        plot_full_protein_png(
            fig_dir + '/D_uptake_EXP.png',
            times,
            T_targets,
            plot_data['D_res'],
            plot_data['D_res_theor'],
            'D',
            'D Uptake Prediction at Optimal T = {:.3f}'.format(T_targets[0]),
        )
        plot_full_protein_png(
            fig_dir + '/D_uptake_norm_EXP.png',
            times,
            T_targets,
            plot_data['D_norm_res'],
            plot_data['D_norm_res_theor'],
            'Normalized %D',
            '%D Uptake Prediction at Optimal T = {:.3f}'.format(T_targets[0]),
        )
    else:
        print('WARNING: No full-protein experimental uptake data available.')

    if has_peptide_data:
        np.save('{}/{}_percentD_peps_EXP.npy'.format(result_dir, pdb_id), peptide_labels)
        np.save('{}/{}_percentD_time_peps_EXP.npy'.format(result_dir, pdb_id), np.log10(times))
        np.save('{}/{}_percentD_d_norm_peps_EXP.npy'.format(result_dir, pdb_id), plot_data['D_norm_peps'][0])
        np.save('{}/{}_percentD_d_norm_theor_EXP.npy'.format(result_dir, pdb_id), plot_data['D_norm_peps_theors'][0])

    if has_full_protein_data:
        with open('{}/{}_percentD_EXP.csv'.format(result_dir, pdb_id), 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(np.log10(times))
            writer.writerow(plot_data['D_norm_res'][0])
            writer.writerow(plot_data['D_norm_res_theor'][0])


legacy_times, legacy_pep_starts, legacy_pep_ends, legacy_pep_labels = load_legacy_inputs()
if len(legacy_times) == 0:
    legacy_times = build_default_time_grid()
    print('WARNING: Legacy HXMS time grid unavailable. Using default simulation-only time grid.')

exp_df = load_experimental_data()
exp_pep_starts, exp_pep_ends, exp_pep_labels = build_experimental_peptides(exp_df)

# ============================================
# mbar
# ============================================

kB = 1.0
T = np.array(T)
beta = kB * T ** (-1)

cE0 = Pot[:, start_frame:]

FN = cE0[0].size
FNs = np.zeros([n_rep], np.int32) + FN
reducedPot0 = np.zeros([n_rep, n_rep, FN], np.float32)
for k in range(n_rep):
    for l in range(n_rep):
        reducedPot0[k, l] = beta[l] * cE0[k]
mbar0 = pymbar.MBAR(reducedPot0, FNs, verbose=True)

# ============================================
# legacy workflow
# ============================================

legacy_valid_T_targets, legacy_plot_data = collect_temperature_slices(
    legacy_T_range,
    legacy_times,
    legacy_pep_starts,
    legacy_pep_ends,
    'legacy',
)
export_legacy_outputs(
    legacy_times,
    legacy_pep_labels,
    legacy_pep_starts,
    legacy_pep_ends,
    legacy_valid_T_targets,
    legacy_plot_data,
)

# ============================================
# experimental optimization workflow
# ============================================

exp_times = build_default_time_grid()


def calculate_total_error(T_target):
    hdx_data = calculate_all_hdx_data(T_target, exp_times, exp_pep_starts, exp_pep_ends, 'exp')
    if hdx_data is None:
        return 1e99

    total_squared_error = 0.0
    D_norm_peps = hdx_data['D_norm_peps']

    for p in range(len(exp_pep_labels)):
        if exp_pep_starts[p] <= 11:
            continue

        peptide_exp_data = exp_df[(exp_df['Start'] == exp_pep_starts[p]) & (exp_df['End'] == exp_pep_ends[p])]
        if peptide_exp_data.empty:
            continue

        try:
            sim_curve_interp = interp1d(
                exp_times,
                D_norm_peps[p],
                kind='linear',
                bounds_error=False,
                fill_value=(D_norm_peps[p][0], D_norm_peps[p][-1]),
            )
        except Exception as exc:
            print(
                'WARNING: Interpolation failed for peptide {}-{} at T={}: {}'.format(
                    exp_pep_starts[p],
                    exp_pep_ends[p],
                    T_target,
                    exc,
                )
            )
            continue

        exp_times_sec = peptide_exp_data['Deut Time (sec)'].values
        exp_uptake_pct = peptide_exp_data['%D'].values
        sim_uptake_at_exp_times = sim_curve_interp(exp_times_sec)
        total_squared_error += np.sum((sim_uptake_at_exp_times - exp_uptake_pct) ** 2)

    print('Testing T = {:.5f}, Total Squared Error = {:.2f}'.format(T_target, total_squared_error))
    return total_squared_error


if not exp_df.empty and len(exp_pep_labels) > 0:
    print('--- Starting Temperature Optimization ---')
    print('Minimizing error, ignoring peptides overlapping with res 1-11.')

    opt_result = minimize_scalar(
        calculate_total_error,
        bounds=T_search_range,
        method='bounded',
        options={'xatol': 1e-3},
    )

    T_optimal = float(opt_result.x)
    print('--- Optimization Finished ---')
    print('Optimal T found: {:.5f}'.format(T_optimal))
    print('Minimum Error: {:.2f}'.format(opt_result.fun))

    exp_valid_T_targets, exp_plot_data = collect_temperature_slices(
        np.array([T_optimal]),
        exp_times,
        exp_pep_starts,
        exp_pep_ends,
        'exp',
    )
    export_experimental_outputs(
        exp_times,
        exp_pep_labels,
        exp_pep_starts,
        exp_pep_ends,
        exp_valid_T_targets,
        exp_plot_data,
        exp_df,
    )
else:
    print('WARNING: No experimental data or peptides found. Skipping experimental optimization workflow.')

print('Script finished.')
