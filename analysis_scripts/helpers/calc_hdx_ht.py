import os

import matplotlib
if 'MPLBACKEND' not in os.environ:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pymbar
from matplotlib.font_manager import FontProperties
from scipy.stats import ks_2samp, linregress

from helpers.advanced_analysis_utils import (
    build_paths,
    float_list_env,
    load_optional_numeric_csv,
    load_replica_results,
    resolve_start_frame,
    save_csv_rows,
)


pdb_id = os.environ.get('pdb_id', 'Pab1_RRM1')  # CHECKME
sim_id = os.environ.get('sim_id', 'REMD')  # CHECKME
n_rep = int(os.environ.get('n_rep', '48'))  # CHECKME
target_temps = float_list_env('T_targets', [0.80, 0.85, 0.90])  # CHECKME
plot_target_temps = float_list_env('plot_T_targets', [0.70, 0.75, 0.80, 0.85, 0.90])
target_T = float(os.environ.get('target_T', '0.8'))  # CHECKME
m_sens = float(os.environ.get('m_sens', '0.04'))  # CHECKME
reference_temperature = float(os.environ.get('reference_temperature', '0.85'))  # CHECKME

paths = build_paths(pdb_id, sim_id)
results = load_replica_results(pdb_id, sim_id, paths['result_dir'], n_rep)

Pot = results['Pot']
Rg = results['Rg']
Hb = results['Hb']
Rmsd = results['Rmsd']
T = results['T']
PS = results['PS']
res = results['res']
n_res = res.size

if os.environ.get('start_frame') is None:
    start_frame = max(0, min(100, max(Rmsd.shape[1] - 1, 0)))
else:
    start_frame = resolve_start_frame(Rmsd.shape[1])
end_frame1 = max(start_frame + 1, int(round(Rmsd.shape[1] * 0.60)))
end_frame2 = Rmsd.shape[1]

dGhx_HXMS = load_optional_numeric_csv('{}/{}_HXMS.csv'.format(paths['pdb_dir'], pdb_id))
dGhx_NMR = load_optional_numeric_csv('{}/{}_NMR.csv'.format(paths['pdb_dir'], pdb_id), skip_header=1)
dGhx_NMR_MS = load_optional_numeric_csv('{}/{}_NMR_MS.csv'.format(paths['pdb_dir'], pdb_id), skip_header=1)

kB = 1.0
beta = T ** -1
cE0 = Pot[:, start_frame:]
FN = cE0[0].size
FNs = np.zeros([n_rep], np.int32) + FN
reduced_pot = np.zeros([n_rep, n_rep, FN], np.float32)
for k in range(n_rep):
    for l in range(n_rep):
        reduced_pot[k, l] = beta[l] * cE0[k]
mbar0 = pymbar.MBAR(reduced_pot, FNs, verbose=True)


def reweight(T_target_value):
    u_n = (cE0 / (T_target_value * kB)).flatten()
    log_w = mbar0._computeUnnormalizedLogWeights(u_n)
    weights = np.exp(log_w)
    weights /= np.sum(weights)
    return weights


ps_conv = PS[:, start_frame:, :]


def residue_dg_from_pf_jscripts(mean_pf, temp_scale):
    # Match the legacy jscripts residue-dG mean calculation used for DG_res_T_slice.
    if mean_pf >= 1.0:
        return 1000.0
    with np.errstate(divide='ignore', invalid='ignore'):
        return 0.001987 * temp_scale * np.log(mean_pf / (1.0 - mean_pf))


def residue_dg_from_pf_step6_plot(mean_pf, temp_scale):
    # Match the deleted step-6 residue plot sentinels.
    if mean_pf >= 0.99999:
        return 1000.0
    if mean_pf <= 0.00001:
        return -100.0
    with np.errstate(divide='ignore', invalid='ignore'):
        return 0.001987 * temp_scale * np.log(mean_pf / (1.0 - mean_pf))


def residue_dg_profile_with_error(weights, temp_scale, plot_mode=False):
    weight_matrix = np.asarray(weights, dtype=float).reshape(n_rep, -1)
    sum_w_per_traj = np.sum(weight_matrix, axis=1)
    total_den = np.sum(sum_w_per_traj)

    dG_mean = np.zeros(n_res)
    dG_err = np.zeros(n_res)

    for residue_index in range(n_res):
        residue_pf = ps_conv[:, :, residue_index]
        sum_wp_per_traj = np.sum(residue_pf * weight_matrix, axis=1)
        total_num = np.sum(sum_wp_per_traj)
        mean_pf = total_num / total_den
        if plot_mode:
            dG_mean[residue_index] = residue_dg_from_pf_step6_plot(mean_pf, temp_scale)
            if n_rep < 2 or mean_pf <= 0.00001 or mean_pf >= 0.99999:
                continue
        else:
            dG_mean[residue_index] = residue_dg_from_pf_jscripts(mean_pf, temp_scale)
            if n_rep < 2 or mean_pf <= 0.0 or mean_pf >= 1.0:
                continue

        if n_rep < 2:
            continue

        jk_dg = np.empty(n_rep)
        for replica_index in range(n_rep):
            jk_den = total_den - sum_w_per_traj[replica_index]
            if jk_den <= 0.0:
                jk_dg[replica_index] = dG_mean[residue_index]
                continue

            jk_pf = (total_num - sum_wp_per_traj[replica_index]) / jk_den
            jk_pf = np.clip(jk_pf, 1.0e-6, 1.0 - 1.0e-6)
            jk_dg[replica_index] = 0.001987 * temp_scale * np.log(jk_pf / (1.0 - jk_pf))

        dG_err[residue_index] = np.sqrt(max((n_rep - 1) * np.var(jk_dg), 0.0))

    return dG_mean, dG_err


def calculate_residue_dg_profiles(temp_values, plot_mode=False):
    profile_temps = np.asarray(temp_values, dtype=float)
    dG_profiles = []
    dG_errors = []
    for temp_value in profile_temps:
        temp_scale = temp_value / 0.85 * 298.0
        dG_mean, dG_err = residue_dg_profile_with_error(reweight(temp_value), temp_scale, plot_mode=plot_mode)
        dG_profiles.append(dG_mean)
        dG_errors.append(dG_err)
    return profile_temps, np.asarray(dG_profiles), np.asarray(dG_errors)


def residue_plot_limits(values):
    return (-10.0, 30.0)


def plot_residue_profile_with_extremes(ax, x_values, y_values, y_errors, y_limits, label, color):
    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    y_errors = np.asarray(y_errors, dtype=float)

    finite_mask = np.isfinite(y_values) & np.isfinite(y_errors)
    if np.any(finite_mask):
        ax.errorbar(
            x_values[finite_mask],
            y_values[finite_mask],
            yerr=y_errors[finite_mask],
            fmt='o-',
            color=color,
            linewidth=1.0,
            markersize=3.5,
            elinewidth=0.6,
            capsize=1.5,
            label=label,
        )
    else:
        ax.plot([], [], 'o-', color=color, linewidth=1.0, markersize=3.5, label=label)

    y_min, y_max = y_limits
    high_mask = np.isposinf(y_values) | (y_values >= 100.0) | (y_values > y_max)
    low_mask = np.isneginf(y_values) | (y_values <= -100.0) | (y_values < y_min)

    if np.any(high_mask):
        ax.scatter(
            x_values[high_mask],
            np.full(np.count_nonzero(high_mask), y_max),
            marker='^',
            s=28,
            color=color,
            clip_on=False,
            zorder=4,
        )
    if np.any(low_mask):
        ax.scatter(
            x_values[low_mask],
            np.full(np.count_nonzero(low_mask), y_min),
            marker='v',
            s=28,
            color=color,
            clip_on=False,
            zorder=4,
        )


Rg_conv = Rg[:, start_frame:]
Hb_conv = Hb[:, start_frame:]
Rmsd_conv = Rmsd[:, start_frame:]

for T_target_value in target_temps:
    mbar_weights = reweight(T_target_value).reshape(n_rep, -1)

    Rg_mean = [np.average(Rg_conv[:, k], weights=mbar_weights[:, k]) for k in range(end_frame2 - start_frame)]
    Hb_mean = [np.average(Hb_conv[:, k], weights=mbar_weights[:, k]) for k in range(end_frame2 - start_frame)]
    Rmsd_mean = [np.average(Rmsd_conv[:, k], weights=mbar_weights[:, k]) for k in range(end_frame2 - start_frame)]

    split = end_frame1 - start_frame
    first_half = {
        'Rg': Rg_mean[:split],
        'Hb': Hb_mean[:split],
        'Rmsd': Rmsd_mean[:split],
    }
    second_half = {
        'Rg': Rg_mean[split:],
        'Hb': Hb_mean[split:],
        'Rmsd': Rmsd_mean[split:],
    }

    def conv_plot1(param, mean1, mean2, name):
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial']
        axes[0].set_xlabel('Frames')
        axes[0].set_ylabel('Average {}'.format(param))
        axes[0].plot(range(start_frame, end_frame1), mean1)
        axes[0].set_title('Weighted Average {} at Upside T={:.2f}'.format(param, T_target_value))
        axes[1].set_xlabel('Frames')
        axes[1].set_ylabel('Average {}'.format(param))
        axes[1].plot(range(end_frame1, end_frame2), mean2)
        axes[1].set_title('Weighted Average {} at Upside T={:.2f}'.format(param, T_target_value))
        plt.tight_layout()
        plt.savefig('{}/{}_{}_convergence_{:.2f}.png'.format(paths['result_dir'], pdb_id, name, T_target_value), dpi=300, bbox_inches='tight')
        plt.close()

    def conv_plot2(param, mean1, mean2, name):
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial']
        axes[0].set_xlabel('Average {}'.format(param))
        axes[0].set_ylabel('Density')
        axes[0].hist(mean1, bins=10, density=True)
        axes[0].set_title('Weighted Average {} at Upside T={:.2f}'.format(param, T_target_value))
        axes[1].set_xlabel('Average {}'.format(param))
        axes[1].set_ylabel('Density')
        axes[1].hist(mean2, bins=10, density=True)
        axes[1].set_title('Weighted Average {} at Upside T={:.2f}'.format(param, T_target_value))
        plt.tight_layout()
        plt.savefig('{}/{}_{}_convergence_hist_{:.2f}.png'.format(paths['result_dir'], pdb_id, name, T_target_value), dpi=300, bbox_inches='tight')
        plt.close()

    conv_plot1('Rg ($\\mathrm{\\AA}$)', first_half['Rg'], second_half['Rg'], 'Rg')
    conv_plot2('Rg ($\\mathrm{\\AA}$)', first_half['Rg'], second_half['Rg'], 'Rg')
    ks_Rg = ks_2samp(first_half['Rg'], second_half['Rg'])

    conv_plot1('# H-bond', first_half['Hb'], second_half['Hb'], 'Hb')
    conv_plot2('# H-bond', first_half['Hb'], second_half['Hb'], 'Hb')
    ks_Hb = ks_2samp(first_half['Hb'], second_half['Hb'])

    conv_plot1('Rmsd ($\\mathrm{\\AA}$)', first_half['Rmsd'], second_half['Rmsd'], 'Rmsd')
    conv_plot2('Rmsd ($\\mathrm{\\AA}$)', first_half['Rmsd'], second_half['Rmsd'], 'Rmsd')
    ks_Rmsd = ks_2samp(first_half['Rmsd'], second_half['Rmsd'])

    save_csv_rows(
        '{}/{}_{:.2f}.csv'.format(paths['result_dir'], pdb_id, T_target_value),
        [ks_Rg, ks_Hb, ks_Rmsd],
    )

Rg_flatten = Rg[:, start_frame:].flatten()
Hb_flatten = Hb[:, start_frame:].flatten()
Rg_mean = []
Hb_mean = []
Rg_reweight = []
Hb_reweight = []
dG_hbond = []

for t in T:
    weights = reweight(t)
    Rg_mean.append(np.mean(Rg[T == t, start_frame:]))
    Hb_mean.append(np.mean(Hb[T == t, start_frame:]))
    Rg_reweight.append(np.sum(Rg_flatten * weights))
    Hb_reweight.append(np.sum(Hb_flatten * weights))

    hist, bins = np.histogram(Hb_flatten, 20, weights=weights)
    hist = np.clip(hist, 1e-12, None)
    bc = (bins[1:] + bins[:-1]) * 0.5
    g = -0.593 * np.log(hist) * t
    g -= np.min(g)
    dG_hbond.append(g)

dG_hbond = np.asarray(dG_hbond)
dG_slice_temps, dGhx_slice, dGhx_slice_err = calculate_residue_dg_profiles(plot_target_temps, plot_mode=True)
dG_profile_temps, dGhx_T, dGhx_T_err = calculate_residue_dg_profiles(T)

pf_frame = np.sum(PS[:, start_frame:, :], axis=2).flatten()
den = np.linspace(0.0, 15.0, 151)
TT = target_T / 0.85 * 298.0
m = m_sens * -pf_frame
w1 = reweight(target_T)

dGhx_D = []
for d in den:
    weights = np.exp(m * d / target_T) * w1
    weights /= np.sum(weights)
    probp = np.zeros(n_res)
    for j in range(n_res):
        pf_i = PS[:, start_frame:, j].flatten()
        probp[j] = np.sum(pf_i * weights)
    probp = np.clip(probp, 1e-8, 1.0 - 1e-8)
    dGhx_D.append(TT * np.log(probp / (1.0 - probp)) * 0.001987)
dGhx_D = np.asarray(dGhx_D)
mValue = np.diff(dGhx_D, axis=0) / (den[0] - den[1])

fig, ax = plt.subplots()
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
ax.set_xlabel('T (Upside)')
ax.set_ylabel('Rg (A)', color='tab:blue')
rg_line_1 = ax.plot(T, Rg_mean, label='Rg direct average', color='tab:blue')
rg_line_2 = ax.plot(T, Rg_reweight, label='Rg MBAR reweighting', color='tab:green')
ax.tick_params(axis='y', labelcolor='tab:blue')
twin = ax.twinx()
twin.set_ylabel('#H-bond', color='tab:orange')
hb_line_1 = twin.plot(T, Hb_mean, label='# H-bonds direct average', color='tab:orange')
hb_line_2 = twin.plot(T, Hb_reweight, label='# H-bonds MBAR reweighting', color='tab:red')
twin.tick_params(axis='y', labelcolor='tab:orange')
twin.invert_yaxis()
legend_lines = rg_line_1 + rg_line_2 + hb_line_1 + hb_line_2
ax.legend(legend_lines, [line.get_label() for line in legend_lines], loc=0)
plt.title('Melting Temperature Curves')
plt.savefig('{}/{}_Tm_curve.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure()
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.hist(Hb_flatten, bins=10, density=True)
plt.xlabel('#H-bond')
plt.ylabel('Density')
plt.title('Native Fraction of #H-bond')
plt.savefig('{}/{}_histogram_Hbond.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure()
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
selected_replica = min(8, len(T) - 1)
plt.plot(bc, dG_hbond[selected_replica])
plt.title('Upside T={:.2f}'.format(T[selected_replica]))
plt.xlabel('#H-bond')
plt.ylabel('∆G (kcal/mol)')
plt.savefig('{}/{}_DG_Hbond.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plot_limits = residue_plot_limits(dGhx_slice)
ax.set_ylim(*plot_limits)
color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['tab:blue'])
for i in range(dG_slice_temps.size):
    plot_residue_profile_with_extremes(
        ax,
        res,
        dGhx_slice[i],
        dGhx_slice_err[i],
        plot_limits,
        'T = {:.2f}'.format(dG_slice_temps[i]),
        color_cycle[i % len(color_cycle)],
    )
if dGhx_NMR_MS is not None and getattr(dGhx_NMR_MS, 'ndim', 0) == 2 and dGhx_NMR_MS.shape[1] >= 2:
    x, y = dGhx_NMR_MS.T[:2]
    ax.plot(x, y, color='k', marker='o')
ax.set_xlabel('Residue Number (from N-term)')
ax.set_ylabel('∆G (kcal/mol)')
ax.set_title('∆G for Upside T Slices')
fontP = FontProperties()
fontP.set_size('x-small')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2, prop=fontP, columnspacing=0.8)
plt.savefig('{}/{}_DG_res_T_slice.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
for i in range(dG_profile_temps.size):
    ax.plot(np.sort(dGhx_T[i])[::-1], label='T = {:.2f}'.format(dG_profile_temps[i]))
ranked_hxms = None
if dGhx_HXMS is not None:
    ranked_hxms = np.sort(np.asarray(dGhx_HXMS, dtype=float))[::-1]
    ax.plot(ranked_hxms, color='k')
ax.set_xlabel('Residue (ranked by ∆G)')
ax.set_ylabel('∆G (kcal/mol)')
ax.set_title('∆G Ranked for Upside T Slices')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2, prop=fontP, columnspacing=0.8)
plt.savefig('{}/{}_DG_res_rank_T_slice.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure()
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
for i, residue_id in enumerate(res):
    plt.plot(den, dGhx_D[:, i])
plt.title('∆G Upside vs. [den]')
plt.xlabel('[den] (M)')
plt.ylabel('∆G (kcal/mol)')
plt.savefig('{}/{}_DG_den.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure()
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.bar(res, dGhx_D[0, :])
plt.title('∆G Upside vs. Sequential Residues')
plt.xlabel('Residue Number (from N-term)')
plt.ylabel('∆G (kcal/mol)')
plt.savefig('{}/{}_DG_res.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
plt.close()

width = 0.5
fig, ax1 = plt.subplots()
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
ax1.set_xlabel('Residue Number (from N-term)')
ax1.set_ylabel('∆G (kcal/mol)', color='tab:blue')
ax1.bar(res - width, dGhx_D[0, :], color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')
ax2 = ax1.twinx()
ax2.set_ylabel('m (kcal/mol*M)', color='tab:red')
ax2.bar(res + width, mValue[0, :], color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')
fig.tight_layout()
plt.title('∆G Upside and m-Value vs. Residue')
plt.savefig('{}/{}_DG_mvalue.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
plt.close()

if ranked_hxms is not None:
    diff_sorts = []
    for i in range(T.size):
        ranked_sim = np.sort(dGhx_T[i])[::-1]
        diff_sorts.append(np.abs(ranked_sim[0] - ranked_hxms[0]))
    selected_index = int(np.argmin(diff_sorts))
else:
    selected_index = int(np.argmin(np.abs(T - reference_temperature)))

selected_T = T[selected_index]
selected_profile = dGhx_T[selected_index]
selected_ranked = np.sort(selected_profile)[::-1]

if dGhx_HXMS is not None or dGhx_NMR is not None:
    fig = plt.figure()
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.plot(selected_ranked, label='Upside Simulation')
    if ranked_hxms is not None:
        plt.plot(ranked_hxms, label='HX-MS Experiment')
    if dGhx_NMR is not None:
        plt.plot(np.asarray(dGhx_NMR, dtype=float), label='NMR Experiment')
    plt.legend()
    plt.xlabel('Residue (ranked by ∆G)')
    plt.ylabel('∆G (kcal/mol)')
    plt.title('∆G Ranked at Upside T = {:.2f}'.format(selected_T))
    plt.savefig('{}/{}_DG_res_rankorder.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
    plt.close()

if dGhx_NMR_MS is not None and getattr(dGhx_NMR_MS, 'ndim', 0) == 2 and dGhx_NMR_MS.shape[1] >= 2:
    nmr_ms = np.asarray(dGhx_NMR_MS, dtype=float)
    nmr_res = nmr_ms[:, 0].astype(int)
    nmr_vals = nmr_ms[:, 1].astype(float)
    res_int = np.intersect1d(res.astype(int), nmr_res)

    if res_int.size >= 2:
        sorter_res = np.argsort(res.astype(int))
        sorter_nmr = np.argsort(nmr_res)
        sim_idx = sorter_res[np.searchsorted(res.astype(int), res_int, sorter=sorter_res)]
        nmr_idx = sorter_nmr[np.searchsorted(nmr_res, res_int, sorter=sorter_nmr)]
        exp_vals = nmr_vals[nmr_idx]

        best_index = None
        best_r2 = None
        best_profile = None
        best_mask = None
        for i in range(T.size):
            sim_vals = dGhx_T[i, sim_idx]
            finite_mask = np.isfinite(sim_vals) & np.isfinite(exp_vals)
            if np.count_nonzero(finite_mask) < 2:
                continue

            if np.allclose(sim_vals[finite_mask], sim_vals[finite_mask][0]) or np.allclose(exp_vals[finite_mask], exp_vals[finite_mask][0]):
                r_square = 0.0
            else:
                r_square = linregress(sim_vals[finite_mask], exp_vals[finite_mask]).rvalue ** 2

            if best_r2 is None or r_square > best_r2:
                best_index = i
                best_r2 = r_square
                best_profile = sim_vals
                best_mask = finite_mask

        if best_index is not None:
            fig = plt.figure()
            plt.rcParams['font.family'] = 'sans-serif'
            plt.rcParams['font.sans-serif'] = ['Arial']
            plt.plot(res, dGhx_T[best_index], '-o', label='Upside Simulation', color='tab:blue')
            plt.scatter(res_int, exp_vals, label='NMR Experiment', color='tab:orange')
            plt.legend()
            plt.xlabel('Residue Number (from N-term)')
            plt.ylabel('∆G (kcal/mol)')
            plt.title('∆G for Upside and NMR at Upside T = {:.2f}'.format(T[best_index]))
            plt.savefig('{}/{}_DG_res_NMR.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
            plt.close()

            fig = plt.figure()
            plt.rcParams['font.family'] = 'sans-serif'
            plt.rcParams['font.sans-serif'] = ['Arial']
            sim_vals = best_profile[best_mask]
            exp_plot_vals = exp_vals[best_mask]
            plt.scatter(sim_vals, exp_plot_vals, c='b')
            plt.xlabel('∆G Upside (kcal/mol)')
            plt.ylabel('∆G NMR (kcal/mol)')
            plt.title('∆G for Upside vs. NMR at Upside T = {:.2f}'.format(T[best_index]))

            if np.count_nonzero(best_mask) >= 2 and not np.allclose(sim_vals, sim_vals[0]):
                slope, intercept = np.polyfit(sim_vals, exp_plot_vals, 1)
                line = np.linspace(np.nanmin(sim_vals), np.nanmax(sim_vals), 100)
                plt.plot(line, slope * line + intercept, color='k', lw=2.5)

            plt.text(0.05, 0.95, '$R^2$ = {:.8f}'.format(best_r2), transform=plt.gca().transAxes, va='top')
            plt.savefig('{}/{}_DG_res_NMR_R2.png'.format(paths['result_dir'], pdb_id), dpi=300, bbox_inches='tight')
            plt.close()

            save_csv_rows(
                '{}/{}_NMR_compare.csv'.format(paths['result_dir'], pdb_id),
                [res_int, best_profile, exp_vals],
            )

save_csv_rows('{}/{}.csv'.format(paths['result_dir'], pdb_id), [res, selected_profile])
rank_rows = [selected_ranked]
if ranked_hxms is not None:
    rank_rows.append(ranked_hxms)
save_csv_rows('{}/{}_ranked.csv'.format(paths['result_dir'], pdb_id), rank_rows)
