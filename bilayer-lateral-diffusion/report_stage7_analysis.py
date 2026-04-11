#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


MPLCONFIGDIR = Path(tempfile.gettempdir()) / "upside-matplotlib"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPLCONFIGDIR))

import matplotlib

matplotlib.use("Agg")
from matplotlib import colors, pyplot as plt


REPORT_DIRNAME = "report"
TASK_KEYS = ["temperature", "damping", "mass_scale"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate a downloaded stage7-analysis tree and generate plots plus a recommendation summary."
    )
    parser.add_argument(
        "--analysis-dir",
        type=Path,
        required=True,
        help="Path to the downloaded bilayer-lateral-diffusion/stage7-analysis directory.",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=None,
        help="Optional output directory. Defaults to <analysis-dir>/report.",
    )
    parser.add_argument(
        "--upside-timestep-ps",
        type=float,
        default=40.0,
        help="Physical time represented by one Upside integration step in picoseconds.",
    )
    parser.add_argument(
        "--integration-step",
        type=float,
        default=0.01,
        help="Integrator time step used in the bilayer scan, in native workflow time units.",
    )
    parser.add_argument(
        "--fit-r2-threshold",
        type=float,
        default=0.95,
        help="Threshold used to flag weak MSD fits in validation and candidate ranking.",
    )
    parser.add_argument(
        "--target-um2-s",
        type=float,
        nargs="*",
        default=[0.5, 1.0, 2.0],
        help="Target DOPC diffusion coefficients in um^2/s for candidate ranking.",
    )
    parser.add_argument(
        "--max-cv-threshold",
        type=float,
        default=0.25,
        help="Maximum replicate coefficient of variation allowed in the stable analysis subset.",
    )
    parser.add_argument(
        "--mass-jump-factor",
        type=float,
        default=4.0,
        help="Maximum allowed jump in diffusion between adjacent mass points at fixed temperature and damping before the lighter-mass branch is treated as unstable.",
    )
    return parser.parse_args()


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def load_inputs(analysis_dir: Path) -> tuple[dict[str, Any], dict[str, Any], pd.DataFrame, pd.DataFrame]:
    manifest_path = require_file(analysis_dir / "analysis_manifest.json")
    summary_path = require_file(analysis_dir / "assembled" / "summary.json")
    task_csv = require_file(analysis_dir / "assembled" / "task_results.csv")
    condition_csv = require_file(analysis_dir / "assembled" / "condition_summary.csv")
    manifest = json.loads(manifest_path.read_text())
    summary = json.loads(summary_path.read_text())
    task = pd.read_csv(task_csv)
    condition = pd.read_csv(condition_csv)
    return manifest, summary, task, condition


def validate_inputs(
    analysis_dir: Path,
    manifest: dict[str, Any],
    summary: dict[str, Any],
    task: pd.DataFrame,
    condition: pd.DataFrame,
    fit_r2_threshold: float,
) -> dict[str, Any]:
    result_json_count = len(list((analysis_dir / "results" / "tasks").glob("*.json")))
    expected_conditions = (
        task[TASK_KEYS].drop_duplicates().shape[0]
        if all(key in task.columns for key in TASK_KEYS)
        else int(condition.shape[0])
    )
    required_task_columns = {
        "task_id",
        "temperature",
        "damping",
        "mass_scale",
        "replicate",
        "diffusion_mean_nm2_per_time",
        "fit_r2_mean",
    }
    required_condition_columns = {
        "temperature",
        "damping",
        "mass_scale",
        "n_replicates_expected",
        "n_replicates_completed",
        "diffusion_mean_nm2_per_time_mean",
    }

    errors: list[str] = []
    warnings: list[str] = []

    missing_task = sorted(required_task_columns - set(task.columns))
    missing_condition = sorted(required_condition_columns - set(condition.columns))
    if missing_task:
        errors.append(f"task_results.csv is missing required columns: {', '.join(missing_task)}")
    if missing_condition:
        errors.append(f"condition_summary.csv is missing required columns: {', '.join(missing_condition)}")

    manifest_total = manifest.get("n_stage7_files")
    summary_total = summary.get("n_tasks_total")
    summary_success = summary.get("n_tasks_completed_successfully")

    if manifest_total != len(task):
        errors.append(
            f"Manifest stage-7 count ({manifest_total}) does not match task_results.csv rows ({len(task)})."
        )
    if summary_total != len(task):
        errors.append(f"Summary total tasks ({summary_total}) does not match task_results.csv rows ({len(task)}).")
    if summary_success != len(task):
        errors.append(
            f"Summary successful tasks ({summary_success}) does not match task_results.csv rows ({len(task)})."
        )
    if result_json_count != len(task):
        errors.append(f"Found {result_json_count} result JSONs but {len(task)} task rows.")
    if len(condition) != expected_conditions:
        errors.append(
            f"condition_summary.csv rows ({len(condition)}) do not match unique task conditions ({expected_conditions})."
        )

    if task["task_id"].duplicated().any():
        errors.append("task_results.csv contains duplicate task_id values.")

    critical_task_columns = ["diffusion_mean_nm2_per_time", "fit_r2_mean", "temperature", "damping", "mass_scale"]
    critical_condition_columns = ["diffusion_mean_nm2_per_time_mean", "temperature", "damping", "mass_scale"]
    if task[critical_task_columns].isna().any().any():
        errors.append("task_results.csv contains NaN values in critical columns.")
    if condition[critical_condition_columns].isna().any().any():
        errors.append("condition_summary.csv contains NaN values in critical columns.")

    incomplete_conditions = condition[condition["n_replicates_completed"] != condition["n_replicates_expected"]]
    if not incomplete_conditions.empty:
        errors.append(f"{len(incomplete_conditions)} conditions have incomplete replicate counts.")

    low_r2_count = int((task["fit_r2_mean"] < fit_r2_threshold).sum())
    very_low_r2_count = int((task["fit_r2_mean"] < 0.9).sum())
    severe_low_r2_count = int((task["fit_r2_mean"] < 0.8).sum())
    if severe_low_r2_count:
        warnings.append(f"{severe_low_r2_count} task has fit_r2_mean below 0.8.")
    if very_low_r2_count:
        warnings.append(f"{very_low_r2_count} tasks have fit_r2_mean below 0.9.")
    if low_r2_count:
        warnings.append(f"{low_r2_count} tasks have fit_r2_mean below {fit_r2_threshold:.2f}.")

    if "thickness_angstrom_mean" in condition.columns:
        thickness = condition["thickness_angstrom_mean"].to_numpy(dtype=float)
        if np.any(thickness < 0.0) and np.any(thickness > 0.0):
            warnings.append(
                "thickness_angstrom_mean spans both signs in the assembled outputs; the report does not use thickness for ranking."
            )

    status = "valid" if not errors else "invalid"
    return {
        "status": status,
        "errors": errors,
        "warnings": warnings,
        "counts": {
            "manifest_n_stage7_files": manifest_total,
            "summary_n_tasks_total": summary_total,
            "summary_n_tasks_completed_successfully": summary_success,
            "task_results_rows": int(len(task)),
            "condition_summary_rows": int(len(condition)),
            "result_json_count": int(result_json_count),
            "expected_condition_rows": int(expected_conditions),
        },
        "fit_quality": {
            "fit_r2_mean_min": float(task["fit_r2_mean"].min()),
            "fit_r2_mean_median": float(task["fit_r2_mean"].median()),
            "fit_r2_mean_below_threshold_count": low_r2_count,
            "fit_r2_mean_below_0p9_count": very_low_r2_count,
            "fit_r2_mean_below_0p8_count": severe_low_r2_count,
        },
    }


def add_physical_units(df: pd.DataFrame, time_unit_ns: float, mean_column: str, std_column: str | None = None) -> pd.DataFrame:
    out = df.copy()
    out[f"{mean_column}_nm2_per_ns"] = out[mean_column] / time_unit_ns
    out[f"{mean_column}_um2_per_s"] = out[f"{mean_column}_nm2_per_ns"] * 1000.0
    if std_column and std_column in out.columns:
        out[f"{std_column}_nm2_per_ns"] = out[std_column] / time_unit_ns
        out[f"{std_column}_um2_per_s"] = out[f"{std_column}_nm2_per_ns"] * 1000.0
    return out


def geometric_mean(values: pd.Series | np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr) & (arr > 0.0)]
    if arr.size == 0:
        return float("nan")
    return float(np.exp(np.log(arr).mean()))


def build_condition_metrics(condition: pd.DataFrame, task: pd.DataFrame, fit_r2_threshold: float) -> pd.DataFrame:
    task_metrics = (
        task.groupby(TASK_KEYS)
        .agg(
            task_count=("task_id", "count"),
            fit_r2_mean_mean=("fit_r2_mean", "mean"),
            fit_r2_mean_min=("fit_r2_mean", "min"),
            fit_r2_below_threshold_count=("fit_r2_mean", lambda x: int((x < fit_r2_threshold).sum())),
            diffusion_nm2_per_time_cv=("diffusion_mean_nm2_per_time", lambda x: float(x.std(ddof=1) / x.mean()) if x.mean() > 0 else float("nan")),
        )
        .reset_index()
    )
    return condition.merge(task_metrics, on=TASK_KEYS, how="left")


def mark_mass_continuity_stability(
    condition: pd.DataFrame,
    diffusion_column: str,
    jump_factor: float,
) -> pd.Series:
    stable = pd.Series(True, index=condition.index, dtype=bool)
    diffusion = pd.to_numeric(condition[diffusion_column], errors="coerce")
    stable &= np.isfinite(diffusion.to_numpy(dtype=float))

    for _, group in condition.groupby(["temperature", "damping"]):
        ordered = group.sort_values("mass_scale", ascending=False)
        idx = list(ordered.index)
        values = pd.to_numeric(ordered[diffusion_column], errors="coerce").to_numpy(dtype=float)
        branch_is_stable = True
        for i, value in enumerate(values):
            if not branch_is_stable or not (math.isfinite(value) and value > 0.0):
                stable.loc[idx[i]] = False
                branch_is_stable = False
                continue
            if i == 0:
                continue
            heavy = values[i - 1]
            ratio = value / heavy
            if ratio > jump_factor or ratio < (1.0 / jump_factor):
                stable.loc[idx[i]] = False
                branch_is_stable = False
    return stable


def build_stable_condition_subset(
    condition: pd.DataFrame,
    fit_r2_threshold: float,
    max_cv_threshold: float,
    mass_jump_factor: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    working = condition.copy()
    working["continuity_stable"] = mark_mass_continuity_stability(
        working,
        "diffusion_mean_nm2_per_time_mean_um2_per_s",
        mass_jump_factor,
    )
    finite_diffusion = np.isfinite(pd.to_numeric(working["diffusion_mean_nm2_per_time_mean_um2_per_s"], errors="coerce"))
    finite_cv = np.isfinite(pd.to_numeric(working["diffusion_nm2_per_time_cv"], errors="coerce"))
    stable_mask = (
        finite_diffusion
        & finite_cv
        & (working["n_replicates_completed"] == working["n_replicates_expected"])
        & (working["fit_r2_mean_min"] >= fit_r2_threshold)
        & (working["diffusion_nm2_per_time_cv"] <= max_cv_threshold)
        & working["continuity_stable"]
    )
    stable = working.loc[stable_mask].copy()
    summary = {
        "n_total_conditions": int(len(working)),
        "n_stable_conditions": int(len(stable)),
        "n_nonfinite_conditions": int((~finite_diffusion).sum()),
        "n_high_cv_conditions": int((finite_cv & (working["diffusion_nm2_per_time_cv"] > max_cv_threshold)).sum()),
        "n_discontinuous_conditions": int((~working["continuity_stable"]).sum()),
        "max_cv_threshold": float(max_cv_threshold),
        "mass_jump_factor": float(mass_jump_factor),
    }
    return working, stable, summary


def build_main_effect_table(condition: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for factor in ["temperature", "mass_scale", "damping"]:
        for value, group in condition.groupby(factor):
            rows.append(
                {
                    "factor": factor,
                    "value": float(value),
                    "geometric_mean_um2_per_s": geometric_mean(group["diffusion_mean_nm2_per_time_mean_um2_per_s"]),
                    "arithmetic_mean_um2_per_s": float(group["diffusion_mean_nm2_per_time_mean_um2_per_s"].mean()),
                    "median_um2_per_s": float(group["diffusion_mean_nm2_per_time_mean_um2_per_s"].median()),
                    "n_conditions": int(len(group)),
                }
            )
    return pd.DataFrame(rows)


def select_target_candidates(condition: pd.DataFrame, targets: list[float]) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    robust = condition[
        (condition["n_replicates_completed"] == condition["n_replicates_expected"])
        & (condition["fit_r2_mean_min"] >= 0.95)
        & (condition["diffusion_nm2_per_time_cv"] <= 0.25)
    ].copy()
    fallback = condition.copy()
    for target in targets:
        pool = robust if not robust.empty else fallback
        ranked = pool.copy()
        ranked["target_um2_per_s"] = target
        ranked["abs_error_um2_per_s"] = (ranked["diffusion_mean_nm2_per_time_mean_um2_per_s"] - target).abs()
        ranked = ranked.sort_values(
            ["abs_error_um2_per_s", "diffusion_nm2_per_time_cv", "fit_r2_mean_min", "temperature"],
            ascending=[True, True, False, True],
        ).head(5)
        ranked["used_robust_filter"] = bool(pool is robust)
        rows.append(ranked)
    return pd.concat(rows, ignore_index=True)


def plot_heatmaps(condition: pd.DataFrame, output_path: Path) -> None:
    temperatures = sorted(condition["temperature"].unique())
    damping_values = sorted(condition["damping"].unique())
    mass_values = sorted(condition["mass_scale"].unique(), reverse=True)
    positive = condition["diffusion_mean_nm2_per_time_mean_um2_per_s"]
    vmin = float(positive[positive > 0.0].min())
    vmax = float(positive.max())
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(1, len(temperatures), figsize=(4.0 * len(temperatures), 4.4), constrained_layout=True)
    if len(temperatures) == 1:
        axes = [axes]
    image = None
    for axis, temperature in zip(axes, temperatures):
        subset = condition[condition["temperature"] == temperature]
        pivot = (
            subset.pivot(index="mass_scale", columns="damping", values="diffusion_mean_nm2_per_time_mean_um2_per_s")
            .reindex(index=mass_values, columns=damping_values)
        )
        image = axis.imshow(pivot.to_numpy(dtype=float), aspect="auto", interpolation="nearest", norm=norm, cmap="viridis")
        axis.set_title(f"T = {temperature:.1f}")
        axis.set_xticks(range(len(damping_values)))
        axis.set_xticklabels([f"{value:g}" for value in damping_values], rotation=45, ha="right")
        axis.set_xlabel("Damping timescale")
        axis.set_yticks(range(len(mass_values)))
        axis.set_yticklabels([f"{value:g}" for value in mass_values])
        axis.set_ylabel("Mass scale")
    assert image is not None
    cbar = fig.colorbar(image, ax=axes, shrink=0.82)
    cbar.set_label("Diffusion (um^2/s)")
    fig.suptitle("DOPC lateral diffusion by temperature, damping, and mass scale", fontsize=13)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_main_effects(main_effects: pd.DataFrame, output_path: Path) -> None:
    order = {
        "temperature": sorted(main_effects.loc[main_effects["factor"] == "temperature", "value"].unique()),
        "mass_scale": sorted(main_effects.loc[main_effects["factor"] == "mass_scale", "value"].unique()),
        "damping": sorted(main_effects.loc[main_effects["factor"] == "damping", "value"].unique()),
    }
    titles = {
        "temperature": "Temperature effect",
        "mass_scale": "Mass-scale effect",
        "damping": "Damping effect",
    }

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.8), constrained_layout=True)
    for axis, factor in zip(axes, ["temperature", "mass_scale", "damping"]):
        subset = main_effects[main_effects["factor"] == factor].copy()
        subset["value"] = pd.Categorical(subset["value"], categories=order[factor], ordered=True)
        subset = subset.sort_values("value")
        axis.plot(
            subset["value"].astype(float),
            subset["geometric_mean_um2_per_s"],
            marker="o",
            linewidth=2.0,
            color="#1f4e79",
        )
        axis.set_yscale("log")
        axis.set_title(titles[factor])
        axis.set_xlabel(factor.replace("_", " "))
        axis.set_ylabel("Geometric mean diffusion (um^2/s)")
        axis.grid(alpha=0.25, which="both")
    fig.suptitle("Main effects across the scanned dryMARTINI controls", fontsize=13)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_interactions(condition: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.0), constrained_layout=True)
    palette = ["#0b3c5d", "#328cc1", "#d9b310", "#b82601", "#7d1935"]
    mass_values = sorted(condition["mass_scale"].unique())

    for color, mass_scale in zip(palette, mass_values):
        subset = condition[condition["mass_scale"] == mass_scale]
        temp_curve = (
            subset.groupby("temperature")["diffusion_mean_nm2_per_time_mean_um2_per_s"].apply(geometric_mean).reset_index()
        )
        damp_curve = subset.groupby("damping")["diffusion_mean_nm2_per_time_mean_um2_per_s"].apply(geometric_mean).reset_index()
        axes[0].plot(temp_curve["temperature"], temp_curve["diffusion_mean_nm2_per_time_mean_um2_per_s"], marker="o", label=f"mass {mass_scale:g}", color=color)
        axes[1].plot(damp_curve["damping"], damp_curve["diffusion_mean_nm2_per_time_mean_um2_per_s"], marker="o", label=f"mass {mass_scale:g}", color=color)

    axes[0].set_title("Temperature trends by mass scale")
    axes[0].set_xlabel("Temperature")
    axes[0].set_ylabel("Geometric mean diffusion (um^2/s)")
    axes[0].set_yscale("log")
    axes[0].grid(alpha=0.25, which="both")

    axes[1].set_title("Damping trends by mass scale")
    axes[1].set_xlabel("Damping timescale")
    axes[1].set_ylabel("Geometric mean diffusion (um^2/s)")
    axes[1].set_yscale("log")
    axes[1].grid(alpha=0.25, which="both")
    axes[1].legend(loc="best", fontsize=8)

    fig.suptitle("Mass scale dominates the diffusion response", fontsize=13)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_fit_quality(task: pd.DataFrame, fit_r2_threshold: float, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0), constrained_layout=True)
    axes[0].hist(task["fit_r2_mean"], bins=30, color="#1f4e79", edgecolor="white")
    axes[0].axvline(fit_r2_threshold, color="#b82601", linestyle="--", linewidth=1.5)
    axes[0].set_xlabel("Fit R^2")
    axes[0].set_ylabel("Task count")
    axes[0].set_title("MSD fit-quality distribution")

    axes[1].scatter(
        task["diffusion_mean_nm2_per_time_um2_per_s"],
        task["fit_r2_mean"],
        s=20,
        alpha=0.7,
        color="#328cc1",
        edgecolors="none",
    )
    axes[1].axhline(fit_r2_threshold, color="#b82601", linestyle="--", linewidth=1.5)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Diffusion (um^2/s)")
    axes[1].set_ylabel("Fit R^2")
    axes[1].set_title("Fit quality versus measured diffusion")

    fig.suptitle("Stage-7 MSD fit checks", fontsize=13)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def build_report_markdown(
    validation: dict[str, Any],
    time_unit_ns: float,
    main_effects: pd.DataFrame,
    target_candidates: pd.DataFrame,
    condition: pd.DataFrame,
    stability_summary: dict[str, Any],
    full_condition: pd.DataFrame,
) -> str:
    temperature_values = sorted(float(value) for value in condition["temperature"].unique())
    mass_values = sorted(float(value) for value in condition["mass_scale"].unique())
    damping_values = sorted(float(value) for value in condition["damping"].unique())
    mass_effect = main_effects[main_effects["factor"] == "mass_scale"]["geometric_mean_um2_per_s"]
    temp_effect = main_effects[main_effects["factor"] == "temperature"]["geometric_mean_um2_per_s"]
    damping_effect = main_effects[main_effects["factor"] == "damping"]["geometric_mean_um2_per_s"]
    mass_fold = float(mass_effect.max() / mass_effect.min())
    temp_fold = float(temp_effect.max() / temp_effect.min())
    damping_fold = float(damping_effect.max() / damping_effect.min())
    fastest = condition.sort_values("diffusion_mean_nm2_per_time_mean_um2_per_s", ascending=False).iloc[0]
    fastest_full = full_condition.sort_values("diffusion_mean_nm2_per_time_mean_um2_per_s", ascending=False).iloc[0]
    robust = condition[
        (condition["n_replicates_completed"] == condition["n_replicates_expected"])
        & (condition["fit_r2_mean_min"] >= 0.95)
        & (condition["diffusion_nm2_per_time_cv"] <= 0.25)
    ].copy()
    best_robust = (
        robust.sort_values("diffusion_mean_nm2_per_time_mean_um2_per_s", ascending=False).iloc[0]
        if not robust.empty
        else None
    )
    damping_cv = (
        condition.groupby("damping")["diffusion_nm2_per_time_cv"].mean().sort_values()
        if "diffusion_nm2_per_time_cv" in condition.columns
        else pd.Series(dtype=float)
    )
    noisiest_damping = float(damping_cv.idxmax()) if not damping_cv.empty else float("nan")
    noisiest_damping_cv = float(damping_cv.max()) if not damping_cv.empty else float("nan")
    lowest_mass = mass_values[0]
    lowest_mass_subset = condition[condition["mass_scale"] == lowest_mass]
    noisy_lowest_mass = lowest_mass_subset[lowest_mass_subset["diffusion_nm2_per_time_cv"] > 0.25]

    lines = [
        "# Stage-7 Diffusion Report",
        "",
        "## Validation",
        f"- Status: {validation['status']}",
        f"- Stage-7 files / task rows / task JSONs: {validation['counts']['manifest_n_stage7_files']} / {validation['counts']['task_results_rows']} / {validation['counts']['result_json_count']}",
        f"- Condition rows: {validation['counts']['condition_summary_rows']}",
        f"- Fit quality: median R^2 = {validation['fit_quality']['fit_r2_mean_median']:.4f}, minimum R^2 = {validation['fit_quality']['fit_r2_mean_min']:.4f}",
        f"- Tasks below R^2 < 0.90: {validation['fit_quality']['fit_r2_mean_below_0p9_count']}",
        f"- Tasks below R^2 < 0.95: {validation['fit_quality']['fit_r2_mean_below_threshold_count']}",
        "",
        "## Calibration",
        f"- Upside physical calibration: 40 ps per integration step.",
        f"- Workflow integration step: 0.01 native time units per step.",
        f"- Derived conversion: 1 workflow time unit = {time_unit_ns:.3f} ns.",
        "- Diffusion conversion:",
        "  - nm^2/ns = diffusion_mean_nm2_per_time / 4",
        "  - um^2/s = diffusion_mean_nm2_per_time * 250",
        "",
        "## Main Findings",
        f"- Stable analysis subset: {stability_summary['n_stable_conditions']} / {stability_summary['n_total_conditions']} conditions passed finite-value, fit, CV, and continuity filters.",
        f"- Mass scale changes geometric-mean diffusion by {mass_fold:.1f}x across mass = {mass_values[-1]:g} down to {mass_values[0]:g}.",
        f"- Temperature changes geometric-mean diffusion by {temp_fold:.1f}x across T = {temperature_values[0]:.1f} to {temperature_values[-1]:.1f}.",
        f"- Fastest condition that passed the robustness filters: T = {fastest['temperature']:.1f}, tau = {fastest['damping']:.0f}, mass = {fastest['mass_scale']:.2g}, diffusion = {fastest['diffusion_mean_nm2_per_time_mean_um2_per_s']:.3f} um^2/s.",
        "",
        "## Target-Matching Candidates",
    ]
    if len(damping_values) == 1:
        lines.insert(lines.index("## Target-Matching Candidates"), f"- Damping was fixed at tau = {damping_values[0]:g} in this scan, so there is no in-grid damping comparison.")
    else:
        lines.insert(
            lines.index("## Target-Matching Candidates"),
            f"- Damping changes geometric-mean diffusion by {damping_fold:.2f}x across tau = {damping_values[0]:g} to {damping_values[-1]:g}.",
        )

    for target in sorted(target_candidates["target_um2_per_s"].unique()):
        lines.append(f"- Target {target:.2f} um^2/s:")
        subset = target_candidates[target_candidates["target_um2_per_s"] == target].head(3)
        for row in subset.itertuples(index=False):
            lines.append(
                f"  - T = {row.temperature:.1f}, tau = {row.damping:.0f}, mass = {row.mass_scale:.2g}, diffusion = {row.diffusion_mean_nm2_per_time_mean_um2_per_s:.3f} um^2/s, CV = {row.diffusion_nm2_per_time_cv:.3f}, min R^2 = {row.fit_r2_mean_min:.3f}"
            )

    lines.extend(
        [
            "",
            "## Recommendation",
        ]
    )

    if mass_fold >= temp_fold:
        lines.extend(
            [
                "- Lower mass does accelerate the bilayer inside the robust subset, but the gain is still modest compared with the remaining timescale mismatch.",
            ]
        )
    else:
        lines.extend(
            [
                f"- Temperature still changes diffusion more strongly than mass within the robust subset of the tested `mass = {mass_values[-1]:g} -> {mass_values[0]:g}` window.",
                "- The joint mass+damping scan improves bilayer fluidity, but the trustworthy part of the grid is still too slow to match the physical DOPC target at the same temperature.",
            ]
        )
    if best_robust is not None:
        lines.append(
            f"- Best robust point in the current grid: T = {best_robust['temperature']:.1f}, tau = {best_robust['damping']:.0f}, mass = {best_robust['mass_scale']:.3g}, diffusion = {best_robust['diffusion_mean_nm2_per_time_mean_um2_per_s']:.3f} um^2/s, CV = {best_robust['diffusion_nm2_per_time_cv']:.3f}."
        )
    if fastest_full["diffusion_mean_nm2_per_time_mean_um2_per_s"] > fastest["diffusion_mean_nm2_per_time_mean_um2_per_s"] * 2.0:
        if fastest_full["diffusion_mean_nm2_per_time_mean_um2_per_s"] > 100.0:
            lines.append(
                f"- The full grid reaches pathologically large diffusion on excluded branches at T = {fastest_full['temperature']:.1f}, tau = {fastest_full['damping']:.0f}, mass = {fastest_full['mass_scale']:.3g}; those runs are numerical-instability frontier points, not calibration candidates."
            )
        else:
            lines.append(
                f"- The full grid reaches much faster diffusion on excluded branches, up to {fastest_full['diffusion_mean_nm2_per_time_mean_um2_per_s']:.3f} um^2/s at T = {fastest_full['temperature']:.1f}, tau = {fastest_full['damping']:.0f}, mass = {fastest_full['mass_scale']:.3g}, but those branches fail the robustness filters."
            )
    if not noisy_lowest_mass.empty:
        noisy_temps = ", ".join(f"{value:.1f}" for value in sorted(set(noisy_lowest_mass["temperature"].tolist())))
        lines.append(
            f"- The lowest tested mass `{lowest_mass:g}` is already noisy at T = {noisy_temps}; treat those highest-diffusion points as exploratory rather than calibration anchors."
        )
    if math.isfinite(noisiest_damping) and noisiest_damping_cv > 0.3:
        lines.append(
            f"- The `tau = {noisiest_damping:g}` branch is much noisier than the rest of the grid (mean CV = {noisiest_damping_cv:.3f}) and should not be used as the main calibration anchor."
        )
    if condition["diffusion_mean_nm2_per_time_mean_um2_per_s"].max() < 1.0:
        lines.append(
            "- No condition that passes the robustness filters reaches 1.0 um^2/s; target-like diffusion appears only on excluded high-variance or discontinuous branches."
        )
    lines.extend(
        [
            "- Within the robust low-CV subset, the closest observed points are still well below a physical DOPC diffusion target at the same temperature.",
            "- Practical implication: the joint scan shows where the fast-but-noisy frontier begins, but it does not yet deliver a robust dryMARTINI calibration point that matches the Upside timescale proxy.",
        ]
    )

    if validation["warnings"]:
        lines.extend(["", "## Warnings"])
        for warning in validation["warnings"]:
            lines.append(f"- {warning}")
    lines.extend(
        [
            f"- {stability_summary['n_nonfinite_conditions']} conditions produced non-finite diffusion outputs.",
            f"- {stability_summary['n_high_cv_conditions']} conditions exceeded CV > {stability_summary['max_cv_threshold']:.2f}.",
            f"- {stability_summary['n_discontinuous_conditions']} conditions were excluded by the adjacent-mass jump filter (factor > {stability_summary['mass_jump_factor']:.0f}).",
        ]
    )

    return "\n".join(lines) + "\n"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> int:
    args = parse_args()
    analysis_dir = args.analysis_dir.resolve()
    report_dir = args.report_dir.resolve() if args.report_dir else (analysis_dir / REPORT_DIRNAME)
    report_dir.mkdir(parents=True, exist_ok=True)

    manifest, summary, task, condition = load_inputs(analysis_dir)
    validation = validate_inputs(analysis_dir, manifest, summary, task, condition, args.fit_r2_threshold)

    time_unit_ns = args.upside_timestep_ps / (1000.0 * args.integration_step)

    task_physical = add_physical_units(task, time_unit_ns, "diffusion_mean_nm2_per_time")
    condition_physical = add_physical_units(
        condition,
        time_unit_ns,
        "diffusion_mean_nm2_per_time_mean",
        std_column="diffusion_mean_nm2_per_time_std",
    )
    condition_metrics = build_condition_metrics(condition_physical, task_physical, args.fit_r2_threshold)
    condition_metrics, stable_condition_metrics, stability_summary = build_stable_condition_subset(
        condition_metrics,
        fit_r2_threshold=args.fit_r2_threshold,
        max_cv_threshold=args.max_cv_threshold,
        mass_jump_factor=args.mass_jump_factor,
    )
    analysis_condition = stable_condition_metrics if not stable_condition_metrics.empty else condition_metrics
    main_effects = build_main_effect_table(analysis_condition)
    target_candidates = select_target_candidates(analysis_condition, sorted(args.target_um2_s))

    plot_heatmaps(analysis_condition, report_dir / "diffusion_heatmaps_um2_s.png")
    plot_main_effects(main_effects, report_dir / "main_effects_um2_s.png")
    plot_interactions(analysis_condition, report_dir / "interaction_trends_um2_s.png")
    plot_fit_quality(task_physical, args.fit_r2_threshold, report_dir / "fit_quality.png")

    task_physical.to_csv(report_dir / "task_results_physical_units.csv", index=False)
    condition_metrics.to_csv(report_dir / "condition_summary_physical_units.csv", index=False)
    stable_condition_metrics.to_csv(report_dir / "stable_condition_summary_physical_units.csv", index=False)
    main_effects.to_csv(report_dir / "main_effects_um2_s.csv", index=False)
    target_candidates.to_csv(report_dir / "target_candidates_um2_s.csv", index=False)

    write_json(report_dir / "validation_summary.json", validation)
    recommendation = {
        "upside_timestep_ps": args.upside_timestep_ps,
        "integration_step": args.integration_step,
        "time_unit_ns": time_unit_ns,
        "targets_um2_per_s": sorted(args.target_um2_s),
        "stable_grid_max_um2_per_s": float(analysis_condition["diffusion_mean_nm2_per_time_mean_um2_per_s"].max()),
        "full_grid_max_um2_per_s": float(condition_metrics["diffusion_mean_nm2_per_time_mean_um2_per_s"].max()),
        "stable_subset_summary": stability_summary,
        "main_effects_fold_change": {
            "temperature": float(
                main_effects[main_effects["factor"] == "temperature"]["geometric_mean_um2_per_s"].max()
                / main_effects[main_effects["factor"] == "temperature"]["geometric_mean_um2_per_s"].min()
            ),
            "mass_scale": float(
                main_effects[main_effects["factor"] == "mass_scale"]["geometric_mean_um2_per_s"].max()
                / main_effects[main_effects["factor"] == "mass_scale"]["geometric_mean_um2_per_s"].min()
            ),
            "damping": float(
                main_effects[main_effects["factor"] == "damping"]["geometric_mean_um2_per_s"].max()
                / main_effects[main_effects["factor"] == "damping"]["geometric_mean_um2_per_s"].min()
            ),
        },
        "target_candidates": json.loads(target_candidates.to_json(orient="records")),
    }
    write_json(report_dir / "recommendation_summary.json", recommendation)
    report_dir.joinpath("report.md").write_text(
        build_report_markdown(
            validation,
            time_unit_ns,
            main_effects,
            target_candidates,
            analysis_condition,
            stability_summary,
            condition_metrics,
        )
    )

    print(f"Validation status: {validation['status']}")
    print(f"Report directory: {report_dir}")
    print(f"Time conversion: 1 workflow time unit = {time_unit_ns:.3f} ns")
    print("Generated plots:")
    for name in [
        "diffusion_heatmaps_um2_s.png",
        "main_effects_um2_s.png",
        "interaction_trends_um2_s.png",
        "fit_quality.png",
    ]:
        print(f"- {report_dir / name}")
    return 0 if validation["status"] == "valid" else 1


if __name__ == "__main__":
    raise SystemExit(main())
