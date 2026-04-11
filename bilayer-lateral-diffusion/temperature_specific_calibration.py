#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


MPLCONFIGDIR = Path(tempfile.gettempdir()) / "upside-matplotlib"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPLCONFIGDIR))

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt


BAG2014_TEMPERATURES_K = np.array([298.0, 303.0, 308.0, 313.0], dtype=float)
BAG2014_DIFFUSION_UM2_S = np.array([2.46, 2.77, 3.07, 3.40], dtype=float)
BAG2014_REPORTED_EA_KJ_MOL = 17.66
BAG2014_CITATION = (
    "Bag N, Yap DHX, Wohland T. Temperature dependence of diffusion in model and live cell membranes "
    "characterized by imaging fluorescence correlation spectroscopy. Biochim Biophys Acta Biomembr. 2014;1838(3):802-813. "
    "doi:10.1016/j.bbamem.2013.10.009"
)
R_GAS_J_MOL_K = 8.31446261815324


@dataclass
class ArrheniusModel:
    intercept: float
    slope: float
    activation_energy_kj_mol: float

    def predict_um2_s(self, temperature_k: float) -> float:
        return float(math.exp(self.intercept + self.slope / temperature_k))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calibrate per-temperature bilayer mass/damping corrections from the stage7-analysis report outputs."
    )
    parser.add_argument(
        "--analysis-dir",
        type=Path,
        required=True,
        help="Path to bilayer-lateral-diffusion/stage7-analysis.",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=None,
        help="Optional output directory. Defaults to <analysis-dir>/report.",
    )
    parser.add_argument(
        "--upside-temperature-k",
        type=float,
        default=350.588235,
        help="Kelvin represented by 1.0 Upside temperature unit.",
    )
    parser.add_argument(
        "--stable-damping-values",
        type=str,
        default="5,8,12,16",
        help="Comma-separated damping values considered stable for correction selection.",
    )
    parser.add_argument(
        "--fit-r2-threshold",
        type=float,
        default=0.95,
        help="Minimum per-condition fit quality allowed in the correction analysis.",
    )
    parser.add_argument(
        "--max-cv-threshold",
        type=float,
        default=0.25,
        help="Maximum replicate coefficient of variation allowed in the robust correction subset.",
    )
    parser.add_argument(
        "--minimum-reportable-mass-scale",
        type=float,
        default=0.01,
        help="Do not present extrapolated mass corrections smaller than this as resolved recommendations.",
    )
    parser.add_argument(
        "--minimum-useful-mass-slope",
        type=float,
        default=0.05,
        help="Minimum absolute value of the fitted mass exponent required before reporting an extrapolated mass correction.",
    )
    return parser.parse_args()


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def fit_arrhenius() -> ArrheniusModel:
    x = 1.0 / BAG2014_TEMPERATURES_K
    y = np.log(BAG2014_DIFFUSION_UM2_S)
    slope, intercept = np.polyfit(x, y, 1)
    activation_energy_kj_mol = float(-slope * R_GAS_J_MOL_K / 1000.0)
    return ArrheniusModel(intercept=float(intercept), slope=float(slope), activation_energy_kj_mol=activation_energy_kj_mol)


def load_condition_summary(report_dir: Path) -> pd.DataFrame:
    path = report_dir / "condition_summary_physical_units.csv"
    require_file(path)
    return pd.read_csv(path)


def resolve_stable_dampings(condition: pd.DataFrame, requested: set[float]) -> list[float]:
    available = sorted({float(value) for value in condition["damping"].to_numpy(dtype=float)})
    if not available:
        raise RuntimeError("No damping values available in condition summary.")
    if not requested:
        return available
    selected = sorted(value for value in available if value in requested)
    return selected if selected else available


def fit_mass_curve(subset: pd.DataFrame) -> tuple[float, float, float]:
    x = np.log(subset["mass_scale"].to_numpy(dtype=float))
    y = np.log(subset["diffusion_mean_nm2_per_time_mean_um2_per_s"].to_numpy(dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    pred = intercept + slope * x
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    return float(slope), float(intercept), float(r2)


def build_robust_subset(
    condition: pd.DataFrame,
    fit_r2_threshold: float,
    max_cv_threshold: float,
) -> pd.DataFrame:
    robust = condition.copy()
    robust = robust[np.isfinite(robust["diffusion_mean_nm2_per_time_mean_um2_per_s"].to_numpy(dtype=float))]
    robust = robust[np.isfinite(robust["diffusion_nm2_per_time_cv"].to_numpy(dtype=float))]
    robust = robust[robust["n_replicates_completed"] == robust["n_replicates_expected"]]
    robust = robust[robust["fit_r2_mean_min"] >= fit_r2_threshold]
    robust = robust[robust["diffusion_nm2_per_time_cv"] <= max_cv_threshold]
    if "continuity_stable" in robust.columns:
        robust = robust[robust["continuity_stable"].fillna(False)]
    return robust.copy()


def solve_mass_target(target_um2_s: float, slope: float, intercept: float) -> tuple[float, float]:
    if not math.isfinite(target_um2_s) or not math.isfinite(slope) or not math.isfinite(intercept):
        return float("nan"), float("nan")
    if target_um2_s <= 0.0 or slope >= 0.0:
        return float("nan"), float("nan")
    log_mass = (math.log(target_um2_s) - intercept) / slope
    if log_mass < -700.0:
        return 0.0, log_mass
    if log_mass > 700.0:
        return float("inf"), log_mass
    return float(math.exp(log_mass)), log_mass


def temperature_status(temperature_k: float) -> str:
    low = float(BAG2014_TEMPERATURES_K.min())
    high = float(BAG2014_TEMPERATURES_K.max())
    if low <= temperature_k <= high:
        return "within_source_range"
    if temperature_k < low:
        return "low_temperature_extrapolation"
    return "high_temperature_extrapolation"


def choose_recommended_tau(
    subset: pd.DataFrame,
    stable_dampings: set[float],
    target_um2_s: float,
    minimum_reportable_mass_scale: float,
    minimum_useful_mass_slope: float,
) -> dict[str, float] | None:
    candidates: list[dict[str, float]] = []
    for damping, group in subset.groupby("damping"):
        damping_value = float(damping)
        if damping_value not in stable_dampings:
            continue
        if len(group) < 3:
            continue
        slope, intercept, mass_fit_r2 = fit_mass_curve(group)
        mass_target, log_mass_target = solve_mass_target(target_um2_s, slope, intercept)
        is_reliable = (
            math.isfinite(mass_target)
            and mass_target >= minimum_reportable_mass_scale
            and slope <= -minimum_useful_mass_slope
        )
        candidates.append(
            {
                "damping": damping_value,
                "mass_slope": slope,
                "mass_intercept": intercept,
                "mass_fit_r2": mass_fit_r2,
                "mean_cv": float(group["diffusion_nm2_per_time_cv"].mean()),
                "pred_um2_s_at_0p1": float(math.exp(intercept + slope * math.log(0.1))),
                "required_mass_scale": mass_target,
                "required_log_mass_scale": log_mass_target,
                "is_reliable": float(1.0 if is_reliable else 0.0),
            }
        )
    if not candidates:
        raise RuntimeError("No stable damping candidates available for temperature-specific calibration.")
    reliable_candidates = [row for row in candidates if bool(row["is_reliable"])]
    if not reliable_candidates:
        return None
    ranked = sorted(
        reliable_candidates,
        key=lambda row: (
            row["mean_cv"],
            -row["required_mass_scale"],
            -row["mass_fit_r2"],
            row["damping"],
        ),
    )
    return ranked[0]


def choose_nearest_condition(
    subset: pd.DataFrame,
    target_um2_s: float,
    stable_dampings: set[float],
) -> pd.Series:
    stable_subset = subset[subset["damping"].isin(stable_dampings)].copy()
    pool = stable_subset if not stable_subset.empty else subset.copy()
    pool["log_error"] = np.abs(
        np.log(pool["diffusion_mean_nm2_per_time_mean_um2_per_s"].to_numpy(dtype=float) / target_um2_s)
    )
    pool = pool.sort_values(
        ["log_error", "diffusion_nm2_per_time_cv", "mass_scale", "damping"],
        ascending=[True, True, True, True],
    )
    return pool.iloc[0]


def build_correction_table(
    condition: pd.DataFrame,
    upside_temperature_k: float,
    stable_dampings: set[float],
    fit_r2_threshold: float,
    max_cv_threshold: float,
    minimum_reportable_mass_scale: float,
    minimum_useful_mass_slope: float,
) -> tuple[pd.DataFrame, ArrheniusModel]:
    arrhenius = fit_arrhenius()
    rows: list[dict[str, float | str]] = []

    robust = build_robust_subset(condition, fit_r2_threshold=fit_r2_threshold, max_cv_threshold=max_cv_threshold)
    if robust.empty:
        raise RuntimeError("No conditions satisfy the requested robust-subset filters.")

    for temperature_up, subset in robust.groupby("temperature"):
        temperature_up = float(temperature_up)
        temperature_k = temperature_up * upside_temperature_k
        temperature_c = temperature_k - 273.15
        target_um2_s = arrhenius.predict_um2_s(temperature_k)
        status = temperature_status(temperature_k)

        nearest = choose_nearest_condition(subset, target_um2_s, stable_dampings)

        recommended_tau = choose_recommended_tau(
            subset,
            stable_dampings,
            target_um2_s,
            minimum_reportable_mass_scale=minimum_reportable_mass_scale,
            minimum_useful_mass_slope=minimum_useful_mass_slope,
        )
        if recommended_tau is None:
            note = "current scan too slow to resolve a reliable mass correction; need lower masses or a wider scan"
            recommended_damping = float(nearest["damping"])
            recommended_mass_scale = float("nan")
            recommended_mass_fit_r2 = float("nan")
            recommended_mass_slope = float("nan")
            recommended_pred_um2_s_at_mass_0p1 = float("nan")
            estimated_mass_scale = float("nan")
        else:
            if recommended_tau["required_mass_scale"] >= subset["mass_scale"].min():
                note = "in-grid"
            elif recommended_tau["required_mass_scale"] >= minimum_reportable_mass_scale:
                note = "requires mass below scan minimum"
            else:
                note = "current scan too slow to resolve a reliable mass correction; need lower masses or a wider scan"
            recommended_damping = float(recommended_tau["damping"])
            recommended_mass_scale = (
                float(recommended_tau["required_mass_scale"])
                if recommended_tau["required_mass_scale"] >= minimum_reportable_mass_scale
                else float("nan")
            )
            recommended_mass_fit_r2 = float(recommended_tau["mass_fit_r2"])
            recommended_mass_slope = float(recommended_tau["mass_slope"])
            recommended_pred_um2_s_at_mass_0p1 = float(recommended_tau["pred_um2_s_at_0p1"])
            estimated_mass_scale = float(recommended_tau["required_mass_scale"])

        rows.append(
            {
                "temperature_up": temperature_up,
                "temperature_k": temperature_k,
                "temperature_c": temperature_c,
                "target_um2_s": target_um2_s,
                "target_status": status,
                "nearest_stable_mass_scale": float(nearest["mass_scale"]),
                "nearest_stable_damping": float(nearest["damping"]),
                "nearest_stable_um2_s": float(nearest["diffusion_mean_nm2_per_time_mean_um2_per_s"]),
                "nearest_stable_ratio_sim_over_target": float(
                    nearest["diffusion_mean_nm2_per_time_mean_um2_per_s"] / target_um2_s
                ),
                "nearest_stable_cv": float(nearest["diffusion_nm2_per_time_cv"]),
                "recommended_damping": recommended_damping,
                "recommended_mass_scale": recommended_mass_scale,
                "estimated_mass_scale_raw": estimated_mass_scale,
                "recommended_mass_fit_r2": recommended_mass_fit_r2,
                "recommended_mass_slope": recommended_mass_slope,
                "recommended_pred_um2_s_at_mass_0p1": recommended_pred_um2_s_at_mass_0p1,
                "note": note,
            }
        )
    return pd.DataFrame(rows).sort_values("temperature_up").reset_index(drop=True), arrhenius


def plot_corrections(table: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)

    axes[0].plot(
        table["temperature_up"],
        table["target_um2_s"],
        marker="o",
        linewidth=2.0,
        color="#0b3c5d",
        label="Target DOPC diffusion",
    )
    axes[0].plot(
        table["temperature_up"],
        table["nearest_stable_um2_s"],
        marker="o",
        linewidth=2.0,
        color="#328cc1",
        label="Closest stable scanned condition",
    )
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Upside temperature")
    axes[0].set_ylabel("Diffusion (um^2/s)")
    axes[0].set_title("Target diffusion versus current scan")
    axes[0].grid(alpha=0.25, which="both")
    axes[0].legend(loc="upper left", fontsize=8)

    finite_mask = np.isfinite(table["recommended_mass_scale"].to_numpy(dtype=float))
    if finite_mask.any():
        axes[1].plot(
            table.loc[finite_mask, "temperature_up"],
            table.loc[finite_mask, "recommended_mass_scale"],
            marker="o",
            linewidth=2.0,
            color="#b82601",
            label="Resolved mass correction",
        )
    unresolved_mask = ~finite_mask
    if unresolved_mask.any():
        axes[1].scatter(
            table.loc[unresolved_mask, "temperature_up"],
            np.full(int(unresolved_mask.sum()), 0.1),
            marker="x",
            s=70,
            linewidths=2.0,
            color="#b82601",
            label="Unresolved below scan range",
        )
    for row in table.itertuples(index=False):
        if math.isfinite(row.recommended_mass_scale):
            axes[1].annotate(
                f"tau={row.recommended_damping:g}",
                (row.temperature_up, row.recommended_mass_scale),
                textcoords="offset points",
                xytext=(0, 6),
                ha="center",
                fontsize=8,
            )
        else:
            axes[1].annotate(
                f"tau={row.recommended_damping:g}\nunresolved",
                (row.temperature_up, 0.1),
                textcoords="offset points",
                xytext=(0, 8),
                ha="center",
                fontsize=8,
            )
    axes[1].axhline(0.1, color="#666666", linestyle="--", linewidth=1.2, label="Current scan minimum mass")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Upside temperature")
    axes[1].set_ylabel("Recommended mass scale")
    axes[1].set_title("Mass correction needed at stable damping")
    axes[1].grid(alpha=0.25, which="both")
    axes[1].legend(loc="lower left", fontsize=8)

    fig.suptitle("Temperature-specific dryMARTINI correction at fixed Upside temperature", fontsize=13)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def write_markdown(
    path: Path,
    table: pd.DataFrame,
    arrhenius: ArrheniusModel,
    stable_dampings: list[float],
    minimum_reportable_mass_scale: float,
    max_cv_threshold: float,
) -> None:
    lines = [
        "# Temperature-Specific Bilayer Correction",
        "",
        "## Target Model",
        f"- Upside temperature conversion: `1.0 T_up = 350.588235 K` from `AGENTS.md`.",
        f"- DOPC diffusion target source: {BAG2014_CITATION}",
        "- Source table values used:",
        "  - 298 K -> 2.46 um^2/s",
        "  - 303 K -> 2.77 um^2/s",
        "  - 308 K -> 3.07 um^2/s",
        "  - 313 K -> 3.40 um^2/s",
        f"- Inference from those source values: Arrhenius fit gives `E_A = {arrhenius.activation_energy_kj_mol:.2f} kJ/mol`.",
        f"- Stable damping window used for correction selection: {', '.join(f'{value:g}' for value in stable_dampings)}.",
        "",
        "## Recommendation",
        f"- Robust subset used for correction selection: completed replicates, `R^2 >= 0.95`, `CV <= {max_cv_threshold:.2f}`, and continuity-stable mass progression.",
        "- Damping is not the main correction knob in this scan; use it only to stay in a stable window.",
        "- Mass scaling is the main correction. The table below shows the mass scale needed to hit the DOPC target at each fixed Upside temperature.",
        "",
        "| T_up | K | target D (um^2/s) | nearest robust in-grid (m, tau, D) | recommended tau | recommended mass | status |",
        "| --- | ---: | ---: | --- | ---: | ---: | --- |",
    ]
    for row in table.itertuples(index=False):
        recommended_mass_label = (
            f"{row.recommended_mass_scale:.4f}"
            if math.isfinite(row.recommended_mass_scale)
            else f"unresolved (< {minimum_reportable_mass_scale:g})"
        )
        lines.append(
            f"| {row.temperature_up:.1f} | {row.temperature_k:.1f} | {row.target_um2_s:.3f} | "
            f"({row.nearest_stable_mass_scale:.2g}, {row.nearest_stable_damping:g}, {row.nearest_stable_um2_s:.3f}) | "
            f"{row.recommended_damping:g} | {recommended_mass_label} | {row.target_status} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "- If the reported mass is `unresolved`, the refreshed scan does not contain a robust branch that reaches the physical target at that temperature.",
            "- If the recommended mass is below `0.1`, the current scan did not reach the physical target at that temperature.",
            "- `low_temperature_extrapolation` and `high_temperature_extrapolation` mean the target comes from extrapolating the experimental Arrhenius fit outside the source measurement window `298-313 K`.",
            "- In this joint scan, target-like diffusion appears only on higher-variance or discontinuous branches; those are intentionally excluded from the recommendation table.",
            "- The recommended damping values should be read as a stable working choice, not as a strong physical calibration knob by themselves.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    analysis_dir = args.analysis_dir.resolve()
    report_dir = args.report_dir.resolve() if args.report_dir else (analysis_dir / "report")
    report_dir.mkdir(parents=True, exist_ok=True)

    condition = load_condition_summary(report_dir)
    requested_stable_dampings = {float(value) for value in args.stable_damping_values.split(",") if value.strip()}
    stable_dampings = resolve_stable_dampings(condition, requested_stable_dampings)
    table, arrhenius = build_correction_table(
        condition,
        upside_temperature_k=args.upside_temperature_k,
        stable_dampings=set(stable_dampings),
        fit_r2_threshold=args.fit_r2_threshold,
        max_cv_threshold=args.max_cv_threshold,
        minimum_reportable_mass_scale=args.minimum_reportable_mass_scale,
        minimum_useful_mass_slope=args.minimum_useful_mass_slope,
    )

    csv_path = report_dir / "temperature_specific_corrections.csv"
    md_path = report_dir / "temperature_specific_correction.md"
    plot_path = report_dir / "temperature_specific_correction.png"
    table.to_csv(csv_path, index=False)
    plot_corrections(table, plot_path)
    write_markdown(
        md_path,
        table,
        arrhenius,
        stable_dampings,
        minimum_reportable_mass_scale=args.minimum_reportable_mass_scale,
        max_cv_threshold=args.max_cv_threshold,
    )

    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")
    print(f"Wrote {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
