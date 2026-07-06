#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import h5py
import numpy as np

import py.martini_build_tables as martini_build_tables
from martini_itp_reader import parse_dopc_from_itp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--representatives", type=int, default=1)
    args = parser.parse_args()

    base = Path(args.base).expanduser().resolve()
    out = Path(args.out).expanduser().resolve()
    representative_count = max(1, int(args.representatives))

    shutil.copy2(base, out)

    corrections = martini_build_tables._build_single_cgl_compaction_corrections_from_base_h5(
        base_h5_path=out,
        include_sc=False,
        include_target=False,
        include_pair_compaction=False,
    )

    defaults = martini_build_tables._default_dry_martini_repo_paths()
    dry_ff_path = defaults["dry_ff_path"]
    lipids_itp_path = defaults["lipids_itp_path"]
    _, pair_params = martini_build_tables.parse_dry_forcefield(dry_ff_path)
    dopc = parse_dopc_from_itp(lipids_itp_path)
    compaction_pool_conformers = martini_build_tables._positive_int_env(
        "UPSIDE_CGL_COMPACTION_POOL_CONFORMERS", 32
    )
    compaction_burnin_steps = martini_build_tables._positive_int_env(
        "UPSIDE_CGL_COMPACTION_BURNIN_STEPS", 20000
    )
    compaction_steps_per_conf = martini_build_tables._positive_int_env(
        "UPSIDE_CGL_COMPACTION_STEPS_PER_CONFORMER", 500
    )

    def load_or_build_isolated_refs() -> np.ndarray:
        cache_path = Path(
            f"/private/tmp/cgl_isolated_compaction_refs_{compaction_pool_conformers}_"
            f"{compaction_burnin_steps}_{compaction_steps_per_conf}.npy"
        )
        if cache_path.exists():
            return np.load(cache_path)
        refs = martini_build_tables.sample_isolated_dopc_bonded_conformers(
            dopc,
            lipids_itp_path=lipids_itp_path,
            pair_params=pair_params,
            conformer_count=compaction_pool_conformers,
            temperature_upside=martini_build_tables.DEFAULT_PRODUCTION_TEMP_UPSIDE,
            seed=1777,
            mc_burnin_steps=compaction_burnin_steps,
            mc_steps_per_conformer=compaction_steps_per_conf,
        )
        np.save(cache_path, refs)
        return refs

    with h5py.File(out, "r+") as h5:
        cg_grp = h5["cg_lipid_table"]
        martini_build_tables._apply_single_cgl_compaction_corrections_to_h5(
            cg_grp,
            corrections,
        )

        ref_dataset = martini_build_tables._select_cgl_compaction_reference_dataset_name(cg_grp, None)
        ref_ensemble_nm = np.asarray(cg_grp[ref_dataset][:], dtype=np.float64)
        base_states = martini_build_tables._reference_compaction_state_metadata_from_ensemble(
            ref_ensemble_nm,
            representative_count=representative_count,
        )
        compact_center_ang = float(base_states["compact_center_ang"])

        ref_compaction_values_ang = martini_build_tables._dopc_tail_extension_series_ang(ref_ensemble_nm)
        supplemental_refs_nm = None
        sample_compaction_support_ang = None
        comp_grp = cg_grp["cg_lipid_compaction"]
        if "target_compaction_ang" in comp_grp:
            sample_compaction_support_ang = np.asarray(
                comp_grp["target_compaction_ang"][:],
                dtype=np.float64,
            )
        if int(np.count_nonzero(ref_compaction_values_ang > compact_center_ang + 1.0e-6)) < max(
            2, representative_count
        ):
            extra_ref_pools = []
            for dataset_name in ("sc_interface_ref_bead_positions_nm", "ref_bead_positions_nm"):
                if dataset_name != ref_dataset and dataset_name in cg_grp:
                    extra_ref_pools.append(np.asarray(cg_grp[dataset_name][:], dtype=np.float64))
            extra_ref_pools.append(np.asarray(load_or_build_isolated_refs(), dtype=np.float64))
            if extra_ref_pools:
                supplemental_refs_nm = np.concatenate(extra_ref_pools, axis=0)

        compaction_states = martini_build_tables._augment_compaction_states_with_compressed_branch(
            ref_ensemble_nm,
            base_states,
            representative_count=representative_count,
            supplemental_refs_nm=supplemental_refs_nm,
            sample_compaction_ang=sample_compaction_support_ang,
        )
        compressed_center_ang = float(compaction_states["compressed_center_ang"])
        if not compressed_center_ang > compact_center_ang + 1.0e-6:
            raise RuntimeError("Failed to build a valid compressed-state representative set")

        pair_grp = cg_grp["cg_lipid_pair"]
        pair_energy_grid = np.asarray(pair_grp["energy_grid_raw_kj_mol"][:], dtype=np.float64)
        n_radial, n_angle1, n_angle2 = pair_energy_grid.shape
        if n_angle1 != n_angle2:
            raise RuntimeError("CGL pair retrofit requires a square angular grid")
        pair_result_cg = {
            "energy_grid_raw": pair_energy_grid,
            "reference_energy_eup": float(pair_grp["reference_energy_eup"][0, 0]),
            "r_values_nm": np.linspace(
                float(pair_grp.attrs["fit_r_min_nm"]),
                float(pair_grp.attrs["fit_r_max_nm"]),
                int(n_radial),
                dtype=np.float64,
            ),
            "cos_theta_grid": np.linspace(-1.0, 1.0, int(n_angle1), dtype=np.float64),
            "n_radial": int(pair_grp.attrs["n_radial"]),
            "n_angular": int(pair_grp.attrs["n_angular"]),
            "knot_spacing_ang": float(pair_grp.attrs["knot_spacing_ang"]),
            "fit_smooth": float(pair_grp.attrs.get("fit_smooth", 0.1)),
            "azimuthal_count": int(pair_grp.attrs.get("azimuthal_count", 1)),
            "bead_frame_count": int(pair_grp.attrs.get("cgl_bead_frame_count", 1)),
        }

        def sample_state_grid(refs_i_nm: np.ndarray, refs_j_nm: np.ndarray) -> np.ndarray:
            return martini_build_tables._sample_cg_pair_energy_grid(
                refs_i_nm,
                refs_j_nm,
                bead_types=list(dopc["bead_types"]),
                bead_charges=list(dopc["bead_charges"]),
                pair_params=pair_params,
                r_values_nm=pair_result_cg["r_values_nm"],
                cos_theta_grid=pair_result_cg["cos_theta_grid"],
                azimuthal_count=int(pair_result_cg["azimuthal_count"]),
                bead_frame_count=int(pair_result_cg["bead_frame_count"]),
                dist_min_nm=1.0e-6,
                average_temperature=martini_build_tables.DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
            )

        def sample_symmetric_state_grid(refs_i_nm: np.ndarray, refs_j_nm: np.ndarray) -> np.ndarray:
            forward = sample_state_grid(refs_i_nm, refs_j_nm)
            reverse = sample_state_grid(refs_j_nm, refs_i_nm)
            return 0.5 * (forward + np.swapaxes(reverse, 1, 2))

        def sample_relax_correction(refs_i_nm: np.ndarray, refs_j_nm: np.ndarray | None) -> np.ndarray:
            result = martini_build_tables._pair_conditioned_tail_relaxation_correction_grid(
                refs_i_nm,
                refs_j_nm,
                bead_types=list(dopc["bead_types"]),
                bead_charges=list(dopc["bead_charges"]),
                pair_params=pair_params,
                lipids_itp_path=lipids_itp_path,
                r_values_nm=pair_result_cg["r_values_nm"],
                cos_theta_grid=pair_result_cg["cos_theta_grid"],
                average_temperature=martini_build_tables.DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
                azimuthal_count=int(pair_result_cg["azimuthal_count"]),
                bead_frame_count=int(pair_result_cg["bead_frame_count"]),
                face_cos_min=0.5,
                radial_cutoff_nm=2.5,
                dist_min_nm=1.0e-6,
            )
            return np.asarray(result["correction_grid_kj_mol"], dtype=np.float64)

        def symmetric_relax_correction(refs_i_nm: np.ndarray, refs_j_nm: np.ndarray) -> np.ndarray:
            forward = sample_relax_correction(refs_i_nm, refs_j_nm)
            reverse = sample_relax_correction(refs_j_nm, refs_i_nm)
            return 0.5 * (forward + np.swapaxes(reverse, 1, 2))

        state_grid_xx = sample_state_grid(
            compaction_states["compressed_refs_nm"],
            compaction_states["compressed_refs_nm"],
        ) + sample_relax_correction(compaction_states["compressed_refs_nm"], None)
        state_grid_ex = sample_symmetric_state_grid(
            compaction_states["extended_refs_nm"],
            compaction_states["compressed_refs_nm"],
        ) + symmetric_relax_correction(
            compaction_states["extended_refs_nm"],
            compaction_states["compressed_refs_nm"],
        )
        state_grid_cx = sample_symmetric_state_grid(
            compaction_states["compact_refs_nm"],
            compaction_states["compressed_refs_nm"],
        ) + symmetric_relax_correction(
            compaction_states["compact_refs_nm"],
            compaction_states["compressed_refs_nm"],
        )

        reference_energy_kj_mol = (
            float(pair_result_cg["reference_energy_eup"])
            * martini_build_tables.ENERGY_CONVERSION_KJ_PER_EUP
        )
        control_kbt_kj_mol = (
            martini_build_tables.DEFAULT_PRODUCTION_TEMP_UPSIDE
            * martini_build_tables.ENERGY_CONVERSION_KJ_PER_EUP
        )
        if control_kbt_kj_mol <= 0.0:
            raise RuntimeError("Compressed pair retrofit requires positive control temperature")

        def to_log1p_control(grid_kj_mol: np.ndarray) -> np.ndarray:
            reduced = (np.asarray(grid_kj_mol, dtype=np.float64) - reference_energy_kj_mol) / control_kbt_kj_mol
            return np.log1p(np.maximum(reduced, 0.0))

        base_control_grid = to_log1p_control(pair_result_cg["energy_grid_raw"])
        fit_smooth = martini_build_tables._float_env(
            "UPSIDE_CGL_COMPACTION_FIT_SMOOTH",
            float(pair_result_cg.get("fit_smooth", 0.1)),
        )

        def fit_correction_tensor(correction_grid_control: np.ndarray) -> np.ndarray:
            tensor = martini_build_tables._fit_radial_angular_angular_tensor_bspline(
                pair_result_cg["r_values_nm"],
                pair_result_cg["cos_theta_grid"],
                correction_grid_control,
                n_knot_radial=int(pair_result_cg["n_radial"]),
                n_knot_angular=int(pair_result_cg["n_angular"]),
                knot_spacing_ang=float(pair_result_cg["knot_spacing_ang"]),
                energy_conversion=1.0,
                smooth=fit_smooth,
            )
            corr_min = float(np.min(correction_grid_control))
            corr_max = float(np.max(correction_grid_control))
            return np.clip(tensor, corr_min, corr_max).reshape(-1)

        correction_ex_control = to_log1p_control(state_grid_ex) - base_control_grid
        correction_cx_control = to_log1p_control(state_grid_cx) - base_control_grid
        correction_xx_control = to_log1p_control(state_grid_xx) - base_control_grid

        comp_grp = cg_grp["cg_lipid_compaction"]
        for dataset_name in (
            "delta_extended_compressed",
            "delta_compact_compressed",
            "delta_compressed_compressed",
            "grid_extended_compressed_kj_mol",
            "grid_compact_compressed_kj_mol",
            "grid_compressed_compressed_kj_mol",
        ):
            if dataset_name in comp_grp:
                del comp_grp[dataset_name]
        comp_grp.attrs["pair_reference_compressed_center_ang"] = np.float32(compressed_center_ang)
        comp_grp.attrs["correction_center_mode"] = "base"
        comp_grp.attrs["pair_state_model"] = "bilinear"
        comp_grp.create_dataset(
            "delta_extended_compressed",
            data=np.asarray(fit_correction_tensor(correction_ex_control), dtype=np.float32),
        )
        comp_grp.create_dataset(
            "delta_compact_compressed",
            data=np.asarray(fit_correction_tensor(correction_cx_control), dtype=np.float32),
        )
        comp_grp.create_dataset(
            "delta_compressed_compressed",
            data=np.asarray(fit_correction_tensor(correction_xx_control), dtype=np.float32),
        )
        comp_grp.create_dataset(
            "grid_extended_compressed_kj_mol",
            data=np.asarray(state_grid_ex, dtype=np.float32),
        )
        comp_grp.create_dataset(
            "grid_compact_compressed_kj_mol",
            data=np.asarray(state_grid_cx, dtype=np.float32),
        )
        comp_grp.create_dataset(
            "grid_compressed_compressed_kj_mol",
            data=np.asarray(state_grid_xx, dtype=np.float32),
        )

    print(
        f"Built exact compressed pair candidate at {out} "
        f"with representative_count={representative_count} "
        f"and compressed_center_ang={compressed_center_ang:.6f}"
    )


if __name__ == "__main__":
    main()
