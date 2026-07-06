#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import py.martini_prepare_system as martini_prepare_system
from py.martini_prepare_system_lib import inject_cg_lipid_nodes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--stock-run-dir", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--dopc-h5", required=True)
    parser.add_argument("--seed", required=True, type=int)
    args = parser.parse_args()

    project_root = Path(args.project_root).expanduser().resolve()
    stock_run_dir = Path(args.stock_run_dir).expanduser().resolve()
    run_dir = Path(args.run_dir).expanduser().resolve()
    dopc_h5 = Path(args.dopc_h5).expanduser().resolve()
    pdb_id = args.pdb_id
    seed = int(args.seed)

    workflow_args = martini_prepare_system.parse_hybrid_workflow_args(
        [
            "--pdb-id",
            pdb_id,
            "--runtime-pdb-id",
            f"{pdb_id}_hybrid",
            "--upside-home",
            str(project_root),
            "--run-dir",
            str(run_dir),
            "--protein-aa-pdb",
            str(project_root / "example" / "16.MARTINI" / "pdb" / f"{pdb_id}.pdb"),
            "--bilayer-pdb",
            str(project_root / "parameters" / "dryMARTINI" / "DOPC.pdb"),
            "--extract-vtf-script",
            str(project_root / "py" / "martini_extract_vtf.py"),
            "--lipid-resolution",
            "coarse",
            "--dopc-h5",
            str(dopc_h5),
            "--seed",
            str(seed),
            "--prep-seed",
            str(seed),
        ]
    )
    martini_prepare_system.validate_hybrid_workflow_args(workflow_args)
    files = martini_prepare_system.workflow_stage_files(workflow_args)

    stock_prepared = stock_run_dir / "checkpoints" / f"{pdb_id}.stage_7.0.prepared.up"
    stock_source = stock_run_dir / "checkpoints" / f"{pdb_id}.stage_6.0.up"
    shutil.copy2(stock_prepared, files["prepared_70"])
    inject_cg_lipid_nodes(files["prepared_70"], dopc_h5)

    shutil.copy2(files["prepared_70"], files["stage_70"])
    martini_prepare_system.handoff_initial_position(
        workflow_args,
        stock_source,
        files["stage_70"],
        "production_hybrid",
    )
    martini_prepare_system.assert_hybrid_stage_active(
        files["stage_70"],
        "production",
        "production",
    )
    martini_prepare_system.run_minimization_stage(
        workflow_args,
        "7.0",
        files["stage_70"],
        workflow_args.min_70_max_iter,
        preserve_stage=True,
    )
    if float(workflow_args.stage_70_burnin_protein_restraint_spring) > 0.0:
        martini_prepare_system.inject_protein_position_restraints(
            files["stage_70"],
            spring_const=float(workflow_args.stage_70_burnin_protein_restraint_spring),
        )
    martini_prepare_system.run_stage70_burnin(workflow_args, files["stage_70"])
    if (
        int(workflow_args.prod_70_burnin_nsteps) > 0
        and float(workflow_args.stage_70_burnin_protein_restraint_spring) > 0.0
        and (
            int(workflow_args.stage_70_release_sc_env_backbone_hold_steps) > 0
            or int(workflow_args.stage_70_release_sc_env_po4_z_hold_steps) > 0
        )
    ):
        martini_prepare_system.reset_stage70_release_hybrid_transition(
            files["stage_70"],
            workflow_args,
        )
    martini_prepare_system.remove_protein_position_restraints(files["stage_70"])
    martini_prepare_system.run_md_stage(
        workflow_args,
        "7.0",
        files["stage_70"],
        files["stage_70"],
        workflow_args.prod_70_nsteps,
        workflow_args.prod_time_step,
        workflow_args.prod_frame_steps,
    )
    martini_prepare_system.extract_stage_vtf(
        workflow_args,
        "7.0",
        files["stage_70"],
        "2",
    )


if __name__ == "__main__":
    main()
