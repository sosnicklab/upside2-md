#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def parse_args():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=(
            "Batch-run rigid dry-MARTINI relaxation for systems listed in a "
            "training manifest, and emit a relaxed manifest for RBM training."
        )
    )
    parser.add_argument(
        "--manifest",
        default="train-data/training_manifest.json",
        help="Input manifest from download_and_prepare_training_inputs.py",
    )
    parser.add_argument(
        "--run-script",
        default=str(script_dir / "run_relax_6x_rigid_dry.sh"),
        help="Path to standalone rigid dry relaxation script (stages 6.0-6.6).",
    )
    parser.add_argument(
        "--output-root",
        default="train-data/relax_runs",
        help="Root directory for per-system relaxation outputs.",
    )
    parser.add_argument(
        "--relaxed-manifest-out",
        default=None,
        help="Output JSON manifest for relaxed systems. Defaults to <output-root>/relaxed_training_manifest.json",
    )
    parser.add_argument(
        "--target-stage",
        default="6.6",
        choices=["6.6"],
        help="Stage file to expose as training input in relaxed manifest.",
    )
    parser.add_argument(
        "--skip-existing",
        type=int,
        choices=[0, 1],
        default=1,
        help="Reuse existing target-stage outputs if present.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Process at most this many systems (0 = all).",
    )
    parser.add_argument(
        "--only-pdb",
        nargs="*",
        default=[],
        help="Optional PDB IDs to include.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=7090685331,
        help="Seed passed to workflow script.",
    )
    parser.add_argument(
        "--extra-env",
        action="append",
        default=[],
        help="Extra KEY=VALUE env override; may be passed multiple times.",
    )
    return parser.parse_args()


def parse_key_value_overrides(items):
    out = {}
    for raw in items:
        if "=" not in raw:
            raise ValueError(f"Invalid --extra-env value (expected KEY=VALUE): {raw}")
        k, v = raw.split("=", 1)
        k = k.strip()
        if not k:
            raise ValueError(f"Invalid --extra-env key: {raw}")
        out[k] = v
    return out


def load_manifest(path):
    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"Manifest must be a JSON object: {path}")
    systems = data.get("systems", [])
    if not isinstance(systems, list):
        raise ValueError(f"Manifest field 'systems' must be a list: {path}")
    return data, systems


def require_path(entry, key):
    raw = entry.get(key)
    if not raw:
        raise ValueError(f"Manifest entry missing '{key}'")
    p = Path(raw).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"Missing path for '{key}': {p}")
    return p


def stage_file_path(output_root, pdb_id, stage):
    return (output_root / pdb_id / "checkpoints" / f"{pdb_id}.stage_{stage}.up").resolve()


def build_env(base_env, entry, output_root, args, extra_env):
    pdb_id = entry["pdb_id"]
    aa_pdb = Path(entry.get("aa_pdb_effective") or entry.get("aa_pdb_downloaded", "")).expanduser().resolve()
    if not aa_pdb.exists():
        raise FileNotFoundError(f"{pdb_id}: AA PDB not found: {aa_pdb}")
    cg_pdb = require_path(entry, "cg_pdb")
    itp_file = require_path(entry, "itp_file")

    run_dir = (output_root / pdb_id).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    env = dict(base_env)
    env["RUN_DIR"] = str(run_dir)
    env["PROTEIN_AA_PDB"] = str(aa_pdb)
    env["PROTEIN_CG_PDB"] = str(cg_pdb)
    env["PROTEIN_ITP"] = str(itp_file)
    env["MARTINIZE_ENABLE"] = "0"
    env["SEED"] = str(int(args.seed))
    env.update(extra_env)
    return env, run_dir


def run_one(entry, run_script, output_root, args, extra_env):
    pdb_id = str(entry["pdb_id"]).strip().lower()
    if not pdb_id:
        raise ValueError("Empty pdb_id in manifest entry")

    target_file = stage_file_path(output_root, pdb_id, args.target_stage)
    run_dir = (output_root / pdb_id).resolve()
    batch_log = run_dir / "batch_relax.log"
    if args.skip_existing and target_file.exists():
        return {
            "pdb_id": pdb_id,
            "status": "skipped_existing",
            "run_dir": str(run_dir),
            "relaxed_stage": args.target_stage,
            "relaxed_up_file": str(target_file),
            "batch_log": str(batch_log),
        }

    env, run_dir = build_env(os.environ, entry, output_root, args, extra_env)
    batch_log.parent.mkdir(parents=True, exist_ok=True)

    cmd = ["bash", str(run_script), f"PDB_ID={pdb_id}"]
    with batch_log.open("w", encoding="utf-8") as log_fh:
        log_fh.write(f"# cmd: {' '.join(cmd)}\n")
        log_fh.write(f"# run_at: {datetime.now(timezone.utc).isoformat()}\n")
        proc = subprocess.run(
            cmd,
            cwd=str(run_script.parent),
            env=env,
            text=True,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            check=False,
        )

    if proc.returncode != 0:
        raise RuntimeError(
            f"{pdb_id}: rigid dry run failed with exit {proc.returncode}; see {batch_log}"
        )
    if not target_file.exists():
        raise FileNotFoundError(
            f"{pdb_id}: target relaxed stage file missing after run: {target_file}"
        )

    return {
        "pdb_id": pdb_id,
        "status": "ok",
        "run_dir": str(run_dir),
        "relaxed_stage": args.target_stage,
        "relaxed_up_file": str(target_file),
        "batch_log": str(batch_log),
    }


def main():
    args = parse_args()
    manifest_path = Path(args.manifest).expanduser().resolve()
    run_script = Path(args.run_script).expanduser().resolve()
    output_root = Path(args.output_root).expanduser().resolve()
    relaxed_manifest = (
        Path(args.relaxed_manifest_out).expanduser().resolve()
        if args.relaxed_manifest_out
        else (output_root / "relaxed_training_manifest.json")
    )

    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")
    if not run_script.exists():
        raise FileNotFoundError(f"Run script not found: {run_script}")

    input_manifest, systems = load_manifest(manifest_path)
    only = {x.strip().lower() for x in args.only_pdb if x.strip()}
    if only:
        systems = [s for s in systems if str(s.get("pdb_id", "")).strip().lower() in only]
    if args.limit > 0:
        systems = systems[: args.limit]

    output_root.mkdir(parents=True, exist_ok=True)
    extra_env = parse_key_value_overrides(args.extra_env)

    relaxed_systems = []
    failed = []

    total = len(systems)
    print(f"Systems selected: {total}")
    for i, entry in enumerate(systems, start=1):
        pdb_id = str(entry.get("pdb_id", "")).strip().lower()
        if not pdb_id:
            failed.append({"pdb_id": "", "reason": "missing pdb_id"})
            print(f"[{i}/{total}] FAIL <missing>: missing pdb_id")
            continue
        try:
            result = run_one(entry, run_script, output_root, args, extra_env)
            merged = dict(entry)
            if "up_file" in merged and "prepared_up_file" not in merged:
                merged["prepared_up_file"] = merged["up_file"]
            merged["up_file"] = result["relaxed_up_file"]
            merged.update(result)
            relaxed_systems.append(merged)
            print(f"[{i}/{total}] {result['status'].upper()} {pdb_id}: {result['relaxed_up_file']}")
        except Exception as exc:
            failed.append({"pdb_id": pdb_id, "reason": str(exc)})
            print(f"[{i}/{total}] FAIL {pdb_id}: {exc}")

    out = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_manifest": str(manifest_path),
        "run_script": str(run_script),
        "output_root": str(output_root),
        "target_stage": args.target_stage,
        "seed": int(args.seed),
        "input_system_count": int(len(input_manifest.get("systems", []))),
        "selected_system_count": int(total),
        "systems": relaxed_systems,
        "failed": failed,
    }

    relaxed_manifest.parent.mkdir(parents=True, exist_ok=True)
    with relaxed_manifest.open("w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(f"Relaxed manifest written: {relaxed_manifest}")
    print(f"Relaxed systems: {len(relaxed_systems)}")
    print(f"Failed systems: {len(failed)}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
