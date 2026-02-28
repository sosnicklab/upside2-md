#!/usr/bin/env python3

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
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
        help=(
            "Output JSON manifest for relaxed systems. Defaults to "
            "<output-root>/relaxed_training_manifest.json"
        ),
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
    parser.add_argument(
        "--runner",
        choices=["local", "slurm"],
        default="local",
        help="Execution backend. 'slurm' submits one job per system.",
    )
    parser.add_argument(
        "--dry-run",
        type=int,
        choices=[0, 1],
        default=0,
        help="In slurm mode, write wrapper scripts and print sbatch commands without submitting.",
    )
    parser.add_argument(
        "--slurm-wait",
        type=int,
        choices=[0, 1],
        default=1,
        help="In slurm mode, wait for all submitted jobs before finalizing the manifest.",
    )
    parser.add_argument(
        "--slurm-poll-seconds",
        type=int,
        default=30,
        help="Polling interval while waiting for Slurm jobs to finish.",
    )
    parser.add_argument(
        "--slurm-account",
        default="",
        help="Optional Slurm account passed to sbatch.",
    )
    parser.add_argument(
        "--slurm-partition",
        default="",
        help="Optional Slurm partition passed to sbatch.",
    )
    parser.add_argument(
        "--slurm-time",
        default="24:00:00",
        help="Slurm walltime passed to sbatch.",
    )
    parser.add_argument(
        "--slurm-cpus-per-task",
        type=int,
        default=1,
        help="Slurm cpus-per-task value for each submitted run.",
    )
    parser.add_argument(
        "--slurm-mem",
        default="",
        help="Optional Slurm memory request, for example 8G.",
    )
    parser.add_argument(
        "--slurm-job-name-prefix",
        default="martini-relax",
        help="Prefix used for per-system Slurm job names.",
    )
    parser.add_argument(
        "--slurm-extra-arg",
        action="append",
        default=[],
        help="Additional raw sbatch argument. May be passed multiple times.",
    )
    return parser.parse_args()


def parse_key_value_overrides(items):
    out = {}
    for raw in items:
        if "=" not in raw:
            raise ValueError(f"Invalid --extra-env value (expected KEY=VALUE): {raw}")
        key, value = raw.split("=", 1)
        key = key.strip()
        if not key:
            raise ValueError(f"Invalid --extra-env key: {raw}")
        out[key] = value
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
    path = Path(raw).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Missing path for '{key}': {path}")
    return path


def stage_file_path(output_root, pdb_id, stage):
    return (output_root / pdb_id / "checkpoints" / f"{pdb_id}.stage_{stage}.up").resolve()


def normalize_pdb_id(entry):
    pdb_id = str(entry.get("pdb_id", "")).strip().lower()
    if not pdb_id:
        raise ValueError("missing pdb_id")
    return pdb_id


def build_env(base_env, entry, output_root, args, extra_env):
    pdb_id = normalize_pdb_id(entry)
    aa_pdb = Path(
        entry.get("aa_pdb_effective") or entry.get("aa_pdb_downloaded", "")
    ).expanduser().resolve()
    if not aa_pdb.exists():
        raise FileNotFoundError(f"{pdb_id}: AA PDB not found: {aa_pdb}")
    cg_pdb = require_path(entry, "cg_pdb")
    itp_file = require_path(entry, "itp_file")

    run_dir = (output_root / pdb_id).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    env_overrides = {
        "RUN_DIR": str(run_dir),
        "PROTEIN_AA_PDB": str(aa_pdb),
        "PROTEIN_CG_PDB": str(cg_pdb),
        "PROTEIN_ITP": str(itp_file),
        "MARTINIZE_ENABLE": "0",
        "SEED": str(int(args.seed)),
    }
    env_overrides.update(extra_env)

    env = dict(base_env)
    env.update(env_overrides)
    return env, env_overrides, run_dir


def build_skip_result(pdb_id, run_dir, target_file, batch_log, args):
    return {
        "pdb_id": pdb_id,
        "status": "skipped_existing",
        "run_dir": str(run_dir),
        "relaxed_stage": args.target_stage,
        "relaxed_up_file": str(target_file),
        "batch_log": str(batch_log),
    }


def prepare_run(entry, run_script, output_root, args, extra_env):
    pdb_id = normalize_pdb_id(entry)
    target_file = stage_file_path(output_root, pdb_id, args.target_stage)
    run_dir = (output_root / pdb_id).resolve()
    batch_log = run_dir / "batch_relax.log"
    if args.skip_existing and target_file.exists():
        return build_skip_result(pdb_id, run_dir, target_file, batch_log, args), None

    env, env_overrides, run_dir = build_env(os.environ, entry, output_root, args, extra_env)
    run_dir.mkdir(parents=True, exist_ok=True)
    return None, {
        "pdb_id": pdb_id,
        "env": env,
        "env_overrides": env_overrides,
        "run_dir": run_dir,
        "batch_log": batch_log,
        "target_file": target_file,
        "run_script": run_script,
    }


def merge_result(entry, result):
    merged = dict(entry)
    if "up_file" in merged and "prepared_up_file" not in merged:
        merged["prepared_up_file"] = merged["up_file"]
    merged["up_file"] = result["relaxed_up_file"]
    merged.update(result)
    return merged


def run_one_local(prepared, args):
    cmd = ["bash", str(prepared["run_script"]), f"PDB_ID={prepared['pdb_id']}"]
    with prepared["batch_log"].open("w", encoding="utf-8") as log_fh:
        log_fh.write(f"# runner: local\n")
        log_fh.write(f"# cmd: {' '.join(cmd)}\n")
        log_fh.write(f"# run_at: {datetime.now(timezone.utc).isoformat()}\n")
        proc = subprocess.run(
            cmd,
            cwd=str(prepared["run_script"].parent),
            env=prepared["env"],
            text=True,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            check=False,
        )

    if proc.returncode != 0:
        raise RuntimeError(
            f"{prepared['pdb_id']}: rigid dry run failed with exit {proc.returncode}; "
            f"see {prepared['batch_log']}"
        )
    if not prepared["target_file"].exists():
        raise FileNotFoundError(
            f"{prepared['pdb_id']}: target relaxed stage file missing after run: "
            f"{prepared['target_file']}"
        )

    return {
        "pdb_id": prepared["pdb_id"],
        "status": "ok",
        "run_dir": str(prepared["run_dir"]),
        "relaxed_stage": args.target_stage,
        "relaxed_up_file": str(prepared["target_file"]),
        "batch_log": str(prepared["batch_log"]),
    }


def require_command(name):
    path = shutil.which(name)
    if not path:
        raise FileNotFoundError(f"Required command not found in PATH: {name}")
    return path


def slurm_job_name(pdb_id, prefix):
    token = "".join(ch if ch.isalnum() or ch in "-_" else "-" for ch in pdb_id)
    return f"{prefix}-{token}"[:128]


def build_slurm_script(prepared):
    slurm_script = prepared["run_dir"] / "submit_relax.slurm.sh"
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        f"cd {shlex.quote(str(prepared['run_script'].parent))}",
    ]
    for key in sorted(prepared["env_overrides"]):
        lines.append(f"export {key}={shlex.quote(prepared['env_overrides'][key])}")
    lines.append(
        f"bash {shlex.quote(str(prepared['run_script']))} "
        f"PDB_ID={shlex.quote(prepared['pdb_id'])}"
    )
    slurm_script.write_text("\n".join(lines) + "\n", encoding="utf-8")
    slurm_script.chmod(0o755)
    return slurm_script


def build_sbatch_command(prepared, slurm_script, args):
    slurm_output = prepared["run_dir"] / "slurm-%j.out"
    cmd = [
        "sbatch",
        "--parsable",
        "--job-name",
        slurm_job_name(prepared["pdb_id"], args.slurm_job_name_prefix),
        "--output",
        str(slurm_output),
        "--chdir",
        str(prepared["run_script"].parent),
        "--nodes",
        "1",
        "--ntasks",
        "1",
        "--cpus-per-task",
        str(int(args.slurm_cpus_per_task)),
        "--time",
        str(args.slurm_time),
    ]
    if args.slurm_account:
        cmd.extend(["--account", args.slurm_account])
    if args.slurm_partition:
        cmd.extend(["--partition", args.slurm_partition])
    if args.slurm_mem:
        cmd.extend(["--mem", args.slurm_mem])
    for raw in args.slurm_extra_arg:
        pieces = shlex.split(raw)
        if not pieces:
            continue
        cmd.extend(pieces)
    cmd.append(str(slurm_script))
    return cmd, slurm_output


def submit_one_slurm(prepared, args):
    slurm_script = build_slurm_script(prepared)
    cmd, slurm_output = build_sbatch_command(prepared, slurm_script, args)
    with prepared["batch_log"].open("w", encoding="utf-8") as log_fh:
        log_fh.write("# runner: slurm\n")
        log_fh.write(f"# submit_at: {datetime.now(timezone.utc).isoformat()}\n")
        log_fh.write(f"# cmd: {' '.join(shlex.quote(x) for x in cmd)}\n")
        log_fh.write(f"# slurm_script: {slurm_script}\n")
        log_fh.write(f"# slurm_output: {slurm_output}\n")

    if args.dry_run:
        return {
            "pdb_id": prepared["pdb_id"],
            "status": "dry_run",
            "run_dir": str(prepared["run_dir"]),
            "relaxed_stage": args.target_stage,
            "relaxed_up_file": str(prepared["target_file"]),
            "batch_log": str(prepared["batch_log"]),
            "slurm_script": str(slurm_script),
            "slurm_output": str(slurm_output),
            "slurm_submit_cmd": cmd,
            "slurm_job_id": None,
        }

    require_command("sbatch")
    proc = subprocess.run(
        cmd,
        cwd=str(prepared["run_script"].parent),
        text=True,
        capture_output=True,
        check=False,
    )
    with prepared["batch_log"].open("a", encoding="utf-8") as log_fh:
        if proc.stdout:
            log_fh.write("# sbatch_stdout:\n")
            log_fh.write(proc.stdout)
            if not proc.stdout.endswith("\n"):
                log_fh.write("\n")
        if proc.stderr:
            log_fh.write("# sbatch_stderr:\n")
            log_fh.write(proc.stderr)
            if not proc.stderr.endswith("\n"):
                log_fh.write("\n")

    if proc.returncode != 0:
        raise RuntimeError(
            f"{prepared['pdb_id']}: sbatch submission failed with exit {proc.returncode}; "
            f"see {prepared['batch_log']}"
        )
    job_id = proc.stdout.strip().split(";", 1)[0]
    if not job_id:
        raise RuntimeError(
            f"{prepared['pdb_id']}: sbatch did not return a job id; see {prepared['batch_log']}"
        )

    return {
        "pdb_id": prepared["pdb_id"],
        "status": "submitted",
        "run_dir": str(prepared["run_dir"]),
        "relaxed_stage": args.target_stage,
        "relaxed_up_file": str(prepared["target_file"]),
        "batch_log": str(prepared["batch_log"]),
        "slurm_script": str(slurm_script),
        "slurm_output": str(slurm_output),
        "slurm_submit_cmd": cmd,
        "slurm_job_id": job_id,
    }


def summarize_active_states(status_by_job):
    counts = {}
    for state in status_by_job.values():
        counts[state] = counts.get(state, 0) + 1
    parts = [f"{state}={counts[state]}" for state in sorted(counts)]
    return ", ".join(parts)


def wait_for_slurm_jobs(submitted, poll_seconds):
    require_command("squeue")
    remaining = {str(item["slurm_job_id"]): item for item in submitted if item.get("slurm_job_id")}
    if not remaining:
        return

    while remaining:
        query_ids = ",".join(sorted(remaining))
        proc = subprocess.run(
            ["squeue", "-h", "-j", query_ids, "-o", "%A|%T"],
            text=True,
            capture_output=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"squeue failed with exit {proc.returncode}: {proc.stderr.strip()}"
            )

        active = {}
        for raw in proc.stdout.splitlines():
            line = raw.strip()
            if not line:
                continue
            parts = line.split("|", 1)
            if len(parts) != 2:
                continue
            job_id, state = parts[0].strip(), parts[1].strip()
            if job_id:
                active[job_id] = state or "UNKNOWN"

        if not active:
            break

        print(
            f"Waiting on {len(active)} Slurm job(s): {summarize_active_states(active)}"
        )
        time.sleep(max(1, int(poll_seconds)))
        remaining = {job_id: remaining[job_id] for job_id in active}


def finalize_slurm_submission(submission, args):
    target_file = Path(submission["relaxed_up_file"]).expanduser().resolve()
    if not target_file.exists():
        raise FileNotFoundError(
            f"{submission['pdb_id']}: target relaxed stage file missing after Slurm run: "
            f"{target_file}"
        )
    result = dict(submission)
    result["status"] = "ok"
    result["relaxed_stage"] = args.target_stage
    result["relaxed_up_file"] = str(target_file)
    return result


def write_manifest(relaxed_manifest, payload):
    relaxed_manifest.parent.mkdir(parents=True, exist_ok=True)
    with relaxed_manifest.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")


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
    if args.runner == "slurm":
        if not args.dry_run:
            require_command("sbatch")
        if args.slurm_wait and not args.dry_run:
            require_command("squeue")

    input_manifest, systems = load_manifest(manifest_path)
    only = {x.strip().lower() for x in args.only_pdb if x.strip()}
    if only:
        systems = [
            entry
            for entry in systems
            if str(entry.get("pdb_id", "")).strip().lower() in only
        ]
    if args.limit > 0:
        systems = systems[: args.limit]

    output_root.mkdir(parents=True, exist_ok=True)
    extra_env = parse_key_value_overrides(args.extra_env)

    relaxed_systems = []
    submitted = []
    failed = []

    total = len(systems)
    print(f"Systems selected: {total}")
    print(f"Runner: {args.runner}")
    for index, entry in enumerate(systems, start=1):
        try:
            pdb_id = normalize_pdb_id(entry)
        except Exception as exc:
            failed.append({"pdb_id": "", "reason": str(exc)})
            print(f"[{index}/{total}] FAIL <missing>: {exc}")
            continue

        try:
            skipped, prepared = prepare_run(entry, run_script, output_root, args, extra_env)
            if skipped is not None:
                relaxed_systems.append(merge_result(entry, skipped))
                print(
                    f"[{index}/{total}] {skipped['status'].upper()} {pdb_id}: "
                    f"{skipped['relaxed_up_file']}"
                )
                continue

            if args.runner == "local":
                result = run_one_local(prepared, args)
                relaxed_systems.append(merge_result(entry, result))
                print(
                    f"[{index}/{total}] {result['status'].upper()} {pdb_id}: "
                    f"{result['relaxed_up_file']}"
                )
                continue

            submission = submit_one_slurm(prepared, args)
            submitted.append({"entry": entry, "submission": submission})
            if submission["status"] == "dry_run":
                print(
                    f"[{index}/{total}] DRY_RUN {pdb_id}: "
                    f"{' '.join(shlex.quote(x) for x in submission['slurm_submit_cmd'])}"
                )
            else:
                print(
                    f"[{index}/{total}] SUBMITTED {pdb_id}: job {submission['slurm_job_id']}"
                )
        except Exception as exc:
            failed.append({"pdb_id": pdb_id, "reason": str(exc)})
            print(f"[{index}/{total}] FAIL {pdb_id}: {exc}")

    if args.runner == "slurm" and args.slurm_wait and not args.dry_run and submitted:
        wait_for_slurm_jobs(
            [item["submission"] for item in submitted],
            args.slurm_poll_seconds,
        )
        completed = []
        for item in submitted:
            try:
                result = finalize_slurm_submission(item["submission"], args)
                relaxed_systems.append(merge_result(item["entry"], result))
                completed.append(result["pdb_id"])
            except Exception as exc:
                failed.append({"pdb_id": item["submission"]["pdb_id"], "reason": str(exc)})
        for pdb_id in completed:
            print(f"[DONE] OK {pdb_id}")

    payload = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_manifest": str(manifest_path),
        "run_script": str(run_script),
        "output_root": str(output_root),
        "target_stage": args.target_stage,
        "seed": int(args.seed),
        "runner": args.runner,
        "dry_run": int(args.dry_run),
        "slurm_wait": int(args.slurm_wait),
        "input_system_count": int(len(input_manifest.get("systems", []))),
        "selected_system_count": int(total),
        "systems": relaxed_systems,
        "submitted": [item["submission"] for item in submitted],
        "failed": failed,
    }

    write_manifest(relaxed_manifest, payload)
    print(f"Relaxed manifest written: {relaxed_manifest}")
    print(f"Relaxed systems: {len(relaxed_systems)}")
    print(f"Submitted jobs: {len(submitted)}")
    print(f"Failed systems: {len(failed)}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
