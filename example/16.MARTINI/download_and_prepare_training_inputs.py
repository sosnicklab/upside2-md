#!/usr/bin/env python3

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def parse_args():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=(
            "Download protein inputs and run MARTINI preparation to generate "
            "training-ready .up files plus a manifest for RBM training."
        )
    )
    parser.add_argument(
        "--list-file",
        default="/Users/yinhan/Documents/Train/upside_input/list",
        help="File containing one PDB ID per line.",
    )
    parser.add_argument("--output-root", default="outputs/rbm_training_inputs")
    parser.add_argument("--manifest-out", default=None)
    parser.add_argument("--prepare-script", default=str(script_dir / "prepare_system.py"))
    parser.add_argument(
        "--prepare-python",
        default=None,
        help="Python executable for running prepare_system.py (auto-resolved if omitted).",
    )
    parser.add_argument("--bilayer-pdb", default=str(script_dir / "pdb" / "bilayer.MARTINI.pdb"))
    parser.add_argument("--protein-cg-template", default=str(script_dir / "pdb" / "{pdb_id}.MARTINI.pdb"))
    parser.add_argument("--protein-itp-template", default=str(script_dir / "pdb" / "{pdb_id}_proa.itp"))
    parser.add_argument("--auto-martinize", type=int, choices=[0, 1], default=1)
    parser.add_argument("--martinize-script", default=str(script_dir / "martinize.py"))
    parser.add_argument("--martinize-ff", default="martini22")
    parser.add_argument("--martinize-molname", default="PROA")
    parser.add_argument("--stage", default="minimization")
    parser.add_argument("--xy-scale", type=float, default=1.0)
    parser.add_argument("--box-padding-xy", type=float, default=0.0)
    parser.add_argument("--box-padding-z", type=float, default=20.0)
    parser.add_argument("--salt-molar", type=float, default=0.15)
    parser.add_argument("--protein-lipid-cutoff", type=float, default=3.0)
    parser.add_argument("--ion-cutoff", type=float, default=4.0)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--timeout", type=float, default=20.0)
    parser.add_argument("--skip-existing", type=int, choices=[0, 1], default=1)
    parser.add_argument("--allow-missing-cg", type=int, choices=[0, 1], default=1)
    parser.add_argument(
        "--aa-clean-mode",
        choices=["largest_complete_chain", "none"],
        default="largest_complete_chain",
        help="How to sanitize AA PDB before martinize.",
    )
    parser.add_argument(
        "--min-complete-residues",
        type=int,
        default=40,
        help="Minimum complete residues required after AA cleaning.",
    )
    return parser.parse_args()


def read_pdb_list(path):
    pdb_ids = []
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            token = line.split()[0].strip().lower()
            if not token:
                continue
            pdb_ids.append(token)
    return pdb_ids


def render_template(path_template, pdb_id):
    return Path(path_template.format(pdb_id=pdb_id, PDB_ID=pdb_id.upper())).expanduser().resolve()


def try_download_file(url, out_path, timeout):
    with urlopen(url, timeout=timeout) as response:
        payload = response.read()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_bytes(payload)


def download_aa_pdb(pdb_id, out_path, timeout):
    candidates = [
        f"https://opm-assets.storage.googleapis.com/pdb/{pdb_id.lower()}.pdb",
        f"https://opm-assets.storage.googleapis.com/pdb/{pdb_id.upper()}.pdb",
        f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb",
    ]
    errors = []
    for url in candidates:
        try:
            try_download_file(url, out_path, timeout)
            return url
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            errors.append(f"{url}: {exc}")
    raise RuntimeError("Failed download attempts:\n" + "\n".join(errors))


def parse_pdb_atom_line(line):
    return {
        "record": line[:6].strip(),
        "serial": int(line[6:11]),
        "name": line[12:16].strip(),
        "altloc": line[16:17],
        "resname": line[17:21].strip().upper(),
        "chain": line[21:22],
        "resseq": int(line[22:26]),
        "icode": line[26:27],
        "line": line.rstrip("\n"),
    }


def clean_aa_for_martinize(args, pdb_id, aa_pdb, output_root):
    if args.aa_clean_mode == "none":
        return aa_pdb.resolve(), {"mode": "none", "kept_residues": None}

    aa_clean_dir = (output_root / "aa_pdb_clean").resolve()
    aa_clean_dir.mkdir(parents=True, exist_ok=True)
    out_path = (aa_clean_dir / f"{pdb_id}.cleaned.pdb").resolve()

    aa_residues = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "HID",
        "HIE",
        "HIP",
        "HSD",
        "HSE",
        "HSP",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "CYX",
    }
    backbone_required = {"N", "CA", "C", "O"}
    residue_heavy_requirements = {
        "ALA": {"CB"},
        "ARG": {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
        "ASN": {"CB", "CG", "OD1", "ND2"},
        "ASP": {"CB", "CG", "OD1", "OD2"},
        "CYS": {"CB", "SG"},
        "CYX": {"CB", "SG"},
        "GLN": {"CB", "CG", "CD", "OE1", "NE2"},
        "GLU": {"CB", "CG", "CD", "OE1", "OE2"},
        "GLY": set(),
        "HIS": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
        "HID": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
        "HIE": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
        "HIP": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
        "HSD": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
        "HSE": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
        "HSP": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
        "ILE": {"CB", "CG1", "CG2", "CD1"},
        "LEU": {"CB", "CG", "CD1", "CD2"},
        "LYS": {"CB", "CG", "CD", "CE", "NZ"},
        "MET": {"CB", "CG", "SD", "CE"},
        "PHE": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
        "PRO": {"CB", "CG", "CD"},
        "SER": {"CB", "OG"},
        "THR": {"CB", "OG1", "CG2"},
        "TRP": {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
        "TYR": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
        "VAL": {"CB", "CG1", "CG2"},
    }

    residue_atom_alias = {
        "ILE": {"CD": "CD1"},
    }

    cryst1 = None
    segment_order = []
    segment_residue_order = defaultdict(list)
    residue_atoms = {}
    residue_names = {}

    seg_idx = 0
    with aa_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            if raw.startswith("CRYST1") and cryst1 is None:
                cryst1 = raw.rstrip("\n")
            if raw.startswith("TER"):
                seg_idx += 1
                continue
            if not raw.startswith(("ATOM", "HETATM")):
                continue

            atom = parse_pdb_atom_line(raw)
            if atom["record"] != "ATOM":
                continue
            if atom["resname"] not in aa_residues:
                continue
            if atom["altloc"].strip() not in {"", "A"}:
                continue

            chain = atom["chain"]
            res_key = (chain, atom["resseq"], atom["icode"], atom["resname"], seg_idx)
            if res_key not in residue_atoms:
                segment_order.append(res_key)
                segment_residue_order[(chain, seg_idx)].append(res_key)
                residue_atoms[res_key] = {}
                residue_names[res_key] = atom["resname"]
            residue_atoms[res_key][atom["name"].upper()] = atom

    best_chain_seg = None
    best_complete = []
    for chain_seg, res_keys in segment_residue_order.items():
        complete = []
        for rk in res_keys:
            names = set(residue_atoms[rk].keys())
            resname = residue_names[rk]
            for src, dst in residue_atom_alias.get(resname, {}).items():
                if src in names:
                    names.add(dst)
            if not backbone_required.issubset(names):
                continue
            needed = residue_heavy_requirements.get(resname, {"CB"})
            if not needed.issubset(names):
                continue
            complete.append(rk)
        if len(complete) > len(best_complete):
            best_complete = complete
            best_chain_seg = chain_seg

    if len(best_complete) < int(args.min_complete_residues):
        raise RuntimeError(
            f"{pdb_id}: cleaned AA chain has {len(best_complete)} complete residues, "
            f"below --min-complete-residues={args.min_complete_residues}"
        )

    # Rewrite as a single clean chain A with sequential residue numbering.
    serial = 1
    res_map = {rk: i + 1 for i, rk in enumerate(best_complete)}
    with out_path.open("w", encoding="utf-8") as out:
        if cryst1:
            out.write(cryst1 + "\n")
        for rk in best_complete:
            atoms_by_name = residue_atoms[rk]
            ordered_names = ["N", "CA", "C", "O"]
            ordered_names += sorted(n for n in atoms_by_name.keys() if n not in {"N", "CA", "C", "O"})
            for name in ordered_names:
                atom = atoms_by_name[name]
                line = atom["line"]
                new_line = (
                    f"{line[:6]}{serial:5d}{line[11:21]}A{res_map[rk]:4d}{line[26:]}"
                )
                out.write(new_line + "\n")
                serial += 1
        out.write("TER\nEND\n")

    meta = {
        "mode": "largest_complete_chain",
        "selected_chain": best_chain_seg[0] if best_chain_seg else "",
        "selected_segment_index": int(best_chain_seg[1]) if best_chain_seg else -1,
        "kept_residues": len(best_complete),
        "total_candidate_residues": len(segment_order),
    }
    return out_path, meta


def python_has_prepare_deps(python_exe):
    try:
        probe = subprocess.run(
            [python_exe, "-c", "import numpy, h5py"],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return False
    return probe.returncode == 0


def resolve_prepare_python(args):
    if args.prepare_python:
        exe = str(Path(args.prepare_python).expanduser())
        if not python_has_prepare_deps(exe):
            raise RuntimeError(
                f"--prepare-python does not have required deps (numpy,h5py): {exe}"
            )
        return exe

    candidates = [sys.executable, "/opt/homebrew/bin/python3.10", "python3.10", "python3"]
    for exe in candidates:
        if python_has_prepare_deps(exe):
            return exe
    raise RuntimeError(
        "Unable to find a Python executable with numpy+h5py for prepare_system.py. "
        "Pass --prepare-python explicitly."
    )


def has_atoms_section(itp_path):
    with itp_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].strip().lower()
            if line in {"[ atoms ]", "[atoms]"}:
                return True
    return False


def resolve_itp_from_top(top_file):
    top_file = Path(top_file).expanduser().resolve()
    base = top_file.parent
    include_candidates = []
    with top_file.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            m = re.search(r'#include\s+"([^"]+\.itp)"', raw)
            if not m:
                continue
            p = (base / m.group(1)).resolve()
            if p.exists():
                include_candidates.append(p)

    for p in include_candidates:
        if has_atoms_section(p):
            return p

    for p in sorted(base.glob("*.itp")):
        if has_atoms_section(p):
            return p.resolve()
    return None


def top_was_generated_from(top_file, input_pdb):
    top_file = Path(top_file).expanduser().resolve()
    if not top_file.exists():
        return False
    marker = str(Path(input_pdb).expanduser().resolve())
    with top_file.open("r", encoding="utf-8", errors="ignore") as fh:
        head = "".join([next(fh, "") for _ in range(40)])
    return marker in head


def run_martinize(args, python_exe, pdb_id, aa_pdb_effective, output_root):
    martinize_script = Path(args.martinize_script).expanduser().resolve()
    if not martinize_script.exists():
        raise FileNotFoundError(f"martinize.py not found: {martinize_script}")

    martinize_dir = output_root / "martinize" / pdb_id
    martinize_dir.mkdir(parents=True, exist_ok=True)
    cg_pdb = (martinize_dir / f"{pdb_id}.MARTINI.pdb").resolve()
    top_file = (martinize_dir / f"{pdb_id}.top").resolve()

    if (
        args.skip_existing
        and cg_pdb.exists()
        and top_file.exists()
        and top_was_generated_from(top_file, aa_pdb_effective)
    ):
        itp_file = resolve_itp_from_top(top_file)
        if itp_file is not None and itp_file.exists():
            return cg_pdb, itp_file.resolve(), "martinize_cached"

    cmd = [
        python_exe,
        str(martinize_script),
        "-f",
        str(aa_pdb_effective.resolve()),
        "-x",
        str(cg_pdb),
        "-o",
        str(top_file),
        "-ff",
        args.martinize_ff,
        "-name",
        args.martinize_molname,
    ]
    proc = subprocess.run(
        cmd,
        cwd=str(martinize_dir),
        check=False,
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"martinize.py failed for {pdb_id} (exit {proc.returncode})\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )
    if not cg_pdb.exists():
        raise FileNotFoundError(f"martinize.py did not produce CG PDB: {cg_pdb}")
    if not top_file.exists():
        raise FileNotFoundError(f"martinize.py did not produce topology file: {top_file}")

    itp_file = resolve_itp_from_top(top_file)
    if itp_file is None or not itp_file.exists():
        raise FileNotFoundError(f"Unable to resolve protein ITP from top file: {top_file}")
    return cg_pdb, itp_file.resolve(), "martinize"


def resolve_cg_and_itp(args, python_exe, pdb_id, aa_pdb_effective, output_root):
    if args.auto_martinize:
        return run_martinize(args, python_exe, pdb_id, aa_pdb_effective, output_root)

    templ_cg = render_template(args.protein_cg_template, pdb_id)
    templ_itp = render_template(args.protein_itp_template, pdb_id)
    if templ_cg.exists() and templ_itp.exists():
        return templ_cg, templ_itp, "template"

    msg = f"Missing CG/ITP for {pdb_id}: cg={templ_cg.exists()} itp={templ_itp.exists()}"
    if args.allow_missing_cg:
        return None, None, f"skip:{msg}"
    raise FileNotFoundError(msg)


def run_prepare_system(args, prepare_python, pdb_id, aa_pdb, cg_pdb, itp_file, output_root):
    prep_script = Path(args.prepare_script).expanduser().resolve()
    if not prep_script.exists():
        raise FileNotFoundError(f"prepare_system.py not found: {prep_script}")

    runtime_dir = output_root / "runtime"
    run_dir = output_root / "runs" / pdb_id
    mapping_dir = output_root / "mapping"
    summary_dir = output_root / "summary"

    runtime_pdb = runtime_dir / f"{pdb_id}.MARTINI.pdb"
    runtime_itp = runtime_dir / f"{pdb_id}_proa.itp"
    mapping_h5 = mapping_dir / f"{pdb_id}.hybrid_mapping.h5"
    mapping_json = mapping_dir / f"{pdb_id}.hybrid_bb_map.json"
    prep_summary = summary_dir / f"{pdb_id}.prepare_summary.json"

    runtime_dir.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)
    mapping_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        prepare_python,
        str(prep_script),
        "--mode",
        "both",
        "--pdb-id",
        pdb_id,
        "--runtime-pdb-output",
        str(runtime_pdb),
        "--runtime-itp-output",
        str(runtime_itp),
        "--prepare-structure",
        "1",
        "--stage",
        args.stage,
        "--run-dir",
        str(run_dir),
        "--bilayer-pdb",
        str(Path(args.bilayer_pdb).expanduser().resolve()),
        "--protein-aa-pdb",
        str(aa_pdb),
        "--protein-cg-pdb",
        str(cg_pdb),
        "--protein-itp",
        str(itp_file),
        "--hybrid-mapping-output",
        str(mapping_h5),
        "--hybrid-bb-map-json-output",
        str(mapping_json),
        "--xy-scale",
        str(args.xy_scale),
        "--box-padding-xy",
        str(args.box_padding_xy),
        "--box-padding-z",
        str(args.box_padding_z),
        "--salt-molar",
        str(args.salt_molar),
        "--protein-lipid-cutoff",
        str(args.protein_lipid_cutoff),
        "--ion-cutoff",
        str(args.ion_cutoff),
        "--seed",
        str(args.seed),
        "--summary-json",
        str(prep_summary),
    ]

    # Pass through current environment and set UPSIDE_HOME if absent.
    env = os.environ.copy()
    if not env.get("UPSIDE_HOME"):
        env["UPSIDE_HOME"] = str(script_dir_from_prepare(prep_script))

    proc = subprocess.run(cmd, check=False, text=True, capture_output=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(
            f"prepare_system.py failed for {pdb_id} (exit {proc.returncode})\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )

    up_file = run_dir / "test.input.up"
    if not up_file.exists():
        raise FileNotFoundError(f"Expected UPSIDE input missing after prep: {up_file}")

    return {
        "up_file": manifest_path_value(up_file),
        "runtime_pdb": manifest_path_value(runtime_pdb),
        "runtime_itp": manifest_path_value(runtime_itp),
        "mapping_h5": manifest_path_value(mapping_h5),
        "mapping_json": manifest_path_value(mapping_json),
        "prep_summary": manifest_path_value(prep_summary),
        "run_dir": manifest_path_value(run_dir),
    }


def script_dir_from_prepare(prepare_script):
    # /.../example/16.MARTINI/prepare_system.py -> repo root is ../../
    return prepare_script.resolve().parents[2]


def manifest_path_value(path):
    path = Path(path).expanduser().resolve()
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def main():
    args = parse_args()
    prepare_python = resolve_prepare_python(args)

    list_file = Path(args.list_file).expanduser().resolve()
    if not list_file.exists():
        raise FileNotFoundError(f"PDB list file not found: {list_file}")

    output_root = Path(args.output_root).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    manifest_out = (
        Path(args.manifest_out).expanduser().resolve()
        if args.manifest_out
        else (output_root / "training_manifest.json")
    )

    pdb_ids = read_pdb_list(list_file)
    if not pdb_ids:
        raise ValueError(f"No PDB IDs found in list file: {list_file}")

    systems = []
    failed = []

    for pdb_id in pdb_ids:
        aa_dir = output_root / "aa_pdb"
        aa_pdb = aa_dir / f"{pdb_id}.pdb"

        try:
            if args.skip_existing and aa_pdb.exists():
                source_url = "existing_local"
            else:
                source_url = download_aa_pdb(pdb_id, aa_pdb, timeout=args.timeout)

            aa_pdb_effective, aa_clean_meta = clean_aa_for_martinize(
                args=args,
                pdb_id=pdb_id,
                aa_pdb=aa_pdb,
                output_root=output_root,
            )

            cg_pdb, itp_file, cg_source = resolve_cg_and_itp(
                args=args,
                python_exe=prepare_python,
                pdb_id=pdb_id,
                aa_pdb_effective=aa_pdb_effective,
                output_root=output_root,
            )
            if cg_pdb is None or itp_file is None:
                msg = cg_source[5:] if cg_source.startswith("skip:") else "missing CG/ITP"
                failed.append({"pdb_id": pdb_id, "reason": msg})
                print(f"SKIP {pdb_id}: {msg}")
                continue

            prep = run_prepare_system(
                args=args,
                prepare_python=prepare_python,
                pdb_id=pdb_id,
                aa_pdb=aa_pdb_effective,
                cg_pdb=cg_pdb,
                itp_file=itp_file,
                output_root=output_root,
            )
            entry = {
                "pdb_id": pdb_id,
                "aa_pdb_downloaded": manifest_path_value(aa_pdb),
                "aa_pdb_effective": manifest_path_value(aa_pdb_effective),
                "cg_pdb": manifest_path_value(cg_pdb),
                "itp_file": manifest_path_value(itp_file),
                "cg_source": cg_source,
                "download_url": source_url,
            }
            entry["aa_clean"] = aa_clean_meta
            entry.update(prep)
            systems.append(entry)
            print(f"OK {pdb_id}: {entry['up_file']}")

        except Exception as exc:
            failed.append({"pdb_id": pdb_id, "reason": str(exc)})
            print(f"FAIL {pdb_id}: {exc}")

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "list_file": manifest_path_value(list_file),
        "systems": systems,
        "failed": failed,
    }

    manifest_out.parent.mkdir(parents=True, exist_ok=True)
    with manifest_out.open("w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(f"Manifest written: {manifest_out}")
    print(f"Prepared systems: {len(systems)}")
    print(f"Failed systems: {len(failed)}")


if __name__ == "__main__":
    main()
