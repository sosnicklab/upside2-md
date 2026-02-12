#!/usr/bin/env python3
import argparse
import csv
import json
import pickle
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

XRAY_METHODS = ("X-RAY DIFFRACTION",)
NMR_METHODS = ("SOLUTION NMR", "SOLID-STATE NMR")
HEME_RESNAMES = {
    "HEM",
    "HEA",
    "HEB",
    "HEC",
    "HED",
    "HEE",
    "HEG",
    "HEO",
    "HEV",
}


def parse_args():
    script_dir = Path(__file__).resolve().parent
    default_converter = (script_dir / ".." / ".." / "py" / "PDB_to_initial_structure.py").resolve()

    p = argparse.ArgumentParser(
        description="Filter/download RCSB structures and prepare Upside input files."
    )
    p.add_argument("--output-root", default="rcsb_dataset", help="Output root directory")
    p.add_argument("--min-residues", type=int, default=50, help="Minimum sequence length")
    p.add_argument("--max-residues", type=int, default=151, help="Maximum sequence length")
    p.add_argument(
        "--max-candidates-per-class",
        type=int,
        default=1000,
        help="Max RCSB candidates to inspect for each class (xray/nmr)",
    )
    p.add_argument(
        "--max-keep-per-class",
        type=int,
        default=250,
        help="Max accepted structures to keep for each class",
    )
    p.add_argument(
        "--xray-resolution-max",
        type=float,
        default=2.5,
        help="Max allowed X-ray resolution in angstrom",
    )
    p.add_argument(
        "--nmr-model",
        default="1",
        help="Model index to extract from NMR PDB files (passed to PDB_to_initial_structure.py --model)",
    )
    p.add_argument(
        "--converter",
        default=str(default_converter),
        help="Path to PDB_to_initial_structure.py",
    )
    p.add_argument("--timeout", type=int, default=30, help="Network timeout in seconds")
    return p.parse_args()


def post_json(url, payload, timeout):
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.load(resp)


def get_json(url, timeout):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.load(resp)


def fetch_candidate_ids(methods, max_candidates, timeout):
    ids = []
    page_size = min(1000, max_candidates)
    start = 0

    while len(ids) < max_candidates:
        payload = {
            "query": {
                "type": "group",
                "logical_operator": "or",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "exptl.method",
                            "operator": "exact_match",
                            "value": method,
                        },
                    }
                    for method in methods
                ],
            },
            "return_type": "entry",
            "request_options": {"paginate": {"start": start, "rows": page_size}},
        }

        data = post_json(RCSB_SEARCH_URL, payload, timeout=timeout)
        result_set = data.get("result_set", [])
        if not result_set:
            break

        for item in result_set:
            ident = item.get("identifier", "").lower()
            if ident:
                ids.append(ident)
                if len(ids) >= max_candidates:
                    break

        start += page_size

    # Preserve order, remove duplicates
    seen = set()
    uniq = []
    for pdb_id in ids:
        if pdb_id not in seen:
            seen.add(pdb_id)
            uniq.append(pdb_id)
    return uniq


def extract_entry_filters(entry):
    info = entry.get("rcsb_entry_info", {})
    accession = entry.get("rcsb_accession_info", {})
    methods = [x.get("method", "") for x in entry.get("exptl", []) if isinstance(x, dict)]
    method_set = set(methods)
    res_list = info.get("resolution_combined", [])
    resolution = None
    if isinstance(res_list, list) and res_list:
        try:
            resolution = float(res_list[0])
        except (TypeError, ValueError):
            resolution = None

    length = info.get("deposited_polymer_monomer_count")
    try:
        length = int(length)
    except (TypeError, ValueError):
        length = None

    poly_types = info.get("selected_polymer_entity_types", "")
    if isinstance(poly_types, list):
        poly_types = " ".join(str(x) for x in poly_types)
    poly_types = str(poly_types)

    deposit_date = accession.get("deposit_date")
    if not deposit_date:
        deposit_date = accession.get("initial_release_date", "")
    deposit_date = str(deposit_date) if deposit_date is not None else ""

    deposit_year = ""
    if len(deposit_date) >= 4 and deposit_date[:4].isdigit():
        deposit_year = deposit_date[:4]

    return method_set, resolution, length, poly_types, deposit_date, deposit_year


def is_disallowed_identifier(pdb_id):
    up = pdb_id.upper()
    return up.startswith("AF_") or up.startswith("MA_")


def download_pdb_file(pdb_id, out_path, timeout):
    url = RCSB_PDB_URL.format(pdb_id=pdb_id.upper())
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        content = resp.read()
    out_path.write_bytes(content)


def protein_chain_ids_from_pdb(pdb_file):
    chains = set()
    with open(pdb_file, "r", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM  "):
                chain_id = line[21].strip()
                chains.add(chain_id if chain_id else "_")
    return chains


def pdb_has_heme(pdb_file):
    with open(pdb_file, "r", errors="replace") as f:
        for line in f:
            if line.startswith("HETATM"):
                resname = line[17:20].strip().upper()
                if resname in HEME_RESNAMES:
                    return True
    return False


def run_converter(converter, pdb_file, out_base, nmr_model=None):
    cmd = [sys.executable, str(converter), str(pdb_file), str(out_base)]
    cmd.append("--record-chain-breaks")
    if nmr_model is not None:
        cmd.extend(["--model", str(nmr_model)])
    subprocess.run(cmd, check=True)


def convert_initial_npy_to_pkl(base):
    npy_path = Path(str(base) + ".initial.npy")
    pkl_path = Path(str(base) + ".initial.pkl")
    if not npy_path.exists():
        return False

    coords = np.load(npy_path)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"Unexpected .initial.npy shape: {coords.shape}")

    with open(pkl_path, "wb") as f:
        pickle.dump(coords[:, :, None], f, -1)
    return True


def chain_break_ok(base):
    pkl_path = Path(str(base) + ".initial.pkl")
    npy_path = Path(str(base) + ".initial.npy")
    if pkl_path.exists():
        with open(pkl_path, "rb") as f:
            native_pos = pickle.load(f, encoding="latin1")[:, :, 0]
    elif npy_path.exists():
        native_pos = np.load(npy_path)
    else:
        raise FileNotFoundError(f"Missing initial structure for {base}")

    max_sep = np.sqrt(np.sum(np.diff(native_pos, axis=0) ** 2, axis=-1)).max()
    return max_sep < 2.0


def derived_n_res(base):
    pkl_path = Path(str(base) + ".initial.pkl")
    npy_path = Path(str(base) + ".initial.npy")
    if pkl_path.exists():
        with open(pkl_path, "rb") as f:
            native_pos = pickle.load(f, encoding="latin1")[:, :, 0]
    elif npy_path.exists():
        native_pos = np.load(npy_path)
    else:
        raise FileNotFoundError(f"Missing initial structure for {base}")

    n_atom = len(native_pos)
    if n_atom % 3 != 0:
        return None
    return n_atom // 3


def cleanup_generated(base):
    for ext in (".fasta", ".initial.pkl", ".initial.npy", ".chi", ".chain_breaks"):
        p = Path(str(base) + ext)
        if p.exists():
            p.unlink()


def ensure_dirs(root):
    root = Path(root)
    dirs = {
        "xray_raw": root / "xray" / "raw_pdb",
        "xray_upside": root / "xray" / "upside_input",
        "nmr_raw": root / "nmr" / "raw_pdb",
        "nmr_upside": root / "nmr" / "upside_input",
        "lists": root / "lists",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def write_pdb_list(path, ids):
    with open(path, "w") as f:
        f.write("prot\n")
        for x in sorted(ids):
            f.write(f"{x}\n")


def write_index_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["pdb_id", "class", "deposit_date", "deposit_year", "method", "resolution", "length"],
        )
        w.writeheader()
        w.writerows(rows)


def main():
    args = parse_args()
    converter = Path(args.converter).resolve()
    if not converter.exists():
        raise FileNotFoundError(f"Converter script not found: {converter}")

    dirs = ensure_dirs(args.output_root)
    manifest_path = Path(args.output_root) / "manifest.csv"

    try:
        xray_ids = fetch_candidate_ids(XRAY_METHODS, args.max_candidates_per_class, args.timeout)
        nmr_ids = fetch_candidate_ids(NMR_METHODS, args.max_candidates_per_class, args.timeout)
    except urllib.error.URLError as e:
        raise RuntimeError(
            f"Failed to query RCSB search API ({RCSB_SEARCH_URL}). "
            "Check network/DNS access from this machine."
        ) from e

    accepted = {"xray": [], "nmr": []}
    accepted_rows = {"xray": [], "nmr": []}
    rows = []

    for class_name, ids in (("xray", xray_ids), ("nmr", nmr_ids)):
        for pdb_id in ids:
            if len(accepted[class_name]) >= args.max_keep_per_class:
                break

            row = {
                "pdb_id": pdb_id,
                "class": class_name,
                "status": "skip",
                "reason": "",
                "method": "",
                "resolution": "",
                "length": "",
                "deposit_date": "",
                "deposit_year": "",
            }

            if is_disallowed_identifier(pdb_id):
                row["reason"] = "disallowed_identifier_prefix"
                rows.append(row)
                continue

            try:
                entry = get_json(RCSB_ENTRY_URL.format(pdb_id=pdb_id.upper()), timeout=args.timeout)
            except urllib.error.URLError as e:
                row["reason"] = f"entry_fetch_failed:{e.reason}"
                rows.append(row)
                continue
            except urllib.error.HTTPError as e:
                row["reason"] = f"entry_http_error:{e.code}"
                rows.append(row)
                continue

            methods, resolution, length, poly_types, deposit_date, deposit_year = extract_entry_filters(entry)
            row["method"] = ";".join(sorted(methods))
            row["resolution"] = "" if resolution is None else f"{resolution:.3f}"
            row["length"] = "" if length is None else str(length)
            row["deposit_date"] = deposit_date
            row["deposit_year"] = deposit_year

            if "protein" not in poly_types.lower():
                row["reason"] = "not_protein_entry"
                rows.append(row)
                continue
            if "dna" in poly_types.lower() or "rna" in poly_types.lower():
                row["reason"] = "contains_nucleic_acid"
                rows.append(row)
                continue
            if length is None or not (args.min_residues <= length <= args.max_residues):
                row["reason"] = "length_out_of_range"
                rows.append(row)
                continue
            if class_name == "xray":
                if not methods.intersection(XRAY_METHODS):
                    row["reason"] = "not_xray_method"
                    rows.append(row)
                    continue
                if resolution is None or resolution > args.xray_resolution_max:
                    row["reason"] = "xray_resolution_too_low"
                    rows.append(row)
                    continue
            if class_name == "nmr":
                if not methods.intersection(NMR_METHODS):
                    row["reason"] = "not_nmr_method"
                    rows.append(row)
                    continue

            raw_dir = dirs["xray_raw"] if class_name == "xray" else dirs["nmr_raw"]
            upside_dir = dirs["xray_upside"] if class_name == "xray" else dirs["nmr_upside"]
            pdb_file = raw_dir / f"{pdb_id}.pdb"
            out_base = upside_dir / pdb_id

            try:
                download_pdb_file(pdb_id, pdb_file, timeout=args.timeout)
            except urllib.error.HTTPError as e:
                row["reason"] = f"download_http_error:{e.code}"
                rows.append(row)
                continue
            except urllib.error.URLError as e:
                row["reason"] = f"download_failed:{e.reason}"
                rows.append(row)
                continue

            chain_ids = protein_chain_ids_from_pdb(pdb_file)
            if len(chain_ids) != 1:
                row["reason"] = "not_single_chain"
                rows.append(row)
                continue

            if pdb_has_heme(pdb_file):
                row["reason"] = "contains_heme"
                rows.append(row)
                continue

            try:
                run_converter(
                    converter,
                    pdb_file,
                    out_base,
                    nmr_model=args.nmr_model if class_name == "nmr" else None,
                )
            except subprocess.CalledProcessError:
                cleanup_generated(out_base)
                row["reason"] = "converter_failed"
                rows.append(row)
                continue

            try:
                convert_initial_npy_to_pkl(out_base)
                ok = chain_break_ok(out_base)
                n_res = derived_n_res(out_base)
            except Exception:
                cleanup_generated(out_base)
                row["reason"] = "initial_parse_failed"
                rows.append(row)
                continue

            if n_res is None:
                cleanup_generated(out_base)
                row["reason"] = "invalid_backbone_atom_count"
                rows.append(row)
                continue

            if not (args.min_residues <= n_res <= args.max_residues):
                cleanup_generated(out_base)
                row["reason"] = "derived_length_out_of_range"
                rows.append(row)
                continue

            if not ok:
                cleanup_generated(out_base)
                row["reason"] = "chain_break_filter_failed"
                rows.append(row)
                continue

            row["status"] = "keep"
            row["reason"] = ""
            rows.append(row)
            accepted[class_name].append(pdb_id)
            accepted_rows[class_name].append(
                {
                    "pdb_id": pdb_id,
                    "class": class_name,
                    "deposit_date": row["deposit_date"],
                    "deposit_year": row["deposit_year"],
                    "method": row["method"],
                    "resolution": row["resolution"],
                    "length": str(n_res),
                }
            )

    with open(manifest_path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "pdb_id",
                "class",
                "status",
                "reason",
                "method",
                "resolution",
                "length",
                "deposit_date",
                "deposit_year",
            ],
        )
        w.writeheader()
        w.writerows(rows)

    write_pdb_list(dirs["lists"] / "pdb_list_xray", accepted["xray"])
    write_pdb_list(dirs["lists"] / "pdb_list_nmr", accepted["nmr"])
    write_pdb_list(dirs["lists"] / "pdb_list_combined", accepted["xray"] + accepted["nmr"])
    write_index_csv(dirs["lists"] / "structure_index_xray.csv", accepted_rows["xray"])
    write_index_csv(dirs["lists"] / "structure_index_nmr.csv", accepted_rows["nmr"])
    write_index_csv(
        dirs["lists"] / "structure_index_combined.csv",
        accepted_rows["xray"] + accepted_rows["nmr"],
    )

    print(f"Accepted X-ray: {len(accepted['xray'])}")
    print(f"Accepted NMR:   {len(accepted['nmr'])}")
    print(f"Manifest:       {manifest_path}")
    print(f"Lists:          {dirs['lists']}")


if __name__ == "__main__":
    main()
