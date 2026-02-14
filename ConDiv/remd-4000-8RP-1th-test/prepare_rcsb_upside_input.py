#!/usr/bin/env python3
import argparse
import csv
import difflib
import http.client
import json
import math
import pickle
import random
import socket
import subprocess
import sys
import time
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
MEMBRANE_TERMS = (
    "membrane",
    "transmembrane",
    "lipid bilayer",
    "micelle",
    "nanodisc",
)
STANDARD_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
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
}


def parse_args():
    script_dir = Path(__file__).resolve().parent
    default_converter = (script_dir / ".." / ".." / "py" / "PDB_to_initial_structure.py").resolve()

    p = argparse.ArgumentParser(
        description="Filter/download RCSB structures and prepare Upside input files."
    )
    p.add_argument(
        "--filter-profile",
        choices=("contrastive", "sidechain"),
        default="contrastive",
        help="contrastive: 50-100 residues, sidechain: 50-500 residues.",
    )
    p.add_argument("--output-root", default="rcsb_dataset", help="Output root directory")
    p.add_argument("--min-residues", type=int, default=None, help="Minimum sequence length")
    p.add_argument("--max-residues", type=int, default=None, help="Maximum sequence length")
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
        default=2.2,
        help="Max allowed X-ray resolution in angstrom",
    )
    p.add_argument(
        "--seqid-threshold",
        type=float,
        default=0.30,
        help="Maximum pairwise sequence similarity among accepted proteins.",
    )
    p.add_argument(
        "--globularity-ransac-iter",
        type=int,
        default=400,
        help="RANSAC iterations for globularity outlier filtering.",
    )
    p.add_argument(
        "--globularity-threshold",
        type=float,
        default=0.12,
        help="Residual threshold in log(Nres)-log(Rg) for globularity filtering.",
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
    p.add_argument(
        "--network-retries",
        type=int,
        default=3,
        help="Retry count for transient network failures (default: 3).",
    )
    p.add_argument(
        "--network-retry-backoff",
        type=float,
        default=1.0,
        help="Base backoff seconds between network retries (default: 1.0).",
    )
    p.add_argument(
        "--exclude-membrane",
        dest="exclude_membrane",
        action="store_true",
        help="Exclude membrane-related proteins based on RCSB metadata text fields (default: enabled).",
    )
    p.add_argument(
        "--include-membrane",
        dest="exclude_membrane",
        action="store_false",
        help="Disable membrane-related exclusion.",
    )
    p.set_defaults(exclude_membrane=True)
    return p.parse_args()


def _is_retryable_network_error(exc):
    if isinstance(exc, urllib.error.HTTPError):
        return exc.code in {408, 425, 429, 500, 502, 503, 504}
    if isinstance(exc, urllib.error.URLError):
        reason = getattr(exc, "reason", None)
        if isinstance(reason, socket.timeout):
            return True
        if isinstance(reason, ConnectionResetError):
            return True
        if isinstance(reason, OSError):
            return True
        return False
    return isinstance(
        exc,
        (http.client.IncompleteRead, http.client.RemoteDisconnected, TimeoutError, socket.timeout),
    )


def _read_url_bytes(req, timeout, retries, backoff):
    attempt = 0
    while True:
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:
            if attempt >= retries or not _is_retryable_network_error(exc):
                raise
            sleep_s = backoff * (2**attempt)
            time.sleep(sleep_s)
            attempt += 1


def _network_error_reason(exc):
    if isinstance(exc, urllib.error.HTTPError):
        return f"http_{exc.code}"
    if isinstance(exc, urllib.error.URLError):
        reason = getattr(exc, "reason", None)
        return str(reason if reason is not None else exc)
    if isinstance(exc, http.client.IncompleteRead):
        return f"incomplete_read:{len(exc.partial)}"
    return exc.__class__.__name__


def post_json(url, payload, timeout, retries, backoff):
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    content = _read_url_bytes(req, timeout=timeout, retries=retries, backoff=backoff)
    return json.loads(content.decode("utf-8"))


def get_json(url, timeout, retries, backoff):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    content = _read_url_bytes(req, timeout=timeout, retries=retries, backoff=backoff)
    return json.loads(content.decode("utf-8"))


def fetch_candidate_ids(methods, max_candidates, timeout, retries, backoff):
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

        data = post_json(
            RCSB_SEARCH_URL,
            payload,
            timeout=timeout,
            retries=retries,
            backoff=backoff,
        )
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


def is_membrane_related_entry(entry):
    fields = []
    struct = entry.get("struct", {})
    struct_keywords = entry.get("struct_keywords", {})
    exptl = entry.get("exptl", [])

    if isinstance(struct, dict):
        fields.append(str(struct.get("title", "")))
    if isinstance(struct_keywords, dict):
        fields.append(str(struct_keywords.get("pdbx_keywords", "")))
        text = struct_keywords.get("text")
        if isinstance(text, list):
            fields.extend(str(x) for x in text)
        else:
            fields.append(str(text or ""))
    if isinstance(exptl, list):
        for e in exptl:
            if isinstance(e, dict):
                fields.append(str(e.get("details", "")))

    blob = " ".join(fields).lower()
    return any(term in blob for term in MEMBRANE_TERMS)


def is_disallowed_identifier(pdb_id):
    up = pdb_id.upper()
    return up.startswith("AF_") or up.startswith("MA_")


def download_pdb_file(pdb_id, out_path, timeout, retries, backoff):
    url = RCSB_PDB_URL.format(pdb_id=pdb_id.upper())
    req = urllib.request.Request(url)
    content = _read_url_bytes(req, timeout=timeout, retries=retries, backoff=backoff)
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


def pdb_has_nonstandard_residue(pdb_file):
    with open(pdb_file, "r", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM  "):
                resname = line[17:20].strip().upper()
                if resname not in STANDARD_RESNAMES:
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


def read_fasta_sequence(fasta_path):
    seq = []
    with open(fasta_path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq)


def seq_similarity(seq_a, seq_b):
    return difflib.SequenceMatcher(a=seq_a, b=seq_b).ratio()


def compute_rg_from_initial(base):
    pkl_path = Path(str(base) + ".initial.pkl")
    npy_path = Path(str(base) + ".initial.npy")
    if pkl_path.exists():
        with open(pkl_path, "rb") as f:
            native_pos = pickle.load(f, encoding="latin1")[:, :, 0]
    elif npy_path.exists():
        native_pos = np.load(npy_path)
    else:
        raise FileNotFoundError(f"Missing initial structure for {base}")

    ca = native_pos[1::3]
    center = ca.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((ca - center) ** 2, axis=1))))


def globularity_keep_ids(records, n_iter=400, threshold=0.12):
    # records contain: pdb_id, n_res, rg
    if len(records) < 12:
        return set(r["pdb_id"] for r in records)

    x = np.array([math.log(float(r["n_res"])) for r in records], dtype=float)
    y = np.array([math.log(float(r["rg"])) for r in records], dtype=float)
    idx_all = list(range(len(records)))

    best_inliers = None
    best_count = -1
    for _ in range(max(1, n_iter)):
        i, j = random.sample(idx_all, 2)
        if abs(x[j] - x[i]) < 1e-12:
            continue
        a = (y[j] - y[i]) / (x[j] - x[i])
        b = y[i] - a * x[i]
        resid = np.abs(y - (a * x + b))
        inliers = resid <= threshold
        count = int(np.sum(inliers))
        if count > best_count:
            best_count = count
            best_inliers = inliers

    if best_inliers is None or int(np.sum(best_inliers)) < 3:
        return set(r["pdb_id"] for r in records)

    x_in = x[best_inliers]
    y_in = y[best_inliers]
    A = np.vstack([x_in, np.ones_like(x_in)]).T
    a_refit, b_refit = np.linalg.lstsq(A, y_in, rcond=None)[0]

    resid_final = np.abs(y - (a_refit * x + b_refit))
    keep = resid_final <= threshold
    return set(records[i]["pdb_id"] for i in range(len(records)) if keep[i])


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
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(rows)


def main():
    args = parse_args()
    if args.min_residues is None or args.max_residues is None:
        if args.filter_profile == "contrastive":
            args.min_residues = 50
            args.max_residues = 100
        else:
            args.min_residues = 50
            args.max_residues = 500

    converter = Path(args.converter).resolve()
    if not converter.exists():
        raise FileNotFoundError(f"Converter script not found: {converter}")

    dirs = ensure_dirs(args.output_root)
    manifest_path = Path(args.output_root) / "manifest.csv"

    try:
        xray_ids = fetch_candidate_ids(
            XRAY_METHODS,
            args.max_candidates_per_class,
            args.timeout,
            args.network_retries,
            args.network_retry_backoff,
        )
        nmr_ids = fetch_candidate_ids(
            NMR_METHODS,
            args.max_candidates_per_class,
            args.timeout,
            args.network_retries,
            args.network_retry_backoff,
        )
    except (urllib.error.URLError, urllib.error.HTTPError, http.client.IncompleteRead) as e:
        raise RuntimeError(
            f"Failed to query RCSB search API ({RCSB_SEARCH_URL}). "
            "Check network/DNS access from this machine."
        ) from e

    accepted = {"xray": [], "nmr": []}
    accepted_rows = {"xray": [], "nmr": []}
    accepted_sequences = []
    rows = []
    kept_row_by_id = {}

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
                entry = get_json(
                    RCSB_ENTRY_URL.format(pdb_id=pdb_id.upper()),
                    timeout=args.timeout,
                    retries=args.network_retries,
                    backoff=args.network_retry_backoff,
                )
            except urllib.error.HTTPError as e:
                row["reason"] = f"entry_http_error:{e.code}"
                rows.append(row)
                continue
            except (urllib.error.URLError, http.client.IncompleteRead) as e:
                row["reason"] = f"entry_fetch_failed:{_network_error_reason(e)}"
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
            if args.exclude_membrane and is_membrane_related_entry(entry):
                row["reason"] = "membrane_related_entry"
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
                download_pdb_file(
                    pdb_id,
                    pdb_file,
                    timeout=args.timeout,
                    retries=args.network_retries,
                    backoff=args.network_retry_backoff,
                )
            except urllib.error.HTTPError as e:
                row["reason"] = f"download_http_error:{e.code}"
                rows.append(row)
                continue
            except (urllib.error.URLError, http.client.IncompleteRead) as e:
                row["reason"] = f"download_failed:{_network_error_reason(e)}"
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

            if pdb_has_nonstandard_residue(pdb_file):
                row["reason"] = "contains_nonstandard_residue"
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

            seq = read_fasta_sequence(Path(str(out_base) + ".fasta"))
            too_similar = False
            for _, ref_seq in accepted_sequences:
                if seq_similarity(seq, ref_seq) > args.seqid_threshold:
                    too_similar = True
                    break
            if too_similar:
                cleanup_generated(out_base)
                row["reason"] = "seq_similarity_too_high"
                rows.append(row)
                continue

            rg = compute_rg_from_initial(out_base)

            row["status"] = "keep"
            row["reason"] = ""
            rows.append(row)
            kept_row_by_id[pdb_id] = row
            accepted[class_name].append(pdb_id)
            accepted_sequences.append((pdb_id, seq))
            accepted_rows[class_name].append(
                {
                    "pdb_id": pdb_id,
                    "class": class_name,
                    "deposit_date": row["deposit_date"],
                    "deposit_year": row["deposit_year"],
                    "method": row["method"],
                    "resolution": row["resolution"],
                    "length": str(n_res),
                    "n_res": int(n_res),
                    "rg": float(rg),
                    "out_base": str(out_base),
                }
            )

    combined_records = accepted_rows["xray"] + accepted_rows["nmr"]
    keep_ids = globularity_keep_ids(
        combined_records,
        n_iter=args.globularity_ransac_iter,
        threshold=args.globularity_threshold,
    )
    for class_name in ("xray", "nmr"):
        kept_ids = []
        kept_meta = []
        for r in accepted_rows[class_name]:
            if r["pdb_id"] in keep_ids:
                kept_ids.append(r["pdb_id"])
                kept_meta.append(r)
            else:
                cleanup_generated(Path(r["out_base"]))
                row = kept_row_by_id.get(r["pdb_id"])
                if row is not None:
                    row["status"] = "skip"
                    row["reason"] = "globularity_outlier"
        accepted[class_name] = kept_ids
        accepted_rows[class_name] = kept_meta

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
