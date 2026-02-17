#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../source.sh"
source "${SCRIPT_DIR}/../../.venv/bin/activate"
set -euo pipefail
cd "${SCRIPT_DIR}"

# Hybrid 1RKL workflow:
# 0) Hybrid preparation (packed MARTINI + hybrid mapping export)
# 1) Stage input generation (dry MARTINI)
# 2) Inject hybrid mapping into each stage .up file
# 3) Run 6.0 -> 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5 -> 6.6 -> 7.0
# 4) Extract key VTF trajectories

# =============================================================================
# USER CONFIGURATION
# =============================================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        PDB_ID=*)
            PDB_ID="${1#*=}"
            shift
            ;;
        *)
            echo "Unknown parameter $1"
            exit 1
            ;;
    esac
done

PDB_ID="${PDB_ID:-1rkl}"
RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-${PDB_ID}_hybrid}"

INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="${RUN_DIR:-outputs/martini_test_1rkl_hybrid}"
CHECKPOINT_DIR="${CHECKPOINT_DIR:-${RUN_DIR}/checkpoints}"
LOG_DIR="${LOG_DIR:-${RUN_DIR}/logs}"
HYBRID_PREP_DIR="${HYBRID_PREP_DIR:-${RUN_DIR}/hybrid_prep}"

PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/${PDB_ID}.pdb}"
PROTEIN_CG_PDB="${PROTEIN_CG_PDB:-pdb/${PDB_ID}.MARTINI.pdb}"
PROTEIN_ITP="${PROTEIN_ITP:-pdb/${PDB_ID}_proa.itp}"
BILAYER_PDB="${BILAYER_PDB:-pdb/bilayer.MARTINI.pdb}"

RUNTIME_PDB_FILE="${SCRIPT_DIR}/pdb/${RUNTIME_PDB_ID}.MARTINI.pdb"
RUNTIME_ITP_FILE="${SCRIPT_DIR}/pdb/${RUNTIME_PDB_ID}_proa.itp"
HYBRID_MAPPING_FILE="${HYBRID_MAPPING_FILE:-${HYBRID_PREP_DIR}/hybrid_mapping.h5}"
HYBRID_PACKED_PDB="${HYBRID_PACKED_PDB:-${HYBRID_PREP_DIR}/hybrid_packed.MARTINI.pdb}"
HYBRID_VALIDATE="${HYBRID_VALIDATE:-1}"

# martinize controls (MARTINI 2.x generation for dry-MARTINI compatibility)
MARTINIZE_ENABLE="${MARTINIZE_ENABLE:-1}"
MARTINIZE_FF="${MARTINIZE_FF:-martini22}"
MARTINIZE_MOLNAME="${MARTINIZE_MOLNAME:-PROA}"
MARTINIZE_SCRIPT="${MARTINIZE_SCRIPT:-${SCRIPT_DIR}/martinize.py}"

PREPARED_60_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.0.prepared.up"
STAGE_60_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.0.up"
PREPARED_61_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.1.prepared.up"
STAGE_61_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.1.up"
STAGE_62_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.2.up"
PREPARED_63_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.3.prepared.up"
STAGE_63_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.3.up"
PREPARED_64_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.4.prepared.up"
STAGE_64_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.4.up"
PREPARED_65_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.5.prepared.up"
STAGE_65_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.5.up"
PREPARED_66_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.6.prepared.up"
STAGE_66_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.6.up"
PREPARED_70_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.prepared.up"
STAGE_70_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.up"

# prepare_hybrid_system.py controls
SALT_MOLAR="${SALT_MOLAR:-0.15}"
PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-3.0}"
ION_CUTOFF="${ION_CUTOFF:-4.0}"
BOX_PADDING_XY="${BOX_PADDING_XY:-0.0}"
BOX_PADDING_Z="${BOX_PADDING_Z:-20.0}"
PREP_SEED="${PREP_SEED:-2026}"

# Simulation controls
TEMPERATURE="${TEMPERATURE:-0.8647}"
THERMOSTAT_TIMESCALE="${THERMOSTAT_TIMESCALE:-4.0}"
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
SEED="${SEED:-7090685331}"

MIN_60_MAX_ITER="${MIN_60_MAX_ITER:-500}"
MIN_61_MAX_ITER="${MIN_61_MAX_ITER:-500}"
EQ_62_NSTEPS="${EQ_62_NSTEPS:-500}"
EQ_63_NSTEPS="${EQ_63_NSTEPS:-500}"
EQ_64_NSTEPS="${EQ_64_NSTEPS:-500}"
EQ_65_NSTEPS="${EQ_65_NSTEPS:-500}"
EQ_66_NSTEPS="${EQ_66_NSTEPS:-500}"
PROD_70_NSTEPS="${PROD_70_NSTEPS:-5000}"

EQ_TIME_STEP="${EQ_TIME_STEP:-0.010}"
PROD_TIME_STEP="${PROD_TIME_STEP:-0.020}"
MIN_TIME_STEP="${MIN_TIME_STEP:-0.010}"

EQ_FRAME_STEPS="${EQ_FRAME_STEPS:-1000}"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-5000}"

COMP_3E4_BAR_INV_TO_A3_PER_EUP="${COMP_3E4_BAR_INV_TO_A3_PER_EUP:-14.521180763676}"
export UPSIDE_NPT_TARGET_PXY="${UPSIDE_NPT_TARGET_PXY:-0.0}"
export UPSIDE_NPT_TARGET_PZ="${UPSIDE_NPT_TARGET_PZ:-0.0}"
export UPSIDE_NPT_TAU="${UPSIDE_NPT_TAU:-4.0}"
export UPSIDE_NPT_COMPRESSIBILITY="${UPSIDE_NPT_COMPRESSIBILITY:-$COMP_3E4_BAR_INV_TO_A3_PER_EUP}"
export UPSIDE_NPT_COMPRESSIBILITY_XY="${UPSIDE_NPT_COMPRESSIBILITY_XY:-$COMP_3E4_BAR_INV_TO_A3_PER_EUP}"
export UPSIDE_NPT_COMPRESSIBILITY_Z="${UPSIDE_NPT_COMPRESSIBILITY_Z:-0.0}"
export UPSIDE_NPT_INTERVAL="${UPSIDE_NPT_INTERVAL:-10}"
export UPSIDE_NPT_SEMI="${UPSIDE_NPT_SEMI:-1}"
export UPSIDE_NPT_DEBUG="${UPSIDE_NPT_DEBUG:-1}"
PROD_70_NPT_ENABLE="${PROD_70_NPT_ENABLE:-0}"
PROD_70_BAROSTAT_TYPE="${PROD_70_BAROSTAT_TYPE:-1}"

export UPSIDE_OVERWRITE_SPLINES="${UPSIDE_OVERWRITE_SPLINES:-1}"
export UPSIDE_EWALD_ENABLE="${UPSIDE_EWALD_ENABLE:-1}"
export UPSIDE_EWALD_ALPHA="${UPSIDE_EWALD_ALPHA:-0.2}"
export UPSIDE_EWALD_KMAX="${UPSIDE_EWALD_KMAX:-5}"
export UPSIDE_MARTINI_FF_DIR="${UPSIDE_MARTINI_FF_DIR:-ff_dry}"
MASS_FF_FILE="${MASS_FF_FILE:-${UPSIDE_MARTINI_FF_DIR}/dry_martini_v2.1.itp}"
SIDECHAIN_LIBRARY="${SIDECHAIN_LIBRARY:-${UPSIDE_HOME}/parameters/ff_2.1/sidechain.h5}"

# =============================================================================
# VALIDATION
# =============================================================================
if [ -z "${UPSIDE_HOME:-}" ]; then
    echo "ERROR: UPSIDE_HOME environment variable is not set"
    exit 1
fi

UPSIDE_EXECUTABLE="${UPSIDE_HOME}/obj/upside"
if [ ! -f "$UPSIDE_EXECUTABLE" ]; then
    echo "ERROR: UPSIDE executable not found: $UPSIDE_EXECUTABLE"
    exit 1
fi

if [ ! -f "${PROTEIN_AA_PDB}" ]; then
    echo "ERROR: protein AA PDB not found: ${PROTEIN_AA_PDB}"
    exit 1
fi
if [ ! -f "${BILAYER_PDB}" ]; then
    echo "ERROR: bilayer PDB not found: ${BILAYER_PDB}"
    exit 1
fi
if [ "${MARTINIZE_ENABLE}" = "1" ]; then
    if [ ! -f "${MARTINIZE_SCRIPT}" ]; then
        echo "ERROR: martinize.py script not found: ${MARTINIZE_SCRIPT}"
        exit 1
    fi
fi
if [ "${MARTINIZE_ENABLE}" = "0" ]; then
    if [ ! -f "${PROTEIN_CG_PDB}" ]; then
        echo "ERROR: protein CG PDB not found: ${PROTEIN_CG_PDB}"
        exit 1
    fi
    if [ ! -f "${PROTEIN_ITP}" ]; then
        echo "ERROR: protein ITP not found: ${PROTEIN_ITP}"
        exit 1
    fi
fi

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR" "$HYBRID_PREP_DIR"

PROTEIN_CG_EFFECTIVE="${PROTEIN_CG_PDB}"
PROTEIN_ITP_EFFECTIVE="${PROTEIN_ITP}"

run_martinize() {
    local input_aa="$1"
    local out_cg_pdb="$2"
    local out_top="$3"

    local script_abs input_abs cg_abs top_abs workdir
    script_abs="$(python3 - "$MARTINIZE_SCRIPT" << 'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)"
    input_abs="$(python3 - "$input_aa" << 'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)"
    cg_abs="$(python3 - "$out_cg_pdb" << 'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)"
    top_abs="$(python3 - "$out_top" << 'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)"
    workdir="$(dirname "$top_abs")"
    mkdir -p "$workdir"

    echo "Running martinize.py command:"
    echo "  python3 ${script_abs} -f ${input_abs} -x ${cg_abs} -o ${top_abs} -ff ${MARTINIZE_FF} -name ${MARTINIZE_MOLNAME}"
    if ! (cd "$workdir" && python3 "$script_abs" -f "$input_abs" -x "$cg_abs" -o "$top_abs" -ff "$MARTINIZE_FF" -name "$MARTINIZE_MOLNAME"); then
        echo "ERROR: martinize.py execution failed."
        exit 1
    fi
}

resolve_itp_from_top() {
    local top_file="$1"
    local out_path
    out_path="$(python3 - "$top_file" << 'PY'
import os, re, sys, glob
top = os.path.abspath(sys.argv[1])
base = os.path.dirname(top)
cands = []
with open(top, "r", encoding="utf-8", errors="ignore") as fh:
    for ln in fh:
        m = re.search(r'#include\s+"([^"]+\.itp)"', ln)
        if not m:
            continue
        p = os.path.join(base, m.group(1))
        if os.path.exists(p):
            cands.append(p)

def has_atoms(path):
    txt = open(path, "r", encoding="utf-8", errors="ignore").read().lower()
    return ("[ atoms ]" in txt) or ("[atoms]" in txt)

for p in cands:
    if has_atoms(p):
        print(p)
        sys.exit(0)

for p in sorted(glob.glob(os.path.join(base, "*.itp"))):
    if has_atoms(p):
        print(p)
        sys.exit(0)

print("")
PY
)"
    echo "${out_path}"
}

assert_itp_types_have_masses() {
    local itp_file="$1"
    local ff_file="$2"

    python3 - "$itp_file" "$ff_file" << 'PY'
import os, sys

itp_file, ff_file = sys.argv[1], sys.argv[2]
if not os.path.exists(itp_file):
    raise SystemExit(f"ERROR: protein ITP not found: {itp_file}")
if not os.path.exists(ff_file):
    raise SystemExit(f"ERROR: mass force-field file not found: {ff_file}")

masses = set()
in_atomtypes = False
with open(ff_file, "r", encoding="utf-8", errors="ignore") as fh:
    for raw in fh:
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        low = line.lower()
        if low == "[ atomtypes ]" or low == "[atomtypes]":
            in_atomtypes = True
            continue
        if line.startswith("[") and line.endswith("]"):
            in_atomtypes = False
            continue
        if in_atomtypes:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    float(parts[1])
                except ValueError:
                    continue
                masses.add(parts[0])

types = set()
in_atoms = False
with open(itp_file, "r", encoding="utf-8", errors="ignore") as fh:
    for raw in fh:
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        low = line.lower()
        if low == "[ atoms ]" or low == "[atoms]":
            in_atoms = True
            continue
        if line.startswith("[") and line.endswith("]"):
            in_atoms = False
            continue
        if in_atoms:
            parts = line.split()
            if len(parts) >= 2:
                types.add(parts[1])

missing = sorted(t for t in types if t not in masses)
if missing:
    raise SystemExit(
        "ERROR: protein ITP contains bead types missing in dry-MARTINI mass table.\n"
        f"  ITP: {itp_file}\n"
        f"  FF:  {ff_file}\n"
        f"  Missing types: {missing}\n"
        "Use MARTINI2-compatible martinize.py settings (e.g., MARTINIZE_FF=martini22)."
    )

print(f"Protein ITP mass compatibility OK: {itp_file}")
PY
}

prepare_protein_inputs() {
    local mass_ff_path="${SCRIPT_DIR}/${MASS_FF_FILE}"
    if [ "${MARTINIZE_ENABLE}" = "1" ]; then
        local martinize_dir="${RUN_DIR}/martinize"
        mkdir -p "${martinize_dir}"
        local cg_pdb="${martinize_dir}/${RUNTIME_PDB_ID}.MARTINI.pdb"
        local top_file="${martinize_dir}/${RUNTIME_PDB_ID}.top"

        run_martinize "${PROTEIN_AA_PDB}" "${cg_pdb}" "${top_file}"
        if [ ! -f "${cg_pdb}" ]; then
            echo "ERROR: martinize.py did not produce CG PDB: ${cg_pdb}"
            exit 1
        fi
        if [ ! -f "${top_file}" ]; then
            echo "ERROR: martinize.py did not produce topology file: ${top_file}"
            exit 1
        fi

        local generated_itp
        generated_itp="$(resolve_itp_from_top "${top_file}")"
        if [ -z "${generated_itp}" ] || [ ! -f "${generated_itp}" ]; then
            echo "ERROR: unable to resolve protein ITP from martinize.py output top: ${top_file}"
            exit 1
        fi

        assert_itp_types_have_masses "${generated_itp}" "${mass_ff_path}"

        PROTEIN_CG_EFFECTIVE="${cg_pdb}"
        PROTEIN_ITP_EFFECTIVE="${generated_itp}"
    else
        if [ ! -f "${PROTEIN_CG_PDB}" ]; then
            echo "ERROR: protein CG PDB not found: ${PROTEIN_CG_PDB}"
            exit 1
        fi
        if [ ! -f "${PROTEIN_ITP}" ]; then
            echo "ERROR: protein ITP not found: ${PROTEIN_ITP}"
            exit 1
        fi
        assert_itp_types_have_masses "${PROTEIN_ITP}" "${mass_ff_path}"
        PROTEIN_CG_EFFECTIVE="${PROTEIN_CG_PDB}"
        PROTEIN_ITP_EFFECTIVE="${PROTEIN_ITP}"
    fi
}

set_barostat_type() {
    local up_file="$1"
    local barostat_type="$2"
    python3 - "$up_file" "$barostat_type" << 'PY'
import sys
import tables as tb
up_file = sys.argv[1]
barostat_type = int(sys.argv[2])
with tb.open_file(up_file, 'r+') as t:
    if '/input/barostat' in t:
        t.root.input.barostat._v_attrs.type = barostat_type
PY
}

set_stage_label() {
    local up_file="$1"
    local stage_label="$2"
    python3 - "$up_file" "$stage_label" << 'PY'
import sys
import h5py
import numpy as np
up_file = sys.argv[1]
stage_label = sys.argv[2]
with h5py.File(up_file, "r+") as h5:
    inp = h5.require_group("input")
    grp = inp.require_group("stage_parameters")
    grp.attrs["enable"] = np.int8(1)
    grp.attrs["current_stage"] = np.bytes_(stage_label)
PY
}

inject_hybrid_mapping() {
    local up_file="$1"
    local mapping_file="$2"
    python3 - "$up_file" "$mapping_file" << 'PY'
import sys
import h5py

up_file = sys.argv[1]
mapping_file = sys.argv[2]
groups = [
    "hybrid_control",
    "hybrid_bb_map",
    "hybrid_sc_map",
    "hybrid_env_topology",
]

with h5py.File(mapping_file, "r") as src, h5py.File(up_file, "r+") as dst:
    src_inp = src["/input"]
    dst_inp = dst.require_group("input")

    n_atom = int(dst["/input/pos"].shape[0])
    mem_len = int(src["/input/hybrid_env_topology/protein_membership"].shape[0])
    if n_atom != mem_len:
        raise ValueError(
            f"Hybrid mapping n_atom mismatch for {up_file}: up has {n_atom}, mapping has {mem_len}"
        )

    for g in groups:
        if g not in src_inp:
            raise ValueError(f"Missing mapping group in {mapping_file}: /input/{g}")
        if g in dst_inp:
            del dst_inp[g]
        src.copy(src_inp[g], dst_inp, name=g)
PY
}

augment_production_rotamer_nodes() {
    local up_file="$1"
    local protein_itp="$2"
    local sidechain_lib="$3"
    python3 - "$up_file" "$protein_itp" "$sidechain_lib" << 'PY'
import os
import sys
import h5py
import numpy as np

up_file, protein_itp, sidechain_lib = sys.argv[1:4]

if not os.path.exists(up_file):
    raise SystemExit(f"ERROR: stage file not found: {up_file}")
if not os.path.exists(protein_itp):
    raise SystemExit(f"ERROR: protein ITP not found: {protein_itp}")
if not os.path.exists(sidechain_lib):
    raise SystemExit(f"ERROR: sidechain library not found: {sidechain_lib}")


def parse_itp_residue_names(path):
    resnames = []
    seen = set()
    in_atoms = False
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            low = line.lower()
            if low == "[ atoms ]" or low == "[atoms]":
                in_atoms = True
                continue
            if line.startswith("[") and line.endswith("]"):
                in_atoms = False
                continue
            if not in_atoms:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                resnr = int(parts[2])
            except ValueError:
                continue
            role = parts[4].strip().upper()
            if role != "BB":
                continue
            if resnr in seen:
                continue
            seen.add(resnr)
            resnames.append(parts[3].strip().upper())
    return resnames


def normalize_resname(name):
    name = name.upper()
    aliases = {
        "HSD": "HIS",
        "HSE": "HIS",
        "HSP": "HIS",
        "HID": "HIS",
        "HIE": "HIS",
        "HIP": "HIS",
        "CYX": "CYS",
    }
    return aliases.get(name, name)


with h5py.File(up_file, "r+") as up, h5py.File(sidechain_lib, "r") as sclib:
    inp = up["/input"]
    pot = inp["potential"]

    if "hybrid_bb_map" not in inp:
        raise ValueError("Missing /input/hybrid_bb_map in stage file")
    if "hybrid_sc_map" not in inp:
        raise ValueError("Missing /input/hybrid_sc_map in stage file")

    bb_residue_raw = inp["hybrid_bb_map"]["bb_residue_index"][:].astype(np.int32)
    bb_atom_idx_raw = inp["hybrid_bb_map"]["bb_atom_index"][:].astype(np.int32)
    if bb_residue_raw.size == 0:
        raise ValueError("hybrid_bb_map is empty")

    residue_ids = []
    residue_to_bb = {}
    seen = set()
    for resid, bb_idx in zip(bb_residue_raw.tolist(), bb_atom_idx_raw.tolist()):
        rid = int(resid)
        bi = int(bb_idx)
        if rid not in seen:
            residue_ids.append(rid)
            seen.add(rid)
        if rid not in residue_to_bb and bi >= 0:
            residue_to_bb[rid] = bi

    if not residue_ids:
        raise ValueError("No protein residues found in hybrid_bb_map")

    residue_to_sc = {}
    sc_res = inp["hybrid_sc_map"]["residue_index"][:].astype(np.int32)
    sc_proxy = inp["hybrid_sc_map"]["proxy_atom_index"][:].astype(np.int32)
    for resid, proxy in zip(sc_res.tolist(), sc_proxy.tolist()):
        rid = int(resid)
        pi = int(proxy)
        if rid not in residue_to_sc and pi >= 0:
            residue_to_sc[rid] = pi

    itp_resnames = [normalize_resname(x) for x in parse_itp_residue_names(protein_itp)]
    if len(itp_resnames) != len(residue_ids):
        raise ValueError(
            f"ITP/BB residue count mismatch: ITP has {len(itp_resnames)} residues, "
            f"hybrid_bb_map has {len(residue_ids)} residues"
        )

    restype_order = [x.decode("ascii") if isinstance(x, bytes) else str(x) for x in sclib["restype_order"][:]]
    restype_num = {x: i for i, x in enumerate(restype_order)}
    start_stop_bead = sclib["rotamer_start_stop_bead"][:].astype(np.int32)
    rot_center_fixed = sclib["rotamer_center_fixed"][:].astype(np.float32)
    if "rotamer_prob_fixed" in sclib:
        rot_energy_fixed = sclib["rotamer_prob_fixed"][:].astype(np.float32).reshape(-1)
    elif "rotamer_prob" in sclib:
        rot_prob = sclib["rotamer_prob"][:].astype(np.float64)
        if rot_prob.ndim != 3:
            raise ValueError(
                f"Unsupported rotamer_prob shape in {sidechain_lib}: {rot_prob.shape}"
            )
        # Convert phi/psi-dependent probabilities into a fixed one-body energy per layer.
        rot_prob_mean = np.clip(rot_prob.mean(axis=(0, 1)), 1.0e-12, None)
        rot_energy_fixed = (-np.log(rot_prob_mean)).astype(np.float32)
    else:
        raise ValueError(
            f"Missing rotamer probability tables in {sidechain_lib}: "
            "need rotamer_prob_fixed or rotamer_prob"
        )
    bead_order = [x.decode("ascii") if isinstance(x, bytes) else str(x) for x in sclib["bead_order"][:]]
    bead_num = {x: i for i, x in enumerate(bead_order)}
    pair_interaction = sclib["pair_interaction"][:].astype(np.float32)

    n_bit_rotamer = 4
    count_by_n_rot = {}
    affine_residue = []
    layer_index = []
    beadtype_seq = []
    id_seq = []

    for rid, resname in zip(residue_ids, itp_resnames):
        if resname not in restype_num:
            raise ValueError(f"Residue '{resname}' not found in sidechain library")
        rt = restype_num[resname]
        start, stop, n_bead = [int(x) for x in start_stop_bead[rt]]
        if n_bead <= 0 or stop <= start:
            continue
        n_rot = (stop - start) // n_bead
        if n_rot <= 0:
            continue

        if n_rot not in count_by_n_rot:
            count_by_n_rot[n_rot] = 0
        base_id = (count_by_n_rot[n_rot] << n_bit_rotamer) + n_rot
        count_by_n_rot[n_rot] += 1

        n_rows = stop - start
        for rel in range(n_rows):
            lid = start + rel
            layer_index.append(lid)
            affine_residue.append(rid)
            beadtype_seq.append(f"{resname}_{rel % n_bead}")
            id_seq.append((rel // n_bead) + (base_id << n_bit_rotamer))

    if not layer_index:
        raise ValueError("No placement rows generated from sidechain library")

    layer_index_arr = np.asarray(layer_index, dtype=np.int32)
    affine_residue_arr = np.asarray(affine_residue, dtype=np.int32)
    id_seq_arr = np.asarray(id_seq, dtype=np.int32)
    beadtype_seq_arr = np.asarray([np.bytes_(x) for x in beadtype_seq], dtype="S16")
    placement_scalar_data = rot_energy_fixed[layer_index_arr][:, None].astype(np.float32)

    pos = inp["pos"][:, :, 0].astype(np.float32)
    n_atom = int(pos.shape[0])

    residue_to_order = {rid: i for i, rid in enumerate(residue_ids)}
    fallback_rid = residue_ids[0]
    bb_list_ordered = [int(residue_to_bb[rid]) for rid in residue_ids]

    def unique_preserve(seq):
        out = []
        seen_local = set()
        for x in seq:
            xi = int(x)
            if xi in seen_local:
                continue
            seen_local.add(xi)
            out.append(xi)
        return out

    def pick_triplet(rid):
        i = residue_to_order[rid]
        bb = int(residue_to_bb[rid])
        scp = int(residue_to_sc.get(rid, -1))

        candidates = []
        for off in (-2, -1, 1, 2, -3, 3):
            j = i + off
            if 0 <= j < len(residue_ids):
                candidates.append(int(residue_to_bb[residue_ids[j]]))
        if scp >= 0:
            candidates.append(scp)
        candidates.extend([bb_list_ordered[0], bb_list_ordered[-1]])
        candidates = [x for x in unique_preserve(candidates) if 0 <= x < n_atom and x != bb]

        if len(candidates) < 2:
            for x in bb_list_ordered:
                xi = int(x)
                if 0 <= xi < n_atom and xi != bb and xi not in candidates:
                    candidates.append(xi)
                if len(candidates) >= 2:
                    break
        if len(candidates) < 2:
            for x in range(n_atom):
                if x != bb and x not in candidates:
                    candidates.append(x)
                if len(candidates) >= 2:
                    break
        if len(candidates) < 2:
            raise ValueError(f"Unable to build affine triplet for residue {rid}")

        best_pair = None
        best_area = -1.0
        for a in candidates:
            for c in candidates:
                if a == c:
                    continue
                area = np.linalg.norm(np.cross(pos[a] - pos[bb], pos[c] - pos[bb]))
                if np.isfinite(area) and area > best_area:
                    best_area = float(area)
                    best_pair = (a, c)

        if best_pair is None:
            return int(candidates[0]), bb, int(candidates[1])
        return int(best_pair[0]), bb, int(best_pair[1])

    n_affine = int(max(residue_ids) + 1)
    affine_atoms = np.zeros((n_affine, 3), dtype=np.int32)
    affine_ref_geom = np.zeros((n_affine, 3, 3), dtype=np.float32)

    fa, fb, fc = pick_triplet(fallback_rid)
    for rid in range(n_affine):
        affine_atoms[rid, :] = [fa, fb, fc]
        ref = pos[[fa, fb, fc], :] - pos[[fa, fb, fc], :].mean(axis=0, keepdims=True)
        affine_ref_geom[rid, :, :] = ref.astype(np.float32)

    for rid in residue_ids:
        a, b, c = pick_triplet(rid)
        for idx in (a, b, c):
            if idx < 0 or idx >= n_atom:
                raise ValueError(f"Affine atom index out of bounds: residue {rid}, index {idx}, n_atom {n_atom}")
        affine_atoms[rid, :] = [a, b, c]
        ref = pos[[a, b, c], :] - pos[[a, b, c], :].mean(axis=0, keepdims=True)
        affine_ref_geom[rid, :, :] = ref.astype(np.float32)

    for node_name in [
        "rotamer",
        "placement_fixed_scalar",
        "placement_fixed_point_vector_only",
        "affine_alignment",
    ]:
        if node_name in pot:
            del pot[node_name]

    g_aff = pot.create_group("affine_alignment")
    g_aff.attrs["arguments"] = np.array([b"pos"])
    g_aff.create_dataset("atoms", data=affine_atoms, dtype=np.int32)
    g_aff.create_dataset("ref_geom", data=affine_ref_geom, dtype=np.float32)

    g_sc = pot.create_group("placement_fixed_point_vector_only")
    g_sc.attrs["arguments"] = np.array([b"affine_alignment"])
    g_sc.create_dataset("rama_residue", data=affine_residue_arr, dtype=np.int32)
    g_sc.create_dataset("affine_residue", data=affine_residue_arr, dtype=np.int32)
    g_sc.create_dataset("layer_index", data=layer_index_arr, dtype=np.int32)
    g_sc.create_dataset("placement_data", data=rot_center_fixed[:, :6].astype(np.float32), dtype=np.float32)
    g_sc.create_dataset("beadtype_seq", data=beadtype_seq_arr)
    g_sc.create_dataset("id_seq", data=id_seq_arr, dtype=np.int32)
    g_sc.create_dataset("fix_rotamer", data=np.zeros((0, 2), dtype=np.int32), dtype=np.int32)

    g_pl = pot.create_group("placement_fixed_scalar")
    g_pl.attrs["arguments"] = np.array([b"affine_alignment"])
    g_pl.create_dataset("rama_residue", data=affine_residue_arr, dtype=np.int32)
    g_pl.create_dataset("affine_residue", data=affine_residue_arr, dtype=np.int32)
    g_pl.create_dataset("layer_index", data=layer_index_arr, dtype=np.int32)
    g_pl.create_dataset("placement_data", data=placement_scalar_data, dtype=np.float32)

    g_rot = pot.create_group("rotamer")
    g_rot.attrs["arguments"] = np.array([b"placement_fixed_point_vector_only", b"placement_fixed_scalar"])
    g_rot.attrs["integrator_level"] = np.int32(1)
    g_rot.attrs["max_iter"] = np.int32(1000)
    g_rot.attrs["tol"] = np.float32(1.0e-3)
    g_rot.attrs["damping"] = np.float32(0.4)
    g_rot.attrs["iteration_chunk_size"] = np.int32(2)

    g_pair = g_rot.create_group("pair_interaction")
    g_pair.create_dataset("interaction_param", data=pair_interaction, dtype=np.float32)
    g_pair.create_dataset("index", data=np.arange(len(beadtype_seq), dtype=np.int32), dtype=np.int32)
    g_pair.create_dataset(
        "type",
        data=np.asarray([bead_num[x] for x in beadtype_seq], dtype=np.int32),
        dtype=np.int32,
    )
    g_pair.create_dataset("id", data=id_seq_arr, dtype=np.int32)

    print(
        f"Injected production rotamer nodes into {up_file}: "
        f"n_res={len(residue_ids)} n_sc_rows={len(beadtype_seq)} n_affine={n_affine}"
    )
PY
}

prepare_hybrid_artifacts() {
    echo "=== Stage 0: Hybrid Packing + Mapping Export ==="
    prepare_protein_inputs

    python3 prepare_hybrid_system.py \
        --protein-pdb "${PROTEIN_AA_PDB}" \
        --protein-cg-pdb "${PROTEIN_CG_EFFECTIVE}" \
        --protein-itp "${PROTEIN_ITP_EFFECTIVE}" \
        --bilayer-pdb "${BILAYER_PDB}" \
        --output-dir "${HYBRID_PREP_DIR}" \
        --salt-molar "${SALT_MOLAR}" \
        --protein-lipid-cutoff "${PROTEIN_LIPID_CUTOFF}" \
        --ion-cutoff "${ION_CUTOFF}" \
        --box-padding-xy "${BOX_PADDING_XY}" \
        --box-padding-z "${BOX_PADDING_Z}" \
        --seed "${PREP_SEED}"

    if [ ! -f "${HYBRID_PACKED_PDB}" ]; then
        echo "ERROR: hybrid packed PDB not found: ${HYBRID_PACKED_PDB}"
        exit 1
    fi
    if [ ! -f "${HYBRID_MAPPING_FILE}" ]; then
        echo "ERROR: hybrid mapping file not found: ${HYBRID_MAPPING_FILE}"
        exit 1
    fi

    if [ "${HYBRID_VALIDATE}" = "1" ]; then
        python3 validate_hybrid_mapping.py "${HYBRID_MAPPING_FILE}"
    fi

    cp -f "${HYBRID_PACKED_PDB}" "${RUNTIME_PDB_FILE}"
    cp -f "${PROTEIN_ITP_EFFECTIVE}" "${RUNTIME_ITP_FILE}"
    echo "Runtime MARTINI PDB: ${RUNTIME_PDB_FILE}"
    echo "Runtime protein ITP: ${RUNTIME_ITP_FILE}"
}

prepare_stage_file() {
    local target_file="$1"
    local prepare_stage="$2"
    local npt_enable="$3"
    local barostat_type="$4"
    local lipidhead_fc="${5:-0}"
    local stage_label="${6:-minimization}"

    export UPSIDE_SIMULATION_STAGE="$prepare_stage"
    export UPSIDE_NPT_ENABLE="$npt_enable"
    export UPSIDE_BILAYER_LIPIDHEAD_FC="$lipidhead_fc"

    python3 prepare_martini.py "${RUNTIME_PDB_ID}" --stage "$prepare_stage" "$RUN_DIR"

    local prepared_tmp="${RUN_DIR}/test.input.up"
    if [ ! -f "$prepared_tmp" ]; then
        echo "ERROR: preparation failed for stage ${prepare_stage}: $prepared_tmp not found"
        exit 1
    fi

    mv -f "$prepared_tmp" "$target_file"
    inject_hybrid_mapping "$target_file" "${HYBRID_MAPPING_FILE}"
    set_stage_label "$target_file" "$stage_label"
    if [ "$stage_label" = "production" ]; then
        augment_production_rotamer_nodes "$target_file" "${PROTEIN_ITP_EFFECTIVE}" "${SIDECHAIN_LIBRARY}"
    fi

    if [ "$npt_enable" = "1" ]; then
        set_barostat_type "$target_file" "$barostat_type"
    fi
}

run_minimization_stage() {
    local stage_label="$1"
    local up_file="$2"
    local max_iter="$3"

    local log_file="${LOG_DIR}/stage_${stage_label}.log"
    echo "=== Stage ${stage_label}: Minimization ==="
    echo "Input/Output: $up_file"

    local cmd=(
        "$UPSIDE_EXECUTABLE"
        "$up_file"
        "--duration" "0"
        "--frame-interval" "1"
        "--temperature" "$TEMPERATURE"
        "--time-step" "$MIN_TIME_STEP"
        "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
        "--thermostat-interval" "$THERMOSTAT_INTERVAL"
        "--seed" "$SEED"
        "--integrator" "v"
        "--disable-recentering"
        "--minimize"
        "--min-max-iter" "$max_iter"
        "--min-energy-tol" "1e-6"
        "--min-force-tol" "1e-3"
        "--min-step" "0.01"
    )

    if ! "${cmd[@]}" 2>&1 | tee "$log_file"; then
        echo "ERROR: Stage ${stage_label} failed"
        exit 1
    fi
}

run_md_stage() {
    local stage_label="$1"
    local input_file="$2"
    local output_file="$3"
    local nsteps="$4"
    local dt="$5"
    local frame_steps="$6"

    local effective_frame_steps="$frame_steps"

    if (( effective_frame_steps >= nsteps )); then
        effective_frame_steps=$(( nsteps / 10 ))
        if (( effective_frame_steps < 1 )); then
            effective_frame_steps=1
        fi
        echo "NOTICE: frame_steps (${frame_steps}) >= nsteps (${nsteps}); using frame_steps=${effective_frame_steps}"
    fi

    # Frame interval remains time-based in CLI. Duration uses step-count override.
    local frame_interval
    frame_interval="$(awk -v n="$effective_frame_steps" -v dt="$dt" 'BEGIN{printf "%.10g", n*dt}')"

    if [ "$input_file" != "$output_file" ]; then
        cp -f "$input_file" "$output_file"
        python3 set_initial_position.py "$input_file" "$output_file"
    fi

    local log_file="${LOG_DIR}/stage_${stage_label}.log"
    echo "=== Stage ${stage_label}: MD ==="
    echo "Input:  $input_file"
    echo "Output: $output_file"
    echo "nsteps=${nsteps}, dt=${dt}, duration(steps)=${nsteps}, frame_steps=${effective_frame_steps}, frame_interval(time)=${frame_interval}"

    local cmd=(
        "$UPSIDE_EXECUTABLE"
        "$output_file"
        "--duration-steps" "$nsteps"
        "--frame-interval" "$frame_interval"
        "--temperature" "$TEMPERATURE"
        "--time-step" "$dt"
        "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
        "--thermostat-interval" "$THERMOSTAT_INTERVAL"
        "--seed" "$SEED"
        "--integrator" "v"
        "--disable-recentering"
    )

    if ! "${cmd[@]}" 2>&1 | tee "$log_file"; then
        echo "ERROR: Stage ${stage_label} failed"
        exit 1
    fi
}

extract_stage_vtf() {
    local stage_label="$1"
    local stage_file="$2"
    local mode="$3"
    local vtf_file="${RUN_DIR}/${PDB_ID}.stage_${stage_label}.vtf"

    echo "=== Stage ${stage_label}: VTF Extraction (mode ${mode}) ==="
    echo "Input:  $stage_file"
    echo "Output: $vtf_file"
    python3 extract_martini_vtf.py "$stage_file" "$vtf_file" "$stage_file" "$RUNTIME_PDB_ID" --mode "$mode"
}

echo "=== Hybrid 1RKL Dry MARTINI Workflow ==="
echo "Protein ID: $PDB_ID"
echo "Runtime PDB ID: $RUNTIME_PDB_ID"
echo "Hybrid prep: $HYBRID_PREP_DIR"
echo "Simulation stages: 6.0 -> 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5 -> 6.6 -> 7.0"
echo

prepare_hybrid_artifacts

# 6.0: soft-core minimization (pre-production / hybrid inactive)
prepare_stage_file "$PREPARED_60_FILE" "minimization" "1" "0" "0" "minimization"
cp -f "$PREPARED_60_FILE" "$STAGE_60_FILE"
run_minimization_stage "6.0" "$STAGE_60_FILE" "$MIN_60_MAX_ITER"
extract_stage_vtf "6.0" "$STAGE_60_FILE" "1"

# 6.1: hard minimization (pre-production / hybrid inactive)
prepare_stage_file "$PREPARED_61_FILE" "npt_prod" "1" "0" "0" "minimization"
cp -f "$PREPARED_61_FILE" "$STAGE_61_FILE"
python3 set_initial_position.py "$STAGE_60_FILE" "$STAGE_61_FILE"
run_minimization_stage "6.1" "$STAGE_61_FILE" "$MIN_61_MAX_ITER"
extract_stage_vtf "6.1" "$STAGE_61_FILE" "1"

# 6.2: soft equilibration
prepare_stage_file "$STAGE_62_FILE" "npt_equil" "1" "0" "200" "minimization"
python3 set_initial_position.py "$STAGE_61_FILE" "$STAGE_62_FILE"
run_md_stage "6.2" "$STAGE_62_FILE" "$STAGE_62_FILE" "$EQ_62_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.2" "$STAGE_62_FILE" "1"

# 6.3: reduced softening equilibration
prepare_stage_file "$PREPARED_63_FILE" "npt_equil_reduced" "1" "0" "100" "minimization"
cp -f "$PREPARED_63_FILE" "$STAGE_63_FILE"
python3 set_initial_position.py "$STAGE_62_FILE" "$STAGE_63_FILE"
run_md_stage "6.3" "$STAGE_63_FILE" "$STAGE_63_FILE" "$EQ_63_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.3" "$STAGE_63_FILE" "1"

# 6.4-6.6: hard equilibration with restraint ramp
prepare_stage_file "$PREPARED_64_FILE" "npt_prod" "1" "0" "50" "minimization"
cp -f "$PREPARED_64_FILE" "$STAGE_64_FILE"
python3 set_initial_position.py "$STAGE_63_FILE" "$STAGE_64_FILE"
run_md_stage "6.4" "$STAGE_64_FILE" "$STAGE_64_FILE" "$EQ_64_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.4" "$STAGE_64_FILE" "1"

prepare_stage_file "$PREPARED_65_FILE" "npt_prod" "1" "0" "20" "minimization"
cp -f "$PREPARED_65_FILE" "$STAGE_65_FILE"
python3 set_initial_position.py "$STAGE_64_FILE" "$STAGE_65_FILE"
run_md_stage "6.5" "$STAGE_65_FILE" "$STAGE_65_FILE" "$EQ_65_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.5" "$STAGE_65_FILE" "1"

prepare_stage_file "$PREPARED_66_FILE" "npt_prod" "1" "0" "10" "minimization"
cp -f "$PREPARED_66_FILE" "$STAGE_66_FILE"
python3 set_initial_position.py "$STAGE_65_FILE" "$STAGE_66_FILE"
run_md_stage "6.6" "$STAGE_66_FILE" "$STAGE_66_FILE" "$EQ_66_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.6" "$STAGE_66_FILE" "1"

# 7.0: production (hybrid active)
prepare_stage_file "$PREPARED_70_FILE" "npt_prod" "$PROD_70_NPT_ENABLE" "$PROD_70_BAROSTAT_TYPE" "0" "production"
cp -f "$PREPARED_70_FILE" "$STAGE_70_FILE"
python3 set_initial_position.py "$STAGE_66_FILE" "$STAGE_70_FILE"
run_md_stage "7.0" "$STAGE_70_FILE" "$STAGE_70_FILE" "$PROD_70_NSTEPS" "$PROD_TIME_STEP" "$PROD_FRAME_STEPS"
extract_stage_vtf "7.0" "$STAGE_70_FILE" "2"

echo
echo "=== Workflow Complete ==="
echo "Hybrid prep:"
echo "  Packed PDB: ${HYBRID_PACKED_PDB}"
echo "  Mapping:    ${HYBRID_MAPPING_FILE}"
echo "Checkpoints:"
echo "  6.0: $STAGE_60_FILE"
echo "  6.1: $STAGE_61_FILE"
echo "  6.2: $STAGE_62_FILE"
echo "  6.3: $STAGE_63_FILE"
echo "  6.4: $STAGE_64_FILE"
echo "  6.5: $STAGE_65_FILE"
echo "  6.6: $STAGE_66_FILE"
echo "  7.0: $STAGE_70_FILE"
echo "VTF:"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.0.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.1.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.2.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.3.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.4.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.5.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.6.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_7.0.vtf"
