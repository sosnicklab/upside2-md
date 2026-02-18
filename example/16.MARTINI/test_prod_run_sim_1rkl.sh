#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../source.sh"
source "${SCRIPT_DIR}/../../.venv/bin/activate"
set -euo pipefail
cd "${SCRIPT_DIR}"

# Production-only runner for 1RKL hybrid workflow.
# Assumes stages 6.0-6.6 and hybrid prep artifacts already exist.

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

RUNTIME_PDB_FILE="${RUNTIME_PDB_FILE:-${SCRIPT_DIR}/pdb/${RUNTIME_PDB_ID}.MARTINI.pdb}"
RUNTIME_ITP_FILE="${RUNTIME_ITP_FILE:-${SCRIPT_DIR}/pdb/${RUNTIME_PDB_ID}_proa.itp}"
HYBRID_MAPPING_FILE="${HYBRID_MAPPING_FILE:-${HYBRID_PREP_DIR}/hybrid_mapping.h5}"

STAGE_66_FILE="${STAGE_66_FILE:-${CHECKPOINT_DIR}/${PDB_ID}.stage_6.6.up}"
PREPARED_70_FILE="${PREPARED_70_FILE:-${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.prepared.up}"
STAGE_70_FILE="${STAGE_70_FILE:-${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.up}"

# Simulation controls
TEMPERATURE="${TEMPERATURE:-0.8647}"
THERMOSTAT_TIMESCALE="${THERMOSTAT_TIMESCALE:-4.0}"
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
SEED="${SEED:-7090685331}"
PROD_70_NSTEPS="${PROD_70_NSTEPS:-5000}"
PROD_TIME_STEP="${PROD_TIME_STEP:-0.020}"
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
SIDECHAIN_LIBRARY="${SIDECHAIN_LIBRARY:-${UPSIDE_HOME}/parameters/ff_2.1/sidechain.h5}"
SC_ENV_LJ_FORCE_CAP="${SC_ENV_LJ_FORCE_CAP:-25.0}"
SC_ENV_COUL_FORCE_CAP="${SC_ENV_COUL_FORCE_CAP:-25.0}"
SC_ENV_RELAX_STEPS="${SC_ENV_RELAX_STEPS:-200}"
SC_ENV_RELAX_DT="${SC_ENV_RELAX_DT:-0.002}"
SC_ENV_RESTRAINT_K="${SC_ENV_RESTRAINT_K:-5.0}"
SC_ENV_MAX_DISPLACEMENT="${SC_ENV_MAX_DISPLACEMENT:-2.0}"

if [ -z "${UPSIDE_HOME:-}" ]; then
    echo "ERROR: UPSIDE_HOME environment variable is not set"
    exit 1
fi

UPSIDE_EXECUTABLE="${UPSIDE_HOME}/obj/upside"
if [ ! -f "$UPSIDE_EXECUTABLE" ]; then
    echo "ERROR: UPSIDE executable not found: $UPSIDE_EXECUTABLE"
    exit 1
fi

if [ ! -f "${RUNTIME_PDB_FILE}" ]; then
    echo "ERROR: runtime MARTINI PDB not found: ${RUNTIME_PDB_FILE}"
    echo "Run full run_sim_1rkl.sh once to generate hybrid runtime assets."
    exit 1
fi
if [ ! -f "${RUNTIME_ITP_FILE}" ]; then
    echo "ERROR: runtime protein ITP not found: ${RUNTIME_ITP_FILE}"
    echo "Run full run_sim_1rkl.sh once to generate hybrid runtime assets."
    exit 1
fi
if [ ! -f "${HYBRID_MAPPING_FILE}" ]; then
    echo "ERROR: hybrid mapping file not found: ${HYBRID_MAPPING_FILE}"
    echo "Run full run_sim_1rkl.sh once to generate hybrid prep outputs."
    exit 1
fi
if [ ! -f "${STAGE_66_FILE}" ]; then
    echo "ERROR: stage 6.6 checkpoint not found: ${STAGE_66_FILE}"
    echo "This script assumes stages 6.0-6.6 have already completed."
    exit 1
fi

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR"

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

set_hybrid_sc_controls() {
    local up_file="$1"
    local lj_cap="$2"
    local coul_cap="$3"
    local relax_steps="$4"
    local relax_dt="$5"
    local rest_k="$6"
    local max_disp="$7"
    python3 - "$up_file" "$lj_cap" "$coul_cap" "$relax_steps" "$relax_dt" "$rest_k" "$max_disp" << 'PY'
import sys
import h5py
import numpy as np

up_file = sys.argv[1]
lj_cap = float(sys.argv[2])
coul_cap = float(sys.argv[3])
relax_steps = int(sys.argv[4])
relax_dt = float(sys.argv[5])
rest_k = float(sys.argv[6])
max_disp = float(sys.argv[7])

with h5py.File(up_file, "r+") as h5:
    grp = h5.require_group("input").require_group("hybrid_control")
    grp.attrs["sc_env_lj_force_cap"] = np.float32(lj_cap)
    grp.attrs["sc_env_coul_force_cap"] = np.float32(coul_cap)
    grp.attrs["sc_env_relax_steps"] = np.int32(relax_steps)
    grp.attrs["sc_env_relax_dt"] = np.float32(relax_dt)
    grp.attrs["sc_env_restraint_k"] = np.float32(rest_k)
    grp.attrs["sc_env_max_displacement"] = np.float32(max_disp)
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
        for rel in range(stop - start):
            lid = start + rel
            layer_index.append(lid)
            affine_residue.append(rid)
            beadtype_seq.append(f"{resname}_{rel % n_bead}")
            id_seq.append((rel // n_bead) + (base_id << n_bit_rotamer))

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

prepare_stage_file() {
    local target_file="$1"
    local prepare_stage="$2"
    local npt_enable="$3"
    local barostat_type="$4"
    local lipidhead_fc="${5:-0}"
    local stage_label="${6:-production}"

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
        set_hybrid_sc_controls \
            "$target_file" \
            "$SC_ENV_LJ_FORCE_CAP" \
            "$SC_ENV_COUL_FORCE_CAP" \
            "$SC_ENV_RELAX_STEPS" \
            "$SC_ENV_RELAX_DT" \
            "$SC_ENV_RESTRAINT_K" \
            "$SC_ENV_MAX_DISPLACEMENT"
        augment_production_rotamer_nodes "$target_file" "${RUNTIME_ITP_FILE}" "${SIDECHAIN_LIBRARY}"
    fi

    if [ "$npt_enable" = "1" ]; then
        set_barostat_type "$target_file" "$barostat_type"
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

echo "=== Hybrid 1RKL Production-Only Runner ==="
echo "Protein ID: $PDB_ID"
echo "Runtime PDB ID: $RUNTIME_PDB_ID"
echo "Assumed completed checkpoint: ${STAGE_66_FILE}"
echo

prepare_stage_file "$PREPARED_70_FILE" "npt_prod" "$PROD_70_NPT_ENABLE" "$PROD_70_BAROSTAT_TYPE" "0" "production"
cp -f "$PREPARED_70_FILE" "$STAGE_70_FILE"
python3 set_initial_position.py "$STAGE_66_FILE" "$STAGE_70_FILE"
run_md_stage "7.0" "$STAGE_70_FILE" "$STAGE_70_FILE" "$PROD_70_NSTEPS" "$PROD_TIME_STEP" "$PROD_FRAME_STEPS"
extract_stage_vtf "7.0" "$STAGE_70_FILE" "2"

echo
echo "=== Production-Only Run Complete ==="
echo "Checkpoint:"
echo "  7.0: $STAGE_70_FILE"
echo "VTF:"
echo "  ${RUN_DIR}/${PDB_ID}.stage_7.0.vtf"
