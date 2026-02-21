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
# Production runs with injected Upside all-atom backbone nodes; use a smaller
# stable timestep than MARTINI-only stages.
PROD_TIME_STEP="${PROD_TIME_STEP:-0.002}"
MIN_TIME_STEP="${MIN_TIME_STEP:-0.010}"

EQ_FRAME_STEPS="${EQ_FRAME_STEPS:-1000}"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-5000}"

COMP_3E4_BAR_INV_TO_A3_PER_EUP="${COMP_3E4_BAR_INV_TO_A3_PER_EUP:-14.521180763676}"
# Unit conversion from AGENTS.md:
#   1 bar = 0.000020659477 E_up / Angstrom^3
BAR_1_TO_EUP_PER_A3="${BAR_1_TO_EUP_PER_A3:-0.000020659477}"
NPT_REF_P_ZERO="${NPT_REF_P_ZERO:-0.0}"
NPT_REF_P_ONE_BAR="${NPT_REF_P_ONE_BAR:-$BAR_1_TO_EUP_PER_A3}"
export UPSIDE_NPT_TARGET_PXY="${UPSIDE_NPT_TARGET_PXY:-$NPT_REF_P_ZERO}"
export UPSIDE_NPT_TARGET_PZ="${UPSIDE_NPT_TARGET_PZ:-$NPT_REF_P_ZERO}"
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
UPSIDE_RAMA_LIBRARY="${UPSIDE_RAMA_LIBRARY:-${UPSIDE_HOME}/parameters/common/rama.dat}"
UPSIDE_RAMA_SHEET_MIXING="${UPSIDE_RAMA_SHEET_MIXING:-${UPSIDE_HOME}/parameters/ff_2.1/sheet}"
UPSIDE_HBOND_ENERGY="${UPSIDE_HBOND_ENERGY:-${UPSIDE_HOME}/parameters/ff_2.1/hbond.h5}"
UPSIDE_REFERENCE_STATE_RAMA="${UPSIDE_REFERENCE_STATE_RAMA:-${UPSIDE_HOME}/parameters/common/rama_reference.pkl}"
SC_ENV_LJ_FORCE_CAP="${SC_ENV_LJ_FORCE_CAP:-25.0}"
SC_ENV_COUL_FORCE_CAP="${SC_ENV_COUL_FORCE_CAP:-25.0}"
SC_ENV_RELAX_STEPS="${SC_ENV_RELAX_STEPS:-200}"
SC_ENV_RELAX_DT="${SC_ENV_RELAX_DT:-0.002}"
SC_ENV_RESTRAINT_K="${SC_ENV_RESTRAINT_K:-5.0}"
SC_ENV_MAX_DISPLACEMENT="${SC_ENV_MAX_DISPLACEMENT:-2.0}"
NONPROTEIN_HS_FORCE_CAP="${NONPROTEIN_HS_FORCE_CAP:-100.0}"
NONPROTEIN_HS_POTENTIAL_CAP="${NONPROTEIN_HS_POTENTIAL_CAP:-5000.0}"
INTEGRATION_RMSD_ALIGN_ENABLE="${INTEGRATION_RMSD_ALIGN_ENABLE:-1}"

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
if [ ! -f "${UPSIDE_RAMA_LIBRARY}" ]; then
    echo "ERROR: Upside rama library not found: ${UPSIDE_RAMA_LIBRARY}"
    exit 1
fi
if [ ! -f "${UPSIDE_HBOND_ENERGY}" ]; then
    echo "ERROR: Upside hbond energy file not found: ${UPSIDE_HBOND_ENERGY}"
    exit 1
fi
if [ ! -f "${UPSIDE_REFERENCE_STATE_RAMA}" ]; then
    echo "ERROR: Upside reference-state rama file not found: ${UPSIDE_REFERENCE_STATE_RAMA}"
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

set_hybrid_sc_controls() {
    local up_file="$1"
    local lj_cap="$2"
    local coul_cap="$3"
    local relax_steps="$4"
    local relax_dt="$5"
    local rest_k="$6"
    local max_disp="$7"
    local hs_force_cap="$8"
    local hs_pot_cap="$9"
    local rmsd_align_enable="${10}"
    python3 - "$up_file" "$lj_cap" "$coul_cap" "$relax_steps" "$relax_dt" "$rest_k" "$max_disp" "$hs_force_cap" "$hs_pot_cap" "$rmsd_align_enable" << 'PY'
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
hs_force_cap = float(sys.argv[8])
hs_pot_cap = float(sys.argv[9])
rmsd_align_enable = int(sys.argv[10])

with h5py.File(up_file, "r+") as h5:
    grp = h5.require_group("input").require_group("hybrid_control")
    grp.attrs["sc_env_lj_force_cap"] = np.float32(lj_cap)
    grp.attrs["sc_env_coul_force_cap"] = np.float32(coul_cap)
    grp.attrs["sc_env_relax_steps"] = np.int32(relax_steps)
    grp.attrs["sc_env_relax_dt"] = np.float32(relax_dt)
    grp.attrs["sc_env_restraint_k"] = np.float32(rest_k)
    grp.attrs["sc_env_max_displacement"] = np.float32(max_disp)
    grp.attrs["nonprotein_hs_force_cap"] = np.float32(hs_force_cap)
    grp.attrs["nonprotein_hs_potential_cap"] = np.float32(hs_pot_cap)
    grp.attrs["integration_rmsd_align_enable"] = np.int8(1 if rmsd_align_enable else 0)
PY
}

inject_hybrid_mapping() {
    local up_file="$1"
    local mapping_file="$2"
    python3 - "$up_file" "$mapping_file" << 'PY'
import sys
import h5py
import numpy as np

up_file = sys.argv[1]
mapping_file = sys.argv[2]
groups = [
    "hybrid_control",
    "hybrid_bb_map",
    "hybrid_sc_map",
    "hybrid_env_topology",
]

BB_COMPONENT_MASS = {
    "N": 14.0,
    "CA": 12.0,
    "C": 12.0,
    "O": 16.0,
}
DEFAULT_COMPONENT_NAMES = ("N", "CA", "C", "O")

def as_text(x):
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode("utf-8", errors="ignore")
    return str(x)

def pad_bytes(value, dtype):
    n = np.dtype(dtype).itemsize
    raw = str(value).encode("ascii", errors="ignore")[:n]
    return raw.ljust(n, b" ")

def replace_dataset(parent, name, data):
    if name in parent:
        del parent[name]
    parent.create_dataset(name, data=data)

with h5py.File(mapping_file, "r") as src, h5py.File(up_file, "r+") as dst:
    src_inp = src["/input"]
    dst_inp = dst.require_group("input")

    base_n_atom = int(dst["/input/pos"].shape[0])
    src_mem = src["/input/hybrid_env_topology/protein_membership"][:].astype(np.int32)
    mem_len = int(src_mem.shape[0])
    if base_n_atom != mem_len:
        raise ValueError(
            f"Hybrid mapping n_atom mismatch for {up_file}: up has {base_n_atom}, mapping has {mem_len}"
        )

    for g in groups:
        if g not in src_inp:
            raise ValueError(f"Missing mapping group in {mapping_file}: /input/{g}")
        if g in dst_inp:
            del dst_inp[g]
        src.copy(src_inp[g], dst_inp, name=g)

    bb_grp = dst_inp["hybrid_bb_map"]
    env_grp = dst_inp["hybrid_env_topology"]
    pot_grp = dst_inp["potential"]["martini_potential"]

    bb_residue = bb_grp["bb_residue_index"][:].astype(np.int32)
    bb_ref_idx = (
        bb_grp["reference_atom_indices"][:].astype(np.int32)
        if "reference_atom_indices" in bb_grp
        else bb_grp["atom_indices"][:].astype(np.int32)
    )
    n_bb = int(bb_ref_idx.shape[0])
    if bb_ref_idx.ndim != 2 or bb_ref_idx.shape[1] != 4:
        raise ValueError("hybrid_bb_map reference/active indices must have shape (n_bb,4)")

    if "reference_atom_coords" in bb_grp:
        bb_ref_xyz = bb_grp["reference_atom_coords"][:].astype(np.float32)
        if bb_ref_xyz.shape != (n_bb, 4, 3):
            raise ValueError("hybrid_bb_map/reference_atom_coords must have shape (n_bb,4,3)")
    else:
        bb_ref_xyz = np.zeros((n_bb, 4, 3), dtype=np.float32)

    comp_names = list(DEFAULT_COMPONENT_NAMES)
    if "reference_atom_names" in bb_grp:
        names_raw = bb_grp["reference_atom_names"][:]
        if names_raw.shape == (4,):
            comp_names = [as_text(x).strip().upper() or DEFAULT_COMPONENT_NAMES[i] for i, x in enumerate(names_raw)]

    valid_ref = bb_ref_idx >= 0
    max_ref_idx = int(np.max(bb_ref_idx[valid_ref])) if np.any(valid_ref) else -1
    n_ref_index = max_ref_idx + 1
    ref_offset = base_n_atom

    if n_ref_index > 0:
        pos_dtype = dst_inp["pos"].dtype
        mom_dtype = dst_inp["mom"].dtype
        vel_dtype = dst_inp["vel"].dtype
        mass_dtype = dst_inp["mass"].dtype
        charge_dtype = dst_inp["charges"].dtype
        type_dtype = dst_inp["type"].dtype
        name_dtype = dst_inp["atom_names"].dtype if "atom_names" in dst_inp else type_dtype
        role_dtype = dst_inp["atom_roles"].dtype if "atom_roles" in dst_inp else type_dtype
        residue_dtype = dst_inp["residue_ids"].dtype
        molecule_dtype = dst_inp["molecule_ids"].dtype
        pclass_dtype = dst_inp["particle_class"].dtype

        ref_pos = np.zeros((n_ref_index, 3, 1), dtype=pos_dtype)
        ref_mom = np.zeros((n_ref_index, 3, 1), dtype=mom_dtype)
        ref_vel = np.zeros((n_ref_index, 3), dtype=vel_dtype)
        ref_mass = np.ones((n_ref_index,), dtype=mass_dtype)
        ref_charge = np.zeros((n_ref_index,), dtype=charge_dtype)
        ref_type = np.full((n_ref_index,), pad_bytes("AA", type_dtype), dtype=type_dtype)
        ref_name = np.full((n_ref_index,), pad_bytes("AA", name_dtype), dtype=name_dtype)
        ref_role = np.full((n_ref_index,), pad_bytes("AA", role_dtype), dtype=role_dtype)
        ref_residue = np.full((n_ref_index,), -1, dtype=residue_dtype)
        ref_molecule = np.full((n_ref_index,), 0, dtype=molecule_dtype)
        ref_pclass = np.full((n_ref_index,), pad_bytes("PROTEINAA", pclass_dtype), dtype=pclass_dtype)

        for k in range(n_bb):
            resid = int(bb_residue[k])
            for d in range(4):
                ref_idx = int(bb_ref_idx[k, d])
                if ref_idx < 0:
                    continue
                ref_pos[ref_idx, :, 0] = bb_ref_xyz[k, d, :]
                cname = comp_names[d] if d < len(comp_names) else DEFAULT_COMPONENT_NAMES[d]
                ref_name[ref_idx] = pad_bytes(cname, name_dtype)
                ref_role[ref_idx] = pad_bytes(cname, role_dtype)
                ref_residue[ref_idx] = resid
                ref_mass[ref_idx] = mass_dtype.type(BB_COMPONENT_MASS.get(cname, 12.0) / 12.0)

        def append_input_dataset(name, appendix):
            if name not in dst_inp:
                return
            base = dst_inp[name][:]
            if base.shape[0] != base_n_atom:
                raise ValueError(f"{up_file}: /input/{name} length changed unexpectedly before augmentation")
            merged = np.concatenate([base, appendix.astype(base.dtype, copy=False)], axis=0)
            replace_dataset(dst_inp, name, merged)

        append_input_dataset("pos", ref_pos)
        append_input_dataset("mom", ref_mom)
        append_input_dataset("vel", ref_vel)
        append_input_dataset("mass", ref_mass)
        append_input_dataset("charges", ref_charge)
        append_input_dataset("type", ref_type)
        append_input_dataset("atom_names", ref_name)
        append_input_dataset("atom_roles", ref_role)
        append_input_dataset("residue_ids", ref_residue)
        append_input_dataset("molecule_ids", ref_molecule)
        append_input_dataset("particle_class", ref_pclass)

        pot_atom_indices = pot_grp["atom_indices"][:]
        if pot_atom_indices.shape[0] != base_n_atom:
            raise ValueError(f"{up_file}: /input/potential/martini_potential/atom_indices length mismatch")
        pot_ai_append = np.arange(ref_offset, ref_offset + n_ref_index, dtype=pot_atom_indices.dtype)
        replace_dataset(pot_grp, "atom_indices", np.concatenate([pot_atom_indices, pot_ai_append], axis=0))

        pot_charges = pot_grp["charges"][:]
        if pot_charges.shape[0] != base_n_atom:
            raise ValueError(f"{up_file}: /input/potential/martini_potential/charges length mismatch")
        pot_q_append = np.zeros((n_ref_index,), dtype=pot_charges.dtype)
        replace_dataset(pot_grp, "charges", np.concatenate([pot_charges, pot_q_append], axis=0))

    n_atom_aug = int(dst_inp["pos"].shape[0])
    if n_atom_aug != base_n_atom + max(0, n_ref_index):
        raise ValueError(f"{up_file}: unexpected n_atom after reference index augmentation")

    bb_runtime_idx = np.full((n_bb, 4), -1, dtype=np.int32)
    bb_runtime_mask = np.zeros((n_bb, 4), dtype=np.int8)
    bb_runtime_w = np.zeros((n_bb, 4), dtype=np.float32)
    for k in range(n_bb):
        raw_w = np.zeros((4,), dtype=np.float32)
        for d in range(4):
            ref_idx = int(bb_ref_idx[k, d])
            if ref_idx < 0:
                continue
            run_idx = ref_offset + ref_idx
            bb_runtime_idx[k, d] = run_idx
            bb_runtime_mask[k, d] = 1
            cname = comp_names[d] if d < len(comp_names) else DEFAULT_COMPONENT_NAMES[d]
            raw_w[d] = np.float32(BB_COMPONENT_MASS.get(cname, 12.0))
        wsum = float(raw_w.sum())
        if wsum > 0.0:
            bb_runtime_w[k, :] = raw_w / wsum

    replace_dataset(bb_grp, "atom_indices", bb_runtime_idx)
    replace_dataset(bb_grp, "atom_mask", bb_runtime_mask)
    replace_dataset(bb_grp, "weights", bb_runtime_w)
    bb_grp.attrs["atom_index_space"] = np.bytes_("stage_runtime")
    bb_grp.attrs["reference_index_space"] = np.bytes_("protein_aa_pdb_0based")
    bb_grp.attrs["reference_index_offset"] = np.int32(ref_offset)
    bb_grp.attrs["reference_index_count"] = np.int32(max(0, n_ref_index))

    # Rewrite SC projection targets to AA carrier indices (no MARTINI BB fallback).
    sc_grp = dst_inp["hybrid_sc_map"]
    sc_residue = sc_grp["residue_index"][:].astype(np.int32)
    n_sc = int(sc_residue.shape[0])
    sc_runtime_idx = np.full((n_sc, 4), -1, dtype=np.int32)
    sc_runtime_w = np.zeros((n_sc, 4), dtype=np.float32)

    residue_to_bb_row = {}
    for k, resid in enumerate(bb_residue.tolist()):
        residue_to_bb_row.setdefault(int(resid), int(k))

    for r in range(n_sc):
        resid = int(sc_residue[r])
        if resid not in residue_to_bb_row:
            raise ValueError(
                f"{up_file}: hybrid_sc_map residue {resid} has no matching hybrid_bb_map row"
            )
        k = residue_to_bb_row[resid]
        idx_row = bb_runtime_idx[k].astype(np.int32, copy=True)
        w_row = bb_runtime_w[k].astype(np.float32, copy=True)
        active = (idx_row >= 0) & (w_row > 0.0)
        if not np.any(active):
            raise ValueError(
                f"{up_file}: hybrid_sc_map residue {resid} maps to empty AA carrier targets"
            )
        sc_runtime_idx[r] = idx_row
        sc_runtime_w[r] = w_row

    replace_dataset(sc_grp, "proj_target_indices", sc_runtime_idx)
    replace_dataset(sc_grp, "proj_weights", sc_runtime_w)
    sc_grp.attrs["target_index_space"] = np.bytes_("stage_runtime")
    sc_grp.attrs["target_projection"] = np.bytes_("bb_component_carriers")

    membership = np.full((n_atom_aug,), -1, dtype=np.int32)
    membership[:base_n_atom] = src_mem
    if n_ref_index > 0:
        membership[ref_offset:ref_offset + n_ref_index] = 0
    replace_dataset(env_grp, "protein_membership", membership)

    if env_grp["protein_membership"].shape[0] != dst_inp["pos"].shape[0]:
        raise ValueError(f"{up_file}: hybrid_env_topology/protein_membership length mismatch after augmentation")
PY
}

augment_production_rotamer_nodes() {
    local up_file="$1"
    local protein_itp="$2"
    local sidechain_lib="$3"
    local upside_home="$4"
    local rama_library="$5"
    local rama_sheet_mixing="$6"
    local hbond_energy="$7"
    local reference_state_rama="$8"
    python3 - "$up_file" "$protein_itp" "$sidechain_lib" "$upside_home" "$rama_library" "$rama_sheet_mixing" "$hbond_energy" "$reference_state_rama" << 'PY'
import os
import sys
import h5py
import numpy as np
import types
import tables as tb
import importlib.util
import _pickle as cPickle

up_file, protein_itp, sidechain_lib, upside_home, rama_library, rama_sheet_mixing, hbond_energy, reference_state_rama = sys.argv[1:9]

if not os.path.exists(up_file):
    raise SystemExit(f"ERROR: stage file not found: {up_file}")
if not os.path.exists(protein_itp):
    raise SystemExit(f"ERROR: protein ITP not found: {protein_itp}")
if not os.path.exists(sidechain_lib):
    raise SystemExit(f"ERROR: sidechain library not found: {sidechain_lib}")
if not os.path.exists(rama_library):
    raise SystemExit(f"ERROR: Upside rama library not found: {rama_library}")
if not os.path.exists(rama_sheet_mixing):
    raise SystemExit(f"ERROR: Upside rama sheet-mixing file not found: {rama_sheet_mixing}")
if not os.path.exists(hbond_energy):
    raise SystemExit(f"ERROR: Upside hbond energy file not found: {hbond_energy}")
if not os.path.exists(reference_state_rama):
    raise SystemExit(f"ERROR: Upside reference-state rama file not found: {reference_state_rama}")

upside_config_py = os.path.join(upside_home, "py", "upside_config.py")
if not os.path.exists(upside_config_py):
    raise SystemExit(f"ERROR: upside_config.py not found: {upside_config_py}")

spec = importlib.util.spec_from_file_location("upside_config_runtime", upside_config_py)
uc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(uc)


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

def remap_atom_index_array(arr, atom_map):
    out = arr.copy()
    if out.size == 0:
        return out
    finite_mask = np.isfinite(out)
    if not np.any(finite_mask):
        return out
    nonneg_mask = finite_mask & (out >= 0)
    if not np.any(nonneg_mask):
        return out
    idx = out[nonneg_mask].astype(np.int64)
    if np.any(idx >= atom_map.shape[0]):
        raise ValueError(
            f"Backbone atom remap index exceeds reference map size: "
            f"max={int(idx.max())}, map_size={atom_map.shape[0]}"
        )
    mapped = atom_map[idx]
    if np.any(mapped < 0):
        bad = np.where(mapped < 0)[0][0]
        raise ValueError(f"Backbone atom remap produced negative target for reference index {int(idx[bad])}")
    out[nonneg_mask] = mapped.astype(out.dtype, copy=False)
    return out

BACKBONE_NODES = [
    "Distance3D",
    "Angle",
    "Dihedral_omega",
    "Dihedral_phi",
    "Dihedral_psi",
    "Spring_bond",
    "Spring_angle",
    "Spring_omega",
    "rama_coord",
    "rama_map_pot",
    "rama_map_pot_ref",
    "infer_H_O",
    "protein_hbond",
    "hbond_energy",
    "backbone_pairs",
]

CANONICAL_AFFINE_REF = np.array(
    [
        [-1.19280531, -0.83127186, 0.0],  # N
        [0.0, 0.0, 0.0],                  # CA
        [1.25222632, -0.87268266, 0.0],   # C
    ],
    dtype=np.float32,
)
CANONICAL_AFFINE_REF -= CANONICAL_AFFINE_REF.mean(axis=0, keepdims=True)


with h5py.File(up_file, "r+") as up, h5py.File(sidechain_lib, "r") as sclib:
    inp = up["/input"]
    pot = inp["potential"]

    if "hybrid_bb_map" not in inp:
        raise ValueError("Missing /input/hybrid_bb_map in stage file")
    if "hybrid_sc_map" not in inp:
        raise ValueError("Missing /input/hybrid_sc_map in stage file")

    bb_grp = inp["hybrid_bb_map"]
    bb_residue_raw = bb_grp["bb_residue_index"][:].astype(np.int32)
    bb_atom_idx_raw = bb_grp["bb_atom_index"][:].astype(np.int32)
    bb_ref_idx = (
        bb_grp["reference_atom_indices"][:].astype(np.int32)
        if "reference_atom_indices" in bb_grp
        else None
    )
    ref_offset = int(bb_grp.attrs.get("reference_index_offset", -1))
    if bb_residue_raw.size == 0:
        raise ValueError("hybrid_bb_map is empty")

    pos = inp["pos"][:, :, 0].astype(np.float32)
    n_atom = int(pos.shape[0])

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

    residue_to_ncac = {}
    if ref_offset >= 0 and bb_ref_idx is not None:
        for resid, ref_row in zip(bb_residue_raw.tolist(), bb_ref_idx.tolist()):
            rid = int(resid)
            n_idx, ca_idx, c_idx = [int(ref_row[0]), int(ref_row[1]), int(ref_row[2])]
            if n_idx < 0 or ca_idx < 0 or c_idx < 0:
                continue
            n_rt = ref_offset + n_idx
            ca_rt = ref_offset + ca_idx
            c_rt = ref_offset + c_idx
            for idx in (n_rt, ca_rt, c_rt):
                if idx < 0 or idx >= n_atom:
                    raise ValueError(
                        f"Backbone carrier index out of bounds for residue {rid}: "
                        f"idx={idx}, n_atom={n_atom}, offset={ref_offset}"
                    )
            if rid not in residue_to_ncac:
                residue_to_ncac[rid] = (n_rt, ca_rt, c_rt)
    missing_backbone = [rid for rid in residue_ids if rid not in residue_to_ncac]
    if missing_backbone:
        raise ValueError(
            f"Missing N/CA/C runtime mapping for residues: {missing_backbone[:8]} "
            f"(total {len(missing_backbone)})."
        )

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

    n_affine = int(max(residue_ids) + 1)
    affine_atoms = np.zeros((n_affine, 3), dtype=np.int32)
    affine_ref_geom = np.zeros((n_affine, 3, 3), dtype=np.float32)

    fa, fb, fc = residue_to_ncac[residue_ids[0]]
    for rid in range(n_affine):
        affine_atoms[rid, :] = [fa, fb, fc]
        affine_ref_geom[rid, :, :] = CANONICAL_AFFINE_REF

    for rid in residue_ids:
        a, b, c = residue_to_ncac[rid]
        affine_atoms[rid, :] = [a, b, c]
        affine_ref_geom[rid, :, :] = CANONICAL_AFFINE_REF

    for node_name in [
        "rotamer",
        "placement_fixed_scalar",
        "placement_fixed_point_vector_only",
        "affine_alignment",
        *BACKBONE_NODES,
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

    if "sequence" in inp:
        del inp["sequence"]
    seq_data = np.asarray([np.bytes_(x) for x in itp_resnames], dtype="S3")
    inp.create_dataset("sequence", data=seq_data)

    ref_n_atom = 3 * len(residue_ids)
    atom_map = np.full((ref_n_atom,), -1, dtype=np.int64)
    for i, rid in enumerate(residue_ids):
        n_idx, ca_idx, c_idx = residue_to_ncac[rid]
        atom_map[3 * i + 0] = n_idx
        atom_map[3 * i + 1] = ca_idx
        atom_map[3 * i + 2] = c_idx

fasta_seq = np.array(itp_resnames)
spring_args = types.SimpleNamespace(
    bond_stiffness=48.0,
    angle_stiffness=175.0,
    omega_stiffness=30.0,
)

with tb.open_file(up_file, mode="a") as tf:
    uc.t = tf
    uc.potential = tf.root.input.potential
    uc.n_atom = ref_n_atom
    uc.n_chains = 1
    uc.chain_starts = np.array([0], dtype=np.int32)
    uc.use_intensive_memory = False

    uc.write_dist_spring(spring_args)
    uc.write_angle_spring(spring_args)
    uc.write_omega_spring1(spring_args, fasta_seq)
    uc.write_rama_map_pot(
        fasta_seq,
        rama_library,
        rama_sheet_mixing,
        secstr_bias='',
        mode='mixture',
        param_deriv=False,
    )

    ref_state_cor = np.log(cPickle.load(open(reference_state_rama, "rb"), encoding="latin1"))
    ref_state_cor -= ref_state_cor.mean()
    grp = tf.create_group(tf.root.input.potential, 'rama_map_pot_ref')
    grp._v_attrs.arguments = np.array([b'rama_coord'])
    grp._v_attrs.log_pot = 0
    uc.create_array(grp, 'residue_id', obj=np.arange(len(fasta_seq), dtype=np.int32))
    uc.create_array(grp, 'rama_map_id', obj=np.zeros(len(fasta_seq), dtype=np.int32))
    uc.create_array(grp, 'rama_pot', obj=ref_state_cor[None])

    uc.write_infer_H_O(fasta_seq, np.array([], dtype=np.int32))
    uc.write_count_hbond(fasta_seq, False)
    uc.write_short_hbond(fasta_seq, hbond_energy)
    uc.write_rama_coord2()
    uc.write_backbone_pair(fasta_seq)

with h5py.File(up_file, "r+") as up:
    pot = up["/input/potential"]
    for node_name in BACKBONE_NODES:
        if node_name not in pot:
            raise ValueError(f"Missing generated backbone node in stage file {up_file}: {node_name}")

    remap_datasets = [
        ("Distance3D", "id"),
        ("Angle", "id"),
        ("Dihedral_omega", "id"),
        ("Dihedral_phi", "id"),
        ("Dihedral_psi", "id"),
        ("rama_coord", "id"),
        ("infer_H_O/donors", "id"),
        ("infer_H_O/acceptors", "id"),
    ]
    for grp_name, dset_name in remap_datasets:
        ds = pot[grp_name][dset_name]
        ds[...] = remap_atom_index_array(ds[:], atom_map).astype(ds.dtype, copy=False)

    print(
        f"Injected production rotamer nodes into {up_file}: "
            f"n_res={len(residue_ids)} n_sc_rows={len(beadtype_seq)} "
            f"n_affine={n_affine} backbone_nodes={len(BACKBONE_NODES)}"
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
        set_hybrid_sc_controls \
            "$target_file" \
            "$SC_ENV_LJ_FORCE_CAP" \
            "$SC_ENV_COUL_FORCE_CAP" \
            "$SC_ENV_RELAX_STEPS" \
            "$SC_ENV_RELAX_DT" \
            "$SC_ENV_RESTRAINT_K" \
            "$SC_ENV_MAX_DISPLACEMENT" \
            "$NONPROTEIN_HS_FORCE_CAP" \
            "$NONPROTEIN_HS_POTENTIAL_CAP" \
            "$INTEGRATION_RMSD_ALIGN_ENABLE"
        augment_production_rotamer_nodes \
            "$target_file" \
            "${PROTEIN_ITP_EFFECTIVE}" \
            "${SIDECHAIN_LIBRARY}" \
            "${UPSIDE_HOME}" \
            "${UPSIDE_RAMA_LIBRARY}" \
            "${UPSIDE_RAMA_SHEET_MIXING}" \
            "${UPSIDE_HBOND_ENERGY}" \
            "${UPSIDE_REFERENCE_STATE_RAMA}"
    fi

    if [ "$npt_enable" = "1" ]; then
        set_barostat_type "$target_file" "$barostat_type"
    fi
}

set_stage_npt_targets() {
    local stage_label="$1"
    case "$stage_label" in
        6.0|6.1)
            # CHARMM-GUI step6.0/6.1 minimization.mdp: ref_p = 0.0 0.0
            export UPSIDE_NPT_TARGET_PXY="$NPT_REF_P_ZERO"
            export UPSIDE_NPT_TARGET_PZ="$NPT_REF_P_ZERO"
            ;;
        6.2|6.3|6.4|6.5|6.6)
            # CHARMM-GUI step6.2-6.6 equilibration.mdp: ref_p defaults to 1 bar.
            # semi-isotropic coupling keeps z compressibility at 0.0, so z scaling remains disabled.
            export UPSIDE_NPT_TARGET_PXY="$NPT_REF_P_ONE_BAR"
            export UPSIDE_NPT_TARGET_PZ="$NPT_REF_P_ONE_BAR"
            ;;
        *)
            return
            ;;
    esac
    echo "NPT targets for stage ${stage_label}: Pxy=${UPSIDE_NPT_TARGET_PXY}, Pz=${UPSIDE_NPT_TARGET_PZ} (E_up/Angstrom^3)"
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
set_stage_npt_targets "6.0"
prepare_stage_file "$PREPARED_60_FILE" "minimization" "1" "0" "0" "minimization"
cp -f "$PREPARED_60_FILE" "$STAGE_60_FILE"
run_minimization_stage "6.0" "$STAGE_60_FILE" "$MIN_60_MAX_ITER"
extract_stage_vtf "6.0" "$STAGE_60_FILE" "1"

# 6.1: hard minimization (pre-production / hybrid inactive)
set_stage_npt_targets "6.1"
prepare_stage_file "$PREPARED_61_FILE" "npt_prod" "1" "0" "0" "minimization"
cp -f "$PREPARED_61_FILE" "$STAGE_61_FILE"
python3 set_initial_position.py "$STAGE_60_FILE" "$STAGE_61_FILE"
run_minimization_stage "6.1" "$STAGE_61_FILE" "$MIN_61_MAX_ITER"
extract_stage_vtf "6.1" "$STAGE_61_FILE" "1"

# 6.2: soft equilibration
set_stage_npt_targets "6.2"
prepare_stage_file "$STAGE_62_FILE" "npt_equil" "1" "0" "200" "minimization"
python3 set_initial_position.py "$STAGE_61_FILE" "$STAGE_62_FILE"
run_md_stage "6.2" "$STAGE_62_FILE" "$STAGE_62_FILE" "$EQ_62_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.2" "$STAGE_62_FILE" "1"

# 6.3: reduced softening equilibration
set_stage_npt_targets "6.3"
prepare_stage_file "$PREPARED_63_FILE" "npt_equil_reduced" "1" "0" "100" "minimization"
cp -f "$PREPARED_63_FILE" "$STAGE_63_FILE"
python3 set_initial_position.py "$STAGE_62_FILE" "$STAGE_63_FILE"
run_md_stage "6.3" "$STAGE_63_FILE" "$STAGE_63_FILE" "$EQ_63_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.3" "$STAGE_63_FILE" "1"

# 6.4-6.6: hard equilibration with restraint ramp
set_stage_npt_targets "6.4"
prepare_stage_file "$PREPARED_64_FILE" "npt_prod" "1" "0" "50" "minimization"
cp -f "$PREPARED_64_FILE" "$STAGE_64_FILE"
python3 set_initial_position.py "$STAGE_63_FILE" "$STAGE_64_FILE"
run_md_stage "6.4" "$STAGE_64_FILE" "$STAGE_64_FILE" "$EQ_64_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.4" "$STAGE_64_FILE" "1"

set_stage_npt_targets "6.5"
prepare_stage_file "$PREPARED_65_FILE" "npt_prod" "1" "0" "20" "minimization"
cp -f "$PREPARED_65_FILE" "$STAGE_65_FILE"
python3 set_initial_position.py "$STAGE_64_FILE" "$STAGE_65_FILE"
run_md_stage "6.5" "$STAGE_65_FILE" "$STAGE_65_FILE" "$EQ_65_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.5" "$STAGE_65_FILE" "1"

set_stage_npt_targets "6.6"
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
