#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CALLER_UPSIDE_HOME="${UPSIDE_HOME:-}"
if [ "${UPSIDE_SKIP_SOURCE_SH:-0}" != "1" ]; then
    source "${PROJECT_ROOT}/source.sh"
    if [ -n "${CALLER_UPSIDE_HOME}" ]; then
        export UPSIDE_HOME="${CALLER_UPSIDE_HOME}"
    fi
elif [ -n "${CALLER_UPSIDE_HOME}" ]; then
    export UPSIDE_HOME="${CALLER_UPSIDE_HOME}"
else
    export UPSIDE_HOME="${PROJECT_ROOT}"
fi
if [ -d "${PROJECT_ROOT}/.venv/bin" ]; then
    export PATH="${PROJECT_ROOT}/.venv/bin:$PATH"
fi
export PATH="${PROJECT_ROOT}/obj:$PATH"
export PYTHONPATH="${PROJECT_ROOT}/py${PYTHONPATH:+:$PYTHONPATH}"
set -euo pipefail
cd "${SCRIPT_DIR}"
PYTHON_WORKFLOW_DIR="${PYTHON_WORKFLOW_DIR:-${PROJECT_ROOT}/py}"

generate_random_seed() {
    local seed=""
    if [ -r /dev/urandom ] && command -v od >/dev/null 2>&1; then
        seed="$(od -An -N4 -tu4 /dev/urandom 2>/dev/null || true)"
        seed="${seed//[[:space:]]/}"
    fi
    if [ -z "$seed" ] || [ "$seed" = "0" ]; then
        seed="$(( ((RANDOM & 32767) << 17) ^ ((RANDOM & 32767) << 2) ^ (RANDOM & 3) ))"
    fi
    if [ -z "$seed" ] || [ "$seed" = "0" ]; then
        seed="1"
    fi
    printf '%s\n' "$seed"
}

if [ -z "${PREP_SEED:-}" ]; then
    PREP_SEED="$(generate_random_seed)"
fi
if [ -z "${SEED:-}" ]; then
    SEED="$(generate_random_seed)"
    if [ "$SEED" = "$PREP_SEED" ]; then
        SEED="$(generate_random_seed)"
    fi
fi

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
PROTEIN_CG_PDB="${PROTEIN_CG_PDB:-}"
PROTEIN_ITP="${PROTEIN_ITP:-}"
BILAYER_PDB="${BILAYER_PDB:-${UPSIDE_HOME}/parameters/dryMARTINI/DOPC.pdb}"
UNIVERSAL_PREP_SCRIPT="${UNIVERSAL_PREP_SCRIPT:-${PYTHON_WORKFLOW_DIR}/martini_prepare_system.py}"
UNIVERSAL_PREP_MODE="${UNIVERSAL_PREP_MODE:-both}"

RUNTIME_PDB_FILE="${RUNTIME_PDB_FILE:-${HYBRID_PREP_DIR}/${RUNTIME_PDB_ID}.MARTINI.pdb}"
RUNTIME_ITP_FILE="${RUNTIME_ITP_FILE:-${HYBRID_PREP_DIR}/${RUNTIME_PDB_ID}_proa.itp}"
HYBRID_MAPPING_FILE="${HYBRID_MAPPING_FILE:-${HYBRID_PREP_DIR}/hybrid_mapping.h5}"
HYBRID_PACKED_PDB="${HYBRID_PACKED_PDB:-${HYBRID_PREP_DIR}/hybrid_packed.MARTINI.pdb}"
HYBRID_VALIDATE="${HYBRID_VALIDATE:-1}"
HYBRID_PREPROD_ACTIVATION_STAGE="${HYBRID_PREPROD_ACTIVATION_STAGE:-__hybrid_disabled__}"

# martinize controls (MARTINI 2.x generation for dry-MARTINI compatibility)
MARTINIZE_ENABLE="${MARTINIZE_ENABLE:-1}"
MARTINIZE_FF="${MARTINIZE_FF:-martini22}"
MARTINIZE_MOLNAME="${MARTINIZE_MOLNAME:-PROA}"
MARTINIZE_SCRIPT="${MARTINIZE_SCRIPT:-${PYTHON_WORKFLOW_DIR}/martini_martinize.py}"
EXTRACT_VTF_SCRIPT="${EXTRACT_VTF_SCRIPT:-${PYTHON_WORKFLOW_DIR}/martini_extract_vtf.py}"

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
CONTINUE_STAGE_70_FROM="${CONTINUE_STAGE_70_FROM:-}"
CONTINUE_STAGE_70_OUTPUT="${CONTINUE_STAGE_70_OUTPUT:-${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.continue.up}"
CONTINUE_STAGE_70_LABEL="${CONTINUE_STAGE_70_LABEL:-7.0_continue}"

# Stage-0 hybrid structure prep controls
SALT_MOLAR="${SALT_MOLAR:-0.15}"
# More conservative default exclusion to prevent lipids packed too close to protein surface.
PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-4.5}"
ION_CUTOFF="${ION_CUTOFF:-4.0}"
XY_SCALE="${XY_SCALE:-1.0}"
BOX_PADDING_XY="${BOX_PADDING_XY:-0.0}"
BOX_PADDING_Z="${BOX_PADDING_Z:-20.0}"
PROTEIN_PLACEMENT_MODE="${PROTEIN_PLACEMENT_MODE:-embed}"
PROTEIN_ORIENTATION_MODE="${PROTEIN_ORIENTATION_MODE:-input}"
PROTEIN_SURFACE_GAP="${PROTEIN_SURFACE_GAP:-6.0}"
PROTEIN_LIPID_MIN_GAP="${PROTEIN_LIPID_MIN_GAP:-4.5}"
PROTEIN_LIPID_CUTOFF_STEP="${PROTEIN_LIPID_CUTOFF_STEP:-0.5}"
PROTEIN_LIPID_CUTOFF_MAX="${PROTEIN_LIPID_CUTOFF_MAX:-8.0}"
BB_AA_MIN_MATCHED_RESIDUES="${BB_AA_MIN_MATCHED_RESIDUES:-8}"
BB_AA_MAX_RIGID_RMSD="${BB_AA_MAX_RIGID_RMSD:-1.5}"

# Simulation controls
TEMPERATURE="${TEMPERATURE:-0.8647}"
# Match the standard Upside examples unless explicitly overridden.
THERMOSTAT_TIMESCALE="${THERMOSTAT_TIMESCALE:-5.0}"
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
STRICT_STAGE_HANDOFF="${STRICT_STAGE_HANDOFF:-1}"

MIN_60_MAX_ITER="${MIN_60_MAX_ITER:-500}"
MIN_61_MAX_ITER="${MIN_61_MAX_ITER:-500}"
EQ_62_NSTEPS="${EQ_62_NSTEPS:-500}"
EQ_63_NSTEPS="${EQ_63_NSTEPS:-500}"
EQ_64_NSTEPS="${EQ_64_NSTEPS:-500}"
EQ_65_NSTEPS="${EQ_65_NSTEPS:-500}"
EQ_66_NSTEPS="${EQ_66_NSTEPS:-500}"
PROD_70_NSTEPS="${PROD_70_NSTEPS:-10000}"

EQ_TIME_STEP="${EQ_TIME_STEP:-0.010}"
# Production runs with injected Upside all-atom backbone nodes; use a smaller
# stable timestep than MARTINI-only stages.
PROD_TIME_STEP="${PROD_TIME_STEP:-0.002}"
# Leave the stage-7 protein backbone non-rigid by default in the hybrid
# workflow. Set this to 1 to re-enable the production fix_rigid mask over
# dry-MARTINI BB beads and AA backbone carrier roles BB/N/CA/C/O.
PROD_70_BACKBONE_FIX_RIGID_ENABLE="${PROD_70_BACKBONE_FIX_RIGID_ENABLE:-0}"
MIN_TIME_STEP="${MIN_TIME_STEP:-0.010}"

EQ_FRAME_STEPS="${EQ_FRAME_STEPS:-1000}"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-50}"

# AGENTS.md: 1 E_up = 2.914952774272 kJ/mol
export UPSIDE_MARTINI_ENERGY_CONVERSION="${UPSIDE_MARTINI_ENERGY_CONVERSION:-2.914952774272}"
# Dry-MARTINI native lengths are in nm; simulation inputs use Angstrom.
export UPSIDE_MARTINI_LENGTH_CONVERSION="${UPSIDE_MARTINI_LENGTH_CONVERSION:-10}"

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

export UPSIDE_EWALD_ENABLE="${UPSIDE_EWALD_ENABLE:-1}"
export UPSIDE_EWALD_ALPHA="${UPSIDE_EWALD_ALPHA:-0.2}"
export UPSIDE_EWALD_KMAX="${UPSIDE_EWALD_KMAX:-5}"
export UPSIDE_MARTINI_FF_DIR="${UPSIDE_MARTINI_FF_DIR:-${UPSIDE_HOME}/parameters/dryMARTINI}"
MASS_FF_FILE="${MASS_FF_FILE:-${UPSIDE_MARTINI_FF_DIR}/dry_martini_v2.1.itp}"
SC_MARTINI_LIBRARY="${SC_MARTINI_LIBRARY:-${UPSIDE_HOME}/parameters/ff_2.1/martini.h5}"
SC_MARTINI_TABLE_JSON_DEFAULT="${SC_MARTINI_TABLE_JSON_DEFAULT:-${UPSIDE_HOME}/SC-training/runs/default/results/assembled/sc_table.json}"
SC_MARTINI_TABLE_JSON="${SC_MARTINI_TABLE_JSON:-${SC_MARTINI_TABLE_JSON_DEFAULT}}"
UPSIDE_RAMA_LIBRARY="${UPSIDE_RAMA_LIBRARY:-${UPSIDE_HOME}/parameters/common/rama.dat}"
UPSIDE_RAMA_SHEET_MIXING="${UPSIDE_RAMA_SHEET_MIXING:-${UPSIDE_HOME}/parameters/ff_2.1/sheet}"
UPSIDE_HBOND_ENERGY="${UPSIDE_HBOND_ENERGY:-${UPSIDE_HOME}/parameters/ff_2.1/hbond.h5}"
UPSIDE_REFERENCE_STATE_RAMA="${UPSIDE_REFERENCE_STATE_RAMA:-${UPSIDE_HOME}/parameters/common/rama_reference.pkl}"
# Initial SC pair-force caps applied when production hybrid coupling turns on.
SC_ENV_LJ_FORCE_CAP="${SC_ENV_LJ_FORCE_CAP:-25.0}"
SC_ENV_COUL_FORCE_CAP="${SC_ENV_COUL_FORCE_CAP:-25.0}"
# Number of production steps used to ramp startup SC-env / SC-BB / env-BB
# LJ+Coulomb forces from capped to uncapped, leaving the last 50 steps of the
# 200-step coupling window fully regular while PO4 z remains held for 150 steps.
SC_ENV_RELAX_STEPS="${SC_ENV_RELAX_STEPS:-150}"
SC_ENV_BACKBONE_HOLD_STEPS="${SC_ENV_BACKBONE_HOLD_STEPS:-200}"
SC_ENV_PO4_Z_HOLD_STEPS="${SC_ENV_PO4_Z_HOLD_STEPS:-150}"
SC_ENV_RELAX_DT="${SC_ENV_RELAX_DT:-0.002}"
SC_ENV_RESTRAINT_K="${SC_ENV_RESTRAINT_K:-5.0}"
SC_ENV_MAX_DISPLACEMENT="${SC_ENV_MAX_DISPLACEMENT:-2.0}"
SC_ENV_PO4_Z_CLAMP_ENABLE="${SC_ENV_PO4_Z_CLAMP_ENABLE:-1}"
SC_ENV_PO4_Z_CLAMP_MODE="${SC_ENV_PO4_Z_CLAMP_MODE:-initial}"
SC_ENV_ENERGY_DUMP_ENABLE="${SC_ENV_ENERGY_DUMP_ENABLE:-1}"
SC_ENV_ENERGY_DUMP_STRIDE="${SC_ENV_ENERGY_DUMP_STRIDE:-1}"
# Keep the environment on the dry-MARTINI LJ/Coulomb Hamiltonian in hybrid
# production. The repulsive-only WCA replacement destabilizes stage 7.
PRODUCTION_NONPROTEIN_HARD_SPHERE="${PRODUCTION_NONPROTEIN_HARD_SPHERE:-0}"
PROTEIN_ENV_INTERFACE_SCALE="${PROTEIN_ENV_INTERFACE_SCALE:-1.00}"
NONPROTEIN_HS_FORCE_CAP="${NONPROTEIN_HS_FORCE_CAP:-100.0}"
NONPROTEIN_HS_POTENTIAL_CAP="${NONPROTEIN_HS_POTENTIAL_CAP:-5000.0}"

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
        echo "ERROR: MARTINI workflow martinize script not found: ${MARTINIZE_SCRIPT}"
        exit 1
    fi
fi
if [ "${MARTINIZE_ENABLE}" = "0" ]; then
    if [ -z "${PROTEIN_CG_PDB}" ]; then
        echo "ERROR: MARTINIZE_ENABLE=0 requires PROTEIN_CG_PDB to be set explicitly"
        exit 1
    fi
    if [ -z "${PROTEIN_ITP}" ]; then
        echo "ERROR: MARTINIZE_ENABLE=0 requires PROTEIN_ITP to be set explicitly"
        exit 1
    fi
    if [ ! -f "${PROTEIN_CG_PDB}" ]; then
        echo "ERROR: protein CG PDB not found: ${PROTEIN_CG_PDB}"
        exit 1
    fi
    if [ ! -f "${PROTEIN_ITP}" ]; then
        echo "ERROR: protein ITP not found: ${PROTEIN_ITP}"
        exit 1
    fi
fi
if [ ! -f "${UNIVERSAL_PREP_SCRIPT}" ]; then
    echo "ERROR: universal prep script not found: ${UNIVERSAL_PREP_SCRIPT}"
    exit 1
fi
if [ ! -f "${EXTRACT_VTF_SCRIPT}" ]; then
    echo "ERROR: MARTINI VTF extractor script not found: ${EXTRACT_VTF_SCRIPT}"
    exit 1
fi

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR" "$HYBRID_PREP_DIR"

PROTEIN_CG_EFFECTIVE="${PROTEIN_CG_PDB}"
PROTEIN_ITP_EFFECTIVE="${PROTEIN_ITP}"

ensure_sc_martini_library() {
    local need_build="0"
    if [ -f "${SC_MARTINI_LIBRARY}" ]; then
        if python3 - "${SC_MARTINI_LIBRARY}" << 'PY'
import sys
import h5py

required = [
    "grid_nm",
    "cos_theta_grid",
    "rotamer_count",
    "rotamer_probability_fixed",
    "rotamer_radial_energy_kj_mol",
    "rotamer_angular_energy_kj_mol",
    "rotamer_angular_profile",
]

path = sys.argv[1]
with h5py.File(path, "r") as h5:
    missing = [name for name in required if name not in h5]
    if missing:
        raise SystemExit(1)
PY
        then
            return
        fi
        echo "NOTICE: ${SC_MARTINI_LIBRARY} is missing rotamer-resolved SC datasets; rebuilding it."
        need_build="1"
    else
        need_build="1"
    fi

    if [ ! -f "${SC_MARTINI_TABLE_JSON}" ]; then
        echo "ERROR: SC MARTINI table JSON not found: ${SC_MARTINI_TABLE_JSON}"
        echo "The stage-7 SC injector now requires the rotamer-resolved martini.h5 schema."
        exit 1
    fi
    if [ "${need_build}" = "1" ]; then
        echo "Building ${SC_MARTINI_LIBRARY} from ${SC_MARTINI_TABLE_JSON}"
    fi
    python3 "${UNIVERSAL_PREP_SCRIPT}" build-sc-martini-h5 \
        --sc-table-json "${SC_MARTINI_TABLE_JSON}" \
        --output-h5 "${SC_MARTINI_LIBRARY}"
    if ! python3 - "${SC_MARTINI_LIBRARY}" << 'PY'
import sys
import h5py

required = [
    "grid_nm",
    "cos_theta_grid",
    "rotamer_count",
    "rotamer_probability_fixed",
    "rotamer_radial_energy_kj_mol",
    "rotamer_angular_energy_kj_mol",
    "rotamer_angular_profile",
]

path = sys.argv[1]
with h5py.File(path, "r") as h5:
    missing = [name for name in required if name not in h5]
    if missing:
        raise SystemExit(
            "rebuilt SC MARTINI library is still missing required datasets: "
            + ",".join(missing)
        )
PY
    then
        echo "ERROR: rebuilt SC MARTINI library failed schema validation: ${SC_MARTINI_LIBRARY}"
        exit 1
    fi
}

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

    echo "Running MARTINI martinize command:"
    echo "  python3 ${script_abs} -f ${input_abs} -x ${cg_abs} -o ${top_abs} -ff ${MARTINIZE_FF} -name ${MARTINIZE_MOLNAME}"
    if ! (cd "$workdir" && python3 "$script_abs" -f "$input_abs" -x "$cg_abs" -o "$top_abs" -ff "$MARTINIZE_FF" -name "$MARTINIZE_MOLNAME"); then
        echo "ERROR: MARTINI martinize execution failed."
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
        "Use MARTINI2-compatible MARTINI martinize settings (e.g., MARTINIZE_FF=martini22)."
    )

print(f"Protein ITP mass compatibility OK: {itp_file}")
PY
}

prepare_protein_inputs() {
    local mass_ff_path
    mass_ff_path="$(python3 - "$MASS_FF_FILE" << 'PY'
import os, sys
print(os.path.abspath(os.path.expanduser(sys.argv[1])))
PY
)"
    if [ "${MARTINIZE_ENABLE}" = "1" ]; then
        local martinize_dir="${RUN_DIR}/martinize"
        mkdir -p "${martinize_dir}"
        local cg_pdb="${martinize_dir}/${RUNTIME_PDB_ID}.MARTINI.pdb"
        local top_file="${martinize_dir}/${RUNTIME_PDB_ID}.top"

        run_martinize "${PROTEIN_AA_PDB}" "${cg_pdb}" "${top_file}"
        if [ ! -f "${cg_pdb}" ]; then
            echo "ERROR: MARTINI martinize did not produce CG PDB: ${cg_pdb}"
            exit 1
        fi
        if [ ! -f "${top_file}" ]; then
            echo "ERROR: MARTINI martinize did not produce topology file: ${top_file}"
            exit 1
        fi

        local generated_itp
        generated_itp="$(resolve_itp_from_top "${top_file}")"
        if [ -z "${generated_itp}" ] || [ ! -f "${generated_itp}" ]; then
            echo "ERROR: unable to resolve protein ITP from MARTINI martinize output top: ${top_file}"
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

set_hybrid_control_mode() {
    local up_file="$1"
    local activation_stage="$2"
    local preprod_mode="${3:-rigid}"
    python3 - "$up_file" "$activation_stage" "$preprod_mode" << 'PY'
import sys
import h5py
import numpy as np
up_file = sys.argv[1]
activation_stage = sys.argv[2]
preprod_mode = sys.argv[3]
with h5py.File(up_file, "r+") as h5:
    grp = h5.require_group("input").require_group("hybrid_control")
    grp.attrs["enable"] = np.int8(1)
    grp.attrs["activation_stage"] = np.bytes_(activation_stage)
    grp.attrs["preprod_protein_mode"] = np.bytes_(preprod_mode)
PY
}

set_hybrid_activation_stage() {
    local up_file="$1"
    local activation_stage="$2"
    set_hybrid_control_mode "$up_file" "$activation_stage" "rigid"
}

assert_hybrid_stage_active() {
    local up_file="$1"
    local expected_stage="$2"
    local expected_activation="$3"
    python3 - "$up_file" "$expected_stage" "$expected_activation" << 'PY'
import sys
import h5py
import numpy as np

def as_text(v):
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode("utf-8", errors="ignore")
    return str(v)

up_file = sys.argv[1]
expected_stage = sys.argv[2]
expected_activation = sys.argv[3]
with h5py.File(up_file, "r") as h5:
    inp = h5["/input"]
    stage = as_text(inp["stage_parameters"].attrs.get("current_stage", b"")).strip()
    hy = inp["hybrid_control"].attrs
    enable = int(hy.get("enable", 0))
    activation = as_text(hy.get("activation_stage", b"")).strip()
    if stage != expected_stage:
        raise SystemExit(
            f"ERROR: {up_file} current_stage={stage!r} expected {expected_stage!r}"
        )
    if enable != 1:
        raise SystemExit(
            f"ERROR: {up_file} hybrid_control.enable={enable} expected 1"
        )
    if activation != expected_activation:
        raise SystemExit(
            f"ERROR: {up_file} activation_stage={activation!r} expected {expected_activation!r}"
        )
print(
    f"Hybrid activation verified: stage={expected_stage}, activation_stage={expected_activation}, enable=1"
)
PY
}

set_hybrid_production_controls() {
    local up_file="$1"
    local nonprotein_hard_sphere="$2"
    local protein_env_interface_scale="$3"
    python3 - "$up_file" "$nonprotein_hard_sphere" "$protein_env_interface_scale" << 'PY'
import sys
import math
import h5py
import numpy as np

up_file = sys.argv[1]
nonprotein_hard_sphere = int(sys.argv[2])
protein_env_interface_scale = float(sys.argv[3])

if not math.isfinite(protein_env_interface_scale) or protein_env_interface_scale <= 0.0:
    raise SystemExit(
        f"ERROR: protein_env_interface_scale must be finite and > 0, got {protein_env_interface_scale!r}"
    )

with h5py.File(up_file, "r+") as h5:
    grp = h5.require_group("input").require_group("hybrid_control")
    grp.attrs["production_nonprotein_hard_sphere"] = np.int8(1 if nonprotein_hard_sphere else 0)
    grp.attrs["protein_env_interface_scale"] = np.float32(protein_env_interface_scale)
PY
}

set_production_backbone_fix_rigid() {
    local up_file="$1"
    python3 - "$up_file" << 'PY'
import sys
import h5py
import numpy as np

up_file = sys.argv[1]
required_roles = ("BB", "N", "CA", "C", "O")

def as_text(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)

with h5py.File(up_file, "r+") as h5:
    inp = h5.require_group("input")
    if "hybrid_env_topology" not in inp:
        raise SystemExit(f"ERROR: missing /input/hybrid_env_topology in {up_file}")
    if "pos" not in inp:
        raise SystemExit(f"ERROR: missing /input/pos in {up_file}")

    membership = inp["hybrid_env_topology"]["protein_membership"][:].astype(np.int32, copy=False)
    n_atom = int(inp["pos"].shape[0])
    if membership.shape[0] != n_atom:
        raise SystemExit(
            f"ERROR: {up_file} protein_membership length {membership.shape[0]} does not match n_atom {n_atom}"
        )

    role_source = ""
    if "atom_roles" in inp:
        role_source = "atom_roles"
    elif "atom_names" in inp:
        role_source = "atom_names"
    else:
        raise SystemExit(f"ERROR: missing /input/atom_roles or /input/atom_names in {up_file}")

    role_values = [as_text(x).strip().upper() for x in inp[role_source][:]]
    if len(role_values) != n_atom:
        raise SystemExit(
            f"ERROR: {up_file} {role_source} length {len(role_values)} does not match n_atom {n_atom}"
        )

    role_counts = {role: 0 for role in required_roles}
    selected = []
    for atom_idx, role in enumerate(role_values):
        if membership[atom_idx] < 0 or role not in role_counts:
            continue
        selected.append(atom_idx)
        role_counts[role] += 1

    missing_roles = [role for role, count in role_counts.items() if count == 0]
    if missing_roles:
        raise SystemExit(
            f"ERROR: {up_file} missing protein backbone roles for fix_rigid selection: "
            f"{','.join(missing_roles)}"
        )

    bb_count = role_counts["BB"]
    mismatched_roles = [role for role in required_roles[1:] if role_counts[role] != bb_count]
    if mismatched_roles:
        detail = ", ".join(f"{role}={role_counts[role]}" for role in required_roles)
        raise SystemExit(
            f"ERROR: {up_file} protein backbone role counts are inconsistent for fix_rigid selection ({detail})"
        )

    existing = np.zeros((0,), dtype=np.int32)
    if "fix_rigid" in inp and "atom_indices" in inp["fix_rigid"]:
        existing = inp["fix_rigid"]["atom_indices"][:].astype(np.int32, copy=False)

    merged = np.unique(np.concatenate([existing, np.asarray(selected, dtype=np.int32)]))

    if "fix_rigid" in inp:
        del inp["fix_rigid"]
    grp = inp.require_group("fix_rigid")
    grp.attrs["enable"] = np.int8(1)
    grp.attrs["selection"] = np.bytes_("protein_backbone_roles_bb_n_ca_c_o")
    grp.attrs["selection_source"] = np.bytes_(role_source)
    grp.attrs["selection_count"] = np.int32(merged.shape[0])
    for role, count in role_counts.items():
        grp.attrs[f"count_{role}"] = np.int32(count)
    grp.create_dataset("atom_indices", data=merged.astype(np.int32, copy=False))
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
    "hybrid_env_topology",
]

BB_COMPONENT_MASS = {
    "N": 14.0,
    "CA": 12.0,
    "C": 12.0,
    "O": 16.0,
}
BB_BEAD_MASS = 72.0
BB_COMPONENT_SUM = float(sum(BB_COMPONENT_MASS.values()))
BB_COMPONENT_SCALE = (BB_BEAD_MASS / BB_COMPONENT_SUM) if BB_COMPONENT_SUM > 0.0 else 1.0
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
                comp_mass = float(BB_COMPONENT_MASS.get(cname, 12.0))
                ref_mass[ref_idx] = mass_dtype.type((BB_COMPONENT_SCALE * comp_mass) / 12.0)

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

    membership = np.full((n_atom_aug,), -1, dtype=np.int32)
    membership[:base_n_atom] = src_mem
    if n_ref_index > 0:
        membership[ref_offset:ref_offset + n_ref_index] = 0
    replace_dataset(env_grp, "protein_membership", membership)

    if env_grp["protein_membership"].shape[0] != dst_inp["pos"].shape[0]:
        raise ValueError(f"{up_file}: hybrid_env_topology/protein_membership length mismatch after augmentation")
PY
}

inject_stage7_sc_table_nodes() {
    local up_file="$1"
    local sc_library="$2"
    python3 "${UNIVERSAL_PREP_SCRIPT}" inject-stage7-sc \
        "$up_file" \
        "$sc_library" \
        "${UPSIDE_HOME}" \
        "${UPSIDE_RAMA_LIBRARY}" \
        "${UPSIDE_RAMA_SHEET_MIXING}" \
        "${UPSIDE_HBOND_ENERGY}" \
        "${UPSIDE_REFERENCE_STATE_RAMA}" \
        --protein-itp "${PROTEIN_ITP_EFFECTIVE}"
}

prepare_hybrid_artifacts() {
    echo "=== Stage 0: Hybrid Packing + Mapping Export ==="
    prepare_protein_inputs

    local packing_cutoff="${PROTEIN_LIPID_CUTOFF}"
    local packing_min_gap="nan"

    while true; do
        echo "Hybrid packing attempt with protein-lipid cutoff: ${packing_cutoff} Å"

        python3 "${UNIVERSAL_PREP_SCRIPT}" \
            --mode "both" \
            --pdb-id "${RUNTIME_PDB_ID}" \
            --runtime-pdb-output "${HYBRID_PACKED_PDB}" \
            --runtime-itp-output "${RUNTIME_ITP_FILE}" \
            --prepare-structure 1 \
            --protein-aa-pdb "${PROTEIN_AA_PDB}" \
            --protein-cg-pdb "${PROTEIN_CG_EFFECTIVE}" \
            --protein-itp "${PROTEIN_ITP_EFFECTIVE}" \
            --hybrid-mapping-output "${HYBRID_MAPPING_FILE}" \
            --hybrid-bb-map-json-output "${HYBRID_PREP_DIR}/hybrid_bb_map.json" \
            --bilayer-pdb "${BILAYER_PDB}" \
            --salt-molar "${SALT_MOLAR}" \
            --protein-lipid-cutoff "${packing_cutoff}" \
            --ion-cutoff "${ION_CUTOFF}" \
            --xy-scale "${XY_SCALE}" \
            --box-padding-xy "${BOX_PADDING_XY}" \
            --box-padding-z "${BOX_PADDING_Z}" \
            --protein-placement-mode "${PROTEIN_PLACEMENT_MODE}" \
            --protein-orientation-mode "${PROTEIN_ORIENTATION_MODE}" \
            --protein-surface-gap "${PROTEIN_SURFACE_GAP}" \
            --bb-aa-min-matched-residues "${BB_AA_MIN_MATCHED_RESIDUES}" \
            --bb-aa-max-rigid-rmsd "${BB_AA_MAX_RIGID_RMSD}" \
            --seed "${PREP_SEED}" \
            --summary-json "${HYBRID_PREP_DIR}/hybrid_prep_summary.json"

        if [ ! -f "${HYBRID_PACKED_PDB}" ]; then
            echo "ERROR: hybrid packed PDB not found: ${HYBRID_PACKED_PDB}"
            exit 1
        fi

        packing_min_gap="$(python3 - "${HYBRID_PACKED_PDB}" << 'PY'
import sys
import numpy as np

pdb_file = sys.argv[1]
protein_residues = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX"
}
lipid_residues = {"DOP", "DOPC"}

protein_xyz = []
lipid_xyz = []
with open(pdb_file, "r", encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:21].strip().upper()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        if resname in protein_residues:
            protein_xyz.append((x, y, z))
        elif resname in lipid_residues:
            lipid_xyz.append((x, y, z))

if not protein_xyz or not lipid_xyz:
    raise SystemExit("nan")

protein_xyz = np.asarray(protein_xyz, dtype=float)
lipid_xyz = np.asarray(lipid_xyz, dtype=float)

min_d2 = float("inf")
chunk = 2000
for i in range(0, lipid_xyz.shape[0], chunk):
    block = lipid_xyz[i:i + chunk]
    d = block[:, None, :] - protein_xyz[None, :, :]
    d2 = np.einsum("ijk,ijk->ij", d, d)
    block_min = float(d2.min())
    if block_min < min_d2:
        min_d2 = block_min

print(f"{min_d2 ** 0.5:.6f}")
PY
)"

        if python3 - "${packing_min_gap}" "${PROTEIN_LIPID_MIN_GAP}" << 'PY'
import sys
observed = float(sys.argv[1])
target = float(sys.argv[2])
raise SystemExit(0 if observed >= target else 1)
PY
        then
            echo "Hybrid packing accepted: min protein-lipid distance ${packing_min_gap} Å (target >= ${PROTEIN_LIPID_MIN_GAP} Å)"
            break
        fi

        if python3 - "${packing_cutoff}" "${PROTEIN_LIPID_CUTOFF_MAX}" << 'PY'
import sys
cutoff = float(sys.argv[1])
limit = float(sys.argv[2])
raise SystemExit(0 if cutoff >= limit else 1)
PY
        then
            echo "ERROR: hybrid packing still overpacked near protein."
            echo "  Observed min protein-lipid distance: ${packing_min_gap} Å"
            echo "  Target min distance: ${PROTEIN_LIPID_MIN_GAP} Å"
            echo "  Reached cutoff limit: ${PROTEIN_LIPID_CUTOFF_MAX} Å"
            echo "  Try increasing PROTEIN_LIPID_CUTOFF_MAX and/or BOX_PADDING_XY."
            exit 1
        fi

        packing_cutoff="$(python3 - "${packing_cutoff}" "${PROTEIN_LIPID_CUTOFF_STEP}" << 'PY'
import sys
cutoff = float(sys.argv[1])
step = float(sys.argv[2])
print(f"{cutoff + step:.3f}")
PY
)"
        echo "Hybrid packing too tight near protein (min ${packing_min_gap} Å). Retrying with cutoff ${packing_cutoff} Å."
    done

    if [ ! -f "${HYBRID_MAPPING_FILE}" ]; then
        echo "ERROR: hybrid mapping file not found: ${HYBRID_MAPPING_FILE}"
        exit 1
    fi

    if [ "${HYBRID_VALIDATE}" = "1" ]; then
        python3 "${UNIVERSAL_PREP_SCRIPT}" validate-hybrid-mapping "${HYBRID_MAPPING_FILE}"
    fi

    if [ "$(python3 - "${HYBRID_PACKED_PDB}" "${RUNTIME_PDB_FILE}" << 'PY'
import os
import sys
print(int(os.path.abspath(sys.argv[1]) == os.path.abspath(sys.argv[2])))
PY
)" = "0" ]; then
        cp -f "${HYBRID_PACKED_PDB}" "${RUNTIME_PDB_FILE}"
    fi
    if [ "$(python3 - "${PROTEIN_ITP_EFFECTIVE}" "${RUNTIME_ITP_FILE}" << 'PY'
import os
import sys
print(int(os.path.abspath(sys.argv[1]) == os.path.abspath(sys.argv[2])))
PY
)" = "0" ]; then
        cp -f "${PROTEIN_ITP_EFFECTIVE}" "${RUNTIME_ITP_FILE}"
    fi
    echo "Runtime MARTINI PDB: ${RUNTIME_PDB_FILE}"
    echo "Runtime protein ITP: ${RUNTIME_ITP_FILE}"
}

prepare_stage_file() {
    local target_file="$1"
    local prepare_stage="$2"
    local npt_enable="$3"
    local lipidhead_fc="${4:-0}"
    local stage_label="${5:-minimization}"

    export UPSIDE_SIMULATION_STAGE="$prepare_stage"
    export UPSIDE_NPT_ENABLE="$npt_enable"
    export UPSIDE_BILAYER_LIPIDHEAD_FC="$lipidhead_fc"

    python3 "${UNIVERSAL_PREP_SCRIPT}" \
        --mode "${UNIVERSAL_PREP_MODE}" \
        --pdb-id "${RUNTIME_PDB_ID}" \
        --runtime-pdb-output "${RUNTIME_PDB_FILE}" \
        --runtime-itp-output "${RUNTIME_ITP_FILE}" \
        --prepare-structure 0 \
        --stage "$prepare_stage" \
        --run-dir "$RUN_DIR" \
        --summary-json "${HYBRID_PREP_DIR}/stage_${prepare_stage}.summary.json"

    local prepared_tmp="${RUN_DIR}/test.input.up"
    if [ ! -f "$prepared_tmp" ]; then
        echo "ERROR: preparation failed for stage ${prepare_stage}: $prepared_tmp not found"
        exit 1
    fi

    mv -f "$prepared_tmp" "$target_file"
    inject_hybrid_mapping "$target_file" "${HYBRID_MAPPING_FILE}"
    set_hybrid_activation_stage "$target_file" "${HYBRID_PREPROD_ACTIVATION_STAGE}"
    set_stage_label "$target_file" "$stage_label"
    if [ "$stage_label" = "production" ]; then
        ensure_sc_martini_library
        set_hybrid_activation_stage "$target_file" "production"
        set_hybrid_production_controls \
            "$target_file" \
            "$PRODUCTION_NONPROTEIN_HARD_SPHERE" \
            "$PROTEIN_ENV_INTERFACE_SCALE"
        inject_stage7_sc_table_nodes \
            "$target_file" \
            "${SC_MARTINI_LIBRARY}"
        if [ "${PROD_70_BACKBONE_FIX_RIGID_ENABLE}" = "1" ]; then
            set_production_backbone_fix_rigid "$target_file"
        fi
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
        handoff_initial_position "$input_file" "$output_file"
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

handoff_initial_position() {
    local input_file="$1"
    local output_file="$2"
    local mode="${3:-default}"

    local refresh_hybrid="0"
    local recenter_production="0"
    if [ "$mode" = "production_hybrid" ]; then
        # Hybrid activates at stage 7.0; refresh AA carrier coordinates onto
        # current BB anchors to avoid a first-step force/energy spike.
        refresh_hybrid="1"
        recenter_production="0"
    fi

    UPSIDE_SET_INITIAL_STRICT_COPY="$STRICT_STAGE_HANDOFF" \
    UPSIDE_SET_INITIAL_REFRESH_HYBRID_CARRIERS="$refresh_hybrid" \
    UPSIDE_SET_INITIAL_RECENTER_PRODUCTION="$recenter_production" \
        python3 "${UNIVERSAL_PREP_SCRIPT}" set-initial-position "$input_file" "$output_file"
}

list_stage_vtf_outputs() {
    local stage_label="$1"
    local stage_file="$2"

    python3 - "$RUN_DIR" "$PDB_ID" "$stage_label" "$stage_file" << 'PY'
import os
import re
import sys

import h5py

run_dir, pdb_id, stage_label, stage_file = sys.argv[1:5]
prefix = os.path.join(run_dir, f"{pdb_id}.stage_{stage_label}")
pattern = re.compile(r"^output_previous_(\d+)$")

with h5py.File(stage_file, "r") as f:
    groups = []
    for key in f.keys():
        match = pattern.fullmatch(str(key))
        if match and isinstance(f[key], h5py.Group):
            groups.append((int(match.group(1)), str(key)))
    groups.sort()
    output_groups = [name for _, name in groups]
    if "output" in f and isinstance(f["output"], h5py.Group):
        output_groups.append("output")

if len(output_groups) <= 1:
    print(f"{prefix}.vtf")
else:
    for segment_index, _ in enumerate(output_groups):
        print(f"{prefix}.segment_{segment_index}.vtf")
PY
}

print_stage_vtf_outputs() {
    local stage_label="$1"
    local stage_file="$2"
    while IFS= read -r vtf_file; do
        [ -n "$vtf_file" ] || continue
        echo "  ${vtf_file}"
    done < <(list_stage_vtf_outputs "$stage_label" "$stage_file")
}

extract_stage_vtf() {
    local stage_label="$1"
    local stage_file="$2"
    local mode="$3"
    local vtf_file="${RUN_DIR}/${PDB_ID}.stage_${stage_label}.vtf"

    echo "=== Stage ${stage_label}: VTF Extraction (mode ${mode}) ==="
    echo "Input:  $stage_file"
    echo "Output prefix: $vtf_file"
    python3 "${EXTRACT_VTF_SCRIPT}" "$stage_file" "$vtf_file" "$stage_file" "$RUNTIME_PDB_ID" --mode "$mode" --split-segments
}

run_stage70_continuation() {
    local source_file="$1"
    local output_file="$2"
    local stage_label="${3:-7.0_continue}"

    if [ ! -f "$source_file" ]; then
        echo "ERROR: continuation source not found: $source_file"
        exit 1
    fi

    assert_hybrid_stage_active "$source_file" "production" "production"

    mkdir -p "$(dirname "$output_file")"

    if [ "$(python3 - "$source_file" "$output_file" << 'PY'
import os
import sys
print(int(os.path.abspath(sys.argv[1]) == os.path.abspath(sys.argv[2])))
PY
)" = "0" ]; then
        cp -f "$source_file" "$output_file"
    fi

    handoff_initial_position "$source_file" "$output_file" "production_hybrid"
    assert_hybrid_stage_active "$output_file" "production" "production"
    run_md_stage "$stage_label" "$output_file" "$output_file" "$PROD_70_NSTEPS" "$PROD_TIME_STEP" "$PROD_FRAME_STEPS"
    extract_stage_vtf "$stage_label" "$output_file" "2"
}

echo "=== Hybrid 1RKL Dry MARTINI Workflow ==="
echo "Protein ID: $PDB_ID"
echo "Runtime PDB ID: $RUNTIME_PDB_ID"
echo "Preparation seed: ${PREP_SEED}"
echo "Simulation seed: ${SEED}"
echo "Universal prep: ${UNIVERSAL_PREP_SCRIPT} (mode=${UNIVERSAL_PREP_MODE})"
echo "VTF extractor: ${EXTRACT_VTF_SCRIPT}"
echo "Hybrid prep: $HYBRID_PREP_DIR"
if [ -n "$CONTINUE_STAGE_70_FROM" ]; then
    echo "Continuation mode: production stage 7.0 only"
    echo "Continuation source: ${CONTINUE_STAGE_70_FROM}"
    echo "Continuation output: ${CONTINUE_STAGE_70_OUTPUT}"
else
    echo "Simulation stages: 6.0 -> 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5 -> 6.6 -> 7.0"
fi
echo

if [ -n "$CONTINUE_STAGE_70_FROM" ]; then
    run_stage70_continuation "$CONTINUE_STAGE_70_FROM" "$CONTINUE_STAGE_70_OUTPUT" "$CONTINUE_STAGE_70_LABEL"
else
    prepare_hybrid_artifacts

    # 6.0: soft-core minimization (pre-production / hybrid inactive)
    set_stage_npt_targets "6.0"
    prepare_stage_file "$PREPARED_60_FILE" "minimization" "1" "0" "minimization"
    cp -f "$PREPARED_60_FILE" "$STAGE_60_FILE"
    run_minimization_stage "6.0" "$STAGE_60_FILE" "$MIN_60_MAX_ITER"
    extract_stage_vtf "6.0" "$STAGE_60_FILE" "1"

    # 6.1: hard minimization (pre-production / hybrid inactive)
    set_stage_npt_targets "6.1"
    prepare_stage_file "$PREPARED_61_FILE" "npt_prod" "1" "0" "minimization"
    cp -f "$PREPARED_61_FILE" "$STAGE_61_FILE"
    handoff_initial_position "$STAGE_60_FILE" "$STAGE_61_FILE"
    run_minimization_stage "6.1" "$STAGE_61_FILE" "$MIN_61_MAX_ITER"
    extract_stage_vtf "6.1" "$STAGE_61_FILE" "1"

    # 6.2: soft equilibration
    set_stage_npt_targets "6.2"
    prepare_stage_file "$STAGE_62_FILE" "npt_equil" "1" "200" "minimization"
    handoff_initial_position "$STAGE_61_FILE" "$STAGE_62_FILE"
    run_md_stage "6.2" "$STAGE_62_FILE" "$STAGE_62_FILE" "$EQ_62_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
    extract_stage_vtf "6.2" "$STAGE_62_FILE" "1"

    # 6.3: reduced softening equilibration
    set_stage_npt_targets "6.3"
    prepare_stage_file "$PREPARED_63_FILE" "npt_equil_reduced" "1" "100" "minimization"
    cp -f "$PREPARED_63_FILE" "$STAGE_63_FILE"
    handoff_initial_position "$STAGE_62_FILE" "$STAGE_63_FILE"
    run_md_stage "6.3" "$STAGE_63_FILE" "$STAGE_63_FILE" "$EQ_63_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
    extract_stage_vtf "6.3" "$STAGE_63_FILE" "1"

    # 6.4-6.6: hard equilibration with restraint ramp
    set_stage_npt_targets "6.4"
    prepare_stage_file "$PREPARED_64_FILE" "npt_prod" "1" "50" "minimization"
    cp -f "$PREPARED_64_FILE" "$STAGE_64_FILE"
    handoff_initial_position "$STAGE_63_FILE" "$STAGE_64_FILE"
    run_md_stage "6.4" "$STAGE_64_FILE" "$STAGE_64_FILE" "$EQ_64_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
    extract_stage_vtf "6.4" "$STAGE_64_FILE" "1"

    set_stage_npt_targets "6.5"
    prepare_stage_file "$PREPARED_65_FILE" "npt_prod" "1" "20" "minimization"
    cp -f "$PREPARED_65_FILE" "$STAGE_65_FILE"
    handoff_initial_position "$STAGE_64_FILE" "$STAGE_65_FILE"
    run_md_stage "6.5" "$STAGE_65_FILE" "$STAGE_65_FILE" "$EQ_65_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
    extract_stage_vtf "6.5" "$STAGE_65_FILE" "1"

    set_stage_npt_targets "6.6"
    prepare_stage_file "$PREPARED_66_FILE" "npt_prod" "1" "10" "minimization"
    cp -f "$PREPARED_66_FILE" "$STAGE_66_FILE"
    handoff_initial_position "$STAGE_65_FILE" "$STAGE_66_FILE"
    run_md_stage "6.6" "$STAGE_66_FILE" "$STAGE_66_FILE" "$EQ_66_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
    extract_stage_vtf "6.6" "$STAGE_66_FILE" "1"

    # 7.0: production (hybrid active)
    prepare_stage_file "$PREPARED_70_FILE" "npt_prod" "$PROD_70_NPT_ENABLE" "0" "production"
    cp -f "$PREPARED_70_FILE" "$STAGE_70_FILE"
    handoff_initial_position "$STAGE_66_FILE" "$STAGE_70_FILE" "production_hybrid"
    assert_hybrid_stage_active "$STAGE_70_FILE" "production" "production"
    run_md_stage "7.0" "$STAGE_70_FILE" "$STAGE_70_FILE" "$PROD_70_NSTEPS" "$PROD_TIME_STEP" "$PROD_FRAME_STEPS"
    extract_stage_vtf "7.0" "$STAGE_70_FILE" "2"
fi

echo
echo "=== Workflow Complete ==="
if [ -n "$CONTINUE_STAGE_70_FROM" ]; then
    echo "Continuation source: ${CONTINUE_STAGE_70_FROM}"
    echo "Continuation checkpoint: ${CONTINUE_STAGE_70_OUTPUT}"
    echo "Continuation VTF:"
    print_stage_vtf_outputs "$CONTINUE_STAGE_70_LABEL" "$CONTINUE_STAGE_70_OUTPUT"
else
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
    print_stage_vtf_outputs "6.0" "$STAGE_60_FILE"
    print_stage_vtf_outputs "6.1" "$STAGE_61_FILE"
    print_stage_vtf_outputs "6.2" "$STAGE_62_FILE"
    print_stage_vtf_outputs "6.3" "$STAGE_63_FILE"
    print_stage_vtf_outputs "6.4" "$STAGE_64_FILE"
    print_stage_vtf_outputs "6.5" "$STAGE_65_FILE"
    print_stage_vtf_outputs "6.6" "$STAGE_66_FILE"
    print_stage_vtf_outputs "7.0" "$STAGE_70_FILE"
fi
