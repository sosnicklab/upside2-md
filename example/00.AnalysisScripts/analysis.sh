#!/bin/bash

set -euo pipefail

# ============================================
# User settings
# ============================================
#

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SELF_PATH="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"

# CHECKME: update these for the target run mode, workflow directory, and repository install.
RUNNER="local" #local or slurm
WORK_DIR="./"
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/../.." && pwd) # Repository root containing source.sh, py/, and obj/.
PYTHON_ENV_CMD="" # Leave empty to use PROJECT_ROOT/.venv when present, or set e.g. "module load anaconda3".

# CHECKME: update these for the target system and experimental condition.
PDB_ID="glpG-RKRK-79HIS"
SIM_ID="memb_test"
N_REP="48"
START_FRAME="100"

# CHECKME: set SKIP_EXPERIMENT_DATA="true" if no experimental data available, otherwise set all parameters for experimental data
SKIP_EXPERIMENT_DATA="false"
HXMS_METHOD="stretch_exp"
PROTEIN_STATE="pd9"
EXP_DATA_FILE="GlpG psWT Sub final peptides up sum 11192024.csv"

# CHECKME: optional legacy HX plot step. Use auto to run only when the legacy inputs are present.
HX_PLOT_ENABLED="auto" # auto, true, or false
HX_PLOT_PREFIX="glpg"
HX_PLOT_STATE="${PROTEIN_STATE}"
HX_PLOT_DFOUT_FILE=""
HX_PLOT_FITDATA_FILE=""
HX_PLOT_DG_FILE=""
HX_PLOT_RESID_FILE=""
HX_PLOT_OUTPUT_DIR=""

# ============================================
# Slurm settings
# ============================================

# CHECKME: adjust these only when RUNNER=slurm.
SLURM_JOB_NAME="analysis_${PDB_ID}"
SLURM_TIME="08:00:00"
SLURM_CPUS_PER_TASK="4"
SLURM_MEM="16G"
SLURM_PARTITION=""
SLURM_ACCOUNT=""




log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}


require_dir() {
    if [[ ! -d "$1" ]]; then
        echo "Missing directory: $1" >&2
        exit 1
    fi
}


require_file() {
    if [[ ! -f "$1" ]]; then
        echo "Missing file: $1" >&2
        exit 1
    fi
}


is_true() {
    case "${1}" in
        1|true|TRUE|True|yes|YES|Yes|y|Y|on|ON|On)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}


is_false() {
    case "${1}" in
        0|false|FALSE|False|no|NO|No|n|N|off|OFF|Off)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}


is_auto() {
    case "${1}" in
        auto|AUTO|Auto)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}


validate_hx_plot_mode() {
    if is_true "${HX_PLOT_ENABLED}" || is_false "${HX_PLOT_ENABLED}" || is_auto "${HX_PLOT_ENABLED}"; then
        return 0
    fi

    echo "Unsupported HX_PLOT_ENABLED: ${HX_PLOT_ENABLED}. Use auto, true, or false." >&2
    exit 1
}


WORK_DIR=$(cd "${WORK_DIR}" && pwd)
RESULT_DIR="${WORK_DIR}/results"
mkdir -p "${RESULT_DIR}"


workflow_env=(
    "pdb_id=${PDB_ID}"
    "sim_id=${SIM_ID}"
    "n_rep=${N_REP}"
    "start_frame=${START_FRAME}"
    "skip_experiment_data=${SKIP_EXPERIMENT_DATA}"
    "HXMS_method=${HXMS_METHOD}"
    "protein_state=${PROTEIN_STATE}"
    "exp_data_file=${EXP_DATA_FILE}"
)


resolve_work_path() {
    local configured_path="$1"
    local default_relative_path="$2"

    if [[ -n "${configured_path}" ]]; then
        if [[ "${configured_path}" = /* ]]; then
            printf '%s\n' "${configured_path}"
        else
            printf '%s\n' "${WORK_DIR}/${configured_path}"
        fi
        return
    fi

    printf '%s\n' "${WORK_DIR}/${default_relative_path}"
}


hx_plot_dfout_path() {
    resolve_work_path "${HX_PLOT_DFOUT_FILE}" "${HX_PLOT_PREFIX}-dfout-${HX_PLOT_STATE}.csv"
}


hx_plot_fitdata_path() {
    resolve_work_path "${HX_PLOT_FITDATA_FILE}" "${HX_PLOT_PREFIX}-fitdata-${HX_PLOT_STATE}"
}


hx_plot_dg_path() {
    resolve_work_path "${HX_PLOT_DG_FILE}" "${HX_PLOT_PREFIX}-dg.csv"
}


hx_plot_resid_path() {
    resolve_work_path "${HX_PLOT_RESID_FILE}" "results/${PDB_ID}.resid"
}


hx_plot_output_dir_path() {
    resolve_work_path "${HX_PLOT_OUTPUT_DIR}" "."
}


submit_to_slurm() {
    require_dir "${WORK_DIR}"
    command -v sbatch >/dev/null 2>&1 || {
        echo "RUNNER=slurm requested, but sbatch was not found in PATH." >&2
        exit 1
    }

    local sbatch_args=(
        --job-name "${SLURM_JOB_NAME}"
        --time "${SLURM_TIME}"
        --cpus-per-task "${SLURM_CPUS_PER_TASK}"
        --mem "${SLURM_MEM}"
        --chdir "${WORK_DIR}"
    )

    if [[ -n "${SLURM_PARTITION}" ]]; then
        sbatch_args+=(--partition "${SLURM_PARTITION}")
    fi
    if [[ -n "${SLURM_ACCOUNT}" ]]; then
        sbatch_args+=(--account "${SLURM_ACCOUNT}")
    fi

    log "Submitting analysis workflow to Slurm from ${WORK_DIR}"
    env \
        RUNNER="${RUNNER}" \
        ANALYSIS_DRIVER_INNER=1 \
        WORK_DIR="${WORK_DIR}" \
        PDB_ID="${PDB_ID}" \
        SIM_ID="${SIM_ID}" \
        N_REP="${N_REP}" \
        START_FRAME="${START_FRAME}" \
        SKIP_EXPERIMENT_DATA="${SKIP_EXPERIMENT_DATA}" \
        HXMS_METHOD="${HXMS_METHOD}" \
        PROTEIN_STATE="${PROTEIN_STATE}" \
        EXP_DATA_FILE="${EXP_DATA_FILE}" \
        SLURM_JOB_NAME="${SLURM_JOB_NAME}" \
        SLURM_TIME="${SLURM_TIME}" \
        SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK}" \
        SLURM_MEM="${SLURM_MEM}" \
        SLURM_PARTITION="${SLURM_PARTITION}" \
        SLURM_ACCOUNT="${SLURM_ACCOUNT}" \
        sbatch "${sbatch_args[@]}" "${SELF_PATH}"
}


activate_runtime() {
    local venv_activate

    require_dir "${PROJECT_ROOT}"
    PROJECT_ROOT=$(cd "${PROJECT_ROOT}" && pwd)
    require_file "${PROJECT_ROOT}/source.sh"
    venv_activate="${PROJECT_ROOT}/.venv/bin/activate"

    set +u
    source "${PROJECT_ROOT}/source.sh"

    if [[ -n "${PYTHON_ENV_CMD}" ]]; then
        log "Activating Python runtime with PYTHON_ENV_CMD"
        eval "${PYTHON_ENV_CMD}"
    elif [[ -f "${venv_activate}" ]]; then
        log "Activating Python runtime from ${venv_activate}"
        source "${venv_activate}"
    else
        log "No ${venv_activate} found. Using the current python from PATH."
    fi

    set -u
}


verify_python_environment() {
    local module_list="numpy pandas scipy matplotlib pymbar mdtraj_upside"
    local missing_modules

    if ! is_true "${SKIP_EXPERIMENT_DATA}"; then
        module_list="${module_list} mdtraj"
    fi
    command -v python >/dev/null 2>&1 || {
        echo "Python was not found after runtime activation. Set PYTHON_ENV_CMD or create ${PROJECT_ROOT}/.venv." >&2
        exit 1
    }

    missing_modules=$(ANALYSIS_REQUIRED_MODULES="${module_list}" python -c 'import importlib.util, os; modules = os.environ["ANALYSIS_REQUIRED_MODULES"].split(); missing = [name for name in modules if importlib.util.find_spec(name) is None]; print(",".join(missing))')
    if [[ -n "${missing_modules}" ]]; then
        echo "Python environment is missing required modules: ${missing_modules}. Set PYTHON_ENV_CMD to load a configured environment or install the dependencies into ${PROJECT_ROOT}/.venv." >&2
        exit 1
    fi

    log "Using python executable: $(command -v python)"
}


run_python_step() {
    local label="$1"
    local mode="$2"
    local script_path="$3"
    shift 3

    log "Running ${label}"
    env "${workflow_env[@]}" "analysis_mode=${mode}" python "${script_path}" "$@"
}


run_config_step() {
    local action="$1"
    log "Running 1.config.py (${action})"
    env "${workflow_env[@]}" "workflow_action=${action}" python "${SCRIPT_DIR}/1.config.py"
}


run_shell_step() {
    local label="$1"
    local script_path="$2"
    log "Running ${label}"
    env "${workflow_env[@]}" bash "${script_path}"
}


run_hx_plot_step() {
    local dfout_path fitdata_path dg_path resid_path output_dir
    local missing_inputs=()

    if is_false "${HX_PLOT_ENABLED}"; then
        log "Skipping 6.generate_hx_plots.py because HX_PLOT_ENABLED=${HX_PLOT_ENABLED}"
        return
    fi

    dfout_path=$(hx_plot_dfout_path)
    fitdata_path=$(hx_plot_fitdata_path)
    dg_path=$(hx_plot_dg_path)
    resid_path=$(hx_plot_resid_path)
    output_dir=$(hx_plot_output_dir_path)

    [[ -f "${dfout_path}" ]] || missing_inputs+=("${dfout_path}")
    [[ -f "${fitdata_path}" ]] || missing_inputs+=("${fitdata_path}")
    [[ -f "${dg_path}" ]] || missing_inputs+=("${dg_path}")
    [[ -f "${resid_path}" ]] || missing_inputs+=("${resid_path}")

    if (( ${#missing_inputs[@]} > 0 )); then
        if is_true "${HX_PLOT_ENABLED}"; then
            printf 'HX plot step requested but required inputs are missing:\n' >&2
            printf '  %s\n' "${missing_inputs[@]}" >&2
            exit 1
        fi

        log "Skipping 6.generate_hx_plots.py because required HX plot inputs are unavailable"
        return
    fi

    log "Running 6.generate_hx_plots.py"
    env \
        "${workflow_env[@]}" \
        "hx_plot_work_dir=${WORK_DIR}" \
        "hx_plot_results_dir=${RESULT_DIR}" \
        "hx_plot_output_dir=${output_dir}" \
        "hx_plot_prefix=${HX_PLOT_PREFIX}" \
        "hx_plot_state=${HX_PLOT_STATE}" \
        "hx_plot_dfout_path=${dfout_path}" \
        "hx_plot_fitdata_path=${fitdata_path}" \
        "hx_plot_dg_path=${dg_path}" \
        "hx_plot_resid_path=${resid_path}" \
        python "${SCRIPT_DIR}/6.generate_hx_plots.py"
}


run_workflow() {
    require_dir "${WORK_DIR}"
    require_dir "${WORK_DIR}/inputs"
    require_dir "${WORK_DIR}/pdb"

    local current_mplconfigdir
    local current_xdg_cache

    current_mplconfigdir=$(printenv MPLCONFIGDIR 2>/dev/null || true)
    current_xdg_cache=$(printenv XDG_CACHE_HOME 2>/dev/null || true)

    if [[ -n "${current_mplconfigdir}" ]]; then
        export MPLCONFIGDIR="${current_mplconfigdir}"
    else
        export MPLCONFIGDIR="${RESULT_DIR}/.mpl-cache"
    fi

    if [[ -n "${current_xdg_cache}" ]]; then
        export XDG_CACHE_HOME="${current_xdg_cache}"
    else
        export XDG_CACHE_HOME="${RESULT_DIR}/.cache"
    fi

    mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"

    activate_runtime
    validate_hx_plot_mode
    verify_python_environment
    cd "${WORK_DIR}"

    if is_true "${SKIP_EXPERIMENT_DATA}"; then
        log "Skipping 0.run_HXMS.py because SKIP_EXPERIMENT_DATA=${SKIP_EXPERIMENT_DATA}"
    else
        log "Running 0.run_HXMS.py"
        env "${workflow_env[@]}" python "${SCRIPT_DIR}/0.run_HXMS.py"
    fi

    run_config_step "config"
    run_shell_step "2.traj_ana.sh" "${SCRIPT_DIR}/2.traj_ana.sh"
    run_shell_step "3.get_protaction_states.sh" "${SCRIPT_DIR}/3.get_protaction_states.sh"
    run_python_step "4.calc_D_uptake.py (uptake)" "uptake" "${SCRIPT_DIR}/4.calc_D_uptake.py"
    run_python_step "4.calc_D_uptake.py (stability)" "stability" "${SCRIPT_DIR}/4.calc_D_uptake.py"
    if is_true "${SKIP_EXPERIMENT_DATA}"; then
        log "Skipping 5.analyze_D_uptake.py (uptake) because SKIP_EXPERIMENT_DATA=${SKIP_EXPERIMENT_DATA}"
    else
        run_python_step "5.analyze_D_uptake.py (uptake)" "uptake" "${SCRIPT_DIR}/5.analyze_D_uptake.py"
    fi
    run_python_step "5.analyze_D_uptake.py (dg_summary)" "dg_summary" "${SCRIPT_DIR}/5.analyze_D_uptake.py"
    run_hx_plot_step
}


case "${RUNNER}" in
    local)
        run_workflow
        ;;
    slurm)
        if printenv ANALYSIS_DRIVER_INNER >/dev/null 2>&1; then
            run_workflow
        else
            submit_to_slurm
        fi
        ;;
    *)
        echo "Unsupported RUNNER: ${RUNNER}. Use local or slurm." >&2
        exit 1
        ;;
esac
