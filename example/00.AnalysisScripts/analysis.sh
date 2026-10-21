#!/bin/bash
set -eo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../../source.sh"
set -u

#----------------------------------------------------------------------
## User settings
#----------------------------------------------------------------------

self_path=$script_dir/$(basename $0)

# checkme: update these for the target run mode, workflow directory, and repository install.
runner=${runner:-local} # local or slurm
work_dir=${work_dir:-./}
python_env_cmd=${python_env_cmd:-} # Leave empty to use UPSIDE_HOME/.venv when present.

# checkme: update these for the target system and experimental condition.
pdb_id=${pdb_id:-glpG-RKRK-79HIS}
sim_id=${sim_id:-memb_test}
n_rep=${n_rep:-48}
start_frame=${start_frame:-100}

# Use trajectory_mode=martini_hybrid to project full hybrid trajectories into the
# ordinary protein-only HDX trajectory contract before running steps 2--6.
trajectory_mode=${trajectory_mode:-upside} # upside or martini_hybrid
hybrid_source_dir=${hybrid_source_dir:-}
if [ -z "${hybrid_source_pattern:-}" ]; then
    hybrid_source_pattern="$pdb_id.stage_7.{replica}.up"
fi
hybrid_project_overwrite=${hybrid_project_overwrite:-true}

# Optional calibrated water-accessibility arrays, shaped exactly like PS.npy.
# Leave the pattern empty to retain the stock protein-only protection state.
water_accessibility_dir=${water_accessibility_dir:-}
water_accessibility_pattern=${water_accessibility_pattern:-}

# checkme: set skip_experiment_data=true if no experimental data available, otherwise set all parameters for experimental data
skip_experiment_data=${skip_experiment_data:-false}
hxms_method=${hxms_method:-stretch_exp}
protein_state=${protein_state:-pd9}
exp_data_file=${exp_data_file:-GlpG psWT Sub final peptides up sum 11192024.csv}

# checkme: optional legacy HX plot step. Use auto to run only when the legacy inputs are present.
hx_plot_enabled=${hx_plot_enabled:-auto} # auto, true, or false
hx_plot_prefix=${hx_plot_prefix:-glpg}
hx_plot_state=${hx_plot_state:-$protein_state}
hx_plot_dfout_file=${hx_plot_dfout_file:-}
hx_plot_fitdata_file=${hx_plot_fitdata_file:-}
hx_plot_dg_file=${hx_plot_dg_file:-}
hx_plot_resid_file=${hx_plot_resid_file:-}
hx_plot_output_dir=${hx_plot_output_dir:-}

#----------------------------------------------------------------------
## Slurm settings
#----------------------------------------------------------------------

# checkme: adjust these only when runner=slurm.
slurm_job_name=${slurm_job_name:-analysis_$pdb_id}
slurm_time=${slurm_time:-08:00:00}
slurm_cpus_per_task=${slurm_cpus_per_task:-4}
slurm_mem=${slurm_mem:-16G}
slurm_partition=${slurm_partition:-}
slurm_account=${slurm_account:-}

work_dir=$(cd $work_dir && pwd)
result_dir=$work_dir/results
mkdir -p $result_dir

workflow_env=(
    "pdb_id=$pdb_id"
    "sim_id=$sim_id"
    "n_rep=$n_rep"
    "start_frame=$start_frame"
    "skip_experiment_data=$skip_experiment_data"
    "HXMS_method=$hxms_method"
    "protein_state=$protein_state"
    "exp_data_file=$exp_data_file"
    "trajectory_mode=$trajectory_mode"
    "hybrid_source_dir=$hybrid_source_dir"
    "hybrid_source_pattern=$hybrid_source_pattern"
    "hybrid_project_overwrite=$hybrid_project_overwrite"
    "water_accessibility_dir=$water_accessibility_dir"
    "water_accessibility_pattern=$water_accessibility_pattern"
)

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

require_dir() {
    if [ ! -d $1 ]; then
        echo "Missing directory: $1" >&2
        exit 1
    fi
}

require_file() {
    if [ ! -f $1 ]; then
        echo "Missing file: $1" >&2
        exit 1
    fi
}

is_true() {
    case $1 in
        1|true|TRUE|True|yes|YES|Yes|y|Y|on|ON|On)
            return 0 ;;
        *)
            return 1 ;;
    esac
}

is_false() {
    case $1 in
        0|false|FALSE|False|no|NO|No|n|N|off|OFF|Off)
            return 0 ;;
        *)
            return 1 ;;
    esac
}

is_auto() {
    case $1 in
        auto|AUTO|Auto)
            return 0 ;;
        *)
            return 1 ;;
    esac
}

validate_hx_plot_mode() {
    if is_true $hx_plot_enabled || is_false $hx_plot_enabled || is_auto $hx_plot_enabled; then
        return 0
    fi
    echo "Unsupported hx_plot_enabled: $hx_plot_enabled. Use auto, true, or false." >&2
    exit 1
}

resolve_work_path() {
    local configured_path=$1
    local default_relative_path=$2

    if [ -n "$configured_path" ]; then
        case $configured_path in
            /*)
                echo $configured_path ;;
            *)
                echo $work_dir/$configured_path ;;
        esac
        return
    fi

    echo $work_dir/$default_relative_path
}

hx_plot_dfout_path() {
    resolve_work_path "$hx_plot_dfout_file" "$hx_plot_prefix-dfout-$hx_plot_state.csv"
}

hx_plot_fitdata_path() {
    resolve_work_path "$hx_plot_fitdata_file" "$hx_plot_prefix-fitdata-$hx_plot_state"
}

hx_plot_dg_path() {
    resolve_work_path "$hx_plot_dg_file" "$hx_plot_prefix-dg.csv"
}

hx_plot_resid_path() {
    resolve_work_path "$hx_plot_resid_file" "results/$pdb_id.resid"
}

hx_plot_output_dir_path() {
    resolve_work_path "$hx_plot_output_dir" "."
}

submit_to_slurm() {
    require_dir $work_dir
    command -v sbatch >/dev/null 2>&1 || {
        echo "runner=slurm requested, but sbatch was not found in PATH." >&2
        exit 1
    }

    local sbatch_args=(
        --job-name $slurm_job_name
        --time $slurm_time
        --cpus-per-task $slurm_cpus_per_task
        --mem $slurm_mem
        --chdir $work_dir
    )

    if [ -n "$slurm_partition" ]; then
        sbatch_args+=(--partition $slurm_partition)
    fi
    if [ -n "$slurm_account" ]; then
        sbatch_args+=(--account $slurm_account)
    fi

    log "Submitting analysis workflow to Slurm from $work_dir"
    env \
        runner=$runner \
        ANALYSIS_DRIVER_INNER=1 \
        work_dir=$work_dir \
        pdb_id=$pdb_id \
        sim_id=$sim_id \
        n_rep=$n_rep \
        start_frame=$start_frame \
        skip_experiment_data=$skip_experiment_data \
        hxms_method=$hxms_method \
        protein_state=$protein_state \
        exp_data_file=$exp_data_file \
        trajectory_mode=$trajectory_mode \
        hybrid_source_dir=$hybrid_source_dir \
        hybrid_source_pattern=$hybrid_source_pattern \
        hybrid_project_overwrite=$hybrid_project_overwrite \
        water_accessibility_dir=$water_accessibility_dir \
        water_accessibility_pattern=$water_accessibility_pattern \
        slurm_job_name=$slurm_job_name \
        slurm_time=$slurm_time \
        slurm_cpus_per_task=$slurm_cpus_per_task \
        slurm_mem=$slurm_mem \
        slurm_partition=$slurm_partition \
        slurm_account=$slurm_account \
        sbatch ${sbatch_args[@]} $self_path
}

activate_runtime() {
    local venv_activate

    require_dir $UPSIDE_HOME
    venv_activate=$UPSIDE_HOME/.venv/bin/activate

    if [ -n "$python_env_cmd" ]; then
        log "Activating Python runtime with python_env_cmd"
        eval "$python_env_cmd"
    elif [ -f $venv_activate ]; then
        log "Activating Python runtime from $venv_activate"
        source $venv_activate
    else
        log "No $venv_activate found. Using the current python from PATH."
    fi
}

verify_python_environment() {
    local module_list="numpy pandas scipy matplotlib pymbar mdtraj_upside tables"
    local missing_modules

    if ! is_true $skip_experiment_data; then
        module_list="$module_list mdtraj"
    fi
    command -v python >/dev/null 2>&1 || {
        echo "Python was not found after runtime activation. Set python_env_cmd or create $UPSIDE_HOME/.venv." >&2
        exit 1
    }

    missing_modules=$(ANALYSIS_REQUIRED_MODULES="$module_list" python -c 'import importlib.util, os; modules = os.environ["ANALYSIS_REQUIRED_MODULES"].split(); missing = [name for name in modules if importlib.util.find_spec(name) is None]; print(",".join(missing))')
    if [ -n "$missing_modules" ]; then
        echo "Python environment is missing required modules: $missing_modules. Set python_env_cmd to load a configured environment or install the dependencies into $UPSIDE_HOME/.venv." >&2
        exit 1
    fi

    log "Using python executable: $(command -v python)"
}

run_python_step() {
    local label=$1
    local mode=$2
    local script_path=$3
    shift 3

    log "Running $label"
    env "${workflow_env[@]}" analysis_mode=$mode python $script_path "$@"
}

run_config_step() {
    local action=$1
    log "Running 1.config.py ($action)"
    env "${workflow_env[@]}" workflow_action=$action python $script_dir/1.config.py
}

run_shell_step() {
    local label=$1
    local script_path=$2
    log "Running $label"
    env "${workflow_env[@]}" bash $script_path
}

run_hx_plot_step() {
    local dfout_path fitdata_path dg_path resid_path output_dir
    local missing_inputs=()

    if is_false $hx_plot_enabled; then
        log "Skipping 6.generate_hx_plots.py because hx_plot_enabled=$hx_plot_enabled"
        return
    fi

    dfout_path=$(hx_plot_dfout_path)
    fitdata_path=$(hx_plot_fitdata_path)
    dg_path=$(hx_plot_dg_path)
    resid_path=$(hx_plot_resid_path)
    output_dir=$(hx_plot_output_dir_path)

    [ -f $dfout_path ] || missing_inputs+=($dfout_path)
    [ -f $fitdata_path ] || missing_inputs+=($fitdata_path)
    [ -f $dg_path ] || missing_inputs+=($dg_path)
    [ -f $resid_path ] || missing_inputs+=($resid_path)

    if [ ${#missing_inputs[@]} -gt 0 ]; then
        if is_true $hx_plot_enabled; then
            printf 'HX plot step requested but required inputs are missing:\n' >&2
            printf '  %s\n' ${missing_inputs[@]} >&2
            exit 1
        fi

        log "Skipping 6.generate_hx_plots.py because required HX plot inputs are unavailable"
        return
    fi

    log "Running 6.generate_hx_plots.py"
    env \
        "${workflow_env[@]}" \
        hx_plot_work_dir=$work_dir \
        hx_plot_results_dir=$result_dir \
        hx_plot_output_dir=$output_dir \
        hx_plot_prefix=$hx_plot_prefix \
        hx_plot_state=$hx_plot_state \
        hx_plot_dfout_path=$dfout_path \
        hx_plot_fitdata_path=$fitdata_path \
        hx_plot_dg_path=$dg_path \
        hx_plot_resid_path=$resid_path \
        python $script_dir/6.generate_hx_plots.py
}

run_workflow() {
    require_dir $work_dir
    require_dir $work_dir/inputs
    require_dir $work_dir/pdb

    local current_mplconfigdir
    local current_xdg_cache

    current_mplconfigdir=$(printenv MPLCONFIGDIR 2>/dev/null || true)
    current_xdg_cache=$(printenv XDG_CACHE_HOME 2>/dev/null || true)

    if [ -n "$current_mplconfigdir" ]; then
        export MPLCONFIGDIR=$current_mplconfigdir
    else
        export MPLCONFIGDIR=$result_dir/.mpl-cache
    fi

    if [ -n "$current_xdg_cache" ]; then
        export XDG_CACHE_HOME=$current_xdg_cache
    else
        export XDG_CACHE_HOME=$result_dir/.cache
    fi

    mkdir -p $MPLCONFIGDIR $XDG_CACHE_HOME

    activate_runtime
    validate_hx_plot_mode
    case $trajectory_mode in
        upside)
            ;;
        martini_hybrid)
            if [ -z "$hybrid_source_dir" ]; then
                echo "trajectory_mode=martini_hybrid requires hybrid_source_dir" >&2
                exit 1
            fi
            require_dir $hybrid_source_dir
            ;;
        *)
            echo "Unsupported trajectory_mode: $trajectory_mode" >&2
            exit 1
            ;;
    esac
    verify_python_environment
    cd $work_dir

    if is_true $skip_experiment_data; then
        log "Skipping 0.run_HXMS.py because skip_experiment_data=$skip_experiment_data"
    else
        log "Running 0.run_HXMS.py"
        env "${workflow_env[@]}" python $script_dir/0.run_HXMS.py
    fi

    run_config_step config
    run_shell_step "2.traj_ana.sh" $script_dir/2.traj_ana.sh
    run_shell_step "3.get_protaction_states.sh" $script_dir/3.get_protaction_states.sh
    run_python_step "4.calc_D_uptake.py (uptake)" uptake $script_dir/4.calc_D_uptake.py
    run_python_step "4.calc_D_uptake.py (stability)" stability $script_dir/4.calc_D_uptake.py
    if is_true $skip_experiment_data; then
        log "Skipping 5.analyze_D_uptake.py (uptake) because skip_experiment_data=$skip_experiment_data"
    else
        run_python_step "5.analyze_D_uptake.py (uptake)" uptake $script_dir/5.analyze_D_uptake.py
    fi
    run_python_step "5.analyze_D_uptake.py (dg_summary)" dg_summary $script_dir/5.analyze_D_uptake.py
    run_hx_plot_step
}

case $runner in
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
        echo "Unsupported runner: $runner. Use local or slurm." >&2
        exit 1
        ;;
esac
