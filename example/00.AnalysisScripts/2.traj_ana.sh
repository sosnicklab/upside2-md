#!/bin/bash
set -eo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../../source.sh"
set -u

pdb_id=${pdb_id:-glpG-RKRK-79HIS} # checkme
sim_id=${sim_id:-memb_test} # checkme
n_rep=${n_rep:-48} # checkme
trajectory_mode=${trajectory_mode:-upside}
hybrid_source_dir=${hybrid_source_dir:-}
if [ -z "${hybrid_source_pattern:-}" ]; then
    hybrid_source_pattern="$pdb_id.stage_7.{replica}.up"
fi
hybrid_project_overwrite=${hybrid_project_overwrite:-true}
replica_token='{replica}'

work_dir=./
input_dir=$work_dir/inputs
output_dir=$work_dir/outputs
run_dir=$output_dir/$sim_id

mkdir -p results "$run_dir"

for i in $(seq 0 $((n_rep-1)))
do
    traj=$output_dir/$sim_id/$pdb_id.run.$i.up
    if [ "$trajectory_mode" = martini_hybrid ]; then
        if [ -z "$hybrid_source_dir" ]; then
            echo "trajectory_mode=martini_hybrid requires hybrid_source_dir" >&2
            exit 1
        fi
        source_name=${hybrid_source_pattern//$replica_token/$i}
        source_traj=$hybrid_source_dir/$source_name
        if [ ! -f "$source_traj" ]; then
            echo "Missing hybrid trajectory: $source_traj" >&2
            exit 1
        fi
        case $hybrid_project_overwrite in
            1|true|TRUE|True|yes|YES|Yes|y|Y|on|ON|On)
                python "$UPSIDE_HOME/py/martini_hdx_project.py" \
                    "$input_dir/$pdb_id-HDX.up" "$source_traj" "$traj" --overwrite ;;
            *)
                python "$UPSIDE_HOME/py/martini_hdx_project.py" \
                    "$input_dir/$pdb_id-HDX.up" "$source_traj" "$traj" ;;
        esac
    elif [ "$trajectory_mode" != upside ]; then
        echo "Unsupported trajectory_mode: $trajectory_mode" >&2
        exit 1
    fi

    python "$script_dir/helpers/get_info_from_upside_traj.py" "$traj" results/${pdb_id}_${sim_id}_$i
    python "$UPSIDE_HOME/py/extract_vtf.py" "$traj" results/${pdb_id}_${sim_id}_$i.vtf
done
