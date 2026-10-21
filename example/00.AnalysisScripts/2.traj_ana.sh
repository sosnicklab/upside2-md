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
# {replica} numbers the files of a REMD ladder; {segment} numbers the consecutive continuations of one run
# (a production stage_7.0, stage_7.1, ...), which are joined in order into that replica's trajectory.
if [ -z "${hybrid_source_pattern:-}" ]; then
    hybrid_source_pattern="$pdb_id.stage_7.{segment}.up"
fi
hybrid_project_overwrite=${hybrid_project_overwrite:-true}
replica_token='{replica}'
segment_token='{segment}'

work_dir=./
input_dir=$work_dir/inputs
output_dir=$work_dir/outputs
run_dir=$output_dir/$sim_id

mkdir -p results "$run_dir"

if [ "$trajectory_mode" = martini_hybrid ]; then
    if [ -z "$hybrid_source_dir" ]; then
        echo "trajectory_mode=martini_hybrid requires hybrid_source_dir" >&2
        exit 1
    fi
    if [[ $hybrid_source_pattern != *"$replica_token"* ]] && [ "$n_rep" -ne 1 ]; then
        echo "hybrid_source_pattern $hybrid_source_pattern names one run, not a ladder: n_rep must be 1" >&2
        exit 1
    fi
    case $hybrid_project_overwrite in
        1|true|TRUE|True|yes|YES|Yes|y|Y|on|ON|On)
            overwrite=--overwrite ;;
        *)
            overwrite= ;;
    esac
elif [ "$trajectory_mode" != upside ]; then
    echo "Unsupported trajectory_mode: $trajectory_mode" >&2
    exit 1
fi

# The segments of one replica, numbered 0, 1, ... with none missing.
replica_segments() {
    local replica_pattern=$1
    local prefix suffix path name middle n_numbered
    segments=()
    if [[ $replica_pattern != *"$segment_token"* ]]; then
        segments=("$hybrid_source_dir/$replica_pattern")
        return
    fi
    while [ -f "$hybrid_source_dir/${replica_pattern//$segment_token/${#segments[@]}}" ]; do
        segments+=("$hybrid_source_dir/${replica_pattern//$segment_token/${#segments[@]}}")
    done
    prefix=${replica_pattern%%"$segment_token"*}
    suffix=${replica_pattern#*"$segment_token"}
    n_numbered=0
    for path in "$hybrid_source_dir"/*; do
        name=$(basename "$path")
        case $name in
            "$prefix"*"$suffix")
                middle=${name#"$prefix"}
                middle=${middle%"$suffix"}
                if [[ $middle =~ ^[0-9]+$ ]]; then
                    n_numbered=$((n_numbered+1))
                fi ;;
        esac
    done
    if [ ${#segments[@]} -eq 0 ] || [ $n_numbered -ne ${#segments[@]} ]; then
        echo "segments of $hybrid_source_dir/$replica_pattern are not numbered 0..N-1 without a gap" >&2
        exit 1
    fi
}

for i in $(seq 0 $((n_rep-1)))
do
    traj=$output_dir/$sim_id/$pdb_id.run.$i.up
    if [ "$trajectory_mode" = martini_hybrid ]; then
        replica_segments "${hybrid_source_pattern//$replica_token/$i}"
        for segment in "${segments[@]}"; do
            if [ ! -f "$segment" ]; then
                echo "Missing hybrid trajectory: $segment" >&2
                exit 1
            fi
        done
        python "$UPSIDE_HOME/py/martini_hdx_project.py" "$input_dir/$pdb_id-HDX.up" "$traj" "${segments[@]}" \
            --pdb-numbering "$input_dir/$pdb_id.chi" $overwrite
    fi

    python "$UPSIDE_HOME/py/get_info_from_upside_traj.py" "$traj" results/${pdb_id}_${sim_id}_$i
    python "$UPSIDE_HOME/py/extract_vtf.py" "$traj" results/${pdb_id}_${sim_id}_$i.vtf
done
