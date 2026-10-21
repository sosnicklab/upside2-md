#!/bin/bash
set -eo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../../source.sh"
set -u

pdb_id=${pdb_id:-glpG-RKRK-79HIS} # checkme
sim_id=${sim_id:-memb_test} # checkme
n_rep=${n_rep:-48} # checkme
water_accessibility_dir=${water_accessibility_dir:-}
water_accessibility_pattern=${water_accessibility_pattern:-}
replica_token='{replica}'

work_dir=./
input_dir=$work_dir/inputs
output_dir=$work_dir/outputs
run_dir=$output_dir/$sim_id

mkdir -p results

for i in $(seq 0 $((n_rep-1)))
do
    traj=$output_dir/$sim_id/$pdb_id.run.$i.up
    protein_ps=results/${pdb_id}_${sim_id}_${i}_PS_protein.npy
    output_ps=results/${pdb_id}_${sim_id}_${i}_PS.npy
    python "$UPSIDE_HOME/py/get_protection_state.py" \
        "$input_dir/$pdb_id-HDX.up" "$traj" "$protein_ps" --residue results/$pdb_id.resid

    if [ -n "$water_accessibility_pattern" ]; then
        accessibility_name=${water_accessibility_pattern//$replica_token/$i}
        accessibility_path=$accessibility_name
        if [ -n "$water_accessibility_dir" ]; then
            accessibility_path=$water_accessibility_dir/$accessibility_name
        fi
        if [ ! -f "$accessibility_path" ]; then
            echo "Missing water-accessibility array: $accessibility_path" >&2
            exit 1
        fi
        python "$UPSIDE_HOME/py/combine_hdx_protection.py" \
            "$protein_ps" "$output_ps" --water-accessibility "$accessibility_path"
    else
        python "$UPSIDE_HOME/py/combine_hdx_protection.py" "$protein_ps" "$output_ps"
    fi
done
