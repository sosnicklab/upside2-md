#!/bin/bash
set -eo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../../source.sh"
set -u

pdb_id=${pdb_id:-glpG-RKRK-79HIS} # checkme
sim_id=${sim_id:-memb_test} # checkme
n_rep=${n_rep:-48} # checkme
trajectory_mode=${trajectory_mode:-upside}
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
    python "$UPSIDE_HOME/py/get_protection_state.py" \
        "$input_dir/$pdb_id-HDX.up" "$output_dir/$sim_id/$pdb_id.run.$i.up" \
        results/${pdb_id}_${sim_id}_${i}_PS_protein.npy --residue results/$pdb_id.resid
done

# A hybrid trajectory still has its lipid, so the membrane term of the protection state is measured from it:
# one contact radius from the tail-tail g(r) pooled over every replica, applied to all of them.
if [ -z "$water_accessibility_pattern" ] && [ "$trajectory_mode" = martini_hybrid ]; then
    views=()
    for i in $(seq 0 $((n_rep-1))); do
        views+=("$output_dir/$sim_id/$pdb_id.run.$i.up")
    done
    python "$UPSIDE_HOME/py/martini_hdx_membrane_accessibility.py" "$input_dir/$pdb_id-HDX.up" \
        "results/${pdb_id}_${sim_id}_${replica_token}_ACC.npy" "${views[@]}"
    water_accessibility_dir=results
    water_accessibility_pattern="${pdb_id}_${sim_id}_${replica_token}_ACC.npy"
fi

for i in $(seq 0 $((n_rep-1)))
do
    protein_ps=results/${pdb_id}_${sim_id}_${i}_PS_protein.npy
    output_ps=results/${pdb_id}_${sim_id}_${i}_PS.npy
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
