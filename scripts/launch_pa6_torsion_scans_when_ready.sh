#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 REPOSITORY_ROOT GEOMETRY_DIR SCAN_DIR" >&2
    exit 2
fi

repository_root=$1
geometry_dir=$2
scan_dir=$3
uv_bin=/home/jpc/.local/bin/uv
scratch_root=/data/scratch/gaussian_scratch/jpc

optimization_sessions=(
    pa6_torsion_internal_neutral_opt_20261028
    pa6_torsion_n_terminus_cation_opt_20261028
    pa6_torsion_c_terminus_anion_opt_20261028
)

while true; do
    active=0
    for session in "${optimization_sessions[@]}"; do
        if tmux has-session -t "$session" 2>/dev/null; then
            active=1
        fi
    done
    if [[ $active -eq 0 ]]; then
        break
    fi
    sleep 30
done

for fragment in internal_neutral n_terminus_cation c_terminus_anion; do
    status_path="${geometry_dir}/${fragment}/job_control/status.json"
    if ! jq -e '.exit_status == 0' "$status_path" >/dev/null; then
        echo "Fragment optimization did not exit successfully: ${fragment}" >&2
        exit 1
    fi
done

cd -- "$repository_root"
"$uv_bin" run python scripts/prepare_pa6_torsion_scans.py \
    --geometry-dir "$geometry_dir" \
    --output-dir "$scan_dir"

manifest_path="${scan_dir}/manifest.json"
while IFS=$'\t' read -r fragment label; do
    pair_dir="${scan_dir}/${fragment}/${label}"
    forward_input="${pair_dir}/forward/scan.gjf"
    reverse_input="${pair_dir}/reverse/scan.gjf"
    session="pa6_scan_${fragment}_${label}_20261028"
    "$uv_bin" run python scripts/run_nyl12_job.py launch \
        --session "$session" \
        --output-dir "${pair_dir}/job_control" \
        --cwd "$pair_dir" \
        --input "$forward_input" \
        --input "$reverse_input" \
        --input "$manifest_path" \
        --input "${geometry_dir}/${fragment}/optimization.log" \
        --input "${repository_root}/scripts/run_gaussian16_file.sh" \
        --input "${repository_root}/scripts/run_gaussian16_bidirectional_scan.sh" \
        --input "${repository_root}/scripts/prepare_pa6_torsion_scans.py" \
        --input "${repository_root}/src/nyl12_covalent/pa6_torsion_scans.py" \
        -- /usr/bin/sg g16 -c \
        "/bin/bash ${repository_root}/scripts/run_gaussian16_bidirectional_scan.sh ${forward_input} ${pair_dir}/forward/scan.log ${scratch_root}/${session}_forward ${reverse_input} ${pair_dir}/reverse/scan.log ${scratch_root}/${session}_reverse"
done < <(jq -r '.scans[] | [.fragment, .label] | @tsv' "$manifest_path")
