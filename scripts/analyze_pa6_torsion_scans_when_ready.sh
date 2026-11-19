#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 REPOSITORY_ROOT SCAN_DIR ANALYSIS_DIR" >&2
    exit 2
fi

repository_root=$1
scan_dir=$2
analysis_dir=$3
manifest_path="${scan_dir}/manifest.json"
uv_bin=/home/jpc/.local/bin/uv

while true; do
    all_finished=1
    while IFS=$'\t' read -r fragment label; do
        status_path="${scan_dir}/${fragment}/${label}/job_control/status.json"
        exit_status=$(jq -r '.exit_status' "$status_path")
        if [[ "$exit_status" == "null" ]]; then
            all_finished=0
        elif [[ "$exit_status" != "0" ]]; then
            echo "Scan failed: ${fragment}/${label}, exit ${exit_status}" >&2
            exit 1
        fi
    done < <(jq -r '.scans[] | [.fragment, .label] | @tsv' "$manifest_path")
    if [[ $all_finished -eq 1 ]]; then
        break
    fi
    sleep 60
done

cd -- "$repository_root"
"$uv_bin" run python scripts/analyze_pa6_torsion_scans.py \
    --manifest "$manifest_path" \
    --output-dir "$analysis_dir"
