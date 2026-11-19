#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 INPUT_GJF OUTPUT_LOG SCRATCH_DIR" >&2
    exit 2
fi

input_gjf=$1
output_log=$2
scratch_dir=$3

if [[ ! -f "$input_gjf" ]]; then
    echo "Gaussian input does not exist: $input_gjf" >&2
    exit 2
fi
if [[ -e "$output_log" ]]; then
    echo "Refusing to overwrite Gaussian output: $output_log" >&2
    exit 2
fi

mkdir -p "$scratch_dir"
export g16root=/usr/local
export GAUSS_SCRDIR="$scratch_dir"
# Gaussian's vendor profile reads optional variables before assigning them.
# Restrict nounset relaxation to that sourced file, then restore fail-fast mode.
set +u
source /usr/local/g16/bsd/g16.profile
set -u

/usr/local/g16/g16 < "$input_gjf" > "$output_log" 2>&1
