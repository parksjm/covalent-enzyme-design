#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 INPUT_ROOT TOPOLOGY GAUSSIAN_SCRATCH_ROOT" >&2
    exit 2
fi

input_root=$1
topology=$2
scratch_root=$3
project_root=${NYL12_PROJECT_ROOT:-/data/jpc/covalent-enzyme-design}
cli=${NYL12_BENCHMARK_CLI:-$project_root/scripts/enzyme_path_method_benchmark.py}
sander=/data/jpc/code/ambertools26/venv/bin/sander
python=/data/jpc/code/ambertools26/venv/bin/python
xtb=/data/jpc/code/enzyme-substrate-modeling/.pixi/envs/analysis/bin/xtb
g16=/usr/local/g16/g16
g16_profile=/usr/local/g16/bsd/g16.profile
worker_count=${BENCHMARK_WORKERS:-3}

for required in \
    "$input_root/benchmark_manifest.json" "$topology" "$cli" \
    "$sander" "$python" "$xtb" "$g16" "$g16_profile"; do
    [[ -e "$required" ]] || { echo "missing prerequisite: $required" >&2; exit 2; }
done
[[ -d "$scratch_root" && -w "$scratch_root" ]] || {
    echo "Gaussian scratch is not writable: $scratch_root" >&2
    exit 2
}
[[ "$worker_count" =~ ^[1-9][0-9]*$ ]] || {
    echo "BENCHMARK_WORKERS must be a positive integer" >&2
    exit 2
}

expected_topology_sha256=$(
    "$python" -c \
        'import json,sys; print(json.load(open(sys.argv[1]))["topology_sha256"])' \
        "$input_root/benchmark_manifest.json"
)
actual_topology_sha256=$(sha256sum "$topology" | awk '{print $1}')
[[ "$actual_topology_sha256" == "$expected_topology_sha256" ]] || {
    echo "topology hash mismatch" >&2
    exit 2
}

export g16root=/usr/local
export AMBERHOME=/data/jpc/code/ambertools26/venv
export LD_LIBRARY_PATH="$AMBERHOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
set +u
source "$g16_profile"
set -u
export input_root topology scratch_root sander python xtb cli

run_window() {
    local window=$1
    local window_dir=$input_root/windows/$window
    local dftb_dir=$window_dir/dftb3_mm
    local dft_dir=$window_dir/wb97xd_mm
    local xtb_dir=$window_dir/gfn2_xtb
    local scratch=$scratch_root/$window
    local restart=$window_dir/input.ncrst
    local output

    for output in \
        "$dftb_dir/single_point.mdout" \
        "$dft_dir/single_point.mdout" \
        "$xtb_dir/gfn2.out"; do
        [[ ! -e "$output" ]] || { echo "refusing to overwrite: $output" >&2; return 2; }
    done
    mkdir -p "$scratch"
    chmod 700 "$scratch"
    export GAUSS_SCRDIR=$scratch

    cd "$dftb_dir"
    "$sander" -O -i minimize.mdin -p "$topology" -c "$restart" \
        -o single_point.mdout -r single_point.ncrst -x single_point.nc \
        -inf single_point.mdinfo > launcher.log 2>&1
    grep -q 'DFTBESCF=' single_point.mdout
    grep -q 'Energy calculation finished' single_point.mdout

    cd "$dft_dir"
    "$sander" -O -i minimize.mdin -p "$topology" -c "$restart" \
        -o single_point.mdout -r single_point.ncrst -x single_point.nc \
        -inf single_point.mdinfo > launcher.log 2>&1
    grep -q 'EXTERNESCF=' single_point.mdout
    grep -q 'Energy calculation finished' single_point.mdout
    grep -q 'Normal termination of Gaussian' gau_job.log.*
    [[ -s old.gau_job.inp ]]

    "$python" "$cli" convert \
        --gaussian-input "$dft_dir/old.gau_job.inp" \
        --output "$xtb_dir"
    cd "$xtb_dir"
    OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 "$xtb" qm_link.xyz --gfn 2 \
        --chrg -2 --uhf 0 --input xtb.inp --namespace gfn2 \
        > gfn2.out 2>&1
    grep -q 'normal termination of xtb' gfn2.out
}
export -f run_window

mapfile -t windows < <(
    "$python" -c \
        'import json,sys; print("\n".join(json.load(open(sys.argv[1]))["windows"]))' \
        "$input_root/benchmark_manifest.json"
)
printf '%s\0' "${windows[@]}" | \
    xargs -0 -n 1 -P "$worker_count" bash -c 'run_window "$1"' _

"$python" "$cli" analyze \
    --input-root "$input_root" \
    --output "$input_root/analysis/result"

"$python" - "$input_root" <<'PY'
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
files = sorted(
    path for path in root.glob("**/*")
    if path.is_file() and path.name not in {"process_acceptance.json"}
)
required = [
    root / "analysis/result/matched_energies.json",
    root / "analysis/result/method_errors.json",
    root / "analysis/result/embedded_method_comparison.png",
]
if not all(path.is_file() and path.stat().st_size > 0 for path in required):
    raise RuntimeError("analysis acceptance failed")
record = {
    "process_complete": True,
    "scientific_acceptance_pending_interpretation": True,
    "sha256": {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in files
    },
}
(root / "process_acceptance.json").write_text(json.dumps(record, indent=2) + "\n")
PY
