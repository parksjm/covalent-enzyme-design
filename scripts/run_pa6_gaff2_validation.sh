#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 VALIDATION_DIRECTORY" >&2
    exit 2
fi
project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
validation_dir=$(cd "$1" && pwd)
common_input_dir="$project_root/parameters/pa6_l4_zwitterion_resp_gaff2/04_gaff2_validation"

set +u
source /data/jpc/code/ambertools26/venv/amber.sh
set -u
export PATH="/data/jpc/code/pmemd26/bin:$PATH"

cd "$validation_dir"
for path in isolated.prmtop solvated.prmtop gb_minimized.rst7 solvated_minimized.rst7; do
    if [[ -e "$path" ]]; then
        echo "Refusing to overwrite existing validation artifact: $path" >&2
        exit 2
    fi
done

tleap -f tleap.in > tleap.log 2>&1
if ! grep -q "Errors = 0" tleap.log; then
    echo "LEaP validation reported errors; inspect tleap.log" >&2
    exit 2
fi

sander -O \
    -i "$common_input_dir/gb_minimize.mdin" \
    -o gb_minimize.mdout \
    -p isolated.prmtop \
    -c isolated.inpcrd \
    -r gb_minimized.rst7 \
    -x gb_minimize.nc \
    -inf gb_minimize.mdinfo

CUDA_VISIBLE_DEVICES=0 pmemd.cuda -O \
    -i "$common_input_dir/solvated_minimize.mdin" \
    -o solvated_minimize.mdout \
    -p solvated.prmtop \
    -c solvated.inpcrd \
    -r solvated_minimized.rst7 \
    -x solvated_minimize.nc \
    -inf solvated_minimize.mdinfo

cpptraj -p solvated.prmtop -i "$common_input_dir/strip_minimization.cpptraj" > cpptraj_strip.log 2>&1
