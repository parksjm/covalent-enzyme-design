#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 VALIDATION_DIRECTORY" >&2
    exit 2
fi
validation_dir=$(cd "$1" && pwd)
set +u
source /data/jpc/code/ambertools26/venv/amber.sh
set -u
export PATH="/data/jpc/code/pmemd26/bin:$PATH"
cd "$validation_dir"

for path in md_heat.rst7 md_density.rst7 md_unrestrained.rst7; do
    if [[ -e "$path" ]]; then
        echo "Refusing to overwrite existing MD artifact: $path" >&2
        exit 2
    fi
done

CUDA_VISIBLE_DEVICES=0 pmemd.cuda -O \
    -i md_heat.mdin -o md_heat.mdout \
    -p solvated.prmtop -c solvated_xmin_finish.rst7 -ref solvated_xmin_finish.rst7 \
    -r md_heat.rst7 -x md_heat.nc -inf md_heat.mdinfo

CUDA_VISIBLE_DEVICES=0 pmemd.cuda -O \
    -i md_density.mdin -o md_density.mdout \
    -p solvated.prmtop -c md_heat.rst7 -ref solvated_xmin_finish.rst7 \
    -r md_density.rst7 -x md_density.nc -inf md_density.mdinfo

CUDA_VISIBLE_DEVICES=0 pmemd.cuda -O \
    -i md_unrestrained.mdin -o md_unrestrained.mdout \
    -p solvated.prmtop -c md_density.rst7 \
    -r md_unrestrained.rst7 -x md_unrestrained.nc -inf md_unrestrained.mdinfo

cpptraj -p solvated.prmtop -i short_md_strip.cpptraj > cpptraj_short_md_strip.log 2>&1
