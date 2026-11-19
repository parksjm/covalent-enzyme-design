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
cd "$validation_dir"

for path in gb_xmin.rst7 solvated_xmin.rst7; do
    if [[ -e "$path" ]]; then
        echo "Refusing to overwrite existing XMIN artifact: $path" >&2
        exit 2
    fi
done

sander -O \
    -i gb_xmin_continue.mdin \
    -o gb_xmin.mdout \
    -p isolated.prmtop \
    -c gb_minimized.rst7 \
    -r gb_xmin.rst7 \
    -x gb_xmin.nc \
    -inf gb_xmin.mdinfo

sander -O \
    -i solvated_xmin_continue.mdin \
    -o solvated_xmin.mdout \
    -p solvated.prmtop \
    -c solvated_minimized.rst7 \
    -r solvated_xmin.rst7 \
    -x solvated_xmin.nc \
    -inf solvated_xmin.mdinfo

cpptraj -p solvated.prmtop -i strip_xmin.cpptraj > cpptraj_xmin_strip.log 2>&1
