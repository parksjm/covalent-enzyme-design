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

for path in solvated_xmin_finish.rst7 solvated_xmin_finish.mdout; do
    if [[ -e "$path" ]]; then
        echo "Refusing to overwrite existing XMIN finish artifact: $path" >&2
        exit 2
    fi
done

sander -O \
    -i solvated_xmin_finish.mdin \
    -o solvated_xmin_finish.mdout \
    -p solvated.prmtop \
    -c solvated_xmin.rst7 \
    -r solvated_xmin_finish.rst7 \
    -x solvated_xmin_finish.nc \
    -inf solvated_xmin_finish.mdinfo

if grep -q "ERROR in" solvated_xmin_finish.mdout; then
    echo "XMIN reported a scientific error; inspect solvated_xmin_finish.mdout" >&2
    exit 2
fi
cpptraj -p solvated.prmtop -i strip_xmin_finish.cpptraj > cpptraj_xmin_finish.log 2>&1
