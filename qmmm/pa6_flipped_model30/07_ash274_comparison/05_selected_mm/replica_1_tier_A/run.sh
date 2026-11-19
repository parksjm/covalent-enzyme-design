#!/usr/bin/env bash
set -euo pipefail
case_dir="$(cd "$(dirname "$0")" && pwd)"
export PYTHONPATH="/data/jpc/covalent-enzyme-design/src"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
export AMBERHOME="/data/jpc/code/ambertools26/venv"
export LD_LIBRARY_PATH="/data/jpc/code/ambertools26/venv/lib"
exec /data/jpc/code/enzyme-substrate-modeling/.pixi/envs/analysis/bin/python \
  "$case_dir/prepare_ash274_mm.py"
