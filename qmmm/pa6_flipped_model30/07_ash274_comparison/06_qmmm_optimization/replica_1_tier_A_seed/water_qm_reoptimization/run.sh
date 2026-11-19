#!/usr/bin/env bash
set -euo pipefail
case_dir="$(cd "$(dirname "$0")" && pwd)"
if [[ "${WATER_QM_LAUNCH_CONFIRMED:-}" != "1" ]]; then
  echo "Set WATER_QM_LAUNCH_CONFIRMED=1 only for the authorized single launch." >&2
  exit 2
fi
export PYTHONPATH="/data/jpc/covalent-enzyme-design/src"
export PYTHONDONTWRITEBYTECODE=1
export AMBERHOME="/data/jpc/code/ambertools26/venv"
export LD_LIBRARY_PATH="/data/jpc/code/ambertools26/venv/lib"
export OMP_NUM_THREADS=1
exec /data/jpc/code/enzyme-substrate-modeling/.pixi/envs/analysis/bin/python \
  "$case_dir/prepare_and_run_water_qm_qmmm.py"
