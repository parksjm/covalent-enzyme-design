#!/usr/bin/env bash
set -euo pipefail

project_root=/data/jpc/covalent-enzyme-design
attempt_root="$project_root/qmmm/pa6_resp_gaff2_wt/06_replica03_native_thr_screen_20261106_attempt1"
md_root=/data/jpc/covalent-enzyme-design-pa6-md/models/pa6_modular_md_20261031/campaign_azalea_blackwell_attempt1
amber_python=/data/jpc/code/ambertools26/venv/bin/python
cpptraj=/data/jpc/code/ambertools26/venv/bin/cpptraj

export PYTHONPATH="$project_root/src"
export OMP_NUM_THREADS=1

"$amber_python" "$project_root/scripts/screen_modular_wt_thr_torsions.py" \
  --topology "$md_root/systems/wt/md_solvated.prmtop" \
  --trajectory "$md_root/systems/wt/replica_03/07_prod_100ns.nc" \
  --geometry-csv "$md_root/analysis/full_100ns_attempt1/wt_replica_03/geometry_timeseries.csv" \
  --cpptraj "$cpptraj" \
  --output "$attempt_root/torsion_screen/results"

"$amber_python" "$project_root/scripts/screen_modular_wt_native_thr.py" \
  --config "$project_root/configs/pa6_modular_wt_replica03_native_thr_screen_20261106.json" \
  --output "$attempt_root/native_screen/results" \
  --cpptraj "$cpptraj" \
  --top-k 10 \
  --verbose
