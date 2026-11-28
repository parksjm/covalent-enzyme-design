#!/usr/bin/env bash
set -euo pipefail
: "${GERANIUM_BOUNDARY_RELEASE_APPROVED:?Set only after PI shared-host resource approval and final QC}"
test "$(hostname -s | tr "[:upper:]" "[:lower:]")" = geranium
cd /data/jpc/covalent-enzyme-design
export PYTHONPATH=/data/jpc/covalent-enzyme-design/src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES=
job_root=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/activation_support_boundary_20261128/job
mkdir "$job_root"
pids=()
/data/jpc/code/ambertools26/venv/bin/python /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/reactant_coverage_plan_20261128/preparation/run_boundary.py --prepared /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/activation_support_boundary_20261128/prepared.json --label activation_boundary_01_original_es --sha256 ea35a86d025ef106c7397fad832a158e9bb5e77740dfad69f39d88e53dad6d4d > "$job_root/activation_boundary_01_original_es.log" 2>&1 &
pids+=("$!")
/data/jpc/code/ambertools26/venv/bin/python /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/reactant_coverage_plan_20261128/preparation/run_boundary.py --prepared /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/activation_support_boundary_20261128/prepared.json --label activation_boundary_01_dominant_es --sha256 ea35a86d025ef106c7397fad832a158e9bb5e77740dfad69f39d88e53dad6d4d > "$job_root/activation_boundary_01_dominant_es.log" 2>&1 &
pids+=("$!")
/data/jpc/code/ambertools26/venv/bin/python /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/reactant_coverage_plan_20261128/preparation/run_boundary.py --prepared /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/activation_support_boundary_20261128/prepared.json --label activation_boundary_02_original_es --sha256 ea35a86d025ef106c7397fad832a158e9bb5e77740dfad69f39d88e53dad6d4d > "$job_root/activation_boundary_02_original_es.log" 2>&1 &
pids+=("$!")
/data/jpc/code/ambertools26/venv/bin/python /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/reactant_coverage_plan_20261128/preparation/run_boundary.py --prepared /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/activation_support_boundary_20261128/prepared.json --label activation_boundary_02_dominant_es --sha256 ea35a86d025ef106c7397fad832a158e9bb5e77740dfad69f39d88e53dad6d4d > "$job_root/activation_boundary_02_dominant_es.log" 2>&1 &
pids+=("$!")
/data/jpc/code/ambertools26/venv/bin/python /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/reactant_coverage_plan_20261128/preparation/run_boundary.py --prepared /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/activation_support_boundary_20261128/prepared.json --label activation_boundary_03_original_es --sha256 ea35a86d025ef106c7397fad832a158e9bb5e77740dfad69f39d88e53dad6d4d > "$job_root/activation_boundary_03_original_es.log" 2>&1 &
pids+=("$!")
/data/jpc/code/ambertools26/venv/bin/python /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/reactant_coverage_plan_20261128/preparation/run_boundary.py --prepared /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/activation_support_boundary_20261128/prepared.json --label activation_boundary_03_dominant_es --sha256 ea35a86d025ef106c7397fad832a158e9bb5e77740dfad69f39d88e53dad6d4d > "$job_root/activation_boundary_03_dominant_es.log" 2>&1 &
pids+=("$!")
result=0
for pid in "${pids[@]}"; do
  if wait "$pid"; then :; else result=1; fi
done
exit "$result"
