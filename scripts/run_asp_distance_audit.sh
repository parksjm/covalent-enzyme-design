#!/usr/bin/env bash
# Reproduce the four-site Asp272/Asp274 distance audit on original full MD.
set -euo pipefail
cd /data/jpc/covalent-enzyme-design
root=/data/jpc/covalent-enzyme-design
model="$root/qmmm/pa6_flipped_model30"
args=(--topology "$model/00_system/md_solvated.prmtop"
      --restart "initial=$model/00_system/md_solvated.inpcrd"
      --restart "minimized=$model/01_mm_minimization/min.rst7"
      --output "$root/analysis/negoro_2023_protonation/nyl12_distances")
for replica in 01 02; do
    for segment in 05_production_10ns 06_production_90ns; do
        args+=(--trajectory "neutral${replica}_${segment}=$model/03_md/neutral_replica_${replica}/${segment}.nc")
    done
done
exec /data/jpc/enzyme-substrate-modeling/.pixi/envs/analysis/bin/python \
    scripts/audit_nyl12_asp_pairs.py "${args[@]}"
