#!/bin/bash
set -euo pipefail
cd /data/jpc/covalent-enzyme-design
trial_root=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/04_neb/spacing_tests
for trial_name in 01_uniform_24_lbfgs 02_uniform_48_lbfgs 03_adaptive_24_lbfgs 04_uniform_24_cg; do
    /data/jpc/code/ambertools26/venv/bin/python scripts/run_nyl12_job.py launch \
        --session "nyl12_neb_${trial_name}_view_20261119" \
        --output-dir "$trial_root/$trial_name/progress_job" \
        --cwd /data/jpc/covalent-enzyme-design \
        --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
        --env PYTHONPATH=/data/jpc/covalent-enzyme-design:/data/jpc/covalent-enzyme-design/src \
        --input scripts/watch_neb_spacing_trial.py \
        --input src/nyl12_covalent/relay_neb_progress.py \
        --input src/nyl12_covalent/neb_xyz_snapshot.py \
        --input src/nyl12_covalent/neb_snapshot.py \
        --input src/nyl12_covalent/neb_progress.py \
        --input "$trial_root/$trial_name/config.json" \
        -- /data/jpc/code/enzyme-substrate-modeling/.pixi/envs/analysis/bin/python \
        scripts/watch_neb_spacing_trial.py --config "$trial_root/$trial_name/config.json" \
        --watch --interval 60 --hours 72
done
