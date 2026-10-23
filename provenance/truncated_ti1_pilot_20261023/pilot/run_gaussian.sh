#!/usr/bin/env bash
set -eo pipefail
export g16root=/usr/local
source /usr/local/g16/bsd/g16.profile
export GAUSS_SCRDIR=/data/scratch/gaussian_scratch/jpc/truncated_ti1_20261023_pilot01
mkdir -- "$GAUSS_SCRDIR"
cd -- /data/jpc/covalent-enzyme-design/output_nyl12_truncated_reference_20261023/truncated_pilot_01
/usr/bin/taskset -c 0-15 /usr/local/g16/g16 < input.gjf > optimization.log 2>&1
PYTHONPATH=/data/jpc/covalent-enzyme-design/src /data/jpc/covalent-enzyme-design/.venv/bin/python /data/jpc/covalent-enzyme-design/output_nyl12_truncated_reference_20261023/truncated_pilot_01/assess_pilot.py
