#!/usr/bin/env bash
set -eo pipefail
export g16root=/usr/local
source /usr/local/g16/bsd/g16.profile
export GAUSS_SCRDIR=/data/scratch/gaussian_scratch/jpc/oxyanion63_svp_20261023_pilot01
mkdir -- "$GAUSS_SCRDIR"
cd -- /data/jpc/covalent-enzyme-design/output_nyl12_oxyanion_cluster_20261023/svp_pilot_01
engine_status=0
if /usr/bin/taskset -c 0-15 /usr/local/g16/g16 < input.gjf > optimization.log 2>&1; then
    engine_status=0
else
    engine_status=$?
fi
assessment_status=0
if PYTHONPATH=/data/jpc/covalent-enzyme-design/src /data/jpc/covalent-enzyme-design/.venv/bin/python /data/jpc/covalent-enzyme-design/output_nyl12_oxyanion_cluster_20261023/svp_pilot_01/assess_pilot.py; then
    assessment_status=0
else
    assessment_status=$?
fi
if [[ $engine_status -ne 0 ]]; then
    exit "$engine_status"
fi
exit "$assessment_status"
