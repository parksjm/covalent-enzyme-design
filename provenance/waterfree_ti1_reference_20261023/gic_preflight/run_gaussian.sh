#!/usr/bin/env bash
set -eo pipefail
export g16root=/usr/local
source /usr/local/g16/bsd/g16.profile
export GAUSS_SCRDIR=/data/scratch/gaussian_scratch/jpc/waterfree39_gic_preflight01
mkdir -- "$GAUSS_SCRDIR"
cd -- /data/jpc/covalent-enzyme-design/output_nyl12_waterfree_reference_20261023/gic_preflight_01
/usr/local/g16/g16 < input.gjf > optimization.log 2>&1
