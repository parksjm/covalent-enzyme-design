# cwd: /data/jpc/covalent-enzyme-design; shell: bash; login=false; permission: escalated
# Numerical library threads limited to 1; GPU visibility empty; uv builds limited to 2.
git clone --no-hardlinks /tmp/nyl12_ipsae2_review_20261127 environments/ipsae2_20261128/source
git -C environments/ipsae2_20261128/source checkout fe560f4c2461eeb91aebded61ac8c1e83b264f60
uv add --project environments/ipsae2_20261128 --python 3.12 'ipsae2 @ git+https://github.com/DunbrackLab/IPSAE2.git@fe560f4c2461eeb91aebded61ac8c1e83b264f60'
uv add --project environments/ipsae2_20261128 --dev pytest
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=''
export MPLCONFIGDIR=/tmp/nyl12_ipsae2_qualification_mpl
environments/ipsae2_20261128/.venv/bin/ipsae --version
environments/ipsae2_20261128/.venv/bin/python -m pytest -q environments/ipsae2_20261128/source/tests
environments/ipsae2_20261128/.venv/bin/ipsae environments/ipsae2_20261128/source/tests/data/model.pdb environments/ipsae2_20261128/source/tests/data/pae.json environments/ipsae2_20261128/source/tests/data/summary.json --output /data/jpc/covalent-enzyme-design/analysis/ipsae2_qualification_20261128/environment/official_example.tsv
