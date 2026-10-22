#!/usr/bin/env bash
export STATUS_PATH=/data/jpc/covalent-enzyme-design/output_nyl12_parameter_validation_20261022/curvature_waiter_verification_01/job/status.json
exec >>/data/jpc/covalent-enzyme-design/output_nyl12_parameter_validation_20261022/curvature_waiter_verification_01/job/job.log 2>&1
cd -- /data/jpc/covalent-enzyme-design || exit $?
export OPENBLAS_NUM_THREADS=1
/data/jpc/covalent-enzyme-design/.venv/bin/python -c 'import datetime,json,os; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'started_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); json.dump(d,open(p,'"'"'w'"'"'),indent=2,sort_keys=True); open(p,'"'"'a'"'"').write('"'"'\n'"'"')' || exit $?
finish_job() { JOB_EXIT=$?; export JOB_EXIT; trap - EXIT HUP INT TERM; /data/jpc/covalent-enzyme-design/.venv/bin/python -c 'import datetime,json,os,tempfile; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'ended_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); d['"'"'exit_status'"'"']=int(os.environ['"'"'JOB_EXIT'"'"']); fd,tmp=tempfile.mkstemp(dir=os.path.dirname(p)); f=os.fdopen(fd,'"'"'w'"'"'); json.dump(d,f,indent=2,sort_keys=True); f.write('"'"'\n'"'"'); f.close(); os.replace(tmp,p)' || exit 125; exit "$JOB_EXIT"; }
trap finish_job EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
/usr/bin/bash -c 'set -e; uv run ruff check scripts/wait_nyl12_ti1_curvature.py tests/test_wait_ti1_curvature.py; uv run pytest tests/test_wait_ti1_curvature.py --cov=scripts.wait_nyl12_ti1_curvature --cov-report=term-missing --cov-report=json:output_nyl12_parameter_validation_20261022/curvature_waiter_verification_01/coverage.json'
