#!/usr/bin/env bash
export STATUS_PATH=/data/jpc/covalent-enzyme-design/design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation/reporting_repair/job_final/status.json
exec >>/data/jpc/covalent-enzyme-design/design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation/reporting_repair/job_final/job.log 2>&1
cd -- /data/jpc/covalent-enzyme-design || exit $?
export PYTHONPATH=/data/jpc/covalent-enzyme-design:/data/jpc/covalent-enzyme-design/src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPLBACKEND=Agg
export LD_LIBRARY_PATH=/data/jpc/code/pmemd26/lib:
/data/jpc/covalent-enzyme-design/.venv/bin/python3 -c 'import datetime,json,os; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'started_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); json.dump(d,open(p,'"'"'w'"'"'),indent=2,sort_keys=True); open(p,'"'"'a'"'"').write('"'"'\n'"'"')' || exit $?
finish_job() { JOB_EXIT=$?; export JOB_EXIT; trap - EXIT HUP INT TERM; /data/jpc/covalent-enzyme-design/.venv/bin/python3 -c 'import datetime,json,os,tempfile; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'ended_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); d['"'"'exit_status'"'"']=int(os.environ['"'"'JOB_EXIT'"'"']); fd,tmp=tempfile.mkstemp(dir=os.path.dirname(p)); f=os.fdopen(fd,'"'"'w'"'"'); json.dump(d,f,indent=2,sort_keys=True); f.write('"'"'\n'"'"'); f.close(); os.replace(tmp,p)' || exit 125; exit "$JOB_EXIT"; }
trap finish_job EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
/data/jpc/covalent-enzyme-design/.venv/bin/python /data/jpc/covalent-enzyme-design/design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation/reporting_repair/publish_existing_comparison.py
