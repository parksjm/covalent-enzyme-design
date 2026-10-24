#!/usr/bin/env bash
export STATUS_PATH=/data/jpc/covalent-enzyme-design/output_nyl12_pa6_hydration_20261024/jobs/positive_control_attempt2/status.json
exec >>/data/jpc/covalent-enzyme-design/output_nyl12_pa6_hydration_20261024/jobs/positive_control_attempt2/job.log 2>&1
cd -- /data/jpc/covalent-enzyme-design || exit $?
export OMP_NUM_THREADS=1
export AMBERHOME=/data/jpc/code/ambertools26/venv
export LD_LIBRARY_PATH=/data/jpc/code/ambertools26/venv/lib
/data/jpc/covalent-enzyme-design/.venv/bin/python -c 'import datetime,json,os; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'started_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); json.dump(d,open(p,'"'"'w'"'"'),indent=2,sort_keys=True); open(p,'"'"'a'"'"').write('"'"'\n'"'"')' || exit $?
finish_job() { JOB_EXIT=$?; export JOB_EXIT; trap - EXIT HUP INT TERM; /data/jpc/covalent-enzyme-design/.venv/bin/python -c 'import datetime,json,os,tempfile; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'ended_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); d['"'"'exit_status'"'"']=int(os.environ['"'"'JOB_EXIT'"'"']); fd,tmp=tempfile.mkstemp(dir=os.path.dirname(p)); f=os.fdopen(fd,'"'"'w'"'"'); json.dump(d,f,indent=2,sort_keys=True); f.write('"'"'\n'"'"'); f.close(); os.replace(tmp,p)' || exit 125; exit "$JOB_EXIT"; }
trap finish_job EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
/data/jpc/covalent-enzyme-design/.venv/bin/python /data/jpc/covalent-enzyme-design/output_nyl12_pa6_hydration_20261024/source/positive_control_triclinic.py
