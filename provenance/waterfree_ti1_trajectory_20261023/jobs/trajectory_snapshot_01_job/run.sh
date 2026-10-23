#!/usr/bin/env bash
export STATUS_PATH=/data/jpc/covalent-enzyme-design/output_nyl12_waterfree_reference_20261023/trajectory_snapshot_01_job/status.json
exec >>/data/jpc/covalent-enzyme-design/output_nyl12_waterfree_reference_20261023/trajectory_snapshot_01_job/job.log 2>&1
cd -- /data/jpc/covalent-enzyme-design || exit $?
export PYTHONPATH=/data/jpc/covalent-enzyme-design/src:/data/jpc/covalent-enzyme-design
export OPENBLAS_NUM_THREADS=1
/data/jpc/covalent-enzyme-design/.venv/bin/python -c 'import datetime,json,os; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'started_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); json.dump(d,open(p,'"'"'w'"'"'),indent=2,sort_keys=True); open(p,'"'"'a'"'"').write('"'"'\n'"'"')' || exit $?
finish_job() { JOB_EXIT=$?; export JOB_EXIT; trap - EXIT HUP INT TERM; /data/jpc/covalent-enzyme-design/.venv/bin/python -c 'import datetime,json,os,tempfile; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'ended_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); d['"'"'exit_status'"'"']=int(os.environ['"'"'JOB_EXIT'"'"']); fd,tmp=tempfile.mkstemp(dir=os.path.dirname(p)); f=os.fdopen(fd,'"'"'w'"'"'); json.dump(d,f,indent=2,sort_keys=True); f.write('"'"'\n'"'"'); f.close(); os.replace(tmp,p)' || exit 125; exit "$JOB_EXIT"; }
trap finish_job EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
/data/jpc/covalent-enzyme-design/.venv/bin/python output_nyl12_waterfree_reference_20261023/export_optimization_trajectory.py --log output_nyl12_waterfree_reference_20261023/waterfree_pilot_01/optimization.log --graph output_nyl12_waterfree_reference_20261023/waterfree_pilot_01/native_seed.sdf --atom-map output_nyl12_waterfree_reference_20261023/waterfree_pilot_01/atom_map.json --policy output_nyl12_waterfree_reference_20261023/waterfree_pilot_01/pilot_policy.json --output output_nyl12_waterfree_reference_20261023/trajectory_snapshot_01
