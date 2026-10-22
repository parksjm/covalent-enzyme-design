#!/usr/bin/env bash
export STATUS_PATH=/data/jpc/covalent-enzyme-design/output_nyl12_ti1_parameters_20261022/pa6_local_bonded/attack_R_cap_phi_minus60/conditional_curvature_controller_01/job/status.json
exec >>/data/jpc/covalent-enzyme-design/output_nyl12_ti1_parameters_20261022/pa6_local_bonded/attack_R_cap_phi_minus60/conditional_curvature_controller_01/job/job.log 2>&1
cd -- /data/jpc/covalent-enzyme-design || exit $?
export OPENBLAS_NUM_THREADS=1
/data/jpc/covalent-enzyme-design/.venv/bin/python -c 'import datetime,json,os; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'started_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); json.dump(d,open(p,'"'"'w'"'"'),indent=2,sort_keys=True); open(p,'"'"'a'"'"').write('"'"'\n'"'"')' || exit $?
finish_job() { JOB_EXIT=$?; export JOB_EXIT; trap - EXIT HUP INT TERM; /data/jpc/covalent-enzyme-design/.venv/bin/python -c 'import datetime,json,os,tempfile; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'ended_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); d['"'"'exit_status'"'"']=int(os.environ['"'"'JOB_EXIT'"'"']); fd,tmp=tempfile.mkstemp(dir=os.path.dirname(p)); f=os.fdopen(fd,'"'"'w'"'"'); json.dump(d,f,indent=2,sort_keys=True); f.write('"'"'\n'"'"'); f.close(); os.replace(tmp,p)' || exit 125; exit "$JOB_EXIT"; }
trap finish_job EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
/data/jpc/covalent-enzyme-design/.venv/bin/python /data/jpc/covalent-enzyme-design/output_nyl12_ti1_parameters_20261022/pa6_local_bonded/conditional_curvature_source/scripts/wait_nyl12_ti1_curvature.py --case /data/jpc/covalent-enzyme-design/output_nyl12_ti1_parameters_20261022/pa6_local_bonded/attack_R_cap_phi_minus60 --snapshot /data/jpc/covalent-enzyme-design/output_nyl12_ti1_parameters_20261022/pa6_local_bonded/conditional_curvature_source --journal /data/jpc/covalent-enzyme-design/output_agent_actions/20261021.jsonl --wait-hours 48 --poll-seconds 60
