#!/usr/bin/env bash
export STATUS_PATH=/data/jpc/covalent-enzyme-design/analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/repaired_comparison_20261120/5XYG/rerun_20261123/job/status.json
exec >>/data/jpc/covalent-enzyme-design/analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/repaired_comparison_20261120/5XYG/rerun_20261123/job/job.log 2>&1
cd -- /data/jpc/covalent-enzyme-design || exit $?
export LD_LIBRARY_PATH=/data/jpc/code/pmemd26/lib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=/data/jpc/covalent-enzyme-design/analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/repaired_comparison_20261120/source_metadata_fix_20261123
/data/jpc/covalent-enzyme-design/.venv/bin/python3 -c 'import datetime,json,os; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'started_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); json.dump(d,open(p,'"'"'w'"'"'),indent=2,sort_keys=True); open(p,'"'"'a'"'"').write('"'"'\n'"'"')' || exit $?
finish_job() { JOB_EXIT=$?; export JOB_EXIT; trap - EXIT HUP INT TERM; /data/jpc/covalent-enzyme-design/.venv/bin/python3 -c 'import datetime,json,os,tempfile; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'ended_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); d['"'"'exit_status'"'"']=int(os.environ['"'"'JOB_EXIT'"'"']); fd,tmp=tempfile.mkstemp(dir=os.path.dirname(p)); f=os.fdopen(fd,'"'"'w'"'"'); json.dump(d,f,indent=2,sort_keys=True); f.write('"'"'\n'"'"'); f.close(); os.replace(tmp,p)' || exit 125; exit "$JOB_EXIT"; }
trap finish_job EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
/usr/bin/nice -n 10 /data/jpc/enzyme-substrate-modeling/.pixi/envs/analysis/bin/python /data/jpc/covalent-enzyme-design/analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/repaired_comparison_20261120/source_metadata_fix_20261123/run_comparison.py --profile /data/jpc/covalent-enzyme-design/analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/repaired_comparison_20261120/5XYG/rerun_20261123/profile.json
