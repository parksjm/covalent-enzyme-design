#!/usr/bin/env bash
export STATUS_PATH=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/reaction_distances_20261128/local_plot_job/status.json
exec >>/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/reaction_distances_20261128/local_plot_job/job.log 2>&1
cd -- /data/jpc/covalent-enzyme-design || exit $?
export PYTHONPATH=/data/jpc/covalent-enzyme-design/src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=/data/jpc/code/pmemd26/lib:
/data/jpc/covalent-enzyme-design/.venv/bin/python -c 'import datetime,json,os; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'started_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); json.dump(d,open(p,'"'"'w'"'"'),indent=2,sort_keys=True); open(p,'"'"'a'"'"').write('"'"'\n'"'"')' || exit $?
finish_job() { JOB_EXIT=$?; export JOB_EXIT; trap - EXIT HUP INT TERM; /data/jpc/covalent-enzyme-design/.venv/bin/python -c 'import datetime,json,os,tempfile; p=os.environ['"'"'STATUS_PATH'"'"']; d=json.load(open(p)); d['"'"'ended_utc'"'"']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('"'"'+00:00'"'"','"'"'Z'"'"'); d['"'"'exit_status'"'"']=int(os.environ['"'"'JOB_EXIT'"'"']); fd,tmp=tempfile.mkstemp(dir=os.path.dirname(p)); f=os.fdopen(fd,'"'"'w'"'"'); json.dump(d,f,indent=2,sort_keys=True); f.write('"'"'\n'"'"'); f.close(); os.replace(tmp,p)' || exit 125; exit "$JOB_EXIT"; }
trap finish_job EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
/usr/bin/nice -n 15 taskset -c 0 /data/jpc/covalent-enzyme-design/.venv/bin/python /data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/reaction_distances_20261128/render_local.py
