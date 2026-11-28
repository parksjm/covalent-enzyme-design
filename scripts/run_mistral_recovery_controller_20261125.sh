#!/usr/bin/env bash
# Lavender-only serving/evaluation lifecycle; never signal the Nemotron host.
set -euo pipefail
root=/data/jpc/mechanism-agent/benchmarks/nemotron_mistral_nano_20261125
attempt="$root/mistral_immediate_compat_retry01"
remote=/data/jpc/projects/mechanism-agent/benchmarks/nemotron_mistral_nano_20261125/mistral_immediate_compat_retry01
job="$attempt/controller_job"
ssh_options=(-T -i /home/jpc/.ssh/id_ed25519_aster -o BatchMode=yes
 -o IdentitiesOnly=yes -o ControlPath=none -o ConnectTimeout=10
 -o StrictHostKeyChecking=yes -o "UserKnownHostsFile=$root/deployment/nano_known_hosts")
tunnel_pid=''
finish() {
 code=$?
 trap - EXIT HUP INT TERM
 if ! ssh "${ssh_options[@]}" jpc@lavender \
  "python3 $remote/stop_owned_server.py --root $remote --model mistral" \
  > "$job/server_shutdown.log" 2>&1; then
  code=1
 fi
 if [[ -n "$tunnel_pid" ]]; then
  if ! python3 - "$tunnel_pid" <<'PY'
import logging
import os
import signal
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO)
pid = int(sys.argv[1])
proc = Path('/proc') / str(pid)
if proc.exists():
    args = (proc / 'cmdline').read_bytes().split(b'\0')
    if (proc.stat().st_uid != os.getuid() or b'-N' not in args
            or b'127.0.0.1:18112:127.0.0.1:18000' not in args
            or b'jpc@lavender' not in args):
        raise RuntimeError('Refuse to stop an unverified recovery tunnel')
    os.kill(pid, signal.SIGTERM)
else:
    logging.info('Owned recovery tunnel already exited')
PY
  then code=1; fi
 fi
 event=error
 message='Mistral-only recovery stopped with an error; inspect retained controller, canary and native server records. Nemotron was not controlled by this recovery.'
 if [[ "$code" == 0 ]]; then
  event=completed
  message='Mistral-only benchmark resume completed; scientific grading remains pending. Its owned server and tunnel were stopped; Nemotron was not controlled by this recovery.'
 fi
 /data/jpc/mechanism-agent/.venv/bin/mechanism-agent notify \
  --notifications-config /data/jpc/mechanism-agent/deployment/notifications.toml \
  --event "$event" --run "$job" --message "$message"
 exit "$code"
}
trap finish EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
ssh "${ssh_options[@]}" jpc@lavender \
 "python3 -c 'import hashlib,json,pathlib; r=json.load(open(\"$remote/runtime_patch.json\")); assert hashlib.sha256(pathlib.Path(r[\"target\"]).read_bytes()).hexdigest()==r[\"after_sha256\"]'"
ssh "${ssh_options[@]}" -o ExitOnForwardFailure=yes -o ServerAliveInterval=30 \
 -o ServerAliveCountMax=3 -N -L 127.0.0.1:18112:127.0.0.1:18000 jpc@lavender \
 > "$job/tunnel.log" 2>&1 &
tunnel_pid=$!
deadline=$((SECONDS+3600))
ready=0
while (( SECONDS < deadline )); do
 kill -0 "$tunnel_pid"
 if curl --noproxy '*' --silent --fail --max-time 4 \
  http://127.0.0.1:18112/v1/models > "$attempt/endpoint.json"; then
  ready=1
  break
 fi
 status=$(ssh "${ssh_options[@]}" jpc@lavender \
  "python3 -c 'import json; s=json.load(open(\"$remote/model_server_job/status.json\")); print(\"running\" if s[\"exit_status\"] is None else s[\"exit_status\"])'")
 [[ "$status" == running ]] || exit 1
 sleep 10
done
[[ "$ready" == 1 ]] || exit 1
cd "$attempt"
/usr/bin/timeout --signal=TERM --kill-after=60s 12h \
 /data/jpc/mechanism-agent/.venv/bin/python "$attempt/run_mistral_only_resume.py" \
 --root "$root" --output "$attempt/benchmark_resume" --port 18112
