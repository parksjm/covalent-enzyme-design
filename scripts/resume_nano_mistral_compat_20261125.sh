#!/usr/bin/env bash
# One bounded infrastructure repair; leave the current Nemotron run untouched.
set -euo pipefail
cd /data/jpc/covalent-enzyme-design
python3 - <<'PY'
import ast
import hashlib
import json
import logging
import shlex
import subprocess
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO)
root = Path('/data/jpc/mechanism-agent/benchmarks/nemotron_mistral_nano_20261125')
name = 'mistral_grammar_compat_retry01'
attempt = root / name
previous = root / 'startup_autotune_off_campaignroot_retry01/controller_job/status.json'
deadline = time.monotonic() + 43200
while time.monotonic() < deadline:
    if json.loads(previous.read_text())['exit_status'] is not None:
        break
    logging.info('Waiting for the original benchmark to finish before runtime repair')
    time.sleep(60)
else:
    raise TimeoutError('Original benchmark exceeded the repair queue deadline')

remote = '/data/jpc/projects/mechanism-agent/benchmarks/nemotron_mistral_nano_20261125'
options = ['-i', '/home/jpc/.ssh/id_ed25519_aster', '-oBatchMode=yes',
           '-oIdentitiesOnly=yes', '-oControlPath=none', '-oConnectTimeout=10',
           '-oStrictHostKeyChecking=yes',
           '-oUserKnownHostsFile=' + str(root / 'deployment/nano_known_hosts')]
for host, model in [('marigold', 'nemotron'), ('lavender', 'mistral')]:
    code = '''import ast,hashlib,json,socket,subprocess
from pathlib import Path
r=Path(REMOTE_ROOT);a=r/ATTEMPT_NAME;a.mkdir(exist_ok=False)
active=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name','--format=csv,noheader'],text=True).strip()
if active: raise RuntimeError('Existing GPU work: '+active)
s=socket.socket();s.bind(('127.0.0.1',18000));s.close()
for n in ('.venv','weights'): (a/n).symlink_to(r/n,target_is_directory=True)
(a/(MODEL+'_spec.json')).write_bytes((r/(MODEL+'_spec.json')).read_bytes())
if MODEL=='mistral':
 p=r/'.venv/lib/python3.12/site-packages/vllm/parser/mistral.py'
 factory=ast.parse((r/'.venv/lib/python3.12/site-packages/mistral_common/guidance/grammar_factory.py').read_text())
 funcs=[n for n in ast.walk(factory) if isinstance(n,ast.FunctionDef) and n.name=='get_lark_from_jinja']
 if len(funcs)!=1 or 'force_tool_names' not in [x.arg for x in funcs[0].args.args]: raise ValueError('Unexpected mistral-common API; no patch applied')
 original=p.read_text()
 old='                    parallel_tool_calls=request.parallel_tool_calls,\\n                    json_only=False,\\n'
 if original.count(old)!=1: raise ValueError('Unexpected vLLM grammar call; no patch applied')
 new=old+'                    force_tool_names=False,\\n'
 patched=original.replace(old,new)
 tree=ast.parse(patched)
 calls=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='get_lark_from_jinja']
 if len(calls)!=1 or not any(k.arg=='force_tool_names' and isinstance(k.value,ast.Constant) and k.value.value is False for k in calls[0].keywords): raise ValueError('Invalid compatibility patch')
 (a/'vllm_mistral_original.py').write_text(original)
 (a/'vllm_mistral_patched.py').write_text(patched)
 (a/'runtime_patch.json').write_text(json.dumps({'target':str(p),'before_sha256':hashlib.sha256(original.encode()).hexdigest(),'after_sha256':hashlib.sha256(patched.encode()).hexdigest(),'change':'Explicitforce_tool_names=False; no tools, no dependency changes'},indent=2))
 p.write_text(patched)
'''.replace('REMOTE_ROOT', repr(remote)).replace('ATTEMPT_NAME', repr(name)).replace('MODEL', repr(model))
    # Validate the exact remote program, including its literal newline handling.
    ast.parse(code)
    subprocess.run(['ssh', '-T', *options, 'jpc@' + host,
                    'python3 -c ' + shlex.quote(code)], check=True, timeout=45)
    subprocess.run(['scp', *options, str(attempt / 'serve_when_ready.sh'),
                    str(attempt / 'stop_owned_server.py'),
                    'jpc@' + host + ':' + remote + '/' + name + '/'],
                   check=True, timeout=45)

command = ['/data/jpc/covalent-enzyme-design/.venv/bin/python',
           'scripts/run_nyl12_job.py', 'launch', '--session',
           'mechanism_nano_grammar_resume_20261125', '--output-dir',
           str(attempt / 'controller_job'), '--cwd', str(attempt)]
for filename in ['controller.sh', 'serve_when_ready.sh', 'source_manifest.json']:
    command += ['--input', str(attempt / filename)]
command += ['--', '/usr/bin/bash', str(attempt / 'controller.sh'), 'controller_job']
subprocess.run(command, check=True)
pointer = json.loads((root / 'current_run.json').read_text())
pointer.update(controller_job=name + '/controller_job',
               session='mechanism_nano_grammar_resume_20261125',
               phase='grammar_compatibility_resume_submitted',
               runtime_change='autotuneoff; Lavenderexplicitforce_tool_names=False',
               cases_unchanged=True, checkpoints_unchanged=True)
(root / 'current_run.json').write_text(json.dumps(pointer, indent=2) + '\n')
logging.info('Exact compatibility repair applied; hash-preserving benchmark resume submitted')
PY
