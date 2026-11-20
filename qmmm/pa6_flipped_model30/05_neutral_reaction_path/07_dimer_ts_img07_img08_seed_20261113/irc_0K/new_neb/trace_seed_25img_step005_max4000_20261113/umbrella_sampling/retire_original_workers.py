#!/usr/bin/env python3
"""Retire only exact original queue-controller processes; preserve native jobs."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import signal
import time


def identity(pid: int) -> dict:
    """Read a live PID identity including its kernel start time and exact argv."""
    proc = Path('/proc') / str(pid)
    return {'pid': pid, 'argv': [x.decode() for x in (proc / 'cmdline').read_bytes().split(b'\0') if x],
            'cwd': os.readlink(proc / 'cwd'),
            'start_ticks': (proc / 'stat').read_text().rsplit(')', 1)[1].split()[19]}


def native_jobs(root: Path) -> dict:
    """Identify only native Sander jobs inside this declared campaign."""
    result = {}
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit():
            continue
        try:
            exe = os.readlink(proc / 'exe')
            cwd = Path(os.readlink(proc / 'cwd'))
        except (FileNotFoundError, PermissionError):
            continue  # Process exited or is outside this user's readable proc records.
        if Path(exe).name == 'sander' and root in cwd.parents:
            result[proc.name] = identity(int(proc.name))
    return result


def main() -> None:
    """Require a dry receipt and an explicit execute action for controller retirement."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    root = Path(__file__).resolve().parent
    prepared = root / 'prepared.json'
    cfg = json.loads(prepared.read_text())
    records = []
    for i in range(8):
        worker = root / 'workers' / f'worker_{i + 1:02d}'
        pid = int((worker / 'claim').read_text())
        rec = identity(pid)
        expected = [cfg['python'], cfg['runtime_script'], 'worker', '--config', str(prepared), '--worker', str(i)]
        if rec['argv'] != expected or rec['cwd'] != str(root):
            raise ValueError(f'Controller identity mismatch: {pid}')
        records.append(rec)
    native = native_jobs(root)
    if len(native) != 8:
        raise ValueError('Expected eight independent native windows before takeover')
    receipt = root / 'controller_retirement_plan.json'
    plan = {'created_utc': datetime.now(timezone.utc).isoformat(),
            'prepared_sha256': hashlib.sha256(prepared.read_bytes()).hexdigest(),
            'controllers': records, 'protected_native_jobs': native}
    if not args.execute:
        with receipt.open('x') as stream:
            json.dump(plan, stream, indent=2)
            stream.write('\n')
        logging.info('Dry-run validated eight exact controllers and eight native jobs')
        return
    old = json.loads(receipt.read_text())
    for key in ['prepared_sha256', 'controllers', 'protected_native_jobs']:
        if plan[key] != old[key]:
            raise ValueError(f'Identity changed since dry-run: {key}')
    for rec in records:
        if identity(rec['pid']) != rec:
            raise ValueError('Controller identity changed immediately before signaling')
        os.kill(rec['pid'], signal.SIGTERM)
    for _ in range(30):
        if not any(Path('/proc', str(rec['pid']), 'cmdline').exists() for rec in records):
            break
        time.sleep(0.1)
    if native_jobs(root) != native:
        raise ValueError('Native window identity changed during takeover; inspect before launching')
    if any(Path('/proc', str(rec['pid'])).exists() for rec in records):
        raise ValueError('Old controller remains present; inspect before launching')
    plan['retired_utc'] = datetime.now(timezone.utc).isoformat()
    plan['only_original_controllers_signaled'] = True
    plan['native_jobs_preserved'] = True
    target = root / 'controller_retirement_receipt.json'
    with target.open('x') as stream:
        json.dump(plan, stream, indent=2)
        stream.write('\n')
    logging.info('Retired eight original controllers; eight native jobs unchanged')


if __name__ == '__main__':
    main()
