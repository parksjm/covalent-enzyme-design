"""Intentionally stop only the verified local direct-PT umbrella campaign."""
from __future__ import annotations

import hashlib
import json
import logging
import os
from pathlib import Path
import re
import signal
import time
from datetime import datetime, timezone

ROOT = Path('/data/jpc/covalent-enzyme-design')
CAMPAIGN = ROOT / 'qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/07_direct_pt_umbrella_30ps_violet'
EVIDENCE = ROOT / 'provenance/direct_pt_umbrella_stop_20261110'
CONTROLLER = 2008323
WRAPPER = 2008314


def inventory(pid: int) -> dict[str, object]:
    """Read exact process identity from procfs."""
    proc = Path('/proc') / str(pid)
    stat = (proc / 'stat').read_text().split(') ', 1)[1].split()
    return {'pid': pid, 'ppid': int(stat[1]), 'state': stat[0],
            'starttime': stat[19], 'cwd': os.readlink(proc / 'cwd'),
            'command': (proc / 'cmdline').read_bytes().decode().rstrip('\0').split('\0')}


def save(name: str, value: object) -> None:
    """Write new evidence without overwriting existing records."""
    with (EVIDENCE / name).open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def main() -> None:
    """Verify controller and children, terminate controller first, retain state."""
    controller = inventory(CONTROLLER)
    command = controller['command']
    assert controller['cwd'] == str(ROOT)
    assert controller['ppid'] == WRAPPER
    assert command[1] == str(ROOT / 'scripts/run_qmmm_umbrella_extension.py')
    assert command[command.index('--output') + 1] == str(CAMPAIGN)
    assert command[command.index('--coordinate') + 1] == 'direct_pt'
    children = []
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit():
            continue
        try:
            row = inventory(int(proc.name))
        except (FileNotFoundError, ProcessLookupError, PermissionError) as error:
            logging.debug('Process inventory raced or inaccessible: %s', error)
            continue
        if row['ppid'] != CONTROLLER:
            continue
        args = row['command']
        assert args[0] == '/data/jpc/code/ambertools26/venv/bin/sander', row
        assert Path(row['cwd']).is_relative_to(CAMPAIGN), row
        assert args[args.index('-i') + 1] == str(Path(row['cwd']) / 'run.mdin')
        children.append(row)
    assert len(children) == 20, len(children)
    save('process_inventory_before.json', {'controller': controller,
         'wrapper': inventory(WRAPPER), 'children': children,
         'utc': datetime.now(timezone.utc).isoformat()})
    save('campaign_status_before.json', json.loads((CAMPAIGN / 'status.json').read_text()))
    # No controller signal handler: TERM prevents any subsequent block launch.
    assert inventory(CONTROLLER)['starttime'] == controller['starttime']
    os.kill(CONTROLLER, signal.SIGTERM)
    time.sleep(2)
    if (Path('/proc') / str(CONTROLLER)).exists():
        assert inventory(CONTROLLER)['state'] == 'Z', 'Controller did not stop'
    signals = []
    for row in children:
        pid = row['pid']
        if not (Path('/proc') / str(pid)).exists():
            signals.append({'pid': pid, 'signal': 'none; already exited'})
            continue
        current = inventory(pid)
        assert current['starttime'] == row['starttime']
        assert current['command'] == row['command']
        assert current['cwd'] == row['cwd']
        os.kill(pid, signal.SIGTERM)
        signals.append({'pid': pid, 'signal': 'SIGTERM'})
    time.sleep(4)
    remaining = []
    for row in children:
        proc = Path('/proc') / str(row['pid'])
        if proc.exists():
            current = inventory(row['pid'])
            if current['state'] != 'Z':
                remaining.append(current)
    save('signal_outcomes.json', {'signals': signals, 'remaining': remaining,
         'utc': datetime.now(timezone.utc).isoformat()})
    assert not remaining, remaining
    progress = []
    for window in sorted(CAMPAIGN.glob('window_*')):
        completed = []
        latest = None
        for stage in sorted(window.iterdir()):
            if not stage.is_dir():
                continue
            process = stage / 'process.json'
            if process.exists() and json.loads(process.read_text()).get('exit_status') == 0:
                completed.append(stage.name)
            mdout = stage / 'run.mdout'
            if mdout.exists():
                steps = re.findall(r'NSTEP\s*=\s*(\d+)', mdout.read_text(errors='replace'))
                latest = {'stage': stage.name, 'last_output_step': int(steps[-1]) if steps else None}
        last_complete = window / completed[-1] / 'restart.ncrst'
        progress.append({'window': window.name, 'completed_stages': completed,
                         'latest': latest, 'last_completed_restart': str(last_complete),
                         'restart_sha256': hashlib.sha256(last_complete.read_bytes()).hexdigest()})
    save('window_progress.json', progress)
    save('intentional_stop.json', {
        'state': 'intentionally_stopped', 'reason': 'Standalone direct OG1-HG1 to ThrN PT is uphill without a product basin; no current justified convergence target for further refinement. Coupled acylation and alternative pose are separate questions.',
        'utc': datetime.now(timezone.utc).isoformat(),
        'process_acceptance': 'Controller then all 20 verified Sander children terminated; no SIGKILL needed.',
        'scientific_acceptance': 'Not a converged PMF or accepted barrier; interruption is intentional, not spontaneous failure.',
        'outputs': 'All original outputs, status records, partial trajectories and restarts retained unchanged by this script.'})
    logging.info('Intentionally stopped controller and %s children; retained all artifacts.', len(children))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
