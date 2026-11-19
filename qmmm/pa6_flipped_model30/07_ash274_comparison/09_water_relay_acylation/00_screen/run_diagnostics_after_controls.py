"""Wait for bounded controls then execute read-only scientific diagnostics."""
import json
from pathlib import Path
import subprocess
import sys
import time


def main() -> None:
    """Require controls process success before running dependent diagnostics."""
    base = Path(__file__).resolve().parent
    status = base / "controls_job/status.json"
    deadline = time.monotonic() + 1800
    while time.monotonic() < deadline:
        record = json.loads(status.read_text())
        if record["ended_utc"] is not None:
            if record["exit_status"] != 0:
                raise RuntimeError("Controls failed; preserve evidence and inspect before diagnosis")
            break
        time.sleep(15)
    else:
        raise TimeoutError("Controls not complete within bounded diagnostic wait")
    root = Path('/data/jpc/covalent-enzyme-design')
    subprocess.run([sys.executable, str(root/'scripts/diagnose_nalpha_relay_screen.py'),
                    '--config', str(base/'primary_config.json'),
                    '--primary', str(base/'primary_attempt1'),
                    '--controls', str(base/'controls_attempt1'),
                    '--output', str(base/'diagnostics_attempt1')], check=True)


if __name__ == '__main__':
    main()
