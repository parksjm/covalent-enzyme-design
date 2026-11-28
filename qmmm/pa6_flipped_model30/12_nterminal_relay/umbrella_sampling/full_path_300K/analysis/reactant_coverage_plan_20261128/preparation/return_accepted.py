"""Return accepted Geranium stages before publishing their acceptance markers."""

from __future__ import annotations

import hashlib
import json
import logging
import subprocess
import time
from pathlib import Path
from typing import Any

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/activation_support_boundary_20261128"
)
WORK = BASE.parent / "analysis/reactant_coverage_plan_20261128/preparation"
SSH = [
    "ssh",
    "-i",
    "/home/jpc/.ssh/id_ed25519_aster",
    "-o",
    "BatchMode=yes",
    "-o",
    "ConnectTimeout=5",
]
REMOTE = r"""
from pathlib import Path
import json
b=Path('/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/activation_support_boundary_20261128')
cfg=json.loads((b/'prepared.json').read_text());stages=[];live=[]
for row in cfg['windows']:
 folder=b/row['label']
 for spec in row['stages']:
  stage=b/Path(spec['input']).parent
  marker=stage/'analysis_acceptance.json'
  if not marker.exists():continue
  receipt=json.loads(marker.read_text())
  if not receipt['accepted']:raise ValueError(str(marker))
  process=json.loads((stage/'process.json').read_text())
  if process.get('exit_status')!=0:raise ValueError(str(stage))
  files=[str(f.relative_to(b)) for f in stage.rglob('*') if f.is_file() and not f.is_symlink() and f.name!='analysis_acceptance.json']
  stages.append({'stage':str(stage.relative_to(b)),'marker':str(marker.relative_to(b)),'files':files,'csv_sha256':receipt['csv_sha256']})
 for name in ['window_status.json','diagnostic_checkpoint_10ps.json']:
  f=folder/name
  if f.exists():live.append(str(f.relative_to(b)))
for directory in [b/'controller_job',b/'job']:
 for f in directory.glob('*'):
  if f.is_file() and not f.is_symlink():live.append(str(f.relative_to(b)))
controller=json.loads((b/'controller_job/status.json').read_text())
print(json.dumps({'stages':stages,'live':live,'controller':controller}))
"""


def journal(phase: str, action: str, result: str = "") -> None:
    """Append one required collector action record or fail immediately."""
    command = [
        "python3",
        str(ROOT / "scripts/record_agent_action.py"),
        "--log",
        str(ROOT / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "reactant_coverage_return",
        "--phase",
        phase,
        "--action",
        action,
        "--permission-mode",
        "escalated",
        "--reason",
        "Return accepted boundary stages for frozen10/15ps analysis",
        "--rollback",
        "Do not delete remote originals; local accepted outputs retain source receipts",
    ]
    if result:
        command += ["--result", result]
    subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)


def sync_files(files: list[str]) -> None:
    """Copy exact listed files without deletion or partial scientific acceptance."""
    if not files:
        return
    command = [
        "rsync",
        "-a",
        "--files-from=-",
        "-e",
        " ".join(SSH),
        "geranium:" + str(BASE) + "/",
        str(BASE) + "/",
    ]
    subprocess.run(
        command,
        input="".join(name + "\n" for name in files),
        text=True,
        check=True,
        capture_output=True,
    )


def poll() -> bool:
    """Return completed native stages; publish acceptance only after local checks."""
    journal(
        "intent",
        "Read Geranium accepted-stage inventory and return new immutable stages",
    )
    remote = subprocess.run(
        [*SSH, "geranium", "python3", "-"],
        input=REMOTE,
        text=True,
        check=True,
        capture_output=True,
    )
    inventory = json.loads(remote.stdout)
    transferred = []
    for row in inventory["stages"]:
        marker = BASE / row["marker"]
        if marker.exists():
            continue
        sync_files(row["files"])
        csv_path = BASE / row["stage"] / "raw_features.csv"
        if hashlib.sha256(csv_path.read_bytes()).hexdigest() != row["csv_sha256"]:
            raise ValueError(f"CSV identity mismatch: {csv_path}")
        process = json.loads((csv_path.parent / "process.json").read_text())
        if process["exit_status"] != 0:
            raise ValueError(f"Native process did not pass: {csv_path.parent}")
        sync_files([row["marker"]])
        receipt = json.loads(marker.read_text())
        if receipt["csv_sha256"] != row["csv_sha256"] or not receipt["accepted"]:
            raise ValueError(f"Acceptance changed during return: {marker}")
        transferred.append(row["stage"])
    sync_files(inventory["live"])
    status = {
        "remote_controller": inventory["controller"],
        "returned_stages": len(inventory["stages"]),
        "new_stages": transferred,
        "poll_seconds": 300,
        "publication_policy": "Native completed-stage payload first; analysis_acceptance marker last after CSV/process checks",
        "MBAR": "Separate frozen boundary_mbar_queue owns10/15ps analysis; return does not claim scientific convergence",
    }
    temp = WORK / "return_state.tmp"
    temp.write_text(json.dumps(status, indent=2) + "\n")
    temp.replace(WORK / "return_state.json")
    journal(
        "outcome",
        "Returned accepted stages and live process records",
        json.dumps(status),
    )
    logging.info(
        "Returned %d new stages; %d total accepted",
        len(transferred),
        len(inventory["stages"]),
    )
    exit_status = inventory["controller"]["exit_status"]
    if exit_status is not None and exit_status != 0:
        raise RuntimeError(
            f"Boundary controller failed: {exit_status}; native artifacts retained"
        )
    return exit_status == 0


def main() -> None:
    """Poll every five minutes until native completion, logging failures explicitly."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    try:
        while not poll():
            time.sleep(300)
    except Exception as error:
        journal(
            "outcome",
            "Boundary return collector failed; no native jobs stopped",
            repr(error),
        )
        raise


if __name__ == "__main__":
    main()
