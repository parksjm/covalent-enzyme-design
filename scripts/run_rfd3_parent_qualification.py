"""Run one immutable native RFD3 attempt, then independent derived acceptance."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import subprocess
import time

LOGGER = logging.getLogger(__name__)


def save(path: Path, data: object) -> None:
    """Create a retained process or command record."""
    with path.open("x") as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write("\n")


def journal(manifest: dict, phase: str, action: str, result: str = "") -> None:
    """Record exact native and audit execution intent and outcome."""
    subprocess.run(
        [
            "/usr/bin/python3",
            manifest["journal_script"],
            "--log",
            manifest["journal"],
            "--agent",
            "rfd3_export_and_parent03",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "PI authorized bounded second original TS parent qualification",
            "--rollback",
            "Retain native outputs; no automatic retries or other job changes",
            "--result",
            result,
        ],
        check=True,
    )


def main() -> None:
    """Execute one hash-bound native job and preserve separate scientific status."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    manifest = json.loads(args.manifest.read_text())
    for path, expected in manifest["source_sha256"].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Frozen prerequisite changed: {path}")
    if len(manifest["jobs"]) != 1:
        raise ValueError("This launcher authorizes exactly one native attempt")
    job = manifest["jobs"][0]
    output = Path(job["output"])
    journal(manifest, "intent", f"Run exact native argv {job['argv']}")
    output.mkdir()
    save(output / "command.json", job["argv"])
    started = time.time()
    with (output / "native.log").open("x") as handle:
        native = subprocess.run(
            job["argv"], stdout=handle, stderr=subprocess.STDOUT, check=False
        )
    save(
        output / "process_acceptance.json",
        {
            "exit_code": native.returncode,
            "elapsed_seconds": time.time() - started,
            "argv": job["argv"],
        },
    )
    journal(
        manifest,
        "outcome",
        "Native parent03 qualification completed",
        f"exit={native.returncode}",
    )
    if native.returncode:
        raise subprocess.CalledProcessError(native.returncode, job["argv"])
    command = [
        manifest["python"],
        manifest["audit_script"],
        str(args.manifest),
        job["job_id"],
        manifest["derived_export"],
    ]
    journal(
        manifest,
        "intent",
        f"Run source-graph and unchanged scientific export gates: {command}",
    )
    with (output / "derived_audit.log").open("x") as handle:
        audit = subprocess.run(
            command, stdout=handle, stderr=subprocess.STDOUT, check=False
        )
    save(
        output / "derived_audit_process.json",
        {"exit_code": audit.returncode, "argv": command},
    )
    journal(
        manifest,
        "outcome",
        "Derived scientific audit completed",
        f"exit={audit.returncode}; receipt={manifest['derived_export']}/scientific_acceptance.json",
    )
    if audit.returncode:
        raise subprocess.CalledProcessError(audit.returncode, command)
    LOGGER.info("Native and independent scientific acceptance completed")


if __name__ == "__main__":
    main()
