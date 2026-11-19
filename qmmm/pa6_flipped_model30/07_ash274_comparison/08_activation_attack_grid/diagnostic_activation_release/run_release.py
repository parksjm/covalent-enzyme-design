#!/usr/bin/env python3
"""Run the prepared single DFTB3/MM diagnostic release."""

from __future__ import annotations
import hashlib
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

FOLDER = Path(__file__).resolve().parent
TOPOLOGY = Path(
    "/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/07_ash274_comparison/00_system/ash274_solvated.prmtop"
)
PROTOCOL = FOLDER / "source_grid_protocol.json"
SANDER = Path("/data/jpc/code/ambertools26/venv/bin/sander")


def digest(path: Path) -> str:
    """Return a SHA-256 digest for an input file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    """Check immutable inputs, run Sander, and record process evidence."""
    manifest = json.loads((FOLDER / "release.json").read_text())
    mdin = FOLDER / "unrestrained.mdin"
    text = mdin.read_text()
    qm = json.loads(PROTOCOL.read_text())["qm"]
    for key, path in (
        ("starting_restart_sha256", FOLDER / "starting.ncrst"),
        ("prepared_mdin_sha256", mdin),
        ("topology_sha256", TOPOLOGY),
        ("source_protocol_sha256", PROTOCOL),
    ):
        if digest(path) != manifest[key]:
            raise ValueError(f"Input hash mismatch: {path}")
    if (
        qm["charge"] != -1
        or qm["real_atoms"] != 102
        or qm["mask"] not in text
        or "qmcharge=-1" not in text
        or "qm_theory='DFTB3'" not in text
        or "ntr=0" not in text
        or "nmropt=0" not in text
        or re.search(r"(?im)^\s*(?:DISANG|&nmr|&wt)\b", text)
    ):
        raise ValueError("QM mask/charge or unrestrained mdin preflight failed")
    if (FOLDER / "process.json").exists():
        raise FileExistsError("process.json exists; create a new attempt for reruns")
    env = os.environ.copy() | {
        "AMBERHOME": "/data/jpc/code/ambertools26/venv",
        "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
    command = [
        str(SANDER),
        "-O",
        "-i",
        mdin.name,
        "-o",
        "minimized.mdout",
        "-p",
        str(TOPOLOGY),
        "-c",
        "starting.ncrst",
        "-r",
        "optimized.ncrst",
        "-x",
        "optimization.nc",
        "-inf",
        "mdinfo",
    ]
    record = {
        "status": "running",
        "command": command,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "permission_mode": "escalated",
        "source_restart_sha256": manifest["source_restart_sha256"],
        "source_qm_mask": qm["mask"],
        "source_qm_charge": qm["charge"],
        "environment": {
            k: env[k]
            for k in (
                "AMBERHOME",
                "LD_LIBRARY_PATH",
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
            )
        },
        "input_sha256": {
            p.name: digest(p)
            for p in (FOLDER / "starting.ncrst", mdin, TOPOLOGY, PROTOCOL)
        },
    }
    (FOLDER / "process.json").write_text(json.dumps(record, indent=2) + "\n")
    with (FOLDER / "engine.log").open("w") as log:
        result = subprocess.run(
            command,
            cwd=FOLDER,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    out = (
        (FOLDER / "minimized.mdout").read_text(errors="replace")
        if (FOLDER / "minimized.mdout").exists()
        else ""
    )
    record.update(
        {
            "status": "finished" if result.returncode == 0 else "failed",
            "exit_status": result.returncode,
            "force_converged": "Optimisation converged" in out
            or "Optimization converged" in out,
            "scc_converged": "DFTBESCF=" in out
            and not re.search(
                r"(?:SCC|SCF).{0,100}(?:fail|not converg|maximum iteration)", out, re.I
            ),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "output_sha256": {
                p.name: digest(p)
                for p in FOLDER.iterdir()
                if p.is_file() and p.name not in {"process.json", "release.json"}
            },
        }
    )
    (FOLDER / "process.json").write_text(json.dumps(record, indent=2) + "\n")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
