"""Refine the unchanged restrained ES model before estimating a PES barrier."""

from __future__ import annotations
import hashlib
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
from datetime import datetime, timezone
import numpy as np
from parmed.amber import AmberParm

ROOT = Path("/data/jpc/covalent-enzyme-design/output_nyl12_pa6_dftb_pilot_20261024")
BASE = ROOT / "analytic_relaxation_continuation_attempt1"
OUT = ROOT / "es_tight_convergence_attempt1"


def sha(path: Path) -> str:
    """Return a source digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path: Path, data: object) -> None:
    """Retain process and acceptance records separately."""
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def main() -> None:
    """Execute one fresh stricter minimization with unchanged chemistry."""
    logging.basicConfig(level=logging.INFO)
    OUT.mkdir(exist_ok=False)
    source = OUT / "source"
    source.mkdir()
    metadata = json.loads((BASE / "preparation.json").read_text())
    text = (BASE / "continuation.mdin").read_text()
    changes = {"maxcyc": 3000, "drms": 0.0005, "ntpr": 25, "ntwx": 25, "ntwr": 100}
    for key, value in changes.items():
        text, count = re.subn(
            rf"(\b{key}\s*=\s*)[\d.Ee+-]+",
            lambda match: match[1] + str(value),
            text,
            count=1,
        )
        if count != 1:
            raise ValueError("Missing control " + key)
    mdin = OUT / "minimize.mdin"
    mdin.write_text(text)
    command = metadata["command"].copy()
    initial = BASE / "optimized.ncrst"
    reference = Path(metadata["restraint_reference"])
    topology = Path(command[command.index("-p") + 1])
    for flag, path in {
        "-i": mdin,
        "-c": initial,
        "-ref": reference,
        "-o": OUT / "minimized.mdout",
        "-r": OUT / "optimized.ncrst",
        "-x": OUT / "optimization.nc",
        "-inf": OUT / "mdinfo",
    }.items():
        command[command.index(flag) + 1] = str(path)
    inputs = [
        Path(__file__),
        BASE / "preparation.json",
        BASE / "continuation.mdin",
        mdin,
        initial,
        reference,
        topology,
    ]
    records = []
    for index, path in enumerate(inputs):
        copied = source / f"{index:02d}_{path.name}"
        shutil.copy2(path, copied)
        if sha(path) != sha(copied):
            raise ValueError("Snapshot changed")
        records.append(
            {"path": str(path), "sha256": sha(path), "snapshot": str(copied)}
        )
    write(
        OUT / "preparation.json",
        {
            "command": command,
            "environment": metadata["environment"],
            "sources": records,
            "changes_only": changes,
            "reference": str(reference),
            "scope": "Strict ES minimization, not TS or barrier",
        },
    )
    record = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "exit_status": None,
    }
    write(OUT / "process.json", record)
    environment = os.environ.copy()
    environment.update(metadata["environment"])
    logging.info("Starting tighter ES convergence; original MM reference retained")
    with (OUT / "engine.log").open("w") as handle:
        result = subprocess.run(
            command,
            cwd=OUT,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    record.update(
        exit_status=result.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    write(OUT / "process.json", record)
    if result.returncode:
        raise RuntimeError("Sander failed")
    rows = re.findall(
        r"^\s*(\d+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+(\S+)\s+(\d+)\s*$",
        (OUT / "minimized.mdout").read_text(),
        re.M,
    )
    if not rows:
        raise ValueError("No minimization record")
    last = rows[-1]
    parm = AmberParm(str(topology), xyz=str(OUT / "optimized.ncrst"))
    write(
        OUT / "scientific_acceptance.json",
        {
            "cycles": int(last[0]),
            "RMS": float(last[2]),
            "GMAX": float(last[3]),
            "RMS_target": 0.0005,
            "RMS_target_met": float(last[2]) <= 0.0005,
            "finite_coordinates": bool(np.isfinite(parm.coordinates).all()),
            "box_preserved": bool(
                np.allclose(parm.box, metadata["source_box"], rtol=0, atol=1e-10)
            ),
            "chemical_state_accepted": False,
            "scope": "Endpoint proton/bond review still required; no barrier established",
        },
    )


if __name__ == "__main__":
    main()
