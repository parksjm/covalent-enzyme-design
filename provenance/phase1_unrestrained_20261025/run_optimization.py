"""Release all phase-I driving and positional restraints with unchanged QM/MM."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from parmed.amber import AmberParm
from parmed.amber.mask import AmberMask


def digest(path: Path) -> str:
    """Return a file's SHA-256 identity."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, data: object) -> None:
    """Write auditable finite-valued metadata."""
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def unrestrained_input(original: str) -> str:
    """Remove only declared biases and retain the native optimization controls."""
    parts = original.split("&wt type='END' /")
    if len(parts) != 2 or "DISANG=restraints.RST" not in parts[1]:
        raise ValueError("Expected the known phase-I driving input")
    text = parts[0]
    for key in ("nmropt", "ntr"):
        text, count = re.subn(rf"\b{key}=1(?=,)", f"{key}=0", text)
        if count != 1:
            raise ValueError(f"Expected exactly one active {key}")
    for pattern in (
        r"\s*restraint_wt=10\.0,",
        r"\s*restraintmask='[^']+',",
    ):
        text, count = re.subn(pattern, "", text)
        if count != 1:
            raise ValueError("Unexpected positional-restraint definition")
    dlfind = re.findall(r"&dlfind\b.*?/", parts[1], flags=re.S)
    if len(dlfind) != 1 or "active='@*'" not in dlfind[0]:
        raise ValueError("Expected native DL-FIND with all atoms active")
    if "maxcycle=1000" not in dlfind[0]:
        raise ValueError("Unexpected original optimization budget")
    control = dlfind[0].replace("maxcycle=1000", "maxcycle=2000")
    result = "Unrestrained phase-I endpoint: same QM/MM, fixed cell\n"
    result += text[text.index("&cntrl") :] + "\n" + control + "\n"
    if any(word in result.lower() for word in ("restraint", "disang", "dumpave")):
        # The title contains 'Unrestrained', so check namelists only.
        body = result[result.index("&cntrl") :].lower()
        if any(word in body for word in ("restraint", "disang", "dumpave")):
            raise ValueError("Artificial restraint controls remain")
    return result


def geometry(parm: AmberParm) -> dict[str, float]:
    """Measure the catalytic distances without assuming an endpoint state."""
    pairs = {
        "OG1_C": (13916, 19903),
        "OG1_HG1": (13916, 13917),
        "Asp272_HG1": (14559, 13917),
        "N_HG1": (13905, 13917),
        "C_O49": (19903, 19904),
        "C_N59": (19903, 19905),
    }
    xyz = parm.coordinates
    return {
        name: float(np.linalg.norm(xyz[a - 1] - xyz[b - 1]))
        for name, (a, b) in pairs.items()
    }


def prepare(run: Path) -> Path:
    """Freeze the exact phase-I endpoint and prepare a separate relaxation."""
    output = run / "phase1_unrestrained"
    inputs = json.loads((run / "inputs.json").read_text())
    source = run / "stage1_16"
    accepted = json.loads((source / "acceptance.json").read_text())
    if not accepted["native_converged"]:
        raise ValueError("Phase-I endpoint has not converged")
    topology = Path(inputs["topology"])
    start = source / "optimized.ncrst"
    parm = AmberParm(str(topology), xyz=str(start))
    if len(parm.atoms) != 103251 or not np.isfinite(parm.coordinates).all():
        raise ValueError("Expected finite full-system phase-I endpoint")
    text = unrestrained_input((source / "minimize.mdin").read_text())
    mask = re.search(r"qmmask='([^']+)'", text)
    if mask is None or sum(AmberMask(parm, mask[1]).Selection()) != 98:
        raise ValueError("Expected the unchanged 98-atom QM region")
    output.mkdir(exist_ok=False)
    frozen = output / "source"
    frozen.mkdir()
    paths = [
        topology,
        start,
        source / "minimize.mdin",
        source / "acceptance.json",
        run / "inputs.json",
        Path(__file__),
    ]
    sources = []
    for number, path in enumerate(paths):
        target = frozen / f"{number:02d}_{path.name}"
        shutil.copy2(path, target)
        if digest(path) != digest(target):
            raise ValueError("Frozen source differs")
        sources.append(
            {"path": str(path), "sha256": digest(path), "frozen_path": str(target)}
        )
    shutil.copy2(start, output / "starting.ncrst")
    (output / "minimize.mdin").write_text(text)
    command = [
        inputs["command"][0],
        "-O",
        "-i",
        str(output / "minimize.mdin"),
        "-p",
        str(frozen / f"00_{topology.name}"),
        "-c",
        str(output / "starting.ncrst"),
        "-o",
        str(output / "minimized.mdout"),
        "-r",
        str(output / "optimized.ncrst"),
        "-x",
        str(output / "optimization.nc"),
        "-inf",
        str(output / "mdinfo"),
    ]
    write_json(
        output / "preparation.json",
        {
            "command": command,
            "environment": inputs["environment"],
            "topology": command[command.index("-p") + 1],
            "start_coordinates": str(output / "starting.ncrst"),
            "sources": sources,
            "input_sha256": digest(output / "minimize.mdin"),
            "initial_geometry_A": geometry(parm),
            "box": parm.box.tolist(),
            "restraints": "None: ntr=0, nmropt=0, no DISANG, no frozen atoms",
            "unchanged": "DFTB3/3ob-3-1, QM98/link8/charge-2, MM Hamiltonian, cell",
            "trajectory": "ntwx=5 force-evaluation frames; not physical MD time",
        },
    )
    return output


def execute(output: Path) -> None:
    """Run the frozen input and distinguish process, convergence and chemistry."""
    metadata = json.loads((output / "preparation.json").read_text())
    if (output / "process.json").exists():
        raise FileExistsError("Never overwrite an optimization attempt")
    for row in metadata["sources"]:
        if digest(Path(row["frozen_path"])) != row["sha256"]:
            raise ValueError("Frozen source changed")
    if digest(output / "minimize.mdin") != metadata["input_sha256"]:
        raise ValueError("Optimization input changed")
    record = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "exit_status": None,
        "command": metadata["command"],
    }
    write_json(output / "process.json", record)
    environment = os.environ.copy()
    environment.update(metadata["environment"])
    logging.info("Starting fully unrestrained optimization in %s", output)
    with (output / "engine.log").open("w") as handle:
        result = subprocess.run(
            metadata["command"],
            cwd=output,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    record.update(
        exit_status=result.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    write_json(output / "process.json", record)
    if result.returncode:
        raise RuntimeError(
            f"sander exited {result.returncode}; retain failed artifacts"
        )
    text = (output / "minimized.mdout").read_text()
    converged = (
        "Optimisation converged" in text
        and "WARNING: OPTIMISATION HAS NOT CONVERGED" not in text
    )
    final = AmberParm(metadata["topology"], xyz=str(output / "optimized.ncrst"))
    finite = bool(np.isfinite(final.coordinates).all())
    box_ok = bool(np.allclose(final.box, metadata["box"], rtol=0, atol=1e-10))
    write_json(
        output / "acceptance.json",
        {
            "process_completed": True,
            "native_converged": converged,
            "finite_coordinates": finite,
            "fixed_box_preserved": box_ok,
            "initial_geometry_A": metadata["initial_geometry_A"],
            "final_geometry_A": geometry(final),
            "chemical_state_validated": False,
            "scope": "Unrestrained rearrangement diagnostic; endpoint review required",
        },
    )
    if not (converged and finite and box_ok):
        raise ValueError("Optimization did not meet numerical acceptance")
    logging.info("Unrestrained optimization converged; chemistry awaits review")


def main() -> None:
    """Prepare once, then run the retained input through the tmux launcher."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "run"))
    parser.add_argument("--path", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args.path.resolve())
    else:
        execute(args.path.resolve())


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
