"""Isolated TI1 attempts with immutable inputs and checkpoint compatibility."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any

from .qm_region import QM_MASK, audit_qm_region
from .ti1_checks import check_convergence, check_geometry

SPECIFICATION = {
    "qmmask": QM_MASK,
    "qmcharge": 0,
    "qm_theory": "DFTB3",
    "scfconv": 1e-8,
    "igb": 8,
    "drms": 0.05,
    "restraint_wt": 5.0,
    "restraintmask": ":1-928,930-1312@N,CA,C,O",
    "schema": 1,
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(path: Path, value: Any) -> None:
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _snapshot(output: Path) -> dict[str, str]:
    root = Path(__file__).resolve().parents[2]
    paths = sorted((root / "src/nyl12_covalent").glob("*.py"))
    paths += sorted((root / "scripts").glob("*.py"))
    paths += [root / "pyproject.toml", root / "uv.lock"]
    hashes = {}
    for source in paths:
        relative = source.relative_to(root)
        destination = output / "source" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        hashes[str(relative)] = _sha(destination)
    return hashes


def _checkpoint(
    source: Path,
    coordinates: Path,
    identities: dict[str, str],
    target: str,
    executable: Path,
) -> dict[str, str]:
    source = source.resolve(strict=True)
    manifest_path = source / "attempt.json"
    manifest = json.loads(manifest_path.read_text())
    process = json.loads((source / "process.json").read_text())
    if manifest["sander_sha256"] != _sha(executable):
        raise ValueError("Checkpoint Sander executable mismatch")
    if process["ended_utc"] is None:
        raise ValueError("Checkpoint source has no terminal process record")
    if manifest["specification"] != SPECIFICATION or manifest["target"] != target:
        raise ValueError("Checkpoint scientific specification or target mismatch")
    for name in ("topology", "reference"):
        if manifest["inputs"][name]["sha256"] != identities[name]:
            raise ValueError(f"Checkpoint {name} mismatch")
        if _sha(source / manifest["inputs"][name]["snapshot"]) != identities[name]:
            raise ValueError(f"Checkpoint source {name} snapshot is corrupt")
    restart = source / "ti1_dftb3.rst7"
    if coordinates != restart.resolve(strict=True):
        raise ValueError("Coordinates must be the explicit source attempt checkpoint")
    digest = _sha(coordinates)
    if process["restart_sha256"] != digest:
        raise ValueError("Checkpoint hash differs from terminal process record")
    return {
        "attempt": str(source),
        "manifest_sha256": _sha(manifest_path),
        "checkpoint_sha256": digest,
    }


def run_attempt(
    *,
    topology: Path,
    coordinates: Path,
    reference: Path,
    output: Path,
    amber_bin: Path,
    target: str,
    max_cycles: int = 300,
    checkpoint_attempt: Path | None = None,
) -> Path:
    """Run a screened TI1 attempt; callers must execute through the tmux launcher.

    Args:
        topology: Product-gated source topology.
        coordinates: Starting coordinates or explicit terminal-attempt restart.
        reference: Original restraint reference, retained across continuation.
        output: New attempt directory.
        amber_bin: Explicit AmberTools executable directory.
        target: PA6 L2 or poly(X,AMA) MAXA register identifier.
        max_cycles: Positive minimization cycle budget.
        checkpoint_attempt: Prior attempt when continuing from its restart.

    Returns:
        Sander output path after convergence and chemistry screens pass.

    Raises:
        ValueError: If input validation or scientific acceptance fails.
        FileExistsError: If the attempt directory already exists.
        RuntimeError: If the engine fails.
    """
    from .qmmm import write_ti1_dftb3_input

    if target not in {"pa6_l2", "poly_xama_maxa"} or max_cycles < 1:
        raise ValueError("Explicit target and positive max_cycles are required")
    inputs = {
        name: path.resolve(strict=True)
        for name, path in {
            "topology": topology,
            "coordinates": coordinates,
            "reference": reference,
        }.items()
    }
    identities = {name: _sha(path) for name, path in inputs.items()}
    executable = (amber_bin / "sander").absolute()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise FileNotFoundError(f"Sander is not executable: {executable}")
    continuation = (
        None
        if checkpoint_attempt is None
        else _checkpoint(
            checkpoint_attempt,
            inputs["coordinates"],
            identities,
            target,
            executable,
        )
    )
    output = output.absolute()
    output.mkdir(parents=True, exist_ok=False)
    snapshots = {
        "topology": "ti1.prmtop",
        "coordinates": "start.rst7",
        "reference": "reference.rst7",
    }
    for name, source in inputs.items():
        shutil.copyfile(source, output / snapshots[name])
        if _sha(output / snapshots[name]) != identities[name]:
            raise ValueError(f"Input changed during snapshot: {source}")
    manifest = {
        "created_utc": _now(),
        "target": target,
        "specification": SPECIFICATION,
        "max_cycles": max_cycles,
        "checkpoint_source": continuation,
        "inputs": {
            name: {
                "path": str(path),
                "sha256": identities[name],
                "snapshot": snapshots[name],
            }
            for name, path in inputs.items()
        },
        "source_sha256": _snapshot(output),
        "python": sys.executable,
        "python_version": sys.version,
        "software_versions": {
            name: version(name) for name in ("numpy", "ParmEd", "rdkit")
        },
        "sander": str(executable),
        "sander_sha256": _sha(executable),
    }
    _json(output / "attempt.json", manifest)
    acceptance: dict[str, Any] = {"accepted": False, "stage": "preflight"}
    process: dict[str, Any] = {
        "started_utc": None,
        "ended_utc": None,
        "exit_status": None,
        "restart_sha256": None,
    }
    restart = output / "ti1_dftb3.rst7"
    try:
        _json(output / "qm_region.json", audit_qm_region(output / "ti1.prmtop"))
        acceptance["starting_geometry"] = check_geometry(
            output / "ti1.prmtop",
            output / "start.rst7",
            target,
        )
        # Parse the reference too, so an incompatible coordinate count fails early.
        check_geometry(output / "ti1.prmtop", output / "reference.rst7", target)
        if not acceptance["starting_geometry"]["accepted"]:
            raise ValueError(
                "Starting geometry failed TI1/relay screening; repair is required"
            )
        mdin = write_ti1_dftb3_input(output, max_cycles)
        slko = amber_bin.parent / "dat/slko"
        parameter_files = sorted(p for p in slko.rglob("*") if p.is_file())
        if not parameter_files:
            raise FileNotFoundError(f"No installed DFTB parameter files: {slko}")
        parameter_hashes = {str(p): _sha(p) for p in parameter_files}
        if checkpoint_attempt is not None:
            previous = json.loads((checkpoint_attempt / "engine_data.json").read_text())
            if previous != parameter_hashes:
                raise ValueError("Checkpoint DFTB parameter data mismatch")
        shutil.copytree(slko, output / "engine_data_snapshot")
        _json(output / "engine_data.json", parameter_hashes)
        mdout = output / "ti1_dftb3.mdout"
        command = [
            str(executable),
            "-i",
            str(mdin),
            "-o",
            str(mdout),
            "-p",
            str(output / "ti1.prmtop"),
            "-c",
            str(output / "start.rst7"),
            "-ref",
            str(output / "reference.rst7"),
            "-r",
            str(restart),
            "-inf",
            str(output / "ti1_dftb3.mdinfo"),
        ]
        environment = os.environ.copy()
        overrides = {
            "AMBERHOME": str(amber_bin.absolute().parent),
            "OMP_NUM_THREADS": "4",
            "LD_LIBRARY_PATH": str(amber_bin.absolute().parent / "lib"),
        }
        if environment.get("LD_LIBRARY_PATH"):
            overrides["LD_LIBRARY_PATH"] += ":" + environment["LD_LIBRARY_PATH"]
        environment.update(overrides)
        _json(output / "command.json", {"argv": command, "environment": overrides})
        acceptance["stage"] = "engine"
        process["started_utc"] = _now()
        with (output / "sander.log").open("x") as log:
            result = subprocess.run(
                command,
                cwd=output,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        process["exit_status"] = result.returncode
        if result.returncode:
            raise RuntimeError(f"Sander exited {result.returncode}; inspect {output}")
        acceptance["stage"] = "postflight"
        acceptance["convergence"] = check_convergence(mdout)
        acceptance["final_geometry"] = check_geometry(
            output / "ti1.prmtop", restart, target
        )
        acceptance["accepted"] = (
            acceptance["convergence"]["accepted"]
            and acceptance["final_geometry"]["accepted"]
        )
        if not acceptance["accepted"]:
            raise ValueError(
                "Process completed but convergence/chemistry acceptance failed"
            )
        return mdout
    except Exception as exc:
        acceptance["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        process["ended_utc"] = _now()
        if restart.is_file():
            process["restart_sha256"] = _sha(restart)
        _json(output / "process.json", process)
        _json(output / "acceptance.json", acceptance)
