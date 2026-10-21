#!/usr/bin/env python3
"""Snapshot and launch a relay-water repair attempt in a detached tmux job."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.jobs import launch_job  # noqa: E402

LOGGER = logging.getLogger("launch_nyl12_relay_repair")
TARGETS = ("pa6_l4", "maxa_maxa")
THREAD_ENV = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
}
VERSION_QUERY = (
    "import importlib.metadata as m, json, sys; "
    "print(json.dumps({'python': sys.version, 'executable': sys.executable, "
    "'numpy': m.version('numpy'), 'parmed': m.version('ParmEd')}))"
)


def launch_repair_attempt(
    target: str,
    topology: Path,
    coordinates: Path,
    output: Path,
    session: str,
) -> Path:
    """Snapshot required files and launch one fresh relay repair attempt.

    Args:
        target: Canonical substrate identifier (``pa6_l4`` or ``maxa_maxa``).
        topology: Existing Amber topology file.
        coordinates: Existing Amber coordinate file.
        output: New attempt directory; it must not already exist.
        session: Unique tmux session name.

    Returns:
        Path to the generated tmux run script.

    Raises:
        ValueError: If the target or output directory is invalid.
        FileNotFoundError: If an input, source file, or runtime is missing.
    """
    if target not in TARGETS:
        raise ValueError(f"target must be one of {TARGETS}, got {target!r}")
    attempt = output.expanduser().absolute()
    if attempt.exists():
        raise ValueError(f"attempt directory already exists: {attempt}")

    topology_source = _required_file(topology, "topology")
    coordinates_source = _required_file(coordinates, "coordinates")
    python = ROOT / ".venv/bin/python"
    if not python.is_file() or not os.access(python, os.X_OK):
        raise FileNotFoundError(f"repository Python is missing or not executable: {python}")
    cli = _required_file(ROOT / "scripts/repair_nyl12_relay_water.py", "repair CLI")
    module_dir = ROOT / "src/nyl12_covalent"
    _required_file(module_dir / "relay_repair_geometry.py", "repair kernel")
    _required_file(module_dir / "relay_repair.py", "repair module")
    modules = sorted(module_dir.glob("*.py"))
    replay_sources = [
        cli,
        ROOT / "scripts/launch_nyl12_relay_repair.py",
        *modules,
        _required_file(ROOT / "pyproject.toml", "project configuration"),
        _required_file(ROOT / "uv.lock", "dependency lockfile"),
    ]
    versions = _runtime_versions(python)
    git = _git_snapshot(ROOT)

    attempt.mkdir(parents=True)
    inputs_dir = attempt / "inputs"
    source_dir = attempt / "source"
    input_records = [
        _copy_record(topology_source, inputs_dir / "ti1.prmtop", attempt),
        _copy_record(coordinates_source, inputs_dir / "ti1.inpcrd", attempt),
    ]
    source_records = []
    for source in replay_sources:
        destination = source_dir / source.relative_to(ROOT)
        source_records.append(_copy_record(source, destination, attempt))

    if git["diff"]:
        diff_path = source_dir / "working_tree.diff"
        diff_path.parent.mkdir(parents=True, exist_ok=True)
        diff_path.write_bytes(git["diff"])
        git["dirty_diff"] = {
            "path": str(diff_path.relative_to(attempt)),
            "sha256": _sha256(diff_path),
        }
    git.pop("diff", None)
    if git["status"]:
        status_path = source_dir / "working_tree_status.txt"
        status_path.write_text(git["status"], encoding="utf-8")
        git["status_snapshot"] = {
            "path": str(status_path.relative_to(attempt)),
            "sha256": _sha256(status_path),
        }
    git.pop("status", None)

    manifest = {
        "target": target,
        "session": session,
        "attempt_directory": str(attempt),
        "python_executable": str(python),
        "software_versions": versions,
        "git": git,
        "inputs": input_records,
        "source_snapshot": source_records,
        "thread_environment": THREAD_ENV,
    }
    manifest_path = attempt / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    command = [
        str(python),
        str(source_dir / "scripts/repair_nyl12_relay_water.py"),
        "--topology",
        str(inputs_dir / "ti1.prmtop"),
        "--coordinates",
        str(inputs_dir / "ti1.inpcrd"),
        "--target",
        target,
        "--output",
        str(attempt / "result"),
        "--step-degrees",
        "15",
    ]
    hashed_inputs = [topology_source, coordinates_source, manifest_path]
    hashed_inputs.extend(Path(record["snapshot_path"]) for record in [*input_records, *source_records])
    hashed_inputs.extend(
        path
        for path in (source_dir / "working_tree.diff", source_dir / "working_tree_status.txt")
        if path.is_file()
    )
    run_script = launch_job(
        session,
        attempt / "job",
        attempt,
        command,
        env_overrides={**THREAD_ENV, "PYTHONPATH": str(source_dir / "src")},
        input_paths=hashed_inputs,
    )
    LOGGER.info("Launched %s relay repair in tmux session %s", target, session)
    return run_script


def _required_file(path: Path, label: str) -> Path:
    """Resolve a required input and ensure it is a regular file."""
    resolved = path.expanduser().resolve(strict=True)
    if not resolved.is_file():
        raise ValueError(f"{label} is not a regular file: {resolved}")
    return resolved


def _copy_record(source: Path, destination: Path, attempt: Path) -> dict[str, str]:
    """Copy a file and record original and snapshot identities."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    original_hash = _sha256(source)
    snapshot_hash = _sha256(destination)
    if original_hash != snapshot_hash:
        raise RuntimeError(f"Snapshot copy digest mismatch for {source}")
    return {
        "original_path": str(source),
        "original_sha256": original_hash,
        "snapshot_path": str(destination),
        "snapshot_relative_path": str(destination.relative_to(attempt)),
        "snapshot_sha256": snapshot_hash,
    }


def _sha256(path: Path) -> str:
    """Return a streamed SHA-256 digest for a regular file."""
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _runtime_versions(python: Path) -> dict[str, str]:
    """Read Python, NumPy, and ParmEd versions from the job interpreter."""
    result = subprocess.run(
        [str(python), "-c", VERSION_QUERY],
        check=True,
        capture_output=True,
        text=True,
    )
    versions: Any = json.loads(result.stdout)
    required = {"python", "executable", "numpy", "parmed"}
    if not isinstance(versions, dict) or not required <= versions.keys():
        raise ValueError("repository Python returned incomplete runtime versions")
    return {key: str(value) for key, value in versions.items()}


def _git_snapshot(repo: Path) -> dict[str, Any]:
    """Capture Git HEAD, dirty status, and diff when Git is available."""
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repo, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"],
        cwd=repo, check=True, capture_output=True, text=True,
    ).stdout
    diff = subprocess.run(
        ["git", "diff", "--binary", "HEAD"], cwd=repo, check=True,
        capture_output=True,
    ).stdout
    return {"available": True, "head": head, "status": status, "diff": diff}


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=TARGETS, required=True)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--session", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Parse required arguments and launch a relay repair attempt."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _parser().parse_args(argv)
    run_script = launch_repair_attempt(
        args.target, args.topology, args.coordinates, args.output, args.session
    )
    LOGGER.info("Run script: %s", run_script)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
