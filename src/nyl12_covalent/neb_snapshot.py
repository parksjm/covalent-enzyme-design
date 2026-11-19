"""Capture coherent native Amber NEB bands using ordered unbuffered output."""

from __future__ import annotations

import hashlib
import json
import math
import tempfile
import os
import re
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from nyl12_covalent.neb_progress import NebReport, parse_last_complete_neb_report

_SAVE = re.compile(rb"Saving img\s+(\d+) to restart:\s*(\S+)")
_CYCLE = re.compile(rb"Testing convergence\s+in cycle\s+(\d+)")


def _require_safe_writer(writer_pid: int) -> str:
    if writer_pid <= 0:
        raise ValueError("A positive sander writer PID is required")
    process = Path("/proc") / str(writer_pid)
    executable = (process / "exe").resolve(strict=True)
    if executable.name not in {"sander", "sander.LES"}:
        raise ValueError("Writer PID is not a supported serial sander executable")
    environment = (process / "environ").read_bytes().split(b"\0")
    enabled = {b"1", b"Y", b"YES", b"TRUE"}
    for entry in environment:
        key, _, value = entry.partition(b"=")
        if key == b"GFORTRAN_UNBUFFERED_ALL" and value.upper() in enabled:
            return "live sander with verified GFORTRAN_UNBUFFERED_ALL"
    raise ValueError(
        "Cannot prove output/write ordering without unbuffered sander output"
    )


def _require_finished_run(status_path: Path, working_directory: Path) -> str:
    """Prove a native run is terminal and has no remaining Sander writer."""
    status = json.loads(status_path.read_text())
    if status.get("ended_utc") is None or not isinstance(
        status.get("exit_status"), int
    ):
        raise ValueError("Terminal job status is required for offline NEB capture")
    expected_cwd = working_directory.resolve()
    for process in (Path("/proc")).iterdir():
        if not process.name.isdigit():
            continue
        try:
            executable = (process / "exe").resolve(strict=True)
            cwd = (process / "cwd").resolve(strict=True)
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if executable.name in {"sander", "sander.LES"} and cwd == expected_cwd:
            raise ValueError("A Sander writer remains after the job status ended")
    return (
        f"terminal job status exit={status['exit_status']}; "
        "no live Sander writer in the calculation directory"
    )


def _signature(path: Path) -> tuple[int, int, int, int, int]:
    status = path.stat()
    return (
        status.st_dev,
        status.st_ino,
        status.st_size,
        status.st_mtime_ns,
        status.st_ctime_ns,
    )


def _restart_natom(payload: bytes) -> int:
    if payload.startswith((b"CDF", b"\x89HDF")):
        import numpy as np
        from parmed.amber import Rst7

        with tempfile.TemporaryDirectory(prefix="neb_restart_check_") as directory:
            path = Path(directory) / "snapshot.ncrst"
            path.write_bytes(payload)
            restart = Rst7.open(str(path))
            xyz = restart.coordinates.reshape(-1, 3)
            if len(xyz) != restart.natom or not np.isfinite(xyz).all():
                raise ValueError("Invalid NetCDF restart coordinates")
            return restart.natom
    lines = payload.splitlines()
    if len(lines) < 3 or not payload.endswith(b"\n"):
        raise ValueError("Incomplete ASCII restart")
    natom = int(lines[1].split()[0])
    values = []
    for line in lines[2:]:
        if len(line) % 12:
            raise ValueError("Restart coordinates are not native 12-column records")
        values.extend(float(line[i : i + 12]) for i in range(0, len(line), 12))
    if natom < 1 or len(values) not in {3 * natom, 3 * natom + 6}:
        raise ValueError("Restart atom/coordinate count is invalid")
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Restart contains nonfinite coordinates")
    return natom


def snapshot_complete_band(
    mdout: Path,
    restart_base: Path,
    destination: Path,
    expected_images: int,
    *,
    writer_pid: int,
) -> NebReport:
    """Copy a band only after proving a stable completed native save loop.

    Args:
        mdout: Native Amber output from the verified live serial writer.
        restart_base: Value passed to sander -r, for example neb.rst7.
        destination: New immutable snapshot directory; must not already exist.
        expected_images: Regular image count, excluding climbing slot.
        writer_pid: Live sander PID, launched with GFORTRAN_UNBUFFERED_ALL=1.

    Returns:
        The profile associated with copied full-system restarts. JSON manifest
        records hashes, writer proof, report bytes and absolute/relative energies.

    Raises:
        ValueError: If writer ordering or complete image association is unproven.
        RuntimeError: If any input changes during the snapshot window.
        FileExistsError: If destination already exists.
        OSError: If prerequisites or source files cannot be read.
    """
    writer_proof = _require_safe_writer(writer_pid)
    paths = [
        restart_base.with_name(
            f"{restart_base.stem}.{image:04d}{restart_base.suffix or '.rst7'}"
        )
        for image in range(1, expected_images + 1)
    ]
    inputs = [mdout, *paths]
    before = [_signature(path) for path in inputs]
    payload = mdout.read_bytes()
    report = parse_last_complete_neb_report(payload, expected_images)
    tail = payload[report.report_end_byte :]
    if b"NEB Report" in tail or report.cycle is None:
        raise ValueError("Latest report has no closed native save loop")
    cycle_match = _CYCLE.search(tail)
    if cycle_match is None or int(cycle_match[1]) != report.cycle:
        raise ValueError("Following convergence marker is unavailable")
    saves = list(_SAVE.finditer(tail[: cycle_match.start()]))
    if [int(match[1]) for match in saves] != list(range(1, expected_images + 1)):
        raise ValueError("Report does not have a complete ordered restart save loop")
    for match, path in zip(saves, paths, strict=True):
        logged = Path(os.fsdecode(match[2]))
        if logged.is_absolute():
            matches = logged.resolve() == path.resolve()
        else:
            # Native sander logs its supplied restart path, relative to its cwd.
            process_cwd = (Path("/proc") / str(writer_pid) / "cwd").resolve(strict=True)
            matches = (process_cwd / logged).resolve() == path.resolve()
        if not matches:
            raise ValueError("Logged restart path does not match requested input")
    coordinates = [path.read_bytes() for path in paths]
    atom_counts = [_restart_natom(data) for data in coordinates]
    if len(set(atom_counts)) != 1:
        raise ValueError("Band restarts have inconsistent atom counts")
    after = [_signature(path) for path in inputs]
    if before != after or mdout.read_bytes() != payload:
        raise RuntimeError("NEB inputs changed during band snapshot")
    if _require_safe_writer(writer_pid) != writer_proof:
        raise RuntimeError("Writer identity changed during band snapshot")
    if [_signature(path) for path in inputs] != before:
        raise RuntimeError("NEB inputs changed before snapshot publication")
    metadata = {
        "report": asdict(report),
        "source": str(mdout.resolve()),
        "snapshot_utc": datetime.now(timezone.utc).isoformat(),
        "writer_pid": writer_pid,
        "writer_proof": writer_proof,
        "coordinate_association": "closed save loop and stable cross-file snapshot",
        "natom": atom_counts[0],
        "scientific_status": "NEB band; minima unverified",
        "restarts": [
            {
                "image": index,
                "source": str(path.resolve()),
                "name": path.name,
                "sha256": hashlib.sha256(data).hexdigest(),
            }
            for index, (path, data) in enumerate(
                zip(paths, coordinates, strict=True), 1
            )
        ],
    }
    destination.mkdir(parents=True, exist_ok=False)
    for path, data in zip(paths, coordinates, strict=True):
        target = destination / path.name
        target.write_bytes(data)
        if (
            hashlib.sha256(target.read_bytes()).hexdigest()
            != hashlib.sha256(data).hexdigest()
        ):
            raise RuntimeError("Copied restart hash verification failed")
    (destination / "source_output_snapshot.mdout").write_bytes(payload)
    (destination / "band_snapshot.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return report
