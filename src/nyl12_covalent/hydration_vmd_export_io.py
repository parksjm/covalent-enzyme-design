"""Shared validated input mapping and cpptraj helpers for hydration VMD export."""

from __future__ import annotations

import hashlib
import json
import logging
import struct
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm, NetCDFTraj

LOGGER = logging.getLogger(__name__)
QM_RANGES = (
    (1621, 1655),
    (2727, 2734),
    (13905, 13931),
    (14555, 14560),
    (14574, 14579),
    (19914, 19917),
    (19925, 19931),
    (19937, 19941),
)
ATOM_MAP = {
    "C199": 19916,
    "O199": 19917,
    "OG1": 13916,
    "Nalpha": 13905,
    "Nleave": 19929,
    "HG1": 13917,
}
EXPECTED_ATOM_NAMES = {
    "C199": "C",
    "O199": "O",
    "OG1": "OG1",
    "Nalpha": "N",
    "Nleave": "N",
    "HG1": "HG1",
}


def sha(path: Path) -> str:
    """Return SHA256 for a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def water_residues(summary: Path, parm: AmberParm) -> list[int]:
    """Read the fixed five closest-approaching complete water residue IDs."""
    payload = json.loads(summary.read_text(encoding="utf-8"))
    waters = payload["waters"][:5]
    if len(waters) != 5:
        raise ValueError("Water analysis must provide at least five candidates")
    ids = [int(item["water_residue_index_1based"]) for item in waters]
    if len(set(ids)) != 5:
        raise ValueError("Fixed hydration water residue IDs must be unique")
    for residue_id in ids:
        residue = parm.residues[residue_id - 1]
        if residue.name != "WAT" or len(residue.atoms) != 3:
            raise ValueError(f"Residue {residue_id} is not a complete three-atom WAT")
    return ids


def selected_indices(parm: AmberParm, water_ids: list[int]) -> np.ndarray:
    """Return original zero-based atom indices retained in the movie."""
    keep_residues = set(range(1, 1317)) | set(water_ids)
    kept = np.array(
        [atom.idx for atom in parm.atoms if atom.residue.idx + 1 in keep_residues],
        dtype=int,
    )
    qm = [idx for lo, hi in QM_RANGES for idx in range(lo - 1, hi)]
    if len(qm) != 98 or not set(qm) <= set(kept):
        raise ValueError("Expected all 98 physical QM atoms in the stripped topology")
    for name, index in ATOM_MAP.items():
        atom = parm.atoms[index - 1]
        expected = EXPECTED_ATOM_NAMES[name]
        if atom.name.strip() != expected:
            raise ValueError(f"Atom mapping mismatch for {name}: {atom.name}")
    return kept


def cpptraj_input(
    topology: Path,
    restart: Path,
    heat: Path,
    equil: Path,
    water_ids: list[int],
    extensions: list[Path],
    dcd: Path,
) -> str:
    """Build trajectory-only cpptraj input with PBC imaging and coordinate stripping."""
    residues = "1-1316," + ",".join(str(value) for value in water_ids)
    lines = [
        f"parm {topology}",
        f"trajin {restart}",
        f"trajin {heat}",
        f"trajin {equil}",
    ]
    lines.extend(f"trajin {segment / 'trajectory.nc'}" for segment in extensions)
    lines.extend(
        (
            "autoimage anchor :1313",
            f"strip !:{residues}",
            f"trajout {dcd} dcd namdcell",
            "run",
            "quit",
            "",
        )
    )
    return "\n".join(lines)


def topology_cpptraj_input(topology: Path, water_ids: list[int], stripped: Path) -> str:
    """Build a topology-only cpptraj pass to write the matching stripped parm7."""
    residues = "1-1316," + ",".join(str(value) for value in water_ids)
    return "\n".join(
        (
            f"parm {topology}",
            f"parmstrip !:{residues}",
            f"parmwrite out {stripped} parm7",
            "quit",
            "",
        )
    )


def accepted_extensions(source: Path) -> list[Path]:
    """Return extension trajectory folders whose scientific acceptance gates pass."""
    accepted = []
    for segment in sorted(source.glob("extend_[0-9][0-9]")):
        status_path = segment / "status.json"
        acceptance_path = segment / "scientific_acceptance.json"
        if not status_path.exists() or not acceptance_path.exists():
            continue
        status = json.loads(status_path.read_text(encoding="utf-8"))
        result = json.loads(acceptance_path.read_text(encoding="utf-8"))
        if status.get("status") != "accepted":
            continue
        checks = (
            status.get("exit_code") == 0,
            result["numeric"]["electronic_history"]["verified_scc_history"] is True,
            result["boundary"]["box_unchanged"] is True,
            float(result["boundary"]["frozen_max_displacement_A"]) <= 1e-6,
            result["acyl_retained"] is True,
            result["scc_verified"] is True,
            result["sampling_converged"] is False,
            (segment / "trajectory.nc").is_file(),
        )
        if not all(checks):
            raise ValueError(
                f"Extension acceptance record failed validation: {segment}"
            )
        accepted.append(segment)
    return accepted


def dcd_dimensions(path: Path) -> tuple[int, int]:
    """Read DCD frame and atom counts from its Fortran records."""
    with path.open("rb") as handle:
        header = _dcd_record(handle)
        if header[:4] != b"CORD":
            raise ValueError("cpptraj did not write a coordinate DCD")
        frame_count = struct.unpack_from("<i", header, 4)[0]
        _dcd_record(handle)
        atom_count = struct.unpack("<i", _dcd_record(handle))[0]
    if frame_count < 1 or atom_count < 1:
        raise ValueError("DCD contains zero frames or atoms")
    return frame_count, atom_count


def _dcd_record(handle: Any) -> bytes:
    marker = handle.read(4)
    if len(marker) != 4:
        raise ValueError("Incomplete DCD record")
    size = struct.unpack("<i", marker)[0]
    data = handle.read(size)
    if len(data) != size or handle.read(4) != marker:
        raise ValueError("Invalid DCD Fortran record")
    return data


def trajectory_times(path: Path) -> np.ndarray:
    """Read validated native NetCDF frame times in picoseconds."""
    trajectory = NetCDFTraj.open_old(str(path))
    variable = trajectory._ncfile.variables["time"]
    units = (
        variable.units.decode() if isinstance(variable.units, bytes) else variable.units
    )
    times = np.asarray(variable[:], dtype=float)
    trajectory.close()
    if units != "picosecond" or times.ndim != 1 or not np.isfinite(times).all():
        raise ValueError(f"Invalid trajectory times or units: {path}")
    if len(times) > 1 and not np.all(np.diff(times) > 0):
        raise ValueError(f"Trajectory times are not strictly increasing: {path}")
    return times
