"""Geometry, topology, and trajectory readers for hydration analysis."""

from __future__ import annotations

import hashlib
import json
import logging
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf_file

from nyl12_covalent.design_md_geometry import image_near

LOGGER = logging.getLogger(__name__)
BASE = Path("qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance")
ATTACK_C = 19916
CARBONYL_O = 19917
NALPHA = 13905
OG1 = 13916
LEAVING_N = 19929
ASP_RESIDUES = {968: (14559, 14560), 970: (14578, 14579)}


def sha256(path: Path) -> str:
    """Return the SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def angle(a: np.ndarray, vertex: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Calculate angles A-vertex-C in degrees for matching coordinate arrays."""
    if not (
        np.isfinite(a).all() and np.isfinite(vertex).all() and np.isfinite(c).all()
    ):
        raise ValueError("Angle coordinates must be finite")
    va, vc = a - vertex, c - vertex
    denom = np.linalg.norm(va, axis=-1) * np.linalg.norm(vc, axis=-1)
    if (denom == 0).any():
        raise ValueError("Cannot calculate angle with zero-length vector")
    cosine = np.sum(va * vc, axis=-1) / denom
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def _topology(path: Path) -> tuple[Any, list[dict[str, Any]]]:
    """Load topology and identify each water oxygen and bonded hydrogens."""
    import parmed

    topology = parmed.load_file(str(path))
    waters = []
    for residue in topology.residues:
        if residue.name not in {"WAT", "HOH", "TIP3"}:
            continue
        oxygens = [atom for atom in residue.atoms if atom.name in {"O", "OW"}]
        if len(oxygens) != 1:
            raise ValueError(f"Water residue {residue.idx + 1} lacks unique O")
        oxygen = oxygens[0]
        hydrogens = sorted(
            atom.idx
            for atom in oxygen.bond_partners
            if atom.element == 1 or atom.name.startswith("H")
        )
        if len(hydrogens) != 2:
            raise ValueError(f"Water O atom {oxygen.idx + 1} lacks two H partners")
        waters.append(
            {
                "residue": residue.idx + 1,
                "name": residue.name,
                "o": oxygen.idx,
                "h": hydrogens,
            }
        )
    return topology, waters


def _read_frames(path: Path, stage: str) -> list[dict[str, Any]]:
    """Read coordinate frames and retain source times and periodic cells."""
    frames: list[dict[str, Any]] = []
    with netcdf_file(path, "r", mmap=False) as nc:
        coords = np.asarray(nc.variables["coordinates"][:])
        rawtime = np.asarray(nc.variables["time"][:]).reshape(-1)
        if coords.ndim == 2:
            coords = coords[None, :, :]
            rawtime = np.asarray([float(np.asarray(nc.variables["time"][:]))])
            lengths = np.asarray(nc.variables["cell_lengths"][:])[None, :]
            angles = np.asarray(nc.variables["cell_angles"][:])[None, :]
        else:
            lengths = np.asarray(nc.variables["cell_lengths"][:])
            angles = np.asarray(nc.variables["cell_angles"][:])
        if not np.isfinite(coords).all() or not np.isfinite(rawtime).all():
            raise ValueError(f"Nonfinite coordinates or times in {path}")
        if len(rawtime) > 1 and not np.all(np.diff(rawtime) > 0):
            raise ValueError(f"Trajectory times are not strictly increasing in {path}")
        if not np.isfinite(lengths).all() or not np.isfinite(angles).all():
            raise ValueError(f"Nonfinite unit cell in {path}")
        if (lengths <= 0).any() or (angles <= 0).any() or (angles >= 180).any():
            raise ValueError(f"Invalid unit cell in {path}")
        for i, xyz in enumerate(coords):
            frames.append(
                {
                    "stage": stage,
                    "time_ps": float(rawtime[i]),
                    "xyz": xyz,
                    "box": np.concatenate((lengths[i], angles[i])),
                }
            )
    return frames


def _hbond(
    oxygen: np.ndarray,
    hydrogens: np.ndarray,
    acceptors: np.ndarray,
    acceptor_indices: list[int],
    box: np.ndarray,
) -> tuple[float, float, float, int]:
    """Return nearest donor-acceptor geometry over all H and acceptor atoms."""
    oh = image_near(hydrogens, oxygen, box)
    acc = image_near(acceptors, oxygen, box)
    choices = []
    for h_index, h in enumerate(oh):
        for a_index, a in enumerate(acc):
            ha = float(np.linalg.norm(h - a))
            da = float(np.linalg.norm(oxygen - a))
            dha = float(angle(oxygen[None, :], h[None, :], a[None, :])[0])
            choices.append((ha, da, dha, acceptor_indices[a_index]))
    return min(choices, key=lambda item: item[0])


def _acyl_state(acceptance: dict[str, Any]) -> tuple[float, bool]:
    """Read acyl-retention fields from a recognized acceptance schema.

    Args:
        acceptance: Scientific acceptance record using flat producer fields or
            the historical nested ``chemistry`` object.

    Returns:
        The carbonyl C--OG distance in angstroms and acyl-retention flag.

    Raises:
        KeyError: If required schema fields are missing.
        TypeError: If the chemistry object or values have invalid types.
        ValueError: If fields are incomplete, nonfinite, or conflicting.
    """

    def read_pair(record: dict[str, Any], label: str) -> tuple[float, bool]:
        missing = {"C_OG_A", "acyl_retained"} - record.keys()
        if missing:
            raise KeyError(f"{label} lacks required fields: {sorted(missing)}")
        distance = record["C_OG_A"]
        retained = record["acyl_retained"]
        if isinstance(distance, bool) or not isinstance(distance, (int, float)):
            raise TypeError(f"{label}.C_OG_A must be numeric")
        if not math.isfinite(float(distance)):
            raise ValueError(f"{label}.C_OG_A must be finite")
        if type(retained) is not bool:
            raise TypeError(f"{label}.acyl_retained must be boolean")
        return float(distance), retained

    flat_present = "C_OG_A" in acceptance or "acyl_retained" in acceptance
    nested_present = "chemistry" in acceptance
    if not flat_present and not nested_present:
        raise KeyError("Acceptance record lacks flat or nested acyl chemistry fields")
    flat = read_pair(acceptance, "acceptance") if flat_present else None
    nested: tuple[float, bool] | None = None
    if nested_present:
        chemistry = acceptance["chemistry"]
        if not isinstance(chemistry, dict):
            raise TypeError("acceptance.chemistry must be an object")
        nested = read_pair(chemistry, "acceptance.chemistry")
    if flat is not None and nested is not None and flat != nested:
        raise ValueError("Flat and nested acyl chemistry fields conflict")
    if flat is not None:
        return flat
    if nested is None:
        raise ValueError("Acceptance record has no recognized acyl chemistry")
    return nested


def _extension_sources(base: Path) -> list[tuple[Path, str]]:
    """Return accepted extension trajectories after verifying their gates."""
    result = []
    accepted_indices = []
    directories = []
    for directory in base.glob("extend_*"):
        suffix = directory.name.removeprefix("extend_")
        if not directory.is_dir() or not suffix.isdecimal():
            continue
        index = int(suffix)
        if index < 1:
            raise ValueError(f"Extension index must be positive: {directory}")
        directories.append((index, directory))
    indices = [index for index, _ in directories]
    if len(indices) != len(set(indices)):
        raise ValueError("Duplicate numeric extension indices")
    for index, directory in sorted(directories):
        status_path = directory / "status.json"
        acceptance_path = directory / "scientific_acceptance.json"
        if not status_path.exists():
            continue
        status = json.loads(status_path.read_text())
        if status.get("status") != "accepted":
            continue
        if not acceptance_path.is_file():
            raise FileNotFoundError(acceptance_path)
        acceptance = json.loads(acceptance_path.read_text())
        numeric = acceptance["numeric"]
        boundary = acceptance["boundary"]
        acyl_distance, acyl_retained = _acyl_state(acceptance)
        temp = float(numeric["final_temperature_K"])
        gates = {
            "steps": numeric["completed_steps"] == 10000,
            "temperature": math.isfinite(temp) and 250 <= temp <= 350,
            "scc_history": numeric["electronic_history"]["verified_scc_history"]
            is True,
            "frozen_atoms": boundary["frozen_max_displacement_A"] <= 1e-6,
            "box": boundary["box_unchanged"] is True,
            "finite_temperature": acceptance["finite_temperature_passed"] is True,
            "native_process": acceptance["native_process_passed"] is True,
            "acyl_bond": acyl_distance < 1.8 and acyl_retained,
            "scc": acceptance["scc_verified"] is True,
        }
        failed = [name for name, passed in gates.items() if not passed]
        if failed:
            raise ValueError(f"Accepted {directory} failed gates: {failed}")
        trajectory = directory / "trajectory.nc"
        if not trajectory.is_file():
            raise FileNotFoundError(trajectory)
        result.append((trajectory, f"extend_{index:02d}"))
        accepted_indices.append(index)
    if accepted_indices and accepted_indices != list(
        range(1, accepted_indices[-1] + 1)
    ):
        raise ValueError("Accepted extension sequence has a missing segment")
    return result


def _plot(output: Path, rows: list[tuple[Any, ...]]) -> None:
    """Write stable distance and angle time-series plots."""
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    times = [r[2] for r in rows]
    axes[0].plot(times, [r[5] for r in rows])
    axes[0].axhline(3.5, color="r", ls="--")
    axes[0].set_ylabel("Nearest water O-C (Å)")
    axes[1].plot(times, [r[6] for r in rows])
    axes[1].axhspan(95, 125, color="orange", alpha=0.2)
    axes[1].set_ylabel("O-C-Ocarbonyl (°)")
    axes[1].set_xlabel("Time (ps)")
    fig.tight_layout()
    fig.savefig(output / "nearest_water_geometry.png", dpi=160)
    plt.close(fig)
