#!/usr/bin/env python3
"""Prepare continuation restarts from a completed dimer optimization."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import shutil
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

LOGGER = logging.getLogger(__name__)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_option(text: str, name: str) -> str:
    """Return one quoted DL-FIND namelist option."""
    matches = re.findall(rf"\b{name}\s*=\s*'([^']*)'", text, re.IGNORECASE)
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one quoted {name} option")
    return matches[0]


def read_xyz_frames(path: Path, natom: int) -> list[np.ndarray]:
    """Read a multi-frame XYZ file and validate every frame's atom count."""
    lines = path.read_text().splitlines()
    frames: list[np.ndarray] = []
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        try:
            count = int(lines[cursor].strip())
        except ValueError as error:
            raise ValueError(f"Invalid XYZ frame count at line {cursor + 1}") from error
        if count != natom:
            raise ValueError(f"XYZ has {count} atoms; expected {natom}")
        if cursor + count + 1 >= len(lines):
            raise ValueError("Truncated XYZ frame")
        cursor += 2  # atom count and comment
        coords = []
        for line in lines[cursor : cursor + count]:
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"Malformed XYZ coordinate at line {cursor + 1}")
            coords.append([float(value) for value in fields[1:4]])
            cursor += 1
        frames.append(np.asarray(coords, dtype=float))
    if not frames:
        raise ValueError(f"No coordinate frames in {path}")
    return frames


def write_restart(path: Path, coordinates: np.ndarray, reference: Rst7) -> None:
    """Write coordinates as an Amber NetCDF restart with the source box."""
    restart = Rst7(natom=len(coordinates))
    restart.coordinates = coordinates
    restart.box = reference.box
    restart.write(str(path), netcdf=True)


def prepare(base: Path, source_run: Path, output: Path) -> dict[str, object]:
    """Validate cycle60 geometry and write continuation input files."""
    base = base.resolve(strict=True)
    source_run = source_run.resolve(strict=True)
    output = output.resolve()
    if output.exists():
        raise FileExistsError(output)

    mdin_path = source_run / "dimer.mdin"
    optimized_path = source_run / "optimized.ncrst"
    path_path = source_run / "path.xyz"
    mode_path = source_run / "path_tsmode.xyz"
    topology_path = base / "inputs/system.prmtop"
    reference_path = base / "inputs/reference.ncrst"
    sources = [
        mdin_path,
        optimized_path,
        path_path,
        mode_path,
        topology_path,
        reference_path,
    ]
    for path in sources:
        if not path.is_file():
            raise FileNotFoundError(path)
    source_hashes = {str(path): sha256(path) for path in sources}

    mdin = mdin_path.read_text()
    active_mask = parse_option(mdin, "active")
    topology = AmberParm(str(topology_path))
    reference = Rst7.open(str(reference_path))
    optimized = Rst7.open(str(optimized_path))
    full_xyz = optimized.coordinates.reshape(-1, 3).copy()
    if len(full_xyz) != len(topology.atoms) or len(full_xyz) != len(
        reference.coordinates.reshape(-1, 3)
    ):
        raise ValueError(
            "Topology, reference, and optimized restart atom counts differ"
        )
    selected = np.asarray(AmberMask(topology, active_mask).Selection(), dtype=bool)
    active_indices = np.flatnonzero(selected)
    if len(active_indices) != 4165:
        raise ValueError(
            f"Active mask selected {len(active_indices)} atoms; expected 4165"
        )

    path_frames = read_xyz_frames(path_path, len(active_indices))
    mode_frames = read_xyz_frames(mode_path, len(active_indices))
    optimized_active = full_xyz[active_indices]
    path_error = float(np.max(np.abs(path_frames[-1] - optimized_active)))
    if path_error > 0.0025:
        raise ValueError(
            f"Last path frame differs from optimized restart by {path_error:.6g} A"
        )
    mode = mode_frames[-1]
    mode_norm = float(np.linalg.norm(mode))
    if not np.isfinite(mode_norm) or mode_norm <= 1e-12:
        raise ValueError("Last path mode frame has invalid norm")
    continuation = full_xyz.copy()
    continuation[active_indices] += 0.01 * mode / mode_norm
    exterior_change = float(
        np.max(np.abs(continuation[~selected] - full_xyz[~selected]))
    )
    if exterior_change != 0.0:
        raise ValueError("Continuation mode changed exterior coordinates")
    if source_hashes != {str(path): sha256(path) for path in sources}:
        raise RuntimeError("A source file changed during preparation")

    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(optimized_path, output / "start.ncrst")
    shutil.copy2(reference_path, output / "reference.ncrst")
    shutil.copy2(topology_path, output / "system.prmtop")
    write_restart(output / "ts_mode_guess.ncrst", continuation, optimized)
    manifest: dict[str, object] = {
        "source_sha256": source_hashes,
        "active_mask": active_mask,
        "active_atom_count": len(active_indices),
        "active_atom_indices_zero_based": active_indices.tolist(),
        "path_frame_count": len(path_frames),
        "mode_frame_count": len(mode_frames),
        "last_path_max_abs_error_angstrom": path_error,
        "last_mode_norm_angstrom": mode_norm,
        "displacement_norm_angstrom": 0.01,
        "exterior_max_abs_change_angstrom": exterior_change,
        "mode_geometry": "optimized full coordinates plus normalized final path_tsmode vector on active atoms",
        "prepared_sha256": {
            name: sha256(output / name)
            for name in (
                "start.ncrst",
                "reference.ncrst",
                "system.prmtop",
                "ts_mode_guess.ncrst",
            )
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Prepared continuation inputs in %s", output)
    return manifest


def main() -> None:
    """Parse CLI arguments and prepare validated continuation inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True, help="06_dimer_ts root")
    parser.add_argument(
        "--source-run", type=Path, required=True, help="cycle60 calculation directory"
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="New continuation input directory"
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    prepare(args.base, args.source_run, args.output)


if __name__ == "__main__":
    main()
