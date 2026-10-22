"""Read nonperiodic Amber ASCII coordinate frames without inventing snapshots."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Iterator

import numpy as np
from parmed.amber import AmberParm

from .local_repair_model import identity


def read_frames(path: Path, atoms: int) -> Iterator[np.ndarray]:
    """Yield complete nonperiodic Amber ASCII frames; reject invalid data.

    Args:
        path: ASCII mdcrd file with a title and fixed-width F8.3 coordinates.
        atoms: Exact topology atom count; box coordinates are not expected.

    Yields:
        Full atom coordinates, preserving order and repeated frames.

    Raises:
        ValueError: If atom count, numeric fields or frame length are invalid.
    """
    if atoms < 1:
        raise ValueError("Positive atom count required")
    buffer = []
    count = 0
    with path.open() as stream:
        if not stream.readline():
            raise ValueError("Missing trajectory title")
        for line in stream:
            line = line.rstrip("\r\n")
            if not line or len(line) % 8:
                raise ValueError("Invalid F8.3 coordinate line")
            for offset in range(0, len(line), 8):
                value = float(line[offset : offset + 8])
                if not np.isfinite(value):
                    raise ValueError("Nonfinite trajectory coordinate")
                buffer.append(value)
                if len(buffer) == atoms * 3:
                    count += 1
                    yield np.array(buffer).reshape(atoms, 3)
                    buffer = []
    if buffer or not count:
        raise ValueError("Incomplete or empty trajectory")


def export_trajectory(topology: Path, trajectory: Path, output: Path) -> dict:
    """Write every full-system frame to compressed multi-model PDB.

    Args:
        topology: Matching nonperiodic 330-residue alpha/beta Amber topology.
        trajectory: Actual Sander minimization ASCII coordinate snapshots.
        output: New export directory.

    Returns:
        Frame count, atom count, source hashes and sampling limitations.
    """
    model = AmberParm(str(topology))
    if len(model.residues) != 330:
        raise ValueError("Expected alpha/beta diagnostic topology")
    if [model.residues[i].name for i in (232, 328, 329)] != ["THR", "LIG", "WAT"]:
        raise ValueError("Expected canonical dimer atom order; map temporary water-first trajectories back first")
    output.mkdir(parents=True, exist_ok=False)
    target = output / "optimization_trajectory.pdb.gz"
    elements = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}
    prefixes = []
    suffixes = []
    for atom in model.atoms:
        chain, residue = identity(atom.residue.idx)
        prefixes.append(
            f"ATOM  {atom.idx + 1:5d} {atom.name:>4s} {atom.residue.name:>3s} {chain}{residue:4d}    "
        )
        suffixes.append(f"  1.00  0.00          {elements[atom.atomic_number]:>2s}\n")
    count = 0
    with gzip.open(target, "wt", compresslevel=6) as stream:
        stream.write(
            "REMARK Actual minimizer snapshots; not MD time or guaranteed accepted iterates\n"
        )
        for count, xyz in enumerate(read_frames(trajectory, len(model.atoms)), 1):
            if count > 9999:
                raise ValueError("PDB MODEL numbering limit exceeded")
            stream.write(f"MODEL     {count:4d}\n")
            for prefix, suffix, point in zip(prefixes, suffixes, xyz, strict=True):
                if np.any(point < -999.999) or np.any(point > 9999.999):
                    raise ValueError("Coordinate exceeds PDB field width")
                stream.write(
                    f"{prefix}{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}{suffix}"
                )
            stream.write("ENDMDL\n")
        stream.write("END\n")
    report = {
        "frames": count,
        "atoms_per_frame": len(model.atoms),
        "full_system": True,
        "coordinate_precision_A": 0.001,
        "sampling": "Actual Sander ntwx=1 minimizer snapshots, keyed to force calls; may include trial/line-search and repeated points; not MD time",
        "source_sha256": {
            str(p.absolute()): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (topology, trajectory)
        },
        "output_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        "output_bytes": target.stat().st_size,
    }
    (output / "trajectory_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
