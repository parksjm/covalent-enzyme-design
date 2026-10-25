"""Verified canonical PDB output helpers for the analytic ES pilot."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from typing import Any, Callable

import numpy as np


class _Buffer:
    def __init__(self) -> None:
        self.text = ""

    def write_text(self, text: str) -> int:
        self.text = text
        return len(text)

    def __str__(self) -> str:
        return "PDB memory buffer"


def record(path: Path) -> dict[str, Any]:
    """Hash an input and reject changes during the read."""
    before = path.stat()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Source changed while hashing: {path}")
    return {
        "path": str(path),
        "sha256": digest,
        "bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
    }


def render(
    parm: Any,
    residues: set[int],
    waters: list[int],
    writer: Callable[..., dict],
) -> tuple[str, list[dict[str, Any]], float]:
    """Reuse the canonical writer and verify coordinates and retained atoms."""
    buffer = _Buffer()
    mapping = writer(parm, residues, waters, buffer)["atoms"]
    lines = [
        line
        for line in buffer.text.splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    indices = [row["original_atom_index"] - 1 for row in mapping]
    expected_indices = [
        atom.idx for atom in parm.atoms if atom.residue.idx + 1 in residues
    ]
    if indices != expected_indices or len(lines) != len(indices) or len(lines) > 99999:
        raise ValueError("Selected atom/H inventory or PDB serial failure")
    observed = np.array(
        [[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in lines]
    )
    error = float(np.max(np.abs(observed - np.asarray(parm.coordinates)[indices])))
    if not np.isfinite(observed).all() or error > 0.00051:
        raise ValueError(f"PDB component error exceeds .00051 A: {error}")
    for row, index in zip(mapping, indices, strict=True):
        row.update(
            atomic_number=parm.atoms[index].atomic_number,
            original_residue_name=parm.atoms[index].residue.name,
        )
    return buffer.text, mapping, error


def zipped(path: Path, text: str) -> None:
    """Write deterministic gzip without retaining duplicate raw PDBs."""
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(text.encode())
