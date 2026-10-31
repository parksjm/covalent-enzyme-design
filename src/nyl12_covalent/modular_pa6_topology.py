"""Identity and topology gates for the modular PA6 L4 MD release."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import parmed


PA6_SEQUENCE = ("P6N", "P6I", "P6I", "P6C")


def validate_release(release_dir: Path) -> dict[str, Any]:
    """Validate every hashed artifact in an MD-ready PA6 L4 release.

    Args:
        release_dir: Release directory containing ``manifest.json``.

    Returns:
        Parsed, identity-verified release manifest.
    """
    manifest_path = release_dir / "manifest.json"
    release = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not release.get("md_ready") or release.get("chemistry", {}).get(
        "sequence"
    ) != list(PA6_SEQUENCE):
        raise ValueError("Release is not accepted for modular PA6 L4 MD")
    for name, digest in release["artifact_sha256"].items():
        if sha256(release_dir / name) != digest:
            raise ValueError(f"Release artifact hash mismatch: {name}")
    return release


def validate_topology(
    path: Path, release_dir: Path, *, solvated: bool
) -> dict[str, Any]:
    """Validate PA6 identities, links, charges, and radii in an Amber topology.

    Args:
        path: Amber topology to inspect.
        release_dir: Identity-verified PA6 L4 release directory.
        solvated: Whether the topology must have zero total charge.

    Returns:
        Compact topology identity and size report.
    """
    topology = parmed.load_file(str(path))
    pa6_residues = list(topology.residues[1312:1316])
    if [residue.name for residue in pa6_residues] != list(PA6_SEQUENCE):
        raise ValueError("PA6 residue sequence is not P6N-P6I-P6I-P6C at 1313-1316")
    for residue, template in zip(pa6_residues, PA6_SEQUENCE, strict=True):
        expected = _mol2_atoms(release_dir / f"{template}.mol2")
        actual = [(atom.name, atom.type, float(atom.charge)) for atom in residue.atoms]
        if len(actual) != len(expected):
            raise ValueError(f"Atom count mismatch for {template}")
        for observed, reference in zip(actual, expected, strict=True):
            if (
                observed[:2] != reference[:2]
                or abs(observed[2] - reference[2]) > 1.0e-5
            ):
                raise ValueError(
                    f"Released atom mismatch for {template}: {observed} != {reference}"
                )
    residue_ids = {id(residue): index for index, residue in enumerate(pa6_residues)}
    links = {
        tuple(
            sorted(
                (
                    residue_ids[id(bond.atom1.residue)],
                    residue_ids[id(bond.atom2.residue)],
                )
            )
        )
        for bond in topology.bonds
        if id(bond.atom1.residue) in residue_ids
        and id(bond.atom2.residue) in residue_ids
        and bond.atom1.residue is not bond.atom2.residue
    }
    if links != {(0, 1), (1, 2), (2, 3)}:
        raise ValueError(f"Unexpected PA6 inter-residue bonds: {links}")
    radius_text = " ".join(
        str(value) for value in topology.parm_data.get("RADIUS_SET", [])
    )
    if "mbondi3" not in radius_text.lower():
        raise ValueError(f"Topology does not identify mbondi3 radii: {radius_text!r}")
    charge = sum(float(atom.charge) for atom in topology.atoms)
    if solvated and abs(charge) > 1.0e-4:
        raise ValueError(f"Solvated topology is not neutral: {charge:.8f}")
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "atom_count": len(topology.atoms),
        "residue_count": len(topology.residues),
        "total_charge": charge,
        "radius_set": radius_text,
    }


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of one file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _mol2_atoms(path: Path) -> list[tuple[str, str, float]]:
    lines = path.read_text(encoding="ascii").splitlines()
    start = lines.index("@<TRIPOS>ATOM") + 1
    end = lines.index("@<TRIPOS>BOND")
    atoms = []
    for line in lines[start:end]:
        fields = line.split()
        if fields:
            atoms.append((fields[1], fields[5], float(fields[8])))
    return atoms
