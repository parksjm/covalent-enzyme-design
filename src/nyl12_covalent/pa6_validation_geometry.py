"""Place accepted PA6 QM coordinates onto a RESP/GAFF2 MOL2 candidate."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

from nyl12_covalent.chemistry import load_pa6
from nyl12_covalent.gaussian_checks import final_orientation
from nyl12_covalent.ti1_reference_graph import validate_pa6_zwitterion

LOGGER = logging.getLogger(__name__)


def prepare_pa6_validation_geometry(
    source_mol2: Path,
    typed_mol2: Path,
    gaussian_log: Path,
    acceptance_json: Path,
    output_mol2: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    """Write a typed MOL2 using one scientifically accepted DFT geometry.

    Args:
        source_mol2: Canonical named PA6 L4 zwitterion graph.
        typed_mol2: RESP-charged, GAFF2-typed MOL2 in canonical atom order.
        gaussian_log: Completed fixed-order DFT optimization log.
        acceptance_json: Existing geometry scientific-acceptance record.
        output_mol2: New coordinate-updated MOL2 path.
        manifest_path: New provenance manifest path.

    Returns:
        Coordinate-transfer provenance and geometry summary.

    Raises:
        FileExistsError: If either output already exists.
        RuntimeError: If the source geometry was not accepted.
        ValueError: If atom identity/order or MOL2 structure is inconsistent.
    """
    if output_mol2.exists() or manifest_path.exists():
        raise FileExistsError("Validation geometry outputs already exist")
    acceptance = json.loads(acceptance_json.read_text())
    if not acceptance.get("accepted"):
        raise RuntimeError("Gaussian geometry lacks positive scientific acceptance")
    molecule = load_pa6(source_mol2).mol
    validate_pa6_zwitterion(molecule)
    expected_numbers = [atom.GetAtomicNum() for atom in molecule.GetAtoms()]
    expected_names = [atom.GetProp("name") for atom in molecule.GetAtoms()]
    elements, xyz = final_orientation(gaussian_log.read_text())
    if elements != expected_numbers or xyz.shape != (len(expected_names), 3):
        raise ValueError("Gaussian atom identity/order differs from canonical PA6")
    rendered = replace_mol2_coordinates(typed_mol2.read_text(), expected_names, xyz)
    output_mol2.parent.mkdir(parents=True, exist_ok=True)
    output_mol2.write_text(rendered)
    report = {
        "schema_version": 1,
        "scope": "Accepted DFT geometry transferred to PA6 RESP/GAFF2 candidate",
        "conformer": acceptance.get("conformer"),
        "atom_count": len(expected_names),
        "h33_h36_distance_angstrom": float(np.linalg.norm(xyz[65] - xyz[68])),
        "scientific_acceptance_sha256": _sha256(acceptance_json),
        "source_sha256": {
            "canonical_mol2": _sha256(source_mol2),
            "typed_mol2": _sha256(typed_mol2),
            "gaussian_log": _sha256(gaussian_log),
        },
        "output_mol2_sha256": _sha256(output_mol2),
    }
    manifest_path.write_text(json.dumps(report, indent=2) + "\n")
    LOGGER.info("Prepared accepted-coordinate GAFF2 MOL2 at %s", output_mol2)
    return report


def replace_mol2_coordinates(
    mol2_text: str, expected_names: list[str], xyz: np.ndarray
) -> str:
    """Replace only MOL2 atom coordinates while preserving types and charges.

    Args:
        mol2_text: Complete typed MOL2 text.
        expected_names: Required atom names in canonical order.
        xyz: Replacement coordinates in Angstrom.

    Returns:
        MOL2 text with updated atom coordinates.

    Raises:
        ValueError: If the atom section, names, or coordinate array is invalid.
    """
    coordinates = np.asarray(xyz, dtype=float)
    if (
        coordinates.shape != (len(expected_names), 3)
        or not np.isfinite(coordinates).all()
    ):
        raise ValueError("Invalid replacement MOL2 coordinates")
    lines = mol2_text.splitlines()
    try:
        start = lines.index("@<TRIPOS>ATOM") + 1
        end = lines.index("@<TRIPOS>BOND", start)
    except ValueError as exc:
        raise ValueError("MOL2 lacks one ordered ATOM/BOND section") from exc
    atom_lines = lines[start:end]
    if len(atom_lines) != len(expected_names):
        raise ValueError("MOL2 atom count differs from canonical PA6")
    replacement = []
    for index, (line, name, point) in enumerate(
        zip(atom_lines, expected_names, coordinates, strict=True), start=1
    ):
        fields = line.split()
        if len(fields) < 9 or int(fields[0]) != index or fields[1] != name:
            raise ValueError(f"MOL2 atom identity/order mismatch at atom {index}")
        fields[2:5] = [f"{value:.10f}" for value in point]
        replacement.append(" ".join(fields))
    return "\n".join([*lines[:start], *replacement, *lines[end:], ""])


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
