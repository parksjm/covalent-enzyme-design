"""Validate preservation and completion of an MCCE step 1 preparation."""

from __future__ import annotations

import hashlib
import logging
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)


def audit_preparation(
    folder: Path,
    input_pdb: Path,
    expected_termini: dict[str, dict[str, int]],
) -> dict[str, Any]:
    """Audit MCCE step1 outputs against their canonical protein input.

    Args:
        folder: Directory containing ``step1_out.pdb``, ``step1.log`` and
            integer ``step1.exit``.
        input_pdb: Validated canonical heavy-protein PDB used as MCCE input.
        expected_termini: Mapping of every chain to exact ``n`` and ``c``
            terminal residue numbers.

    Returns:
        A JSON-compatible acceptance report. Scientific and process failures
        are reported in ``errors``; missing artifacts and malformed PDB raise.

    Raises:
        FileNotFoundError: A required input or output artifact is missing.
        ValueError: A PDB record or ``step1.exit`` is malformed.
    """
    output_path = folder / "step1_out.pdb"
    log_path = folder / "step1.log"
    exit_path = folder / "step1.exit"
    for path in (input_pdb, output_path, log_path, exit_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if not expected_termini or any(
        set(value) != {"n", "c"}
        or any(type(number) is not int for number in value.values())
        for value in expected_termini.values()
    ):
        raise ValueError("expected_termini must define integer n and c per chain")

    input_atoms = _read_atoms(input_pdb)
    native_atoms = _read_atoms(output_path)
    exit_text = exit_path.read_text(encoding="utf-8").strip()
    try:
        exit_code = int(exit_text)
    except ValueError:
        exit_code = None

    input_chains = {chain for chain, _, _ in input_atoms}
    expected_chains = set(expected_termini)
    native_chains = _read_chain_ids(output_path)
    original_identities = set(input_atoms)
    native_identities = set(native_atoms)
    extras = sorted(native_identities - original_identities)
    minimum_distances: list[float] = []
    missing: list[tuple[str, int, str]] = []
    for identity, original_xyz in input_atoms.items():
        candidates = native_atoms.get(identity, [])
        if not candidates:
            missing.append(identity)
            continue
        distance = min(_distance(original_xyz[0], xyz) for xyz in candidates)
        minimum_distances.append(distance)
        if distance > 0.002:
            missing.append(identity)

    observed_ntr = _terminal_residues("NTR", native_path=output_path)
    observed_ntg = _terminal_residues("NTG", native_path=output_path)
    observed_ctr = _terminal_residues("CTR", native_path=output_path)
    expected_n = {(chain, termini["n"]) for chain, termini in expected_termini.items()}
    expected_c = {(chain, termini["c"]) for chain, termini in expected_termini.items()}
    observed_n = {
        (item["chain"], item["residue_number"]) for item in observed_ntr + observed_ntg
    }
    observed_c = {(item["chain"], item["residue_number"]) for item in observed_ctr}
    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    bad_log_lines = [
        line.strip()
        for line in log_text.splitlines()
        if any(marker in line for marker in ("Error!", "ERROR:", "FATAL"))
    ]
    criteria: dict[str, Any] = {
        "exit_code_is_integer": exit_code is not None,
        "exit_code_zero": exit_code == 0,
        "log_has_no_errors": not bad_log_lines,
        "input_chain_count_matches": input_chains == expected_chains,
        "native_chain_count_matches": native_chains == expected_chains,
        "n_termini_match": observed_n == expected_n,
        "c_termini_match": observed_c == expected_c,
        "all_original_heavy_atoms_preserved": not missing,
        "original_coordinates_within_0_002_angstrom": not missing,
        "no_extra_native_heavy_atoms": not extras,
        "max_minimum_coordinate_displacement_angstrom": max(
            minimum_distances, default=None
        ),
        "missing_or_moved_atoms": [_identity_json(item) for item in missing],
        "extra_native_heavy_atoms": [_identity_json(item) for item in extras],
        "error_log_lines": bad_log_lines,
    }
    errors = [
        key for key, value in criteria.items() if isinstance(value, bool) and not value
    ]
    report: dict[str, Any] = {
        "accepted": not errors,
        "criteria": criteria,
        "expected_termini": expected_termini,
        "observed_NTR": observed_ntr,
        "observed_NTG": observed_ntg,
        "observed_CTR": observed_ctr,
        "counts": {
            "input_chains": len(input_chains),
            "expected_chains": len(expected_chains),
            "input_heavy_atoms": len(input_atoms),
            "native_heavy_atom_identities": len(native_identities),
            "missing_or_moved_heavy_atoms": len(missing),
            "extra_native_heavy_atom_identities": len(extras),
            "error_log_lines": len(bad_log_lines),
        },
        "errors": errors,
        "artifact_sha256": {
            str(path.name if path.parent == folder else path): _sha256(path)
            for path in (input_pdb, output_path, log_path, exit_path)
        },
    }
    LOGGER.info(
        "MCCE step1 acceptance: %s (%d failed criteria)",
        report["accepted"],
        len(errors),
    )
    return report


def _read_atoms(
    path: Path,
) -> dict[tuple[str, int, str], list[tuple[float, float, float]]]:
    """Read heavy ATOM records, retaining coordinate alternatives."""
    atoms: dict[tuple[str, int, str], list[tuple[float, float, float]]] = defaultdict(
        list
    )
    for line_number, line in enumerate(
        path.read_text(encoding="ascii").splitlines(), 1
    ):
        if not line.startswith("ATOM"):
            continue
        if not line.startswith("ATOM  "):
            raise ValueError(f"Malformed ATOM record in {path} line {line_number}")
        if len(line) < 54:
            raise ValueError(f"Truncated ATOM record in {path} line {line_number}")
        atom_name = line[12:16].strip().upper()
        if (
            not atom_name
            or atom_name.startswith("H")
            or (len(atom_name) > 1 and atom_name[0].isdigit() and atom_name[1] == "H")
        ):
            continue
        try:
            chain = line[21:22].strip()
            residue_number = int(line[22:26])
            xyz = tuple(float(line[start : start + 8]) for start in (30, 38, 46))
        except ValueError as error:
            raise ValueError(
                f"Malformed ATOM record in {path} line {line_number}"
            ) from error
        if not chain or not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"Malformed ATOM record in {path} line {line_number}")
        atoms[(chain, residue_number, atom_name)].append(xyz)  # type: ignore[arg-type]
    if not atoms:
        raise ValueError(f"No heavy ATOM records in {path}")
    return dict(atoms)


def _terminal_residues(
    residue_name: str,
    native_path: Path,
) -> list[dict[str, Any]]:
    """Collect unique chain/residue pseudo-terminal records from PDB text."""
    residues: set[tuple[str, int]] = set()
    for line_number, line in enumerate(
        native_path.read_text(encoding="ascii").splitlines(), 1
    ):
        if not line.startswith("ATOM  ") or line[17:20].strip() != residue_name:
            continue
        if len(line) < 54:
            raise ValueError(
                f"Truncated terminal record in {native_path} line {line_number}"
            )
        try:
            residues.add((line[21:22].strip(), int(line[22:26])))
        except ValueError as error:
            raise ValueError(
                f"Malformed terminal record in {native_path} line {line_number}"
            ) from error
    return [
        {"chain": chain, "residue_number": number} for chain, number in sorted(residues)
    ]


def _read_chain_ids(path: Path) -> set[str]:
    """Return chain identifiers from all native ATOM records, including H."""
    chains: set[str] = set()
    for line_number, line in enumerate(
        path.read_text(encoding="ascii").splitlines(), 1
    ):
        if not line.startswith("ATOM"):
            continue
        if not line.startswith("ATOM  ") or len(line) < 22:
            raise ValueError(f"Malformed ATOM record in {path} line {line_number}")
        chain = line[21:22].strip()
        if not chain:
            raise ValueError(f"Blank chain identifier in {path} line {line_number}")
        chains.add(chain)
    return chains


def _distance(
    first: tuple[float, float, float], second: tuple[float, float, float]
) -> float:
    """Return Cartesian distance in angstroms."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second, strict=True)))


def _identity_json(identity: tuple[str, int, str]) -> dict[str, Any]:
    """Convert an atom identity tuple to a report object."""
    chain, residue_number, atom_name = identity
    return {"chain": chain, "residue_number": residue_number, "atom_name": atom_name}


def _sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()
