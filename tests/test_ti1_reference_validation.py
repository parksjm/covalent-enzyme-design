"""Tests for scoped TI1 geometry and Hessian acceptance reports."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from nyl12_covalent.ti1_bonded_model import build_ti1_bonded_model
from nyl12_covalent.ti1_reference_coordinates import (
    TORSION_LABELS,
    embed_bonded_reference,
)
from nyl12_covalent.ti1_reference_validation import (
    assess_reference_hessian,
    assess_reference_optimization,
)

TARGETS = {"phi": -120.0, "psi": 134.110276, "chi1": -58.964418}
ATOMIC_NUMBERS = {"C": 6, "H": 1, "N": 7, "O": 8}


def _reference_directory(
    directory: Path,
    *,
    charge: int = -1,
    coordinate_change: tuple[str, tuple[float, float, float]] | None = None,
    proton_transfer: bool = False,
) -> tuple[dict[str, object], np.ndarray]:
    molecule, metadata = build_ti1_bonded_model("R")
    atom_map = metadata["atom_index_by_origin"]
    molecule = embed_bonded_reference(molecule, atom_map, TARGETS, 721)
    coordinates = np.asarray(molecule.GetConformer().GetPositions(), dtype=float)
    if coordinate_change is not None:
        origin, point = coordinate_change
        coordinates[atom_map[origin]] = point
    if proton_transfer:
        hydrogen = atom_map["THR:H"]
        oxygen = coordinates[atom_map["TIL:O"]]
        direction = coordinates[hydrogen] - oxygen
        direction /= np.linalg.norm(direction)
        coordinates[hydrogen] = oxygen + 0.95 * direction
    constraints = []
    for name, labels in TORSION_LABELS.items():
        indices = [atom_map[label] + 1 for label in labels]
        observed = rdMolTransforms.GetDihedralDeg(
            molecule.GetConformer(), *(index - 1 for index in indices)
        )
        constraints.append(
            {"name": name, "atoms_one_based": indices, "degrees": float(observed)}
        )
    metadata.update(
        {
            "coordinate_reference": {
                "chain": "A",
                "residue_number": 233,
                "psi_degrees": TARGETS["psi"],
                "chi1_degrees": TARGETS["chi1"],
            },
            "constraints": constraints,
            "cap_phi_override_degrees": TARGETS["phi"],
        }
    )
    directory.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(directory / "molecule.sdf")) as writer:
        writer.write(molecule)
    (directory / "atom_map.json").write_text(json.dumps(metadata), encoding="utf-8")
    (directory / "optimization.log").write_text(
        _optimization_log(molecule, coordinates, charge), encoding="utf-8"
    )
    return metadata, coordinates


def _optimization_log(molecule: Chem.Mol, coordinates: np.ndarray, charge: int) -> str:
    rows = "\n".join(
        f"{index + 1} {atom.GetAtomicNum()} 0 "
        f"{point[0]:.8f} {point[1]:.8f} {point[2]:.8f}"
        for index, (atom, point) in enumerate(zip(molecule.GetAtoms(), coordinates))
    )
    return (
        f"Charge = {charge} Multiplicity = 1\n"
        "Input orientation:\n"
        " ---------------------------------------------------------------------\n"
        " Center Atomic Atomic Coordinates (Angstroms)\n"
        " Number Number Type X Y Z\n"
        " ---------------------------------------------------------------------\n"
        f"{rows}\n"
        " ---------------------------------------------------------------------\n"
        "Optimization completed.\n"
        "Normal termination of Gaussian 16 at Wed Oct 21 12:00:00 2026\n"
    )


def _fchk_record(label: str, kind: str, value: str) -> str:
    return f"{label:<40}   {kind}   {value}"


def _fchk_array(label: str, kind: str, values: list[str]) -> str:
    rows = [" ".join(values[start : start + 6]) for start in range(0, len(values), 6)]
    return "\n".join([_fchk_record(label, kind, f"N= {len(values)}"), *rows])


def _fchk_text(metadata: dict[str, object]) -> str:
    atoms = metadata["atoms"]
    numbers = [str(ATOMIC_NUMBERS[atom["element"]]) for atom in atoms]
    n_atoms = len(numbers)
    dimension = 3 * n_atoms
    hessian_count = dimension * (dimension + 1) // 2
    return "\n".join(
        [
            "Synthetic TI1 formatted checkpoint",
            "SP RHF/6-31G(d) FChk",
            _fchk_record("Number of atoms", "I", str(n_atoms)),
            _fchk_record("Charge", "I", "-1"),
            _fchk_record("Multiplicity", "I", "1"),
            _fchk_array("Atomic numbers", "I", numbers),
            _fchk_array("Current cartesian coordinates", "R", ["0.0"] * dimension),
            _fchk_array("Cartesian Force Constants", "R", ["0.0"] * hessian_count),
        ]
    )


def _frequency_log(count: int = 105, normal: bool = True) -> str:
    values = [-12.0, *[float(index + 1) for index in range(count - 1)]]
    lines = [
        "Frequencies -- "
        + " ".join(f"{value:.1f}" for value in values[start : start + 3])
        for start in range(0, len(values), 3)
    ]
    if normal:
        lines.append("Normal termination of Gaussian 16 at Wed Oct 21 12:01:00 2026")
    return "\n".join(lines) + "\n"


def test_accepts_constrained_reference_and_keeps_scope_limited(tmp_path: Path) -> None:
    _reference_directory(tmp_path)

    report = assess_reference_optimization(tmp_path)

    assert report["accepted"]
    assert report["geometry"]["accepted"]
    assert report["reactive_bond_screen_passed"]
    assert all(row["error_degrees"] <= 0.10 for row in report["constraints"])
    assert report["scope"] == "constrained capped QM geometry prerequisite only"
    assert report["forcefield_validated"] is False
    assert report["native_charge_transfer_eligible"] is False
    assert report["unconstrained_minimum_claimed"] is False
    assert isinstance(report["oxyanion_close_hydrogens"], list)


@pytest.mark.parametrize(
    ("charge", "coordinate_change"),
    [
        (0, None),
        (-1, ("TIL:acylCH3", (20.0, 20.0, 20.0))),
    ],
)
def test_rejects_wrong_state_or_reactive_geometry(
    tmp_path: Path,
    charge: int,
    coordinate_change: tuple[str, tuple[float, float, float]] | None,
) -> None:
    _reference_directory(tmp_path, charge=charge, coordinate_change=coordinate_change)

    report = assess_reference_optimization(tmp_path)

    assert not report["accepted"]
    if charge != -1:
        assert not report["electronic_state_passed"]
    if coordinate_change is not None:
        assert not report["reactive_bond_screen_passed"]


def test_empty_optimization_log_fails_as_malformed_input(tmp_path: Path) -> None:
    _reference_directory(tmp_path)
    (tmp_path / "optimization.log").write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="no Cartesian orientation"):
        assess_reference_optimization(tmp_path)


def test_rejects_garbage_after_gaussian_termination(tmp_path: Path) -> None:
    _reference_directory(tmp_path)
    path = tmp_path / "optimization.log"
    path.write_text(
        path.read_text(encoding="utf-8") + "trailing garbage\n", encoding="utf-8"
    )

    report = assess_reference_optimization(tmp_path)

    assert not report["accepted"]
    assert not report["optimization_completed_normally"]


def test_oxyanion_covalent_hydrogen_contact_rejects_proton_transfer(
    tmp_path: Path,
) -> None:
    _reference_directory(tmp_path, proton_transfer=True)

    report = assess_reference_optimization(tmp_path)

    assert not report["accepted"]
    assert report["oxyanion_proton_transfer_detected"]
    assert (
        min(contact["distance_A"] for contact in report["oxyanion_close_hydrogens"])
        < 1.20
    )


def test_rejects_reference_identity_and_constraint_mismatch(tmp_path: Path) -> None:
    _reference_directory(tmp_path)
    metadata_path = tmp_path / "atom_map.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["coordinate_reference"]["chain"] = "B"
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(ValueError, match="chain A"):
        assess_reference_optimization(tmp_path)


def test_accepts_finite_hessian_inventory_and_preserves_negative_modes(
    tmp_path: Path,
) -> None:
    metadata, _ = _reference_directory(tmp_path)
    (tmp_path / "hessian.fchk").write_text(_fchk_text(metadata), encoding="utf-8")
    (tmp_path / "frequency.log").write_text(_frequency_log(), encoding="utf-8")

    report = assess_reference_hessian(tmp_path)

    assert report["accepted"]
    assert report["hessian_shape"] == [111, 111]
    assert report["unprojected_frequencies_cm_inverse"][0] == -12.0
    assert report["expected_unprojected_modes"] == 105
    assert report["conditional_curvature_validation_pending"] is True
    assert report["forcefield_validated"] is False
    assert report["native_charge_transfer_eligible"] is False
    assert report["unconstrained_minimum_claimed"] is False


@pytest.mark.parametrize(
    ("count", "normal"),
    [(104, True), (105, False)],
)
def test_hessian_gate_rejects_incomplete_frequency_evidence(
    tmp_path: Path, count: int, normal: bool
) -> None:
    metadata, _ = _reference_directory(tmp_path)
    (tmp_path / "hessian.fchk").write_text(_fchk_text(metadata), encoding="utf-8")
    (tmp_path / "frequency.log").write_text(
        _frequency_log(count, normal), encoding="utf-8"
    )

    report = assess_reference_hessian(tmp_path)

    assert not report["accepted"]
    assert report["unconstrained_minimum_claimed"] is False


def test_hessian_gate_rejects_fchk_atom_order_mismatch(tmp_path: Path) -> None:
    metadata, _ = _reference_directory(tmp_path)
    text = _fchk_text(metadata).replace("7 6 6 8 6 6", "6 6 6 8 6 6", 1)
    (tmp_path / "hessian.fchk").write_text(text, encoding="utf-8")
    (tmp_path / "frequency.log").write_text(_frequency_log(), encoding="utf-8")

    with pytest.raises(ValueError, match="atomic numbers"):
        assess_reference_hessian(tmp_path)
