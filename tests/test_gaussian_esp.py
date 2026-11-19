"""Exercise scientific and malformed-data rejection for standalone HF ESP."""

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from nyl12_covalent.gaussian_esp import BOHR_ANGSTROM, assess_hf_esp, read_espgen


@pytest.fixture
def inputs(tmp_path: Path) -> tuple[Path, Path, Path, np.ndarray]:
    """Create a small valid water graph and standalone HF/ESP fixtures."""
    molecule = Chem.AddHs(Chem.MolFromSmiles("O"))
    xyz = np.array([[0.0, 0.0, 0.0], [0.9572, 0.0, 0.0], [-0.24, 0.927, 0.0]])
    conformer = Chem.Conformer(3)
    for index, point in enumerate(xyz):
        conformer.SetAtomPosition(index, point)
    molecule.AddConformer(conformer)
    sdf = tmp_path / "molecule.sdf"
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(molecule)
    log = tmp_path / "gaussian.log"
    log.write_text(
        "#p HF/6-31G(d) Pop=MK IOp(6/33=2,6/42=6) NoSymm\n"
        "Standard basis: 6-31G(d) (6D, 7F)\nCharge = 0 Multiplicity = 1\n"
        "Input orientation:\n-----\nCenter Number Type X Y Z\n-----\n"
        + "\n".join(
            f"{i + 1} {z} 0 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}"
            for i, (z, p) in enumerate(zip([8, 1, 1], xyz, strict=True))
        )
        + "\n-----\nSCF Done: E(RHF) = -75.1 A.U.\n"
        "ESP Fit Center 4 is at 1 2 3\nElectrostatic Properties\n"
        "Normal termination of Gaussian 16\n"
    )
    esp = tmp_path / "esp.dat"
    esp.write_text(
        "3 1 0\n"
        + "\n".join(" ".join(f"{v:.9E}" for v in p / BOHR_ANGSTROM) for p in xyz)
        + "\n0.2 2.0 3.0 4.0\n"
    )
    return log, sdf, esp, xyz


def _assess(inputs: tuple) -> dict:
    return assess_hf_esp(*inputs, multiplicity=1, coordinate_tolerance_angstrom=2e-5)


def test_valid_standalone_esp_is_not_a_fit(inputs: tuple) -> None:
    report = _assess(inputs)
    assert report["accepted"]
    assert report["esp_grid_points"] == 1
    assert report["energy_hartree"] == -75.1
    assert not report["charges_fitted"]
    assert not report["native_charge_transfer_eligible"]


def test_fixed_width_espgen_header_with_abutting_counts(tmp_path: Path) -> None:
    esp = tmp_path / "fixed_width.esp"
    esp.write_text(
        f"{79:5d}{18391:5d}{0:5d}\n"
        + " 0.0 0.0 0.0\n" * 79
        + " 0.0 0.0 0.0 0.0\n" * 18391
    )

    centers, points = read_espgen(esp, atom_count=79)

    assert centers.shape == (79, 3)
    assert points.shape == (18391, 4)


@pytest.mark.parametrize(
    ("before", "after", "check"),
    [
        ("HF/6-31G(d)", "B3LYP/6-31G(d)", "method_and_basis"),
        ("Pop=MK", "Pop=Full", "mk_output_requested"),
        ("6/42=6", "6/42=5", "mk_output_requested"),
        ("Charge = 0", "Charge = -1", "electronic_state"),
        ("Multiplicity = 1", "Multiplicity = 3", "electronic_state"),
        ("E(RHF)", "E(UHF)", "finite_restricted_hf_energy"),
        ("ESP Fit Center", "Missing fit data", "esp_log_present"),
        ("Normal termination", "Error termination", "standalone_normal_termination"),
    ],
)
def test_failed_scientific_checks(
    inputs: tuple, before: str, after: str, check: str
) -> None:
    log = inputs[0]
    log.write_text(log.read_text().replace(before, after))
    report = _assess(inputs)
    assert not report["accepted"]
    assert not report["checks"][check]


@pytest.mark.parametrize(
    "contents",
    ["", "3 1\n", "4 1 0\n", "3 0 0\n", "3 1 1\n", "3 1 0\n0 0 0\n"],
)
def test_bad_grid_header_or_record_count(inputs: tuple, contents: str) -> None:
    inputs[2].write_text(contents)
    with pytest.raises(ValueError):
        _assess(inputs)


@pytest.mark.parametrize("record", ["0.2 1 2", "nan 1 2 3", "0.2 1 inf 3"])
def test_bad_potential_record(inputs: tuple, record: str) -> None:
    lines = inputs[2].read_text().splitlines()
    lines[-1] = record
    inputs[2].write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        _assess(inputs)


def test_wrong_grid_coordinate_units(inputs: tuple) -> None:
    lines = inputs[2].read_text().splitlines()
    lines[2] = "0.9572 0 0"
    inputs[2].write_text("\n".join(lines) + "\n")
    assert not _assess(inputs)["checks"]["esp_atom_centers_match_log"]


def test_changed_source_geometry(inputs: tuple) -> None:
    inputs[3][1, 0] += 0.01
    assert not _assess(inputs)["checks"]["source_geometry_unchanged"]


@pytest.mark.parametrize("tolerance", [0, -1, float("nan")])
def test_invalid_tolerance(inputs: tuple, tolerance: float) -> None:
    with pytest.raises(ValueError, match="tolerance"):
        assess_hf_esp(*inputs, multiplicity=1, coordinate_tolerance_angstrom=tolerance)


def test_wrong_spin_or_reference_shape(inputs: tuple) -> None:
    with pytest.raises(ValueError, match="singlet"):
        assess_hf_esp(*inputs, multiplicity=3, coordinate_tolerance_angstrom=2e-5)
    with pytest.raises(ValueError, match="reference"):
        assess_hf_esp(
            *inputs[:3],
            np.zeros((2, 3)),
            multiplicity=1,
            coordinate_tolerance_angstrom=2e-5,
        )


def test_combined_job_is_not_standalone(inputs: tuple) -> None:
    log = inputs[0]
    log.write_text(log.read_text() + "Normal termination of Gaussian 16\n")
    assert not _assess(inputs)["checks"]["standalone_normal_termination"]


def test_missing_graph_hydrogens(inputs: tuple) -> None:
    with Chem.SDWriter(str(inputs[1])) as writer:
        writer.write(Chem.MolFromSmiles("O"))
    with pytest.raises(ValueError, match="separate graph atom"):
        _assess(inputs)


def test_invalid_sdf(inputs: tuple) -> None:
    inputs[1].write_text("invalid\n$$$$\n")
    with pytest.raises(ValueError, match="valid explicit-hydrogen"):
        _assess(inputs)
