"""Check external-charge conversion and reject incomplete reference gradients."""

import importlib
from pathlib import Path
import subprocess

import pytest


@pytest.fixture
def modules(monkeypatch):
    folder = Path(__file__).resolve().parents[1] / (
        "qmmm/pa6_flipped_model30/07_ash274_comparison/"
        "08_activation_attack_grid/reference_checks"
    )
    monkeypatch.syspath_prepend(str(folder))
    return (
        importlib.import_module("reference_embedding"),
        importlib.import_module("reference_gaussian"),
    )


def deck(path: Path) -> Path:
    """Create a small-PC deck with the required real/link atom inventory."""
    path.write_text(
        "%chk=gau_job.chk\n#P wB97XD/def2SVP Force Charge Prop=(Field,Read)\n\n"
        "test\n\n-1 1\n"
        + "H 0.0 0.0 0.0\n" * 110
        + "\n1.0 2.0 3.0 -0.3\n\n1.0 2.0 3.0\n\n"
    )
    return path


def test_pc_units_and_field_query_removal(modules, tmp_path):
    embedding, _ = modules
    source = deck(tmp_path / "source.inp")
    direct = embedding.write_direct_gaussian(source, tmp_path)
    parsed = embedding.parse_external_input(direct)
    assert parsed["point_charges"] == [(1.0, 2.0, 3.0, -0.3)]
    assert "Prop=(Field,Read)" not in direct.read_text()
    assert all(not s.strip() for s in parsed["lines"][parsed["pc_section_end"] :])
    embedding.write_xtb_inputs(source, tmp_path)
    values = list(
        map(float, (tmp_path / "point_charges.pc").read_text().splitlines()[1].split())
    )
    assert values[:4] == pytest.approx(
        [-0.3, 1.88972612546, 3.77945225092, 5.66917837637]
    )


@pytest.mark.parametrize(
    "atoms, count, accepted", [(110, 330, True), (0, 0, False), (110, 329, False)]
)
def test_gaussian_gradient_inventory(
    modules, tmp_path, monkeypatch, atoms, count, accepted
):
    _, gaussian = modules
    source = deck(tmp_path / "source.inp")
    folder = tmp_path / "calculation"

    def fake_run(command, **kwargs):
        assert kwargs["cwd"] == folder
        if command[0].endswith("g16"):
            kwargs["stdout"].write(
                "SCF Done: E(RwB97XD) = -100.0 A.U.\n"
                "Self energy of the charges = -2.0 a.u.\n"
                "Normal termination of Gaussian\n"
            )
        else:
            (folder / "direct_force.fchk").write_text(
                f"Number of atoms I {atoms}\nCartesian Gradient R N= {count}\n"
                + "0.1 " * count
            )
        return subprocess.CompletedProcess(command, 0, stdout=b"")

    monkeypatch.setattr(gaussian.subprocess, "run", fake_run)
    monkeypatch.setattr(gaussian, "journal", lambda *args: None)
    if accepted:
        result = gaussian.run_direct_gaussian(source, folder, source, {})
        assert result["corrected_gaussian_energy_kcal_mol"] == pytest.approx(
            -98 * 627.5094740631
        )
    else:
        with pytest.raises(ValueError, match="110 atoms and330"):
            gaussian.run_direct_gaussian(source, folder, source, {})
