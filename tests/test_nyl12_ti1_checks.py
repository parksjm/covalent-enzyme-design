"""Synthetic tests for TI1 topology and Amber convergence screening."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import parmed as pmd
import pytest

from nyl12_covalent import ti1_checks


def _model(target: str = "pa6_l2") -> pmd.Structure:
    """Construct the ordered TI1 inventory with a valid screening geometry."""
    structure = pmd.Structure()
    xyz: list[tuple[float, float, float]] = []

    residue_names = {
        928: "THR",
        1024: "THR",
        1120: "THR",
        1216: "THR",
        1312: "LIG",
        1313: "WAT",
    }

    def atom(
        res_idx: int,
        name: str,
        element: int,
        position: tuple[float, float, float],
        charge: float = 0.0,
    ) -> pmd.Atom:
        item = pmd.Atom(name=name, atomic_number=element, charge=charge)
        structure.add_atom(item, residue_names.get(res_idx, "ALA"), res_idx + 1)
        xyz.append(position)
        return item

    def bond(first: pmd.Atom, second: pmd.Atom) -> None:
        structure.bonds.append(pmd.Bond(first, second))

    def fill(start: int, stop: int) -> None:
        """Add one distant carbon to each otherwise empty residue in a range."""
        for index in range(start, stop):
            atom(index, "CA", 6, (100.0 + index, 100.0, 100.0))

    fill(0, 928)

    thr_n = atom(928, "N", 7, (0.0, 0.0, 0.0), -0.3)
    thr_h1 = atom(928, "H1", 1, (-0.4, 0.9, 0.0), 0.3)
    thr_h2 = atom(928, "H2", 1, (-0.4, -0.45, 0.78), 0.3)
    thr_h3 = atom(928, "H3", 1, (-0.4, -0.45, -0.78), 0.3)
    og1 = atom(928, "OG1", 8, (1.4, 0.0, 0.0), -0.6)
    for hydrogen in (thr_h1, thr_h2, thr_h3):
        bond(thr_n, hydrogen)
    fill(929, 1024)

    for index in (1024, 1120, 1216):
        if index == 1120:
            fill(1025, 1120)
        elif index == 1216:
            fill(1121, 1216)
        beta_n = atom(index, "N", 7, (300.0 + index, 0.0, 0.0), -0.4)
        beta_h1 = atom(index, "H1", 1, (300.0 + index, 1.0, 0.0), 0.2)
        beta_h2 = atom(index, "H2", 1, (300.0 + index, -1.0, 0.0), 0.2)
        for hydrogen in (beta_h1, beta_h2):
            bond(beta_n, hydrogen)
    fill(1217, 1312)

    carbon_name, oxygen_name, nitrogen_name, count = ti1_checks.TARGET_ATOMS[target]
    carbon = atom(1312, carbon_name, 6, (2.8, 0.0, 0.0))
    oxygen = atom(1312, oxygen_name, 8, (4.2, 0.0, 0.0), -0.5)
    nitrogen = atom(1312, nitrogen_name, 7, (2.8, -1.4, 0.0), 0.2)
    fourth = atom(1312, "H_C", 1, (2.8, 0.0, -1.4), 0.3)
    for partner in (og1, oxygen, nitrogen, fourth):
        bond(carbon, partner)

    ligand_names = {carbon_name, oxygen_name, nitrogen_name, "H_C"}
    for number in range(count - 4):
        name = f"X{number:03d}"
        while name in ligand_names:
            name += "X"
        ligand_names.add(name)
        atom(1312, name, 6, (200.0 + number, 100.0, 100.0))
    bond(oxygen, carbon)

    water_o = atom(1313, "O", 8, (0.7, 2.8, 0.0), -0.834)
    water_h1 = atom(1313, "H1", 1, (1.65, 2.8, 0.0), 0.417)
    water_h2 = atom(1313, "H2", 1, (0.46, 3.71, 0.0), 0.417)
    bond(water_o, water_h1)
    bond(water_o, water_h2)

    structure.coordinates = np.asarray(xyz, dtype=float)
    return structure


def _patch_model(monkeypatch: pytest.MonkeyPatch, model: pmd.Structure) -> None:
    """Replace AmberParm loading while retaining real ParmEd atom behavior."""
    monkeypatch.setattr(ti1_checks, "AmberParm", lambda *args, **kwargs: model)


def test_geometry_accepts_valid_product_gated_inventory(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both substrate registers pass inventory, bond, and water screens."""
    for target in ti1_checks.TARGET_ATOMS:
        model = _model(target)
        _patch_model(monkeypatch, model)
        result = ti1_checks.check_geometry(
            Path("top.parm7"), Path("coords.rst7"), target
        )
        assert result["accepted"] is True, result
        assert result["inventory_passed"] is True
        assert all(bond["accepted"] for bond in result["bonds"])
        assert result["water_maximum_contact_shortfall_A"] == 0.0


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("ligand_count", "Ligand inventory"),
        ("beta_h3", "neutral NH2"),
        ("ti1_h3_bond", "proton connectivity"),
        ("c_n_bond", "required bond"),
        ("water_charge", "neutral H2O"),
    ],
)
def test_geometry_rejects_invalid_inventory_and_connectivity(
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
    message: str,
) -> None:
    """Changed ligand count, proton states, or required bonds fail explicitly."""
    model = _model()
    xyz = np.asarray(model.coordinates).reshape((-1, 3)).copy()
    if mutation == "ligand_count":
        last = model.residues[1312].atoms[-1]
        model.atoms.remove(last)
        model = SimpleNamespace(
            atoms=model.atoms,
            residues=model.residues,
            coordinates=xyz[:-1],
        )
    elif mutation == "beta_h3":
        extra = pmd.Atom(name="H3", atomic_number=1, charge=0.0)
        model.residues[1024].atoms.append(extra)
        model.atoms.append(extra)
        model = SimpleNamespace(
            atoms=model.atoms,
            residues=model.residues,
            coordinates=np.vstack([xyz, [400, 0, 0]]),
        )
    elif mutation == "ti1_h3_bond":
        next(b for b in model.bonds if "H3" in (b.atom1.name, b.atom2.name)).delete()
    elif mutation == "c_n_bond":
        next(
            b for b in model.bonds if {b.atom1.name, b.atom2.name} == {"C54", "N59"}
        ).delete()
    elif mutation == "water_charge":
        next(a for a in model.residues[1313].atoms if a.name == "O").charge = -0.5
    _patch_model(monkeypatch, model)
    with pytest.raises(ValueError, match=message):
        ti1_checks.check_geometry(Path("top"), Path("coords"), "pa6_l2")


def test_geometry_reports_distant_proton_and_water_clashes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Geometric failures remain visible as rejected report entries."""
    model = _model()
    xyz = np.asarray(model.coordinates).reshape((-1, 3)).copy()
    xyz[next(a.idx for a in model.atoms if a.name == "H3" and a.residue.idx == 928)] = (
        50.0,
        0.0,
        0.0,
    )
    model.coordinates = xyz
    _patch_model(monkeypatch, model)
    result = ti1_checks.check_geometry(Path("top"), Path("coords"), "pa6_l2")
    assert result["accepted"] is False
    assert any(not item["accepted"] for item in result["bonds"])

    model = _model()
    water_o = next(a for a in model.residues[1313].atoms if a.name == "O")
    xyz = np.asarray(model.coordinates).reshape((-1, 3)).copy()
    xyz[water_o.idx] = (0.0, 0.2, 0.0)
    model.coordinates = xyz
    _patch_model(monkeypatch, model)
    result = ti1_checks.check_geometry(Path("top"), Path("coords"), "pa6_l2")
    assert result["accepted"] is False
    assert result["water_maximum_contact_shortfall_A"] > 0.10


@pytest.mark.parametrize("bad_coordinates", ["nan", "wrong_shape"])
def test_geometry_rejects_invalid_coordinates(
    monkeypatch: pytest.MonkeyPatch,
    bad_coordinates: str,
) -> None:
    """Nonfinite coordinates and topology-order shape mismatches fail closed."""
    model = _model()
    if bad_coordinates == "nan":
        xyz = np.asarray(model.coordinates).reshape((-1, 3)).copy()
        xyz[0, 0] = np.nan
        model.coordinates = xyz
    else:
        model = SimpleNamespace(
            atoms=model.atoms,
            residues=model.residues,
            coordinates=np.asarray(model.coordinates).reshape((-1, 3))[:-1],
        )
    _patch_model(monkeypatch, model)
    with pytest.raises(ValueError, match="Coordinates must be finite"):
        ti1_checks.check_geometry(Path("top"), Path("coords"), "pa6_l2")


def _mdout(rms: str = "0.010", gmax: str = "0.050", extra: str = "") -> str:
    """Return a minimal final Amber minimization result block."""
    return (
        "Amber output\nFINAL RESULTS\nNSTEP ENERGY RMS GMAX\n"
        f"  100 -12.5 {rms} {gmax}\n{extra}"
    )


def test_convergence_accepts_finite_final_result(tmp_path: Path) -> None:
    """Finite converged values after FINAL RESULTS pass."""
    output = tmp_path / "min.mdout"
    output.write_text(_mdout())
    result = ti1_checks.check_convergence(output)
    assert result["accepted"] is True
    assert result["final_energy_rms_gmax"] == [-12.5, 0.01, 0.05]
    assert result["final_results_present"] is True


@pytest.mark.parametrize(
    ("contents", "accepted"),
    [
        (_mdout(rms="0.051"), False),
        (_mdout(gmax="-0.01"), False),
        (_mdout(rms="********"), False),
        (_mdout(rms="NaN"), False),
        ("Amber output without final block", False),
        (_mdout(extra="SCF calculation failed\n"), False),
        (_mdout(extra="SANDER BOMB: abnormal termination\n"), False),
        (_mdout(extra="QM error in evaluation\n"), False),
    ],
)
def test_convergence_fails_bad_missing_or_engine_error_results(
    tmp_path: Path,
    contents: str,
    accepted: bool,
) -> None:
    """Convergence fails for thresholds, malformed output, and engine failures."""
    output = tmp_path / "min.mdout"
    output.write_text(contents)
    result = ti1_checks.check_convergence(output)
    assert result["accepted"] is accepted


def test_convergence_uses_last_results_block_and_fortran_exponents(
    tmp_path: Path,
) -> None:
    """The final block controls acceptance and Amber D exponents are parsed."""
    output = tmp_path / "min.mdout"
    output.write_text(_mdout(rms="0.20") + _mdout(rms="1.0D-2", gmax="2.0D-1"))
    result: dict[str, Any] = ti1_checks.check_convergence(output, drms=0.02)
    assert result["accepted"] is True
    assert result["final_energy_rms_gmax"] == [-12.5, 0.01, 0.2]
