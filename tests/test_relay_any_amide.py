"""Check coherent chemistry, predictor atom maps and nested-site accounting."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from nyl12_covalent.relay_any_amide import _tier, evaluate_model, load_graph

SOURCE = (
    Path(__file__).resolve().parents[1] / "design/nyl12_ligandmpnn_pilot/"
    "boltz_validation/configuration/whole_l4_source.mol2"
)


@pytest.fixture()
def model() -> tuple:
    """Make four separated occupied sites using the actual named ligand graph."""
    graph = load_graph(SOURCE)
    mol = Chem.MolFromMol2File(str(SOURCE), removeHs=False, sanitize=True)
    conformer = mol.GetConformer()
    xyz = {
        a.GetProp("_TriposAtomName"): np.array(conformer.GetAtomPosition(a.GetIdx()))
        for a in mol.GetAtoms()
        if a.GetAtomicNum() != 1
    }
    atoms, maps = {}, {}
    for i, (alpha, beta, ligand) in enumerate(zip("ABCD", "EFGH", "IJKL", strict=True)):
        offset = np.array([100.0 * i, 0, 0])
        maps[ligand] = {}
        for j, name in enumerate(graph.elements):
            model_name = f"X{j}"
            maps[ligand][model_name] = name
            atoms[ligand, 1, model_name] = graph.elements[name], xyz[name] + offset
        c, o, n = (xyz[k] + offset for k in ("C54", "O49", "N59"))
        atoms[alpha, 113, "N"] = "N", o + [0, 0, 3]
        atoms[alpha, 186, "ND2"] = "N", o + [0, 2.8, 0]
        atoms[beta, 1, "OG1"] = "O", c + [0, 0, 3]
        atoms[beta, 1, "N"] = "N", n + [0, 0, 3]
        atoms[beta, 1, "CB"] = "C", c + [0, 1, 4]
        for residue in (40, 42):
            for oxygen in (1, 2):
                atoms[beta, residue, f"OD{oxygen}"] = "O", n + [residue / 10, oxygen, 1]
    return graph, atoms, maps


def test_all_three_amides_use_carbonyl_side_carbon() -> None:
    graph = load_graph(SOURCE)
    assert len(graph.elements) == 33
    assert {
        (a.carbon, a.oxygen, a.nitrogen, a.adjacent_carbon) for a in graph.amides
    } == {
        ("C52", "O47", "N58", "C60"),
        ("C54", "O49", "N59", "C62"),
        ("C55", "O50", "N57", "C63"),
    }
    assert all(
        a.carbon != "C53" for a in graph.amides
    )  # Terminal carboxyl is not amide.


def test_exhaustive_rows_and_one_to_one_assignment(model: tuple) -> None:
    graph, atoms, maps = model
    result = evaluate_model(atoms, maps, graph)
    rows, assigned = result["all_amide_site_rows"], result["assigned"]
    assert len(rows) == 48
    assert len(assigned) == 4
    assert {r["ligand"] for r in assigned} == set("IJKL")
    assert {r["site"] for r in assigned} == {"AE", "BF", "CG", "DH"}
    for row in rows:
        assert row["amide_geometry_valid"]
        assert -180 <= row["herschlag_phi"] <= 180
        assert np.isfinite(row["nterm_to_leaving_n"])
    assert result["native_thr_hard_gate"] is False


def test_missing_wrong_or_nonfinite_mapped_atoms_fail(model: tuple) -> None:
    graph, atoms, maps = model
    key = ("E", 1, "OG1")
    missing = dict(atoms)
    del missing[key]
    with pytest.raises(ValueError, match="Missing mapped atom"):
        evaluate_model(missing, maps, graph)
    for value, message in [
        (("C", [1, 2, 3]), "Wrong element"),
        (("O", [1, 2]), "Wrong element"),
        (("O", [1, 2, float("nan")]), "Nonfinite"),
    ]:
        bad = {**atoms, key: value}
        with pytest.raises(ValueError, match=message):
            evaluate_model(bad, maps, graph)


def test_nonbijective_and_incomplete_ligand_maps_fail(model: tuple) -> None:
    graph, atoms, maps = model
    for key, replacement in [("X0", "bogus"), ("X1", maps["I"]["X0"])]:
        bad = {**maps, "I": {**maps["I"], key: replacement}}
        with pytest.raises(ValueError, match="Nonbijective"):
            evaluate_model(atoms, bad, graph)
    with pytest.raises(ValueError, match="four L4"):
        evaluate_model(atoms, {"I": maps["I"]}, graph)
    with pytest.raises(ValueError, match="three internal"):
        evaluate_model(atoms, maps, replace(graph, amides=graph.amides[:1]))


def test_bad_bond_geometry_is_not_a_passing_tier(model: tuple) -> None:
    graph, atoms, maps = model
    oxygen = next(k for k, v in maps["I"].items() if v == "O49")
    key = ("I", 1, oxygen)
    bad = {**atoms, key: ("O", atoms[key][1] + [0, 0, 4])}
    rows = evaluate_model(bad, maps, graph)["all_amide_site_rows"]
    affected = [r for r in rows if r["ligand"] == "I" and r["carbonyl_c"] == "C54"]
    assert all(
        not r["amide_geometry_valid"] and r["tier"] == "reject" for r in affected
    )


def test_bad_graph_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(Chem, "MolFromMol2File", lambda *a, **kw: None)
    with pytest.raises(ValueError, match="Invalid MOL2"):
        load_graph(SOURCE)


@pytest.mark.parametrize(
    "d1,label", [(3.0, "A"), (3.4, "B"), (3.7, "C"), (4, "reject")]
)
@pytest.mark.parametrize("phi", [-90.0, 90.0])
def test_legacy_tiers_accept_either_signed_orientation(
    d1: float, label: str, phi: float
) -> None:
    row = {
        "amide_geometry_valid": True,
        "d1": d1,
        "d2": 3,
        "d3": 3,
        "a_bd": 105,
        "herschlag_phi": phi,
    }
    assert _tier(row) == label
