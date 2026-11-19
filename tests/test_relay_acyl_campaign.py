"""Chemistry and conditional queue tests without optimization or engine calls."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from scripts import run_relay_acyl_campaign as campaign
from nyl12_covalent.reaction_bias import HarmonicTarget, evaluate_bias


@pytest.fixture
def model() -> SimpleNamespace:
    mapping = dict(
        thr_n=1,
        thr_h1=2,
        thr_h2=3,
        thr_cb=4,
        thr_og=5,
        proton=6,
        carbon=7,
        oxygen=8,
        leaving_n=9,
        leaving_h=10,
        water_o=11,
        water_h1=12,
        water_h2=13,
        ash_oxygen=14579,
        ash_proton=15,
        asp_o1=16,
        asp_o2=17,
        phe_n=18,
        asn_nd2=19,
    )
    coordinates = np.random.default_rng(10).normal(size=(24, 3))
    hydrogens = np.array([1, 2, 5, 9, 11, 12, 14, 19, 20, 21, 22, 23])
    atoms = [SimpleNamespace(atomic_number=6, bond_partners=[]) for _ in range(24)]
    atoms[18] = SimpleNamespace(
        atomic_number=7, bond_partners=[SimpleNamespace(atomic_number=1)]
    )
    return SimpleNamespace(
        atom_map=mapping,
        hydrogens=hydrogens,
        coordinates=coordinates,
        qm=np.array([18]),
        topology=SimpleNamespace(atoms=atoms),
    )


def product_state(model: SimpleNamespace) -> dict:
    """Return an unambiguous balanced product inventory with 54 audited H."""
    m = model.atom_map
    owners = [
        (m["thr_h1"], m["thr_n"]),
        (m["thr_h2"], m["thr_n"]),
        (m["leaving_h"], m["leaving_n"]),
        (m["proton"], m["leaving_n"]),
        (m["water_h1"], m["water_o"]),
        (m["water_h2"], m["water_o"]),
        (m["ash_proton"], m["ash_oxygen"]),
    ]
    owners += [(100 + i, 1000 + i) for i in range(54 - len(owners))]
    owners[-1] = (owners[-1][0], 19)
    return {
        "rOG1C": 1.4,
        "rCN": 3.0,
        "rCO": 1.23,
        "proton_owners": [
            {
                "hydrogen_1based": h,
                "nearest_qm_heavy_1based": parent,
                "nearest_distance_A": 1.0,
                "classification": "covalent",
            }
            for h, parent in owners
        ],
    }


def test_classify_balanced_product_and_reactant(
    model: SimpleNamespace, monkeypatch: pytest.MonkeyPatch
) -> None:
    state = product_state(model)
    monkeypatch.setattr(campaign, "state_record", lambda *_: deepcopy(state))
    assert campaign.classify(model, model.coordinates)["product_geometry_retained"]
    assert not campaign.classify(model, model.coordinates)["reactant_geometry"]
    state["rOG1C"], state["rCN"] = 3.0, 1.35
    transferred = next(
        row
        for row in state["proton_owners"]
        if row["hydrogen_1based"] == model.atom_map["proton"]
    )
    transferred["nearest_qm_heavy_1based"] = model.atom_map["thr_og"]
    assert campaign.classify(model, model.coordinates)["reactant_geometry"]


@pytest.mark.parametrize(
    "defect",
    [
        "ester",
        "scission",
        "carbonyl",
        "ambiguous",
        "water",
        "asp",
        "ash",
        "unassigned",
        "spectator",
    ],
)
def test_classifier_rejects_bad_product(
    model: SimpleNamespace, monkeypatch: pytest.MonkeyPatch, defect: str
) -> None:
    state = product_state(model)
    if defect == "ester":
        state["rOG1C"] = 2.5
    elif defect == "scission":
        state["rCN"] = 1.5
    elif defect == "carbonyl":
        state["rCO"] = 1.5
    elif defect == "ambiguous":
        state["proton_owners"][-1]["classification"] = "shared_ambiguous"
    elif defect == "unassigned":
        state["proton_owners"][-1]["classification"] = "unassigned"
        state["proton_owners"][-1]["nearest_distance_A"] = 3.0
    elif defect == "spectator":
        state["proton_owners"][-1]["nearest_qm_heavy_1based"] = 1001
    else:
        owner = {"water": "water_o", "asp": "asp_o1", "ash": "ash_oxygen"}[defect]
        state["proton_owners"][-1]["nearest_qm_heavy_1based"] = model.atom_map[owner]
    monkeypatch.setattr(campaign, "state_record", lambda *_: state)
    assert not campaign.classify(model, model.coordinates)["product_geometry_retained"]


def test_asynchronous_coordinate_formula_and_hydrogen_permutation(
    model: SimpleNamespace,
) -> None:
    cv = campaign.asynchronous_coordinate(model)
    xyz = model.coordinates
    value = evaluate_bias(xyz, [HarmonicTarget(cv, 0, 0)]).values[0]
    m = model.atom_map
    from nyl12_covalent.proton_coordination import smooth_coordination

    expected = (
        np.linalg.norm(xyz[m["thr_og"] - 1] - xyz[m["carbon"] - 1])
        - np.linalg.norm(xyz[m["carbon"] - 1] - xyz[m["leaving_n"] - 1])
        + smooth_coordination(xyz, m["thr_og"] - 1, model.hydrogens)[0]
        - smooth_coordination(xyz, m["leaving_n"] - 1, model.hydrogens)[0]
    )
    assert value == pytest.approx(expected)
    swapped = xyz.copy()
    swapped[model.hydrogens] = swapped[model.hydrogens[::-1]]
    assert evaluate_bias(swapped, [HarmonicTarget(cv, 0, 0)]).values[
        0
    ] == pytest.approx(value)


@pytest.mark.parametrize(
    "stop", ["seed", "unbiased_force", "unbiased_chemistry", "window", "complete"]
)
def test_conditional_campaign_queue(
    tmp_path: Path, model: SimpleNamespace, monkeypatch: pytest.MonkeyPatch, stop: str
) -> None:
    calls = []
    assessments = []
    monkeypatch.setattr(campaign, "load_case", lambda _: model)
    monkeypatch.setattr(
        campaign.Rst7, "open", lambda _: SimpleNamespace(coordinates=model.coordinates)
    )

    def classify(*_: object) -> dict:
        assessments.append(True)
        return {
            "product_geometry_retained": stop != "unbiased_chemistry",
            "reactant_geometry": True,
        }

    def run_stage(
        case_path: Path, start: Path, output: Path, targets: list, *, max_steps: int
    ) -> dict:
        output.mkdir(parents=True)
        calls.append((output, targets, max_steps))
        force = not (
            (stop == "seed" and len(calls) == 1)
            or (stop == "unbiased_force" and len(calls) == 2)
            or (stop == "window" and len(calls) == 3)
        )
        return {"endpath": str(output / "optimized.ncrst"), "forceconv": force}

    monkeypatch.setattr(campaign, "classify", classify)
    monkeypatch.setattr(campaign, "run_stage", run_stage)
    campaign.run(tmp_path / "case.json")
    expected = {
        "seed": 1,
        "unbiased_force": 2,
        "unbiased_chemistry": 2,
        "window": 3,
        "complete": 39,
    }[stop]
    assert len(calls) == expected
    if len(calls) >= 2:
        assert calls[1][1] == []
    if stop in ("seed", "unbiased_force", "unbiased_chemistry"):
        assert not (tmp_path / "03_reverse_path").exists()
    if stop == "complete":
        assert calls[-1][1] == []
        assert all(len(terms) == 1 for _, terms, _ in calls[2:-1])
    assert (tmp_path / "campaign_status.json").is_file()


@pytest.mark.parametrize("severe", [0, 1])
def test_preparation_uses_winning_recipe_and_only_publishes_zero_severe(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, severe: int
) -> None:
    import json
    from scripts import prepare_relay_acyl_campaign as preparation

    base = tmp_path / "base"
    source = base / "source"
    folder = base / "campaign"
    source.mkdir(parents=True)
    folder.mkdir()
    (base / "00_system").mkdir()
    (base / "00_system/ash274_solvated.prmtop").write_text("synthetic topology")
    (source / "optimized.ncrst").write_text("synthetic accepted restart")
    (source / "acceptance.json").write_text(
        json.dumps({"status": "accepted", "acceptance": {"accepted": True}})
    )
    (source / "minimize.mdin").write_text(
        "&cntrl maxcyc=2000 /\n&qmmm qmmask='QM',qm_theory='DFTB3', /\n&dlfind active='ACTIVE', /\n"
    )
    monkeypatch.setattr(preparation, "BASE", base)
    monkeypatch.setattr(preparation, "SOURCE", source)
    monkeypatch.setattr(preparation, "CAMPAIGN", folder)
    monkeypatch.setattr(preparation.parmed, "load_file", lambda _: SimpleNamespace())
    monkeypatch.setattr(
        preparation.Rst7,
        "open",
        lambda _: SimpleNamespace(coordinates=np.zeros((6, 3))),
    )
    monkeypatch.setattr(
        preparation,
        "AmberMask",
        lambda _, text: SimpleNamespace(
            Selected=lambda: iter(range(102 if text == "QM" else 4662))
        ),
    )
    monkeypatch.setattr(
        preparation,
        "load_case",
        lambda _: SimpleNamespace(active=np.ones(4662, dtype=bool)),
    )
    calls = []
    guess = np.ones((6, 3))

    def search(*args: object, **kwargs: object) -> tuple:
        calls.append("leaving")
        return guess, {"base": True}

    def refine(*args: object, **kwargs: object) -> tuple:
        calls.append("torsion")
        np.testing.assert_array_equal(args[3], guess)
        assert kwargs["max_heavy_displacement_A"] == 3
        assert kwargs["max_hydrogen_displacement_A"] == 5
        return 2 * guess, {
            "search": {"recipe": "285 then 24"},
            "severe_contact_count": severe,
            "close_contacts": [],
            "ester_angle_degrees": 83.968,
        }

    monkeypatch.setattr(preparation, "search_balanced_acyl_seed", search)
    monkeypatch.setattr(preparation, "refine_acyl_torsion_seed", refine)
    monkeypatch.setattr(
        preparation,
        "write_restart",
        lambda path, *_: path.write_text("synthetic selected seed"),
    )
    if severe:
        with pytest.raises(ValueError, match="severe"):
            preparation.prepare()
        assert not (folder / "01_endpoint_inputs/acyl_guess.ncrst").exists()
    else:
        assert preparation.prepare() == folder / "case.json"
        assert (folder / "01_endpoint_inputs/acyl_guess.ncrst").is_file()
        manifest = json.loads((folder / "01_endpoint_inputs/manifest.json").read_text())
        assert manifest["seed_is_optimized"] is False
        assert manifest["seed_ester_angle_degrees"] == 83.968
        assert any(
            path.endswith("relay_seed_refinement.py")
            for path in manifest["input_sha256"]
        )
    assert calls == ["leaving", "torsion"]
