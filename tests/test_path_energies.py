"""Check native energy selection, bias subtraction and explicit model limits."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent.path_energies import relay_endpoint_energies


def _window(run: Path, label: str, stage: int, total: float) -> dict:
    folder = run / label
    folder.mkdir()
    (folder / "minimized.mdout").write_text(f"Etot = -999\nEtot = {total}\n")
    target = 0.2
    (folder / "process.json").write_text(
        json.dumps(
            {
                "stage": stage,
                "attack_target_A": 1.5,
                "proton_difference_target_A": target,
            }
        )
    )
    (folder / "acceptance.json").write_text(
        json.dumps(
            {"distances_A": {"attack": 1.6, "OG_H": 1.3, "Asp_H": 1.2, "N_HG": 1.0}}
        )
    )
    blocks = []
    for t in (1.5, target):
        blocks.append(
            f"&rst r1={t - 2}, r2={t}, r3={t}, r4={t + 2}, rk2=100, rk3=100, /"
        )
    (folder / "restraints.RST").write_text("\n".join(blocks))
    return {"window": label, "accepted": stage == 1}


def test_native_total_bias_and_reference_are_separate(tmp_path: Path) -> None:
    endpoints = [
        _window(tmp_path, "stage1_00", 1, -100),
        _window(tmp_path, "stage2_00", 2, -95),
    ]
    rows = relay_endpoint_energies(tmp_path, endpoints)
    assert rows[0]["total_kcal_mol"] == -100
    assert rows[0]["driving_bias_kcal_mol"] == pytest.approx(2)
    assert rows[1]["driving_bias_kcal_mol"] == pytest.approx(1)
    assert rows[1]["relative_total_kcal_mol"] == 5
    assert rows[1]["relative_minus_driving_kcal_mol"] == pytest.approx(6)
    assert not rows[1]["accepted"]
    assert len(rows[1]["source_sha256"]) == 4


@pytest.mark.parametrize(
    "invalid", ["energy", "stage", "count", "force", "region", "target"]
)
def test_invalid_or_changed_energy_model_fails(tmp_path: Path, invalid: str) -> None:
    endpoint = _window(tmp_path, "stage1_00", 1, -100)
    folder = tmp_path / "stage1_00"
    if invalid == "energy":
        (folder / "minimized.mdout").write_text("not an energy")
    elif invalid == "stage":
        p = folder / "process.json"
        data = json.loads(p.read_text())
        data["stage"] = 3
        p.write_text(json.dumps(data))
    else:
        p = folder / "restraints.RST"
        text = p.read_text()
        if invalid == "count":
            text = ""
        elif invalid == "force":
            text = text.replace("rk2=100", "rk2=50")
        elif invalid == "region":
            text = text.replace("r4=3.5", "r4=1.55")
        else:
            text = text.replace("r2=1.5", "r2=1.4")
        p.write_text(text)
    with pytest.raises(ValueError):
        relay_endpoint_energies(tmp_path, [endpoint])


def test_no_endpoints_fails(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        relay_endpoint_energies(tmp_path, [])
