"""Check native energy selection, bias subtraction and explicit model limits."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent import path_energies
from nyl12_covalent.path_energies import relay_endpoint_energies
from nyl12_covalent.path_energies import (
    driven_path_energies,
    write_driven_path_energy_data,
)

DISTANCE_KEYS = {
    "attack": "attack",
    "donor_h": "OG_H",
    "relay_h": "Asp_H",
    "acceptor_h": "N_HG",
}


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
                **({"relay_retention_target_A": 0.1} if stage == 2 else {}),
            }
        )
    )
    (folder / "acceptance.json").write_text(
        json.dumps(
            {"distances_A": {"attack": 1.6, "OG_H": 1.3, "Asp_H": 1.2, "N_HG": 1.0}}
        )
    )
    blocks = []
    targets = [1.5, target]
    if stage == 2:
        targets.append(0.1)
    for t in targets:
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
    rows = relay_endpoint_energies(tmp_path, endpoints, DISTANCE_KEYS)
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
        relay_endpoint_energies(tmp_path, [endpoint], DISTANCE_KEYS)


def test_no_endpoints_fails(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        relay_endpoint_energies(tmp_path, [], DISTANCE_KEYS)


def test_write_energy_artifacts_uses_explicit_plotter(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    endpoint = _window(tmp_path, "stage1_00", 1, -100)
    output = tmp_path / "output"
    output.mkdir()
    commands: list[list[str]] = []

    def fake_run(command: list[str], **kwargs: object) -> None:
        commands.append(command)
        assert kwargs == {"env": {"MPLBACKEND": "Agg"}, "check": True}
        Path(command[-1]).write_bytes(b"PNG")

    monkeypatch.setattr(path_energies.subprocess, "run", fake_run)
    path_energies.write_energy_artifacts(
        tmp_path,
        [endpoint],
        output,
        "plot-python",
        Path("plot.py"),
        {"MPLBACKEND": "Agg"},
        DISTANCE_KEYS,
    )
    assert commands == [
        [
            "plot-python",
            "plot.py",
            str(output / "reaction_energies.json"),
            str(output / "reaction_energies.png"),
        ]
    ]
    assert (output / "reaction_energies.csv").read_text().startswith("window,")
    assert (output / "reaction_energies.png").read_bytes() == b"PNG"


def test_incomplete_distance_mapping_fails(tmp_path: Path) -> None:
    endpoint = _window(tmp_path, "stage1_00", 1, -100)
    with pytest.raises(ValueError, match="mapping"):
        relay_endpoint_energies(tmp_path, [endpoint], {"attack": "attack"})


def _driven_window(run: Path, name: str, total: float, restraint: float) -> None:
    folder = run / name
    folder.mkdir()
    (folder / "minimized.mdout").write_text(
        f"Etot = -1\nRESTRAINT = 9\nEtot = {total}\nRESTRAINT = {restraint}\n"
    )
    (folder / "acceptance.json").write_text(json.dumps({"native_converged": True}))
    (folder / "restraints.RST").write_text("&rst /\n")


def test_driven_path_uses_reported_final_restraint_energy(tmp_path: Path) -> None:
    _driven_window(tmp_path, "window_01", -100.0, 2.0)
    _driven_window(tmp_path, "window_02", -95.0, 1.0)
    rows = driven_path_energies(tmp_path, ["window_01", "window_02"])
    assert rows[0]["relative_total_kcal_mol"] == 0.0
    assert rows[1]["relative_total_kcal_mol"] == 5.0
    assert rows[1]["relative_minus_driving_kcal_mol"] == 6.0
    assert len(rows[0]["source_sha256"]) == 3
    output = tmp_path / "plots"
    output.mkdir()
    json_path, csv_path = write_driven_path_energy_data(rows, output)
    assert json.loads(json_path.read_text())[1]["window"] == "window_02"
    assert csv_path.read_text().startswith("window,")


def test_driven_path_rejects_unconverged_window(tmp_path: Path) -> None:
    _driven_window(tmp_path, "window_01", -100.0, 2.0)
    (tmp_path / "window_01" / "acceptance.json").write_text(
        json.dumps({"native_converged": False})
    )
    with pytest.raises(ValueError, match="Unconverged"):
        driven_path_energies(tmp_path, ["window_01"])


def test_guarded_continuation_subtracts_both_biases_and_exports_columns(
    tmp_path: Path,
) -> None:
    from nyl12_covalent.direct_acyl_drive import direct_restraint

    _driven_window(tmp_path, "original", -100.0, 2.0)
    _driven_window(tmp_path, "guarded", -95.0, 3.0)
    folder = tmp_path / "guarded"
    (folder / "restraints.RST").write_text(
        direct_restraint(0.3, substrate_proton_guard=True)
    )
    (folder / "acceptance.json").write_text(
        json.dumps(
            {
                "native_converged": True,
                "geometry_A": {"substrate_N_HG1_A": 1.4},
            }
        )
    )
    rows = driven_path_energies(tmp_path, ["original", "guarded"])
    assert rows[1]["substrate_proton_guard_bias_kcal_mol"] == pytest.approx(1.125)
    assert rows[1]["composite_bias_kcal_mol"] == pytest.approx(1.875)
    assert rows[1]["relative_minus_driving_kcal_mol"] == pytest.approx(4.0)
    _, csv_path = write_driven_path_energy_data(rows, tmp_path)
    assert (
        "substrate_proton_guard_bias_kcal_mol" in csv_path.read_text().splitlines()[0]
    )
    (folder / "acceptance.json").write_text(json.dumps({"native_converged": True}))
    with pytest.raises(ValueError, match="guard energy"):
        driven_path_energies(tmp_path, ["guarded"])
