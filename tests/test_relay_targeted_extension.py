"""Test isolated continuation lineage and fail-closed native/diagnostic ordering."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from nyl12_covalent import relay_targeted_extension as prep
from nyl12_covalent import relay_targeted_extension_runtime as runtime
from nyl12_covalent.neb_umbrella_campaign import digest, write_json

SOURCE = """audited source
&qmmm
 qm_theory='DFTB3', qmcharge=-2, dftb_telec=0.0,
 qmmask='@1-98', scfconv=1.d-12,
/
&dlfind
 active=':1-20,25',
/
"""


@pytest.fixture
def prepared(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Create four accepted two-block prefixes with immutable original statuses."""
    base = tmp_path / "additional_sampling_20261126"
    base.mkdir()
    topology = base / "topology.prmtop"
    topology.write_text("native topology")
    source_mdin = base / "original_source.in"
    source_mdin.write_text(SOURCE)
    rows = []
    for number, label in enumerate(prep.LABELS):
        folder = base / label
        folder.mkdir()
        (folder / "seed.ncrst").write_text("native initial seed")
        (folder / "window_status.json").write_text("ORIGINAL STATUS")
        stages = []
        for index, name in enumerate(
            (
                "01_equilibration_1ps",
                "02_equilibration_4ps",
                "03_production_5ps",
                "04_production_5ps",
            )
        ):
            stage = folder / name
            stage.mkdir()
            (stage / "restart.ncrst").write_text(name)
            (stage / "run.mdin").write_text("generated native mdin without dlfind")
            (stage / "umbrella.RST").write_text(
                f"exact unchanged native restraints {label}"
            )
            steps = (2000, 8000, 10000, 10000)[index]
            write_json(
                stage / "acceptance.json",
                {
                    "accepted": True,
                    "completed_steps": steps,
                    "every_integration_step_SCC_converged": True,
                    "checkpoint_sha256": digest(stage / "restart.ncrst"),
                },
            )
            write_json(stage / "process.json", {"exit_status": 0})
            write_json(stage / "analysis_outcome.json", {"status": "completed"})
            write_json(stage / "input_identity.json", {})
            stages.append(
                {
                    "label": name,
                    "input": f"{label}/{name}/run.mdin",
                    "steps": steps,
                    "seed": 100 + number * 10 + index,
                    "kind": "production" if index >= 2 else "equil",
                }
            )
        rows.append(
            {
                "label": label,
                "seed": f"{label}/seed.ncrst",
                "stages": stages,
                "bias": {
                    "feature_names": ["q_acyl", "q_act"],
                    "centers_A": [0.0, 0.0],
                    "force_constants": [25.0, 25.0],
                },
            }
        )
    cfg = {
        "windows": rows,
        "temperature_K": 300,
        "topology": str(topology),
        "temperature_gate_K": [150.0, 450.0],
        "system": {
            "active_atoms": 4165,
            "full_atoms": 104922,
            "qm_charge": -2,
            "qm_real_atoms": 98,
            "qm_method": "DFTB3/3ob-3-1",
            "active_mask": ":1-20,25",
            "qm_atom_ids_mask": "@1-98",
        },
        "frozen_files_sha256": {},
        "external_files_sha256": {str(topology): digest(topology)},
        "targets_plan": {
            "source_mdin": str(source_mdin),
            "sander": "/never/run/native",
            "environment": {"OMP_NUM_THREADS": "1"},
        },
    }
    path = base / "prepared.json"
    write_json(path, cfg)
    monkeypatch.setattr(
        prep, "verify_target", lambda *args: json.loads(path.read_text())
    )
    return path


def test_prepare_preserves_prefix_and_real_native_restart_input(prepared: Path) -> None:
    """All sixteen stages retain velocity continuation, native QM and exact RST."""
    original = prepared.read_bytes()
    manifest = prep.prepare_extension(prepared, prepared.parent / "extension_20261127")
    plan = prep.verify_extension(manifest, digest(manifest))
    seeds = []
    for row in plan["cfg"]["windows"]:
        assert len(row["stages"]) == 8
        for spec in row["extension_stages"]:
            stage = (prepared.parent / spec["input"]).parent
            text = (stage / "run.mdin").read_text()
            assert "irest=1, ntx=5" in text and "tempi=" not in text
            assert "nstlim=10000" in text and "dt=0.0005000" in text
            assert "qmcharge=-2" in text and "qmmask='@1-98'" in text
            assert (stage / "umbrella.RST").read_bytes() == (
                prepared.parent / row["label"] / "04_production_5ps/umbrella.RST"
            ).read_bytes()
            seeds.append(spec["seed"])
        assert (
            prepared.parent / row["label"] / "window_status.json"
        ).read_text() == "ORIGINAL STATUS"
    assert len(set(seeds)) == 16 and prepared.read_bytes() == original
    with pytest.raises(ValueError, match="hash changed"):
        prep.verify_extension(manifest, "wrong")
    with pytest.raises(ValueError, match="logical child"):
        prep.prepare_extension(prepared, manifest.parent)
    restart = prepared.parent / prep.LABELS[0] / "04_production_5ps/restart.ncrst"
    restart.write_text("altered velocities")
    with pytest.raises(ValueError, match="Pinned file"):
        prep.verify_extension(manifest, digest(manifest))


@pytest.mark.parametrize(
    "fault,match",
    [
        ("receipt", "Unaccepted"),
        ("checkpoint", "checkpoint"),
        ("diagnostics", "diagnostics"),
        ("SCC", "SCC"),
    ],
)
def test_prefix_failures_stop_before_mutation(
    prepared: Path, fault: str, match: str
) -> None:
    """Never continue a changed or scientifically unaccepted source prefix."""
    stage = prepared.parent / prep.LABELS[0] / "04_production_5ps"
    if fault == "checkpoint":
        (stage / "restart.ncrst").write_text("changed")
    else:
        path = stage / (
            "analysis_outcome.json" if fault == "diagnostics" else "acceptance.json"
        )
        m = json.loads(path.read_text())
        m[
            {
                "receipt": "accepted",
                "diagnostics": "status",
                "SCC": "every_integration_step_SCC_converged",
            }[fault]
        ] = "failed" if fault == "diagnostics" else False
        write_json(path, m)
    with pytest.raises(ValueError, match=match):
        prep.prepare_extension(prepared, prepared.parent / "extension_20261127")
    assert not (prepared.parent / "extension_20261127").exists()


def fake_runtime(monkeypatch: pytest.MonkeyPatch) -> Mock:
    """Replace native calculations while retaining actual filesystem lineage checks."""
    monkeypatch.setattr(runtime.parmed, "load_file", lambda p: object())
    monkeypatch.setattr(
        runtime.Rst7,
        "open",
        lambda p: SimpleNamespace(coordinates=np.zeros((104922, 3)), box=np.ones(6)),
    )
    monkeypatch.setattr(
        runtime,
        "AmberMask",
        lambda top, mask: SimpleNamespace(
            Selected=lambda: range(98 if mask.startswith("@") else 4165)
        ),
    )
    monkeypatch.setattr(
        runtime,
        "_checkpoint",
        lambda p, *args: ({"checkpoint_sha256": digest(p)}, np.zeros((104922, 3))),
    )
    monkeypatch.setattr(
        runtime, "native_acceptance", Mock(return_value={"accepted": True})
    )
    monkeypatch.setattr(runtime, "ensure_diagnostics", Mock(return_value=True))

    def native(stage: Path, *args: object) -> Path:
        path = stage / "restart.ncrst"
        path.write_text("continued native coordinates and velocities")
        return path

    mock = Mock(side_effect=native)
    monkeypatch.setattr(runtime, "_run_stage", mock)
    return mock


def test_four_stages_are_one_label_continuation(
    prepared: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every stage consumes its own predecessor and completes full diagnostics first."""
    manifest = prep.prepare_extension(prepared, prepared.parent / "extension_20261127")
    native = fake_runtime(monkeypatch)
    label = prep.LABELS[0]
    runtime.run_extension(manifest, label, digest(manifest))
    calls = native.call_args_list
    assert (
        len(calls) == 4
        and calls[0].args[3]
        == prepared.parent / label / "04_production_5ps/restart.ncrst"
    )
    assert [c.args[3].parent.name for c in calls[1:]] == [
        "05_production_5ps",
        "06_production_5ps",
        "07_production_5ps",
    ]
    result = json.loads((manifest.parent / "status" / f"{label}.json").read_text())
    assert result["status"] == "completed" and result["production_ps_total"] == 30
    assert result["native_and_diagnostics_accepted"] is True
    assert runtime.ensure_diagnostics.call_count == 4
    with pytest.raises(FileExistsError):
        runtime.run_extension(manifest, label, digest(manifest))
    with pytest.raises(ValueError, match="exactly one"):
        runtime.run_extension(manifest, "undeclared", digest(manifest))


def test_diagnostic_failure_stops_next_native_stage(
    prepared: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Native completion alone cannot promote a block with failed bias/proton analysis."""
    manifest = prep.prepare_extension(prepared, prepared.parent / "extension_20261127")
    native = fake_runtime(monkeypatch)
    runtime.ensure_diagnostics.return_value = False
    label = prep.LABELS[0]
    with pytest.raises(RuntimeError, match="bias/proton"):
        runtime.run_extension(manifest, label, digest(manifest))
    assert native.call_count == 1
    result = json.loads((manifest.parent / "status" / f"{label}.json").read_text())
    assert (
        result["status"] == "failed" and result["stages"][0]["native_accepted"] is True
    )
