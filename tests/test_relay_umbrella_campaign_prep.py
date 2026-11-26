"""Artifact-backed checks for immutable relay pilot input preparation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from nyl12_covalent import relay_umbrella_campaign_prep as campaign
from nyl12_covalent.relay_umbrella_prepare import MOBILE_MASK, QM_MASK


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))


def _source_tree(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    base = tmp_path / "relay"
    inputs = base / "model/inputs"
    inputs.mkdir(parents=True)
    (inputs / "system.prmtop").write_text("topology\n")
    (inputs / "dimer.mdin").write_text(
        """source mdin
&cntrl
 imin=1,
/
&qmmm
 qm_theory='DFTB3', qmcharge=-2, dftb_telec=0.0,
 qmmask='"""
        + QM_MASK
        + """', scfconv=1.d-10, qmmm_int=1,
/
&dlfind
 active='"""
        + MOBILE_MASK
        + """', neb='FROZEN',
/
"""
    )
    (inputs / "optimized.ncrst").write_text("source restart\n")
    (inputs / "scientific_acceptance.json").write_text("{}\n")
    (base / "config.json").write_text("{}\n")
    model = base / "model"
    _write_json(model / "manifest.json", {})
    _write_json(model / "release_acceptance.json", {"accepted": True, "force_converged": True})
    topology_hash = _sha(inputs / "system.prmtop")
    _write_json(
        model / "source_identity.json",
        {
            "full_atom_count": 104922,
            "active_atom_count": 4165,
            "qm_atom_count": 98,
            "qm_mm_link_count": 8,
            "qm_charge": -2,
            "topology_sha256": topology_hash,
        },
    )
    fixed = {
        "full_atoms": 104922,
        "active_atoms": 4165,
        "active_mask": MOBILE_MASK[1:],
        "qm_atoms": 98,
        "qm_mask": QM_MASK,
        "qm_mm_link_atoms": 8,
        "qm_charge": -2,
        "qm_method": "DFTB3, 3ob-3-1 parameters; 0 K electronic temperature",
        "protons_tracked": [13906, 13907, 13917, 19937],
    }
    seed_path = inputs / "optimized.ncrst"
    seed_sha = _sha(seed_path)
    for branch in ("01_activation_reverse", "02_acylation_forward"):
        plan = {
            "fixed_scientific_inputs": fixed,
            "source_restart": str(seed_path.resolve()),
            "source_restart_sha256": seed_sha,
            "topology": str((inputs / "system.prmtop").resolve()),
            "topology_sha256": topology_hash,
        }
        _write_json(base / branch / "exploratory_plan_20261124.json", plan)
        evidence = base / branch / "01_coordinate_drive/window_001/acceptance.json"
        restart = evidence.parent / "optimized.ncrst"
        restart.parent.mkdir(parents=True, exist_ok=True)
        restart.write_text(f"accepted {branch}\n")
        _write_json(
            evidence,
            {
                "restart": str(restart.resolve()),
                "native_converged": True,
                "electronic_history": {"passed": True},
                "geometry": {
                    "distances_A": {
                        "C_NL": 1.5,
                        "C_OG": 1.75,
                        "HG_OG": 2.2,
                        "HG_N": 1.1,
                        "HG_NL": 2.0,
                    }
                },
            },
        )

    parameter_release = tmp_path / "parameter_release"
    parameter_copy = parameter_release / "parameters_3ob"
    parameter_copy.mkdir(parents=True)
    amberhome = tmp_path / "amber"
    installed_parameters = amberhome / "dat/slko/3ob-3-1"
    installed_parameters.mkdir(parents=True)
    parameter_hashes = {}
    for first in "HCNO":
        for second in "HCNO":
            name = f"{first}-{second}.skf"
            content = f"parameter {name}\n"
            (parameter_copy / name).write_text(content)
            (installed_parameters / name).write_text(content)
            parameter_hashes[name] = _sha(parameter_copy / name)
    _write_json(parameter_release / "prepared.json", {"parameter_sha256": parameter_hashes})
    sander = amberhome / "bin/sander"
    sander.parent.mkdir(parents=True)
    sander.write_text("validated test runtime placeholder\n")
    sander.chmod(0o755)

    for relative in campaign.SOURCE_SNAPSHOT_FILES:
        source = tmp_path / "repository" / relative
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text(f"source {relative}\n")

    monkeypatch.setattr(campaign, "ROOT", tmp_path / "repository")
    monkeypatch.setattr(campaign, "PARAMETER_RELEASE", parameter_release)
    monkeypatch.setattr(campaign, "AMBERHOME", amberhome)
    monkeypatch.setattr(campaign, "SANDER", sander)
    monkeypatch.setattr(campaign, "SK_PARAMETER_DIR", installed_parameters)
    monkeypatch.setattr(campaign, "EXPECTED_SANDER_SHA256", _sha(sander))
    return base, base / "umbrella_sampling/pilot"


def test_preparation_copies_and_freezes_all_computational_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base, output = _source_tree(tmp_path, monkeypatch)
    prepared = campaign.prepare_pilot(base, output)
    payload = campaign.verify_prepared(output)
    assert prepared.is_file()
    assert len(payload["windows"]) == 18
    assert len(payload["frozen_files_sha256"]) > 60
    assert len(payload["runtime"]["parameter_sha256"]) == 16
    assert payload["native_calculations_launched"] is False
    assert payload["charts"] == ["activation=(q_acyl,q_act)", "delivery=(q_acyl,q_del)"]
    assert len(list(output.glob("*/01_heating_1ps/run.mdin"))) == 18
    assert len(list(output.glob("*/02_equilibration_1ps/run.mdin"))) == 18
    assert len(list(output.glob("*/03_production_2ps/run.mdin"))) == 18


def test_runtime_guard_rejects_changed_stage_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base, output = _source_tree(tmp_path, monkeypatch)
    campaign.prepare_pilot(base, output)
    stage_input = output / "activation_01/01_heating_1ps/run.mdin"
    stage_input.write_text(stage_input.read_text() + "! changed\n")
    with pytest.raises(ValueError, match="Frozen pilot file changed"):
        campaign.verify_prepared(output)


def test_preparation_rejects_existing_output_without_overwrite(tmp_path: Path) -> None:
    base = tmp_path / "relay"
    output = base / "umbrella_sampling/pilot"
    output.mkdir(parents=True)
    with pytest.raises(FileExistsError):
        campaign.prepare_pilot(base, output)
