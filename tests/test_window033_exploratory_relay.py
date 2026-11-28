"""Tests for explicit authorization and force-gated exploratory dispatch."""

from __future__ import annotations

import hashlib
import json
from types import ModuleType
from pathlib import Path
import sys

import pytest

import scripts.run_window033_exploratory_relay as dispatcher
from scripts.run_window033_exploratory_relay import (
    AUTHORIZATION_FILENAME,
    BRANCHES,
    DISPATCH_FILENAME,
    PLAN_FILENAME,
    run_exploratory_branch,
    validate_dispatch,
)


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, sort_keys=True) + "\n")


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _campaign(tmp_path: Path) -> Path:
    campaign = tmp_path / "campaign"
    model = campaign / "model"
    (model / "inputs").mkdir(parents=True)
    (model / "inputs/optimized.ncrst").write_bytes(b"restart")
    (model / "inputs/system.prmtop").write_bytes(b"topology")
    identity = {
        "seed": "ti1_window033",
        "full_atom_count": 104922,
        "active_atom_count": 4165,
        "qm_atom_count": 98,
        "qm_charge": -2,
        "coordinate_sha256": hashlib.sha256(b"unoptimized trajectory seed").hexdigest(),
        "topology_sha256": _digest(model / "inputs/system.prmtop"),
    }
    acceptance = {
        "accepted": True,
        "force_converged": True,
        "native_converged": True,
        "max_gradient_hartree_per_bohr": 8e-6,
    }
    _write_json(model / "source_identity.json", identity)
    _write_json(model / "release_acceptance.json", acceptance)
    model_manifest = {
        "inputs/optimized.ncrst": _digest(model / "inputs/optimized.ncrst"),
        "inputs/system.prmtop": _digest(model / "inputs/system.prmtop"),
        "source_identity.json": _digest(model / "source_identity.json"),
        "release_acceptance.json": _digest(model / "release_acceptance.json"),
    }
    _write_json(model / "manifest.json", model_manifest)
    frozen_source = campaign / "frozen_source.py"
    frozen_source.write_text("approved source\n")
    config = {
        "model": "model",
        "branches": BRANCHES,
        "source_sha256": {str(frozen_source): _digest(frozen_source)},
        "source_topology_sha256": _digest(model / "inputs/system.prmtop"),
    }
    _write_json(campaign / "config.json", config)
    _write_json(
        campaign / "preparation.json",
        {
            "config_sha256": _digest(campaign / "config.json"),
            "source_frozen_manifest_sha256": _digest(model / "manifest.json"),
        },
    )
    dispatcher = campaign / "dispatcher.py"
    dispatcher.write_text("authorized fixture dispatcher\n")
    branch_root = campaign / BRANCHES["forward_acylation"]
    authorization = {
        "status": "approved_exploratory_only",
        "authorized_branches": ["forward_acylation"],
    }
    _write_json(branch_root / AUTHORIZATION_FILENAME, authorization)
    _write_json(
        branch_root / PLAN_FILENAME,
        {
            "branch": "forward_acylation",
            "dispatcher": str(dispatcher),
            "dispatcher_sha256": _digest(dispatcher),
            "config_sha256": _digest(campaign / "config.json"),
            "authorization_sha256": _digest(
                branch_root / AUTHORIZATION_FILENAME
            ),
        },
    )
    return campaign


def test_force_accepted_seed_passes_without_full_minimum_claim(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)

    report = validate_dispatch(campaign, "forward_acylation")

    assert report["status"] == "preflight_passed_exploratory_only"
    assert report["window033_force_accepted"] is True
    assert report["full_minimum_verified"] is False
    assert report["ts_claim_authorized"] is False


def test_source_hash_mismatch_blocks_dispatch(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    (campaign / "model/inputs/optimized.ncrst").write_bytes(b"changed")

    with pytest.raises(ValueError, match="hash mismatch"):
        validate_dispatch(campaign, "forward_acylation")


def test_named_frozen_source_hash_mismatch_blocks_dispatch(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    (campaign / "frozen_source.py").write_text("changed source\n")

    with pytest.raises(ValueError, match="source hash mismatch"):
        validate_dispatch(campaign, "forward_acylation")


def test_nonaccepted_force_state_blocks_dispatch(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    path = campaign / "model/release_acceptance.json"
    acceptance = json.loads(path.read_text())
    acceptance["force_converged"] = False
    _write_json(path, acceptance)
    model = campaign / "model"
    manifest = json.loads((model / "manifest.json").read_text())
    manifest["release_acceptance.json"] = _digest(path)
    _write_json(model / "manifest.json", manifest)
    preparation = json.loads((campaign / "preparation.json").read_text())
    preparation["source_frozen_manifest_sha256"] = _digest(model / "manifest.json")
    _write_json(campaign / "preparation.json", preparation)

    with pytest.raises(ValueError, match="force-acceptance"):
        validate_dispatch(campaign, "forward_acylation")


def test_unapproved_branch_is_rejected(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)

    with pytest.raises(ValueError, match="authorization"):
        validate_dispatch(campaign, "reverse_activation")


def test_unknown_branch_is_rejected(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)

    with pytest.raises(ValueError, match="Unsupported relay branch"):
        validate_dispatch(campaign, "direct_proton_transfer")


def test_config_hash_change_blocks_dispatch(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    config_path = campaign / "config.json"
    config = json.loads(config_path.read_text())
    config["unexpected"] = True
    _write_json(config_path, config)

    with pytest.raises(ValueError, match="config differs"):
        validate_dispatch(campaign, "forward_acylation")


def test_non_object_preparation_record_is_rejected(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    (campaign / "preparation.json").write_text("[]\n")

    with pytest.raises(ValueError, match="Expected JSON object"):
        validate_dispatch(campaign, "forward_acylation")


def test_frozen_manifest_change_blocks_dispatch(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    model_manifest = campaign / "model/manifest.json"
    manifest = json.loads(model_manifest.read_text())
    manifest["unexpected"] = "0" * 64
    _write_json(model_manifest, manifest)

    with pytest.raises(ValueError, match="manifest differs"):
        validate_dispatch(campaign, "forward_acylation")


def test_branch_map_change_blocks_dispatch(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    config_path = campaign / "config.json"
    config = json.loads(config_path.read_text())
    config["branches"] = {**BRANCHES, "unexpected": "extra"}
    _write_json(config_path, config)
    preparation_path = campaign / "preparation.json"
    preparation = json.loads(preparation_path.read_text())
    preparation["config_sha256"] = _digest(config_path)
    _write_json(preparation_path, preparation)

    with pytest.raises(ValueError, match="branch map"):
        validate_dispatch(campaign, "forward_acylation")


def test_plan_hash_mismatch_blocks_dispatch(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    plan_path = campaign / BRANCHES["forward_acylation"] / PLAN_FILENAME
    plan = json.loads(plan_path.read_text())
    plan["branch"] = "reverse_activation"
    _write_json(plan_path, plan)

    with pytest.raises(ValueError, match="authorization or branch plan"):
        validate_dispatch(campaign, "forward_acylation")


def test_identity_mismatch_blocks_dispatch(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    model = campaign / "model"
    identity_path = model / "source_identity.json"
    identity = json.loads(identity_path.read_text())
    identity["qm_charge"] = 0
    _write_json(identity_path, identity)
    manifest_path = model / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["source_identity.json"] = _digest(identity_path)
    _write_json(manifest_path, manifest)
    preparation_path = campaign / "preparation.json"
    preparation = json.loads(preparation_path.read_text())
    preparation["source_frozen_manifest_sha256"] = _digest(manifest_path)
    _write_json(preparation_path, preparation)

    with pytest.raises(ValueError, match="source identity"):
        validate_dispatch(campaign, "forward_acylation")


def test_launch_writes_exploratory_receipt_without_minimum_claim(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = _campaign(tmp_path)
    parent = ModuleType("nyl12_covalent")
    parent.__path__ = []
    implementation = ModuleType("nyl12_covalent.nterminal_relay_run")
    implementation.run_branch = lambda campaign, branch: {"status": "completed"}
    monkeypatch.setitem(sys.modules, "nyl12_covalent", parent)
    monkeypatch.setitem(sys.modules, "nyl12_covalent.nterminal_relay_run", implementation)

    result = run_exploratory_branch(campaign, "forward_acylation")
    receipt = json.loads(
        (campaign / BRANCHES["forward_acylation"] / DISPATCH_FILENAME).read_text()
    )

    assert result["status"] == "completed"
    assert receipt["status"] == "path_implementation_returned"
    assert receipt["full_minimum_verified"] is False
    assert receipt["formal_ts_claim"] is False
    with pytest.raises(FileExistsError, match="Preserve existing"):
        run_exploratory_branch(campaign, "forward_acylation")


def test_failed_native_dispatch_preserves_failure_receipt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = _campaign(tmp_path)
    parent = ModuleType("nyl12_covalent")
    parent.__path__ = []
    implementation = ModuleType("nyl12_covalent.nterminal_relay_run")

    def fail_branch(campaign: Path, branch: str) -> dict:
        raise RuntimeError("native branch failed")

    implementation.run_branch = fail_branch
    monkeypatch.setitem(sys.modules, "nyl12_covalent", parent)
    monkeypatch.setitem(sys.modules, "nyl12_covalent.nterminal_relay_run", implementation)

    with pytest.raises(RuntimeError, match="native branch failed"):
        run_exploratory_branch(campaign, "forward_acylation")

    receipt = json.loads(
        (campaign / BRANCHES["forward_acylation"] / DISPATCH_FILENAME).read_text()
    )
    assert receipt["status"] == "dispatch_failed"
    assert receipt["error"] == "RuntimeError: native branch failed"
    assert receipt["full_minimum_verified"] is False
    assert receipt["formal_ts_claim"] is False


def test_cli_preflight_returns_success_without_launch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = _campaign(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_window033_exploratory_relay.py",
            "--campaign",
            str(campaign),
            "--branch",
            "forward_acylation",
        ],
    )

    assert dispatcher.main() == 0


@pytest.mark.parametrize(("branch_status", "expected_exit"), [("completed", 0), ("failed", 1)])
def test_cli_launch_exit_status(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    branch_status: str,
    expected_exit: int,
) -> None:
    campaign = _campaign(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_window033_exploratory_relay.py",
            "--campaign",
            str(campaign),
            "--branch",
            "forward_acylation",
            "--launch",
        ],
    )
    monkeypatch.setattr(
        dispatcher,
        "run_exploratory_branch",
        lambda campaign, branch: {"status": branch_status},
    )

    assert dispatcher.main() == expected_exit
