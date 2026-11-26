"""Focused safeguards for independently forked hydration branches."""

from __future__ import annotations

from pathlib import Path
import json
from types import SimpleNamespace

import pytest

import nyl12_covalent.hydration_replicas as replicas
from nyl12_covalent.hydration_replicas import (
    SPECS,
    _patch_native,
    _patch_postprocess,
    digest,
    verify_hashes,
)


def test_replica_specs_are_disjoint_and_bounded() -> None:
    """Forks use distinct five-block windows and unique independent seeds."""
    assert [(item["checkpoint"], item["first"], item["last"]) for item in SPECS] == [
        (4, 5, 9),
        (7, 8, 12),
    ]
    seeds = [seed for item in SPECS for seed in item["seeds"]]
    assert len(seeds) == len(set(seeds)) == 10
    assert all(len(item["seeds"]) == 5 for item in SPECS)


def test_hash_guard_accepts_identity_and_rejects_change(tmp_path: Path) -> None:
    """Pinned inputs are verified before any worker action."""
    item = tmp_path / "input.dat"
    item.write_text("accepted restart identity\n")
    verify_hashes({str(item): digest(item)})
    item.write_text("changed restart identity\n")
    with pytest.raises(ValueError, match="identity changed"):
        verify_hashes({str(item): "0" * 64})


def test_native_snapshot_is_specialized_without_touching_source(tmp_path: Path) -> None:
    """Copied native code receives only this branch's bounded range and paths."""
    original = tmp_path / "original"
    source = original / "continuation_101ps_source"
    module = source / "nyl12_covalent"
    module.mkdir(parents=True)
    original_extension = (
        'SEGMENT_NAMES = tuple(f"extend_{index:02d}" for index in range(10, 20))\n'
        "SEEDS = (2026112610, 2026112611, 2026112612, 2026112613, 2026112614, "
        "2026112615, 2026112616, 2026112617, 2026112618, 2026112619)\n"
        'CONFIG = "continuation_101ps_config.json"\n'
        'SOURCE = "continuation_101ps_source"\n'
    )
    (module / "deacylation_hydration_extension.py").write_text(original_extension)
    (module / "hydration_followon_gates.py").write_text(
        'CONFIG = "continuation_101ps_config.json"\nSOURCE = "continuation_101ps_source"\n'
    )
    branch = tmp_path / "branch"
    spec = SPECS[0]
    _patch_native(branch, spec, original)
    copied = (
        branch / "native_source/nyl12_covalent/deacylation_hydration_extension.py"
    ).read_text()
    assert "range(5, 10)" in copied
    assert "2026112721" in copied and "2026112725" in copied
    assert "runtime_config.json" in copied and '"native_source"' in copied
    assert (
        module / "deacylation_hydration_extension.py"
    ).read_text() == original_extension
    gate = (
        branch / "native_source/nyl12_covalent/hydration_followon_gates.py"
    ).read_text()
    assert "runtime_config.json" in gate and "native_source" in gate


def test_postprocess_snapshot_caps_global_prefix(tmp_path: Path) -> None:
    """Branch reports scan only the inherited history through its final block."""
    original = tmp_path / "original"
    source = original / "continuation_101ps_postprocess_source"
    (source / "src/nyl12_covalent").mkdir(parents=True)
    adapter = source / "src/nyl12_covalent/hydration_receipt_adapter.py"
    adapter.write_text(
        "for index in range(1, 20):\n    pass\n# extend_01 through extend_19\n"
    )
    (source / "scripts").mkdir()
    script = source / "scripts/refresh_deacylation_hydration.py"
    script.write_text("raise SystemExit(0)\n")
    branch = tmp_path / "branch"
    (branch / "native_source/postprocess").mkdir(parents=True)
    _patch_postprocess(branch, SPECS[0], original)
    text = (
        branch / "postprocess_source/src/nyl12_covalent/hydration_receipt_adapter.py"
    ).read_text()
    assert "range(1, 10)" in text and "extend_01 through extend_09" in text
    assert (
        branch / "native_source/postprocess/hydration_receipt_adapter.py"
    ).read_text() == text
    assert "range(1, 20)" in adapter.read_text()


def test_prepare_branch_builds_only_fork_local_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Preparation links history and stages five inputs without copying MD data."""
    original = tmp_path / "hydration_reconnaissance"
    original.mkdir()
    for name in (
        "config.json",
        "hydration_endpoint_status_20261124.json",
        "extension_20ps_config.json",
        "followon_runtime_config.json",
    ):
        (original / name).write_text("{}\n")
    (original / "inputs").mkdir()
    (original / "inputs/source.mdin").write_text("source\n")
    (original / "inputs/system.prmtop").write_text("topology\n")
    for name in ("heat", "equil"):
        folder = original / name
        folder.mkdir()
        (folder / "restart.ncrst").write_text(name)
    for index in range(1, 5):
        block = original / f"extend_{index:02d}"
        (block / "inputs").mkdir(parents=True)
        for filename in (
            "status.json",
            "process.json",
            "scientific_acceptance.json",
            "restart.ncrst",
            "inputs/restart.ncrst",
            "inputs/source.mdin",
            "inputs/system.prmtop",
            "inputs/segment.mdin",
        ):
            path = block / filename
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"{index}:{filename}")
    fork = original / "extend_04/restart.ncrst"
    template = {
        "pilot_dir": str(original),
        "active_mask": "@1",
        "segment_names": [],
        "seeds": {},
        "source_snapshot_sha256": {},
        "postprocess_source_sha256": {},
        "segment_input_sha256": {},
    }
    (original / "continuation_101ps_template.json").write_text(json.dumps(template))

    def fake_validate(_base: Path, _index: int) -> Path:
        return fork

    def fake_native(branch: Path, _spec: dict, _source: Path) -> dict[str, str]:
        tree = branch / "native_source"
        tree.mkdir()
        marker = tree / "runner.py"
        marker.write_text("frozen runner")
        return {"runner.py": digest(marker)}

    def fake_post(branch: Path, _spec: dict, _source: Path) -> dict[str, str]:
        tree = branch / "postprocess_source/scripts"
        tree.mkdir(parents=True)
        marker = tree / "refresh.py"
        marker.write_text("frozen analysis")
        return {str(marker): digest(marker)}

    monkeypatch.setattr(replicas, "validate_checkpoint", fake_validate)
    monkeypatch.setattr(replicas, "_patch_native", fake_native)
    monkeypatch.setattr(replicas, "_patch_postprocess", fake_post)
    monkeypatch.setattr(
        replicas, "build_input", lambda *_args: "Acyl enzyme hydration only: equil\n"
    )
    result = replicas.prepare_branch(original, SPECS[0], "python")
    branch = Path(result["branch"])
    assert (branch / "extend_01").is_symlink()
    assert (branch / "extend_04").is_symlink()
    assert (branch / "extend_05").is_dir()
    assert not (branch / "extend_05").is_symlink()
    assert (branch / "extend_05/inputs/restart.ncrst").is_symlink()
    assert not (branch / "extend_06/inputs/restart.ncrst").exists()
    assert (branch / "extend_09/status.json").is_file()
    assert (
        json.loads((branch / "runtime_config.json").read_text())["seeds"]["extend_05"]
        == SPECS[0]["seeds"][0]
    )


def test_candidate_check_runs_gate_in_isolated_child(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Candidate screening delegates to the pinned snapshot in a fresh process."""
    observed: dict[str, object] = {}

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        observed["command"] = command
        observed["env"] = kwargs["env"]
        return SimpleNamespace(stdout="clear\n", returncode=0)

    monkeypatch.setattr(replicas.subprocess, "run", fake_run)
    branch, native = tmp_path / "branch", tmp_path / "native"
    assert replicas.candidate_check(branch, native, 4, "python") is False
    assert "hydration_followon_gates" in str(observed["command"])
    assert observed["env"]["PYTHONPATH"] == str(native)
