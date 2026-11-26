"""Tests for bounded preparation of the PA6 101 ps hydration continuation."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np
import pytest

from nyl12_covalent import hydration_continuation_101ps as gates
from scripts.prepare_hydration_continuation_101ps import (
    ROOT,
    SEEDS,
    SEGMENTS,
    patch_native_sources,
    patch_postprocess_sources,
    sha256,
)
from nyl12_covalent.hydration_continuation_sources_101ps import replace_once

BASE = ROOT / "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance"


def test_plan_is_bounded_and_native_float32_101ps_grid_passes() -> None:
    """Use unique explicit seeds and tolerate only native float32 rounding."""
    assert SEGMENTS == tuple(f"extend_{i:02d}" for i in range(10, 20))
    assert len(SEEDS) == 10 and len(set(SEEDS)) == 10
    assert min(SEEDS) > 2026112410
    native = np.asarray(6.0 + np.arange(1, 1901) * 0.05, dtype=np.float32)
    expected = 6.0 + np.arange(1, 1901) * 0.05
    assert all(
        gates.native_time_matches(float(actual), float(target))
        for actual, target in zip(native, expected, strict=True)
    )
    too_far = float(np.nextafter(np.float32(101.0), np.float32(np.inf)))
    assert not gates.native_time_matches(too_far, 101.0)


def test_patch_anchors_fail_closed_on_repeated_text() -> None:
    with pytest.raises(ValueError, match="matched 2 times"):
        replace_once("anchor anchor", "anchor", "changed", "ambiguous")


def test_native_snapshot_selects_each_config_and_exact_extend09_restart() -> None:
    with TemporaryDirectory() as temporary:
        snapshot = Path(temporary) / "continuation_101ps_source"
        hashes = patch_native_sources(BASE / "followon_source_snapshot", snapshot)
        assert hashes
        extension = (
            snapshot / "nyl12_covalent/deacylation_hydration_extension.py"
        ).read_text()
        gates_source = (
            snapshot / "nyl12_covalent/hydration_followon_gates.py"
        ).read_text()
        queue = (snapshot / "nyl12_covalent/hydration_followon_queue.py").read_text()
        assert "range(10, 20)" in extension
        assert "continuation_101ps_config.json" in extension
        assert '"extension_20ps_config.json" if index <= 4' in gates_source
        assert '"followon_runtime_config.json" if index <= 9' in gates_source
        assert '"continuation_101ps_config.json"' in gates_source
        assert 'completion.get("status") != "completed_51ps_review_required"' in queue
        assert 'cfg = read(base / "continuation_101ps_template.json")' in queue
        assert '(base / "extend_09/restart.ncrst").read_bytes()' in queue
        assert 'save(base, "continuation_101ps_config.json", cfg)' in queue
        assert "hold(plan, 9)" in queue
        assert 'enumerate(plan["segment_names"], start=10)' in queue


def test_postprocess_snapshot_keeps_current_vmd_and_fixes_conditional_proton_wording() -> (
    None
):
    with TemporaryDirectory() as temporary:
        snapshot = Path(temporary) / "continuation_101ps_postprocess_source"
        hashes = patch_postprocess_sources(
            BASE / "followon_postprocess_source", snapshot
        )
        assert hashes
        package = snapshot / "src/nyl12_covalent"
        for filename in (
            "hydration_vmd_exporter.py",
            "hydration_vmd_viewer.py",
            "hydration_vmd_export_io.py",
            "vmd_movie_policy.py",
            "relay_progress_vmd.py",
            "vmd_chain_groups.py",
        ):
            assert sha256(package / filename) == sha256(
                ROOT / "src/nyl12_covalent" / filename
            )
        adapter = (package / "hydration_receipt_adapter.py").read_text()
        assert "range(1, 20)" in adapter and "extend_01 through extend_19" in adapter
        analysis = (package / "deacylation_water_analysis.py").read_text()
        report = (package / "deacylation_water_reporting.py").read_text()
        assert "Nalpha is protonated" not in analysis + report
        assert "accepted segment" in analysis + report


def test_patched_queue_prepares_ten_inputs_from_extend09_without_running_md() -> None:
    """Exercise the real cloned prepare_runtime path in a temporary campaign."""
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        base = root / "campaign"
        base.mkdir()
        snapshot = base / "continuation_101ps_source"
        patch_native_sources(BASE / "followon_source_snapshot", snapshot)
        (base / "inputs").mkdir()
        (base / "inputs/source.mdin").write_bytes(
            (BASE / "inputs/source.mdin").read_bytes()
        )
        (base / "inputs/system.prmtop").write_bytes(
            (BASE / "inputs/system.prmtop").read_bytes()
        )
        predecessor = base / "extend_09"
        predecessor.mkdir()
        (predecessor / "restart.ncrst").write_bytes(
            (BASE / "extend_09/restart.ncrst").read_bytes()
        )
        config = json.loads((BASE / "followon_runtime_config.json").read_text())
        (base / "continuation_101ps_template.json").write_text(json.dumps(config))
        plan = {
            "base": str(base),
            "segment_names": list(SEGMENTS),
            "seeds": dict(zip(SEGMENTS, SEEDS, strict=True)),
            "postprocess_argv": [sys.executable],
            "postprocess_environment": {"OMP_NUM_THREADS": "1"},
            "recovery_source_sha256": {"source": "sha256"},
        }
        plan_path = base / "test_plan.json"
        plan_path.write_text(json.dumps(plan))
        code = (
            "import json, sys; "
            f"sys.path.insert(0, {str(snapshot)!r}); "
            "from nyl12_covalent.hydration_followon_queue import prepare_runtime; "
            f"prepare_runtime(json.load(open({str(plan_path)!r})))"
        )
        environment = dict(os.environ)
        environment["PYTHONPATH"] = os.pathsep.join(
            [str(snapshot), str(ROOT / "src"), environment.get("PYTHONPATH", "")]
        )
        subprocess.run([sys.executable, "-c", code], check=True, env=environment)
        runtime = json.loads((base / "continuation_101ps_config.json").read_text())
        assert runtime["endpoint_restart_path"].endswith("extend_09/restart.ncrst")
        assert runtime["seeds"] == plan["seeds"]
        for name in SEGMENTS:
            assert (
                json.loads((base / name / "status.json").read_text())["status"]
                == "prepared"
            )
            assert (base / name / "inputs/segment.mdin").is_file()
            assert (base / name / "inputs/source.mdin").is_file()
            assert (base / name / "inputs/system.prmtop").is_file()
            assert not (base / name / "run.mdout").exists()
            assert not (base / name / "command.json").exists()
        assert (base / "extend_10/inputs/restart.ncrst").read_bytes() == (
            BASE / "extend_09/restart.ncrst"
        ).read_bytes()
        assert not (base / "extend_11/inputs/restart.ncrst").exists()


def test_attempted_continuation_must_be_a_contiguous_prefix(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    with TemporaryDirectory() as temporary:
        base = Path(temporary)
        monkeypatch.setattr(
            gates, "validate_accepted_segment", lambda _base, _index: None
        )
        for index, status in ((10, "prepared"), (12, "accepted")):
            folder = base / f"extend_{index:02d}"
            folder.mkdir()
            (folder / "status.json").write_text(json.dumps({"status": status}))
        with pytest.raises(ValueError, match="not contiguous"):
            gates.validate_prefix(base, 10, 19)


def test_candidate_screen_stops_and_requires_fresh_manifest() -> None:
    with TemporaryDirectory() as temporary:
        base = Path(temporary)
        analysis = base / "analysis"
        analysis.mkdir()
        marker = base / "accepted_source.json"
        marker.write_text("stable")
        (analysis / "source_manifest.json").write_text(
            json.dumps({"sha256": {str(marker): sha256(marker)}})
        )
        stage_counts = {f"extend_{i:02d}": 100 for i in range(1, 10)}
        (analysis / "report.json").write_text(
            json.dumps(
                {
                    "frame_count_by_stage": stage_counts,
                    "criteria": {
                        "near_attack_O_C_A_max": 3.5,
                        "attack_angle_deg": [95, 125],
                    },
                    "water_count": 1,
                }
            )
        )
        (analysis / "per_water_summary.json").write_text(
            json.dumps(
                {
                    "waters": [
                        {"near_attack_frame_count": 1},
                    ]
                }
            )
        )
        assert gates.candidate_found(base, 9) is True
        marker.write_text("changed")
        with pytest.raises(ValueError, match="Analysis source changed"):
            gates.candidate_found(base, 9)


def _write_accepted_history(base: Path, through: int) -> None:
    """Build a tiny internally consistent native restart chain for gate tests."""
    input_root = base / "equil"
    input_root.mkdir(parents=True, exist_ok=True)
    (input_root / "restart.ncrst").write_bytes(b"accepted-equil-restart")
    for index in range(1, through + 1):
        name = f"extend_{index:02d}"
        block = base / name
        inputs = block / "inputs"
        inputs.mkdir(parents=True)
        previous = base / ("equil" if index == 1 else f"extend_{index - 1:02d}")
        start_restart = previous / "restart.ncrst"
        (inputs / "restart.ncrst").write_bytes(start_restart.read_bytes())
        (inputs / "segment.mdin").write_text(f"native-input-{index}")
        (block / "restart.ncrst").write_bytes(f"accepted-restart-{index}".encode())
        (block / "trajectory.nc").write_bytes(b"synthetic-nonempty-trajectory")
        start_hash = sha256(start_restart)
        (block / "status.json").write_text(
            json.dumps(
                {
                    "status": "accepted",
                    "exit_code": 0,
                    "start_restart_sha256": start_hash,
                    "restart_sha256": sha256(block / "restart.ncrst"),
                }
            )
        )
        (block / "process.json").write_text(
            json.dumps(
                {
                    "status": "process_completed",
                    "exit_code": 0,
                    "start_restart_sha256": start_hash,
                }
            )
        )
        (block / "scientific_acceptance.json").write_text(
            json.dumps(
                {
                    "numeric": {
                        "completed_steps": 10000,
                        "electronic_history": {"verified_scc_history": True},
                        "final_temperature_K": 300.0,
                    },
                    "boundary": {
                        "box_unchanged": True,
                        "frozen_max_displacement_A": 0.0,
                    },
                    "C_OG_A": 1.33,
                    "scc_verified": True,
                    "acyl_retained": True,
                    "finite_temperature_passed": True,
                    "native_process_passed": True,
                }
            )
        )
        config_name = (
            "extension_20ps_config.json"
            if index <= 4
            else "followon_runtime_config.json"
            if index <= 9
            else "continuation_101ps_config.json"
        )
        config_path = base / config_name
        config = json.loads(config_path.read_text()) if config_path.exists() else {}
        input_hashes = {path.name: sha256(path) for path in inputs.iterdir()}
        config.setdefault("segment_input_sha256", {})[name] = input_hashes
        config_path.write_text(json.dumps(config))


def test_native_acceptance_validates_full_historical_and_current_restart_chain() -> (
    None
):
    with TemporaryDirectory() as temporary:
        base = Path(temporary)
        _write_accepted_history(base, 10)
        assert gates.validate_prefix(base, 10, 10) == 1


def test_native_acceptance_rejects_chemistry_and_input_hash_failures() -> None:
    with TemporaryDirectory() as temporary:
        base = Path(temporary)
        _write_accepted_history(base, 10)
        acceptance_path = base / "extend_10/scientific_acceptance.json"
        acceptance = json.loads(acceptance_path.read_text())
        acceptance["acyl_retained"] = False
        acceptance_path.write_text(json.dumps(acceptance))
        with pytest.raises(ValueError, match="Native acceptance or restart lineage"):
            gates.validate_accepted_segment(base, 10)

    with TemporaryDirectory() as temporary:
        base = Path(temporary)
        _write_accepted_history(base, 10)
        input_path = base / "extend_10/inputs/segment.mdin"
        input_path.write_text("changed-input")
        with pytest.raises(ValueError, match="Native acceptance or restart lineage"):
            gates.validate_accepted_segment(base, 10)


def test_prefix_returns_accepted_attempt_count_and_stops_at_nonaccepted_block() -> None:
    with TemporaryDirectory() as temporary:
        base = Path(temporary)
        _write_accepted_history(base, 9)
        assert gates.validate_prefix(base, 10, 12) == 0

    with TemporaryDirectory() as temporary:
        base = Path(temporary)
        _write_accepted_history(base, 10)
        (base / "extend_11").mkdir()
        (base / "extend_11/status.json").write_text(json.dumps({"status": "failed"}))
        assert gates.validate_prefix(base, 10, 12) == 1

    with TemporaryDirectory() as temporary:
        base = Path(temporary)
        _write_accepted_history(base, 12)
        (base / "extend_11/status.json").write_text(json.dumps({"status": "failed"}))
        with pytest.raises(ValueError, match="not contiguous"):
            gates.validate_prefix(base, 10, 12)


def test_candidate_gate_rejects_incomplete_stages_changed_criteria_and_water_count() -> (
    None
):
    with TemporaryDirectory() as temporary:
        base = Path(temporary)
        analysis = base / "analysis"
        analysis.mkdir()
        marker = base / "source.json"
        marker.write_text("fixed")
        (analysis / "source_manifest.json").write_text(
            json.dumps({"sha256": {str(marker): sha256(marker)}})
        )
        report_path = analysis / "report.json"
        summary_path = analysis / "per_water_summary.json"
        report = {
            "frame_count_by_stage": {f"extend_{i:02d}": 100 for i in range(1, 10)},
            "criteria": {"near_attack_O_C_A_max": 3.5, "attack_angle_deg": [95, 125]},
            "water_count": 1,
        }
        summary_path.write_text(
            json.dumps({"waters": [{"near_attack_frame_count": 0}]})
        )
        report_path.write_text(json.dumps(report))
        report["frame_count_by_stage"].pop("extend_09")
        report_path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="does not cover"):
            gates.candidate_found(base, 9)
        report["frame_count_by_stage"]["extend_09"] = 100
        report["criteria"]["attack_angle_deg"] = [90, 130]
        report_path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="screen changed"):
            gates.candidate_found(base, 9)
        report["criteria"]["attack_angle_deg"] = [95, 125]
        report["water_count"] = 2
        report_path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="incomplete"):
            gates.candidate_found(base, 9)
