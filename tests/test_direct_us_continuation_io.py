"""Checks immutable continuation plan and accepted-stage validation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from nyl12_covalent.direct_us_continuation_io import (
    _trajectory_digest,
    requested_block_count,
    validate_accepted_stage,
    validate_continuation_plan,
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _plan_fixture(tmp_path: Path, *, materialize_through: int = 7) -> tuple[Path, Path]:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    windows = []
    stages = []
    for index in range(1, 64):
        label = f"window_{index:03d}"
        window_stages = []
        for block in range(7, 13):
            stage_label = f"03_production_block_{block:02d}_5ps"
            stage = campaign / label / stage_label
            window_stages.append(
                {"label": stage_label, "kind": "production", "block": block}
            )
            row = {
                "window": label,
                "block": block,
                "label": stage_label,
                "stage": str(stage.resolve()),
                "mdin_sha256": "missing",
                "restraint_sha256": "missing",
            }
            if block <= materialize_through:
                stage.mkdir(parents=True, exist_ok=True)
                (stage / "run.mdin").write_text(f"mdin {label} {block}")
                (stage / "umbrella.RST").write_text(f"restraint {label}")
                row["mdin_sha256"] = _sha(stage / "run.mdin")
                row["restraint_sha256"] = _sha(stage / "umbrella.RST")
            stages.append(row)
        windows.append({"label": label, "stages": window_stages})
    prepared = campaign / "prepared.json"
    prepared.write_text(
        json.dumps({"windows": [{"label": w["label"]} for w in windows]})
    )
    plan = {
        "schema_version": 1,
        "campaign_root": str(campaign.resolve()),
        "prepared": str(prepared.resolve()),
        "prepared_sha256": _sha(prepared),
        "source_sha256": {},
        "windows": windows,
        "stages": stages,
        "augmentation_root": str((campaign / "augmented").resolve()),
    }
    plan_dir = campaign / "continuation_60ps_20261126"
    plan_dir.mkdir()
    plan_path = plan_dir / "plan.json"
    plan_path.write_text(json.dumps(plan))
    (plan_dir / "plan.sha256").write_text(_sha(plan_path) + "\n")
    return campaign, plan_path


@pytest.mark.parametrize(
    "horizon,blocks", [(35, 7), (40, 8), (45, 9), (50, 10), (55, 11), (60, 12)]
)
def test_requested_block_count(horizon: float, blocks: int) -> None:
    assert requested_block_count(horizon) == blocks


@pytest.mark.parametrize("horizon", [30, 36, 65, float("nan")])
def test_requested_block_count_rejects_unsupported_horizons(horizon: float) -> None:
    with pytest.raises(ValueError, match="35, 40, 45"):
        requested_block_count(horizon)


def test_trajectory_hash_cache_reuses_only_unchanged_file_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import nyl12_covalent.direct_us_continuation_io as continuation_io

    trajectory = tmp_path / "trajectory.nc"
    trajectory.write_bytes(b"first coordinates")
    real_digest = continuation_io.digest
    calls = 0

    def counted_digest(path: Path) -> str:
        nonlocal calls
        calls += 1
        return real_digest(path)

    monkeypatch.setattr(continuation_io, "digest", counted_digest)
    first = _trajectory_digest(trajectory)
    assert _trajectory_digest(trajectory) == first
    assert calls == 1
    trajectory.write_bytes(b"changed coordinates")
    assert _trajectory_digest(trajectory) != first
    assert calls == 2


def test_plan_validation_checks_only_requested_native_inputs(tmp_path: Path) -> None:
    campaign, plan = _plan_fixture(tmp_path)

    by_window, sources = validate_continuation_plan(campaign, plan, 35)

    assert len(by_window) == 63
    assert all(
        len(records) == 1 and records[0]["block"] == 7 for records in by_window.values()
    )
    assert str(plan.resolve()) in sources
    assert not (campaign / "window_001/03_production_block_08_5ps/run.mdin").exists()


def test_plan_validation_rejects_missing_requested_block(tmp_path: Path) -> None:
    campaign, plan = _plan_fixture(tmp_path)
    plan_data = json.loads(plan.read_text())
    plan_data["stages"] = [row for row in plan_data["stages"] if row["block"] != 7]
    plan.write_text(json.dumps(plan_data))
    (plan.parent / "plan.sha256").write_text(_sha(plan) + "\n")

    with pytest.raises(ValueError, match="skips required production blocks"):
        validate_continuation_plan(campaign, plan, 35)


def test_plan_validation_rejects_requested_input_hash_change(tmp_path: Path) -> None:
    campaign, plan = _plan_fixture(tmp_path)
    (campaign / "window_001/03_production_block_07_5ps/run.mdin").write_text("changed")

    with pytest.raises(ValueError, match="input hash mismatch"):
        validate_continuation_plan(campaign, plan, 35)


def test_accepted_stage_requires_native_and_augmentation_receipts(
    tmp_path: Path,
) -> None:
    campaign = tmp_path / "campaign"
    stage = campaign / "window_001/03_production_block_07_5ps"
    native_csv = stage / "observables/observables.csv"
    native_csv.parent.mkdir(parents=True)
    native_csv.write_text(
        "frame,time_ps\n"
        + "".join(f"{frame},{frame / 20:.2f}\n" for frame in range(100))
    )
    trajectory = stage / "trajectory.nc"
    trajectory.write_bytes(b"synthetic native coordinates")
    paths = {
        "process.json": stage / "process.json",
        "acceptance.json": stage / "acceptance.json",
        "observables/observables.json": stage / "observables/observables.json",
    }
    paths["process.json"].write_text('{"exit_status":0,"ended_utc":"now"}')
    paths["acceptance.json"].write_text(
        '{"completed_steps":10000,"finite_energy_temperature":true,"SCC_error_detected":false}'
    )
    native_receipt = {
        "frames": 100,
        "frames_censored": 0,
        "csv_sha256": _sha(native_csv),
        "source_sha256": {str(trajectory): _sha(trajectory)},
    }
    paths["observables/observables.json"].write_text(json.dumps(native_receipt))
    augmented_root = campaign / "augmented"
    derived = augmented_root / "window_001/03_production_block_07_5ps"
    derived.mkdir(parents=True)
    derived_csv = derived / "observables.csv"
    derived_csv.write_text(
        "frame,q_acyl\n" + "".join(f"{frame},0.0\n" for frame in range(100))
    )
    source_hashes = {
        str(path): _sha(path) for path in (*paths.values(), native_csv, trajectory)
    }
    receipt = {
        "rows_count": 100,
        "frames_censored": 0,
        "csv_sha256": _sha(derived_csv),
        "source_sha256": source_hashes,
        "parity": {
            "all_four_distances_passed": True,
            "q_h_partition_passed": True,
        },
    }
    (derived / "augmentation.json").write_text(json.dumps(receipt))
    record = {
        "window": "window_001",
        "label": "03_production_block_07_5ps",
        "stage": str(stage),
    }

    hashes = validate_accepted_stage(record, augmented_root)

    assert len(hashes) >= 6
    assert str(derived_csv) in hashes

    paths["acceptance.json"].write_text('{"completed_steps":9999}')
    with pytest.raises(ValueError, match="acceptance failed"):
        validate_accepted_stage(record, augmented_root)
    paths["acceptance.json"].write_text(
        '{"completed_steps":10000,"finite_energy_temperature":true,"SCC_error_detected":false}'
    )
    trajectory.write_bytes(b"changed native coordinates")
    with pytest.raises(ValueError, match="raw trajectory hash failed"):
        validate_accepted_stage(record, augmented_root)
