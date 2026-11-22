"""Check the hashed, persistent single-GPU relaxation queue contract."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from scripts import watch_design_relaxation_queue as worker


def setup_queue(tmp_path: Path, ready: bool = True) -> tuple[Path, Path]:
    """Create a minimal one-case queue with a canonical, hashed model request."""
    (tmp_path / "prepared").mkdir()
    (tmp_path / "configuration").mkdir()
    (tmp_path / "prepared/candidate_manifest.json").write_text(
        json.dumps({"jobs": [{"variant": "WT", "model_index": 0, "seed": 1}]})
    )
    (tmp_path / "configuration/relaxation.json").write_text("{}")
    case = tmp_path / "models/WT/seed_01"
    case.mkdir(parents=True)
    source = case / "normalized.pdb"
    source.write_text("UNCHANGED_COORDINATES")
    request = case / "request.json"
    if ready:
        request.write_text(
            json.dumps(
                {
                    "source": str(source),
                    "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                }
            )
        )
    return case, request


def invoke(tmp_path: Path) -> None:
    """Run the queue CLI with one explicitly selected GPU."""
    with patch("sys.argv", ["queue", "--campaign", str(tmp_path), "--gpu", "0"]):
        worker.main()


def test_queue_runs_once_and_resumes_only_completed_state(tmp_path: Path) -> None:
    """A valid source is minimized once; an explicit completed receipt prevents rerun."""
    case, _ = setup_queue(tmp_path)
    with patch.object(
        worker, "relax_model", return_value={"gradient_converged": True}
    ) as relax:
        invoke(tmp_path)
        invoke(tmp_path)
    assert relax.call_count == 1
    assert relax.call_args.args[2]["gpu"] == "0"
    assert (
        json.loads((case / "relaxation_worker_status.json").read_text())["status"]
        == "completed"
    )


def test_queue_rejects_modified_source_before_engine(tmp_path: Path) -> None:
    """Hash mismatch prevents minimization instead of silently accepting changed inputs."""
    case, _ = setup_queue(tmp_path)
    (case / "normalized.pdb").write_text("ALTERED")
    with (
        patch.object(worker, "relax_model") as relax,
        pytest.raises(ValueError, match="changed"),
    ):
        invoke(tmp_path)
    relax.assert_not_called()


def test_queue_records_failure_and_requires_review(tmp_path: Path) -> None:
    """A backend failure is recorded and raised; a later run does not retry it silently."""
    case, _ = setup_queue(tmp_path)
    with (
        patch.object(worker, "relax_model", side_effect=RuntimeError("LEaP failed")),
        pytest.raises(RuntimeError, match="LEaP failed"),
    ):
        invoke(tmp_path)
    assert (
        json.loads((case / "relaxation_worker_status.json").read_text())["status"]
        == "failed"
    )
    with (
        patch.object(worker, "relax_model") as relax,
        pytest.raises(RuntimeError, match="restart review"),
    ):
        invoke(tmp_path)
    relax.assert_not_called()


def test_queue_waits_for_request_without_reserving_gpu_work(tmp_path: Path) -> None:
    """A not-yet-produced model is polled without invoking the scientific backend."""
    case, request = setup_queue(tmp_path, False)

    def publish(_: float) -> None:
        source = case / "normalized.pdb"
        request.write_text(
            json.dumps(
                {
                    "source": str(source),
                    "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                }
            )
        )

    with (
        patch.object(worker.time, "sleep", side_effect=publish) as sleep,
        patch.object(worker, "relax_model", return_value={}),
    ):
        invoke(tmp_path)
    sleep.assert_called_once_with(30)
