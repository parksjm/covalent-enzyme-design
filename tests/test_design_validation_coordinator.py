"""Check distributed normalization, queue-readiness and model-level analysis glue."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from scripts import coordinate_design_validation as control


def job(index: int = 0) -> dict:
    """Return one explicitly seeded planned case."""
    return {
        "variant": "WT",
        "model_index": index,
        "seed": index + 1,
        "expected_alpha_sequence": "A" * 232,
        "expected_beta_sequence": "T" * 96,
    }


def publishable(base: Path, item: dict) -> dict:
    """Create native-complete artifacts for one normalization request."""
    path = base / "predictions/canonical/WT" / f"WT_model_{item['model_index']}.pdb"
    path.parent.mkdir(parents=True)
    path.write_text("RAW_COORDINATES")
    path.with_name(f"confidence_{path.stem}.json").write_text("{}")
    component = control.component_path(base, item)
    component.parent.mkdir(parents=True)
    component.write_bytes(b"TRUSTED_NATIVE_COMPONENT")
    return {"mol2": str(base / "reference.mol2")}


def normalize(source: Path, component: Path, reference: Path, out: Path) -> dict:
    """Mock the separately tested graph normalizer without changing coordinates."""
    out.mkdir()
    (out / "normalized.pdb").write_bytes(source.read_bytes())
    return {"role_mapping": {"C54": "C54", "O49": "O49"}}


def sequences() -> dict:
    """Provide exact expected eight-chain sequences for the glue tests."""
    return {**dict.fromkeys("ABCD", "A" * 232), **dict.fromkeys("EFGH", "T" * 96)}


def test_publish_requires_complete_native_outputs_and_preserves_hash(
    tmp_path: Path,
) -> None:
    """Publish once only after normalization and identity checks; verify later reads."""
    item = job()
    manifest = {"jobs": [item]}
    assert control.publish_requests(tmp_path, manifest, {"mol2": "x"}) == 0
    config = publishable(tmp_path, item)
    with (
        patch.object(control, "normalize_model", side_effect=normalize) as normal,
        patch.object(control, "_pdb_sequences", return_value=sequences()),
    ):
        assert control.publish_requests(tmp_path, manifest, config) == 1
        assert control.publish_requests(tmp_path, manifest, config) == 1
    assert normal.call_count == 1
    request = json.loads((tmp_path / "models/WT/seed_01/request.json").read_text())
    assert request["ligand_count"] == 4
    assert request["sha256"] == hashlib.sha256(b"RAW_COORDINATES").hexdigest()
    Path(request["source"]).write_text("MODIFIED")
    with pytest.raises(ValueError, match="changed"):
        control.publish_requests(tmp_path, manifest, config)


def test_identity_change_blocks_readiness(tmp_path: Path) -> None:
    """A changed protein sequence must never be submitted as the named candidate."""
    item = job()
    config = publishable(tmp_path, item)
    observed = sequences()
    observed["A"] = "C" * 232
    with (
        patch.object(control, "normalize_model", side_effect=normalize),
        patch.object(control, "_pdb_sequences", return_value=observed),
        pytest.raises(ValueError, match="alpha sequence"),
    ):
        control.publish_requests(tmp_path, {"jobs": [item]}, config)
    assert not (tmp_path / "models/WT/seed_01/request.json").exists()


def test_component_paths_distinguish_preserved_sentinel(tmp_path: Path) -> None:
    """The separately submitted first WT case must reuse its actual component."""
    assert "/WT/seed_01/" in str(control.component_path(tmp_path, job()))
    assert "/runs/WT/model_1/output/" in str(control.component_path(tmp_path, job(1)))


def configure_main(base: Path) -> list[dict]:
    """Create45 explicit cases and a minimal controller configuration."""
    (base / "prepared").mkdir()
    (base / "configuration").mkdir()
    jobs = [job(i) for i in range(45)]
    (base / "prepared/candidate_manifest.json").write_text(json.dumps({"jobs": jobs}))
    (base / "configuration/relaxation.json").write_text(
        json.dumps(
            {"ssh_identity": "/authorized/key", "reference_cif": "reference.cif"}
        )
    )
    return jobs


def test_requests_are_transferred_after_referenced_models(tmp_path: Path) -> None:
    """Network ordering must keep remote workers from observing incomplete requests."""
    jobs = configure_main(tmp_path)
    for j in jobs:
        case = tmp_path / "models/WT" / f"seed_{j['model_index'] + 1:02d}"
        case.mkdir(parents=True)
        (case / "relaxation_worker_status.json").write_text('{"status":"completed"}')
    with (
        patch("sys.argv", ["control", "--campaign", str(tmp_path)]),
        patch.object(control, "sync_tree") as sync,
        patch.object(control, "publish_requests", return_value=45),
        patch.object(control, "final_comparison") as final,
    ):
        control.main()
    assert json.loads((tmp_path / "status.json").read_text())["status"] == "completed"
    calls = sync.call_args_list
    assert "--exclude=request.json" in calls[1].args[2]
    assert "--include=request.json" in calls[2].args[2]
    final.assert_called_once()


def test_transfer_failure_is_recorded_and_raised(tmp_path: Path) -> None:
    """A failed transfer is not reported as successful campaign completion."""
    configure_main(tmp_path)
    with (
        patch("sys.argv", ["control", "--campaign", str(tmp_path)]),
        patch.object(control, "sync_tree", side_effect=RuntimeError("transfer failed")),
        pytest.raises(RuntimeError, match="transfer failed"),
    ):
        control.main()
    assert json.loads((tmp_path / "status.json").read_text())["status"] == "failed"


def test_final_comparison_retains_all_models_and_stages(tmp_path: Path) -> None:
    """All45 seeds yield90 stage records; sites are not substituted for replicates."""
    jobs = [job(i) for i in range(45)]
    for item in jobs:
        folder = (
            tmp_path
            / "models/WT"
            / f"seed_{item['model_index'] + 1:02d}"
            / "relaxation"
        )
        folder.mkdir(parents=True)
        (folder / "acceptance.json").write_text(
            json.dumps({"force_converged": True, "final_rms_gradient": 0.04})
        )

    def compare(entries: list, reference: Path, output: Path) -> dict:
        assert len(entries) == 90 and sum(e["stage"] == "raw" for e in entries) == 45
        assert all(e["ligand_count"] == 4 for e in entries)
        output.mkdir()
        return {"model_count": 90}

    with (
        patch.object(control, "compare_models", side_effect=compare),
        patch.object(control, "plot_comparison") as plot,
    ):
        control.final_comparison(
            tmp_path, {"jobs": jobs}, {"reference_cif": "reference.cif"}
        )
    assert (
        len(
            list(
                csv.DictReader((tmp_path / "analysis/relaxation_acceptance.csv").open())
            )
        )
        == 45
    )
    assert (
        json.loads((tmp_path / "completion.json").read_text())["planned_models"] == 45
    )
    plot.assert_called_once()


def test_sync_uses_authorized_identity_and_nonzero_failure_gate() -> None:
    """All transfers use explicit nondeleting argv and checked process outcomes."""
    with patch.object(control.subprocess, "run") as run:
        control.sync_tree("source", "destination", ["--exclude=*"], "/authorized/key")
    argv = run.call_args.args[0]
    assert "--delete" not in argv and "/authorized/key" in argv[4]
    assert run.call_args.kwargs["check"] is True
