"""Contract tests for the isolated upstream scoring adapter and sidecar."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from nyl12_covalent import ipsae2_scoring as scoring
from nyl12_covalent import ipsae2_sidecar as sidecar


@pytest.fixture
def inputs(tmp_path: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Create minimal hashed native files and an isolated executable pin."""
    artifacts = {}
    for name in ("cif", "pae", "confidence"):
        path = tmp_path / name
        path.write_text(name)
        artifacts[name] = {"path": str(path), "sha256": scoring.file_sha256(path)}
    executable = tmp_path / "ipsae"
    executable.write_text("isolated command")
    runtime = {
        "command": [str(executable)],
        "revision": scoring.REVISION,
        "version": "0.1.11",
        "source_sha256": {str(executable): scoring.file_sha256(executable)},
    }
    model = {"case_id": "WT_seed_1", "artifacts": artifacts}
    mapping = {
        "accepted": True,
        "artifacts": artifacts,
        "protein_chains": list("ABCDEFGH"),
        "ligand_chains": list("IJKL"),
        "protomer_groups": [list(pair) for pair in ("AE", "BF", "CG", "DH")],
        "physical_pairs_verified": True,
    }
    return model, mapping, runtime


def fake_run(command: list[str], **kwargs: Any) -> SimpleNamespace:
    """Write realistic opposite-direction ligand rows without upstream math."""
    path = Path(command[command.index("--output-file") + 1])
    columns = [
        "chn1",
        "chn2",
        "ipsae",
        "ipsae12",
        "ipsae21",
        "n_scored_contact12",
        "n_scored_contact21",
        "align_token12",
        "n0_12",
        "d0_12",
    ]
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for index, ligand in enumerate("IJKL"):
            left, right = (ligand, "A") if index == 1 else ("A", ligand)
            writer.writerow(
                dict(
                    zip(
                        columns,
                        [left, right, "0.9", "0.2", "0.9", "0", "0", "A:1", "33", "2"],
                        strict=True,
                    )
                )
            )
        writer.writerow(
            dict(
                zip(
                    columns,
                    ["A", "B", "NA", "NA", "0.4", "NA", "NA", "NA", "NA", "NA"],
                    strict=True,
                )
            )
        )
    return SimpleNamespace(returncode=0)


def test_all_sites_and_explicit_direction(
    inputs: tuple, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    model, mapping, runtime = inputs
    monkeypatch.setattr(scoring.subprocess, "run", fake_run)
    result = scoring.score_model_raw(model, mapping, runtime, tmp_path / "scores")
    assert result["status"] == "completed"
    assert set(result["ligand_profiles"]) == set("IJKL")
    assert all(len(rows) == 3 for rows in result["ligand_profiles"].values())
    first = result["ligand_profiles"]["I"][0]
    assert first["ipsae"] == 0.9
    assert first["protein_aligned_ipsae"] == 0.2
    assert first["no_scored_contact"] is True
    assert result["ligand_profiles"]["J"][0]["protein_aligned_direction"] == "21"
    ambiguous = result["runs"]["chain_pairs"]["rows"][-1]
    assert ambiguous["ipsae"] is None and ambiguous["no_scored_contact"] is None
    assert "not_distinguished" in ambiguous["undefined_reason"]
    process = json.loads((tmp_path / "scores/whole_protein_process.json").read_text())
    assert "--all-pairs" in process["command"]
    assert "A+B+C+D+E+F+G+H" in process["command"]


@pytest.mark.parametrize(
    "target", ["runtime", "native", "mapping", "groups", "chains", "accepted"]
)
def test_fail_closed(inputs: tuple, tmp_path: Path, target: str) -> None:
    model, mapping, runtime = inputs
    if target == "runtime":
        Path(runtime["command"][0]).write_text("changed")
    elif target == "native":
        Path(model["artifacts"]["pae"]["path"]).write_text("changed")
    elif target == "mapping":
        mapping["artifacts"] = {**mapping["artifacts"], "pae": {"sha256": "0" * 64}}
    elif target == "groups":
        mapping["physical_pairs_verified"] = False
    elif target == "chains":
        mapping["protein_chains"] = ["A"]
    else:
        mapping["accepted"] = False
    with pytest.raises(ValueError):
        scoring.score_model_raw(model, mapping, runtime, tmp_path / "scores")
    assert not (tmp_path / "scores").exists()


def test_qualification_gate(inputs: tuple) -> None:
    runtime = inputs[2]
    for bad in (
        {},
        {"status": "accepted"},
        {"status": "accepted", "runtime_sha256": scoring.runtime_digest(runtime)},
    ):
        with pytest.raises(ValueError):
            scoring.require_qualification(bad, runtime)
    scoring.require_qualification(
        {
            "status": "accepted",
            "runtime_sha256": scoring.runtime_digest(runtime),
            "qualified_roles": ["WT", "design"],
        },
        runtime,
    )


def test_failed_process_retained(
    inputs: tuple, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        scoring.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=7)
    )
    with pytest.raises(RuntimeError):
        scoring.score_model_raw(*inputs, tmp_path / "failed")
    result = json.loads((tmp_path / "failed/scores.json").read_text())
    assert result["status"] == "failed"
    assert (tmp_path / "failed/chain_pairs.log").exists()


def test_runtime_validation(inputs: tuple) -> None:
    runtime = inputs[2]
    for patch in (
        {"version": "0.1.12"},
        {"source_sha256": {}},
        {"command": ["/another/executable"]},
    ):
        with pytest.raises(ValueError):
            scoring.verify_runtime({**runtime, **patch})


def test_invalid_upstream_scores(tmp_path: Path) -> None:
    path = tmp_path / "empty.tsv"
    path.write_text("chn1\tchn2\tipsae\n")
    with pytest.raises(ValueError, match="No score rows"):
        scoring.parse_scores(path, list("ABCDEFGH"), list("IJKL"))
    with pytest.raises(ValueError, match="Nonfinite"):
        scoring._number("inf")


@pytest.fixture
def config(inputs: tuple, tmp_path: Path) -> dict[str, Any]:
    model, mapping, runtime = inputs
    for name, value in (("runtime", runtime), ("model", model), ("mapping", mapping)):
        (tmp_path / f"{name}.json").write_text(json.dumps(value))
    qualification = {
        "status": "accepted",
        "runtime_sha256": scoring.runtime_digest(runtime),
        "qualified_roles": ["WT", "design"],
    }
    (tmp_path / "qualification.json").write_text(json.dumps(qualification))
    return {
        "runtime": str(tmp_path / "runtime.json"),
        "qualification": str(tmp_path / "qualification.json"),
        "scoring_root": str(tmp_path / "scoring"),
        "journal_command": ["journal"],
        "cases": [
            {
                "case_id": model["case_id"],
                "native_receipt": str(tmp_path / "native.json"),
                "model_manifest": str(tmp_path / "model.json"),
                "mapping_receipt": str(tmp_path / "mapping.json"),
            }
        ],
    }


def test_pending_and_plan(config: dict, tmp_path: Path) -> None:
    case_id = config["cases"][0]["case_id"]
    assert (
        sidecar.poll_cases(config)["cases"][case_id]["status"]
        == "pending_native_prediction"
    )
    native = tmp_path / "native.json"
    native.write_text(json.dumps({"status": "running"}))
    assert (
        sidecar.poll_cases(config)["cases"][case_id]["status"]
        == "pending_native_prediction"
    )
    native.write_text(
        json.dumps({"status": "native_prediction_completed", "case": case_id})
    )
    assert sidecar.poll_cases(config, plan=True)["cases"][case_id]["status"] == "ready"
    (tmp_path / "mapping.json").unlink()
    assert sidecar.poll_cases(config)["cases"][case_id]["status"] == "pending_mapping"
    (tmp_path / "qualification.json").unlink()
    assert sidecar.poll_cases(config)["status"] == "pending_qualification"
    assert not (tmp_path / "scoring").exists()


def test_sidecar_failure_does_not_block_later_case(
    config: dict, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    case_id = config["cases"][0]["case_id"]
    (tmp_path / "native.json").write_text(
        json.dumps({"status": "native_prediction_completed", "case": case_id})
    )
    config["cases"].append(
        {
            **config["cases"][0],
            "case_id": "design_seed_1",
            "native_receipt": str(tmp_path / "later.json"),
        }
    )
    monkeypatch.setattr(sidecar, "_journal", lambda *a: None)
    monkeypatch.setattr(
        sidecar,
        "score_model_raw",
        lambda *a: (_ for _ in ()).throw(ValueError("bad mapping")),
    )
    result = sidecar.poll_cases(config)
    assert result["cases"][case_id]["status"] == "failed"
    assert result["cases"]["design_seed_1"]["status"] == "pending_native_prediction"


def test_sidecar_success_and_cached_result(
    config: dict, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    case_id = config["cases"][0]["case_id"]
    (tmp_path / "native.json").write_text(
        json.dumps({"status": "native_prediction_completed", "case": case_id})
    )
    monkeypatch.setattr(sidecar, "_journal", lambda *a: None)
    monkeypatch.setattr(scoring.subprocess, "run", fake_run)
    assert sidecar.poll_cases(config)["cases"][case_id]["status"] == "completed"
    assert sidecar.poll_cases(config)["cases"][case_id]["status"] == "completed"


def test_journal_failure_stops_mutation(
    config: dict, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "native.json").write_text(
        json.dumps({"status": "native_prediction_completed"})
    )

    def fail_journal(*args: Any) -> None:
        raise OSError("journal unavailable")

    monkeypatch.setattr(sidecar, "_journal", fail_journal)
    with pytest.raises(OSError, match="journal unavailable"):
        sidecar.poll_cases(config)
    assert not (tmp_path / "scoring").exists()
