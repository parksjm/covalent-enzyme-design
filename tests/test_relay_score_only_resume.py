"""Score continuation preserves the pool and reuses matched WT orders."""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent.relay_broad_runtime import rank_broad

SCRIPT = Path(__file__).resolve().parents[1] / (
    "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/"
    "protection_corrected_pilot_20261127/scoring_recovery_attempt01/resume_scores.py"
)
SPEC = importlib.util.spec_from_file_location("score_only_test", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
_ORIGINAL_PATH = sys.path.copy()
try:
    SPEC.loader.exec_module(MODULE)
finally:
    sys.path[:] = _ORIGINAL_PATH
assert MODULE.rank_broad is rank_broad


def _fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict:
    original = tmp_path / "generation"
    original.mkdir()
    reference = {
        **dict.fromkeys("ABCD", "A" * 155 + "K" + "A" * 76),
        **dict.fromkeys("EFGH", "A" * 96),
    }
    monkeypatch.setattr(MODULE, "_pdb_sequences", lambda _: reference)
    monkeypatch.setattr(
        MODULE, "load_runner", lambda _: SimpleNamespace(flags=lambda *_: [])
    )
    monkeypatch.setattr(MODULE, "native_run", lambda *_: 0)
    monkeypatch.setattr(
        MODULE,
        "describe_stats",
        lambda *_: ({"local_AE_design_nll": 1.0}, np.zeros((2, 1312))),
    )
    candidates = []
    for i in range(24):
        sequence = dict(reference)
        for chain in "ABCD":
            value = list(sequence[chain])
            value[i + 1] = "C"
            value[159] = "E"
            sequence[chain] = "".join(value)
        candidates.append(
            {
                "candidate_id": f"relay_corr01_{i + 1:02d}",
                "sequence": ":".join(sequence.values()),
                "mutation_count": 2,
            }
        )
    (original / "deduplicated_pool.json").write_text(json.dumps(candidates))
    manifest = {
        "native": {"python": "native", "score_engine": "engine", "environment": {}},
        "request": {
            "generation": {"mask": {"alpha": [*range(2, 26), 160], "beta": []}}
        },
        "local_pocket_mask": {"alpha": [160], "beta": []},
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    templates = []
    for i in range(1, 6):
        (original / f"score_WT_{i:02d}.process.json").write_text('{"exit_code":0}')
        templates.append(
            {
                "id": str(i),
                "state": "TS" if i < 5 else "ES",
                "pdb": str(tmp_path / f"template{i}.pdb"),
                "WT_stats": "stats",
            }
        )
    return {
        "sha256": {
            str(manifest_path): hashlib.sha256(manifest_path.read_bytes()).hexdigest()
        },
        "raw_count": 24,
        "unique_count": 24,
        "generation_manifest": str(manifest_path),
        "generation_output": str(original),
        "templates": templates,
        "coherent_score_wrapper": "wrapper",
        "broad_runner": "runner",
    }


def test_complete_resume_reuses_WT_and_scores_all_24(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _fixture(tmp_path, monkeypatch)
    before = (Path(config["generation_output"]) / "deduplicated_pool.json").read_bytes()
    MODULE.run(config, tmp_path / "scores")
    acceptance = json.loads(
        (tmp_path / "scores/scientific_acceptance.json").read_text()
    )
    assert acceptance["native_candidate_scores"] == 120
    assert acceptance["generation_was_not_repeated"] is True
    assert len(list((tmp_path / "scores").glob("completed_*.json"))) == 24
    assert (
        Path(config["generation_output"]) / "deduplicated_pool.json"
    ).read_bytes() == before


@pytest.mark.parametrize(
    "failure",
    ["source", "WT", "native", "order", "count", "states", "ids", "mutations"],
)
def test_invalid_recovery_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    config = _fixture(tmp_path, monkeypatch)
    original = Path(config["generation_output"])
    if failure == "source":
        config["sha256"][config["generation_manifest"]] = "wrong"
    elif failure == "WT":
        (original / "score_WT_01.process.json").write_text('{"exit_code":1}')
    elif failure == "native":
        monkeypatch.setattr(MODULE, "native_run", lambda *_: 1)
    elif failure == "order":
        counter = iter(range(125))
        monkeypatch.setattr(
            MODULE,
            "describe_stats",
            lambda *_: (
                {"local_AE_design_nll": 1.0},
                np.full((2, 1312), next(counter)),
            ),
        )
    elif failure == "count":
        config["raw_count"] = 23
    elif failure == "states":
        config["templates"][0]["state"] = "ES"
    else:
        path = original / "deduplicated_pool.json"
        data = json.loads(path.read_text())
        data[0]["candidate_id" if failure == "ids" else "mutation_count"] = (
            "old" if failure == "ids" else 0
        )
        path.write_text(json.dumps(data))
    with pytest.raises((ValueError, RuntimeError)):
        MODULE.run(config, tmp_path / "scores")
    assert not (tmp_path / "scores/scientific_acceptance.json").exists()
