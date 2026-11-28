"""Movie coverage follows accepted source blocks rather than a fixed label."""

import json
from pathlib import Path
from typing import Any

import pytest

from nyl12_covalent.neb_us_movie_scope import validated_sampling_scope


def _fixture(tmp_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    prepared = {"physical_hamiltonian_sha256": "same", "windows": []}
    manifest = {"windows": ["a", "b"], "frame_plan": [], "frames": 16}
    for label in ("a", "b"):
        prepared["windows"].append(
            {
                "label": label,
                "stages": [{"kind": "production", "label": f"p{i}"} for i in range(4)],
            }
        )
        for i in range(4):
            folder = tmp_path / label / f"p{i}"
            folder.mkdir(parents=True)
            (folder / "analysis_acceptance.json").write_text(
                json.dumps(
                    {
                        "accepted": True,
                        "window": label,
                        "physical_hamiltonian_sha256": "same",
                        "sampled_duration_ps": 5.0,
                    }
                )
            )
            manifest["frame_plan"].extend([{"window": label, "block": f"p{i}"}] * 2)
    return prepared, manifest


def test_derives_twenty_ps_instead_of_hardcoded_five(tmp_path: Path) -> None:
    prepared, manifest = _fixture(tmp_path)
    scope = validated_sampling_scope(tmp_path, prepared, manifest)
    assert scope["sampled_horizon_per_window_ps"] == 20
    assert scope["sampled_frames_per_window"] == 8
    assert scope["production_blocks_per_window"] == 4
    assert len(scope["acceptance_receipt_sha256"]) == 8


@pytest.mark.parametrize(
    "error",
    [
        "order",
        "counts",
        "frames",
        "prefix",
        "accepted",
        "hamiltonian",
        "duration",
        "horizons",
    ],
)
def test_rejects_misleading_movie_scope(tmp_path: Path, error: str) -> None:
    prepared, manifest = _fixture(tmp_path)
    if error == "order":
        manifest["windows"].reverse()
    elif error == "counts":
        manifest["frame_plan"].pop()
    elif error == "frames":
        manifest["frames"] += 1
    elif error == "prefix":
        manifest["frame_plan"][0]["block"] = "p9"
    else:
        path = tmp_path / "a/p0/analysis_acceptance.json"
        receipt = json.loads(path.read_text())
        key, value = {
            "accepted": ("accepted", False),
            "hamiltonian": ("physical_hamiltonian_sha256", "different"),
            "duration": ("sampled_duration_ps", 0),
            "horizons": ("sampled_duration_ps", 3),
        }[error]
        receipt[key] = value
        path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError):
        validated_sampling_scope(tmp_path, prepared, manifest)
