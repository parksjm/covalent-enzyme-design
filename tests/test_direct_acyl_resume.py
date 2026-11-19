"""Check explicit, archive-preserving direct-acyl continuation preparation."""

import argparse
import hashlib
import json
from pathlib import Path

import pytest

from nyl12_covalent.direct_acyl_resume import prepare_resume


def _setup(tmp_path: Path) -> argparse.Namespace:
    output = tmp_path / "work"
    drive = output / "01_coupled_coordinate_drive"
    for label in ("asynchronous_20", "asynchronous_21"):
        folder = drive / label
        folder.mkdir(parents=True)
        (folder / "optimized.ncrst").write_text(label)
        (folder / "acceptance.json").write_text(
            json.dumps(
                {
                    "native_converged": True,
                    "chemistry_retained": label.endswith("20"),
                }
            )
        )
    (tmp_path / "source_snapshot").mkdir()
    (tmp_path / "source_snapshot" / "original.py").write_text("immutable source")
    rows = [
        {"label": f"asynchronous_{i}", "composite_target_A": i / 10}
        for i in (20, 21, 22)
    ]
    (output / "inputs.json").write_text(json.dumps({"schedule": rows}))
    (output / "status.json").write_text(
        json.dumps(
            {
                "state": "stopped",
                "completed_windows": [r["label"] for r in rows[:2]],
            }
        )
    )
    start = drive / "asynchronous_20" / "optimized.ncrst"
    return argparse.Namespace(
        output=output,
        start=start,
        resume_after="asynchronous_20",
        checkpoint_sha256=hashlib.sha256(start.read_bytes()).hexdigest(),
        final_composite=2.2,
        substrate_proton_guard=True,
    )


def test_resume_preserves_prefix_and_archives_rejected_suffix(tmp_path: Path) -> None:
    args = _setup(tmp_path)
    original = args.start.read_bytes()
    rows, status = prepare_resume(args)
    assert [r["label"] for r in rows] == ["asynchronous_21", "asynchronous_22"]
    assert status["completed_windows"] == ["asynchronous_20"]
    assert args.start.read_bytes() == original
    archive = tmp_path / "clutter" / "initial_window21_proton_detour"
    assert (
        archive / "asynchronous_21" / "optimized.ncrst"
    ).read_text() == "asynchronous_21"
    assert (archive / "source_snapshot" / "original.py").is_file()
    assert not (args.start.parent.parent / "asynchronous_21").exists()
    assert status["continuation"]["substrate_proton_guard"]
    assert status["continuation"]["rejected_suffix_excluded_from_continuous_path"]


@pytest.mark.parametrize(
    "invalid", ["digest", "start", "state", "prefix", "gate", "final", "archive", "end"]
)
def test_invalid_resume_does_not_move_files(tmp_path: Path, invalid: str) -> None:
    args = _setup(tmp_path)
    if invalid == "digest":
        args.checkpoint_sha256 = "wrong"
    elif invalid == "start":
        args.start = args.start.parent.parent / "asynchronous_21" / "optimized.ncrst"
    elif invalid == "state":
        (args.output / "status.json").write_text('{"state":"running"}')
    elif invalid == "prefix":
        (args.output / "status.json").write_text(
            '{"state":"stopped","completed_windows":[]}'
        )
    elif invalid == "gate":
        (args.start.parent / "acceptance.json").write_text('{"native_converged":false}')
    elif invalid == "final":
        args.final_composite = 99
    elif invalid == "archive":
        (tmp_path / "clutter" / "initial_window21_proton_detour").mkdir(parents=True)
    else:
        args.resume_after = "asynchronous_22"
        args.start = args.start.parent.parent / args.resume_after / "optimized.ncrst"
        args.start.parent.mkdir()
        args.start.write_text("end")
        args.checkpoint_sha256 = hashlib.sha256(args.start.read_bytes()).hexdigest()
    with pytest.raises((ValueError, FileExistsError)):
        prepare_resume(args)
    assert (args.output / "01_coupled_coordinate_drive" / "asynchronous_21").is_dir()
