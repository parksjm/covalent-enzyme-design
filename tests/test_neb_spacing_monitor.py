"""Complete-report, incremental-stream and isolated-comparison contracts."""

import json
from pathlib import Path

import pytest

from nyl12_covalent.neb_spacing_monitor import chemical_gaps, parse_block, update_trial
from scripts.monitor_neb_spacing_tests import refresh


def native(count: int = 24, cycle: int = 50) -> bytes:
    """Return fixed-width native rows and a complete target-ratio test."""
    lines = ["NEB Report", "             Energy       F_tang    F_perp     Dist"]
    for i in range(1, count + 1):
        energy = -10 + (i - 1) / (count - 1) * 0.01
        lines.append(f"Img {i:4d}{energy:15.7f}{0.001:10.5f}{0.002:10.5f}{1.0:10.5f}")
    lines += [
        "Total path length 24.0",
        f"Testing convergence  in cycle  {cycle}",
        "  Max grad  8.0E-03 Target: 1.0E-04 converged? no",
        "  RMS grad  4.0E-04 Target: 2.0E-04 converged? no",
    ]
    return ("\n".join(lines) + "\n").encode()


def config(root: Path, count: int = 24) -> Path:
    """Create an explicit synthetic trial with fixed case identity."""
    trial = root / f"trial_{count}"
    (trial / "calculation").mkdir(parents=True)
    case = root / "case.json"
    case.write_text("{}")
    path = trial / "config.json"
    path.write_text(
        json.dumps(
            {"regular_images": count, "active_atoms": 4662, "case_path": str(case)}
        )
    )
    return path


@pytest.mark.parametrize("count", [24, 48])
def test_regular_count_physical_units_and_force_ratios(count: int) -> None:
    result = parse_block(native(count), count)
    assert len(result["images"]) == count
    assert result["max_gradient_ratio"] == pytest.approx(80.0)
    assert result["rms_gradient_ratio"] == pytest.approx(2.0)
    assert result["endpoint_kcal_mol"] == pytest.approx(6.275094740631)
    assert result["cycle"] == 50
    with pytest.raises(ValueError):
        parse_block(native(count), 48 if count == 24 else 24)
    with pytest.raises(ValueError):
        parse_block(native(count).replace(b"1.0E-04", b"0.0E-04"), count)


def test_incremental_partial_native_blocks_and_reset_timing(tmp_path: Path) -> None:
    path = config(tmp_path)
    source = path.parent / "calculation/minimized.mdout"
    output = tmp_path / "monitor"
    assert update_trial(path, output)["status"] == "waiting_for_native_output"
    first, second = native(cycle=50), native(cycle=100)
    source.write_bytes(first + b"Resetting L-BFGS memory\n" + second[:-20])
    status = update_trial(path, output, max_bytes=7)
    assert status["status"] == "waiting_for_complete_report"
    for _ in range(1500):
        status = update_trial(path, output, max_bytes=53)
        if status["bytes_consumed"] == source.stat().st_size:
            break
    current = json.loads((output / "current_profile.json").read_text())
    assert current["cycle"] == 50 and current["reset_count"] == 0
    assert len((output / "history.jsonl").read_text().splitlines()) == 1
    with source.open("ab") as stream:
        stream.write(second[-20:])
    status = update_trial(path, output)
    current = json.loads((output / "current_profile.json").read_text())
    assert current["cycle"] == 100 and current["reset_count"] == 1
    assert status["report_count"] == 2
    assert (output / "checkpoint_0050.json").is_file()
    assert (output / "checkpoint_0100.json").is_file()
    update_trial(path, output)
    assert len((output / "history.jsonl").read_text().splitlines()) == 2
    source.write_bytes(b"")
    with pytest.raises(ValueError, match="truncated"):
        update_trial(path, output)


def test_config_mutation_and_invalid_count(tmp_path: Path) -> None:
    path = config(tmp_path)
    source = path.parent / "calculation/minimized.mdout"
    source.write_bytes(native())
    output = tmp_path / "monitor"
    update_trial(path, output)
    with pytest.raises(ValueError, match="positive"):
        update_trial(path, output, max_bytes=0)
    cfg = json.loads(path.read_text())
    cfg["active_atoms"] = 1
    path.write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match="4662"):
        update_trial(path, output)
    cfg["active_atoms"] = 4662
    cfg["new_input"] = True
    path.write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match="changed"):
        update_trial(path, output)


def test_coherent_geometry_is_not_paired_to_another_cycle(tmp_path: Path) -> None:
    assert not chemical_gaps(tmp_path, 50)["available"]
    latest = tmp_path / "progress/latest"
    latest.mkdir(parents=True)
    (latest / "progress.json").write_text('{"cycle":50}')
    (latest / "images.csv").write_text("rOG1C,rCN\n1.0,1.5\n2.0,3.5\n")
    assert not chemical_gaps(tmp_path, 100)["available"]
    result = chemical_gaps(tmp_path, 50)
    assert result["max_adjacent_absolute_changes"]["rCN"] == 2.0


def test_explicit_startup_waiting_and_isolated_comparison_plot(tmp_path: Path) -> None:
    output = tmp_path / "plots"
    assert refresh(tmp_path, output)["status"] == "waiting_for_configs"
    path = config(tmp_path)
    (path.parent / "calculation/minimized.mdout").write_bytes(native())
    result = refresh(tmp_path, output)
    assert len(result["trials"]) == 1
    assert (output / "comparison.png").is_file()
    assert (output / "comparison_history.csv").is_file()
    stamp = (output / "comparison.png").stat().st_mtime_ns
    refresh(tmp_path, output)
    assert (output / "comparison.png").stat().st_mtime_ns == stamp


def test_image_resets_split_reads_and_legacy_cursor_migration(tmp_path: Path) -> None:
    """Count geometric events once across chunks, saved cursors and old schemas."""
    path = config(tmp_path)
    source = path.parent / "calculation/minimized.mdout"
    output = tmp_path / "monitor"
    first = native(cycle=1)
    source.write_bytes(first + b"  Resetting images  5 to 7\nResetting images 8 to")
    status = update_trial(path, output)
    assert status["image_reset_count"] == 1
    assert status["reset_count"] == status["bfgs_memory_reset_count"] == 0
    cursor_path = output / "cursor.json"
    cursor = json.loads(cursor_path.read_text())
    del cursor["image_reset_offsets"]  # An existing pre-patch persisted monitor.
    cursor_path.write_text(json.dumps(cursor))
    status = update_trial(path, output)
    assert status["image_reset_count"] == 1  # Migration excludes incomplete tail.
    with source.open("ab") as stream:
        stream.write(b" 8\nResetting L-BFGS memory\n" + native(cycle=2))
    status = update_trial(path, output, max_bytes=2)
    while status["bytes_consumed"] < source.stat().st_size:
        status = update_trial(path, output, max_bytes=71)
    profile = json.loads((output / "current_profile.json").read_text())
    assert profile["image_reset_count"] == 2  # Events, not3+1affected images.
    assert profile["bfgs_memory_reset_count"] == profile["reset_count"] == 1
    assert update_trial(path, output)["image_reset_count"] == 2
    assert len((output / "history.jsonl").read_text().splitlines()) == 2


def test_old_and_new_history_rows_export_together(tmp_path: Path) -> None:
    """A live monitor can append new counter columns to its legacy CSV history."""
    path = config(tmp_path)
    source = path.parent / "calculation/minimized.mdout"
    source.write_bytes(native(cycle=1))
    output = tmp_path / "plots"
    refresh(tmp_path, output)
    history = output / path.parent.name / "history.jsonl"
    old = json.loads(history.read_text())
    old.pop("image_reset_count")
    old.pop("bfgs_memory_reset_count")
    history.write_text(json.dumps(old) + "\n")
    with source.open("ab") as stream:
        stream.write(b"Resetting images 5 to 5\n" + native(cycle=2))
    refresh(tmp_path, output)
    assert "image_reset_count" in (output / "comparison_history.csv").read_text()
