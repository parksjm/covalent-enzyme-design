"""Verify complete hydration-movie transfer without full-system MD data."""

import hashlib
from pathlib import Path
import shlex
import shutil
import subprocess

import pytest


def test_hydration_download_includes_complete_movie(tmp_path: Path) -> None:
    """Exercise actual rsync rules, paired viewer inputs and checksum validation."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("Actual rsync is required to verify portable downloads")
    text = Path("scripts/pull_nyl12_inspection.sh").read_text()
    block = text.split("# Deacylation hydration:", 1)[1].split("# Flipped model30:", 1)[
        0
    ]
    rules = shlex.split(
        block.split('rsync "${hydration_rsync_args[@]}"', 1)[1]
        .split('"$remote:$hydration_source/"', 1)[0]
        .replace("\\\n", "")
    )
    source, target = tmp_path / "source", tmp_path / "laptop"
    portable = [
        "view_hydration.vmd",
        "hydration.dcd",
        "stripped.parm7",
        "qm_indices.tcl",
        "README.md",
        "frame_manifest.json",
        "export.cpptraj",
    ]
    keep = [f"vmd_progress/{name}" for name in portable] + [
        "analysis/summary.json",
        "analysis/per_frame_metrics.csv",
        "analysis/water_approach.png",
        "README.md",
        "hydration_endpoint_status_20261124.json",
    ]
    reject = [
        "heat/trajectory.nc",
        "equil/trajectory.nc",
        "extend_01/trajectory.nc",
        "inputs/system.prmtop",
        "inputs/endpoint.ncrst",
        "job/job.log",
        "analysis/failed_initial_20261124/report.json",
        "analysis/per_frame_metrics_initial_20261124.csv",
        "vmd_progress/clutter/old_movie.dcd",
    ]
    for name in keep + reject:
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    (source / "vmd_progress/MANIFEST.sha256").write_text(
        "".join(
            f"{hashlib.sha256((source / 'vmd_progress' / name).read_bytes()).hexdigest()}  {name}\n"
            for name in portable
        )
    )
    subprocess.run([rsync, "-rt", *rules, f"{source}/", f"{target}/"], check=True)
    assert all((target / name).is_file() for name in keep)
    assert all(not (target / name).exists() for name in reject)
    for line in (target / "vmd_progress/MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        actual = hashlib.sha256((target / "vmd_progress" / name).read_bytes())
        assert actual.hexdigest() == expected
    assert '"${checksum[@]}" MANIFEST.sha256' in block
    assert "--exclude='/13_deacylation_pilot/'" in text
