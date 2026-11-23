"""Exercise real rsync filters for new progress-deck and relay bundles."""

from pathlib import Path
import shlex
import shutil
import subprocess

import pytest


def _options(start: str, end: str, variable: str, source: str) -> list[str]:
    script = Path("scripts/pull_nyl12_inspection.sh").read_text()
    block = script.split(start, 1)[1].split(end, 1)[0]
    return shlex.split(
        block.split(f'rsync "${{{variable}[@]}}"', 1)[1]
        .split(f'"$remote:${source}/"', 1)[0]
        .replace("\\\n", "")
    )


def test_progress_briefing_preserves_editable_deck_and_sources(tmp_path: Path) -> None:
    """Copy the deck, PDF and figure sources while omitting transient renders."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("rsync required for actual export filter checks")
    options = _options(
        "# Progress briefing:",
        "# Frozen native-Thr calibration:",
        "vmd_rsync_args",
        "briefing_source",
    )
    source, target = tmp_path / "source", tmp_path / "target"
    keep = [
        "progress.pptx",
        "progress.pdf",
        "slides.json",
        "README.md",
        "assets/pmf.png",
        "MANIFEST.sha256",
    ]
    reject = ["rendered/slide01.png", "raw_optimization.nc"]
    for name in keep + reject:
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    subprocess.run([rsync, "-rt", *options, f"{source}/", f"{target}/"], check=True)
    assert all((target / name).is_file() for name in keep)
    assert all(not (target / name).exists() for name in reject)


def test_relay_progress_keeps_complete_bundle_and_excludes_raw_data(
    tmp_path: Path,
) -> None:
    """Retain even source snapshots inside a verified checksummed movie bundle."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("rsync required for actual export filter checks")
    script = Path("scripts/pull_nyl12_inspection.sh").read_text()
    block = script.split("# Unbiased N-terminal relay releases:", 1)[1]
    assert "--max-size=500m" in block
    assert "flock --shared -- $relay_release_source/.vmd_export.lock rsync" in block
    options = _options(
        "# Unbiased N-terminal relay releases:",
        "END_NOT_PRESENT",
        "relay_release_rsync_args",
        "relay_release_source",
    )
    # --rsync-path is a remote-side option; the local filter test needs only rules.
    options = [option for option in options if not option.startswith("--rsync-path=")]
    source, target = tmp_path / "source", tmp_path / "target"
    prefix = "ti1_window032/vmd_progress_verified"
    keep = [
        "prepared.json",
        f"{prefix}/view_relay.vmd",
        f"{prefix}/optimization.dcd",
        f"{prefix}/stripped.parm7",
        f"{prefix}/MANIFEST.sha256",
        f"{prefix}/source_snapshot/frames.nc",
        "ti1_window033/minimum_verification/04_full_hessian/scientific_acceptance.json",
    ]
    reject = [
        "source_snapshot/model.py",
        "parameters_3ob/C-C.skf",
        "ti1_window032/inputs/seed.ncrst",
        "ti1_window032/optimization.nc",
        "ti1_window032/vmd_initial_progress/old.dcd",
        "ti1_window033/minimum_verification/04_full_hessian/columns/column_00001/plus/receipt.json",
    ]
    for name in keep + reject:
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    subprocess.run([rsync, "-rt", *options, f"{source}/", f"{target}/"], check=True)
    assert all((target / name).is_file() for name in keep)
    assert all(not (target / name).exists() for name in reject)
