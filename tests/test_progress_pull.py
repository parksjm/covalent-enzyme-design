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
    """Copy every checked deck asset while omitting obsolete sources and raw data."""
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
        "assets/pes.svg",
        "assets/pymol/render.py",
        "assets/pymol/states.pse",
        "assets/relay_candidate_white_greyH.png",
        "rendered/contact_sheet.png",
        "rendered/slide-01.png",
        "slide_list.txt",
        "MANIFEST.sha256",
    ]
    reject = [
        "source_history_before_mask_update/obsolete.pptx",
        "assets/mol_blender/failed_attempt/render.png",
        "assets/mol_blender/failed_attempt/style.json",
        "raw_optimization.nc",
    ]
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


def test_relay_cumulative_pmfs_include_extensions_without_raw_data(
    tmp_path: Path,
) -> None:
    """Export cumulative PMFs past 30 ps with their sources, but omit raw jobs."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("rsync required for actual export filter checks")
    options = _options(
        "# N-terminal relay: both branch movies,",
        "# Deacylation hydration:",
        "nterminal_relay_rsync_args",
        "nterminal_relay_source",
    )
    prefix = "umbrella_sampling/full_path_300K/analysis"
    progression = f"{prefix}/cumulative_progression_20261125"
    keep = [
        "umbrella_sampling/README.md",
        f"{prefix}/cumulative_mbar_progression_20261125.md",
        f"{prefix}/direct_relay_compatibility_20261125.json",
        f"{progression}/latest_progression.png",
        f"{progression}/latest_checkpoint.json",
        f"{progression}/checkpoint_05ps/analysis.png",
        f"{progression}/checkpoint_05ps/source.json",
        f"{progression}/checkpoint_35ps/analysis.json",
        f"{progression}/checkpoint_35ps/analysis_request.json",
        f"{progression}/checkpoint_35ps/progression.png",
        "umbrella_sampling/full_path_300K/vmd_progress/current/view_path.vmd",
        "umbrella_sampling/full_path_300K/vmd_progress/current/path.dcd",
        "umbrella_sampling/full_path_300K/vmd_progress/current/path.parm7",
        "umbrella_sampling/full_path_300K/vmd_progress/current/manifest.json",
        "umbrella_sampling/full_path_300K/vmd_progress/current/qm_indices.json",
        "umbrella_sampling/full_path_300K/vmd_progress/current/SHA256SUMS.txt",
        "umbrella_sampling/full_path_300K/vmd_progress/relay_us_movie_5ps_20261126/view_path.vmd",
        "umbrella_sampling/full_path_300K/plots/pmf_20ps_readable.png",
        "umbrella_sampling/full_path_300K/plots/mbar_overlap_20ps.pdf",
    ]
    reject = [
        f"{progression}/checkpoint_35ps/raw.nc",
        f"{progression}/checkpoint_35ps/source_snapshot/model.py",
        f"{prefix}/watcher_inputs_20261125_continuation03/source_snapshot/code.py",
        f"{prefix}/watcher_runtime_cumulative_progression_20261125/job.log",
        "umbrella_sampling/full_path_300K/movie_20ps_export_20261126/job/status.json",
        "umbrella_sampling/full_path_300K/activation_full_001/07_production_5ps/trajectory.nc",
    ]
    source, target = tmp_path / "source", tmp_path / "target"
    for name in keep + reject:
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    subprocess.run([rsync, "-rt", *options, f"{source}/", f"{target}/"], check=True)
    assert all((target / name).is_file() for name in keep)
    assert all(not (target / name).exists() for name in reject)
