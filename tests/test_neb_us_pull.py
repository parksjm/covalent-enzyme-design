"""Check the actual new-US rsync rules against compact and raw fixture files."""

import hashlib
from pathlib import Path
import shlex
import shutil
import subprocess

import pytest


def test_neb_us_pull_keeps_current_diagnostics_and_complete_views(
    tmp_path: Path,
) -> None:
    """Keep current plots/portable views while excluding histories and raw inputs."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("rsync required for actual filter-order verification")
    script = Path("scripts/pull_nyl12_inspection.sh").read_text()
    block = script.split("# NEB-seeded two-coordinate US:", 1)[1].split(
        "# Cumulative accepted frames", 1
    )[0]
    size_options = shlex.split(
        block.split("trace_us_rsync_args=(", 1)[1].split(")", 1)[0]
    )
    assert "--max-size=500m" in size_options
    options = shlex.split(
        block.split('rsync "${trace_us_rsync_args[@]}"', 1)[1]
        .split('"$remote:$trace_us_source/"', 1)[0]
        .replace("\\\n", "")
    )
    source, target = tmp_path / "source", tmp_path / "target"
    keep = [
        "prepared.json",
        "vmd_dynamic/path.parm7",
        "vmd_dynamic/path.dcd",
        "vmd_dynamic/view_umbrella.vmd",
        "vmd_dynamic/manifest.json",
        "vmd_dynamic/frame_manifest.json",
        "vmd_dynamic/frame_manifest.csv",
        "vmd_dynamic/request.json",
        "vmd_dynamic/README.md",
        "vmd_dynamic/MANIFEST.sha256",
        "analysis/pmf_latest/pmf_1d.png",
        "analysis/pmf_latest/pmf_1d.csv",
        "analysis/pmf_latest/README.md",
        "source_snapshot/scripts/run_neb_umbrella.py",
        "dispatch/azalea_retry02/source_snapshot/src/nyl12_covalent/jobs.py",
        "dispatch/sync_source_snapshot/scripts/sync_neb_umbrella_results.py",
        "dispatch/sync_source_snapshot/source_sha256.json",
        "dispatch/plan.json",
        "dispatch/azalea_retry02/plan.json",
        "dispatch/host_manifests/azalea_prepared.json",
        "azalea_compatibility_receipt.json",
        "workers/worker_01/status.json",
        "window_001/01_heating_1ps/run.mdin",
        "window_001/03_production_block_01_5ps/acceptance.json",
        "window_001/progress/block_01/view.pdb",
        "window_001/progress/block_01/view.vmd",
        "window_001/progress/block_01/MANIFEST.sha256",
    ]
    reject = [
        "inputs/system.prmtop",
        "vmd_dynamic/work/azalea/part.dcd",
        "vmd_dynamic/work/violet/part.dcd",
        "source_snapshot/__pycache__/cached.pyc",
        "dispatch/sync_source_snapshot/unused.bin",
        "execution_preflight/azalea/native_probe_attempt1/process.json",
        "dispatch/preparation_history/ncurses_attempt1/window_jobs/window_053/status.json",
        "dispatch/incoming/azalea/window_053/observables.csv",
        "dispatch/azalea_retry02/source_snapshot/unused.json",
        "dispatch/sync_source_snapshot/unused.json",
        "dispatch/azalea/worker01/job/command.json",
        "dispatch/sync_job/command.json",
        "source_snapshot/source.json",
        "analysis_job/source_snapshot/state.json",
        "preparation_history/old/prepared.json",
        "native_boundary_probe_attempt2/verification.json",
        "analysis/update_000001/coverage.png",
        "window_001/03_production_block_01_5ps/trajectory.nc",
        "window_001/03_production_block_01_5ps/restart.ncrst",
        "window_001/03_production_block_01_5ps/run.mdout",
    ]
    for name in keep + reject + ["analysis/update_000002/coverage.png"]:
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    large_trajectory = source / "vmd_dynamic/path.dcd"
    with large_trajectory.open("r+b") as trajectory:
        trajectory.truncate(100 * 1024 * 1024 + 1)
    manifest_entries = []
    for filename in (
        "path.parm7",
        "path.dcd",
        "view_umbrella.vmd",
        "manifest.json",
        "frame_manifest.json",
        "frame_manifest.csv",
        "request.json",
        "README.md",
    ):
        digest = hashlib.sha256()
        with (source / "vmd_dynamic" / filename).open("rb") as asset:
            for chunk in iter(lambda: asset.read(1024 * 1024), b""):
                digest.update(chunk)
        manifest_entries.append(f"{digest.hexdigest()}  {filename}")
    (source / "vmd_dynamic/MANIFEST.sha256").write_text(
        "\n".join(manifest_entries) + "\n"
    )
    (source / "analysis/latest").symlink_to("update_000002")
    subprocess.run(
        [rsync, "-rtS", *size_options, *options, f"{source}/", f"{target}/"],
        check=True,
    )
    for name in keep:
        assert (target / name).is_file(), name
    assert (target / "vmd_dynamic/path.dcd").stat().st_size == 100 * 1024 * 1024 + 1
    for entry in (target / "vmd_dynamic/MANIFEST.sha256").read_text().splitlines():
        expected_digest, filename = entry.split("  ", 1)
        digest = hashlib.sha256()
        with (target / "vmd_dynamic" / filename).open("rb") as asset:
            for chunk in iter(lambda: asset.read(1024 * 1024), b""):
                digest.update(chunk)
        assert digest.hexdigest() == expected_digest, filename
    assert (target / "analysis/latest/coverage.png").is_file()
    for name in reject + ["analysis/update_000002/coverage.png"]:
        assert not (target / name).exists(), name


def test_model30_pull_excludes_fine_grained_hessian_columns(
    tmp_path: Path,
) -> None:
    """Retain compact TS validation artifacts while omitting column receipts."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("rsync required for actual filter-order verification")
    script = Path("scripts/pull_nyl12_inspection.sh").read_text()
    block = script.split("# Flipped model30:", 1)[1].split(
        "# One compact 99-frame endpoint-minimization movie", 1
    )[0]
    options = shlex.split(
        block.split('rsync "${vmd_rsync_args[@]}" --copy-links', 1)[1]
        .split(
            '"$remote:/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/"', 1
        )[0]
        .replace("\\\n", "")
    )
    source, target = tmp_path / "source", tmp_path / "target"
    root = "05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/ts_hessian"
    keep = [
        f"{root}/summary.json",
        f"{root}/saddle_verification.json",
        f"{root}/hessian_modes.png",
        f"{root}/view_hessian.vmd",
        f"{root}/run.sh",
        f"{root}/min.in",
    ]
    reject = [
        f"{root}/columns/column_01445/minus/receipt.json",
        f"{root}/columns/column_01445/minus/displacement.xyz",
        f"{root}/columns/column_01446/plus/receipt.json",
    ]
    for name in keep + reject:
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    subprocess.run([rsync, "-rt", *options, f"{source}/", f"{target}/"], check=True)
    for name in keep:
        assert (target / name).is_file(), name
    for name in reject:
        assert not (target / name).exists(), name
