"""Check the actual new-US rsync rules against compact and raw fixture files."""

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
    options = shlex.split(
        block.split('rsync "${vmd_rsync_args[@]}"', 1)[1]
        .split('"$remote:$trace_us_source/"', 1)[0]
        .replace("\\\n", "")
    )
    source, target = tmp_path / "source", tmp_path / "target"
    keep = [
        "prepared.json",
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
    (source / "analysis/latest").symlink_to("update_000002")
    subprocess.run([rsync, "-rt", *options, f"{source}/", f"{target}/"], check=True)
    for name in keep:
        assert (target / name).is_file(), name
    assert (target / "analysis/latest/coverage.png").is_file()
    for name in reject + ["analysis/update_000002/coverage.png"]:
        assert not (target / name).exists(), name
