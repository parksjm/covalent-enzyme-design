"""Check that relay US downloads contain compact results without raw sampling."""

from pathlib import Path
import shlex
import shutil
import subprocess

import pytest


def test_relay_us_export_keeps_summaries_without_raw_windows(tmp_path: Path) -> None:
    """Exercise the actual ordered rsync rules against representative artifacts."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("Actual rsync is needed to check ordered include rules")
    text = Path("scripts/pull_nyl12_inspection.sh").read_text()
    block = text.split("# N-terminal relay:", 1)[1].split(
        "# Deacylation hydration:", 1
    )[0]
    options = shlex.split(
        block.split('rsync "${nterminal_relay_rsync_args[@]}"', 1)[1]
        .split('"$remote:$nterminal_relay_source/"', 1)[0]
        .replace("\\\n", "")
    )
    campaign = "umbrella_sampling/full_path_300K/"
    keep = [
        "prepared.json",
        "prepared.sha256",
        "preflight.json",
        "launch_authorization.json",
        "status.json",
        "README.md",
        "relay_us_azalea_20261125/runtime_summary.json",
        "relay_us_azalea_20261125/plan.json",
        "relay_us_azalea_20261125/launch_argv.json",
        "analysis/readiness_20261125.json",
        "analysis/readiness_20261125.png",
        "analysis/first_coverage_20261125/checkpoint_001/analysis.json",
        "analysis/first_coverage_20261125/checkpoint_001/analysis.png",
        "analysis/first_coverage_20261125/checkpoint_001/analysis_request.json",
    ]
    reject = [
        "inputs/system.prmtop",
        "inputs/seeds/source.ncrst",
        "source_snapshot/runtime.py",
        "activation_full_001/03_production_5ps/trajectory.nc",
        "activation_full_001/03_production_5ps/analysis_acceptance.json",
        "activation_full_001/03_production_5ps/receipt.json",
        "analysis/first_coverage_20261125/checkpoint_001/source_rows.json",
        "analysis/first_coverage_20261125/checkpoint_001/native.mdout",
        "clutter/old.json",
        "relay_us_azalea_20261125/source_snapshot/runtime.py",
        "relay_us_azalea_20261125/workers/worker_00/receipt.json",
    ]
    source, target = tmp_path / "source", tmp_path / "laptop"
    for name in keep + reject:
        path = source / campaign / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    subprocess.run([rsync, "-rt", *options, f"{source}/", f"{target}/"], check=True)
    assert all((target / campaign / name).is_file() for name in keep)
    assert all(not (target / campaign / name).exists() for name in reject)


def test_current_movie_uses_immutable_checksummed_publication() -> None:
    """Use the verified immutable bundle for the mutable current movie path."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("Actual rsync is needed to exercise trajectory transfer args")
    script = Path("scripts/pull_nyl12_inspection.sh").read_text()
    transfer = script.split("relay_current_rsync_args=(", 1)[1].split(")", 1)[0]
    options = shlex.split(transfer)
    assert "--checksum" in options
    assert "current_us_movie_rsync_args" not in script
    assert 'movie_source=$(ssh "$remote" readlink -f "$movie_pointer")' in script
    assert '"$remote:$movie_source/" "$movie_destination/"' in script
