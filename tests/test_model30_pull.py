"""Exercise model30 laptop filters with real local rsync transfers."""

import hashlib
import shlex
import shutil
import subprocess
from pathlib import Path

import pytest

PULL = Path(__file__).resolve().parents[1] / "scripts/pull_nyl12_inspection.sh"


def _model30_options() -> list[str]:
    text = PULL.read_text().split("# Flipped model30: compact coordinates", 1)[1]
    block = text.split('\nrsync "${vmd_rsync_args[@]}"', 1)[1]
    block = block.split(
        '\n    "$remote:/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/"',
        1,
    )[0]
    return shlex.split(block.replace("\\\n", "").rstrip(" \\\n"))


def test_model30_latest_bundle_is_complete_and_light(tmp_path: Path) -> None:
    """Dereference latest, retain its manifest, and skip historical/heavy data."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("Local rsync is required to test its filter semantics")
    source = tmp_path / "source"
    campaign = source / "05_neutral_reaction_path/04_direct_acylenzyme"
    progress = campaign / "work/vmd_progress"
    bundle = progress / "update_0002"
    bundle.mkdir(parents=True)
    names = (
        "manifest.json",
        "potential_energy.csv",
        "potential_energy.json",
        "potential_energy.png",
        "reaction_distances.json",
        "reaction_distances.png",
        "reaction_path.dcd",
        "reaction_path.vmd",
        "stripped.parm7",
        "topology.cpptraj",
        "topology_export.log",
        "trajectory.cpptraj",
        "trajectory_export.log",
    )
    for name in names:
        (bundle / name).write_text(name)
    (bundle / "MANIFEST.sha256").write_text(
        "".join(
            f"{hashlib.sha256(name.encode()).hexdigest()}  {name}\n" for name in names
        )
    )
    (progress / "latest").symlink_to(bundle.name, target_is_directory=True)
    included = ("job/run.sh", "inputs/base.mdin", "work/window/restraints.RST")
    excluded = (
        "source_snapshot/scripts/large.py",
        "clutter/old.json",
        "inputs/system.prmtop",
        "inputs/accepted_es.ncrst",
        "work/window/optimization.nc",
        "work/window/minimized.mdout",
    )
    for name in (*included, *excluded):
        path = campaign / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    destination = tmp_path / "destination"
    command = [rsync, "-rt", "--safe-links", *_model30_options()]
    dry = subprocess.run(
        [*command, "--dry-run", "--out-format=%n", f"{source}/", f"{destination}/"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    assert "vmd_progress/latest/reaction_path.dcd" in dry
    assert "vmd_progress/update_0002/" not in dry
    assert not destination.exists()
    subprocess.run([*command, f"{source}/", f"{destination}/"], check=True)
    copied = destination / campaign.relative_to(source)
    latest = copied / "work/vmd_progress/latest"
    assert latest.is_dir() and not latest.is_symlink()
    assert {p.name for p in latest.iterdir()} == {*names, "MANIFEST.sha256"}
    for line in (latest / "MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split("  ", 1)
        assert hashlib.sha256((latest / name).read_bytes()).hexdigest() == digest
    assert all((copied / name).is_file() for name in included)
    assert not any((copied / name).exists() for name in excluded)
    assert not (copied / "work/vmd_progress/update_0002").exists()


def test_model30_neutral_md_vmd_bundles_are_downloaded(tmp_path: Path) -> None:
    """Copy stripped neutral-replica movies but not raw production trajectories."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("Local rsync is required to test its filter semantics")
    source = tmp_path / "source"
    for replica in ("neutral_replica_01", "neutral_replica_02"):
        root = source / "03_md" / replica
        bundle = root / "vmd_stripped_100ns"
        bundle.mkdir(parents=True)
        names = (
            "README.md",
            "manifest.json",
            "production_100ns_stripped.dcd",
            "stripped.parm7",
            "topology.cpptraj",
            "topology_export.log",
            "trajectory.cpptraj",
            "trajectory_export.log",
            "view_production.vmd",
        )
        for name in names:
            (bundle / name).write_text(name)
        (bundle / "MANIFEST.sha256").write_text("manifest")
        (root / "05_production_10ns.nc").write_text("raw")
        (root / "06_production_90ns.nc").write_text("raw")
        (root / "system.prmtop").write_text("raw")
    destination = tmp_path / "destination"
    command = [rsync, "-rt", "--safe-links", *_model30_options()]
    subprocess.run([*command, f"{source}/", f"{destination}/"], check=True)
    for replica in ("neutral_replica_01", "neutral_replica_02"):
        copied = destination / "03_md" / replica
        bundle = copied / "vmd_stripped_100ns"
        assert (bundle / "production_100ns_stripped.dcd").is_file()
        assert (bundle / "view_production.vmd").is_file()
        assert (bundle / "trajectory_export.log").is_file()
        assert not (copied / "05_production_10ns.nc").exists()
        assert not (copied / "06_production_90ns.nc").exists()
        assert not (copied / "system.prmtop").exists()
