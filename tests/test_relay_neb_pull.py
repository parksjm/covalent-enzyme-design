"""Exercise latest-only relay NEB laptop transfer using real rsync."""

from pathlib import Path
import hashlib
import shlex
import shutil
import subprocess

import pytest

from scripts import watch_relay_acyl_neb as monitor


def relay_options() -> list[str]:
    """Extract actual ordered relay rsync rules from the production pull script."""
    text = (monitor.ROOT / "scripts/pull_nyl12_inspection.sh").read_text()
    block = text.split("# Relay acylation campaign:", 1)[1].split(
        "# Final embedded-method", 1
    )[0]
    return shlex.split(
        block.split('rsync "${vmd_rsync_args[@]}"', 1)[1]
        .split('"$remote:$relay_campaign_source/"', 1)[0]
        .replace("\\\n", "")
    )


def test_relay_neb_pull_copies_only_latest_stripped_bundle(tmp_path: Path) -> None:
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("rsync is required to verify its actual filter ordering")
    options = relay_options()
    source = tmp_path / "source"
    target = tmp_path / "target"
    keep = [
        "04_neb/config.json",
        "04_neb/calculation/minimize.mdin",
        "04_neb/progress/update_000001/vmd/path.dcd",
        "04_neb/progress/update_000001/vmd/path.parm7",
        "04_neb/progress/update_000001/vmd/view_path.vmd",
        "04_neb/progress/update_000001/energy.svg",
        "04_neb/progress/update_000001/images.json",
    ]
    reject = [
        "04_neb/calculation/raw.json",
        "04_neb/progress/update_000001/source_nebpath.xyz.gz",
        "04_neb/progress/update_000001/native_report.mdout",
        "04_neb/progress/update_000001/optimized.0001.ncrst",
        "04_neb/progress/source_milestones/milestone_0001/source.py",
    ]
    for name in keep + reject:
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    (source / "04_neb/progress/latest").symlink_to("update_000001")
    subprocess.run(
        [rsync, "-rt", *options, str(source) + "/", str(target) + "/"], check=True
    )
    assert (target / "04_neb/config.json").is_file()
    for name in keep[2:]:
        assert (target / name.replace("update_000001", "latest")).is_file()
    assert not (target / "04_neb/progress/update_000001").exists()
    assert (target / "04_neb/calculation/minimize.mdin").is_file()
    assert not (target / "04_neb/calculation/raw.json").exists()
    assert not (target / "04_neb/progress/latest/source_nebpath.xyz.gz").exists()
    assert not (target / "04_neb/progress/latest/native_report.mdout").exists()
    assert not (target / "04_neb/progress/source_milestones").exists()


@pytest.mark.parametrize(
    "trial",
    [
        "01_uniform_24_lbfgs",
        "02_uniform_48_lbfgs",
        "03_adaptive_24_lbfgs",
        "04_uniform_24_cg",
    ],
)
def test_spacing_trial_pull_keeps_latest_portable_only(
    tmp_path: Path, trial: str
) -> None:
    """All trial latest manifests survive rsync while archives and raw inputs do not."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("rsync required to verify actual ordered filter semantics")
    source, target = tmp_path / "source", tmp_path / "target"
    relative = Path("04_neb/spacing_tests") / trial
    root = source / relative
    keep = [
        "config.json",
        "preflight.json",
        "run.sh",
        "job/status.json",
        "plots/neb_evolution.png",
        "plots/neb_evolution.svg",
        "plots/neb_evolution.csv",
        "plots/plot_energy_evolution.py",
        "calculation/minimize.mdin",
        "progress/progress.json",
    ]
    reject = [
        "calculation/native.json",
        "calculation/nebpath.xyz",
        "calculation/minimized.mdout",
        "calculation/optimized.0001.ncrst",
        "inputs/seed.ncrst",
        "inputs/manifest.json",
        "source_snapshot/source.py",
        "source_snapshots/manifest.json",
        "job/source_snapshot/source.py",
        "progress/source_milestones/milestone_0001/source.py",
        "history/old.json",
        "progress/history/old.json",
    ]
    for name in keep + reject:
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    bundle = root / "progress/update_000002"
    bundle.mkdir()
    payload = {
        "images.json": b"images",
        "energy.svg": b"svg",
        "vmd/path.dcd": b"dcd",
        "vmd/path.parm7": b"parm",
        "vmd/view_path.vmd": b"vmd",
    }
    for name, content in payload.items():
        path = bundle / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    vmd_names = [
        name.removeprefix("vmd/") for name in payload if name.startswith("vmd/")
    ]
    (bundle / "vmd/MANIFEST.sha256").write_text(
        "".join(
            f"{hashlib.sha256((bundle / 'vmd' / name).read_bytes()).hexdigest()}  {name}\n"
            for name in vmd_names
        )
    )
    portable = [*payload, "vmd/MANIFEST.sha256"]
    (bundle / "MANIFEST.sha256").write_text(
        "".join(
            f"{hashlib.sha256((bundle / name).read_bytes()).hexdigest()}  {name}\n"
            for name in portable
        )
    )
    for name in (
        "native_report.mdout",
        "source_nebpath.xyz.gz",
        "optimized.0001.ncrst",
        "frame_links/accepted.ncrst",
    ):
        path = bundle / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("raw")
    old = root / "progress/update_000001"
    old.mkdir()
    (old / "images.json").write_text("old")
    (root / "progress/latest").symlink_to(bundle.name)
    subprocess.run(
        [rsync, "-rt", *relay_options(), f"{source}/", f"{target}/"], check=True
    )
    copied = target / relative
    assert all((copied / name).is_file() for name in keep)
    assert not any((copied / name).exists() for name in reject)
    assert not list((copied / "progress").glob("update_*"))
    latest = copied / "progress/latest"
    assert latest.is_dir() and not latest.is_symlink()
    for manifest in (latest / "MANIFEST.sha256", latest / "vmd/MANIFEST.sha256"):
        for line in manifest.read_text().splitlines():
            expected, name = line.split("  ", 1)
            assert (
                hashlib.sha256((manifest.parent / name).read_bytes()).hexdigest()
                == expected
            )
    assert not (latest / "native_report.mdout").exists()
    assert not (latest / "source_nebpath.xyz.gz").exists()
    assert not (latest / "frame_links").exists()
