"""Exercise the actual native-Thr overlay rsync filters against fixture files."""

import hashlib
from pathlib import Path
import shlex
import shutil
import subprocess

import pytest


def test_overlay_bundle_is_portable_and_excludes_coordinate_cache(
    tmp_path: Path,
) -> None:
    """Copy complete viewer inputs and manifest, without large raw coordinates."""
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("Actual rsync is required to check filter behavior")
    text = Path("scripts/pull_nyl12_inspection.sh").read_text()
    block = text.split("# Native-Thr MD diagnostic:", 1)[1].split(
        "# Current controlling plans", 1
    )[0]
    rules = shlex.split(
        block.split('rsync "${vmd_rsync_args[@]}"', 1)[1]
        .split('"$remote:$thr_diag_source/"', 1)[0]
        .replace("\\\n", "")
    )
    source, target = tmp_path / "source", tmp_path / "target"
    keep = [
        "reference.pdb",
        "WT_typical.pdb",
        "broader_mask07_typical.pdb",
        "native_thr.pse",
        "view_native_thr.py",
        "native_thr_typical_side_by_side.png",
        "gate_distributions.pdf",
        "snapshot_manifest.json",
        "scalar_gate_report.csv",
        "README.md",
    ]
    reject = [
        "coordinate_evidence/WT_typical.npz",
        "coordinate_evidence/unneeded.json",
        "raw/trajectory.nc",
        "system.prmtop",
    ]
    for name in keep + reject:
        p = source / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(name)
    (source / "MANIFEST.sha256").write_text(
        "".join(
            f"{hashlib.sha256((source / n).read_bytes()).hexdigest()}  {n}\n"
            for n in keep
        )
    )
    subprocess.run([rsync, "-rt", *rules, f"{source}/", f"{target}/"], check=True)
    assert all((target / n).is_file() for n in keep)
    assert all(not (target / n).exists() for n in reject)
    for line in (target / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        assert hashlib.sha256((target / name).read_bytes()).hexdigest() == expected
    assert '"${checksum[@]}" MANIFEST.sha256' in block
