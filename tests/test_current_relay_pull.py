"""Exercise exact rsync rules for current relay movies and latest PMF."""

from pathlib import Path
import hashlib
import os
import shlex
import shutil
import subprocess

import pytest


def _rules(marker: str, variable: str) -> list[str]:
    text = Path("scripts/pull_nyl12_inspection.sh").read_text().split(marker, 1)[1]
    return shlex.split(
        text.split('rsync "${relay_current_rsync_args[@]}"', 1)[1]
        .split('"$remote:$' + variable + '/"', 1)[0]
        .replace("\\\n", "")
    )


def _current_transfer_args() -> list[str]:
    text = Path("scripts/pull_nyl12_inspection.sh").read_text()
    declaration = text.split("relay_current_rsync_args=(", 1)[1].split(")", 1)[0]
    return shlex.split(declaration)


def test_mutable_current_is_excluded_and_immutable_source_is_checksummed() -> None:
    """Keep the post-symlink mutable publication out of the legacy relay pull."""
    text = Path("scripts/pull_nyl12_inspection.sh").read_text()
    legacy = text.split("nterminal_relay_rsync_args=", 1)[1].split(
        "# Deacylation hydration:", 1
    )[0]
    assert (
        "--exclude='/umbrella_sampling/full_path_300K/vmd_progress/current/***'"
        in legacy
    )
    assert (
        "--include='/umbrella_sampling/full_path_300K/vmd_progress/current/***'"
        not in legacy
    )
    assert "current_us_movie_rsync_args" not in legacy
    assert 'movie_source=$(ssh "$remote" readlink -f "$movie_pointer")' in text
    assert "--checksum" in _current_transfer_args()


def test_production_current_transfer_repairs_same_size_same_mtime_corruption(
    tmp_path: Path,
) -> None:
    """Checksum mode repairs a laptop file rsync's size/mtime check would reuse."""
    rsync = shutil.which("rsync")
    if not rsync:
        pytest.skip("rsync unavailable")
    source = tmp_path / "immutable_verified_70ps_57windows_001"
    target = tmp_path / "laptop" / "vmd_progress" / "current"
    source.mkdir()
    target.mkdir(parents=True)
    source_file = source / "path.dcd"
    target_file = target / "path.dcd"
    source_file.write_bytes(b"verified trajectory payload")
    target_file.write_bytes(b"corruptd trajectory payload")
    timestamp = 1_700_000_000
    os.utime(source_file, (timestamp, timestamp))
    os.utime(target_file, (timestamp, timestamp))
    assert source_file.stat().st_size == target_file.stat().st_size
    assert source_file.stat().st_mtime == target_file.stat().st_mtime

    subprocess.run(
        [rsync, *_current_transfer_args(), f"{source}/", f"{target}/"],
        check=True,
    )
    assert target_file.read_bytes() == source_file.read_bytes()


@pytest.mark.parametrize("kind", ["movie", "pmf"])
def test_current_publication_filters(tmp_path: Path, kind: str) -> None:
    """Keep exactly manifest-listed payloads, including an oversized DCD name."""
    rsync = shutil.which("rsync")
    if not rsync:
        pytest.skip("rsync unavailable")
    if kind == "movie":
        rules = _rules(
            "# Current N-terminal relay publications:",
            "$movie_source".removeprefix("$"),
        )
        names = [
            "README.md",
            "coordinate_probe_first.in",
            "coordinate_probe_first.log",
            "coordinate_probe_first.pdb",
            "coordinate_probe_last.in",
            "coordinate_probe_last.log",
            "coordinate_probe_last.pdb",
            "manifest.json",
            "path.dcd",
            "path.parm7",
            "qm_indices.json",
            "view_path.vmd",
            "write_topology.in",
            "write_topology.log",
            "write_trajectory.in",
            "write_trajectory.log",
        ]
        manifest = "SHA256SUMS.txt"
    else:
        rules = _rules("# Latest eligible-data PMF:", "relay_pmf_source")
        names = [
            "latest_pmf_comparison.png",
            "latest.json",
            "point.json",
            "sources.json",
            "weights.npz",
        ]
        manifest = "SHA256SUMS.txt"
    source = tmp_path / "source"
    source.mkdir()
    for name in names:
        (source / name).write_text(name)
    (source / manifest).write_text(
        "".join(
            f"{hashlib.sha256((source / n).read_bytes()).hexdigest()}  {n}\n"
            for n in names
        )
    )
    for name in [
        "private/secret.json",
        "columns/receipt.json",
        "Hessian.json",
        "large.nc",
        "unlisted.bin",
    ]:
        f = source / name
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text("excluded")
    target = tmp_path / "laptop"
    subprocess.run([rsync, "-rt", *rules, f"{source}/", f"{target}/"], check=True)
    assert {
        str(f.relative_to(target)) for f in target.rglob("*") if f.is_file()
    } == set(names) | {manifest}
    for line in (target / manifest).read_text().splitlines():
        expected, name = line.split("  ", 1)
        assert hashlib.sha256((target / name).read_bytes()).hexdigest() == expected
    if kind == "movie":
        with (source / "path.dcd").open("wb") as f:
            f.truncate(1300 * 1024 * 1024)
        result = subprocess.run(
            [
                rsync,
                "-rt",
                "--dry-run",
                "--out-format=%n",
                *rules,
                f"{source}/",
                f"{tmp_path}/large/",
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        assert "path.dcd" in result.stdout
        text = (
            Path("scripts/pull_nyl12_inspection.sh")
            .read_text()
            .split("# Current N-terminal relay publications:", 1)[1]
        )
        assert "--max-size" not in text
        assert "additional_sampling_20261126/vmd_progress/current" in text


def test_removed_legacy_movie_variable_is_not_referenced() -> None:
    """Do not leave an unbound set-u variable after removing the legacy block."""
    script = Path("scripts/pull_nyl12_inspection.sh").read_text()
    assert "$current_us_movie" not in script
    assert 'echo "Current relay VMD: $movie_destination/view_path.vmd"' in script
