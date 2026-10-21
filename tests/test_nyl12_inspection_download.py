"""Check curated downloads, byte integrity, and laptop command safety."""

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from scripts import prepare_nyl12_inspection_download as export

ROOT = Path(__file__).resolve().parents[1]
PULL = ROOT / "scripts/pull_nyl12_inspection.sh"


def test_package_integrity_and_collision(tmp_path: Path) -> None:
    """Copied bytes and both manifests agree, with no silent overwrites."""
    source = tmp_path / "original.xyz"
    source.write_text("1\nexample\nH 0 0 0\n")
    output = tmp_path / "download"
    report = export.build_package({"qm/final.xyz": source}, output)
    assert report["payload_bytes"] == source.stat().st_size
    assert (output / "qm/final.xyz").read_bytes() == source.read_bytes()
    assert json.loads((output / "provenance.json").read_text()) == report
    for line in (output / "MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split("  ", 1)
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == digest
    with pytest.raises(FileExistsError):
        export.build_package({"qm/final.xyz": source}, output)


@pytest.mark.parametrize(
    "name", ["../escape", "/escape", "bad\nname", "MANIFEST.sha256", "provenance.json"]
)
def test_unsafe_paths(tmp_path: Path, name: str) -> None:
    source = tmp_path / "source"
    source.write_text("data")
    with pytest.raises(ValueError):
        export.build_package({name: source}, tmp_path / "download")
    assert not (tmp_path / "download").exists()


def test_limits_and_missing_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    output = tmp_path / "download"
    with pytest.raises(FileNotFoundError):
        export.build_package({"a": source}, output)
    source.write_text("12345")
    monkeypatch.setattr(export, "MAX_FILE_BYTES", 4)
    with pytest.raises(ValueError, match="small regular"):
        export.build_package({"a": source}, output)
    monkeypatch.setattr(export, "MAX_FILE_BYTES", 10)
    monkeypatch.setattr(export, "MAX_TOTAL_BYTES", 4)
    with pytest.raises(ValueError, match="25 MiB"):
        export.build_package({"a": source}, output)
    assert not output.exists()


def test_inventory_and_cli(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Only small inspection types enter the explicit inventory."""
    old = tmp_path / "output_nyl12_resume_20261020/nyl12_pymol_starting_models"
    old.mkdir(parents=True)
    for name in ("model.pdb", "view.pml", "scratch.chk", "raw.log"):
        (old / name).write_text("fixture")
    for folder in ("docs", "manuscript"):
        (tmp_path / folder).mkdir()
        (tmp_path / folder / "note.md").write_text("note")
    inventory = export.selected_files(
        tmp_path, tmp_path / "audit", tmp_path / "optimized", tmp_path / "trajectories"
    )
    assert "protein/model.pdb" in inventory
    assert "protein/scratch.chk" not in inventory
    assert "protein/raw.log" not in inventory
    assert (
        inventory["qm/neutral_thr_nme/conformer_01/optimized.sdf"]
        == tmp_path / "optimized/conformer_01.sdf"
    )
    source = tmp_path / "source"
    source.write_text("source")
    monkeypatch.setattr(export, "selected_files", lambda *args: {"small.txt": source})
    output = tmp_path / "export"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "export",
            "--repo",
            str(tmp_path),
            "--qm-audit",
            str(source),
            "--qm-structures",
            str(tmp_path),
            "--qm-trajectories",
            str(tmp_path),
            "--output",
            str(output),
        ],
    )
    export.main()
    assert (output / "small.txt").read_text() == "source"


def test_rsync_dry_run_and_arguments(tmp_path: Path) -> None:
    """A laptop destination with spaces stays one argument; no deletion flag."""
    binary = tmp_path / "bin"
    binary.mkdir()
    captured = tmp_path / "args"
    for name, body in {
        "rsync": '#!/bin/sh\nprintf "%s\\n" "$@" > "$CAPTURE_ARGS"\n',
        "ssh": "#!/bin/sh\nexit 0\n",
        "shasum": "#!/bin/sh\nexit 0\n",
    }.items():
        path = binary / name
        path.write_text(body)
        path.chmod(0o755)
    environment = {
        **os.environ,
        "PATH": f"{binary}:{os.environ['PATH']}",
        "CAPTURE_ARGS": str(captured),
    }
    destination = tmp_path / "laptop files"
    for dry_run in (True, False):
        argv = [
            "bash",
            str(PULL),
            *(["--dry-run"] if dry_run else []),
            "jpc@violet",
            str(destination),
        ]
        subprocess.run(argv, env=environment, check=True, capture_output=True)
        arguments = captured.read_text().splitlines()
        assert ("--dry-run" in arguments) == dry_run
        assert "--delete" not in arguments
        assert arguments[-1] == f"{destination}/"
        assert arguments[-2].endswith("output_laptop_inspection_20261021/")
        assert destination.exists() is not dry_run


@pytest.mark.parametrize(
    "args", [[], ["-evil"], ["host;touch foo"], ["goodhost", "-bad"]]
)
def test_rsync_rejects_invalid_arguments(args: list[str]) -> None:
    result = subprocess.run(["bash", str(PULL), *args], capture_output=True)
    assert result.returncode == 2
