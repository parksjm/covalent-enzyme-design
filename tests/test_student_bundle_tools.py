"""Synthetic stdlib bundle integrity/profile checks; remote I/O is always stubbed."""

from __future__ import annotations

import hashlib
import importlib
import io
import json
from pathlib import Path
import shutil
import tarfile
from typing import Any

import pytest

TOOLS = Path(__file__).resolve().parents[1] / "student_handoff/nyl12_ash274_md/tools"


@pytest.fixture
def modules(monkeypatch: pytest.MonkeyPatch) -> tuple[Any, Any, Any]:
    """Load isolated student command modules without executing their CLI."""
    monkeypatch.syspath_prepend(str(TOOLS))
    return tuple(
        importlib.import_module(name)
        for name in ("build_bundle", "fetch_bundle", "student_bundle_support")
    )


def entry(
    path: Path, destination: str, profile: str = "core", required: bool = True
) -> dict:
    """Declare an exact synthetic payload with source size and digest."""
    return {
        "source": str(path.absolute()),
        "destination": destination,
        "profile": profile,
        "required": required,
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "category": "test",
    }


def manifest(tmp_path: Path, entries: list[dict]) -> Path:
    """Write one explicit manifest without discovering surrounding files."""
    path = tmp_path / "data_manifest.json"
    path.write_text(json.dumps({"entries": entries}))
    return path


def payload(tmp_path: Path) -> Path:
    """Return a synthetic scientific text file, not production data."""
    path = tmp_path / "small.txt"
    path.write_text("mapped atoms and protocol\n")
    return path


def test_profiles_verified_archive_and_optional_missing(
    tmp_path: Path, modules: tuple
) -> None:
    builder, _, support = modules
    source = payload(tmp_path)
    missing = {
        **entry(source, "optional.txt"),
        "source": str(tmp_path / "absent"),
        "required": False,
    }
    raw = {
        **entry(source, "raw.nc", "raw_reference_only"),
        "sha256": None,
        "source": "/must/not/be/opened/75GB.nc",
    }
    path = manifest(
        tmp_path,
        [
            entry(source, "data/input.txt"),
            entry(source, "movie/view.txt", "inspection"),
            missing,
            raw,
        ],
    )
    for profile, name in (("core", "core_data"), ("inspection", "inspection")):
        output = tmp_path / name
        archive = builder.build(path, profile, output)
        support.verify_contents(output)
        assert (output / "data/input.txt").read_bytes() == source.read_bytes()
        assert (output / "movie/view.txt").exists() == (profile == "inspection")
        assert not (output / "raw.nc").exists()
        exported = json.loads((output / "export_manifest.json").read_text())
        assert exported["skipped_optional_sources"] == [str(tmp_path / "absent")]
        assert "data/input.txt" in (output / "README.md").read_text()
        assert (
            archive.with_name(archive.name + ".sha256").read_text()
            == f"{support.digest(archive)}  {archive.name}\n"
        )
        with tarfile.open(archive) as tar:
            assert all(m.name.startswith(name) for m in tar.getmembers())
        with pytest.raises(FileExistsError):
            builder.build(path, profile, output)


@pytest.mark.parametrize(
    "destination", ["../escape", "/absolute", "a/../escape", "a\\bad", "a//b", "a\nb"]
)
def test_unsafe_destinations_fail_before_output(
    tmp_path: Path, modules: tuple, destination: str
) -> None:
    builder, _, _ = modules
    path = manifest(tmp_path, [entry(payload(tmp_path), destination)])
    with pytest.raises(ValueError):
        builder.build(path, "core", tmp_path / "output")
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("destinations", [("README.md",), ("a", "a"), ("a", "a/b")])
def test_generated_duplicate_and_parent_conflicts(
    tmp_path: Path, modules: tuple, destinations: tuple
) -> None:
    builder, _, _ = modules
    source = payload(tmp_path)
    with pytest.raises(ValueError):
        builder.build(
            manifest(tmp_path, [entry(source, p) for p in destinations]),
            "core",
            tmp_path / "output",
        )


@pytest.mark.parametrize("field,value", [("sha256", "0" * 64), ("bytes", 999)])
def test_changed_payload_fails_closed(
    tmp_path: Path, modules: tuple, field: str, value: Any
) -> None:
    builder, _, _ = modules
    source = payload(tmp_path)
    with pytest.raises(ValueError, match="size/SHA"):
        builder.build(
            manifest(tmp_path, [{**entry(source, "input"), field: value}]),
            "core",
            tmp_path / "output",
        )
    assert not (tmp_path / "output.tar.gz").exists()


def test_missing_required_symlink_and_executable_rejected(
    tmp_path: Path, modules: tuple
) -> None:
    builder, _, _ = modules
    source = payload(tmp_path)
    original = entry(source, "input")
    with pytest.raises(FileNotFoundError):
        builder.selected_entries(
            {"entries": [{**original, "source": str(tmp_path / "missing")}]}, "core"
        )
    link = tmp_path / "link"
    link.symlink_to(source)
    with pytest.raises(ValueError, match="Symlink"):
        builder.selected_entries(
            {"entries": [{**original, "source": str(link)}]}, "core"
        )
    binary = tmp_path / "program.so"
    binary.write_bytes(b"not a scientific payload")
    with pytest.raises(ValueError, match="binaries"):
        builder.selected_entries({"entries": [entry(binary, "program")]}, "core")
    elf = tmp_path / "renamed.dat"
    elf.write_bytes(b"\x7fELFexecutable")
    with pytest.raises(ValueError, match="magic"):
        builder.build(
            manifest(tmp_path, [entry(elf, "program")]), "core", tmp_path / "output"
        )


def test_manifest_schema_empty_and_contents_integrity(
    tmp_path: Path, modules: tuple
) -> None:
    builder, _, support = modules
    source = payload(tmp_path)
    row = entry(source, "input")
    for change in (
        {"required": 1},
        {"bytes": -1},
        {"sha256": None},
        {"category": ""},
        {"profile": "unknown"},
        {"source": "relative"},
    ):
        with pytest.raises(ValueError):
            builder.selected_entries({"entries": [{**row, **change}]}, "core")
    with pytest.raises(ValueError, match="No available"):
        builder.selected_entries({"entries": []}, "core")
    output = tmp_path / "core_data"
    builder.build(manifest(tmp_path, [row]), "core", output)
    (output / "input").write_text("tampered")
    with pytest.raises(ValueError, match="checksum"):
        support.verify_contents(output)
    (output / "input").write_bytes(source.read_bytes())
    (output / "extra").write_text("unlisted")
    with pytest.raises(ValueError, match="unlisted"):
        support.verify_contents(output)


@pytest.mark.parametrize("extract", [False, True])
def test_fetch_downloads_only_selected_archive_and_own_checksum(
    tmp_path: Path, modules: tuple, monkeypatch: pytest.MonkeyPatch, extract: bool
) -> None:
    builder, fetcher, _ = modules
    published = tmp_path / "published"
    published.mkdir()
    source = payload(published)
    archive = builder.build(
        manifest(published, [entry(source, "data/input")]),
        "core",
        published / "core_data",
    )
    calls = []

    def rsync(command: list[str], check: bool) -> None:
        calls.append(command)
        assert check and command[:2] == ["rsync", "-rt"]
        assert all("inspection" not in arg for arg in command)
        assert command[-3:-1] == [
            "user@violet:/published/core_data.tar.gz",
            "user@violet:/published/core_data.tar.gz.sha256",
        ]
        target = Path(command[-1])
        shutil.copy2(archive, target / archive.name)
        shutil.copy2(
            archive.with_name(archive.name + ".sha256"),
            target / (archive.name + ".sha256"),
        )

    monkeypatch.setattr(fetcher.subprocess, "run", rsync)
    destination = tmp_path / "download"
    result = fetcher.fetch(
        "user@violet", "/published", "core_data", destination, extract
    )
    assert len(calls) == 1 and result.exists()
    if extract:
        assert (result / "data/input").read_bytes() == source.read_bytes()
    with pytest.raises(FileExistsError):
        fetcher.fetch("user@violet", "/published", "core_data", destination, extract)


@pytest.mark.parametrize(
    "record",
    [
        "bad checksum record\n",
        "0" * 64 + "  inspection.tar.gz\n",
        "0" * 64 + "  core_data.tar.gz\n",
    ],
)
def test_fetch_rejects_wrong_or_optional_checksum(
    tmp_path: Path, modules: tuple, monkeypatch: pytest.MonkeyPatch, record: str
) -> None:
    _, fetcher, _ = modules

    def rsync(command: list[str], check: bool) -> None:
        target = Path(command[-1])
        (target / "core_data.tar.gz").write_bytes(b"archive")
        (target / "core_data.tar.gz.sha256").write_text(record)

    monkeypatch.setattr(fetcher.subprocess, "run", rsync)
    with pytest.raises(ValueError):
        fetcher.fetch("violet", "/published", "core_data", tmp_path / "out", False)


@pytest.mark.parametrize(
    "name,kind",
    [
        ("core_data/../escape", "file"),
        ("/escape", "file"),
        ("inspection/file", "file"),
        ("core_data/link", "symlink"),
        ("core_data/link", "hardlink"),
        ("core_data/device", "device"),
    ],
)
def test_archive_traversal_links_and_devices_rejected(
    tmp_path: Path, modules: tuple, name: str, kind: str
) -> None:
    _, fetcher, _ = modules
    archive = tmp_path / "bad.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo(name)
        if kind == "symlink":
            info.type = tarfile.SYMTYPE
            info.linkname = "/outside"
        elif kind == "hardlink":
            info.type = tarfile.LNKTYPE
            info.linkname = "/outside"
        elif kind == "device":
            info.type = tarfile.CHRTYPE
        else:
            info.size = 1
        tar.addfile(info, io.BytesIO(b"x") if kind == "file" else None)
    with pytest.raises(ValueError):
        fetcher.extract_verified(archive, tmp_path, "core_data")
    assert not (tmp_path / "core_data").exists()


def test_fetch_validation_and_cli_dispatch(
    tmp_path: Path, modules: tuple, monkeypatch: pytest.MonkeyPatch
) -> None:
    builder, fetcher, _ = modules
    for source, remote, package in [
        ("-host", "/path", "core_data"),
        ("u@-host", "/path", "core_data"),
        ("host", "/a/../b", "core_data"),
        ("host", "relative", "core_data"),
        ("host", "/path", "both"),
    ]:
        with pytest.raises(ValueError):
            fetcher.fetch(source, remote, package, tmp_path / "out", False)
    calls = []
    monkeypatch.setattr(builder, "build", lambda *args: calls.append(args))
    builder.main(
        ["--manifest", "data.json", "--profile", "core", "--output", "core_data"]
    )
    assert calls[-1] == (Path("data.json"), "core", Path("core_data"))
    monkeypatch.setattr(fetcher, "fetch", lambda *args: calls.append(args))
    fetcher.main(
        [
            "--source",
            "violet",
            "--remote-dir",
            "/data",
            "--package",
            "inspection",
            "--destination",
            "out",
            "--extract",
        ]
    )
    assert calls[-1] == ("violet", "/data", "inspection", Path("out"), True)
