"""Tests for source identity, scope and unchanged native release dispatch."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from scripts import run_plateau_relay_release as release


@pytest.fixture
def prepared(tmp_path: Path) -> tuple[Path, dict]:
    """Make a small source-bound plan without running scientific calculations."""
    seed = tmp_path / "start.ncrst"
    seed.write_bytes(b"original coordinate fixture")
    source = tmp_path / "source.py"
    source.write_text("immutable source\n")
    plan = {
        "campaign": str(tmp_path),
        "source_sha256": {str(source): hashlib.sha256(source.read_bytes()).hexdigest()},
        "images": [
            {
                "image": 11,
                "reuse_existing_release": False,
                "native_directory": str(tmp_path / "image_011/native"),
                "start_restart": str(seed),
                "start_sha256": hashlib.sha256(seed.read_bytes()).hexdigest(),
            }
        ],
    }
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    return path, plan


def test_verified_record(prepared: tuple[Path, dict]) -> None:
    _, plan = prepared
    assert release.verified_record(plan, 11) == plan["images"][0]


@pytest.mark.parametrize("image", [10, 12, 22])
def test_reject_unprepared_image(prepared: tuple[Path, dict], image: int) -> None:
    with pytest.raises(ValueError, match="plateau image"):
        release.verified_record(prepared[1], image)


@pytest.mark.parametrize(
    "kind", ["source", "coordinates", "duplicate", "outside", "reuse"]
)
def test_reject_bad_identity_or_scope(prepared: tuple[Path, dict], kind: str) -> None:
    _, plan = prepared
    row = plan["images"][0]
    if kind == "source":
        Path(next(iter(plan["source_sha256"]))).write_text("changed")
    elif kind == "coordinates":
        Path(row["start_restart"]).write_text("changed")
    elif kind == "duplicate":
        plan["images"].append(row.copy())
    elif kind == "outside":
        row["native_directory"] = str(Path(plan["campaign"]).parent / "outside")
    else:
        row["reuse_existing_release"] = True
    with pytest.raises(ValueError):
        release.verified_record(plan, 11)


def test_preserve_existing_attempt(prepared: tuple[Path, dict]) -> None:
    Path(prepared[1]["images"][0]["native_directory"]).mkdir(parents=True)
    with pytest.raises(FileExistsError):
        release.verified_record(prepared[1], 11)


@pytest.mark.parametrize("accepted", [True, False])
def test_native_dispatch_without_bias(
    prepared: tuple[Path, dict], monkeypatch: pytest.MonkeyPatch, accepted: bool
) -> None:
    """Pass the full region through unchanged and require native acceptance."""
    calls = []

    def native(*args: object, **kwargs: object) -> dict:
        calls.append((args, kwargs))
        return {"accepted": accepted}

    monkeypatch.setitem(
        sys.modules,
        "nyl12_covalent.nterminal_relay_native",
        SimpleNamespace(run_native=native),
    )
    monkeypatch.setitem(
        sys.modules,
        "nyl12_covalent.relay_release",
        SimpleNamespace(minimization_input=lambda s, **kw: (s, kw)),
    )
    monkeypatch.setitem(
        sys.modules,
        "nyl12_covalent.stationary_native",
        SimpleNamespace(NativeState=lambda p: SimpleNamespace(mdin_source="frozen")),
    )
    if accepted:
        assert release.run_release(prepared[0], 11) == {"accepted": True}
    else:
        with pytest.raises(RuntimeError, match="acceptance"):
            release.run_release(prepared[0], 11)
    assert calls[0][1] == {"restraints": None}
    assert calls[0][0][3] == ("frozen", {"maxcycles": 10000})


def test_cli_arguments(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    calls = []
    monkeypatch.setattr(
        sys, "argv", ["release", "--plan", str(tmp_path / "plan"), "--image", "11"]
    )
    monkeypatch.setattr(release, "run_release", lambda p, i: calls.append((p, i)))
    release.main()
    assert calls == [(tmp_path / "plan", 11)]
