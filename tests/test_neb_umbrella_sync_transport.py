"""Bounded SSH/rsync contracts and immutable completed-bundle synchronization."""

import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace
import pytest
from nyl12_covalent import neb_umbrella_sync as s
from nyl12_covalent import neb_umbrella_sync_transport as t
from neb_umbrella_sync_helpers import STAGE, WINDOW, block, view
from neb_umbrella_sync_helpers import campaign as campaign


def test_sync_once_and_foreign_inventory(
    campaign: tuple, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    path, _, incoming, values, _ = campaign
    views = tmp_path / "remote_view"
    view_values = view(views)
    data = {
        "blocks": [{"window": WINDOW, "stage": STAGE, "sha256": values}],
        "views": [{"window": WINDOW, "generation": "block_01", "sha256": view_values}],
        "statuses": {f"{WINDOW}/window_status.json": {"status": "running"}},
    }
    monkeypatch.setattr(t, "inventory", lambda *a: data)
    fetched = []

    def fake_fetch(
        remote: str, root: Path, destination: Path, names: list[str]
    ) -> None:
        fetched.extend(names)
        for name in names:
            relative = Path(name)
            payload = (
                incoming / Path(*relative.parts[2:])
                if relative.parts[1] == STAGE
                else views / relative.name
            )
            target = destination / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(payload, target)

    monkeypatch.setattr(t, "fetch", fake_fetch)
    assert t.sync_once(path, "azalea") == {"blocks": 1, "views": 1, "statuses": 1}
    assert not any(
        name.endswith((".nc", ".ncrst", "prepared.json")) for name in fetched
    )
    fetched.clear()
    assert t.sync_once(path, "azalea")["blocks"] == 0
    assert not fetched
    assert not (path.parent.parent / WINDOW / "window_status.json").exists()
    assert (
        path.parent / "remote_status/azalea" / WINDOW / "window_status.json"
    ).exists()
    data["blocks"][0]["window"] = "window_001"
    with pytest.raises(ValueError, match="ownership"):
        t.sync_once(path, "azalea")


@pytest.mark.parametrize("case", ["status", "view", "duplicate"])
def test_inventory_rejection(campaign: tuple, case: str) -> None:
    path, _, _, values, _ = campaign
    _, _, _, blocks = s.context(path, "azalea")
    row = {"window": WINDOW, "stage": STAGE, "sha256": values}
    data = {"blocks": [row], "views": [], "statuses": {}}
    if case == "status":
        data["statuses"]["window_001/window_status.json"] = {}
    elif case == "view":
        data["views"].append(
            {"window": WINDOW, "generation": "../../escape", "sha256": {}}
        )
    else:
        data["blocks"].append(copy.deepcopy(row))
    with pytest.raises(ValueError, match="ownership"):
        t.validate_inventory(data, blocks, "azalea")


def test_frozen_source(campaign: tuple) -> None:
    path, _, _, _, _ = campaign
    script = t.freeze_sources(path, "azalea")
    t.verify_sources(script)
    with pytest.raises(FileExistsError):
        t.freeze_sources(path, "azalea")
    script.write_text("changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        t.verify_sources(script)


def test_ssh_rsync_commands(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    calls = []

    def run(command: list, **kwargs: object) -> SimpleNamespace:
        calls.append((command, kwargs))
        return SimpleNamespace(stdout='{"blocks": [], "views": [], "statuses": {}}')

    monkeypatch.setattr(t.subprocess, "run", run)
    assert t.inventory("azalea", tmp_path, {})["blocks"] == []
    t.fetch(
        "azalea", tmp_path, tmp_path / "incoming", ["window_017/stage/acceptance.json"]
    )
    assert calls[0][0][0] == "ssh"
    assert calls[1][0][0] == "rsync" and "--delete" not in calls[1][0]
    assert b"window_017/stage/acceptance.json\0" == calls[1][1]["input"]
    t.fetch("azalea", tmp_path, tmp_path / "incoming", [])
    assert len(calls) == 2


def test_remote_inventory_local_fixture(campaign: tuple) -> None:
    path, cfg, _, _, _ = campaign
    root = path.parent.parent
    folder = root / WINDOW / STAGE
    block(folder, cfg)
    (root / WINDOW / "progress/block_01").mkdir(parents=True)
    request = {
        "prepared": str(root / "prepared.json"),
        "prepared_sha256": s.digest(root / "prepared.json"),
        "blocks": {WINDOW: [STAGE]},
        "critical": s.CRITICAL,
        "remote": "azalea",
        "view_files": sorted(t.VIEW_FILES),
    }
    result = subprocess.run(
        [sys.executable, "-c", t.REMOTE_INVENTORY, str(root), json.dumps(request)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert len(json.loads(result.stdout)["blocks"]) == 1
    (folder / "analysis_request.json").rename(folder / "analysis_request.pending")
    result = subprocess.run(
        [sys.executable, "-c", t.REMOTE_INVENTORY, str(root), json.dumps(request)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert not json.loads(result.stdout)["blocks"]


def test_complete_remote_view_and_status_inventory(campaign: tuple) -> None:
    path, _, _, _, _ = campaign
    root = path.parent.parent
    generation = root / WINDOW / "progress/block_01"
    view(generation)
    s.write(root / WINDOW / "window_status.json", {"status": "running"})
    worker = root / "dispatch/azalea/worker01"
    worker.mkdir(parents=True)
    s.write(worker / "status.json", {"status": "running"})
    request = {
        "prepared": str(root / "prepared.json"),
        "prepared_sha256": s.digest(root / "prepared.json"),
        "blocks": {WINDOW: [STAGE]},
        "critical": s.CRITICAL,
        "remote": "azalea",
        "view_files": sorted(t.VIEW_FILES),
    }
    command = [sys.executable, "-c", t.REMOTE_INVENTORY, str(root), json.dumps(request)]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    data = json.loads(result.stdout)
    assert len(data["views"]) == 1 and len(data["statuses"]) == 2
    (generation / "MANIFEST.sha256").write_text("partial")
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    assert not json.loads(result.stdout)["views"]
    (generation / "MANIFEST.sha256").write_text("0" * 64 + "  trajectory.nc\n")
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(command, check=True, capture_output=True, text=True)


def test_subprocess_failures_propagate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fail(command: list, **kwargs: object) -> None:
        assert kwargs["check"] is True
        raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(t.subprocess, "run", fail)
    with pytest.raises(subprocess.CalledProcessError):
        t.inventory("azalea", tmp_path, {})
    with pytest.raises(subprocess.CalledProcessError):
        t.fetch("azalea", tmp_path, tmp_path / "incoming", ["acceptance.json"])
