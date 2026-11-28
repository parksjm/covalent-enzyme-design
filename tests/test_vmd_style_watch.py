"""Check bounded display watcher stability and strict topology identity pins."""

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from nyl12_covalent import vmd_style_watch as watch
from nyl12_covalent.vmd_style_manifests import digest


def row(tmp_path: Path) -> dict:
    """Create an explicitly pinned pair of mock portable movie files."""
    viewer, topology = tmp_path / "view.vmd", tmp_path / "stripped.parm7"
    viewer.write_text("original")
    topology.write_text("validated topology")
    return dict(
        viewer=str(viewer),
        topology=str(topology),
        topology_sha256=digest(topology),
        topology_identity_sha256="physical_identity",
        qm_indices=[1],
        qm_evidence="source",
    )


def test_stability_checksum_readiness_styled_and_cheap_unchanged(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    item, state = row(tmp_path), {}
    assert watch.inspect_viewer(item, state) == "unstable"
    assert watch.inspect_viewer(item, state) == "refresh"
    monkeypatch.setattr(watch, "verify_checksums", lambda _: ["not ready"])
    assert watch.inspect_viewer(item, state) == "waiting_checksums"
    Path(item["viewer"]).write_text(watch.BEGIN + "\nold policy\n" + watch.END)
    assert watch.inspect_viewer(item, state) == "unstable"
    assert watch.inspect_viewer(item, state) == "waiting_checksums"
    Path(item["viewer"]).write_text(
        watch.BEGIN
        + f"\n# Display policy version: {watch.POLICY_VERSION}\n"
        + watch.END
    )
    assert watch.inspect_viewer(item, state) == "unstable"
    assert watch.inspect_viewer(item, state) == "styled"
    monkeypatch.setattr(
        watch, "digest", lambda _: pytest.fail("unchanged parsed or hashed")
    )
    assert watch.inspect_viewer(item, state) == "unchanged"
    Path(item["viewer"]).unlink()
    assert watch.inspect_viewer(item, state) == "waiting_files"


def test_equivalent_export_header_is_allowed_but_changed_physics_holds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    item, state = row(tmp_path), {}
    Path(item["topology"]).write_text("new header")
    monkeypatch.setattr(watch, "AmberParm", lambda _: object())
    monkeypatch.setattr(watch, "topology_identity", lambda _: "physical_identity")
    assert watch.inspect_viewer(item, state) == "unstable"
    assert watch.inspect_viewer(item, state) == "refresh"
    assert state["equivalent_topology"]["pinned_raw_sha256"] == item["topology_sha256"]
    monkeypatch.setattr(watch, "topology_identity", lambda _: "changed_physics")
    with pytest.raises(ValueError, match="HOLD"):
        watch.inspect_viewer(item, state)


def test_identity_hash_changes_for_order_graph_and_charge() -> None:
    residue = SimpleNamespace(idx=0, name="ALA", chain="A")
    atoms = [
        SimpleNamespace(
            idx=i,
            name=f"C{i}",
            atomic_number=6,
            residue=residue,
            type="CT",
            charge=0,
            mass=12,
        )
        for i in range(2)
    ]
    parm = SimpleNamespace(atoms=atoms, bonds=[])
    initial = watch.topology_identity(parm)
    atoms[0].charge = 0.1
    assert watch.topology_identity(parm) != initial
    atoms[0].charge = 0
    parm.atoms = list(reversed(atoms))
    assert watch.topology_identity(parm) != initial
    parm.atoms = atoms
    parm.bonds = [SimpleNamespace(atom1=atoms[0], atom2=atoms[1])]
    assert watch.topology_identity(parm) != initial


def test_script_watch_refresh_errors_and_bounds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    spec = importlib.util.spec_from_file_location(
        "watch_cli", Path(__file__).resolve().parents[1] / "scripts/watch_vmd_style.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps(dict(movies=[row(tmp_path)])))
    with pytest.raises(ValueError):
        module.run_watch(plan, tmp_path / "output", 61, 72)
    for mode in ("refresh", "error"):
        clock = iter([0, 0, 0, 10000])
        monkeypatch.setattr(module.time, "monotonic", lambda: next(clock))
        monkeypatch.setattr(module.time, "sleep", lambda _: None)
        calls = []
        monkeypatch.setattr(
            module.subprocess, "run", lambda *a, **k: calls.append(a[0])
        )
        if mode == "error":

            def fail(*args: object) -> str:
                raise ValueError("changed identity")

            monkeypatch.setattr(module, "inspect_viewer", fail)
        else:
            monkeypatch.setattr(module, "inspect_viewer", lambda *a: "refresh")
        module.run_watch(plan, tmp_path / mode, 1, 1)
        assert calls
        assert (tmp_path / mode / "watch_status.json").exists()
