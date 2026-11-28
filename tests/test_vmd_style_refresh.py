"""Verify physical mapping, reversible patches and portable checksum updates."""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import parmed
import pytest

from nyl12_covalent import vmd_style_refresh as style
from nyl12_covalent.vmd_style_manifests import (
    backup_file,
    checksum_updates,
    digest,
    verify_checksums,
)


def topology() -> parmed.Structure:
    """Provide native chain identity and a substrate QM/MM carbon edge."""
    parm = parmed.Structure()
    for name, residue, number in (
        ("CA", "ALA", 1),
        ("CB", "ALA", 1),
        ("C1", "LIG", 2),
        ("C2", "LIG", 2),
    ):
        parm.add_atom(
            parmed.Atom(name=name, atomic_number=6), residue, number, chain="A"
        )
    for first, second in ((0, 1), (2, 3)):
        parm.bonds.append(parmed.Bond(parm.atoms[first], parm.atoms[second]))
    return parm


def test_policy_uses_explicit_substrate_native_edges_and_chain_groups() -> None:
    text = style.topology_policy(topology(), [1, 2])
    assert (
        "::nyl12_movie_policy::apply [molinfo top] {1 2} {2 3} {{0 1} {2 3}} {A {0 1}}"
        in text
    )
    compact = topology()[[False, True, True, True]]
    mapped = style.topology_policy(
        compact, [0, 1], original_topology=topology(), original_ids_1based=[2, 3, 4]
    )
    assert "{0 1} {1 2} {{1 2}} {A {0}}" in mapped


@pytest.mark.parametrize("qm", [[-1], [4], [True]])
def test_invalid_physical_qm_atom_fails(qm: list[int]) -> None:
    with pytest.raises(ValueError, match="physical"):
        style.topology_policy(topology(), qm)


@pytest.mark.parametrize("ids", [[1], [1, 1, 3, 4], [0, 2, 3, 4], [2, 1, 3, 4]])
def test_invalid_source_identity_mapping_fails(ids: list[int]) -> None:
    with pytest.raises(ValueError):
        style.topology_policy(
            topology(), [1], original_topology=topology(), original_ids_1based=ids
        )


def test_mapping_requires_source_and_ids_together() -> None:
    with pytest.raises(ValueError, match="together"):
        style.topology_policy(topology(), [1], original_topology=topology())


def test_correction_preserves_commands_and_is_replaceable() -> None:
    original = "load movie\nlabel existing\ncamera existing\n"
    first = style.corrected_viewer(original, "policy1")
    second = style.corrected_viewer(first, "policy2")
    assert second.startswith(original)
    assert "policy1" not in second
    assert second.count(style.BEGIN) == 1
    assert style.corrected_viewer(second, "policy2") == second
    for malformed in (
        style.BEGIN,
        style.END,
        style.END + style.BEGIN,
        style.BEGIN + style.BEGIN + style.END + style.END,
    ):
        with pytest.raises(ValueError):
            style.corrected_viewer(malformed, "policy")


def test_prepare_and_append_validate_before_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    viewer = tmp_path / "view.vmd"
    viewer.write_text("original\n")
    monkeypatch.setattr(style, "AmberParm", lambda _: topology())
    row = dict(
        viewer=str(viewer),
        topology=str(tmp_path / "parm7"),
        qm_indices=[1, 2],
        qm_evidence="explicit source qmmask",
    )
    path, text = style.prepare_correction(row)
    assert path == viewer and style.BEGIN in text
    assert viewer.read_text() == "original\n"
    style.append_topology_policy(viewer, topology(), [1, 2])
    assert viewer.read_text() == text
    row.update(source_topology="source.parm7", original_ids_1based=[1, 2, 3, 4])
    assert style.prepare_correction(row)[1] == text
    row["qm_evidence"] = ""
    with pytest.raises(ValueError, match="evidence"):
        style.prepare_correction(row)
    row["viewer"] = str(tmp_path / "source_snapshot/view.vmd")
    with pytest.raises(ValueError, match="Immutable"):
        style.prepare_correction(row)


def test_backups_and_derived_manifests_preserve_ancestry(tmp_path: Path) -> None:
    bundle, backups = tmp_path / "bundle", tmp_path / "clutter"
    bundle.mkdir()
    viewer = bundle / "view.vmd"
    viewer.write_text("old")
    old = digest(viewer)
    record = backup_file(viewer, tmp_path, backups)
    assert Path(record["backup"]).read_text() == "old"
    assert backup_file(viewer, tmp_path, backups) == record
    source = dict(
        source_sha256={"view.vmd": old},
        sources={"viewer": dict(sha256=old)},
        files={"view.vmd": old},
    )
    (bundle / "manifest.json").write_text(json.dumps(source))
    (bundle / "export_manifest.json").write_text(
        json.dumps(dict(bundle_files={"view.vmd": dict(sha256=old)}))
    )
    (bundle / "MANIFEST.sha256").write_text(f"{old}  view.vmd\n")
    updates = checksum_updates(bundle, {"view.vmd": "new"})
    for path, text in updates.items():
        path.write_text(text)
    viewer.write_text("new")
    assert not verify_checksums(bundle)
    changed = json.loads((bundle / "manifest.json").read_text())
    assert changed["sources"] == source["sources"]
    assert changed["source_sha256"] == source["source_sha256"]
    assert changed["files"]["view.vmd"] == digest(viewer)
    with pytest.raises(ValueError, match="Conflicting"):
        backup_file(viewer, tmp_path, backups)
    viewer.write_text("mismatch")
    assert len(verify_checksums(bundle)) == 1
    (bundle / "MANIFEST.sha256").write_text("invalid")
    with pytest.raises(ValueError, match="Unsupported"):
        checksum_updates(bundle, {})
    (bundle / "MANIFEST.sha256").write_text("")
    (bundle / "manifest.json").write_text(json.dumps(dict(files={"view.vmd": []})))
    with pytest.raises(ValueError, match="Unsupported"):
        checksum_updates(bundle, {"view.vmd": "new"})


def test_refresh_cli_applies_and_reports_only_derived_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    spec = importlib.util.spec_from_file_location(
        "refresh_cli",
        Path(__file__).resolve().parents[1] / "scripts/refresh_vmd_style.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    viewer = tmp_path / "view.vmd"
    viewer.write_text("original")
    plan, report = tmp_path / "plan.json", tmp_path / "report.json"
    plan.write_text(
        json.dumps(
            dict(
                movies=[
                    dict(viewer=str(viewer), qm_indices=[1], qm_evidence="evidence")
                ]
            )
        )
    )
    monkeypatch.setattr(module, "ROOT", tmp_path)
    monkeypatch.setattr(module, "prepare_correction", lambda row: (viewer, "corrected"))
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    args = [
        "refresh",
        "--plan",
        str(plan),
        "--backup",
        str(tmp_path / "clutter"),
        "--report",
        str(report),
    ]
    monkeypatch.setattr(sys, "argv", args)
    module.main()
    assert viewer.read_text() == "original"
    monkeypatch.setattr(sys, "argv", [*args, "--apply"])
    module.main()
    data = json.loads(report.read_text())
    assert (
        data["movies"][0]["changes"][0]["after_sha256"]
        == hashlib.sha256(b"corrected").hexdigest()
    )
    assert (tmp_path / "clutter/view.vmd").read_text() == "original"
    module.main()  # Exact repeat is a no-op, preserving the original backup.
    monkeypatch.setattr(
        module,
        "prepare_correction",
        lambda row: (_ for _ in ()).throw(ValueError("unmapped")),
    )
    with pytest.raises(ValueError, match="Mapping validation"):
        module.main()
    assert json.loads(report.read_text())["failures"]


def test_correction_rechecks_watcher_topology_pin_in_child(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    viewer, parm = tmp_path / "view.vmd", tmp_path / "topology.parm7"
    viewer.write_text("original")
    parm.write_text("raw topology")
    monkeypatch.setattr(style, "AmberParm", lambda _: topology())
    row = dict(
        viewer=str(viewer),
        topology=str(parm),
        qm_indices=[1],
        qm_evidence="source",
        topology_sha256=digest(parm),
    )
    assert style.prepare_correction(row)[0] == viewer
    from nyl12_covalent.vmd_style_watch import topology_identity

    row["topology_identity_sha256"] = topology_identity(topology())
    parm.write_text("new header")
    assert style.prepare_correction(row)[0] == viewer
    row["topology_identity_sha256"] = "changed identity"
    with pytest.raises(ValueError, match="Pinned"):
        style.prepare_correction(row)
    digests = iter(["before", "after"])
    monkeypatch.setattr(style, "digest", lambda _: next(digests))
    with pytest.raises(ValueError, match="while preparing"):
        style.prepare_correction(row)


def test_bundle_caption_replacement_is_explicit_and_repeatable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    viewer = tmp_path / "view.vmd"
    old = 'puts "Frame0: exact sampled seed. Later frames: optimizer iterations, not MD time or reaction path."'
    caption = 'puts "Combined accepted path endpoints; not continuous MD."'
    viewer.write_text(old)
    monkeypatch.setattr(style, "AmberParm", lambda _: topology())
    row = dict(
        viewer=str(viewer),
        topology="topology",
        qm_indices=[1],
        qm_evidence="mapping",
        display_caption=caption,
    )
    _, updated = style.prepare_correction(row)
    assert old not in updated and caption in updated
    viewer.write_text(updated)
    assert style.prepare_correction(row)[1] == updated
    viewer.write_text("unexpected caption")
    with pytest.raises(ValueError, match="caption"):
        style.prepare_correction(row)
