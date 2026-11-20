"""Scientific payload validation and acceptance-last publication checks."""

from pathlib import Path
import numpy as np
import pytest
from nyl12_covalent import neb_umbrella_sync as s
from neb_umbrella_sync_helpers import STAGE, WINDOW, hashes, view
from neb_umbrella_sync_helpers import campaign as campaign


def test_acceptance_last_and_immutable(
    campaign: tuple, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, cfg, incoming, values, expected = campaign
    target = path.parent.parent / WINDOW / STAGE
    real_copy = s.atomic_copy
    order = []

    def tracked(source: Path, destination: Path) -> None:
        if destination.name != "acceptance.json":
            assert not (target / "acceptance.json").exists()
        else:
            assert (target / "transfer_receipt.json").is_file()
            assert all((target / p).is_file() for p in s.CRITICAL[:-1])
        order.append(destination.name)
        real_copy(source, destination)

    monkeypatch.setattr(s, "atomic_copy", tracked)
    assert s.publish_block(incoming, target, values, cfg, expected)
    assert order[-1] == "acceptance.json"
    assert not s.publish_block(incoming, target, values, cfg, expected)
    (target / "process.json").write_text("changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        s.publish_block(incoming, target, values, cfg, expected)


def test_interrupted_publication_resumes(
    campaign: tuple, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, cfg, incoming, values, expected = campaign
    target = path.parent.parent / WINDOW / STAGE
    real_copy = s.atomic_copy
    calls = []

    def interrupted(source: Path, destination: Path) -> None:
        calls.append(source)
        if len(calls) == 3:
            raise OSError("simulated disk error")
        real_copy(source, destination)

    monkeypatch.setattr(s, "atomic_copy", interrupted)
    with pytest.raises(OSError, match="disk error"):
        s.publish_block(incoming, target, values, cfg, expected)
    assert not (target / "acceptance.json").exists()
    monkeypatch.setattr(s, "atomic_copy", real_copy)
    assert s.publish_block(incoming, target, values, cfg, expected)


@pytest.mark.parametrize(
    "case", ["hash", "frames", "csv", "nonfinite", "npz", "bias", "process"]
)
def test_invalid_payloads(campaign: tuple, case: str) -> None:
    _, cfg, incoming, values, _ = campaign
    if case == "hash":
        (incoming / "process.json").write_text("changed")
    elif case in ("frames", "csv"):
        receipt = s.read(incoming / "observables/observables.json")
        receipt["frames" if case == "frames" else "csv_sha256"] = 999
        s.write(incoming / "observables/observables.json", receipt)
    elif case == "nonfinite":
        csv = incoming / "observables/observables.csv"
        csv.write_text(csv.read_text().replace("1,2,3,4", "nan,2,3,4"))
        receipt = s.read(incoming / "observables/observables.json")
        receipt["csv_sha256"] = s.digest(csv)
        s.write(incoming / "observables/observables.json", receipt)
    elif case in ("npz", "bias"):
        np.savez_compressed(
            incoming / "samples_and_cross_bias.npz",
            cv_values=[[9.0, 2.0], [2.0, 3.0]]
            if case == "npz"
            else [[1.0, 2.0], [2.0, 3.0]],
            bias_kcal_mol=np.full((63, 2), np.nan),
            u_kn_bias_only=np.zeros((63, 2)),
        )
    else:
        s.write(incoming / "process.json", {"exit_status": 1, "ended_utc": "now"})
    if case != "hash":
        values = hashes(incoming)
    with pytest.raises(ValueError):
        s.validate_block(incoming, values, cfg)


def test_partial_payload_and_receipt_guards(campaign: tuple) -> None:
    path, cfg, incoming, values, expected = campaign
    target = path.parent.parent / WINDOW / STAGE
    target.mkdir(parents=True)
    (target / "process.json").write_text("unrelated previous partial content")
    with pytest.raises(ValueError, match="partial block"):
        s.publish_block(incoming, target, values, cfg, expected)
    with pytest.raises(ValueError, match="critical payloads"):
        s.validate_block(incoming, {}, cfg)
    with pytest.raises(ValueError, match="Unsafe"):
        s.verify_hashes(incoming, {"../escape": "a"})
    linked = incoming / "linked"
    linked.symlink_to(incoming / "process.json")
    with pytest.raises(ValueError, match="hash mismatch"):
        s.verify_hashes(incoming, {"linked": s.digest(linked)})
    target2 = path.parent.parent / WINDOW / "other"
    s.publish_block(incoming, target2, values, cfg, expected)
    changed = {**expected, "source_host": "other"}
    with pytest.raises(ValueError, match="transfer identity"):
        s.published(target2, values, changed)


def test_view_publication(campaign: tuple, tmp_path: Path) -> None:
    path, _, _, _, _ = campaign
    incoming = tmp_path / "view"
    values = view(incoming)
    target = path.parent.parent / WINDOW / "progress/block_01"
    expected = s.identity(path, "azalea", "remote/view", values)
    assert s.publish_view(incoming, target, values, expected)
    assert not s.publish_view(incoming, target, values, expected)
    assert (target / "transfer_receipt.json").exists()
    other = target.with_name("block_02")
    other.mkdir()
    with pytest.raises(FileExistsError):
        s.publish_view(incoming, other, values, expected)
    (incoming / "MANIFEST.sha256").write_text("incomplete")
    with pytest.raises(ValueError, match="Incomplete"):
        s.manifest_hashes(incoming)


def test_context_identity_and_ownership(campaign: tuple) -> None:
    path, _, _, _, _ = campaign
    _, _, prepared, blocks = s.context(path, "azalea")
    assert len(blocks) == 47 and "window_001" not in blocks
    with pytest.raises(ValueError, match="explicit remote"):
        s.context(path, "violet")
    prepared.write_text(prepared.read_text() + " ")
    with pytest.raises(ValueError, match="Master prepared"):
        s.context(path, "azalea")
    path.write_text(path.read_text() + " ")
    with pytest.raises(ValueError, match="plan changed"):
        s.context(path, "azalea")


def test_latest_rejects_incomplete_and_preserves_unrelated(tmp_path: Path) -> None:
    generation = tmp_path / "block_01"
    generation.mkdir()
    with pytest.raises(ValueError, match="incomplete"):
        s.publish_latest(generation)
    (generation / "MANIFEST.sha256").write_text("complete")
    s.publish_latest(generation)
    assert (tmp_path / "latest").resolve() == generation
    s.publish_latest(generation)
    other = tmp_path / "other"
    other.mkdir()
    (other / "MANIFEST.sha256").write_text("complete")
    (other.parent / "latest").replace(other.parent / "old_pointer")
    (other.parent / "latest").write_text("unrelated")
    with pytest.raises(ValueError, match="non-symlink"):
        s.publish_latest(other)
