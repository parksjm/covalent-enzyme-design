"""Check canonical native controls against real solved conversion artifacts."""

from __future__ import annotations

import json
from pathlib import Path
import shlex
import shutil
import subprocess

import pytest

from nyl12_covalent.control_inspection import (
    export_control,
    pdb_atoms,
    publish_subset,
    restore_modular,
    sha256,
)
from nyl12_covalent.control_pymol import write_views

BASE = Path(
    "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/dgx_wt_control_qualification_20261127"
)
CONFIG = Path("design/nyl12_ligandmpnn_pilot/boltz_validation/configuration")


@pytest.mark.parametrize(
    "predictor,seed",
    [
        ("af3", "20261126"),
        ("af3", "20261127"),
        ("boltz", "20261126"),
        ("boltz", "20261127"),
    ],
)
def test_real_native_roundtrip(tmp_path: Path, predictor: str, seed: str) -> None:
    """All four actual native CIFs survive solved modular heavy roundtrip."""
    candidates = list(
        (BASE / "comparison_20261127").glob(
            f"{predictor}_*{seed}*_contact_pairs/canonicalization_receipt.json"
        )
    )
    if not candidates:
        pytest.skip("Native prediction artifacts unavailable")
    source = Path(json.loads(candidates[0].read_text())["source"])
    maps = json.loads(
        (BASE.parent / "prediction_preparation/ligand_graph_manifest.json").read_text()
    )["ligand_atom_maps"][predictor]
    config = json.loads((CONFIG / "relaxation.json").read_text())
    result = export_control(
        source,
        maps,
        Path(config["mol2"]),
        Path(config["mapping"]),
        Path(config["release"]),
        tmp_path / "model",
    )
    assert result["raw_atom_count"] == result["roundtrip_atom_count"]
    assert result["roundtrip_max_coordinate_component_error_A"] <= 0.000501
    assert len(result["nearest_d1_per_site"]) == 4
    assert len(result["all48"]) == 48
    assert all(r["tier"] == "reject" for r in result["all48"])
    assert sha256(source) == sha256(tmp_path / "model/raw_prediction.cif")
    mapping = result["atom_mapping"]
    thr = next(r for r in mapping if r["mature_local"] == ["E", 1, "OG1"])
    assert thr["biological_precursor"] == ["E", 233, "OG1"]
    assert thr["amber_sequential"] == [929, "OG1"]
    assert len({tuple(r["amber_sequential"]) for r in mapping}) == len(mapping)
    conversion = json.loads(
        (tmp_path / "model/amber_sequential_input/conversion.json").read_text()
    )
    with pytest.raises(FileExistsError):
        restore_modular(
            tmp_path / "model/amber_sequential_input/complex_heavy_modular.pdb",
            conversion,
            tmp_path / "model/roundtrip_mature_local.pdb",
        )
    with pytest.raises(ValueError, match="132unique"):
        restore_modular(
            tmp_path / "model/amber_sequential_input/complex_heavy_modular.pdb",
            {"atom_mapping": []},
            tmp_path / "bad.pdb",
        )


def test_portable_publication_and_actual_pull_filters(tmp_path: Path) -> None:
    """Published manifests exactly match the subset actually copied by rsync."""
    rsync = shutil.which("rsync")
    if not rsync:
        pytest.skip("rsync unavailable")
    stage = tmp_path / "stage"
    stage.mkdir()
    write_views(stage, ["af3_seed_20261126"], [20, 21, 22])
    (stage / "af3_seed_20261126").mkdir()
    (stage / "af3_seed_20261126/canonical_precursor.pdb").write_text("END\n")
    (stage / "crystal_reference.pdb").write_text("END\n")
    (stage / "README.md").write_text("portable\n")
    published = tmp_path / "published"
    publish_subset(stage, published)
    text = Path("scripts/pull_nyl12_inspection.sh").read_text()
    block = text.split("# Native WT predictor controls:", 1)[1]
    rules = shlex.split(
        block.split('rsync "${vmd_rsync_args[@]}"', 1)[1]
        .split('"$remote:$control_source/"', 1)[0]
        .replace("\\\n", "")
    )
    for name in [
        "columns/receipt.json",
        "private/secret.json",
        "nested/hessian.npy",
        "nested/Hessian.json",
        "scratch.dat",
    ]:
        f = published / name
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text("excluded")
    target = tmp_path / "laptop"
    subprocess.run([rsync, "-rt", *rules, f"{published}/", f"{target}/"], check=True)
    for line in (target / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        assert sha256(target / name) == expected
    assert not (target / "columns").exists() or not list(
        (target / "columns").rglob("*.json")
    )
    assert not (target / "private/secret.json").exists()
    assert not (target / "nested/Hessian.json").exists()
    assert not (target / "scratch.dat").exists()
    assert "run view_controls.py" in (target / "load_all_controls.pml").read_text()
    assert "crystal_reference.pdb" in (target / "view_controls.py").read_text()
    assert str(tmp_path) not in (target / "view_controls.py").read_text()
    with pytest.raises(FileExistsError):
        publish_subset(published, target)
    symlinkstage = tmp_path / "symlinks"
    symlinkstage.mkdir()
    (symlinkstage / "link").symlink_to(target / "README.md")
    with pytest.raises(ValueError, match="symlinks"):
        publish_subset(symlinkstage, tmp_path / "no")


def test_pdb_invalid_inputs(tmp_path: Path) -> None:
    """Empty coordinate files fail rather than generating fallback exports."""
    f = tmp_path / "empty.pdb"
    f.write_text("END\n")
    with pytest.raises(ValueError, match="Empty"):
        pdb_atoms(f)


@pytest.mark.parametrize(
    ("metric", "before", "after", "expected"),
    [
        ("herschlag_phi", 179.9846956945411, -179.98993717721007, 0.025367128248834),
        ("herschlag_phi", -179.9, 179.9, 0.2),
        ("herschlag_phi", 80.0, -80.0, 160.0),
        ("a_bd", 179.9, -179.9, 359.8),
        ("d1", 3.2, 3.2005, 0.0005),
    ],
)
def test_roundtrip_metric_errors_respect_only_phi_periodicity(
    metric: str, before: float, after: float, expected: float
) -> None:
    """A signed-phi boundary crossing does not hide linear/angular changes."""
    from nyl12_covalent.control_inspection import _metric_error

    assert _metric_error(metric, before, after) == pytest.approx(expected, abs=1e-10)
