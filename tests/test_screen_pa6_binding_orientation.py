"""Regression tests for static PA6 orientation screening."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from scripts import screen_pa6_binding_orientation as screen

SEED = (
    screen.ROOT
    / "qmmm/pa6_static_es_seed"
    / "pa6_l4z_l2_nyl12_wt_model_23_canonical.pdb"
)


@pytest.fixture
def atoms() -> dict:
    """Read the retained model23 scientific regression fixture."""
    return screen.read_atoms(SEED)


@pytest.fixture
def reference() -> dict:
    """Read the catalytic reference with the production residue selection."""
    return screen._cif_atoms(
        screen.ROOT / "cif/n12c.cif", "A", {r for _, r, _ in screen.SITE_RESIDUES}
    )


def test_model23_signed_geometry_regression(atoms: dict, reference: dict) -> None:
    """Preserve known signed phi, Tier A and held-out native Thr geometry."""
    site = screen.layout_site(atoms, 0)
    ligand = {name: atoms[("I", 1, name)][1] for name in ("C54", "O49", "N59", "C62")}
    row = screen.measure(site, ligand, reference)
    assert row["herschlag_phi"] == pytest.approx(99.01, abs=0.02)
    assert row["tier"] == "A"
    assert row["orientation"] == "positive_phi"
    assert row["native_pass"] and row["joint_pass"]
    assert row["og1_shift"] == pytest.approx(0.18836, abs=0.001)
    assert row["psi_delta"] == pytest.approx(-30.0904, abs=0.001)
    assert row["selection_eligible"] is True


def test_native_failure_cannot_veto_valid_catalytic_geometry(
    atoms: dict, reference: dict, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Keep the pose while recording native similarity as a failed diagnostic."""
    original = screen.state_metrics

    def changed_native(ref: dict, aligned: dict) -> dict:
        result = original(ref, aligned)
        result["thr233_held_out"]["heavy_atom_displacements_angstrom"]["OG1"] = 3.0
        return result

    monkeypatch.setattr(screen, "state_metrics", changed_native)
    ligand = {name: atoms[("I", 1, name)][1] for name in ("C54", "O49", "N59", "C62")}
    row = screen.measure(screen.layout_site(atoms, 0), ligand, reference)
    assert row["tier"] == "A"
    assert row["native_pass"] is False and row["joint_pass"] is False
    assert row["selection_eligible"] is True


@pytest.mark.parametrize("layout", ["split", "flat_amber"])
@pytest.mark.parametrize("site_index", range(4))
def test_alternative_layouts_preserve_site_atoms(
    atoms: dict, layout: str, site_index: int
) -> None:
    """Map each symmetry site identically across protein numbering schemes."""
    remapped = {}
    for (chain, auth, name), value in atoms.items():
        if chain not in "ABCD":
            continue
        index = "ABCD".index(chain)
        if layout == "split":
            target = (
                chain if auth <= 232 else "EFGH"[index],
                auth if auth <= 232 else auth - 232,
                name,
            )
        else:
            number = (
                auth + 232 * index if auth <= 232 else 929 + 96 * index + auth - 233
            )
            target = ("A", number, name)
        remapped[target] = value
    expected = screen.layout_site(atoms, site_index)
    actual = screen.layout_site(remapped, site_index)
    assert actual.keys() == expected.keys()
    for key in expected:
        np.testing.assert_array_equal(actual[key][1], expected[key][1])


def test_invalid_connectivity_fails(atoms: dict, reference: dict) -> None:
    """Reject an atom mapping that cannot represent the central amide."""
    ligand = {
        name: atoms[("I", 1, name)][1].copy() for name in ("C54", "O49", "N59", "C62")
    }
    ligand["O49"] = ligand["C54"] + np.array([20.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="connectivity inconsistent"):
        screen.measure(screen.layout_site(atoms, 0), ligand, reference)


@pytest.mark.parametrize("kind", ["duplicate", "alternate", "nonfinite"])
def test_ambiguous_or_invalid_coordinates_fail(tmp_path: Path, kind: str) -> None:
    """Reject duplicate keys, unsupported alternate positions and NaN inputs."""
    line = next(
        line for line in SEED.read_text().splitlines() if line.startswith("ATOM")
    )
    if kind == "duplicate":
        content = f"{line}\n{line}\n"
        message = "Ambiguous atom"
    elif kind == "alternate":
        content = line[:16] + "B" + line[17:] + "\n"
        message = "Ambiguous atom"
    else:
        content = line[:30] + "     nan" + line[38:] + "\n"
        message = "Nonfinite atom"
    path = tmp_path / "invalid.pdb"
    path.write_text(content)
    with pytest.raises(ValueError, match=message):
        screen.read_atoms(path)


def test_unknown_layout_fails() -> None:
    """Fail explicitly when no supported numbering marker is present."""
    with pytest.raises(ValueError, match="Unrecognized protein layout"):
        screen.layout_site({}, 0)


def test_manifest_integration(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Screen every site and retain input hashes alongside machine-readable rows."""
    manifest = tmp_path / "manifest.csv"
    with manifest.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "atom_names", "ensemble"])
        writer.writeheader()
        writer.writerow(
            {
                "path": str(SEED),
                "atom_names": "C54;O49;N59;C62",
                "ensemble": "retained_seed",
            }
        )
    output = tmp_path / "screen"
    monkeypatch.setattr("sys.argv", ["screen", str(manifest), str(output)])
    screen.main()
    summary = json.loads((output / "summary.json").read_text())
    rows = list(csv.DictReader((output / "all_sites.csv").open()))
    assert summary["models"] == 1
    assert summary["site_ligand_pairs"] == len(rows) == 4
    assert {row["site"] for row in rows} == set("ABCD")
    assert all(row["ensemble"] == "retained_seed" for row in rows)
    assert any(row["joint_pass"] for row in summary["joint_passes"])
    hashes = {item["path"]: item["sha256"] for item in summary["input_hashes"]}
    assert hashes[str(SEED)] == hashlib.sha256(SEED.read_bytes()).hexdigest()
    assert hashes[str(manifest)] == hashlib.sha256(manifest.read_bytes()).hexdigest()
    with pytest.raises(FileExistsError):
        screen.main()
