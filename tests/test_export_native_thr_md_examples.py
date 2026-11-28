"""Meaningful rigid-fit and portable atom-identity checks for MD overlays."""

from pathlib import Path

import numpy as np
import pytest

from scripts.export_native_thr_md_examples import align_site, write_pdb


def test_fit_excludes_thr_and_preserves_its_deviation() -> None:
    """Changing held-out Thr cannot change scaffold alignment or hide its shift."""
    rng = np.random.default_rng(17)
    ref = {
        (r, a): (a[0], rng.normal(size=3))
        for r in (113, 186, 272, 274)
        for a in ("N", "CA", "C", "O")
    }
    ref[(233, "OG1")] = ("O", np.array([1.0, 2.0, 3.0]))
    rot = np.array([[0.0, -1, 0], [1.0, 0, 0], [0, 0, 1]])
    shift = np.array([7.0, 8.0, -9.0])
    mobile = {key: (e, x @ rot + shift) for key, (e, x) in ref.items()}
    mobile[(233, "OG1")] = ("O", mobile[(233, "OG1")][1] + np.array([2.0, 0, 0]))
    aligned, lig, fit = align_site(
        ref, mobile, {"C54": np.array([1, 2, 3]) @ rot + shift}
    )
    assert fit["fit_RMSD_A"] < 1e-12
    assert np.linalg.norm(
        aligned[(233, "OG1")][1] - ref[(233, "OG1")][1]
    ) == pytest.approx(2)
    assert np.allclose(lig["C54"], [1, 2, 3])
    first = np.array(fit["rotation_row_vector"])
    mobile[(233, "OG1")] = ("O", np.array([100, 200, 300]))
    _, _, second = align_site(ref, mobile, {})
    assert np.array_equal(first, second["rotation_row_vector"])


def test_pdb_has_canonical_identity_and_refuses_overwrite(tmp_path: Path) -> None:
    """PDB viewer files preserve canonical residue/atom identities and originals."""
    path = tmp_path / "site.pdb"
    write_pdb(
        path, {(233, "OG1"): ("O", np.array([1.0, 2.0, 3.0]))}, {"C54": np.zeros(3)}
    )
    rows = path.read_text().splitlines()
    atom = next(row for row in rows if row.startswith("ATOM"))
    assert atom[12:16].strip() == "OG1"
    assert atom[21] == "A" and atom[22:26].strip() == "233"
    assert atom[76:78].strip() == "O"
    assert any(row.startswith("HETATM") and row[12:16].strip() == "C54" for row in rows)
    with pytest.raises(FileExistsError):
        write_pdb(path, {}, {})
