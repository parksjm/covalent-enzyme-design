"""Check native I1–I4 mapping and actual Ile/Thr rather than Val stereochemistry."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.rfd3_output_mapping import canonical_output_keys
from nyl12_covalent.rfd3_qualification_audit_v2 import coordinate_gates


def test_native_ligand_identity_and_rejected_numbering_guesses() -> None:
    """The serialized I1–I4 identifiers stay unchanged, without offset guessing."""
    assert canonical_output_keys(["I"] * 4, [1, 2, 3, 4], ["N"] * 4) == [
        ("I", p, "N") for p in range(1, 5)
    ]
    for positions in ([233, 234, 235, 236], [1, 2, 3, 236]):
        with pytest.raises(ValueError, match="I1-I4"):
            canonical_output_keys(["I"] * 4, positions, ["N"] * 4)
    with pytest.raises(ValueError, match="Duplicate"):
        canonical_output_keys(["I"] * 5, [1, 2, 3, 4, 4], ["N"] * 5)
    for chain, position, message in (("Z", 1, "chain"), ("E", 233, "protein")):
        with pytest.raises(ValueError, match=message):
            canonical_output_keys(
                ["I"] * 4 + [chain], [1, 2, 3, 4, position], ["N"] * 5
            )


def test_val_methyl_swap_allowed_and_true_ile_thr_cb_flips_rejected() -> None:
    """Val CB is nonstereogenic, while Ile/Thr CB sign changes remain flagged."""
    atoms, coordinates = [], []
    residues = (("A", 1, "VAL"), ("A", 2, "ILE"), ("A", 3, "THR"), ("B", 1, "ALA"))
    for ordinal, (chain, position, residue) in enumerate(residues):
        names = ["N", "CA", "C", "CB"]
        if residue in ("VAL", "ILE"):
            names += ["CG1", "CG2"]
        if residue == "THR":
            names += ["OG1", "CG2"]
        points = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 0, 1], [2, 0, 1], [1, -1, 1]]
        for name, point in zip(names, points, strict=False):
            atoms.append(
                {
                    "chain": chain,
                    "residue": position,
                    "name": name,
                    "residue_name": residue,
                }
            )
            coordinates.append(np.array(point, dtype=float) + ordinal * 20)
    source = np.array(coordinates)
    thresholds = {
        "fixed_max_A": 0.05,
        "mobile_CA_rmsd_min_A": 0.0,
        "mobile_CA_rmsd_max_A": 3.0,
        "peptide_min_A": 1.15,
        "peptide_max_A": 1.65,
        "clash_min_A": 0.1,
    }
    for position, expected in ((1, []), (2, ["A2:CB"]), (3, ["A3:CB"])):
        changed = source.copy()
        indices = [
            i
            for i, a in enumerate(atoms)
            if a["chain"] == "A"
            and a["residue"] == position
            and a["name"] in ("CG1", "CG2", "OG1")
        ]
        changed[indices] = changed[indices[::-1]]
        result = coordinate_gates(
            source, changed, atoms, [], ["A1", "A2", "A3"], [0, 1], thresholds
        )
        assert result["stereochemical_inversions"] == expected
