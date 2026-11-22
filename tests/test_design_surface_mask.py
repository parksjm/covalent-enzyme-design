"""Tests for SASA-based design mask eligibility."""

from pathlib import Path

import pytest
from nyl12_covalent.design_surface_mask import (
    calculate_residue_sasa,
    passes_all_copy_sasa,
    proposal_shell_positions,
    protein_molecule_from_pdb,
    read_bound_complex_heavy_coordinates,
    relative_accessibility,
    summarize_sasa,
)


def _atom_line(
    serial: int,
    atom: str,
    residue: str,
    chain: str,
    number: int,
    xyz: tuple[float, float, float],
    element: str,
) -> str:
    """Format one fixed-column PDB record for an isolated test residue."""
    x, y, z = xyz
    return (
        f"ATOM  {serial:5d} {atom:>4s} {residue:>3s} {chain}{number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          "
        f"{element:>2s}  "
    )


def _write_alanine(path: Path) -> None:
    """Write a nonbonded five-heavy-atom alanine coordinate fixture."""
    atoms = [
        ("N", (0.0, 0.0, 0.0), "N"),
        ("CA", (1.45, 0.0, 0.0), "C"),
        ("C", (2.0, 1.4, 0.0), "C"),
        ("O", (1.4, 2.45, 0.0), "O"),
        ("CB", (2.0, -0.8, 1.2), "C"),
    ]
    lines = [
        _atom_line(i, name, "ALA", "A", 1, xyz, element)
        for i, (name, xyz, element) in enumerate(atoms, start=1)
    ]
    lines.append(
        "HETATM    6  C1  LIG I   1       8.000   0.000   0.000  1.00  0.00           C  "
    )
    path.write_text("\n".join(lines) + "\nEND\n")


def test_protein_parser_excludes_ligand_and_does_not_infer_bonds(
    tmp_path: Path,
) -> None:
    """Only protein heavy ATOM records are loaded without inferred bonds."""
    pdb = tmp_path / "alanine.pdb"
    _write_alanine(pdb)

    molecule = protein_molecule_from_pdb(pdb)

    assert molecule.GetNumAtoms() == 5
    assert molecule.GetNumBonds() == 0
    assert all(
        atom.GetPDBResidueInfo().GetResidueName().strip() == "ALA"
        for atom in molecule.GetAtoms()
    )


def test_atom_sasa_and_tien_normalization_are_positive(tmp_path: Path) -> None:
    """FreeSASA writes per-atom areas and residue ASA normalizes to Tien table."""
    pdb = tmp_path / "alanine.pdb"
    _write_alanine(pdb)
    molecule = protein_molecule_from_pdb(pdb)

    residue_asa = calculate_residue_sasa(molecule)
    report = summarize_sasa(molecule, expected_residues=1)

    asa = residue_asa[("A", 1, "ALA")]
    assert asa > 0
    assert report["algorithm"] == "Lee-Richards"
    assert report["classifier"] == "Protor"
    assert report["residues"][0]["rsa"] == pytest.approx(asa / 129.0, abs=1e-6)
    assert relative_accessibility("ALA", 129.0) == pytest.approx(1.0)


def test_sasa_rejects_residue_count_mismatch(tmp_path: Path) -> None:
    """The parent audit can require the exact standard-protein residue count."""
    pdb = tmp_path / "alanine.pdb"
    _write_alanine(pdb)

    with pytest.raises(ValueError, match="Expected 2 standard residues"):
        summarize_sasa(protein_molecule_from_pdb(pdb), expected_residues=2)


def test_all_equivalent_copies_must_pass_surface_gate() -> None:
    """One buried copy blocks a candidate even when the other copies are exposed."""
    assert passes_all_copy_sasa([0.25] * 8)
    assert not passes_all_copy_sasa([0.5] * 7 + [0.2499])
    with pytest.raises(ValueError, match="nonempty"):
        passes_all_copy_sasa([])


def test_bound_complex_coordinate_reader_requires_expected_chains(
    tmp_path: Path,
) -> None:
    """Protein and ligand coordinates are separated by chain and record type."""
    lines = [
        _atom_line(i, "CA", "GLY", chr(64 + i), 1, (float(i), 0.0, 0.0), "C")
        for i in range(1, 9)
    ]
    for serial, chain in enumerate("IJKL", start=9):
        atom = _atom_line(serial, "C1", "LIG", chain, 1, (float(serial), 0.0, 0.0), "C")
        lines.append("HETATM" + atom[6:])
    pdb = tmp_path / "complex.pdb"
    pdb.write_text("\n".join(lines) + "\nEND\n")

    protein, ligands = read_bound_complex_heavy_coordinates(pdb)

    assert set(protein) == set("ABCDEFGH")
    assert set(ligands) == set("IJKL")
    assert protein["A"][1][0] == "GLY"
    assert ligands["I"].shape == (1, 3)


def test_proposal_shell_reader_selects_seed_and_distance_shell() -> None:
    """Shell candidates resolve to numeric residue positions from the proposal."""
    proposal = {
        "mask_proposal": {
            "shell_geometry_and_candidates": {
                "alpha": {
                    "seed_02": {
                        "candidate_6_12_positions": [{"residue": "ALA57"}],
                        "distant_12_16_unselected": ["GLY82"],
                    },
                    "seed_03": {
                        "candidate_6_12_positions": [{"residue": "THR58"}],
                        "distant_12_16_unselected": ["ILE83"],
                    },
                }
            }
        }
    }

    assert proposal_shell_positions(proposal, 2026112202, "alpha", "6-12") == {57}
    assert proposal_shell_positions(proposal, 2026112203, "alpha", "12-16") == {83}
    with pytest.raises(ValueError, match="Unsupported seed/shell"):
        proposal_shell_positions(proposal, 1, "alpha", "6-12")
