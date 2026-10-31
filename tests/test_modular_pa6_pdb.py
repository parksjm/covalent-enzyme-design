"""Tests for graph-verified legacy PA6 PDB conversion."""

from __future__ import annotations

from pathlib import Path

from nyl12_covalent.modular_pa6_pdb import _read_mol2_graph, convert_legacy_pa6_l4_pdb

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "inputs/pa6_modular_md_20261031"
RELEASE = ROOT / "parameters/releases/pa6_l4_modular_resp_gaff2_v1_20261031"


def test_convert_legacy_pdb_preserves_coordinates_and_splits_residues(
    tmp_path: Path,
) -> None:
    """A full synthetic legacy PDB converts to the exact modular inventory."""
    source = INPUTS / "pa6_l4_source.mol2"
    names, bonds = _read_mol2_graph(source)
    atom_lines = source.read_text(encoding="ascii").splitlines()
    start = atom_lines.index("@<TRIPOS>ATOM") + 1
    end = atom_lines.index("@<TRIPOS>BOND")
    rows = [line.split() for line in atom_lines[start:end] if line.strip()]
    pdb_lines = [
        "HEADER    SYNTHETIC MODULAR CONVERSION TEST",
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  ",
        "TER       2      ALA A   1",
    ]
    serial_offset = 19892
    for serial, row in enumerate(rows, start=serial_offset + 1):
        name, x, y, z = row[1], float(row[2]), float(row[3]), float(row[4])
        element = row[5][0].upper()
        pdb_lines.append(
            f"HETATM{serial:5d} {name:<4s} LIG A1313    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00          {element:>2s}  "
        )
    for left, right in sorted(bonds):
        pdb_lines.append(f"CONECT{left + serial_offset:5d}{right + serial_offset:5d}")
    pdb_lines.append("END")
    pdb = tmp_path / "legacy.pdb"
    pdb.write_text("\n".join(pdb_lines) + "\n", encoding="ascii")

    output = tmp_path / "converted"
    report = convert_legacy_pa6_l4_pdb(
        pdb,
        source,
        INPUTS / "whole_l4_context_mapping.json",
        RELEASE,
        output,
    )

    assert report["status"] == "accepted_legacy_whole_l4_to_modular_pdb_conversion"
    assert all(report["checks"].values())
    assert len(report["atom_mapping"]) == len(names) == 79
    converted = (output / "complex_modular_pa6_l4.pdb").read_text(encoding="ascii")
    atom_records = "\n".join(
        line for line in converted.splitlines() if line.startswith("HETATM")
    )
    assert atom_records.count(" P6N A1313") == 21
    assert atom_records.count(" P6I A1314") == 19
    assert atom_records.count(" P6I A1315") == 19
    assert atom_records.count(" P6C A1316") == 20
    assert "CONECT" not in converted
