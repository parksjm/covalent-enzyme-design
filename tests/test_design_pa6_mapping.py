"""Focused tests for four-copy PA6 conversion and canonical export glue."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import numpy as np
import pytest

from nyl12_covalent import design_pa6_mapping as mapping


def _atom_line(
    record: str,
    serial: int,
    name: str,
    residue: str,
    chain: str,
    resnum: int,
    element: str = "C",
    xyz: tuple[float, float, float] = (1.0, 2.0, 3.0),
) -> str:
    x, y, z = xyz
    return (
        f"{record:<6s}{serial:5d} {name:>4s} {residue:>3s} {chain}{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2s}"
    )


def _source(path: Path, include_hydrogen: bool = False) -> None:
    lines: list[str] = []
    serial = 0
    for chain, count in [(chain, 232) for chain in "ABCD"] + [
        (chain, 96) for chain in "EFGH"
    ]:
        for residue in range(1, count + 1):
            serial += 1
            lines.append(_atom_line("ATOM", serial, "CA", "ALA", chain, residue))
    for chain in "IJKL":
        for index in range(33):
            serial += 1
            lines.append(
                _atom_line(
                    "HETATM",
                    serial,
                    f"C{index:03d}",
                    "LIG",
                    chain,
                    9,
                    xyz=(float(index), 2.0, 3.0),
                )
            )
    if include_hydrogen:
        serial += 1
        lines.append(_atom_line("ATOM", serial, "H", "ALA", "A", 1, "H"))
        serial += 1
        lines.append(_atom_line("HETATM", serial, "H1", "LIG", "I", 1, "H"))
    path.write_text("\n".join([*lines, "END"]) + "\n", encoding="utf-8")


def _fake_converter(
    source: Path, mol2: Path, mapping_path: Path, release: Path, output: Path
) -> dict[str, Any]:
    source_lines = [
        line
        for line in source.read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    protein: list[str] = []
    for chain in "ABCDEFGH":
        offset = (
            "ABCD".index(chain) * 232
            if chain in "ABCD"
            else 928 + "EFGH".index(chain) * 96
        )
        for line in source_lines:
            if line[21] == chain:
                protein.append(f"{line[:22]}{offset + int(line[22:26]):4d}{line[26:]}")
    source_ligand = [line for line in source_lines if line[21] == "I"]
    ligand: list[str] = []
    records = []
    for index, line in enumerate(source_ligand):
        residue = 1313 + index % 4
        ligand_line = f"{line[:21]}I{residue:4d}{line[26:]}"
        ligand.append(ligand_line)
        records.append(
            {
                "source_atom": line[12:16].strip(),
                "source_serial": int(line[6:11]),
                "target_atom": line[12:16].strip(),
                "target_residue": residue,
            }
        )
    output.write_text("\n".join([*protein, *ligand, "END"]) + "\n")
    return {"atom_mapping": records}


def _setup_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    source = tmp_path / "source.pdb"
    _source(source, include_hydrogen=True)
    mol2 = tmp_path / "whole.mol2"
    mol2.write_text("pinned\n")
    mapping_path = tmp_path / "mapping.json"
    mapping_path.write_text(json.dumps({"blocks": []}))
    release = tmp_path / "release"
    release.mkdir()
    return source, mol2, mapping_path, release


def test_four_copy_offsets_ter_and_component_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source, mol2, context, release = _setup_inputs(tmp_path)
    monkeypatch.setattr(mapping, "convert_heavy_complex", _fake_converter)
    destination = tmp_path / "converted"

    result = mapping.convert_four_ligands(source, mol2, context, release, destination)

    lines = (destination / "complex_heavy_modular.pdb").read_text().splitlines()
    atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    ligand_atoms = [line for line in atoms if line[21] in "IJKL"]
    assert len(ligand_atoms) == 132
    assert {
        chain: sum(line[21] == chain for line in ligand_atoms) for chain in "IJKL"
    } == {chain: 33 for chain in "IJKL"}
    assert sorted({int(line[22:26]) for line in ligand_atoms}) == list(
        range(1313, 1329)
    )
    assert len({int(line[6:11]) for line in atoms}) == len(atoms)
    ters = [line for line in lines if line.startswith("TER   ")]
    assert [int(line[22:26]) for line in ters] == list(mapping._TER_RESIDUES)
    assert all(line[17:20] != "LIG" for line in atoms if line[21] in "ABCDEFGH")
    assert (destination / "components/I/component_input.pdb").is_file()
    component_input = (destination / "components/I/component_input.pdb").read_text()
    component_atoms = [
        line
        for line in component_input.splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    assert not any(line[76:78].strip() == "H" for line in component_atoms)
    assert {line[22:26] for line in component_atoms if line[21] == "I"} == {"   1"}
    assert len(result["atom_mapping"]) == 132
    assert result["atom_mapping"][0]["source_chain"] == "I"
    assert max(row["target_residue"] for row in result["atom_mapping"]) == 1328
    assert (
        json.loads((destination / "conversion.json").read_text())["ligand_copies"] == 4
    )


def test_strict_converter_failure_propagates_without_publishing_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source, mol2, context, release = _setup_inputs(tmp_path)
    destination = tmp_path / "failed"
    calls = 0

    def failing_converter(*args: Any) -> dict[str, Any]:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise ValueError("heavy graph mismatch")
        return _fake_converter(*args)

    monkeypatch.setattr(mapping, "convert_heavy_complex", failing_converter)
    with pytest.raises(ValueError, match="heavy graph mismatch"):
        mapping.convert_four_ligands(source, mol2, context, release, destination)
    assert not destination.exists()


@pytest.mark.parametrize(
    "mutation, message", [("altloc", "Alternate locations"), ("nan", "non-finite")]
)
def test_source_coordinates_and_altlocs_are_rejected(
    tmp_path: Path, mutation: str, message: str
) -> None:
    source = tmp_path / "invalid.pdb"
    _source(source)
    lines = source.read_text().splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith("HETATM"))
    line = lines[index]
    if mutation == "altloc":
        lines[index] = f"{line[:16]}A{line[17:]}"
    else:
        lines[index] = f"{line[:30]}     nan{line[38:]}"
    source.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match=message):
        mapping._validate_source(source)


def _fake_topology(atom_order: list[tuple[int, str, int]]) -> Any:
    atoms: list[Any] = []
    residues: list[Any] = []
    position = 0
    for res_index in range(1328):
        if res_index < 1312:
            resname = "HIE" if res_index == 0 else "ALA"
            pairs = [(6, "CA")]
        else:
            template = ("P6N", "P6I", "P6I", "P6C")[(res_index - 1312) % 4]
            pairs = [
                (6, name) for at_res, name, _ in atom_order if at_res == res_index + 1
            ]
            resname = template
        residue_atoms = []
        for atomic_number, name in pairs:
            atom = SimpleNamespace(idx=position, atomic_number=atomic_number, name=name)
            atoms.append(atom)
            residue_atoms.append(atom)
            position += 1
        residues.append(SimpleNamespace(name=resname, atoms=residue_atoms))
    topology = SimpleNamespace(
        atoms=atoms,
        residues=residues,
        coordinates=np.arange(position * 3, dtype=float).reshape(-1, 3),
    )
    return topology


def test_export_uses_atom_identity_mapping_not_topology_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    atom_order: list[tuple[int, str, int]] = []
    records: list[dict[str, Any]] = []
    atoms_per_residue = (9, 8, 8, 8)
    for residue in range(1313, 1329):
        local_residue = (residue - 1313) % 4
        atom_names = [
            f"T{index:03d}"
            for index in reversed(range(atoms_per_residue[local_residue]))
        ]
        for atom_name in atom_names:
            atom_order.append((residue, atom_name, 6))
            records.append(
                {
                    "target_residue": residue,
                    "target_atom": atom_name,
                    "source_atom": f"B{residue % 10}{atom_name[-2:]}",
                    "source_chain": "IJKL"[(residue - 1313) // 4],
                }
            )
    topology_object = _fake_topology(atom_order)
    amber = ModuleType("parmed.amber")
    amber.AmberParm = lambda *args, **kwargs: topology_object  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "parmed.amber", amber)
    topology = tmp_path / "system.prmtop"
    restart = tmp_path / "system.rst7"
    topology.write_text("topology")
    restart.write_text("coordinates")
    conversion = tmp_path / "conversion.json"
    conversion.write_text(json.dumps({"atom_mapping": records}))
    output = tmp_path / "canonical.pdb"

    result = mapping.export_canonical_heavy(topology, restart, conversion, output)

    ligand_lines = [
        line for line in output.read_text().splitlines() if line.startswith("HETATM")
    ]
    assert len(ligand_lines) == len(records)
    for line in ligand_lines:
        assert line[12:16].strip() in {
            row["source_atom"] for row in records if row["source_chain"] == line[21]
        }
    protein_first = next(
        line for line in output.read_text().splitlines() if line.startswith("ATOM  ")
    )
    assert protein_first[17:20] == "HIS"
    assert result["ligand_residues"] == 16


def test_export_rejects_wrong_residue_inventory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    amber = ModuleType("parmed.amber")
    amber.AmberParm = lambda *args, **kwargs: SimpleNamespace(  # type: ignore[attr-defined]
        coordinates=np.empty((0, 3)), atoms=[], residues=[]
    )
    monkeypatch.setitem(sys.modules, "parmed.amber", amber)
    topology = tmp_path / "system.prmtop"
    restart = tmp_path / "system.rst7"
    topology.write_text("topology")
    restart.write_text("coordinates")
    conversion = tmp_path / "conversion.json"
    records = []
    for residue in range(1313, 1329):
        for index in range((9, 8, 8, 8)[(residue - 1313) % 4]):
            records.append(
                {
                    "target_residue": residue,
                    "target_atom": f"T{index:03d}",
                    "source_atom": f"S{index:03d}",
                }
            )
    conversion.write_text(json.dumps({"atom_mapping": records}))
    with pytest.raises(ValueError, match="1328 total residues"):
        mapping.export_canonical_heavy(
            topology, restart, conversion, tmp_path / "out.pdb"
        )
