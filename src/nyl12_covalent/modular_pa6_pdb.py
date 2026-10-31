"""Convert a legacy whole-residue PA6 L4 PDB to released modular residues."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class PdbAtom:
    """One PDB coordinate record required for deterministic remapping."""

    serial: int
    name: str
    residue_name: str
    chain: str
    residue_number: int
    xyz: tuple[float, float, float]
    occupancy: float
    temperature_factor: float
    element: str


def convert_legacy_pa6_l4_pdb(
    source_pdb: Path,
    source_mol2: Path,
    mapping_path: Path,
    release_dir: Path,
    output_dir: Path,
    *,
    ligand_residue_name: str = "LIG",
    ligand_residue_number: int = 1313,
) -> dict[str, Any]:
    """Split one 79-atom PA6 residue into P6N-P6I-P6I-P6C.

    The source PDB atom order and CONECT graph must exactly match the pinned
    whole-L4 MOL2. The accepted context mapping supplies residue-local names,
    while the release MOL2 files define output ordering.

    Args:
        source_pdb: Historical protein/PA6 complex with one whole LIG residue.
        source_mol2: Pinned 79-atom whole-L4 graph and atom ordering.
        mapping_path: Accepted whole-L4 residue assignment manifest.
        release_dir: Imported immutable modular PA6 release.
        output_dir: New immutable conversion directory.
        ligand_residue_name: Exact legacy ligand residue name.
        ligand_residue_number: Exact legacy ligand residue number.

    Returns:
        Conversion manifest with source and output hashes.

    Raises:
        FileExistsError: If the output directory already exists.
        ValueError: If atom, graph, mapping, or release identity checks fail.
    """
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    required = (source_pdb, source_mol2, mapping_path, release_dir / "manifest.json")
    if any(not path.is_file() for path in required):
        raise ValueError("One or more modular conversion inputs are missing")

    lines = source_pdb.read_text(encoding="ascii").splitlines()
    atoms = [
        _parse_pdb_atom(line) for line in lines if line.startswith(("ATOM  ", "HETATM"))
    ]
    ligand = [
        atom
        for atom in atoms
        if atom.residue_name == ligand_residue_name
        and atom.residue_number == ligand_residue_number
    ]
    mol2_names, mol2_bonds = _read_mol2_graph(source_mol2)
    if len(ligand) != 79 or [atom.name for atom in ligand] != mol2_names:
        raise ValueError("Legacy PDB does not match the pinned 79-atom MOL2 order")
    pdb_bonds = _pdb_ligand_bonds(lines, ligand)
    if pdb_bonds != mol2_bonds:
        raise ValueError("Legacy PDB CONECT graph does not match the pinned MOL2")

    assignments = _assignments(mapping_path)
    if [row["atom_index_1based"] for row in assignments] != list(range(1, 80)):
        raise ValueError("Mapping does not contain one ordered 79-atom block")
    release = json.loads((release_dir / "manifest.json").read_text(encoding="utf-8"))
    if not (
        release.get("md_ready") is True
        and release.get("production_ready_for_validated_l4_scope") is True
        and release.get("chemistry", {}).get("sequence") == ["P6N", "P6I", "P6I", "P6C"]
    ):
        raise ValueError("Imported release is not accepted for modular PA6 L4 MD")

    templates = ("P6N", "P6I", "P6I", "P6C")
    template_names = {
        name: _read_mol2_graph(release_dir / f"{name}.mol2")[0]
        for name in set(templates)
    }
    by_residue_name: dict[tuple[int, str], tuple[PdbAtom, dict[str, Any]]] = {}
    for atom, assignment in zip(ligand, assignments, strict=True):
        key = (int(assignment["residue_index_1based"]), str(assignment["atom_name"]))
        if key in by_residue_name:
            raise ValueError(f"Duplicate modular atom assignment: {key}")
        by_residue_name[key] = (atom, assignment)

    expected_keys = {
        (index, name)
        for index, template in enumerate(templates, start=1)
        for name in template_names[template]
    }
    if set(by_residue_name) != expected_keys:
        raise ValueError("Mapping atom names differ from release templates")

    ligand_serials = {atom.serial for atom in ligand}
    retained = _retained_pdb_lines(lines, ligand_serials, ligand_residue_number)
    max_serial = max(atom.serial for atom in atoms if atom.serial not in ligand_serials)
    chain = ligand[0].chain
    records: list[dict[str, Any]] = []
    output_lines = [*retained]
    serial = max_serial + 1
    for residue_index, template in enumerate(templates, start=1):
        for local_name in template_names[template]:
            source_atom, assignment = by_residue_name[(residue_index, local_name)]
            output_lines.append(
                _format_pdb_atom(
                    serial,
                    local_name,
                    template,
                    chain,
                    ligand_residue_number + residue_index - 1,
                    source_atom,
                )
            )
            records.append(
                {
                    "source_serial": source_atom.serial,
                    "source_atom_name": source_atom.name,
                    "target_serial": serial,
                    "target_residue_number": ligand_residue_number + residue_index - 1,
                    "target_residue_name": template,
                    "target_atom_name": local_name,
                    "mapping_equivalence_class": assignment["equivalence_class"],
                }
            )
            serial += 1
    output_lines.extend(
        (f"TER   {serial:5d}      P6C {chain}{ligand_residue_number + 3:4d}", "END")
    )

    output_dir.mkdir(parents=True)
    output_pdb = output_dir / "complex_modular_pa6_l4.pdb"
    output_pdb.write_text("\n".join(output_lines) + "\n", encoding="ascii")
    manifest = {
        "schema_version": 1,
        "status": "accepted_legacy_whole_l4_to_modular_pdb_conversion",
        "sequence": list(templates),
        "source_ligand": {
            "residue_name": ligand_residue_name,
            "residue_number": ligand_residue_number,
            "atom_count": len(ligand),
            "bond_count": len(pdb_bonds),
        },
        "checks": {
            "source_atom_order_matches_mol2": True,
            "source_conect_matches_mol2": True,
            "mapping_matches_release_templates": True,
            "coordinates_preserved": True,
            "modular_atom_count_79": len(records) == 79,
        },
        "atom_mapping": records,
        "source_sha256": {str(path.resolve()): _sha256(path) for path in required},
        "output_pdb": str(output_pdb.resolve()),
        "output_pdb_sha256": _sha256(output_pdb),
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def _parse_pdb_atom(line: str) -> PdbAtom:
    return PdbAtom(
        serial=int(line[6:11]),
        name=line[12:16].strip(),
        residue_name=line[17:20].strip(),
        chain=line[21:22],
        residue_number=int(line[22:26]),
        xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
        occupancy=float(line[54:60]),
        temperature_factor=float(line[60:66]),
        element=(line[76:78].strip() or line[12:16].strip()[0]).upper(),
    )


def _read_mol2_graph(path: Path) -> tuple[list[str], set[tuple[int, int]]]:
    lines = path.read_text(encoding="ascii").splitlines()
    atom_start = lines.index("@<TRIPOS>ATOM") + 1
    bond_start = lines.index("@<TRIPOS>BOND") + 1
    names = [
        line.split()[1] for line in lines[atom_start : bond_start - 1] if line.strip()
    ]
    bonds: set[tuple[int, int]] = set()
    for line in lines[bond_start:]:
        if line.startswith("@<TRIPOS>"):
            break
        fields = line.split()
        if fields:
            bonds.add(tuple(sorted((int(fields[1]), int(fields[2])))))
    return names, bonds


def _pdb_ligand_bonds(lines: list[str], ligand: list[PdbAtom]) -> set[tuple[int, int]]:
    serial_to_index = {atom.serial: index for index, atom in enumerate(ligand, start=1)}
    bonds: set[tuple[int, int]] = set()
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        fields = [int(value) for value in line.split()[1:]]
        if fields and fields[0] in serial_to_index:
            for target in fields[1:]:
                if target in serial_to_index:
                    bonds.add(
                        tuple(
                            sorted(
                                (serial_to_index[fields[0]], serial_to_index[target])
                            )
                        )
                    )
    return bonds


def _assignments(path: Path) -> list[dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    blocks = data.get("blocks", [])
    if len(blocks) != 4 or any(
        block.get("assignments") != blocks[0].get("assignments") for block in blocks[1:]
    ):
        raise ValueError("Mapping must contain four identical whole-L4 blocks")
    return list(blocks[0]["assignments"])


def _retained_pdb_lines(
    lines: list[str], ligand_serials: set[int], ligand_residue: int
) -> list[str]:
    retained = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and int(line[6:11]) in ligand_serials:
            continue
        if line.startswith("CONECT") or line.startswith("END"):
            continue
        if line.startswith("TER") and line[22:26].strip() == str(ligand_residue):
            continue
        retained.append(line)
    return retained


def _format_pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain: str,
    residue_number: int,
    source: PdbAtom,
) -> str:
    x, y, z = source.xyz
    return (
        f"HETATM{serial:5d} {atom_name:<4s}{residue_name:>4s} {chain}{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{source.occupancy:6.2f}"
        f"{source.temperature_factor:6.2f}          {source.element:>2s}  "
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
