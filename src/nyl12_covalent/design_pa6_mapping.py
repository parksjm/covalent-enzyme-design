"""Assemble and display four-ligand PA6 designvalidation structures."""

from __future__ import annotations

import hashlib
import json
import math
import tempfile
from pathlib import Path
from typing import Any

from .raw_pa6_modular import convert_heavy_complex

_SOURCE_CHAINS = "IJKL"
_LIGAND_TEMPLATES = ("P6N", "P6I", "P6I", "P6C")
_ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}
_TER_RESIDUES = (232, 464, 696, 928, 1024, 1120, 1216, 1312, 1316, 1320, 1324, 1328)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_source(source: Path) -> tuple[list[str], dict[str, list[str]]]:
    lines = [
        line
        for line in source.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    if not lines:
        raise ValueError("Source contains no PDB atoms")
    if any(line[16].strip() or line[26].strip() for line in lines):
        raise ValueError("Alternate locations or insertion codes are unsupported")
    if {line[21] for line in lines} != set("ABCDEFGHIJKL"):
        raise ValueError("Source must contain exactly chains A-L")
    by_chain: dict[str, list[str]] = {chain: [] for chain in "ABCDEFGHIJKL"}
    for line in lines:
        element = line[76:78].strip().upper()
        if element not in {"H", "C", "N", "O", "S"}:
            raise ValueError(f"Unsupported source element {element!r}")
        xyz = tuple(float(line[start : start + 8]) for start in (30, 38, 46))
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError("Source contains non-finite coordinates")
        by_chain[line[21]].append(line)
    for chain in _SOURCE_CHAINS:
        heavy = [line for line in by_chain[chain] if line[76:78].strip() != "H"]
        residue_ids = {(line[22:26], line[26]) for line in heavy}
        if len(heavy) != 33 or len(residue_ids) != 1:
            raise ValueError(f"Chain {chain} must contain one 33-heavy-atom ligand")
        if len({line[12:16].strip() for line in heavy}) != 33:
            raise ValueError(f"Chain {chain} has duplicate ligand heavy-atom names")
    return lines, by_chain


def _renumber(line: str, serial: int) -> str:
    return f"{line[:6]}{serial:5d}{line[11:]}"


def _ter_line(serial: int, chain: str, residue: int) -> str:
    return f"TER   {serial:5d}      {'':3s} {chain}{residue:4d}"


def convert_four_ligands(
    source: Path,
    mol2: Path,
    mapping: Path,
    release: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Convert a validated A-L complex into one modular four-ligand heavy PDB.

    Args:
        source: Boltz PDB containing protein chains A-H and PA6 chains I-L.
        mol2: Pinned whole-PA6 MOL2 graph used by the strict converter.
        mapping: Validated four-block whole-PA6 context mapping.
        release: Validated modular Amber template release.
        output_dir: Destination directory, which must not already exist.

    Returns:
        Conversion manifest including all 132 PA6 atom correspondences.
    """
    if output_dir.exists():
        raise FileExistsError(output_dir)
    source_lines, by_chain = _validate_source(source)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{output_dir.name}.", dir=output_dir.parent
    ) as tmp:
        stage = Path(tmp)
        components = stage / "components"
        components.mkdir()
        atom_mapping: list[dict[str, Any]] = []
        protein_lines: list[str] | None = None
        ligand_lines_by_chain: dict[str, list[str]] = {}
        for copy_index, source_chain in enumerate(_SOURCE_CHAINS):
            component = components / source_chain
            component.mkdir()
            input_path = component / "component_input.pdb"
            mapped_path = component / "converted.pdb"
            ligand_source = [
                line for line in by_chain[source_chain] if line[76:78].strip() != "H"
            ]
            normalized = [
                f"{line[:21]}I{1:4d}{line[26:]}"
                for line in by_chain[source_chain]
                if line[76:78].strip() != "H"
            ]
            normalized = [
                f"{line[:17]}LIG{line[20:]}" if line[21] == "I" else line
                for line in normalized
            ]
            input_path.write_text(
                "\n".join(
                    [
                        *[
                            line
                            for line in source_lines
                            if line[21] in "ABCDEFGH" and line[76:78].strip() != "H"
                        ],
                        *normalized,
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            (component / "whole_l4_mapping.json").write_bytes(mapping.read_bytes())
            source_name_serial = {
                line[12:16].strip(): int(line[6:11]) for line in ligand_source
            }
            conversion = convert_heavy_complex(
                input_path,
                mol2,
                component / "whole_l4_mapping.json",
                release,
                mapped_path,
            )
            (component / "conversion.json").write_text(
                json.dumps(conversion, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            converted_lines = [
                line
                for line in mapped_path.read_text(encoding="utf-8").splitlines()
                if line.startswith(("ATOM  ", "HETATM"))
            ]
            current_protein = [
                line for line in converted_lines if line[21] in "ABCDEFGH"
            ]
            current_ligand = [line for line in converted_lines if line[21] == "I"]
            if len(current_ligand) != 33:
                raise ValueError(
                    f"Strict converter returned an incomplete {source_chain} copy"
                )
            if protein_lines is None:
                protein_lines = current_protein
            elif current_protein != protein_lines:
                raise ValueError(
                    "Protein coordinates changed across component conversions"
                )
            ligand_lines_by_chain[source_chain] = current_ligand
            for record in conversion["atom_mapping"]:
                name = record["source_atom"]
                if name not in source_name_serial:
                    raise ValueError(f"Converter returned unknown source atom {name}")
                atom_mapping.append(
                    {
                        **record,
                        "source_chain": source_chain,
                        "source_serial": source_name_serial[name],
                        "component_source_serial": record["source_serial"],
                        "target_chain": "IJKL"[copy_index],
                        "target_residue": record["target_residue"] + 4 * copy_index,
                    }
                )
        if protein_lines is None or len(atom_mapping) != 132:
            raise ValueError("Expected four complete PA6 heavy-atom mappings")
        merged: list[str] = []
        serial = 0
        chain_ends = dict(zip("ABCDEFGH", _TER_RESIDUES[:8], strict=True))
        for chain in "ABCDEFGH":
            chain_lines = [line for line in protein_lines if line[21] == chain]
            if not chain_lines:
                raise ValueError(f"Converted protein chain {chain} is empty")
            for line in chain_lines:
                serial += 1
                merged.append(_renumber(line, serial))
            serial += 1
            merged.append(_ter_line(serial, chain, chain_ends[chain]))
        for copy_index, source_chain in enumerate(_SOURCE_CHAINS):
            target_chain = "IJKL"[copy_index]
            for line in ligand_lines_by_chain[source_chain]:
                serial += 1
                global_residue = int(line[22:26]) + 4 * copy_index
                merged.append(
                    _renumber(
                        f"{line[:21]}{target_chain}{global_residue:4d}{line[26:]}",
                        serial,
                    )
                )
            serial += 1
            merged.append(_ter_line(serial, target_chain, 1316 + 4 * copy_index))
        merged.append("END")
        pdb_path = stage / "complex_heavy_modular.pdb"
        pdb_path.write_text("\n".join(merged) + "\n", encoding="utf-8")
        manifest: dict[str, Any] = {
            "status": "accepted_four_ligand_heavy_conversion",
            "output": str(output_dir / pdb_path.name),
            "protein_residues": 1312,
            "ligand_copies": 4,
            "ligand_heavy_atoms_per_copy": 33,
            "atom_mapping": atom_mapping,
            "coordinates_preserved": True,
            "hydrogen_policy": "Source hydrogens stripped; Amber templates are authoritative for regeneration.",
            "histidine_view_policy": "Not applicable to the heavy conversion; canonical viewer export normalizes HIE/HID/HIP names to HIS only.",
            "source_sha256": _sha256(source),
            "input_sha256": {
                str(path): _sha256(path) for path in (source, mol2, mapping)
            },
            "output_sha256": _sha256(pdb_path),
        }
        (stage / "conversion.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        stage.rename(output_dir)
    return manifest


def export_canonical_heavy(
    topology: Path,
    restart: Path,
    conversion_path: Path,
    output: Path,
) -> dict[str, Any]:
    """Export Amber coordinates in canonical protein and Boltz ligand names.

    Args:
        topology: Amber topology for the four-ligand modular system.
        restart: Coordinates matching ``topology``.
        conversion_path: ``conversion.json`` from :func:`convert_four_ligands`.
        output: New heavy-only PDB destination.

    Returns:
        Export inventory and coordinate provenance.
    """
    if output.exists():
        raise FileExistsError(output)
    import numpy as np
    from parmed.amber import AmberParm

    manifest = json.loads(conversion_path.read_text(encoding="utf-8"))
    records = manifest.get("atom_mapping")
    if not isinstance(records, list) or len(records) != 132:
        raise ValueError("Conversion manifest must contain 132 PA6 atom mappings")
    ligand_names = {
        (int(row["target_residue"]), row["target_atom"]): row["source_atom"]
        for row in records
    }
    if len(ligand_names) != 132:
        raise ValueError("PA6 atom mapping contains duplicate target identities")
    model = AmberParm(str(topology), xyz=str(restart))
    coordinates = np.asarray(model.coordinates, dtype=float).reshape(-1, 3)
    if len(coordinates) != len(model.atoms) or not np.isfinite(coordinates).all():
        raise ValueError(
            "Amber topology/restart coordinates are incomplete or non-finite"
        )
    residue_atoms = list(model.residues)
    if len(residue_atoms) != 1328:
        raise ValueError(f"Expected 1328 total residues, found {len(residue_atoms)}")
    if any(len(residue.atoms) == 0 for residue in residue_atoms):
        raise ValueError("Topology contains an empty residue")
    expected_per_copy = ("P6N", "P6I", "P6I", "P6C")
    for index, name in enumerate(expected_per_copy * 4, start=1312):
        if residue_atoms[index].name.upper() != name:
            raise ValueError(f"Unexpected modular ligand residue at {index + 1}")
    protein_residue_counts = [len(residue.atoms) for residue in residue_atoms[:1312]]
    if len(protein_residue_counts) != 1312:
        raise ValueError("Protein residue inventory must contain 1312 residues")
    ligand_residues: dict[int, list[Any]] = {}
    for residue_index, residue in enumerate(residue_atoms[1312:], start=1313):
        ligand_residues[residue_index] = list(residue.atoms)
    ligand_heavy = 0
    for residue_index, atoms in ligand_residues.items():
        mapped_atoms = [atom for atom in atoms if int(atom.atomic_number) != 1]
        if not mapped_atoms:
            raise ValueError(f"Ligand residue {residue_index} has no heavy atoms")
        if len({atom.name for atom in mapped_atoms}) != len(mapped_atoms):
            raise ValueError(
                f"Ligand residue {residue_index} has duplicate Amber atom names"
            )
        for atom in mapped_atoms:
            if (residue_index, atom.name) not in ligand_names:
                raise ValueError(
                    f"Unmapped modular PA6 atom {residue_index}:{atom.name}"
                )
        ligand_heavy += len(mapped_atoms)
    if ligand_heavy != 132:
        raise ValueError(f"Expected 132 ligand heavy atoms, found {ligand_heavy}")
    output_lines: list[str] = []
    serial = 0
    ter_after: dict[str, int] = {}
    for global_residue, residue in enumerate(residue_atoms[:1312], start=1):
        chain, local = _protein_chain_residue(global_residue)
        residue_name = (
            "HIS"
            if residue.name.upper() in {"HIE", "HID", "HIP"}
            else residue.name[:3].upper()
        )
        for atom in residue.atoms:
            if int(atom.atomic_number) == 1:
                continue
            element = _ELEMENTS.get(int(atom.atomic_number))
            if element is None or element == "H":
                raise ValueError(
                    f"Unsupported protein atomic number {atom.atomic_number}"
                )
            serial += 1
            x, y, z = coordinates[atom.idx]
            output_lines.append(
                _pdb_atom_line(
                    "ATOM",
                    serial,
                    atom.name,
                    residue_name,
                    chain,
                    local,
                    x,
                    y,
                    z,
                    element,
                )
            )
        if global_residue in _TER_RESIDUES[:8]:
            ter_after[chain] = local
            serial += 1
            output_lines.append(_ter_line(serial, chain, local))
    for global_residue, atoms in ligand_residues.items():
        copy_index = (global_residue - 1313) // 4
        target_chain = "IJKL"[copy_index]
        for atom in atoms:
            if int(atom.atomic_number) == 1:
                continue
            element = _ELEMENTS.get(int(atom.atomic_number))
            if element is None or element == "H":
                raise ValueError(f"Unsupported PA6 atomic number {atom.atomic_number}")
            serial += 1
            x, y, z = coordinates[atom.idx]
            output_lines.append(
                _pdb_atom_line(
                    "HETATM",
                    serial,
                    ligand_names[(global_residue, atom.name)],
                    "LIG",
                    target_chain,
                    1,
                    x,
                    y,
                    z,
                    element,
                )
            )
        if global_residue in {1316, 1320, 1324, 1328}:
            serial += 1
            output_lines.append(_ter_line(serial, target_chain, 1))
    output_lines.append("END")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    return {
        "status": "accepted_canonical_heavy_export",
        "output": str(output),
        "protein_residues": 1312,
        "ligand_residues": 16,
        "ligand_heavy_atoms": ligand_heavy,
        "histidine_view_normalization": "HIE/HID/HIP displayed as HIS; topology protonation states unchanged.",
        "topology_sha256": _sha256(topology),
        "restart_sha256": _sha256(restart),
        "output_sha256": _sha256(output),
    }


def _protein_chain_residue(global_residue: int) -> tuple[str, int]:
    """Map concatenated Amber protein residue numbers back to source chains."""
    if 1 <= global_residue <= 928:
        index, offset = divmod(global_residue - 1, 232)
        return "ABCD"[index], offset + 1
    if 929 <= global_residue <= 1312:
        index, offset = divmod(global_residue - 929, 96)
        return "EFGH"[index], offset + 1
    raise ValueError(f"Protein residue out of range: {global_residue}")


def _pdb_atom_line(
    record: str,
    serial: int,
    atom: str,
    residue: str,
    chain: str,
    number: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    """Format one atom record using fixed-width PDB fields."""
    if not all(math.isfinite(value) for value in (x, y, z)):
        raise ValueError("Cannot write non-finite PDB coordinates")
    return (
        f"{record:<6s}{serial:5d} {atom:>4s} {residue:>3s} {chain}"
        f"{number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}"
        f"          {element:>2s}"
    )
