"""Strict heavy-atom conversion of canonical Nyl12/PA6 AF3 complexes."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any


def convert_heavy_complex(
    source: Path, mol2: Path, mapping: Path, release: Path, output: Path
) -> dict[str, Any]:
    """Map canonical chain-local protein and whole-L4 heavy atoms to Amber.

    Args:
        source: Canonical heavy-only A-H protein, I1 PA6 complex.
        mol2: Pinned whole-L4 atom names and bond graph.
        mapping: Validated whole-L4 context assignments.
        release: Validated modular release directory.
        output: New converted PDB destination.

    Returns:
        Atom correspondence, graph-distance checks, and input identities.
    """
    if output.exists():
        raise FileExistsError(output)
    lines = [
        x for x in source.read_text().splitlines() if x.startswith(("ATOM  ", "HETATM"))
    ]
    if any(x[16].strip() or x[26].strip() for x in lines):
        raise ValueError("Alternate locations or insertion codes are unsupported")
    if any(x[76:78].strip() not in {"C", "N", "O", "S"} for x in lines):
        raise ValueError("Expected explicitly identified protein/PA6 heavy atoms")
    ligand = [x for x in lines if x[21] == "I"]
    protein = [x for x in lines if x[21] in "ABCDEFGH"]
    if len(ligand) + len(protein) != len(lines):
        raise ValueError("Unexpected chain outside A-I")
    content = mol2.read_text().splitlines()
    atom_start = content.index("@<TRIPOS>ATOM") + 1
    bond_start = content.index("@<TRIPOS>BOND") + 1
    atoms = [x.split() for x in content[atom_start : bond_start - 1] if x.strip()]
    heavy = {int(x[0]): x for x in atoms if x[5].split(".")[0] != "H"}
    by_name = {x[12:16].strip(): x for x in ligand}
    if len(by_name) != 33 or set(by_name) != {x[1] for x in heavy.values()}:
        raise ValueError("PA6 heavy names differ from pinned whole-L4 MOL2")
    if any(x[17:20] != "LIG" or int(x[22:26]) != 1 for x in ligand):
        raise ValueError("Expected exactly LIG I1")
    for row in heavy.values():
        if by_name[row[1]][76:78].strip() != row[5].split(".")[0].upper():
            raise ValueError(f"Element mismatch for {row[1]}")
    expected_bonds: set[tuple[int, int]] = set()
    for line in content[bond_start:]:
        if line.startswith("@<TRIPOS>"):
            break
        fields = line.split()
        if fields and int(fields[1]) in heavy and int(fields[2]) in heavy:
            expected_bonds.add(tuple(sorted((int(fields[1]), int(fields[2])))))
    xyz = {
        index: tuple(
            float(by_name[row[1]][a:b]) for a, b in ((30, 38), (38, 46), (46, 54))
        )
        for index, row in heavy.items()
    }
    observed_bonds: set[tuple[int, int]] = set()
    distances = []
    for left in heavy:
        for right in heavy:
            if left >= right:
                continue
            distance = math.dist(xyz[left], xyz[right])
            cutoff = (
                1.65
                if "O" in {heavy[left][5].split(".")[0], heavy[right][5].split(".")[0]}
                else 1.85
            )
            if distance < cutoff:
                if distance < 1.0:
                    raise ValueError("PA6 has a sub-angstrom heavy contact")
                observed_bonds.add((left, right))
                distances.append(
                    {
                        "names": [heavy[left][1], heavy[right][1]],
                        "distance_angstrom": distance,
                    }
                )
    if observed_bonds != expected_bonds:
        raise ValueError(
            f"PA6 geometric graph mismatch: {observed_bonds ^ expected_bonds}"
        )
    blocks = json.loads(mapping.read_text())["blocks"]
    if len(blocks) != 4 or any(
        x["assignments"] != blocks[0]["assignments"] for x in blocks[1:]
    ):
        raise ValueError("Mapping must contain four identical whole-L4 assignments")
    assignments = {int(x["atom_index_1based"]): x for x in blocks[0]["assignments"]}
    templates = ("P6N", "P6I", "P6I", "P6C")
    mapped = {
        (int(assignments[i]["residue_index_1based"]), assignments[i]["atom_name"]): (
            i,
            by_name[row[1]],
        )
        for i, row in heavy.items()
    }
    converted = []
    records = []
    serial = 0
    for chain in "ABCDEFGH":
        chain_lines = [x for x in protein if x[21] == chain]
        expected = 232 if chain in "ABCD" else 96
        residues = list(dict.fromkeys(int(x[22:26]) for x in chain_lines))
        if residues != list(range(1, expected + 1)):
            raise ValueError(f"Unexpected residue inventory for chain {chain}")
        offset = (
            "ABCD".index(chain) * 232
            if chain in "ABCD"
            else 928 + "EFGH".index(chain) * 96
        )
        for line in chain_lines:
            serial += 1
            residue = offset + int(line[22:26])
            converted.append(
                f"{line[:6]}{serial:5d}{line[11:22]}{residue:4d}{line[26:]}"
            )
    for residue, template in enumerate(templates, 1):
        release_lines = (release / f"{template}.mol2").read_text().splitlines()
        start = release_lines.index("@<TRIPOS>ATOM") + 1
        end = release_lines.index("@<TRIPOS>BOND")
        keys = [
            (residue, x.split()[1])
            for x in release_lines[start:end]
            if x.strip() and not x.split()[1].startswith("H")
        ]
        for key in keys:
            index, line = mapped.pop(key)
            serial += 1
            converted.append(
                f"HETATM{serial:5d} {key[1]:<4s}{template:>4s} I{1312 + residue:4d}{line[26:]}"
            )
            records.append(
                {
                    "source_atom": heavy[index][1],
                    "source_serial": int(line[6:11]),
                    "target_atom": key[1],
                    "target_residue": 1312 + residue,
                }
            )
    if mapped:
        raise ValueError(f"Unmapped PA6 atoms: {mapped}")
    output.write_text("\n".join([*converted, "END"]) + "\n")
    return {
        "status": "accepted_raw_heavy_conversion",
        "protein_residues": 1312,
        "ligand_heavy_atoms": 33,
        "atom_mapping": records,
        "heavy_graph_bonds": distances,
        "coordinates_preserved": True,
        "input_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (source, mol2, mapping)
        },
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }
