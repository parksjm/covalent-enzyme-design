#!/usr/bin/env python3
"""Build intact, protein-only 7YU0/7YU1 precursor candidates with PDBFixer."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable

LOGGER = logging.getLogger(__name__)
SEED = 20261120


def parse_seqres(lines: Iterable[str]) -> dict[str, list[str]]:
    """Read full deposited sequences for chains A and B."""
    sequences: dict[str, list[str]] = defaultdict(list)
    for line in lines:
        if line.startswith("SEQRES"):
            sequences[line[11:12]].extend(line[19:70].split())
    if set(sequences) != {"A", "B"}:
        raise ValueError("Expected SEQRES for exactly chains A and B")
    if any(len(seq) != 355 for seq in sequences.values()):
        raise ValueError("Expected full 355-residue precursor sequence per chain")
    return dict(sequences)


def seqres_records(sequences: dict[str, list[str]]) -> list[str]:
    """Format complete sequences as PDB SEQRES records."""
    rows = []
    for chain, sequence in sequences.items():
        for offset in range(0, len(sequence), 13):
            rows.append(
                f"SEQRES {offset // 13 + 1:3d} {chain} {len(sequence):4d}  "
                + " ".join(sequence[offset : offset + 13])
            )
    return rows


def prepare_original(
    source: Path, destination: Path
) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Write protein-only heavy atoms with altloc A preferred and full SEQRES."""
    lines = source.read_text().splitlines()
    sequences = parse_seqres(lines)
    records: dict[tuple[str, int, str], list[str]] = defaultdict(list)
    for line in lines:
        if not line.startswith("ATOM  ") or line[21:22] not in {"A", "B"}:
            continue
        if line[76:78].strip().upper() in {"H", "D"}:
            continue
        if line[26:27].strip():
            raise ValueError("Insertion-coded residue encountered")
        residue = int(line[22:26])
        chain = line[21:22]
        if not 1 <= residue <= 355:
            raise ValueError(f"Residue outside precursor sequence: {chain}:{residue}")
        if line[17:20].strip() != sequences[chain][residue - 1]:
            raise ValueError(f"ATOM/SEQRES mismatch at {chain}:{residue}")
        records[(chain, residue, line[12:16])].append(line)
    output_atoms: dict[str, list[str]] = {"A": [], "B": []}
    choices: dict[str, str] = {}
    for (chain, residue, atom), candidates in records.items():
        alt_a = [line for line in candidates if line[16:17] == "A"]
        blank = [line for line in candidates if line[16:17] == " "]
        pool = alt_a or blank or candidates
        selected = max(pool, key=lambda line: float(line[54:60] or 0.0))
        label = selected[16:17].strip()
        output_atoms[chain].append(selected[:16] + " " + selected[17:])
        if len(candidates) > 1:
            choices[f"{chain}:{residue}:{atom.strip()}"] = label or "blank"
    for chain, rows in output_atoms.items():
        rows.sort(key=lambda row: (int(row[22:26]), row[12:16]))
        if not rows:
            raise ValueError(f"No protein atoms found for chain {chain}")
    output = ["HEADER    INTACT NYLC PRECURSOR PROTEIN HEAVY ATOMS"]
    output.extend(seqres_records(sequences))
    for chain, rows in output_atoms.items():
        output.extend(rows)
        output.append(f"TER   {1:5d}      LYS {chain}{355:4d}")
    output.append("END")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(output) + "\n")
    return sequences, choices


def atoms(path: Path) -> dict[tuple[str, int, str], tuple[float, float, float]]:
    """Read heavy ATOM coordinates keyed by chain, residue number and atom."""
    found = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        key = (line[21:22], int(line[22:26]), line[12:16].strip())
        if key in found:
            raise ValueError(f"Duplicate atom identity: {key}")
        xyz = tuple(float(line[i : i + 8]) for i in (30, 38, 46))
        if not all(math.isfinite(x) for x in xyz):
            raise ValueError(f"Invalid coordinates at {key}")
        found[key] = xyz
    return found


def run(source: Path, outputdir: Path) -> dict[str, object]:
    """Generate a gap-built candidate; record geometry diagnostics without gating."""
    from openmm import Platform
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    outputdir.mkdir(parents=True, exist_ok=False)
    original = outputdir / "original_protein.pdb"
    candidate = outputdir / "candidate.pdb"
    sequences, altlocs = prepare_original(source, original)
    original_atoms = atoms(original)
    fixer = PDBFixer(filename=str(original))
    fixer.platform = Platform.getPlatformByName("CPU")
    fixer.findMissingResidues()
    missing_residues = {str(k): list(v) for k, v in fixer.missingResidues.items()}
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        raise ValueError(
            f"Unexpected nonstandard residues: {fixer.nonstandardResidues}"
        )
    fixer.findMissingAtoms()
    missing_atoms = sum(len(value) for value in fixer.missingAtoms.values())
    missing_terminal_atoms = sum(
        len(value) for value in fixer.missingTerminals.values()
    )
    fixer.addMissingAtoms(seed=SEED)
    with candidate.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    # PDBFile omits SEQRES; restore the complete precursor declaration.
    candidate_lines = candidate.read_text().splitlines()
    candidate.write_text(
        "\n".join(
            [candidate_lines[0], *seqres_records(sequences), *candidate_lines[1:]]
        )
        + "\n"
    )
    new_atoms = atoms(candidate)
    residue_names: dict[str, dict[int, str]] = defaultdict(dict)
    atom_names: dict[tuple[str, int], set[str]] = defaultdict(set)
    for line in candidate.read_text().splitlines():
        if line.startswith("ATOM  "):
            chain, number = line[21:22], int(line[22:26])
            residue_names[chain][number] = line[17:20].strip()
            atom_names[(chain, number)].add(line[12:16].strip())
    if set(residue_names) != {"A", "B"}:
        raise ValueError(f"Unexpected output chains: {sorted(residue_names)}")
    for chain, sequence in sequences.items():
        if set(residue_names[chain]) != set(range(1, 356)):
            raise ValueError(f"Incomplete 1..355 sequence in chain {chain}")
        if any(residue_names[chain][i] != sequence[i - 1] for i in range(1, 356)):
            raise ValueError(f"Output sequence differs from deposited chain {chain}")
        if "OXT" not in atom_names[(chain, 355)]:
            raise ValueError(f"Missing terminal OXT in chain {chain}")
    displacements = [
        math.dist(xyz, new_atoms[key]) for key, xyz in original_atoms.items()
    ]
    maximum_displacement = max(displacements, default=0.0)
    if len(new_atoms) < len(original_atoms) or maximum_displacement > 0.002:
        raise ValueError(
            f"Original heavy atom retention failed: max drift {maximum_displacement:.4f} A"
        )
    bad_cn = []
    for chain in "AB":
        for number in range(1, 355):
            c = new_atoms[(chain, number, "C")]
            n = new_atoms[(chain, number + 1, "N")]
            distance = math.dist(c, n)
            if not 1.1 <= distance <= 1.8:
                bad_cn.append(
                    {
                        "chain": chain,
                        "residues": [number, number + 1],
                        "distance_A": round(distance, 4),
                    }
                )
    report = {
        "candidate_only": True,
        "state": "candidate_built",
        "source": str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "seed": SEED,
        "platform": fixer.platform.getName(),
        "cpu_threads": 2,
        "altloc_policy": "prefer A per atom, then blank, then highest occupancy; selected labels recorded",
        "selected_altlocs": altlocs,
        "full_precursor_ranges": {chain: [1, 355] for chain in "AB"},
        "preserved_sequence_differences": [
            {
                "position": i,
                "7YU0": sequences["A"][i - 1],
                "7YU1": sequences["B"][i - 1],
            }
            for i in range(1, 356)
            if sequences["A"][i - 1] != sequences["B"][i - 1]
        ],
        "internal_residue_267": {chain: sequences[chain][266] for chain in "AB"},
        "missing_residue_groups": missing_residues,
        "missing_atom_count": missing_atoms,
        "missing_terminal_atom_count": missing_terminal_atoms,
        "original_heavy_atom_count": len(original_atoms),
        "candidate_heavy_atom_count": len(new_atoms),
        "max_original_heavy_atom_drift_A": maximum_displacement,
        "oxt_count": sum("OXT" in atom_names[(chain, 355)] for chain in "AB"),
        "bad_peptide_cn_distances": bad_cn,
        "geometry_acceptance": "deferred to explicit Amber19/GBn2 refinement stage",
    }
    (outputdir / "candidate_metadata.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    LOGGER.info("Built PDBFixer precursor candidate at %s", candidate)
    return report


def main() -> None:
    """Run the precursor candidate builder."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outputdir", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    print(json.dumps(run(args.source, args.outputdir), indent=2))


if __name__ == "__main__":
    main()
