#!/usr/bin/env python3
"""Repair unresolved residues in mature NylC 5XYG after an explicit chain split."""

from __future__ import annotations

import argparse
import json
import logging
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable

LOGGER = logging.getLogger(__name__)
CHAIN_MAP = {
    ("A", "alpha"): "A",
    ("B", "alpha"): "B",
    ("A", "beta"): "C",
    ("B", "beta"): "D",
}
CUT_RESIDUE = 266
EXPECTED_RANGES = {"A": (1, 266), "B": (1, 266), "C": (267, 355), "D": (267, 355)}


def parse_seqres(lines: Iterable[str]) -> dict[str, list[str]]:
    """Return deposited three-letter SEQRES sequences by chain."""
    sequences: dict[str, list[str]] = defaultdict(list)
    for line in lines:
        if line.startswith("SEQRES"):
            chain_id = line[11:12]
            sequences[chain_id].extend(line[19:70].split())
    if set(sequences) != {"A", "B"} or sequences["A"] != sequences["B"]:
        raise ValueError(
            "5XYG must contain identical SEQRES records for chains A and B"
        )
    if len(sequences["A"]) != 355:
        raise ValueError(
            f"Expected 355 deposited residues per chain; got {len(sequences['A'])}"
        )
    return dict(sequences)


def select_protein_atoms(
    lines: Iterable[str],
) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Split protein ATOM records into mature alpha/beta chains and select one altloc."""
    source_lines = list(lines)
    sequences = parse_seqres(source_lines)
    candidates: dict[tuple[str, int, str, str], list[str]] = defaultdict(list)
    for line in source_lines:
        if not line.startswith("ATOM  "):
            continue
        source_chain = line[21:22]
        if source_chain not in {"A", "B"}:
            raise ValueError(f"Unexpected protein chain {source_chain!r}")
        if line[26:27].strip():
            raise ValueError(
                "Insertion-coded residue is not supported by this preparation"
            )
        if line[76:78].strip().upper() in {"H", "D"}:
            continue
        residue_number = int(line[22:26])
        if residue_number < 1 or residue_number > 355:
            raise ValueError(
                f"Residue outside deposited SEQRES: {source_chain}:{residue_number}"
            )
        if line[17:20].strip() != sequences[source_chain][residue_number - 1]:
            raise ValueError(
                f"ATOM residue does not match deposited SEQRES: {source_chain}:{residue_number}"
            )
        key = (source_chain, residue_number, line[26:27], line[12:16])
        candidates[key].append(line)

    selected: dict[str, list[str]] = {chain_id: [] for chain_id in "ABCD"}
    choices: dict[str, str] = {}
    for key, records in candidates.items():
        source_chain, residue_number, _, atom_name = key
        blank = [record for record in records if record[16:17] in {"", " "}]
        alt_a = [record for record in records if record[16:17] == "A"]
        pool = blank or alt_a or records
        chosen = max(
            pool, key=lambda rec: (float(rec[54:60] or 0.0), rec[16:17] == "A")
        )
        new_chain = (
            "AB"[source_chain == "B"]
            if residue_number <= CUT_RESIDUE
            else "CD"[source_chain == "B"]
        )
        # Once selected, remove the alternate-location label to make the chosen model explicit.
        rewritten = chosen[:16] + " " + chosen[17:21] + new_chain + chosen[22:]
        selected[new_chain].append(rewritten)
        if len(records) > 1:
            choices[f"{source_chain}:{residue_number}:{atom_name.strip()}"] = (
                f"altloc={chosen[16:17].strip() or 'blank'} occupancy={float(chosen[54:60] or 0.0):.2f}"
            )
    for chain_id, rows in selected.items():
        rows.sort(key=lambda row: (int(row[22:26]), row[12:16], row[16:17]))
        if not rows:
            raise ValueError(
                f"No protein ATOM records found for output chain {chain_id}"
            )
    return selected, choices


def write_split_input(source: Path, destination: Path) -> dict[str, object]:
    """Write protein-only A/B alpha and C/D beta PDBFixer input with SEQRES."""
    lines = source.read_text().splitlines()
    sequences = parse_seqres(lines)
    sequence = sequences["A"]
    selected, choices = select_protein_atoms(lines)
    output = ["HEADER    MATURE NYLC 5XYG SPLIT FOR RESIDUE REPAIR"]
    for chain_id in "ABCD":
        lo, hi = EXPECTED_RANGES[chain_id]
        chain_sequence = sequence[lo - 1 : hi]
        for offset in range(0, len(chain_sequence), 13):
            serial = offset // 13 + 1
            chunk = chain_sequence[offset : offset + 13]
            output.append(
                f"SEQRES {serial:3d} {chain_id} {len(chain_sequence):4d}  "
                + " ".join(chunk)
            )
        output.extend(selected[chain_id])
        last = hi
        last_resname = sequence[hi - 1]
        output.append(f"TER   {1:5d}      {last_resname:3s} {chain_id}{last:4d}")
    output.append("END")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(output) + "\n")
    return {
        "selected_altlocs": choices,
        "input_atom_counts": {key: len(value) for key, value in selected.items()},
        "expected_sequences": {
            chain_id: sequence[lo - 1 : hi]
            for chain_id, (lo, hi) in EXPECTED_RANGES.items()
        },
    }


def read_atoms(
    path: Path,
) -> dict[tuple[str, int, str, str], tuple[float, float, float]]:
    """Read ATOM records keyed by chain, author residue number, insertion code, and atom."""
    atoms: dict[tuple[str, int, str, str], tuple[float, float, float]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        if line[76:78].strip().upper() in {"H", "D"}:
            raise ValueError(
                "Unexpected hydrogen in protein-only heavy-atom input/output"
            )
        key = (line[21:22], int(line[22:26]), line[26:27], line[12:16].strip())
        if key in atoms:
            raise ValueError(f"Duplicate prepared atom identity: {key}")
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"Nonfinite coordinate in prepared structure: {key}")
        atoms[key] = xyz
    return atoms


def validate_repair(
    source_atoms: dict[tuple[str, int, str, str], tuple[float, float, float]],
    repaired: Path,
    expected_sequences: dict[str, list[str]],
    topology: object,
) -> dict[str, object]:
    """Check chain sequences, terminal atoms, peptide geometry, and coordinate retention."""
    output_atoms = read_atoms(repaired)
    moved = []
    for key, xyz in source_atoms.items():
        if key not in output_atoms:
            raise ValueError(f"Original heavy atom missing after repair: {key}")
        displacement = math.dist(xyz, output_atoms[key])
        if displacement > 0.002:
            moved.append((key, displacement))
    if moved:
        raise ValueError(
            f"Original heavy atom displacement exceeds 0.002 A: {moved[:5]}"
        )

    residues: dict[str, dict[int, str]] = defaultdict(dict)
    atoms_by_residue: dict[tuple[str, int], set[str]] = defaultdict(set)
    for line in repaired.read_text().splitlines():
        if line.startswith("ATOM  "):
            chain_id, residue_number, residue_name, atom_name = (
                line[21:22],
                int(line[22:26]),
                line[17:20],
                line[12:16].strip(),
            )
            residues[chain_id][residue_number] = residue_name
            atoms_by_residue[(chain_id, residue_number)].add(atom_name)
    if set(residues) != set("ABCD"):
        raise ValueError(f"Unexpected repaired chain IDs: {sorted(residues)}")
    for chain_id, (lo, hi) in EXPECTED_RANGES.items():
        expected_ids = set(range(lo, hi + 1))
        if set(residues[chain_id]) != expected_ids:
            raise ValueError(
                f"Chain {chain_id} residue inventory differs from {lo}-{hi}"
            )
        expected = expected_sequences[chain_id]
        actual = [residues[chain_id][i] for i in range(lo, hi + 1)]
        if actual != expected:
            raise ValueError(f"Chain {chain_id} does not match deposited SEQRES")
        if (
            "N" not in atoms_by_residue[(chain_id, lo)]
            or "C" not in atoms_by_residue[(chain_id, hi)]
        ):
            raise ValueError(
                f"Missing declared terminal backbone atoms in chain {chain_id}"
            )
        if "OXT" not in atoms_by_residue[(chain_id, hi)]:
            raise ValueError(f"Missing C-terminal OXT in chain {chain_id}")
        for residue_number in range(lo, hi):
            first = atoms_by_residue[(chain_id, residue_number)]
            second = atoms_by_residue[(chain_id, residue_number + 1)]
            if "C" not in first or "N" not in second:
                raise ValueError(
                    f"Missing peptide atoms at {chain_id}:{residue_number}-{residue_number + 1}"
                )
            c_xyz = output_atoms[(chain_id, residue_number, " ", "C")]
            n_xyz = output_atoms[(chain_id, residue_number + 1, " ", "N")]
            bond = math.dist(c_xyz, n_xyz)
            if not 1.1 <= bond <= 1.8:
                raise ValueError(
                    f"Nonphysical C-N distance {bond:.3f} A at {chain_id}:{residue_number}-{residue_number + 1}"
                )
    cleavage_pairs = {("A", "C", 266, 267), ("B", "D", 266, 267)}
    cleavage_bonds = []
    for atom1, atom2 in topology.bonds():
        for left, right in ((atom1, atom2), (atom2, atom1)):
            if (
                (
                    left.residue.chain.id,
                    right.residue.chain.id,
                    int(left.residue.id),
                    int(right.residue.id),
                )
                in cleavage_pairs
                and left.name == "C"
                and right.name == "N"
            ):
                cleavage_bonds.append((left.residue.chain.id, right.residue.chain.id))
    if cleavage_bonds:
        raise ValueError(
            f"Unexpected peptide bond across the mature cleavage: {cleavage_bonds}"
        )
    return {
        "original_heavy_atom_count": len(source_atoms),
        "repaired_heavy_atom_count": len(output_atoms),
        "maximum_original_heavy_displacement_angstrom": max(
            (math.dist(xyz, output_atoms[key]) for key, xyz in source_atoms.items()),
            default=0.0,
        ),
        "residues_per_chain": {
            chain: len(values) for chain, values in residues.items()
        },
        "n_terminal_group_count": sum(
            "N" in atoms_by_residue[(chain, lo)]
            for chain, (lo, _) in EXPECTED_RANGES.items()
        ),
        "c_terminal_group_count": sum(
            "C" in atoms_by_residue[(chain, hi)]
            for chain, (_, hi) in EXPECTED_RANGES.items()
        ),
        "oxt_count": sum(
            "OXT" in atoms_by_residue[(chain, hi)]
            for chain, (_, hi) in EXPECTED_RANGES.items()
        ),
        "cleavage_interchain_cn_bond_count": len(cleavage_bonds),
        "source_and_repaired_structures_heavy_atom_only": True,
    }


def run_repair(source: Path, output_dir: Path) -> dict[str, object]:
    """Run the explicit 5XYG mature-chain repair and write scientific acceptance."""
    from openmm import Platform
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    output_dir.mkdir(parents=True, exist_ok=False)
    split_path = output_dir / "5XYG_split_protein_input.pdb"
    repaired_path = output_dir / "5XYG_repaired_protein.pdb"
    prep = write_split_input(source, split_path)
    original_atoms = read_atoms(split_path)
    fixer = PDBFixer(filename=str(split_path))
    fixer.platform = Platform.getPlatformByName("CPU")
    fixer.findMissingResidues()
    if len(list(fixer.topology.chains())) != 4:
        raise ValueError("PDBFixer input did not produce four explicitly split chains")
    missing_residue_groups = {str(k): list(v) for k, v in fixer.missingResidues.items()}
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        raise ValueError(
            f"Unexpected nonstandard residues; refusing substitutions: {fixer.nonstandardResidues}"
        )
    fixer.findMissingAtoms()
    missing_atom_count = sum(len(value) for value in fixer.missingAtoms.values())
    missing_terminal_count = sum(
        len(value) for value in fixer.missingTerminals.values()
    )
    fixer.addMissingAtoms(seed=20261120)
    with repaired_path.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    repaired_file = PDBFile(str(repaired_path))
    if {chain.id for chain in repaired_file.topology.chains()} != set("ABCD"):
        raise ValueError("PDB write/read changed the required four-chain identity")
    acceptance = validate_repair(
        original_atoms,
        repaired_path,
        prep["expected_sequences"],
        repaired_file.topology,
    )
    report = {
        "state": "accepted",
        "source": str(source),
        "split_input": str(split_path),
        "output": str(repaired_path),
        "seed": 20261120,
        "platform": fixer.platform.getName(),
        "missing_residue_groups": missing_residue_groups,
        "missing_heavy_atom_count": missing_atom_count,
        "missing_terminal_atom_count": missing_terminal_count,
        **prep,
        **acceptance,
    }
    (output_dir / "scientific_acceptance.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    LOGGER.info("5XYG residue repair passed all declared acceptance checks")
    return report


def main() -> None:
    """Parse arguments and run the 5XYG residue repair."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    report = run_repair(args.source, args.output_dir)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
