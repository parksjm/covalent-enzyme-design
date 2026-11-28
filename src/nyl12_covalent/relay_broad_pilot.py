"""Preparation rules for the broad, symmetry-tied relay ensemble pilot."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .ligandmpnn_pilot_audit import THREE_TO_ONE
from .nyl12_sequence_protection import (
    MINIMUM_PROTECTED_ALPHA,
    MINIMUM_PROTECTED_BETA_LOCAL,
)

ALPHA_CHAINS = "ABCD"
BETA_CHAINS = "EFGH"
PROTEIN_CHAINS = ALPHA_CHAINS + BETA_CHAINS
REPEAT_PAIRS = frozenset({frozenset(pair) for pair in ("AE", "BF", "CG", "DH")})
EXPECTED_LENGTHS = {
    **dict.fromkeys(ALPHA_CHAINS, 232),
    **dict.fromkeys(BETA_CHAINS, 96),
}
CORE_PROTECTED = {
    "alpha": {1, *MINIMUM_PROTECTED_ALPHA},
    "beta": set(MINIMUM_PROTECTED_BETA_LOCAL),
}
INTERFACE_DISTANCE_ANGSTROM = 4.0
INTERFACE_COPY_MINIMUM = 3


@dataclass(frozen=True)
class Atom:
    """One heavy protein atom with canonical chain-local residue mapping."""

    chain: str
    position: int
    atom_name: str
    residue_name: str
    element: str
    xyz: tuple[float, float, float]


@dataclass(frozen=True)
class ParentStructure:
    """Validated protein atoms and per-chain sequences from one parent PDB."""

    path: Path
    atoms: tuple[Atom, ...]
    sequences: dict[str, str]


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file's exact bytes."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _local_position(chain: str, residue_number: int) -> int:
    """Map source precursor numbering or LigandMPNN numbering to local chains."""
    if chain in ALPHA_CHAINS:
        position = residue_number
    elif chain in BETA_CHAINS:
        position = residue_number - 232 if residue_number > 232 else residue_number
    else:
        raise ValueError(f"Unexpected protein chain {chain!r}")
    if not 1 <= position <= EXPECTED_LENGTHS[chain]:
        raise ValueError(f"Residue outside canonical chain: {chain}:{residue_number}")
    return position


def read_parent(path: Path) -> ParentStructure:
    """Read a complete eight-chain tetramer and validate canonical numbering.

    Args:
        path: Full heavy-atom PDB with protein chains A-H. Nonprotein chains
            are ignored for sequence/interface census and retained separately
            by the pipeline manifests.

    Returns:
        Validated protein atom coordinates and identical alpha/beta copy
        sequences.

    Raises:
        ValueError: If chain mapping, residue identity, or heavy coordinates
            are incomplete or inconsistent.
    """
    atoms: list[Atom] = []
    residues: dict[str, dict[int, str]] = {chain: {} for chain in PROTEIN_CHAINS}
    cysteine_sg: set[tuple[str, int]] = set()
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.startswith("ATOM  ") or line[21] not in PROTEIN_CHAINS:
            continue
        chain = line[21]
        position = _local_position(chain, int(line[22:26]))
        residue_name = line[17:20].strip()
        if residue_name not in THREE_TO_ONE:
            raise ValueError(f"Nonstandard residue at {path}:{line_number}")
        residue = residues[chain].get(position)
        if residue is not None and residue != residue_name:
            raise ValueError(f"Conflicting residue identities at {chain}:{position}")
        residues[chain][position] = residue_name
        atom_name = line[12:16].strip()
        element = line[76:78].strip().upper()
        if element in {"H", "D"}:
            continue
        xyz = tuple(float(line[start : start + 8]) for start in (30, 38, 46))
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"Nonfinite coordinates at {path}:{line_number}")
        atoms.append(Atom(chain, position, atom_name, residue_name, element, xyz))
        if residue_name == "CYS" and atom_name == "SG":
            cysteine_sg.add((chain, position))

    sequences: dict[str, str] = {}
    for chain, length in EXPECTED_LENGTHS.items():
        if sorted(residues[chain]) != list(range(1, length + 1)):
            raise ValueError(f"Incomplete canonical chain {chain} in {path}")
        sequences[chain] = "".join(
            THREE_TO_ONE[residues[chain][position]] for position in range(1, length + 1)
        )
    if len({sequences[chain] for chain in ALPHA_CHAINS}) != 1:
        raise ValueError(f"Alpha sequences differ across copies in {path}")
    if len({sequences[chain] for chain in BETA_CHAINS}) != 1:
        raise ValueError(f"Beta sequences differ across copies in {path}")
    for chain in PROTEIN_CHAINS:
        for position, name in residues[chain].items():
            if name == "CYS" and (chain, position) not in cysteine_sg:
                raise ValueError(f"Cysteine lacks SG atom at {chain}:{position}")
    return ParentStructure(path.resolve(), tuple(atoms), sequences)


def _cell(xyz: tuple[float, float, float], width: float) -> tuple[int, int, int]:
    return tuple(math.floor(value / width) for value in xyz)


def inter_repeat_contacts(
    parent: ParentStructure, cutoff: float = INTERFACE_DISTANCE_ANGSTROM
) -> dict[str, set[int]]:
    """Find protein positions contacting a different alpha/beta repeat.

    Contacts within AE, BF, CG, or DH are same-repeat contacts and are omitted.
    A spatial hash keeps this audit bounded without a new spatial-library
    dependency.
    """
    if cutoff <= 0 or not math.isfinite(cutoff):
        raise ValueError("Contact cutoff must be positive and finite")
    by_chain: dict[str, list[Atom]] = {chain: [] for chain in PROTEIN_CHAINS}
    for atom in parent.atoms:
        by_chain[atom.chain].append(atom)
    contacts = {chain: set() for chain in PROTEIN_CHAINS}
    width = cutoff
    cutoff_squared = cutoff * cutoff
    neighbors = range(-1, 2)
    for chain_index, left_chain in enumerate(PROTEIN_CHAINS):
        for right_chain in PROTEIN_CHAINS[chain_index + 1 :]:
            if frozenset((left_chain, right_chain)) in REPEAT_PAIRS:
                continue
            right_grid: dict[tuple[int, int, int], list[Atom]] = {}
            for atom in by_chain[right_chain]:
                right_grid.setdefault(_cell(atom.xyz, width), []).append(atom)
            for left_atom in by_chain[left_chain]:
                cx, cy, cz = _cell(left_atom.xyz, width)
                for dx in neighbors:
                    for dy in neighbors:
                        for dz in neighbors:
                            for right_atom in right_grid.get(
                                (cx + dx, cy + dy, cz + dz), ()
                            ):
                                distance_squared = sum(
                                    (left_atom.xyz[index] - right_atom.xyz[index]) ** 2
                                    for index in range(3)
                                )
                                if distance_squared <= cutoff_squared:
                                    contacts[left_chain].add(left_atom.position)
                                    contacts[right_chain].add(right_atom.position)
    return contacts


def _disulfide_positions(parents: list[ParentStructure]) -> dict[str, set[int]]:
    """Return cysteine residues forming an SG-SG bond in any supplied parent."""
    protected = {"alpha": set(), "beta": set()}
    for parent in parents:
        sg_atoms = [
            atom
            for atom in parent.atoms
            if atom.residue_name == "CYS" and atom.atom_name == "SG"
        ]
        for index, left in enumerate(sg_atoms):
            for right in sg_atoms[index + 1 :]:
                distance_squared = sum(
                    (left.xyz[axis] - right.xyz[axis]) ** 2 for axis in range(3)
                )
                if distance_squared <= 2.3**2:
                    for atom in (left, right):
                        family = "alpha" if atom.chain in ALPHA_CHAINS else "beta"
                        protected[family].add(atom.position)
    return protected


def protection_census(
    parent_paths: list[Path],
    *,
    interface_cutoff: float = INTERFACE_DISTANCE_ANGSTROM,
    interface_copy_minimum: int = INTERFACE_COPY_MINIMUM,
) -> dict[str, Any]:
    """Build the canonical hard-protection and mutable-position receipt.

    Interface positions are hard-protected only if each of three independent
    original TS parents shows contact in at least ``interface_copy_minimum`` of
    that sequence family’s four equivalent copies. Other observed interface
    contacts are retained as mutable interface-periphery diagnostics.
    """
    if len(parent_paths) != 3:
        raise ValueError("Interface census requires exactly three original TS parents")
    if not 1 <= interface_copy_minimum <= 4:
        raise ValueError("Interface-copy threshold must be between one and four")
    parents = [read_parent(path) for path in parent_paths]
    reference_sequences = parents[0].sequences
    if any(parent.sequences != reference_sequences for parent in parents[1:]):
        raise ValueError("Original parent sequences differ")

    observed: dict[str, set[int]] = {"alpha": set(), "beta": set()}
    recurring_by_parent: list[dict[str, set[int]]] = []
    recurring: dict[str, set[int]] = {"alpha": set(), "beta": set()}
    per_parent: list[dict[str, Any]] = []
    for parent in parents:
        contacts = inter_repeat_contacts(parent, interface_cutoff)
        family_data: dict[str, Any] = {}
        family_recurring: dict[str, set[int]] = {}
        for family, chains in (("alpha", ALPHA_CHAINS), ("beta", BETA_CHAINS)):
            length = 232 if family == "alpha" else 96
            observed_family = {
                position for chain in chains for position in contacts[chain]
            }
            recurring_family = {
                position
                for position in range(1, length + 1)
                if sum(position in contacts[chain] for chain in chains)
                >= interface_copy_minimum
            }
            observed[family].update(observed_family)
            family_recurring[family] = recurring_family
            family_data[family] = {
                "observed_any_copy": sorted(observed_family),
                "recurring_copies": sorted(recurring_family),
                "copy_contact_positions": {
                    chain: sorted(contacts[chain]) for chain in chains
                },
            }
        per_parent.append(
            {
                "path": str(parent.path),
                "sha256": sha256_file(parent.path),
                "interface": family_data,
            }
        )
        recurring_by_parent.append(family_recurring)

    for family in ("alpha", "beta"):
        recurring[family] = set.intersection(
            *(item[family] for item in recurring_by_parent)
        )

    disulfides = _disulfide_positions(parents)
    protected = {
        family: set(CORE_PROTECTED[family]) | recurring[family] | disulfides[family]
        for family in ("alpha", "beta")
    }
    observed_mutable = {
        family: observed[family] - recurring[family] for family in ("alpha", "beta")
    }
    mutable = {
        "alpha": sorted(set(range(1, 233)) - protected["alpha"]),
        "beta": sorted(set(range(1, 97)) - protected["beta"]),
    }
    precursor_protected = sorted(
        protected["alpha"] | {position + 232 for position in protected["beta"]}
    )
    if len(precursor_protected) + len(mutable["alpha"]) + len(mutable["beta"]) != 328:
        raise AssertionError(
            "Canonical mask does not partition all 328 precursor positions"
        )
    return {
        "schema_version": 1,
        "status": "prepared",
        "canonical_numbering": {
            "alpha": "precursor 1-232; copies A-D",
            "beta": "precursor 233-328; copies E-H; local positions 1-96",
        },
        "source_parents": per_parent,
        "protection_policy": {
            "hard_core": {
                "alpha": sorted(CORE_PROTECTED["alpha"]),
                "beta_local": sorted(CORE_PROTECTED["beta"]),
                "note": "Includes initiator M1, catalytic/maturation positions, and D272-G273-D274 loop.",
            },
            "disulfides": {
                "alpha": sorted(disulfides["alpha"]),
                "beta_local": sorted(disulfides["beta"]),
            },
            "recurring_inter_repeat_interface": {
                "cutoff_angstrom": interface_cutoff,
                "required_copies": interface_copy_minimum,
                "all_three_original_parents_required": True,
                "same_repeat_pairs_excluded": ["AE", "BF", "CG", "DH"],
                "alpha": sorted(recurring["alpha"]),
                "beta_local": sorted(recurring["beta"]),
                "not_experimentally_proven_essential": True,
            },
            "mutable_interface_periphery": {
                "alpha": sorted(observed_mutable["alpha"] - protected["alpha"]),
                "beta_local": sorted(observed_mutable["beta"] - protected["beta"]),
                "interpretation": "Observed contacts that fail the recurring-interface protection threshold remain mutable and require downstream assembly scrutiny.",
            },
        },
        "fixed": {
            "alpha": sorted(protected["alpha"]),
            "beta_local": sorted(protected["beta"]),
            "precursor_positions": precursor_protected,
        },
        "mutable": mutable,
        "counts": {
            "precursor_positions": 328,
            "fixed_alpha": len(protected["alpha"]),
            "fixed_beta": len(protected["beta"]),
            "fixed_total": len(precursor_protected),
            "mutable_alpha": len(mutable["alpha"]),
            "mutable_beta": len(mutable["beta"]),
            "mutable_total": len(mutable["alpha"]) + len(mutable["beta"]),
        },
        "source_sequences": {
            "alpha": reference_sequences["A"],
            "beta": reference_sequences["E"],
        },
    }


def symmetry_manifest(mask: dict[str, list[int]]) -> dict[str, str]:
    """Return explicit A-D/E-H tied mutation arguments for LigandMPNN."""
    groups = [
        ",".join(f"{chain}{position}" for chain in ALPHA_CHAINS)
        for position in mask["alpha"]
    ] + [
        ",".join(f"{chain}{position}" for chain in BETA_CHAINS)
        for position in mask["beta"]
    ]
    designed = [
        f"{chain}{position}"
        for family_chains, family_positions in (
            (ALPHA_CHAINS, mask["alpha"]),
            (BETA_CHAINS, mask["beta"]),
        )
        for position in family_positions
        for chain in family_chains
    ]
    if not groups:
        raise ValueError("Broad mask is empty")
    return {
        "symmetry_residues": "|".join(groups),
        "symmetry_weights": "|".join("1,0,0,0" for _ in groups),
        "redesigned_residues": " ".join(designed),
    }


def sample_schedule(
    original_parent_ids: list[str], remodeled_parent_ids: list[str]
) -> list[dict[str, Any]]:
    """Build the fixed 6+6+3-per-remodelled-parent sample allocation."""
    if len(original_parent_ids) != 2 or len(set(original_parent_ids)) != 2:
        raise ValueError("Exactly two distinct original parent IDs must be selected")
    if not 0 <= len(remodeled_parent_ids) <= 4:
        raise ValueError("Zero to four accepted remodeled parents are permitted")
    all_ids = original_parent_ids + remodeled_parent_ids
    if len(set(all_ids)) != len(all_ids):
        raise ValueError("Parent IDs must be globally unique")
    schedule = [
        {"parent_id": parent_id, "parent_kind": "original_TS", "raw_samples": 6}
        for parent_id in original_parent_ids
    ]
    schedule.extend(
        {
            "parent_id": parent_id,
            "parent_kind": "remodeled_TS",
            "raw_samples": 3,
        }
        for parent_id in remodeled_parent_ids
    )
    raw_count = sum(item["raw_samples"] for item in schedule)
    if raw_count > 24:
        raise AssertionError("Generation schedule exceeds the 24-sample ceiling")
    return schedule


def sha256_json_file(path: Path) -> tuple[str, dict[str, Any]]:
    """Read and hash a JSON file without normalizing its bytes."""
    data = json.loads(path.read_text())
    return sha256_file(path), data
