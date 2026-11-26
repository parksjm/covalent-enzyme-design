"""Validated, display-only chain selections for stripped Nyl12 movies."""

from collections.abc import Mapping
from typing import Any


_NORMALIZE = {
    "ASH": "ASP",
    "GLH": "GLU",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "CYX": "CYS",
    "CYM": "CYS",
    "LYN": "LYS",
}
_PROTEIN = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
)
_BOUNDS = (0, 232, 464, 696, 928, 1024, 1120, 1216, 1312)


def display_chain_groups(
    topology: Any, mapping: Mapping[int, int]
) -> dict[str, list[int]]:
    """Return physical display indices grouped by explicit or validated chains.

    Args:
        topology: ParmEd topology with residues and native bond partners.
        mapping: Original one-based atom IDs to displayed zero-based indices.

    Returns:
        Chain labels mapped to retained atom indices. Existing complete native
        labels are preserved. Blank Nyl12 labels are reconstructed only after
        verifying four 232-residue alpha and four 96-residue beta sequences,
        their ordering, Thr beta starts and all seven peptide-chain breaks.

    Raises:
        ValueError: If labels are partial or the unlabeled topology does not
            support the documented mature Nyl12 ordering.

    Notes:
        Neither native atom/residue attributes nor topology files are changed.
        A/E, B/F, C/G and D/H are the documented alpha/beta display pairs.
    """
    residues = []
    residue_atoms: dict[int, list[Any]] = {}
    for atom in topology.atoms:
        residue = atom.residue
        if _NORMALIZE.get(residue.name, residue.name) not in _PROTEIN:
            continue
        identity = id(residue)
        if identity not in residue_atoms:
            residues.append(residue)
            residue_atoms[identity] = []
        residue_atoms[identity].append(atom)
    if not residues:
        raise ValueError("No protein residues available for chain display")
    labels = [str(residue.chain).strip() for residue in residues]
    if any(labels):
        if not all(labels):
            raise ValueError("Partial native chain labels require an explicit repair")
        assignments = list(zip(residues, labels, strict=True))
    else:
        if len(residues) != 1312 or [r.idx for r in residues] != list(range(1312)):
            raise ValueError(
                "Blank chains require the documented 1312-residue mature Nyl12 topology"
            )
        blocks = [residues[a:b] for a, b in zip(_BOUNDS, _BOUNDS[1:], strict=False)]
        sequences = [
            tuple(_NORMALIZE.get(r.name, r.name) for r in block) for block in blocks
        ]
        if len(set(sequences[:4])) != 1 or len(set(sequences[4:])) != 1:
            raise ValueError(
                "Mature alpha/beta copies have inconsistent residue sequences"
            )
        if any(block[0].name != "THR" for block in blocks[4:]):
            raise ValueError(
                "Mature beta chains must begin at the four expected Thr residues"
            )
        for boundary in _BOUNDS[1:-1]:
            left, right = residues[boundary - 1], residues[boundary]
            carbon = [atom for atom in residue_atoms[id(left)] if atom.name == "C"]
            nitrogen = [atom for atom in residue_atoms[id(right)] if atom.name == "N"]
            if len(carbon) != 1 or len(nitrogen) != 1:
                raise ValueError(
                    "A chain boundary lacks an unambiguous peptide C/N pair"
                )
            if nitrogen[0] in carbon[0].bond_partners:
                raise ValueError(
                    "Expected mature peptide-chain break is covalently connected"
                )
        assignments = [
            (residue, label)
            for label, block in zip("ABCDEFGH", blocks, strict=True)
            for residue in block
        ]
    groups: dict[str, list[int]] = {}
    for residue, label in assignments:
        kept = [
            int(mapping[a.idx + 1])
            for a in residue_atoms[id(residue)]
            if a.idx + 1 in mapping
        ]
        if kept:
            groups.setdefault(label, []).extend(kept)
    if not groups:
        raise ValueError("No mapped protein atoms available for chain cartoons")
    return {label: sorted(indices) for label, indices in groups.items()}
