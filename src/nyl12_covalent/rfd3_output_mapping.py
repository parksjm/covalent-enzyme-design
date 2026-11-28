"""Strict native RFD3 output correspondence, with exported ligand residues1–4."""

from __future__ import annotations

from collections.abc import Sequence


def canonical_output_keys(
    chains: Sequence[str], residues: Sequence[int], names: Sequence[str]
) -> list[tuple[str, int, str]]:
    """Preserve verified output numbering and reject ambiguous/offset identities."""
    keys = [
        (str(c), int(p), str(n))
        for c, p, n in zip(chains, residues, names, strict=True)
    ]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate native atom identities")
    if {p for c, p, _ in keys if c == "I"} != {1, 2, 3, 4}:
        raise ValueError(
            "Native exported ligand must use exactly I1-I4; no numbering fallback"
        )
    for chain, position, _ in keys:
        if chain not in "ABCDEFGHI":
            raise ValueError("Unexpected native chain")
        if chain != "I" and not 1 <= position <= (232 if chain in "ABCD" else 96):
            raise ValueError("Unexpected native protein residue number")
    return keys
