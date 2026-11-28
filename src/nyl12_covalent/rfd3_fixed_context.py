"""Restore only eight immutable terminal oxygens omitted by native RFD3."""

from __future__ import annotations

import numpy as np

from .rfd3_qualification_audit import align_on_fixed

AtomKey = tuple[str, int, str]


def reconcile_terminal_oxygens(
    source_keys: list[AtomKey],
    source: np.ndarray,
    generated_keys: list[AtomKey],
    generated: np.ndarray,
    mobile_residues: set[str],
    anchor_tolerance: float = 0.05,
) -> tuple[np.ndarray, dict]:
    """Align native output and restore exactly the declared eight fixed OXT atoms.

    Args:
        source_keys: Exact complete source atom identities.
        source: Source coordinates in the same order.
        generated_keys: Native generated identities after reversible label mapping.
        generated: Raw generated coordinates in matching order.
        mobile_residues: Explicit residues allowed to move.
        anchor_tolerance: Maximum terminal C/O displacement after fixed-atom alignment.

    Returns:
        Complete aligned coordinates and explicit omission/restoration evidence.

    Raises:
        ValueError: Omissions, anchors or correspondence violate the approved rule.
    """
    expected_missing = {(c, 232 if c in "ABCD" else 96, "OXT") for c in "ABCDEFGH"}
    if len(set(generated_keys)) != len(generated_keys) or len(set(source_keys)) != len(
        source_keys
    ):
        raise ValueError("Ambiguous native/source atom correspondence")
    if set(source_keys) - set(generated_keys) != expected_missing or set(
        generated_keys
    ) - set(source_keys):
        raise ValueError(
            "Native omissions differ from exactly eight source terminal OXT atoms"
        )
    original = {key: i for i, key in enumerate(source_keys)}
    fixed = np.array([f"{c}{p}" not in mobile_residues for c, p, _ in generated_keys])
    target = source[[original[k] for k in generated_keys]]
    aligned = align_on_fixed(target, generated, fixed)
    x, y = generated[fixed], target[fixed]
    cx, cy = x.mean(axis=0), y.mean(axis=0)
    u, _, vt = np.linalg.svd((x - cx).T @ (y - cy))
    parity = np.eye(3)
    parity[2, 2] = np.linalg.det(u @ vt)
    rotation = u @ parity @ vt
    translation = cy - cx @ rotation
    result = np.array(source, copy=True)
    observed = {key: i for i, key in enumerate(generated_keys)}
    for key, index in observed.items():
        result[original[key]] = aligned[index]
    restored = []
    for chain, position, atom in sorted(expected_missing):
        if f"{chain}{position}" in mobile_residues:
            raise ValueError("OXT belongs to a mobile residue")
        anchors = []
        for name in ("C", "O"):
            key = (chain, position, name)
            if key not in observed:
                raise ValueError("Required retained terminal C/O anchor is absent")
            displacement = float(
                np.linalg.norm(result[original[key]] - source[original[key]])
            )
            if displacement > anchor_tolerance:
                raise ValueError("Immutable terminal C/O anchor moved beyond tolerance")
            anchors.append({"atom": name, "displacement_A": displacement})
        restored.append(
            {
                "chain": chain,
                "position": position,
                "atom": atom,
                "source_index0based": original[(chain, position, atom)],
                "source_xyz": source[original[(chain, position, atom)]].tolist(),
                "fixed_C_O_anchors": anchors,
            }
        )
    return result, {
        "accepted": True,
        "raw_native_atom_count": len(generated_keys),
        "complete_derived_atom_count": len(source_keys),
        "restored_source_atoms": restored,
        "raw_native_coordinates_modified": False,
        "raw_to_source_transform": {
            "rotation": rotation.tolist(),
            "translation": translation.tolist(),
        },
    }
