"""Explicit, lossless-identity adapter for centered BF16 native RFD3 output."""

from __future__ import annotations

from typing import Any

import numpy as np

AtomKey = tuple[str, int, str]


def bf16_coordinates(values: np.ndarray) -> np.ndarray:
    """Round float32 coordinates to BF16 using round-to-nearest-even.

    Args:
        values: Finite coordinate components.

    Returns:
        Rounded values represented as float32.
    """
    array = np.asarray(values, dtype=np.float32)
    if not np.isfinite(array).all():
        raise ValueError("Nonfinite source coordinates")
    bits = array.view(np.uint32)
    rounded = (bits + np.uint32(0x7FFF) + ((bits >> 16) & 1)) & np.uint32(0xFFFF0000)
    return rounded.view(np.float32)


def adapt_coordinates(
    source_keys: list[AtomKey],
    source: np.ndarray,
    native_keys: list[AtomKey],
    native: np.ndarray,
    mobile_residues: set[str],
    center: np.ndarray,
    signature_tolerance: float = 1e-7,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Restore fixed source components only after proving their BF16 signature.

    Args:
        source_keys: Complete immutable source atom identities.
        source: High precision coordinates in source order.
        native_keys: Bijective native identities, excluding exactly eight OXT.
        native: Native coordinates, never modified in place.
        mobile_residues: Explicit residues diffused by the native input.
        center: Exact float32 COM used by the native partial-diffusion parser.
        signature_tolerance: Serialization tolerance for exact BF16 components.

    Returns:
        Source-frame derived coordinates and explicit restoration evidence.

    Raises:
        ValueError: Inventory, coordinates or fixed BF16 signature differ.
    """
    if len(set(source_keys)) != len(source_keys) or len(set(native_keys)) != len(
        native_keys
    ):
        raise ValueError("Nonbijective source/native atom identity")
    missing = {(c, 232 if c in "ABCD" else 96, "OXT") for c in "ABCDEFGH"}
    if set(source_keys) - set(native_keys) != missing or set(native_keys) - set(
        source_keys
    ):
        raise ValueError("Unexpected atom inventory; require exactly eight fixed OXT")
    if source.shape != (len(source_keys), 3) or native.shape != (len(native_keys), 3):
        raise ValueError("Coordinate shape differs from identity inventory")
    if not np.isfinite(native).all() or not np.isfinite(source).all():
        raise ValueError("Nonfinite coordinates")
    lookup = {key: i for i, key in enumerate(source_keys)}
    order = np.array([lookup[key] for key in native_keys])
    fixed = np.array([f"{c}{p}" not in mobile_residues for c, p, _ in native_keys])
    if fixed.sum() < 3 or fixed.all():
        raise ValueError("Invalid fixed/mobile mask")
    expected = bf16_coordinates(np.asarray(source, dtype=np.float32)[order] - center)
    error = float(np.max(np.abs(expected[fixed].astype(float) - native[fixed])))
    if error > signature_tolerance:
        raise ValueError(f"Fixed coordinates lack exact BF16 source signature: {error}")
    for chain, position, _ in missing:
        if f"{chain}{position}" in mobile_residues:
            raise ValueError("Missing OXT is mobile")
        if any((chain, position, name) not in native_keys for name in ("C", "O")):
            raise ValueError("Missing fixed terminal anchors")
    result = np.array(source, dtype=float, copy=True)
    # Translation applies equally to all native atoms; mobile components are never repaired.
    result[order[~fixed]] = native[~fixed].astype(float) + center.astype(float)
    return result, {
        "accepted": True,
        "center_A": center.tolist(),
        "fixed_components_compared": int(fixed.sum() * 3),
        "max_BF16_signature_error_A": error,
        "signature_serialization_tolerance_A": signature_tolerance,
        "mobile_coordinates_preserved_exactly_after_translation": bool(
            np.array_equal(result[order[~fixed]], native[~fixed].astype(float) + center)
        ),
        "restored_fixed_atoms": int(fixed.sum()),
        "restored_source_OXT": [list(key) for key in sorted(missing)],
        "raw_native_modified": False,
    }


def reassociate_typed_graph(
    source_keys: list[AtomKey],
    source_bonds: np.ndarray,
    target_keys: list[AtomKey],
) -> np.ndarray:
    """Map every explicit source typed bond through an exact atom bijection.

    Args:
        source_keys: Immutable graph identity order.
        source_bonds: Explicit source bond index pairs and bond types.
        target_keys: Complete derived output identity order.

    Returns:
        All typed bonds in target index order, without spatial inference.

    Raises:
        ValueError: Atom correspondence or typed source graph is invalid.
    """
    if len(set(source_keys)) != len(source_keys) or len(set(target_keys)) != len(
        target_keys
    ):
        raise ValueError("Nonbijective graph identities")
    if set(source_keys) != set(target_keys):
        raise ValueError("Graph atom identities differ")
    bonds = np.asarray(source_bonds, dtype=int)
    if bonds.ndim != 2 or bonds.shape[1] != 3 or not len(bonds):
        raise ValueError("Missing explicit typed source graph")
    if (
        np.any(bonds[:, :2] < 0)
        or np.any(bonds[:, :2] >= len(source_keys))
        or np.any(bonds[:, 2] <= 0)
    ):
        raise ValueError("Invalid typed source bond")
    lookup = {key: i for i, key in enumerate(target_keys)}
    result = bonds.copy()
    result[:, :2] = [
        [lookup[source_keys[a]], lookup[source_keys[b]]] for a, b in bonds[:, :2]
    ]
    return result
