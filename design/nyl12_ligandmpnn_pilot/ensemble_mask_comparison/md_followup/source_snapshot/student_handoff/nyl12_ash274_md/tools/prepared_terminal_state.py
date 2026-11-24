"""Explicit four-copy terminal protonation contracts; never infer state from atoms."""

from __future__ import annotations

import math
from typing import Any

LEGACY_NEUTRAL_STATES = frozenset({"ASP272/ASP274", "ASP272/ASH274", "ASH272/ASP274"})


def terminal_contract(dataset: dict[str, Any]) -> tuple[str, float, str]:
    """Resolve an explicit terminal declaration or the original neutral schema.

    Args:
        dataset: Declared fixed Asp microstate and optional terminal_protonation.

    Returns:
        NH2/NH3+, expected whole-Thr charge and the declaration contract name.

    Raises:
        ValueError: If the declaration or authorized paired microstate is invalid.
    """
    if "terminal_protonation" not in dataset:
        if dataset["microstate"] not in LEGACY_NEUTRAL_STATES:
            raise ValueError(
                "Missing explicit terminal protonation or recognized legacy-neutral contract"
            )
        return "NH2", -0.001, "legacy_neutral_state_contract_v1"
    state = dataset["terminal_protonation"]
    if state not in ("NH2", "NH3+"):
        raise ValueError("terminal_protonation must explicitly be NH2 or NH3+")
    if state == "NH3+" and dataset["microstate"] != "ASP272/ASP274":
        raise ValueError("Charged Thr NH3+ is authorized only with ASP272/ASP274")
    return state, 1.0 if state == "NH3+" else 0.0, "explicit_terminal_protonation_v2"


def validate_terminal(residue: Any, state: str, expected_charge: float) -> float:
    """Require declared N-bound H inventory and the WHOLE Thr-residue charge.

    Args:
        residue: ParmEd catalytic Thr residue with bonded atoms and charges.
        state: Explicitly resolved NH2 or NH3+ declaration.
        expected_charge: Whole-residue expected sum from the declared recipe.

    Returns:
        Validated whole-residue partial-charge sum in elementary-charge units.

    Raises:
        ValueError: If hydrogen identity/bonding or whole-residue charge disagrees.
    """
    if state not in ("NH2", "NH3+"):
        raise ValueError("Unknown declared terminal state")
    atoms = {atom.name: atom for atom in residue.atoms}
    expected = {"H1", "H2", "H3"} if state == "NH3+" else {"H1", "H2"}
    if (
        "N" not in atoms
        or {name for name in atoms if name in {"H1", "H2", "H3"}} != expected
    ):
        raise ValueError(f"Declared terminal {state} hydrogen inventory mismatch")
    bonded = {atom.name for atom in atoms["N"].bond_partners if atom.atomic_number == 1}
    if bonded != expected:
        raise ValueError(f"Declared terminal {state} N-H connectivity mismatch")
    charge = sum(float(atom.charge) for atom in residue.atoms)
    if not math.isfinite(charge) or not math.isclose(
        charge, expected_charge, abs_tol=1e-5
    ):
        raise ValueError(
            f"Whole Thr-residue charge {charge} differs from declared {expected_charge}"
        )
    return charge
