"""Render fixed-geometry Gaussian ESP inputs for accepted PA6 conformers."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem

from nyl12_covalent.pa6_resp_preparation import CORES, MEMORY_GB


def hf_esp_input(elements: list[int], xyz: NDArray[np.float64], label: str) -> str:
    """Return a nonoptimizing HF/6-31G(d) Merz--Kollman Gaussian input.

    Args:
        elements: Atomic numbers in canonical PA6 atom order.
        xyz: Accepted coordinates in Angstrom and the same atom order.
        label: Provenance-safe conformer label used for the checkpoint name.

    Returns:
        Complete Gaussian input text for a neutral singlet fixed geometry.
    """
    if xyz.shape != (len(elements), 3) or not np.isfinite(xyz).all():
        raise ValueError("Invalid accepted coordinates for ESP input")
    periodic = Chem.GetPeriodicTable()
    lines = [
        f"%nprocshared={CORES}",
        f"%mem={MEMORY_GB}GB",
        f"%chk={label}_hf_mk.chk",
        "#p HF/6-31G(d) Pop=MK IOp(6/33=2,6/42=6) SCF=(Tight,XQC) NoSymm",
        "",
        f"PA6 L4 zwitterion {label}; fixed accepted geometry MK ESP",
        "",
        "0 1",
    ]
    lines.extend(
        f"{periodic.GetElementSymbol(number)} {x:.10f} {y:.10f} {z:.10f}"
        for number, (x, y, z) in zip(elements, xyz, strict=True)
    )
    return "\n".join([*lines, "", ""])
