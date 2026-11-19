"""Scientific gate for a preorganized Thr233-to-Asp272 proton relay."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm

THR_OG1 = 13916
THR_HG1 = 13917
ASP272_OD1 = 14559


def relay_distances(topology: Path, coordinates: Path) -> dict[str, float]:
    """Measure the reviewed Thr233–Asp272 OD1 relay contacts.

    Args:
        topology: Amber topology for the solvated modular-PA6 complex.
        coordinates: Amber coordinate or restart file.

    Returns:
        OG1–OD1 and HG1–OD1 distances in angstrom.
    """
    system = AmberParm(str(topology), xyz=str(coordinates))
    xyz = system.coordinates
    return {
        "OG1_Asp272_OD1_A": float(
            np.linalg.norm(xyz[THR_OG1 - 1] - xyz[ASP272_OD1 - 1])
        ),
        "HG1_Asp272_OD1_A": float(
            np.linalg.norm(xyz[THR_HG1 - 1] - xyz[ASP272_OD1 - 1])
        ),
    }


def relay_contact_gate(distances: dict[str, float]) -> dict[str, Any]:
    """Apply the predeclared relay-ready contact bounds.

    Args:
        distances: Output from :func:`relay_distances`.

    Returns:
        Acceptance flag, measured distances, and immutable bounds.
    """
    og1 = distances["OG1_Asp272_OD1_A"]
    hg1 = distances["HG1_Asp272_OD1_A"]
    bounds = {
        "OG1_Asp272_OD1_A": [2.5, 3.2],
        "HG1_Asp272_OD1_A": [1.5, 2.3],
    }
    return {
        "accepted": bounds["OG1_Asp272_OD1_A"][0]
        <= og1
        <= bounds["OG1_Asp272_OD1_A"][1]
        and bounds["HG1_Asp272_OD1_A"][0] <= hg1 <= bounds["HG1_Asp272_OD1_A"][1],
        "distances_A": distances,
        "bounds_A": bounds,
        "interpretation": (
            "This gate tests relay preorganization, not native Thr233 similarity."
        ),
    }
