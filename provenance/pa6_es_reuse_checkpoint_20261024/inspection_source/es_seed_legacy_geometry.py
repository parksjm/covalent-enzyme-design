"""Validate retained historical PDB phi separately from primary NC geometry."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from nyl12_covalent.active_site_geometry import projected_herschlag_phi


def inspect_legacy_phi(parm: Any, pdb: Path) -> dict:
    """Recompute historical phi and verify PDB quantization against the NC seed.

    Args:
        parm: Original, imaged Amber seed before display alignment.
        pdb: Actual retained five-atom historical CPPTRAJ PDB for this frame.

    Returns:
        Historical phi and explicit per-component coordinate differences.

    Raises:
        ValueError: Atom identity, finiteness or precision bound is invalid.
        FileNotFoundError: The required historical PDB is unavailable.
    """
    requirements = [
        (929, "THR", "OG1"),
        (1313, "LIG", "O49"),
        (1313, "LIG", "C54"),
        (1313, "LIG", "N59"),
        (1313, "LIG", "C62"),
    ]
    pdb_atoms = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = (int(line[22:26]), line[17:20].strip(), line[12:16].strip())
        if key in pdb_atoms:
            raise ValueError("Duplicate historical PDB atom")
        pdb_atoms[key] = np.asarray(
            [float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]]
        )
    if set(pdb_atoms) != set(requirements):
        raise ValueError("Historical PDB must contain exactly the mapped five atoms")
    xyz = parm.coordinates.copy()
    observations = []
    for number, resname, name in requirements:
        residue = parm.residues[number - 1]
        atoms = [atom for atom in residue.atoms if atom.name == name]
        if residue.name != resname or len(atoms) != 1:
            raise ValueError("Historical coordinate mapping is ambiguous")
        current, old = xyz[atoms[0].idx], pdb_atoms[(number, resname, name)]
        if (
            not np.isfinite([current, old]).all()
            or np.max(np.abs([current, old])) > 128
        ):
            raise ValueError(
                "Historical coordinate precision bound requires finite coordinates within128A"
            )
        difference = current - old
        if np.max(np.abs(difference)) > 0.00051:
            raise ValueError(
                "Historical PDB/NC coordinate bound exceeded; no alignment applied"
            )
        observations.append(
            {
                "atom_index": atoms[0].idx + 1,
                "residue": number,
                "atom": name,
                "NC_coordinates_A": current.tolist(),
                "PDB_coordinates_A": old.tolist(),
                "difference_A": difference.tolist(),
            }
        )
    points = [tuple(pdb_atoms[key]) for key in requirements]
    return {
        "phi_degrees": projected_herschlag_phi(*points),
        "coordinate_component_bound_A": 0.00051,
        "bound_reason": "0.0005A PDB rounding plus float32 write allowance, coordinates within128A",
        "coordinate_observations": observations,
        "maximum_component_difference_A": max(
            abs(v) for row in observations for v in row["difference_A"]
        ),
    }


METRICS = dict(
    zip(
        [
            "d1_thr929_og1_scissile_c",
            "d2_phe113_n_scissile_o",
            "d3_asn186_nd2_scissile_o",
            "a_BD_og1_c_o",
            "herschlag_phi_projected",
        ],
        [
            "d1_d_OG1_scissile_C",
            "d2_d_Phe113N_scissile_O",
            "d3_d_Asn186ND2_scissile_O",
            "a1_alpha_Burgi_Dunitz",
            "Herschlag_phi",
        ],
        strict=True,
    )
)


def _compare(metrics: dict, reference: dict, legacy_phi: float) -> dict:
    checks = {}
    for name, column in METRICS.items():
        expected, actual = float(reference[column]), float(metrics[name])
        if not math.isfinite(expected) or not math.isfinite(actual):
            raise ValueError(f"Nonfinite geometry metric {name}")
        raw_actual = actual
        if name == "herschlag_phi_projected":
            actual = float(legacy_phi)
            if not math.isfinite(actual):
                raise ValueError("Nonfinite legacyPDBphi")
        delta = actual - expected
        if name == "herschlag_phi_projected":
            delta = (delta + 180) % 360 - 180
        tolerance = 0.001 if name.startswith("d") else 0.01
        if name == "herschlag_phi_projected":
            tolerance = 0.00001
        checks[name] = {
            "reference": expected,
            "raw_NC_observed": raw_actual,
            "comparison_policy": "legacy3-decimalPDBphi"
            if name == "herschlag_phi_projected"
            else "rawNC",
            "raw_NC_difference": (raw_actual - expected + 180) % 360 - 180
            if name == "herschlag_phi_projected"
            else raw_actual - expected,
            "observed": actual,
            "difference": delta,
            "tolerance": tolerance,
            "passed": abs(delta) <= tolerance,
        }
    return checks
