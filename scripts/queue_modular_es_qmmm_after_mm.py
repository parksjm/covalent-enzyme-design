#!/usr/bin/env python3
"""Conditionally run modular-PA6 DFTB3/MM after an accepted MM endpoint."""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

import numpy as np
import parmed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.active_site_geometry import (  # noqa: E402
    angle_degrees,
    distance,
    projected_herschlag_phi,
)
from nyl12_covalent.catalytic_alignment import align_catalytic_site  # noqa: E402
from nyl12_covalent.es_dftb_unrestrained import prepare_dftb_optimization  # noqa: E402
from nyl12_covalent.modular_es_qm_region import build_modular_es_region  # noqa: E402
from nyl12_covalent.modular_es_snapshot_selection import classify_geometry  # noqa: E402
from nyl12_covalent.ti1_checks import check_convergence  # noqa: E402
from nyl12_covalent.unrestrained_relay import native_converged  # noqa: E402

LOGGER = logging.getLogger(__name__)


def _point(structure: Any, residue: int, atom_name: str) -> tuple[float, float, float]:
    matches = [
        atom for atom in structure.residues[residue - 1].atoms if atom.name == atom_name
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one atom at residue {residue} named {atom_name}")
    xyz = np.asarray(structure.coordinates[matches[0].idx], dtype=float)
    if xyz.shape != (3,) or not np.isfinite(xyz).all():
        raise ValueError("Endpoint coordinates must be finite")
    return tuple(float(value) for value in xyz)


def endpoint_geometry(topology: Path, coordinates: Path) -> dict[str, Any]:
    """Measure the five E-S metrics and Thr attack-orientation diagnostic."""
    structure = parmed.load_file(str(topology), xyz=str(coordinates))
    og1 = _point(structure, 929, "OG1")
    cb = _point(structure, 929, "CB")
    carbon = _point(structure, 1314, "C")
    oxygen = _point(structure, 1314, "O")
    nitrogen = _point(structure, 1315, "N")
    adjacent = _point(structure, 1314, "C5")
    metrics = {
        "d1": distance(og1, carbon),
        "d2": distance(_point(structure, 113, "N"), oxygen),
        "d3": distance(_point(structure, 186, "ND2"), oxygen),
        "attack": angle_degrees(og1, carbon, oxygen),
        "ligand_rmsd": 0.0,
        "projected_phi": projected_herschlag_phi(
            og1, oxygen, carbon, nitrogen, adjacent
        ),
        "cb_og1_carbon_angle": angle_degrees(cb, og1, carbon),
    }
    if not all(math.isfinite(value) for value in metrics.values()):
        raise ValueError("Nonfinite endpoint metric")
    return {"metrics": metrics, "tier": classify_geometry(metrics) or "none"}


def native_gate(alignment: dict[str, Any]) -> dict[str, Any]:
    """Apply the predeclared held-out Thr233 native-like limits."""
    state = alignment.get("states", alignment)["continuation_endpoint"]
    thr = state["thr233_held_out"]
    displacement = thr["heavy_atom_displacements_angstrom"]
    checks = {
        "N_displacement_le_1.25_A": displacement["N"] <= 1.25,
        "OG1_displacement_le_0.75_A": displacement["OG1"] <= 0.75,
        "abs_chi1_delta_le_20_deg": abs(thr["chi1_delta_es_minus_reference_degrees"])
        <= 20.0,
        "abs_psi_delta_le_35_deg": abs(thr["psi_delta_es_minus_reference_degrees"])
        <= 35.0,
    }
    return {"accepted": all(checks.values()), "checks": checks, "metrics": thr}


def validate_override_options(
    og1_displacement_max: float | None, override_reason: str | None
) -> str | None:
    """Validate paired explicit near-pass override options.

    Args:
        og1_displacement_max: Optional maximum allowed OG1 displacement in A.
        override_reason: Optional nonempty explanation for using the override.

    Returns:
        The stripped override reason, or ``None`` when no override was supplied.

    Raises:
        ValueError: If only one option is supplied or either value is invalid.
    """
    if (og1_displacement_max is None) != (override_reason is None):
        raise ValueError(
            "--allow-og1-displacement-max and --override-reason must be supplied together"
        )
    if og1_displacement_max is None:
        return None
    if not math.isfinite(og1_displacement_max) or og1_displacement_max <= 0:
        raise ValueError(
            "OG1 displacement override maximum must be finite and positive"
        )
    reason = override_reason.strip()
    if not reason:
        raise ValueError("--override-reason must be nonempty")
    return reason


def mm_acceptance_decision(
    convergence: dict[str, Any],
    geometry: dict[str, Any],
    native: dict[str, Any],
    og1_displacement_max: float | None,
    override_reason: str | None,
) -> dict[str, Any]:
    """Build standard and explicitly overridden MM-to-QM/MM decisions."""
    reason = validate_override_options(og1_displacement_max, override_reason)
    standard_accepted = (
        convergence["accepted"]
        and geometry["tier"] in {"A", "B"}
        and native["accepted"]
    )
    override_failures: list[str] = []
    override_accepted = False
    if og1_displacement_max is not None:
        if not convergence["accepted"]:
            override_failures.append("convergence")
        if geometry["tier"] not in {"A", "B"}:
            override_failures.append("geometry_tier")
        non_og1_checks = {
            name: passed
            for name, passed in native["checks"].items()
            if name != "OG1_displacement_le_0.75_A"
        }
        if not all(non_og1_checks.values()):
            override_failures.extend(
                name for name, passed in non_og1_checks.items() if not passed
            )
        og1_displacement = native["metrics"]["heavy_atom_displacements_angstrom"]["OG1"]
        if og1_displacement > og1_displacement_max:
            override_failures.append("OG1_displacement_exceeds_override_max")
        override_accepted = not override_failures
    return {
        "standard_accepted_for_qmmm": standard_accepted,
        "override": {
            "requested": og1_displacement_max is not None,
            "accepted": override_accepted,
            "reason": reason,
            "og1_displacement_max_angstrom": og1_displacement_max,
            "failed_conditions": override_failures,
        },
        "accepted_for_qmmm": standard_accepted or override_accepted,
    }


def _wait_for_status(path: Path, poll_seconds: int) -> dict[str, Any]:
    while True:
        if path.is_file():
            status = json.loads(path.read_text(encoding="utf-8"))
            if status.get("ended_utc") is not None:
                return status
        time.sleep(poll_seconds)


def main() -> None:
    """Wait for MM, gate its endpoint, then run and assess one local QM/MM shell."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-status", required=True, type=Path)
    parser.add_argument("--mm-directory", required=True, type=Path)
    parser.add_argument("--initial-coordinates", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--reference-cif", required=True, type=Path)
    parser.add_argument("--region-template", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--active-radius", required=True, type=float)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--allow-og1-displacement-max", type=float)
    parser.add_argument("--override-reason")
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    override_reason = validate_override_options(
        args.allow_og1_displacement_max, args.override_reason
    )
    if args.output.exists() or args.active_radius <= 0 or args.poll_seconds < 1:
        raise ValueError("Fresh output, positive radius, and poll interval required")
    args.output.mkdir(parents=True)
    status = _wait_for_status(args.upstream_status.resolve(), args.poll_seconds)
    if status.get("exit_status") != 0:
        raise RuntimeError("Upstream MM job did not exit successfully")
    mm_coordinates = args.mm_directory / "minimized.ncrst"
    convergence = check_convergence(args.mm_directory / "minimized.mdout", drms=0.01)
    alignment = align_catalytic_site(
        args.reference_cif.resolve(),
        args.topology.resolve(),
        args.initial_coordinates.resolve(),
        mm_coordinates.resolve(),
        args.output / "mm_endpoint_alignment",
    )
    mm_geometry = endpoint_geometry(args.topology.resolve(), mm_coordinates.resolve())
    mm_native = native_gate(alignment)
    mm_acceptance = {
        "convergence": convergence,
        "geometry": mm_geometry,
        "native_thr233": mm_native,
    }
    mm_acceptance.update(
        mm_acceptance_decision(
            convergence,
            mm_geometry,
            mm_native,
            args.allow_og1_displacement_max,
            override_reason,
        )
    )
    (args.output / "mm_acceptance.json").write_text(
        json.dumps(mm_acceptance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if not mm_acceptance["accepted_for_qmmm"]:
        raise RuntimeError("MM endpoint failed standard and explicit QM/MM gates")

    region = args.output / "qm_region.json"
    build_modular_es_region(
        args.region_template.resolve(),
        args.topology.resolve(),
        mm_coordinates.resolve(),
        region,
    )
    active_mask = f"(:113,186,929,930,968,970,1313-1316)<:{args.active_radius:.1f}"
    qmmm_dir = args.output / "dftb3_mm"
    prepared = prepare_dftb_optimization(
        region,
        args.sander.resolve(),
        qmmm_dir,
        active_mask=active_mask,
        case_label=f"WT modular PA6 E-S {args.active_radius:.1f} A movable shell",
    )
    environment = dict(os.environ)
    environment.update(prepared["environment"])
    process = subprocess.run(prepared["command"], env=environment, check=False)
    if process.returncode:
        raise RuntimeError(f"DFTB3/MM engine exited with status {process.returncode}")
    optimized = qmmm_dir / "optimized.ncrst"
    final_alignment = align_catalytic_site(
        args.reference_cif.resolve(),
        args.topology.resolve(),
        mm_coordinates.resolve(),
        optimized,
        args.output / "qmmm_endpoint_alignment",
    )
    final_geometry = endpoint_geometry(args.topology.resolve(), optimized)
    final_native = native_gate(final_alignment)
    acceptance = {
        "native_converged": native_converged(qmmm_dir / "minimized.mdout"),
        "geometry": final_geometry,
        "native_thr233": final_native,
    }
    acceptance["catalytically_competent_es"] = (
        acceptance["native_converged"]
        and final_geometry["tier"] in {"A", "B"}
        and final_native["accepted"]
    )
    (args.output / "qmmm_acceptance.json").write_text(
        json.dumps(acceptance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if not acceptance["catalytically_competent_es"]:
        raise RuntimeError("QM/MM endpoint failed scientific acceptance")


if __name__ == "__main__":
    main()
