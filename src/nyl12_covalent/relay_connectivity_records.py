"""Accepted relay record validation and mapped endpoint observables."""

from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Any

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

from .neb_umbrella_observables import frame_observables, topology_proton_pool
from .relay_connectivity_source import ROLES, TRACKED_H, digest, read_json


def _qm_atom_ids(topology: Path, mdin: Path) -> list[int]:
    parm = AmberParm(str(topology))
    source = mdin.read_text()
    masks = re.findall(r"\bqmmask\s*=\s*'([^']+)'", source, re.I)
    if len(masks) != 1:
        raise ValueError("Frozen dimer input must contain one QM mask")
    selection = AmberMask(parm, masks[0]).Selection()
    return [
        atom.idx + 1
        for atom, selected in zip(parm.atoms, selection, strict=True)
        if selected
    ]


def _frame_observables(
    topology: Path, restart_path: Path, qm_atom_ids: list[int]
) -> dict[str, Any]:
    restart = Rst7.open(str(restart_path))
    xyz = np.asarray(restart.coordinates, dtype=float).reshape((-1, 3))
    if len(xyz) != 104922 or restart.box is None:
        raise ValueError(f"Full restart with periodic cell required: {restart_path}")
    pool, initial_owners, owner_ids = topology_proton_pool(topology, qm_atom_ids)
    result = frame_observables(
        xyz,
        np.asarray(restart.box, dtype=float),
        ROLES,
        pool,
        owner_ids,
        initial_owners,
    )
    result["hydrogen_pool_1based"] = pool
    result["initial_owners_1based"] = initial_owners
    result["owner_atom_ids_1based"] = owner_ids
    result["off_path_changed_hydrogens"] = sorted(
        set(result["changed_hydrogen_ids"]) - TRACKED_H
    )
    return result


def _endpoint_chemistry_passed(
    geometry: dict[str, Any], branch: str, bias_kcal_mol: float
) -> bool:
    """Apply the frozen branch endpoint geometry and tracked-H gates."""
    if abs(bias_kcal_mol) > 1e-8:
        return False
    distances = geometry["distances_A"]
    ownership = geometry["ownership"]
    if any(
        len(ownership.get(hydrogen, [])) != 1 for hydrogen in ("H1", "H2", "HG", "HN")
    ):
        return False
    owners = [ownership[name][0] for name in ("H1", "H2", "HG")]
    if owners.count("N") != 2 or ownership["HN"] != ["NL"]:
        return False
    if not 1.1 <= distances["C_O"] <= 1.32:
        return False
    if branch == "reverse_activation":
        return (
            distances["C_OG"] >= 2.5
            and 1.2 <= distances["C_NL"] <= 1.65
            and owners.count("OG") == 1
        )
    if branch == "forward_acylation":
        return (
            1.25 <= distances["C_OG"] <= 1.65
            and distances["C_NL"] >= 2.2
            and owners.count("NL") == 1
        )
    raise ValueError(f"Unknown relay endpoint branch: {branch}")


def _accepted_rows(
    campaign: Path,
    branch_dir: Path,
    rows: list[dict[str, Any]],
    topology_hash: str,
    qm_atom_ids: list[int],
    chemistry: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    drive = branch_dir / "01_coordinate_drive"
    accepted: list[dict[str, Any]] = []
    receipts: list[dict[str, Any]] = []
    previous_path = campaign / "model/inputs/optimized.ncrst"
    previous_hash = digest(previous_path)
    for row in rows:
        label = str(row.get("label", ""))
        if label == "window033_start":
            if Path(row["restart"]).resolve() != previous_path.resolve():
                raise ValueError("Path baseline is not the declared window033 restart")
            seed_distances = row["geometry"]["distances_A"]
            accepted.append(
                {
                    "label": label,
                    "restart_path": str(previous_path),
                    "restart_sha256": previous_hash,
                    "predecessor_restart_path": None,
                    "predecessor_restart_sha256": None,
                    "geometry_label": row["geometry"].get("geometry_label"),
                    "distances_A": {
                        key: seed_distances[key]
                        for key in (
                            "C_OG",
                            "C_NL",
                            "C_O",
                            "H1_N",
                            "H1_OG",
                            "H1_NL",
                            "H2_N",
                            "H2_OG",
                            "H2_NL",
                            "HG_N",
                            "HG_OG",
                            "HG_NL",
                            "HN_N",
                            "HN_OG",
                            "HN_NL",
                        )
                    },
                    "tracked_hydrogen_ownership_1_30_A": row["geometry"]["ownership"],
                    "all_qm_proton_observables": _frame_observables(
                        campaign / "model/inputs/system.prmtop",
                        previous_path,
                        qm_atom_ids,
                    ),
                    "total_kcal_mol": float(row["total_kcal_mol"]),
                    "bias_kcal_mol": float(row["bias_kcal_mol"]),
                    "physical_kcal_mol": float(row["physical_kcal_mol"]),
                    "process_receipt": None,
                    "acceptance_receipt": str(
                        campaign / "model/inputs/scientific_acceptance.json"
                    ),
                }
            )
            continue
        if label == "unrestrained_endpoint":
            folder = branch_dir / "02_unrestrained_endpoint"
        elif label.startswith("window_"):
            folder = drive / label
        else:
            continue
        process_path, acceptance_path = (
            folder / "process.json",
            folder / "acceptance.json",
        )
        process = read_json(process_path)
        receipts.append(
            {
                "label": label,
                "process_path": str(process_path),
                "process_sha256": digest(process_path),
                "acceptance_path": str(acceptance_path)
                if acceptance_path.exists()
                else None,
                "acceptance_sha256": digest(acceptance_path)
                if acceptance_path.exists()
                else None,
            }
        )
        if row.get("accepted") is not True:
            continue
        if not acceptance_path.is_file():
            raise ValueError(f"Progress marks {label} accepted without acceptance.json")
        acceptance = read_json(acceptance_path)
        restart = Path(row["restart"]).resolve(strict=True)
        actual_restart = digest(restart)
        inputs = process.get("input_sha256", {})
        start_keys = [
            key for key in inputs if Path(key).resolve() == previous_path.resolve()
        ]
        topology_keys = [key for key in inputs if Path(key).name == "system.prmtop"]
        if len(start_keys) != 1 or inputs[start_keys[0]] != previous_hash:
            raise ValueError(f"Predecessor restart hash mismatch at {label}")
        if len(topology_keys) != 1 or inputs[topology_keys[0]] != topology_hash:
            raise ValueError(f"Topology hash mismatch at {label}")
        if process.get("exit_status") != 0:
            raise ValueError(
                f"Accepted window lacks successful process receipt: {label}"
            )
        if (
            acceptance.get("accepted") is not True
            or acceptance.get("native_converged") is not True
        ):
            raise ValueError(
                f"Accepted window lacks native/scientific acceptance: {label}"
            )
        if acceptance.get("electronic_history", {}).get("passed") is not True:
            raise ValueError(
                f"Accepted window lacks electronic-history acceptance: {label}"
            )
        if acceptance.get("nonrelay_hydrogen_audit", {}).get("passed") is not True:
            raise ValueError(f"Nonrelay hydrogen gate failed: {label}")
        if float(acceptance.get("frozen_max_abs_change_A", math.inf)) != 0.0:
            raise ValueError(f"Frozen exterior changed or is unverified: {label}")
        if Path(acceptance["restart"]).resolve() != restart:
            raise ValueError(f"Acceptance restart differs from progress row: {label}")
        if acceptance.get("geometry") != row.get("geometry"):
            raise ValueError(f"Progress geometry differs from acceptance: {label}")
        total = float(acceptance["total_kcal_mol"])
        bias = float(acceptance["bias_kcal_mol"])
        physical = float(acceptance["physical_kcal_mol"])
        if not all(math.isfinite(value) for value in (total, bias, physical)):
            raise ValueError(f"Nonfinite energy record: {label}")
        if abs(total - bias - physical) > 0.002:
            raise ValueError(f"Energy decomposition mismatch: {label}")
        distances = row["geometry"]["distances_A"]
        ownership = row["geometry"]["ownership"]
        observables = _frame_observables(
            campaign / "model/inputs/system.prmtop", restart, qm_atom_ids
        )
        accepted.append(
            {
                "label": label,
                "restart_path": str(restart),
                "restart_sha256": actual_restart,
                "predecessor_restart_path": str(previous_path),
                "predecessor_restart_sha256": previous_hash,
                "geometry_label": row["geometry"].get("geometry_label"),
                "distances_A": {
                    key: distances[key]
                    for key in (
                        "C_OG",
                        "C_NL",
                        "C_O",
                        "H1_N",
                        "H1_OG",
                        "H1_NL",
                        "H2_N",
                        "H2_OG",
                        "H2_NL",
                        "HG_N",
                        "HG_OG",
                        "HG_NL",
                        "HN_N",
                        "HN_OG",
                        "HN_NL",
                    )
                },
                "tracked_hydrogen_ownership_1_30_A": ownership,
                "all_qm_proton_observables": observables,
                "total_kcal_mol": total,
                "bias_kcal_mol": bias,
                "physical_kcal_mol": physical,
                "unbiased_endpoint_candidate_passed": None,
                "independent_energy_hartree": None,
                "maximum_gradient_hartree_per_bohr": None,
                "process_receipt": receipts[-1],
                "acceptance_receipt": receipts[-1]["acceptance_path"],
            }
        )
        if label == "unrestrained_endpoint":
            independent_energy = float(acceptance["independent_energy_hartree"])
            maximum_gradient = float(acceptance["max_gradient_hartree_per_bohr"])
            if not math.isfinite(independent_energy) or not math.isfinite(
                maximum_gradient
            ):
                raise ValueError(
                    "Unbiased endpoint has nonfinite independent force receipt"
                )
            accepted[-1]["independent_energy_hartree"] = independent_energy
            accepted[-1]["maximum_gradient_hartree_per_bohr"] = maximum_gradient
            accepted[-1]["unbiased_endpoint_candidate_passed"] = (
                _endpoint_chemistry_passed(acceptance["geometry"], chemistry, bias)
                and maximum_gradient <= 1e-5
            )
        previous_path, previous_hash = restart, actual_restart
    return accepted, receipts
