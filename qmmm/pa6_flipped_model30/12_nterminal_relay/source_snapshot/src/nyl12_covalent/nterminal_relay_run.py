"""Bounded native relay branches with explicit scientific stop conditions."""

from __future__ import annotations

import json
import logging
from pathlib import Path
import re
from typing import Any

import numpy as np
from parmed.amber import Rst7

from nyl12_covalent.direct_proton_transfer import driven_input
from nyl12_covalent.stationary_native import write_json
from nyl12_covalent.nterminal_relay_coordinates import (
    amber_restraint,
    candidate,
    coordinate,
    geometry,
    schedule,
)
from nyl12_covalent.nterminal_relay_progress import publish
from nyl12_covalent.relay_release import minimization_input
from nyl12_covalent.stationary_native import HARTREE_TO_KCAL_PER_MOL, NativeState

LOGGER = logging.getLogger(__name__)
_BRANCHES = {"reverse_activation", "forward_acylation"}


def _restart_xyz(path: Path) -> np.ndarray:
    return np.asarray(Rst7.open(str(path)).coordinates).reshape((-1, 3))


def _driven_mdin(source: str, maxcycles: int) -> str:
    base = minimization_input(source, maxcycles=maxcycles)
    base, count = re.subn(r"\btol\s*=\s*1\.d-5", "tol=1.d-4", base, count=1, flags=re.I)
    if count != 1 or not re.search(r"\btole\s*=\s*1\.d-7", base, re.I):
        raise ValueError("Unable to establish requested driven-window tolerances")
    return driven_input(base, max_cycles=maxcycles)


def run_branch(campaign: Path, branch: str) -> dict[str, Any]:
    """Drive, release and assess one conditionally authorized relay branch.

    Args:
        campaign: Campaign root containing config.json and frozen source model.
        branch: reverse_activation or forward_acylation.

    Returns:
        Scientific branch status, with TS result only after an accepted release.

    Raises:
        ValueError: If the branch, source or schedule fails its input contract.
        Exception: If a native calculation or publication fails; status is saved.
    """
    if branch not in _BRANCHES:
        raise ValueError(f"Unknown branch: {branch}")
    from nyl12_covalent.nterminal_relay_native import run_native
    from nyl12_covalent.nterminal_relay_ts import run_ts

    campaign = campaign.resolve(strict=True)
    config = json.loads((campaign / "config.json").read_text())
    root = (campaign / config["branches"][branch]).resolve()
    if not root.is_relative_to(campaign) or root == campaign:
        raise ValueError("Branch directory must be inside campaign")
    root.mkdir(parents=True, exist_ok=True)
    status_path = root / "branch_status.json"
    status: dict[str, Any] = {"branch": branch, "status": "running"}
    write_json(status_path, status)
    try:
        model = (campaign / config["model"]).resolve(strict=True)
        state = NativeState(model)
        proton = int(config.get("proton", 13917))
        start = model / "inputs/optimized.ncrst"
        start_geometry = geometry(state.xyz)
        initial = coordinate(state.xyz, branch, proton)
        energy = float(config["reference_energy_hartree"]) * HARTREE_TO_KCAL_PER_MOL
        rows: list[dict[str, Any]] = [
            {
                "label": "window033_start",
                "restart": str(start.resolve(strict=True)),
                "accepted": True,
                "native_converged": True,
                "geometry": start_geometry,
                "total_kcal_mol": energy,
                "bias_kcal_mol": 0.0,
                "physical_kcal_mol": energy,
                "target_A": initial,
                "observed_coordinate_A": initial,
            }
        ]
        windows = schedule(
            initial,
            final=float(config.get("path_final_A", 4.2)),
            maxincrement=float(config.get("path_maxincrement_A", 0.15)),
        )
        mdin = _driven_mdin(state.mdin_source, int(config.get("path_maxcycles", 4000)))
        found = False
        previous = start
        for item in windows:
            label, target = item["window"], item["target"]
            result = run_native(
                campaign,
                root / "01_coordinate_drive" / label,
                previous,
                mdin,
                restraints=amber_restraint(target, branch, proton),
            )
            restart = Path(result["restart"]).resolve(strict=True)
            result.update(
                label=label,
                target_A=target,
                observed_coordinate_A=coordinate(_restart_xyz(restart), branch, proton),
                restart=str(restart),
            )
            rows.append(result)
            publish(campaign, root, rows)
            if not result["accepted"]:
                status.update(
                    status="failed", reason="native_window_rejected", label=label
                )
                write_json(status_path, status)
                return status
            previous = restart
            if candidate(result["geometry"], branch):
                found = True
                break
        if not found:
            status.update(
                status="failed", reason="schedule_exhausted_without_candidate"
            )
            write_json(status_path, status)
            return status
        release = run_native(
            campaign,
            root / "02_unrestrained_endpoint",
            previous,
            minimization_input(state.mdin_source, maxcycles=10000),
            restraints=None,
        )
        release_restart = Path(release["restart"]).resolve(strict=True)
        release.update(
            label="unrestrained_endpoint",
            target_A=None,
            observed_coordinate_A=coordinate(
                _restart_xyz(release_restart), branch, proton
            ),
            restart=str(release_restart),
        )
        rows.append(release)
        publish(campaign, root, rows)
        if not release["accepted"] or not candidate(release["geometry"], branch):
            status.update(status="failed", reason="unrestrained_endpoint_not_retained")
            write_json(status_path, status)
            return status
        ts = run_ts(campaign, root, rows)
        status.update(status="completed", endpoint_restart=str(release_restart), ts=ts)
        write_json(status_path, status)
        return status
    except Exception as exc:
        status.update(
            status="failed", reason="exception", error=f"{type(exc).__name__}: {exc}"
        )
        write_json(status_path, status)
        LOGGER.exception("Relay branch %s failed", branch)
        raise
