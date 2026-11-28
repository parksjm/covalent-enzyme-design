"""Dispatch PI-authorized exploratory relay branches from force-accepted window033."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import sys
from typing import Any

LOGGER = logging.getLogger(__name__)
BRANCHES = {
    "reverse_activation": "01_activation_reverse",
    "forward_acylation": "02_acylation_forward",
}
FORCE_LIMIT_HARTREE_PER_BOHR = 1e-5
AUTHORIZATION_FILENAME = "exploratory_authorization_20261124.json"
PLAN_FILENAME = "exploratory_plan_20261124.json"
DISPATCH_FILENAME = "exploratory_dispatch_20261124.json"


def digest(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Atomically write one dispatch record."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def read_json(path: Path) -> dict[str, Any]:
    """Read a JSON object and fail on non-object content."""
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def _verify_frozen_files(campaign: Path, config: dict[str, Any]) -> None:
    """Verify every content hash named by the campaign and frozen model."""
    preparation = read_json(campaign / "preparation.json")
    config_path = campaign / "config.json"
    if digest(config_path) != preparation["config_sha256"]:
        raise ValueError("Campaign config differs from the prepared campaign")
    model = campaign / config["model"]
    manifest_path = model / "manifest.json"
    if digest(manifest_path) != preparation["source_frozen_manifest_sha256"]:
        raise ValueError("Frozen model manifest differs from preparation record")
    manifest = read_json(manifest_path)
    for relative, expected in manifest.items():
        path = model / relative
        if digest(path) != expected:
            raise ValueError(f"Frozen model file hash mismatch: {path}")
    for name, expected in config["source_sha256"].items():
        path = Path(name)
        if digest(path) != expected:
            raise ValueError(f"Frozen campaign source hash mismatch: {path}")


def validate_dispatch(campaign: Path, branch: str) -> dict[str, Any]:
    """Validate authorization, immutable inputs, and force acceptance.

    Args:
        campaign: N-terminal relay campaign directory.
        branch: Exact branch key from the frozen configuration.

    Returns:
        A compact preflight report for an exploratory branch dispatch.

    Raises:
        ValueError: If the branch, authorization, source hashes, or force gate fail.
    """
    if branch not in BRANCHES:
        raise ValueError(f"Unsupported relay branch: {branch}")
    campaign = campaign.resolve(strict=True)
    config = read_json(campaign / "config.json")
    if config["branches"] != BRANCHES:
        raise ValueError("Campaign branch map differs from frozen definitions")
    _verify_frozen_files(campaign, config)

    branch_root = campaign / config["branches"][branch]
    authorization_path = branch_root / AUTHORIZATION_FILENAME
    plan_path = branch_root / PLAN_FILENAME
    if not authorization_path.is_file() or not plan_path.is_file():
        raise ValueError("Exploratory authorization or branch plan is missing")
    authorization = read_json(authorization_path)
    plan = read_json(plan_path)
    dispatcher_path = Path(plan.get("dispatcher", ""))
    if (
        authorization.get("status") != "approved_exploratory_only"
        or branch not in authorization.get("authorized_branches", [])
        or plan.get("authorization_sha256") != digest(authorization_path)
        or plan.get("branch") != branch
        or plan.get("config_sha256") != digest(campaign / "config.json")
        or not dispatcher_path.is_file()
        or digest(dispatcher_path) != plan.get("dispatcher_sha256")
    ):
        raise ValueError("Exploratory authorization or branch plan is invalid")

    model = campaign / config["model"]
    identity = read_json(model / "source_identity.json")
    acceptance = read_json(model / "release_acceptance.json")
    if (
        identity.get("seed") != "ti1_window033"
        or identity.get("full_atom_count") != 104922
        or identity.get("active_atom_count") != 4165
        or identity.get("qm_atom_count") != 98
        or identity.get("qm_charge") != -2
        or digest(model / "inputs/system.prmtop")
        != identity.get("topology_sha256")
        or digest(model / "inputs/system.prmtop")
        != config["source_topology_sha256"]
    ):
        raise ValueError("Window033 source identity differs from approved model")
    if (
        acceptance.get("accepted") is not True
        or acceptance.get("force_converged") is not True
        or acceptance.get("native_converged") is not True
        or acceptance.get("max_gradient_hartree_per_bohr", float("inf"))
        > FORCE_LIMIT_HARTREE_PER_BOHR
    ):
        raise ValueError("Window033 force-acceptance requirements are not met")
    return {
        "status": "preflight_passed_exploratory_only",
        "branch": branch,
        "branch_directory": str(branch_root),
        "window033_force_accepted": True,
        "full_minimum_verified": False,
        "ts_claim_authorized": False,
        "config_sha256": digest(campaign / "config.json"),
        "restart_sha256": digest(model / "inputs/optimized.ncrst"),
        "topology_sha256": digest(model / "inputs/system.prmtop"),
    }


def run_exploratory_branch(campaign: Path, branch: str) -> dict[str, Any]:
    """Run one authorized branch through the frozen native branch implementation.

    The existing implementation may generate a conditional dimer-search
    candidate after a chemically accepted released endpoint. That search does
    not certify a stationary point, TS, or barrier.

    Args:
        campaign: N-terminal relay campaign directory.
        branch: Exact branch key from the frozen configuration.

    Returns:
        Branch status returned by the frozen implementation.
    """
    report = validate_dispatch(campaign, branch)
    campaign = campaign.resolve(strict=True)
    branch_root = campaign / BRANCHES[branch]
    receipt_path = branch_root / DISPATCH_FILENAME
    if receipt_path.exists():
        raise FileExistsError(f"Preserve existing exploratory dispatch record: {receipt_path}")
    source = campaign.resolve(strict=True) / "source_snapshot"
    sys.path.insert(0, str(source / "src"))
    sys.path.insert(0, str(source))
    from nyl12_covalent.nterminal_relay_run import run_branch

    receipt: dict[str, Any] = {
        **report,
        "status": "dispatch_started",
        "started_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "full_minimum_verified": False,
        "formal_ts_claim": False,
    }
    write_json(receipt_path, receipt)
    LOGGER.info("Dispatching exploratory-only branch %s", branch)
    try:
        result = run_branch(campaign, branch)
    except Exception as error:
        receipt.update(
            status="dispatch_failed",
            ended_utc=datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            error=f"{type(error).__name__}: {error}",
        )
        write_json(receipt_path, receipt)
        raise
    receipt.update(
        status="path_implementation_returned",
        ended_utc=datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        implementation_result=result,
        full_minimum_verified=False,
        formal_ts_claim=False,
    )
    write_json(receipt_path, receipt)
    return result


def main() -> int:
    """Preflight or explicitly launch one approved exploratory branch."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--branch", choices=tuple(BRANCHES), required=True)
    parser.add_argument(
        "--launch",
        action="store_true",
        help="Launch the native pathway after preflight; omit for read-only checks.",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    if not args.launch:
        report = validate_dispatch(args.campaign, args.branch)
        LOGGER.info("%s", json.dumps(report, sort_keys=True))
        return 0
    result = run_exploratory_branch(args.campaign, args.branch)
    LOGGER.info("%s", json.dumps(result, sort_keys=True))
    return 0 if result.get("status") == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
