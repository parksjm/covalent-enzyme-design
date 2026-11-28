"""Six reviewed independent starts at unchanged main-campaign two-CV centers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .relay_activation_support import prepare_reviewed_targets
from .relay_activation_support_runtime import run_reviewed_target
from .neb_umbrella_campaign import digest

TARGETS = {
    "A026_original": "activation_full_026",
    "A026_late": "activation_full_026",
    "A028_original": "activation_full_028",
    "A028_late": "activation_full_028",
    "D008_from_D010": "delivery_full_008",
    "D009_from_D007": "delivery_full_009",
}


def validate_mixing_targets(
    targets: list[dict[str, Any]], parent: dict[str, Any]
) -> None:
    """Require the six approved starts at their exact unchanged native centers."""
    if len(targets) != 6 or {r["label"] for r in targets} != set(TARGETS):
        raise ValueError("Exactly the six reviewed independent starts are required")
    windows = {r["label"]: r for r in parent["windows"]}
    for row in targets:
        native = windows[TARGETS[row["label"]]]["bias"]
        if (
            row["feature_names"] != native["feature_names"]
            or row["centers_A"] != native["centers_A"]
            or row["k"] != 25
            or native["force_constants"] != [25.0, 25.0]
        ):
            raise ValueError("Mixing checks must retain exact main-window native bias")


def prepare_mixing_checks(plan_path: Path) -> Path:
    """Prepare exactly six independently initialized reviewed native states."""
    plan = json.loads(plan_path.read_text())
    parent_path = Path(plan["parent_prepared"])
    parent = json.loads(parent_path.read_text())
    if digest(parent_path) != parent_path.with_suffix(".sha256").read_text().strip():
        raise ValueError("Parent hash mismatch")
    return prepare_reviewed_targets(
        plan_path,
        "mixing_checks_20261127",
        lambda targets: validate_mixing_targets(targets, parent),
    )


def run_mixing_check(prepared: Path, label: str, expected_sha256: str) -> None:
    """Run one frozen independent start with native acceptance at every stage."""
    cfg = json.loads(prepared.read_text())
    parent = json.loads(Path(cfg["parent_prepared"]).read_text())
    validate_mixing_targets(cfg["targets_plan"]["targets"], parent)
    run_reviewed_target(prepared, label, expected_sha256, 6)


def main() -> None:
    """Prepare a reviewed plan or execute one explicitly hash-pinned task."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--prepared", type=Path)
    parser.add_argument("--label")
    parser.add_argument("--sha256")
    args = parser.parse_args()
    if args.plan and not any((args.prepared, args.label, args.sha256)):
        prepare_mixing_checks(args.plan)
    elif args.prepared and args.label and args.sha256 and not args.plan:
        run_mixing_check(args.prepared, args.label, args.sha256)
    else:
        parser.error("Specify --plan or all of --prepared/--label/--sha256")


if __name__ == "__main__":
    main()
