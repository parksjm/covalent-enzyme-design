#!/usr/bin/env python3
"""Run one selected candidate through segmented MM and a gated QM/MM shell."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
import math
from pathlib import Path
import subprocess
import sys
from collections.abc import Callable, Sequence
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.candidate_mm_qmmm_controller import run_candidate  # noqa: E402
from nyl12_covalent.catalytic_alignment import align_catalytic_site  # noqa: E402
from nyl12_covalent.asp_relay_gate import (  # noqa: E402
    relay_contact_gate,
    relay_distances,
)
from scripts.queue_modular_es_qmmm_after_mm import (  # noqa: E402
    endpoint_geometry,
    native_gate,
)

LOGGER = logging.getLogger(__name__)


def validate_active_radii(active_radii: Sequence[float]) -> None:
    """Require at least one unique, finite, positive active-shell radius."""
    if not active_radii:
        raise ValueError("At least one --active-radius is required")
    if any(not math.isfinite(radius) or radius <= 0 for radius in active_radii):
        raise ValueError("Active radii must be finite and positive")
    if len(set(active_radii)) != len(active_radii):
        raise ValueError("Active radii must not contain duplicates")


def launch_qmmm_shells(
    queue_script: Path,
    common_arguments: Sequence[str],
    active_radii: Sequence[float],
    output: Path,
    cwd: Path,
    runner: Callable[[Sequence[str], Path], int] | None = None,
) -> dict[str, Any]:
    """Launch and record one isolated QM/MM queue per declared radius.

    Each queue consumes the same already accepted MM endpoint supplied in
    ``common_arguments``. A failed shell is recorded and stops later launches.
    """
    validate_active_radii(active_radii)
    execute = runner or _run_queue
    output.mkdir(parents=True, exist_ok=False)
    report_path = output / "qmmm_shell_launch_report.json"
    shell_results: list[dict[str, Any]] = []

    for radius in active_radii:
        radius_label = f"{radius:g}".replace(".", "p")
        shell_output = output / f"qmmm_dftb3_shell_{radius_label}A"
        command = [
            sys.executable,
            str(queue_script),
            *common_arguments,
            "--output",
            str(shell_output),
            "--active-radius",
            str(radius),
        ]
        try:
            exit_status = execute(command, cwd)
        except OSError as error:
            shell_results.append(
                {
                    "active_radius_A": radius,
                    "output_directory": str(shell_output),
                    "command": command,
                    "launch_error": str(error),
                }
            )
            _write_shell_report(report_path, shell_results)
            raise
        shell_results.append(
            {
                "active_radius_A": radius,
                "output_directory": str(shell_output),
                "command": command,
                "exit_status": exit_status,
            }
        )
        _write_shell_report(report_path, shell_results)
        if exit_status != 0:
            raise RuntimeError(
                f"DFTB3/MM queue for {radius:g} A shell exited with status "
                f"{exit_status}"
            )

    return {"shells": shell_results, "exit_status": 0}


def _run_queue(command: Sequence[str], cwd: Path) -> int:
    """Run one single-radius QM/MM queue without altering its arguments."""
    return subprocess.run(list(command), cwd=cwd, check=False).returncode


def _write_shell_report(path: Path, shells: list[dict[str, Any]]) -> None:
    """Persist each shell command and process result as it completes."""
    path.write_text(
        json.dumps({"shells": shells}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    """Parse one candidate definition and run its conditional pipeline."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--engine", required=True, type=Path)
    parser.add_argument("--reference-cif", required=True, type=Path)
    parser.add_argument("--region-template", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--case-label", required=True)
    parser.add_argument("--segment-cycles", type=int, default=10_000)
    parser.add_argument("--total-cycle-ceiling", type=int, required=True)
    parser.add_argument("--active-radius", type=float, action="append", required=True)
    parser.add_argument(
        "--endpoint-gate",
        choices=("native_thr", "asp272_relay"),
        default="native_thr",
        help="Scientific gate applied after unrestrained MM convergence.",
    )
    args = parser.parse_args(argv)
    try:
        validate_active_radii(args.active_radius)
    except ValueError as error:
        parser.error(str(error))
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )

    topology = args.topology.resolve(strict=True)
    initial_coordinates = args.coordinates.resolve(strict=True)
    reference_cif = args.reference_cif.resolve(strict=True)
    region_template = args.region_template.resolve(strict=True)
    sander = args.sander.resolve(strict=True)

    def assess(endpoint: Path) -> dict[str, Any]:
        alignment = align_catalytic_site(
            reference_cif,
            topology,
            initial_coordinates,
            endpoint,
            args.output / "mm_endpoint_alignment",
        )
        geometry = endpoint_geometry(topology, endpoint)
        thr_gate = native_gate(alignment)
        relay_gate = relay_contact_gate(relay_distances(topology, endpoint))
        state_gate = thr_gate if args.endpoint_gate == "native_thr" else relay_gate
        return {
            "accepted": geometry["tier"] in {"A", "B"} and state_gate["accepted"],
            "endpoint_gate": args.endpoint_gate,
            "geometry": geometry,
            "native_thr233": thr_gate,
            "asp272_relay": relay_gate,
        }

    def launch_qmmm(mm_directory: Path, qmmm_output: Path) -> dict[str, Any]:
        upstream_status = args.output / "mm_gate_passed_status.json"
        upstream_status.write_text(
            json.dumps(
                {
                    "ended_utc": datetime.now(timezone.utc).isoformat(),
                    "exit_status": 0,
                    "gate": f"Tier A/B plus {args.endpoint_gate} passed",
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        common_arguments = [
            "--upstream-status",
            str(upstream_status),
            "--mm-directory",
            str(mm_directory),
            "--initial-coordinates",
            str(initial_coordinates),
            "--topology",
            str(topology),
            "--reference-cif",
            str(reference_cif),
            "--region-template",
            str(region_template),
            "--sander",
            str(sander),
        ]
        return launch_qmmm_shells(
            queue_script=ROOT / "scripts" / "queue_modular_es_qmmm_after_mm.py",
            common_arguments=common_arguments,
            active_radii=args.active_radius,
            output=qmmm_output,
            cwd=ROOT,
        )

    report = run_candidate(
        topology=topology,
        initial_coordinates=initial_coordinates,
        engine=args.engine.resolve(strict=True),
        output=args.output.resolve(),
        case_label=args.case_label,
        segment_cycles=args.segment_cycles,
        total_cycle_ceiling=args.total_cycle_ceiling,
        endpoint_gate=assess,
        qmmm_launcher=launch_qmmm,
    )
    LOGGER.info("Candidate controller completed: %s", report["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
