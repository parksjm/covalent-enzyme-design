"""Wait for two actual accepted TS parents, then run only bounded broad generation."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import subprocess
import time

from nyl12_covalent.relay_broad_release import (
    FOUR_TEMPLATE_PROTOCOL,
    validate_generation_request,
    validate_idle_gpu,
    validate_original_lineage,
)
from nyl12_covalent.relay_broad_runtime import (
    digest,
    journal,
    native_run,
    verify_parent,
    write_new,
)

LOGGER = logging.getLogger(__name__)


def check_hashes(hashes: dict[str, str]) -> None:
    """Fail if any frozen source or input changed."""
    for path, expected in hashes.items():
        if digest(Path(path)) != expected:
            raise ValueError(f"Frozen generation queue source/input changed: {path}")


def execute(config: dict) -> None:
    """Release one 24-sample pilot after native, graph, source and GPU gates.

    Args:
        config: Explicit immutable queue configuration.
    """
    root = Path(config["repository_root"])
    output = Path(config["output"])
    check_hashes(config["source_sha256"])
    request = json.loads(Path(config["request"]).read_text())
    protection = json.loads(Path(config["protection"]).read_text())
    validate_generation_request(request, protection)
    gate_path, terminal_path = Path(config["gate"]), Path(config["terminal"])
    deadline = time.monotonic() + config["max_wait_seconds"]
    while not terminal_path.exists():
        if time.monotonic() >= deadline:
            raise TimeoutError("Finite native-acceptance wait exhausted")
        LOGGER.info(
            "Waiting at single parent03 native-plus-scientific acceptance boundary"
        )
        time.sleep(60)
    terminal = json.loads(terminal_path.read_text())
    if terminal["accepted"] is not True or terminal["terminal"] is not True:
        raise RuntimeError(
            "Parent03 native/scientific qualification failed; no generation"
        )
    gate = json.loads(gate_path.read_text())
    if gate["accepted"] is not True:
        raise ValueError("Two-original-parent manifest is not accepted")
    validate_original_lineage(gate["parents"])
    for parent in gate["parents"]:
        verify_parent(parent)
        for field in ("authoritative_coordinates", "explicit_graph_cif"):
            if digest(Path(parent[field])) != parent[field + "_sha256"]:
                raise ValueError(
                    "Accepted authoritative coordinate/graph artifact changed"
                )
    check_hashes(config["source_sha256"])
    # Re-run the independent structural→Amber→mature-local audit at release.
    command = [
        config["queue_python"],
        config["source_audit_script"],
        config["configuration_path"],
        str(output / "canonical_source_release_acceptance.json"),
    ]
    if native_run(
        root, command, output / "canonical_source_release", config["queue_environment"]
    ):
        raise RuntimeError("Independent canonical-source release audit failed")
    canonical = json.loads(
        (output / "canonical_source_release_acceptance.json").read_text()
    )
    if canonical["accepted"] is not True:
        raise ValueError("Canonical-source audit did not accept")
    rows = subprocess.check_output(
        [
            "nvidia-smi",
            "--id=0",
            "--query-gpu=memory.used,utilization.gpu",
            "--format=csv,noheader,nounits",
        ],
        text=True,
    )
    compute = subprocess.check_output(
        [
            "nvidia-smi",
            "--id=0",
            "--query-compute-apps=pid",
            "--format=csv,noheader,nounits",
        ],
        text=True,
    )
    validate_idle_gpu(rows, compute)
    write_new(
        output / "GPU_release_acceptance.json",
        {
            "accepted": True,
            "GPU0_memory_utilization": rows.strip(),
            "GPU0_compute_pids": compute.strip(),
            "existing_jobs_signaled": False,
        },
    )
    native = json.loads(Path(config["native_manifest"]).read_text())
    manifest = {
        "allocation_protocol": FOUR_TEMPLATE_PROTOCOL,
        "request": request,
        "native": native,
        "accepted_manifests": [gate],
        "local_pocket_mask": config["local_pocket_mask"],
        "source_sha256": {
            **config["source_sha256"],
            str(gate_path): digest(gate_path),
            str(terminal_path): digest(terminal_path),
            str(output / "canonical_source_release_acceptance.json"): digest(
                output / "canonical_source_release_acceptance.json"
            ),
        },
    }
    for parent in gate["parents"]:
        for field in (
            "pdb_path",
            "canonical_mapping",
            "pa6_graph_receipt",
            "scientific_acceptance",
            "authoritative_coordinates",
            "explicit_graph_cif",
        ):
            manifest["source_sha256"][parent[field]] = parent[field + "_sha256"]
    manifest_path = output / "broad_execution_manifest.json"
    write_new(manifest_path, manifest)
    write_new(
        output / "generation_execution_release.json",
        {
            "execution_release_checked": True,
            "qualification_gate_verified": True,
            "canonical_source_mapping_verified": True,
            "output_scientific_acceptance": False,
            "RFD_expansion_launched": False,
            "predictions_authorized": False,
            "manifest": str(manifest_path),
            "manifest_sha256": digest(manifest_path),
            "raw_budget": 24,
            "templates": 4,
            "samples_per_template": 6,
        },
    )
    journal(
        root,
        "intent",
        "Run generation-only four-template broad pilot",
        str(manifest_path),
    )
    code = native_run(
        root,
        [
            native["python"],
            config["broad_runner"],
            str(manifest_path),
            str(output / "generation_run01"),
        ],
        output / "broad_generation",
        native["environment"],
    )
    receipt_path = output / "generation_run01/scientific_acceptance.json"
    if code or not receipt_path.exists():
        raise RuntimeError(f"Native broad generation/crossscore failed: exit={code}")
    receipt = json.loads(receipt_path.read_text())
    if receipt["accepted"] is not True or receipt["raw_count"] != 24:
        raise ValueError("Broad generation scientific budget/tie audit failed")
    write_new(
        output / "generation_scientific_acceptance.json",
        {
            "accepted": True,
            "process_exit": code,
            "generation_scientific_acceptance": str(receipt_path),
            "generation_scientific_acceptance_sha256": digest(receipt_path),
            "variant_predictions_launched": False,
            "activity_validated": False,
        },
    )
    journal(
        root,
        "outcome",
        "Generation-only four-template broad pilot completed",
        "24 raw samples; full symmetry/protection and matched TS/ES crossscores accepted; predictions remain held",
    )


def main() -> None:
    """Run the finite queue once and preserve explicit terminal success/failure."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configuration", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    config = json.loads(args.configuration.read_text())
    try:
        execute(config)
    except Exception as error:
        journal(
            Path(config["repository_root"]),
            "outcome",
            "Generation-only queue terminated",
            repr(error),
        )
        write_new(
            Path(config["output"]) / "queue_terminal.json",
            {
                "terminal": True,
                "accepted": False,
                "error": repr(error),
                "variant_predictions_launched": False,
            },
        )
        raise
    write_new(
        Path(config["output"]) / "queue_terminal.json",
        {"terminal": True, "accepted": True, "variant_predictions_launched": False},
    )


if __name__ == "__main__":
    main()
