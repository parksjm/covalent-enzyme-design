"""Publish a new two-original-parent gate after actual native and science completion."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import subprocess
import time

LOGGER = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Return the SHA256 of an existing retained artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: object) -> None:
    """Create a new receipt without overwriting historical failures."""
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def main() -> None:
    """Wait for the one new job, verify artifact hashes and publish a new gate."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    manifest = json.loads(args.manifest.read_text())
    root = Path(manifest["campaign_root"])
    process = root / "job_parent03_20261127/status.json"
    while True:
        state = json.loads(process.read_text())
        if state["ended_utc"] is not None:
            break
        LOGGER.info("Waiting for independent parent03 native/export completion")
        time.sleep(60)
    action = [
        "/usr/bin/python3",
        manifest["journal_script"],
        "--log",
        manifest["journal"],
        "--agent",
        "rfd3_export_and_parent03",
        "--permission-mode",
        "escalated",
        "--reason",
        "Publish actual two-original-parent qualification status",
        "--rollback",
        "Preserve historical failed receipts and native outputs",
    ]
    subprocess.run(
        action
        + [
            "--phase",
            "intent",
            "--action",
            "Verify and publish new two-parent accepted manifest after parent03 terminal status",
        ],
        check=True,
    )
    parents = []
    candidates = [
        root / "exports_adapter_20261127_parent01_t5",
        Path(manifest["derived_export"]),
    ]
    for directory in candidates:
        receipt_path = directory / "scientific_acceptance.json"
        if not receipt_path.exists():
            continue
        data = json.loads(receipt_path.read_text())
        if not data["accepted"]:
            continue
        for field in (
            "pdb_path",
            "canonical_mapping",
            "pa6_graph_receipt",
            "authoritative_coordinates",
            "explicit_graph_cif",
        ):
            if digest(Path(data[field])) != data[field + "_sha256"]:
                raise ValueError(f"Accepted export changed: {data[field]}")
        if digest(Path(data["native_cif"])) != data["native_cif_sha256"]:
            raise ValueError("Raw native output changed")
        parents.append(
            {
                "parent_id": data["job_id"],
                "sequence_scope": "native_WT",
                "parent_source_id": data["parent_source_id"],
                "noise_A": 5,
                "seed": next(
                    j["seed"]
                    for j in json.loads(Path(data["manifest"]).read_text())["jobs"]
                    if j["job_id"] == data["job_id"]
                ),
                **{
                    k: data[k]
                    for k in (
                        "pdb_path",
                        "canonical_mapping",
                        "pa6_graph_receipt",
                        "pdb_path_sha256",
                        "canonical_mapping_sha256",
                        "pa6_graph_receipt_sha256",
                        "authoritative_coordinates",
                        "authoritative_coordinates_sha256",
                        "explicit_graph_cif",
                        "explicit_graph_cif_sha256",
                    )
                },
                "scientific_acceptance": str(receipt_path),
                "scientific_acceptance_sha256": digest(receipt_path),
            }
        )
    accepted = (
        state["exit_status"] == 0
        and {p["parent_source_id"] for p in parents} == {"medoid01", "medoid03"}
        and len(parents) == 2
    )
    status = "accepted" if accepted else "scientific_failed"
    result = {
        "accepted": accepted,
        "status": status,
        "operational_status": str(process),
        "parents": parents,
        "required_distinct_original_parents": ["medoid01", "medoid03"],
        "historical_failures_retained": True,
        "rejected_parent01_t2": str(
            root
            / "exports_adapter_20261127_parent01_t2_rejected/scientific_acceptance.json"
        ),
        "generation_launched": False,
    }
    save(root / "accepted_parent_manifest_adapter_20261127.json", result)
    save(
        root / "qualification_terminal_adapter_20261127.json",
        {
            "terminal": True,
            "accepted": accepted,
            "status": status,
            "accepted_parent_manifest": str(
                root / "accepted_parent_manifest_adapter_20261127.json"
            ),
        },
    )
    subprocess.run(
        action
        + [
            "--phase",
            "outcome",
            "--action",
            "Published new two-original-parent gate",
            "--result",
            f"accepted={accepted}; native process exit={state['exit_status']}; broader controller not launched",
        ],
        check=True,
    )
    LOGGER.info("Two-original-parent gate accepted=%s", accepted)
    if not accepted:
        raise RuntimeError("Two-original-parent qualification remains held")


if __name__ == "__main__":
    main()
