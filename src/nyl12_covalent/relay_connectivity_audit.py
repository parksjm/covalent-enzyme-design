"""Read-only audit orchestration for accepted N-terminal relay path records."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .relay_connectivity_records import (
    _accepted_rows,
    _qm_atom_ids,
)
from .relay_connectivity_source import (
    BRANCHES,
    DISCONTINUITY_A,
    SnapshotChangedError,
    _read_snapshot,
    _verify_source,
    digest,
    read_json,
)


def audit_campaign(campaign: Path) -> dict[str, Any]:
    """Audit accepted relay records without launching calculations.

    Args:
        campaign: Existing 12_nterminal_relay campaign directory.

    Returns:
        JSON-safe provisional connectivity audit and receipt references.
    """
    campaign = campaign.resolve(strict=True)
    config = read_json(campaign / "config.json")
    source = _verify_source(campaign, config)
    topology = campaign / config["model"] / "inputs/system.prmtop"
    topology_hash = digest(topology)
    qm_atom_ids = _qm_atom_ids(
        topology, campaign / config["model"] / "inputs/dimer.mdin"
    )
    branch_reports = {}
    common_inputs: dict[str, dict[str, str]] = {}
    live_records: list[tuple[Path, str]] = []
    for directory, chemistry in BRANCHES.items():
        branch = campaign / directory
        progress_path = branch / "progress.json"
        progress, progress_hash = _read_snapshot(progress_path)
        live_records.append((progress_path, progress_hash))
        status_path = branch / "branch_status.json"
        status, status_hash = _read_snapshot(status_path)
        live_records.append((status_path, status_hash))
        rows = progress.get("rows")
        if not isinstance(rows, list):
            raise ValueError(f"Missing path rows: {progress_path}")
        accepted, receipts = _accepted_rows(
            campaign, branch, rows, topology_hash, qm_atom_ids, chemistry
        )
        live_records.extend(
            (Path(receipt["process_path"]), receipt["process_sha256"])
            for receipt in receipts
        )
        live_records.extend(
            (Path(receipt["acceptance_path"]), receipt["acceptance_sha256"])
            for receipt in receipts
            if receipt["acceptance_path"] and receipt["acceptance_sha256"]
        )
        for entry in accepted:
            if entry["process_receipt"] is None:
                continue
            process = read_json(Path(entry["process_receipt"]["process_path"]))
            inputs = process["input_sha256"]
            shared = {
                Path(path).name: value
                for path, value in inputs.items()
                if Path(path).name in {"system.prmtop", "sander", "libdlfind.so"}
            }
            common_inputs.setdefault(directory, {}).update(shared)
        discontinuities = []
        previous_geometry = None
        for row in rows:
            if row.get("accepted") is not True:
                continue
            current = row.get("geometry", {}).get("distances_A", {})
            if previous_geometry is not None:
                for key in sorted(set(current) & set(previous_geometry)):
                    delta = float(current[key]) - float(previous_geometry[key])
                    if abs(delta) > DISCONTINUITY_A:
                        discontinuities.append(
                            {
                                "label": row.get("label"),
                                "distance": key,
                                "signed_change_A": delta,
                            }
                        )
            previous_geometry = current
        branch_reports[chemistry] = {
            "branch_directory": str(branch),
            "status": status.get("status", "unknown"),
            "status_receipt": {"path": str(status_path), "sha256": status_hash},
            "progress_receipt": {"path": str(progress_path), "sha256": progress_hash},
            "accepted_path_rows": accepted,
            "process_receipts": receipts,
            "distance_discontinuities_over_0_35_A": discontinuities,
            "torsion_discontinuity_status": "not available in accepted path records",
        }
    shared_keys = (
        set.intersection(
            *(set(common_inputs.get(directory, {})) for directory in BRANCHES)
        )
        if all(directory in common_inputs for directory in BRANCHES)
        else set()
    )
    shared_identity = {
        key: {directory: common_inputs[directory][key] for directory in BRANCHES}
        for key in sorted(shared_keys)
    }
    if any(len(set(values.values())) != 1 for values in shared_identity.values()):
        raise ValueError("Branch topology or shared native Hamiltonian inputs differ")
    for path, expected in live_records:
        if digest(path) != expected:
            raise SnapshotChangedError(
                f"Live campaign record changed during audit: {path}"
            )
    endpoints_accepted = all(
        directory in common_inputs
        and any(
            row["label"] == "unrestrained_endpoint"
            and row["unbiased_endpoint_candidate_passed"] is True
            for row in report["accepted_path_rows"]
        )
        for report in branch_reports.values()
    )
    return {
        "schema_version": 1,
        "audit_status": "complete" if endpoints_accepted else "provisional_incomplete",
        "connectivity_claim": "not established",
        "required_next_gate": "accepted unbiased endpoints plus exact atom-inventory/reference comparison to E-S and acylenzyme records",
        "campaign": str(campaign),
        "atom_inventory": {"atom_count": 104922, "topology_sha256": topology_hash},
        "source_manifest": source,
        "shared_topology_exterior_and_physical_hamiltonian": {
            "verified_shared_input_hashes": shared_identity,
            "frozen_base_qmmm_input_sha256": source["base_hamiltonian_input_sha256"],
            "exterior_gate": "each accepted window reports exactly zero frozen-coordinate change from the common source",
            "biased_hamiltonians": "branch-specific reaction-coordinate biases differ by design; bias and physical energies are separate",
        },
        "branches": branch_reports,
        "interpretation_limits": [
            "A force-accepted path frame is not a chemically accepted endpoint or verified transition state.",
            "No E-S to relay to acylenzyme connectivity is claimed until unbiased endpoint records pass and exact atom inventories/reference states match.",
            "No unconditional barrier or PMF interpretation is supported until a relay-ready E-S state is related to the dominant bound reactant ensemble and its preorganization free energy.",
            "No torsions are present in the accepted path records; only saved distance changes are assessed.",
        ],
    }
