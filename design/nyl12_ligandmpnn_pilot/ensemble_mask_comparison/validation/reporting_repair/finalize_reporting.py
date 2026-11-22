"""Complete reporting from existing accepted prediction/minimization artifacts."""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

from nyl12_covalent.local_design_validation_inputs import atomic_json
from nyl12_covalent.local_design_validation_phases import _compare

LOGGER = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Return SHA256 of a local file."""
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main() -> None:
    """Republish comparison and complete the existing campaign without calculations."""
    repair = Path(__file__).resolve().parent
    campaign = repair.parent
    parent = campaign.parent
    manifest = json.loads((campaign / "prepared/candidate_manifest.json").read_text())
    jobs = manifest["jobs"]
    reused = manifest["reused_models"]
    config = json.loads((campaign / "configuration/relaxation.json").read_text())
    prior = json.loads((campaign / "status.json").read_text())
    if prior["phase"] != "comparison" or prior["status"] != "failed":
        raise ValueError("Expected retained comparison-only failure")
    if (len(jobs), len(reused), prior["relaxed_models"]) != (52, 3, 52):
        raise ValueError("Unexpected accepted model inventory")
    result = _compare(campaign, jobs, reused, config)
    unchanged = ("model_summary.csv", "occupied_sites.csv", "all_sites.csv", "report.json")
    for name in unchanged:
        if digest(campaign / "analysis" / name) != digest(repair / "original" / name):
            raise ValueError(f"Previously computed geometry changed: {name}")
    accepted = result["acceptance_rows"]
    if len(accepted) != 55 or not all(row["process_completed"] for row in accepted):
        raise ValueError("Expected55 process-completed model seeds")
    completion = {
        "status": "completed",
        "new_models": len(jobs),
        "reused_models": len(reused),
        "total_model_seeds": len(accepted),
        "predicted_models": prior["predicted_models"],
        "normalized_models": prior["normalized_models"],
        "relaxed_models": prior["relaxed_models"],
        "process_completed_models": sum(row["process_completed"] for row in accepted),
        "force_converged_models": sum(row["force_converged"] for row in accepted),
        "sample_unit": "Model seed; four catalytic sites are nested within each model.",
        "interpretation": "Geometry and minimization do not validate catalytic activity.",
        "geometry_report": result["geometry_report"],
        "plots": result["plots"],
        "reporting_repair": str(repair),
        "scientific_calculations_repeated": False,
        "original_workflow_exit_status": 1,
    }
    atomic_json(campaign / "COMPLETED.json", completion)
    prior = {key: value for key, value in prior.items() if key not in {"error", "error_type"}}
    prior.update(status="completed", phase="complete", reporting_repair=str(repair))
    atomic_json(campaign / "status.json", prior)
    atomic_json(parent / "status.json", {
        "status": "completed_after_reporting_repair",
        "completion": str(campaign / "COMPLETED.json"),
        "original_workflow_exit_status": 1,
        "scientific_calculations_repeated": False,
    })
    LOGGER.info("Completed comparison: %s model seeds; %s force-converged", len(accepted), completion["force_converged_models"])


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
