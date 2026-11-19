"""Record explicit startup precision review without rewriting strict evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    """Return SHA256 for retained reproducibility evidence."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def review(trial: Path) -> dict:
    """Apply explicitly authorized1e-6A startup precision review to raw evidence.

    Args:
        trial: Existing verified trial with immutable strict startup report.

    Returns:
        Separate review record; raw1e-7report and engine inputs remain unchanged.
    """
    raw_path = trial / "startup_verification.json"
    raw = json.loads(raw_path.read_text())
    cfg = json.loads((trial / "config.json").read_text())
    manifest_path = Path(cfg["source_seed_manifest"])
    manifest = json.loads(manifest_path.read_text())
    if cfg["images"] != manifest["images"]:
        raise ValueError("Trialcontrols differ from frozen manifest")
    for name in cfg["images"]:
        if (
            raw["input_sha256"][name] != manifest["image_sha256"][name]
            or digest(Path(name)) != manifest["image_sha256"][name]
        ):
            raise ValueError("Frozen manifest/control hash mismatch")
    errors = [
        raw[k]
        for k in (
            "max_component_error_A",
            "max_vector_error_A",
            "qm_max_vector_error_A",
            "H_max_vector_error_A",
        )
    ]
    errors += raw["endpoint_vector_errors_A"]
    accepted = bool(
        max(errors) <= 1e-6
        and raw["element_order_verified"]
        and raw["input_frozen_exterior_bitwise"]
        and raw["periodic_box_bitwise"]
    )
    result = dict(
        status="accepted_startup_geometry_precision_review"
        if accepted
        else "failed_precision_review",
        accepted=accepted,
        raw_report=str(raw_path),
        raw_report_sha256=digest(raw_path),
        original_strict_passed=raw["passed"],
        original_component_tolerance_A=raw["coordinate_component_tolerance_A"],
        reviewed_coordinate_tolerance_A=1e-6,
        reviewed_endpoint_tolerance_A=1e-6,
        max_component_error_A=raw["max_component_error_A"],
        max_vector_error_A=raw["max_vector_error_A"],
        qm_max_vector_error_A=raw["qm_max_vector_error_A"],
        H_max_vector_error_A=raw["H_max_vector_error_A"],
        endpoint_vector_errors_A=raw["endpoint_vector_errors_A"],
        frozen_manifest_path=str(manifest_path),
        frozen_manifest_sha256=digest(manifest_path),
        control_hashes_match_frozen_manifest=True,
        rationale="Explicit rootreview: native write_xyz uses f15.7 and computed ANG_AU; Sander inputs use CODATA08_A_TO_BOHRS. Observed~2e-7A roundtrip/alignment/formatting difference is below declared1e-6A startup precision, not a path/atom-label failure.",
        source_evidence=dict(
            native_writer="dlf_util.F90:806-855 (f15.7; ANG_AU)",
            native_input="sander/dlfind_module.F90:17,2189-2195 (CODATA08_A_TO_BOHRS)",
        ),
        review_script_sha256=digest(Path(__file__)),
        coordinates_modified=False,
        raw_report_preserved=True,
        engines_modified_or_stopped=False,
        force_convergence_assessed=False,
        chemical_stationarity_assessed=False,
    )
    with (trial / "startup_precision_review.json").open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trials", type=Path, nargs="+")
    for path in parser.parse_args().trials:
        result = review(path.resolve())
        print(path.name, result["accepted"], result["max_component_error_A"])
