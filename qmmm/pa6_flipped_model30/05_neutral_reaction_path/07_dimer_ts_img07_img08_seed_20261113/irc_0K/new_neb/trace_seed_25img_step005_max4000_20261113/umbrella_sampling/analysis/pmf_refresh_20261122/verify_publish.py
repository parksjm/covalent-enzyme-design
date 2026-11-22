"""Verify and publish the requested PMF refresh without altering sampling."""
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import sys
import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT / "source_snapshot/scripts"))
from requested_path_pmf import KB, HarmonicWindow, collect, solve
from primary_q_acyl import density_profile, regression


def main() -> None:
    """Check immutable inputs and independent marginal, then switch compact pointer."""
    status = json.loads((OUT / "job/status.json").read_text())
    if status["exit_status"] != 0:
        raise RuntimeError("PMF calculation has not completed successfully")
    manifest = json.loads((OUT / "input_source_manifest.json").read_text())
    for relative, item in manifest["files"].items():
        if "__pycache__" in Path(relative).parts and Path(relative).suffix == ".pyc":
            continue  # Python regenerates derived bytecode on import; .py sources remain checked.
        if hashlib.sha256((OUT / relative).read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Frozen input/source changed: {relative}")
    regression()
    cfg, windows = collect(OUT / "frozen_inputs/prepared.json")
    result = json.loads((OUT / "mbar/latest.json").read_text())
    summary = json.loads((OUT / "primary_q_acyl_summary.json").read_text())
    old = Path(manifest["previous_analysis"])
    previous = json.loads((old / "primary_q_acyl_summary.json").read_text())
    np.testing.assert_allclose(summary["reference_bin_edges_A"], previous["reference_bin_edges_A"], atol=1e-12, rtol=0)
    raw = [np.concatenate([b["samples"] for b in w["blocks"]]) for w in windows]
    selected = [a[::s["decorrelation_stride"]] for a, s in zip(raw, result["window_summaries"], strict=True)]
    beta = 1 / (KB * cfg["temperature_K"])
    weights, _, connected = solve(selected, [HarmonicWindow(**w["bias"]) for w in cfg["windows"]], beta)
    if not connected or summary["bootstrap_connected"] != 50:
        raise ValueError("Disconnected estimate or incomplete block bootstrap")
    profile = np.load(OUT / "primary_q_acyl_bootstrap.npz")
    edges = profile["edges"]
    points = np.concatenate(selected)
    reference = int(np.searchsorted(edges, cfg["windows"][0]["q_acyl_center"], side="right") - 1)
    pmf = density_profile(points[:, 0], weights, edges, beta, reference)
    np.testing.assert_allclose(pmf, profile["pmf"], atol=1e-10, rtol=1e-12)
    h2, _, _ = np.histogram2d(points[:, 0], points[:, 1], bins=[edges, 23], weights=weights)
    direct = np.bincount(np.minimum(np.searchsorted(edges, points[:, 0], side="right") - 1, len(edges)-2), weights=weights, minlength=len(edges)-1)
    np.testing.assert_allclose(h2.sum(axis=1), direct, atol=0, rtol=1e-12)
    counts = Counter(r["hg1_owner"] for w in windows for b in w["blocks"] for r in b["rows"])
    retention = {"raw_biased_frame_counts_hg1_owner_atom1based": dict(counts), "roles_1based": cfg["roles_1based"], "unbiased_populations": False, "chemical_frames_censored": 0}
    (OUT / "proton_state_retention.json").write_text(json.dumps(retention, indent=2)+"\n")
    with (OUT / "per_window_accepted_coverage.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["window", "accepted_production_ps", "completed_5ps_blocks"])
        writer.writerows((w["label"], 5*len(w["blocks"]), len(w["blocks"])) for w in windows)
    export = OUT / "compact_export"
    export.mkdir()
    names = ["README.md", "primary_q_acyl_pmf.png", "primary_q_acyl_pmf.svg", "primary_q_acyl_pmf.csv", "primary_q_acyl_summary.json", "primary_scientific_acceptance.json", "proton_state_retention.json", "per_window_accepted_coverage.csv"]
    for name in names:
        shutil.copy2(OUT / name, export / name)
    for name in ["primary_fes.png", "diagnostic_fes.png", "mbar_overlap.png", "histogram_overlap.png"]:
        shutil.copy2(Path(result["generation"]) / name, export / name)
    shutil.copy2(OUT / "job/status.json", export / "execution_status.json")
    shutil.copy2(OUT / "job/command.json", export / "execution_command.json")
    pointer = OUT.parent / "pmf_latest"
    old_pointer = str(pointer.readlink())
    receipt = {"published_utc": datetime.now(timezone.utc).isoformat(), "analysis": str(OUT), "previous_pointer": old_pointer, "input_source_manifest_sha256": hashlib.sha256((OUT/"input_source_manifest.json").read_bytes()).hexdigest(), "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "checks": ["process_exit_zero", "all_frozen_source_and_input_hashes", "same_ES_density_reference_bin", "fifty_decade_regression", "independent_2D_marginal", "50_connected_block_bootstraps"], "software": {n: importlib.metadata.version(n) for n in ["numpy", "scipy", "pymbar", "matplotlib"]}, "artifact_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in export.iterdir() if f.is_file()}}
    (export / "provenance.json").write_text(json.dumps(receipt, indent=2)+"\n")
    temporary = OUT.parent / "pmf_latest.pending_refresh"
    temporary.symlink_to(export.relative_to(OUT.parent))
    temporary.replace(pointer)
    (OUT / "publication_receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    import logging
    logging.basicConfig(level=logging.INFO)
    logging.info("Verified and published updated PMF: %.6f kcal/mol", summary["maximum_relative_ES_bin_kcal_mol"])


if __name__ == "__main__":
    main()
