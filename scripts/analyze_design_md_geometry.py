#!/usr/bin/env python3
"""Stream fixed-assignment four-site design geometry on completed MD replicas."""

from __future__ import annotations
import argparse
import csv
import json
import logging
from pathlib import Path
import sys
import numpy as np
import parmed
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.design_md_completion import (  # noqa: E402
    sha,
    validate_completion,
    add_selection_eligibility,
    replica_path,
)
from nyl12_covalent.design_md_geometry import (  # noqa: E402
    frame_rows,
    index_mapping,
    native_reference,
)
from nyl12_covalent.native_thr_policy import policy_metadata  # noqa: E402


def main() -> None:
    """Analyze every tenth production frame with original assignment and gates."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset",
        required=True,
        choices=("WT", "broader_mask07", "focused_mask03", "V111C"),
    )
    parser.add_argument(
        "--replica",
        choices=("01", "02"),
        default="01",
        help="Replica number (default: 01).",
    )
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = args.project_root / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison"
    follow = base / "md_followup"
    name = args.dataset
    manifest = follow / "systems" / name / "datasets.json"
    data = json.loads(manifest.read_text())["datasets"][name]
    topology = (manifest.parent / data["topology"]["path"]).resolve(strict=True)
    restart = (manifest.parent / data["restart"]["path"]).resolve(strict=True)
    if (
        sha(topology) != data["topology"]["sha256"]
        or sha(restart) != data["restart"]["sha256"]
    ):
        raise ValueError("Dataset input hash mismatch")
    selection = follow / "selection.json"
    source = Path(json.loads(selection.read_text())["source_selection"][name])
    reference = args.project_root / "cif/n12c.cif"
    parm = parmed.load_file(str(topology), xyz=str(restart))
    mapping = index_mapping(parm, source)
    ref = native_reference(reference)
    initial = frame_rows(
        np.asarray(parm.coordinates), np.asarray(parm.box), mapping, ref
    )
    table = base / "validation/analysis/occupied_sites.csv"
    initial_rows = [
        r
        for r in csv.DictReader(table.open())
        if r["path"] == str(source / "relaxed.pdb") and r["stage"] == "relaxed"
    ]
    if len(initial_rows) != 4:
        raise ValueError("Missing original four-site comparison")
    for i, r in enumerate(initial_rows):
        if (
            r["site"] != "AE BF CG DH".split()[i]
            or abs(initial[i]["d1"] - float(r["d1"])) > 0.01
            or abs(initial[i]["herschlag_phi"] - float(r["herschlag_phi"])) > 0.15
            or initial[i]["tier"] != r["tier"]
            or initial[i]["native_pass"] != (r["native_pass"] == "True")
        ):
            raise ValueError(f"Initial mapping mismatch for {r['site']}")
    replica_dir = replica_path(follow, name, args.replica)
    trajectories = validate_completion(
        replica_dir, data["topology"]["sha256"], data["restart"]["sha256"]
    )
    args.output.mkdir(parents=True, exist_ok=False)
    rows = []
    for trajectory in trajectories:
        with netcdf_file(trajectory, "r", mmap=True) as nc:
            n = nc.variables["coordinates"].shape[0]
            expected_frames = 1000 if trajectory == trajectories[0] else 9000
            if n != expected_frames:
                raise ValueError("Unexpected production trajectory frame count")
            times = np.asarray(nc.variables["time"][:], dtype=float)
            first_time = 510.0 if trajectory == trajectories[0] else 10510.0
            if (
                times.shape != (n,)
                or not np.isfinite(times).all()
                or not np.allclose(
                    times, first_time + 10.0 * np.arange(n), atol=0.01, rtol=0
                )
            ):
                raise ValueError("Unexpected production trajectory timing")
            for i in range(9, n, 10):
                xyz = np.array(nc.variables["coordinates"][i])
                box = np.r_[
                    nc.variables["cell_lengths"][i], nc.variables["cell_angles"][i]
                ]
                if xyz.shape != (len(parm.atoms), 3) or not np.isfinite(xyz).all():
                    raise ValueError("Invalid frame")
                geometry = frame_rows(xyz, box, mapping, ref)
                time = float(nc.variables["time"][i])
                for site, row in zip(("AE", "BF", "CG", "DH"), geometry, strict=True):
                    rows.append(
                        {
                            "time_ps": time,
                            "site": site,
                            **add_selection_eligibility(row),
                        }
                    )
    if len(rows) != 4000:
        raise ValueError(
            f"Expected 1000 sampled frames and four sites, got {len(rows)} records"
        )
    with (args.output / "sites.jsonl").open("w") as stream:
        for row in rows:
            stream.write(json.dumps(row, allow_nan=False) + "\n")
    with (args.output / "sites.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    results = []
    for site in ("AE", "BF", "CG", "DH"):
        r = [x for x in rows if x["site"] == site]
        for block, part in [
            ("all", r),
            ("first_half", r[: len(r) // 2]),
            ("second_half", r[len(r) // 2 :]),
        ]:
            results.append(
                {
                    "site": site,
                    "block": block,
                    "frames": len(part),
                    "tier_a": sum(x["tier"] == "A" for x in part),
                    "tier_ab": sum(x["tier"] in ("A", "B") for x in part),
                    "selection_eligible": sum(x["selection_eligible"] for x in part),
                    "native": sum(x["native_pass"] for x in part),
                    "joint": sum(x["joint_pass"] for x in part),
                    "negative_phi": sum(x["herschlag_phi"] < 0 for x in part),
                    "median_d1_A": float(np.median([x["d1"] for x in part])),
                }
            )
    inputs = [
        topology,
        restart,
        reference,
        manifest,
        selection,
        table,
        replica_dir / "100ns_process_complete.json",
        replica_dir / "dynamics/05_production_10ns.process.json",
        replica_dir / "dynamics/05_production_10ns.acceptance.json",
        replica_dir / "dynamics/06_production_90ns.acceptance.json",
        replica_dir / "dynamics/extension_100ns.json",
        *trajectories,
        *[source / "system/components" / c / "conversion.json" for c in "IJKL"],
    ]
    report = {
        "dataset": name,
        "replica": f"replica_{args.replica}",
        "native_thr_role": "diagnostic_only",
        "native_thr_policy": policy_metadata(),
        "selection_policy": "Tier A/B catalytic geometry only; native and joint metrics retained as historical diagnostics",
        "sampling": "Every tenth 10-ps production frame; 100-ps stride, two halves; four sites nested in one replica.",
        "interpretation": "Preliminary occupancy, not independent replicates or sampling convergence",
        "initial_geometry_validated": True,
        "input_sha256": {str(f): sha(f) for f in inputs},
        "summary": results,
    }
    (args.output / "summary.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    logging.warning("Finished %s frames for %s", len(rows) // 4, name)


if __name__ == "__main__":
    main()
