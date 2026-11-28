#!/usr/bin/env python3
"""Save preliminary final-frame geometry with fixed mapping and provenance."""

from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import parmed
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.design_md_geometry import (  # noqa: E402
    frame_rows,
    index_mapping,
    native_reference,
)


def sha(path: Path) -> str:
    """Return full file SHA-256."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()


def main() -> None:
    """Validate minimized mapping and save one completed production frame."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    p = args.base
    name = args.dataset
    manifest = p / "systems" / name / "datasets.json"
    d = json.loads(manifest.read_text())["datasets"][name]
    selection = p / "selection.json"
    source = Path(json.loads(selection.read_text())["source_selection"][name])
    top = (manifest.parent / d["topology"]["path"]).resolve(strict=True)
    rst = (manifest.parent / d["restart"]["path"]).resolve(strict=True)
    if sha(top) != d["topology"]["sha256"] or sha(rst) != d["restart"]["sha256"]:
        raise ValueError("Topology or restart hash mismatch")
    t = parmed.load_file(str(top), xyz=str(rst))
    mapping = index_mapping(t, source)
    ref = native_reference(args.reference)
    baseline = frame_rows(np.asarray(t.coordinates), np.asarray(t.box), mapping, ref)
    table = p.parent / "validation/analysis/occupied_sites.csv"
    before = [
        r
        for r in csv.DictReader(table.open())
        if r["path"] == str(source / "relaxed.pdb") and r["stage"] == "relaxed"
    ]
    if len(before) != 4:
        raise ValueError("Missing original four-site baseline")
    for i, r in enumerate(before):
        if (
            r["site"] != "AE BF CG DH".split()[i]
            or abs(baseline[i]["d1"] - float(r["d1"])) > 0.01
            or abs(baseline[i]["herschlag_phi"] - float(r["herschlag_phi"])) > 0.15
            or baseline[i]["tier"] != r["tier"]
            or baseline[i]["native_pass"] != (r["native_pass"] == "True")
        ):
            raise ValueError("Initial mapped geometry mismatch")
    trajectory = p / "md" / name / "replica_01/dynamics/06_production_90ns.nc"
    receipt = p / "md" / name / "replica_01/100ns_process_complete.json"
    if not json.loads(receipt.read_text())["completed"]:
        raise ValueError("Replica incomplete")
    with netcdf_file(trajectory, "r", mmap=True) as nc:
        count = nc.variables["coordinates"].shape[0]
        xyz = np.array(nc.variables["coordinates"][-1])
        box = np.r_[nc.variables["cell_lengths"][-1], nc.variables["cell_angles"][-1]]
        time = float(nc.variables["time"][-1])
    if xyz.shape != (len(t.atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Invalid final frame")
    rows = frame_rows(xyz, box, mapping, ref)
    files = [
        top,
        rst,
        manifest,
        selection,
        table,
        args.reference,
        receipt,
        *[source / "system/components" / c / "conversion.json" for c in "IJKL"],
    ]
    record = {
        "dataset": name,
        "replica": "replica_01",
        "frame_index_zero_based": count - 1,
        "frame_count": count,
        "trajectory": str(trajectory),
        "trajectory_size_bytes": trajectory.stat().st_size,
        "frame_coordinate_sha256": hashlib.sha256(xyz.tobytes()).hexdigest(),
        "frame_box": box.tolist(),
        "time_ps": time,
        "baseline_mapping_verified": True,
        "interpretation": "Single terminal frame; four sites nested in one replica; no occupancy or activity claim",
        "inputs_sha256": {str(x): sha(x) for x in files},
        "sites": [
            {"site": site, **row}
            for site, row in zip(("AE", "BF", "CG", "DH"), rows, strict=True)
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    with args.output.with_suffix(".csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(record["sites"][0]))
        w.writeheader()
        w.writerows(record["sites"])
    print(
        json.dumps(
            {
                "dataset": name,
                "time_ps": time,
                "output": str(args.output),
                "sites": [
                    (
                        r["site"],
                        r["tier"],
                        r["native_pass"],
                        round(r["herschlag_phi"], 1),
                        round(r["d1"], 2),
                    )
                    for r in record["sites"]
                ],
            }
        )
    )


if __name__ == "__main__":
    main()
