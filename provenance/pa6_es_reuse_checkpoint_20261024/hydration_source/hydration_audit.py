"""Read-only periodic watershell counts for existing WT PA6 trajectories."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
MANIFEST = PROJECT / "configs/pa6_wt_es_reuse_20261024.json"
PRIOR = (
    PROJECT / "output_nyl12_pa6_es_audit_20261024/reuse_audit_attempt1/reuse_audit.json"
)
INSPECTION = (
    PROJECT
    / "output_nyl12_pa6_es_audit_20261024/seed_inspection_attempt4/inspection.json"
)
LOGGER = logging.getLogger(__name__)


def _hash(path: Path) -> dict[str, Any]:
    before = path.stat()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Source changed during hashing: {path}")
    return {
        "path": str(path),
        "bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "sha256": digest,
    }


def _check_inherited(record: dict[str, Any]) -> dict[str, Any]:
    current = Path(record["path"]).stat()
    if (current.st_size, current.st_mtime_ns) != (record["bytes"], record["mtime_ns"]):
        raise ValueError(f"Original source metadata changed: {record['path']}")
    return {
        **record,
        "integrity_policy": "Inherited SHA256 from original audit; current size and mtime verified; not rehashed",
    }


def _counts(path: Path) -> list[tuple[int, int, int]]:
    rows = []
    headers = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if line.lstrip().startswith("#"):
            headers.append(line)
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"Expected exactly frame/lower/upper columns: {path}")
        frame, lower, upper = map(int, fields)
        if frame != len(rows) + 1 or not 0 <= lower <= upper:
            raise ValueError(f"Invalid frame/count values: {path}: {line}")
        rows.append((frame, lower, upper))
    if len(rows) != 2000 or not headers:
        raise ValueError(f"Expected header and 2000 saved frames: {path}")
    LOGGER.info("Validated %s headers=%s", path.name, headers)
    return rows


def run() -> dict[str, Any]:
    """Run sequential single-thread cpptraj analyses and validate all counts."""
    config = json.loads(MANIFEST.read_text())
    prior = json.loads(PRIOR.read_text())
    inspection = json.loads(INSPECTION.read_text())
    system, amberhome = Path(config["system"]), Path(config["amberhome"])
    executable = amberhome / "bin/cpptraj"
    inherited = {row["path"]: row for row in prior["sources"]}
    topology = system / "md_solvated.prmtop"
    paths = [topology] + [
        system / replica / "07_prod_100ns.nc" for replica in config["replicas"]
    ]
    original_records = [_check_inherited(inherited[str(path)]) for path in paths]
    manifest_record = _hash(MANIFEST)
    if manifest_record["sha256"] != inherited[str(MANIFEST)]["sha256"]:
        raise ValueError("Manifest differs from original reuse audit")
    records = [
        manifest_record,
        _hash(PRIOR),
        _hash(INSPECTION),
        _hash(executable),
        _hash(Path(__file__)),
    ]
    for path in (MANIFEST, PRIOR, INSPECTION):
        shutil.copy2(path, ROOT / "source" / path.name)
    (ROOT / "source_manifest.json").write_text(
        json.dumps(
            {
                "hashed_sources": records,
                "inherited_sources": original_records,
            },
            indent=2,
        )
        + "\n"
    )
    selected = {
        (frame["source"]["replica"], frame["source"]["full_frame"]): frame[
            "water_inventory"
        ]["local_water_count"]
        for frame in inspection["frames"]
    }
    if len(selected) != 9 or any(value != 0 for value in selected.values()):
        raise ValueError("Expected exactly nine accepted seeds with zero local waters")
    production = {Path(row["directory"]).name: row for row in prior["production"]}
    reports = {}
    selected_checks = []
    env = os.environ.copy()
    env.update(
        AMBERHOME=str(amberhome),
        OMP_NUM_THREADS="1",
        LD_LIBRARY_PATH=str(amberhome / "lib"),
    )
    for replica in config["replicas"]:
        started = time.time()
        if production[replica]["frame_dt_ns"] != 0.05:
            raise ValueError("Expected full saved-frame spacing of 0.05 ns")
        directory = ROOT / "attempts" / replica
        directory.mkdir(parents=True, exist_ok=False)
        trajectory = system / replica / "07_prod_100ns.nc"
        script = directory / "watershell.cpptraj"
        script.write_text(
            f"parm {topology}\ntrajin {trajectory} 1 last\n"
            f"watershell :929@N out {directory / 'n35.dat'} lower 3.5 upper 5.0 :WAT@O\n"
            f"watershell :929@OG1 out {directory / 'og35.dat'} lower 3.5 upper 5.0 :WAT@O\n"
            f"watershell :929@N,OG1 out {directory / 'union35.dat'} lower 3.5 upper 5.0 :WAT@O\n"
            f"watershell :929@N,OG1 out {directory / 'union4.dat'} lower 4.0 upper 5.0 :WAT@O\n"
            "run\n"
        )
        LOGGER.info(
            "Starting %s: 2000 full frames, four periodic watershell actions", replica
        )
        with (directory / "cpptraj.log").open("w") as log:
            process = subprocess.run(
                [str(executable), "-i", str(script)],
                cwd=directory,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        status = {
            "exit_status": process.returncode,
            "elapsed_seconds": time.time() - started,
            "command": [str(executable), "-i", str(script)],
        }
        (directory / "process_status.json").write_text(
            json.dumps(status, indent=2) + "\n"
        )
        if process.returncode != 0:
            raise RuntimeError(f"CPPTRAJ failed for {replica}: {process.returncode}")
        data = {
            name: _counts(directory / f"{name}.dat")
            for name in ("n35", "og35", "union35", "union4")
        }
        both_series, any4_series = [], []
        with (directory / "hydration_counts.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "full_frame",
                    "physical_time_ns",
                    "n_lower35",
                    "og_lower35",
                    "union_lower35",
                    "both_lower35",
                    "union_lower4",
                ]
            )
            for index in range(2000):
                frame = index + 1
                n, og, union, union4 = (
                    data[name][index][1]
                    for name in ("n35", "og35", "union35", "union4")
                )
                both = n + og - union
                if not 0 <= both <= min(n, og) or not union <= union4:
                    raise ValueError(
                        f"Invalid inclusion-exclusion counts: {replica} frame {frame}"
                    )
                writer.writerow([frame, frame * 0.05, n, og, union, both, union4])
                both_series.append(both)
                any4_series.append(union4)
                if (replica, frame) in selected:
                    expected = selected[(replica, frame)]
                    selected_checks.append(
                        {
                            "replica": replica,
                            "full_frame": frame,
                            "periodic_union4_count": union4,
                            "accepted_inspection_local_count": expected,
                            "passed": union4 == expected,
                        }
                    )
                    if union4 != expected:
                        raise ValueError(
                            f"Selected seed count mismatch: {replica} {frame}"
                        )
        reports[replica] = {
            "saved_frames": 2000,
            "physical_time_ns": [0.05, 100.0],
            "any_water_within4_frame_count": sum(count > 0 for count in any4_series),
            "any_water_within4_frame_fraction": sum(count > 0 for count in any4_series)
            / 2000,
            "same_water_within35_of_both_frame_count": sum(
                count > 0 for count in both_series
            ),
            "same_water_within35_of_both_frame_fraction": sum(
                count > 0 for count in both_series
            )
            / 2000,
            "maximum_union4_water_count": max(any4_series),
            "process_status": status,
            "derived_counts": _hash(directory / "hydration_counts.csv"),
        }
        _check_inherited(inherited[str(trajectory)])
        LOGGER.info(
            "Completed %s in %.1fs: %s",
            replica,
            status["elapsed_seconds"],
            reports[replica],
        )
    if len(selected_checks) != 9:
        raise ValueError("Did not check all nine accepted seeds")
    for row in original_records:
        _check_inherited(row)
    return {
        "scientific_acceptance": {
            "all_6000_frames_counted": True,
            "selected_nine_seed_counts_match": True,
            "original_source_metadata_unchanged": True,
        },
        "replicas": reports,
        "selected_seed_checks": selected_checks,
        "definitions": {
            "water": "WAT oxygen atoms only; default periodic minimum-image distances retained",
            "lower_counts_only": True,
            "both35": "N(lower3.5)+OG1(lower3.5)-union(lower3.5); same water oxygen within 3.5A of both",
            "time": "full_frame*0.05 ns; historical stride10 CSV starts 0.05ns earlier",
        },
        "limitations": [
            "Descriptive saved-frame fractions; no independence or confidence intervals asserted.",
            "Distance contacts do not establish hydrogen bonding, relay orientation, or mechanism.",
            "No new molecular dynamics was run.",
        ],
    }


def main() -> None:
    """Preserve explicit acceptance/failure separately from launcher exit."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    try:
        result = run()
    except Exception as exc:
        (ROOT / "scientific_acceptance.json").write_text(
            json.dumps(
                {
                    "accepted": False,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                },
                indent=2,
            )
            + "\n"
        )
        raise
    (ROOT / "hydration_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    (ROOT / "scientific_acceptance.json").write_text(
        json.dumps({"accepted": True}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
