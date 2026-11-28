"""Join completed relay screens to every saved-frame native geometry audit."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any

from nyl12_covalent.modular_es_snapshot_selection import classify_geometry

METRICS = ("d1", "d2", "d3", "attack", "projected_phi")
FLAGS = (
    "direct_strict_control",
    "strict_assisted",
    "strict_assisted_tier_ab",
    "strict_assisted_tier_ab_native",
    "nearpass_tier_ab",
)
STATES = {"nh2_asp272_asp274", "nh2_asp272_ash274", "nh3_asp272_asp274"}


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_cases(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Require precisely four matched replicas for each explicit WT state.

    Args:
        config: Relay CLI configuration.

    Returns:
        Cases keyed by unique case name.
    """
    cases = {case["name"]: case for case in config["cases"]}
    expected = {f"{state}__replica_{i:02d}" for state in STATES for i in range(1, 5)}
    if len(config["cases"]) != 12 or set(cases) != expected:
        raise ValueError(
            "Require exactly twelve matched WT cases; legacy cases forbidden"
        )
    for name, case in cases.items():
        if case["state"] != name.split("__")[0]:
            raise ValueError(f"Case/state mapping mismatch: {name}")
    return cases


def join_row(row: dict[str, Any], occupancy: dict[str, str]) -> tuple[dict, str]:
    """Validate frame identity and geometry and classify one route.

    Args:
        row: One relay route.
        occupancy: The corresponding all-frame geometry/native row.

    Returns:
        Frame flags and explicit assisted candidate class.
    """
    if row["case"] != occupancy["case"] or int(row["frame"]) != int(occupancy["frame"]):
        raise ValueError("Case/frame mapping mismatch")
    for key in ("time_ps", *METRICS):
        left, right = float(row[key]), float(occupancy[key])
        tolerance = 0.001 if key == "time_ps" else 1e-6
        if (
            not math.isfinite(left)
            or not math.isfinite(right)
            or abs(left - right) > tolerance
        ):
            raise ValueError(f"Joined {key} mismatch")
    if occupancy["native_pass"] not in {"True", "False"}:
        raise ValueError("native_pass must be True or False")
    tier = classify_geometry({key: float(row[key]) for key in METRICS})
    if occupancy["tier"] != (tier or "none"):
        raise ValueError("Tier mapping mismatch")
    assisted = row["water_count"] + row["internal_count"] > 0
    strict = row["strict"]
    if not isinstance(strict, bool):
        raise ValueError("strict must be boolean")
    native = occupancy["native_pass"] == "True"
    flags = dict(
        zip(
            FLAGS,
            (
                strict and not assisted,
                strict and assisted,
                strict and assisted and tier is not None,
                strict and assisted and tier is not None and native,
                not strict and assisted and tier is not None,
            ),
            strict=True,
        )
    )
    candidate = (
        "strict_joint_native"
        if flags[FLAGS[3]]
        else "strict_joint_native_fail"
        if flags[FLAGS[2]]
        else "nearpass_tier_ab"
        if flags[FLAGS[4]]
        else ""
    )
    return flags, candidate


def summarize(
    config_path: Path, occupancy_dir: Path, screen_dirs: list[Path], output: Path
) -> None:
    """Write unique-frame replica occupancies and across-replica statistics.

    Args:
        config_path: Exact twelve-case relay config.
        occupancy_dir: Complete per-case saved-frame CSV directory.
        screen_dirs: Completed relay output directories, possibly overlapping.
        output: New summary directory, which must not exist.
    """
    config = json.loads(config_path.read_text())
    cases = validate_cases(config)
    hashes = {str(config_path): _digest(config_path)}
    sources: dict[str, list[Path]] = {name: [] for name in cases}
    input_hashes = {}
    for directory in screen_dirs:
        summary_path = directory / "summary.json"
        summary = json.loads(summary_path.read_text())
        hashes[str(summary_path)] = _digest(summary_path)
        declared = {entry["case"]: entry for entry in summary["cases"]}
        for path in sorted(directory.glob("*_paths.jsonl")):
            name = path.name.removesuffix("_paths.jsonl")
            if name not in cases or name not in declared:
                raise ValueError(f"Unexpected or incomplete relay case: {name}")
            receipt = declared[name]
            if receipt["total_frames"] != 10000:
                raise ValueError("Screen total_frames must be10000")
            expected_inputs = {cases[name]["topology"], *cases[name]["trajectories"]}
            if set(receipt["input_sha256"]) != expected_inputs:
                raise ValueError("Screen source mapping mismatch")
            if (
                receipt["input_sha256"][cases[name]["topology"]]
                != cases[name]["topology_sha256"]
            ):
                raise ValueError("Screen topology hash mismatch")
            for source, digest in receipt["input_sha256"].items():
                if source in input_hashes and input_hashes[source] != digest:
                    raise ValueError("Conflicting screen input hashes")
                input_hashes[source] = digest
            sources[name].append(path)
            hashes[str(path)] = _digest(path)
    results, candidates = [], []
    for name, case in cases.items():
        if not sources[name]:
            raise ValueError(f"Missing completed screen: {name}")
        receipt_path = occupancy_dir / f"{name}_summary.json"
        occupancy_receipt = json.loads(receipt_path.read_text())
        hashes[str(receipt_path)] = _digest(receipt_path)
        if (
            occupancy_receipt["name"] != name
            or occupancy_receipt["state"] != case["state"]
            or occupancy_receipt["frames"] != 10000
        ):
            raise ValueError("Occupancy summary identity mismatch")
        expected_inputs = {case["topology"], *case["trajectories"]}
        fresh_hashes = occupancy_receipt["input_sha256"]
        if not expected_inputs.issubset(fresh_hashes):
            raise ValueError("Occupancy summary missing MD input hashes")
        if any(
            fresh_hashes[source] != input_hashes[source] for source in expected_inputs
        ):
            raise ValueError("Relay versus occupancy MD input hash mismatch")
        path = occupancy_dir / f"{name}.csv"
        hashes[str(path)] = _digest(path)
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != 10000 or [int(r["frame"]) for r in rows] != list(
            range(1, 10001)
        ):
            raise ValueError("Occupancy must contain all10000 sequential frames")
        if any(r["case"] != name for r in rows):
            raise ValueError("Occupancy case mapping mismatch")
        frames: dict[int, dict] = {}
        best: dict[int, dict] = {}
        for source in sources[name]:
            with source.open() as handle:
                for line in handle:
                    row = json.loads(line)
                    frame = int(row["frame"])
                    if frame < 1 or frame > 10000:
                        raise ValueError("Relay frame outside occupancy")
                    flags, category = join_row(row, rows[frame - 1])
                    flags_for_frame = frames.setdefault(
                        frame, dict.fromkeys(FLAGS, False)
                    )
                    for key in FLAGS:
                        flags_for_frame[key] |= flags[key]
                    if category:
                        rank = (
                            {
                                "strict_joint_native": 0,
                                "strict_joint_native_fail": 1,
                                "nearpass_tier_ab": 2,
                            }[category],
                            row["violation"],
                        )
                        record = {
                            **row,
                            "candidate_class": category,
                            "native_pass": rows[frame - 1]["native_pass"] == "True",
                            "source_jsonl": str(source),
                            "_rank": rank,
                        }
                        if frame not in best or rank < best[frame]["_rank"]:
                            best[frame] = record
        expected = set(range(1, 10001, 5))
        if not expected.issubset(frames):
            raise ValueError(f"Missing coarse sampled frames: {name}")
        result = {"case": name, "state": case["state"], "sampled_frames": len(frames)}
        for key in FLAGS:
            result[key] = sum(flags[key] for flags in frames.values())
            result[key + "_fraction"] = result[key] / len(frames)
        results.append(result)
        candidates.extend(best.values())
    state_rows = []
    for state in sorted(STATES):
        replicas = [r for r in results if r["state"] == state]
        result = {"state": state, "replicas": 4}
        for key in FLAGS:
            fractions = [r[key + "_fraction"] for r in replicas]
            result[key + "_fraction_mean"] = statistics.mean(fractions)
            result[key + "_fraction_sd"] = statistics.stdev(fractions)
        state_rows.append(result)
    chosen = sorted(candidates, key=lambda r: (r["_rank"], r["case"], r["frame"]))[:10]
    for row in chosen:
        row.pop("_rank")
    payload = {
        "replicas": results,
        "states": state_rows,
        "candidates": chosen,
        "source_sha256": hashes,
        "screen_input_sha256": input_hashes,
        "scope": "Geometry and conditional relay hypotheses only; no mechanism claim. "
        "Fractions use unique sampled frames; SD across four replicas, no frame IID CI. "
        "Recorded screen input hashes retained; trajectories not rehashed by joiner.",
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "summary.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n"
    )
    for filename, records in [
        ("per_replica.csv", results),
        ("per_state.csv", state_rows),
    ]:
        with (output / filename).open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
