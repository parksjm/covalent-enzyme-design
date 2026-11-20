"""Incremental, report-only comparison of explicitly configured NEB trials."""

from __future__ import annotations

import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Any

from .neb_progress import parse_last_complete_neb_report

BLOCK = re.compile(rb"^NEB Report\r?\n.*?^\s*RMS grad[^\n]*\n", re.M | re.S)
GRADIENT = re.compile(r"(Max|RMS) grad\s+([-+\d.EeDd]+)\s+Target:\s*([-+\d.EeDd]+)")
RESET = re.compile(
    r"(?:reset|restart).*(?:bfgs|memory|history)|(?:bfgs|memory|history).*reset", re.I
)
IMAGE_RESET = re.compile(
    r"^[ \t]*Resetting images[ \t]+(\d+)[ \t]+to[ \t]+(\d+)[ \t]*\r?\n?$", re.I
)
CHECKPOINTS = (50, 100, 250)


def digest(path: Path) -> str:
    """Return SHA-256 of a small configuration/source artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, record: Any) -> None:
    """Replace a JSON inspection artifact atomically on the same filesystem."""
    pending = path.with_name(path.name + ".pending")
    pending.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    pending.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Atomically write a compact rectangular inspection table."""
    pending = path.with_name(path.name + ".pending")
    with pending.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(dict.fromkeys(key for row in rows for key in row))
        )
        writer.writeheader()
        writer.writerows(rows)
    pending.replace(path)


def parse_block(payload: bytes, count: int) -> dict[str, Any]:
    """Parse a closed report plus its complete convergence block.

    Args:
        payload: Native report through the newline-terminated RMS-gradient test.
        count: Explicit regular-image count,24or48.

    Returns:
        Physical profile, native gradient values/targets/ratios and spacing summary.

    Raises:
        ValueError: If report, cycle or gradient targets are invalid.
    """
    report = parse_last_complete_neb_report(payload, count)
    values = {
        kind: (
            float(value.replace("D", "E").replace("d", "e")),
            float(target.replace("D", "E").replace("d", "e")),
        )
        for kind, value, target in GRADIENT.findall(payload.decode("ascii"))
    }
    if (
        report.cycle is None
        or set(values) != {"Max", "RMS"}
        or any(t <= 0 for _, t in values.values())
    ):
        raise ValueError(
            "Complete native cycle and positive Max/RMS gradient targets required"
        )
    distances = [
        float(line[43:53])
        for line in payload.decode("ascii").splitlines()
        if line.startswith("Img ")
    ]
    nonzero = [v for v in distances if v > 0]
    peak = max(report.images, key=lambda image: image.energy_hartree)
    return {
        "cycle": report.cycle,
        "regular_images": count,
        "physical_max_kcal_mol": peak.relative_energy_kcal_mol,
        "peak_image": peak.image,
        "endpoint_kcal_mol": report.images[-1].relative_energy_kcal_mol,
        "max_gradient": values["Max"][0],
        "max_gradient_target": values["Max"][1],
        "max_gradient_ratio": values["Max"][0] / values["Max"][1],
        "rms_gradient": values["RMS"][0],
        "rms_gradient_target": values["RMS"][1],
        "rms_gradient_ratio": values["RMS"][0] / values["RMS"][1],
        "native_spacing_max_min_ratio": max(nonzero) / min(nonzero)
        if nonzero
        else None,
        "native_block_sha256": hashlib.sha256(payload).hexdigest(),
        "images": [asdict(image) for image in report.images],
        "climbing_image": asdict(report.climbing_image)
        if report.climbing_image
        else None,
    }


def chemical_gaps(trial: Path, cycle: int) -> dict[str, Any]:
    """Use chemical differences only from an immutable coherent exported same-cycle band."""
    latest = trial / "progress/latest"
    if not (latest / "progress.json").is_file():
        return {"available": False, "reason": "No coherent coordinate export"}
    latest = latest.resolve()
    report = json.loads((latest / "progress.json").read_text())
    if report["cycle"] != cycle:
        return {
            "available": False,
            "reason": "Coherent export cycle differs",
            "export_cycle": report["cycle"],
        }
    with (latest / "images.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    keys = ("rOG1C", "rCN", "rCO", "CN_OG1", "CN_leave", "CN_water")
    gaps = {}
    for key in keys:
        if key in rows[0]:
            values = [float(row[key]) for row in rows]
            gaps[key] = max(abs(b - a) for a, b in zip(values, values[1:]))
    return {
        "available": True,
        "coherent_export": str(latest),
        "progress_sha256": digest(latest / "progress.json"),
        "max_adjacent_absolute_changes": gaps,
    }


def update_trial(
    config_path: Path, output: Path, max_bytes: int = 8 * 1024 * 1024
) -> dict[str, Any]:
    """Read bounded new native bytes and retain complete-report history.

    ``reset_count`` preserves its legacy BFGS-memory meaning and equals
    ``bfgs_memory_reset_count``. ``image_reset_count`` counts native geometric
    "Resetting images ... to ..." events, each once, not the images affected.
    Old cursors migrate with one read of their already-consumed complete lines.

    Args:
        config_path: Explicit trial config with regular_images24or48 and4662active.
        output: This monitor's dedicated per-trial inspection directory.
        max_bytes: Upper bound on newly consumed bytes per refresh.

    Returns:
        Explicit waiting/current status; no barrier acceptance is inferred.

    Raises:
        ValueError: If scientific configuration changes, output truncates or reports fail validation.
    """
    if max_bytes <= 0:
        raise ValueError("max_bytes must be positive")
    cfg = json.loads(config_path.read_text())
    if (
        type(cfg["regular_images"]) is not int
        or cfg["regular_images"] not in (24, 48)
        or cfg["active_atoms"] != 4662
    ):
        raise ValueError(
            "Expected explicitly configured24or48regular images and4662active atoms"
        )
    case = Path(cfg["case_path"])
    output.mkdir(parents=True, exist_ok=True)
    state_path = output / "cursor.json"
    identity = {"config": digest(config_path), "case": digest(case)}
    state = (
        json.loads(state_path.read_text())
        if state_path.exists()
        else {
            "identity": identity,
            "offset": 0,
            "pending": "",
            "event_tail": "",
            "report_count": 0,
            "reset_count": 0,
            "reset_offsets": [],
            "image_reset_offsets": [],
            "first_observed_utc": datetime.now(timezone.utc).isoformat(),
        }
    )
    if state["identity"] != identity:
        raise ValueError(
            "Trial config/case changed; use a new monitor output directory"
        )
    native = config_path.parent / "calculation/minimized.mdout"
    if not native.is_file():
        return {"status": "waiting_for_native_output", "trial": config_path.parent.name}
    stat = native.stat()
    if stat.st_size < state["offset"] or state.get("inode", stat.st_ino) != stat.st_ino:
        raise ValueError(
            "Native output truncated or replaced; use a new trial directory"
        )
    with native.open("rb") as stream:
        if "image_reset_offsets" not in state:
            complete_end = state["offset"] - len(state["event_tail"].encode("ascii"))
            prefix = stream.read(complete_end).decode("ascii")
            state["image_reset_offsets"] = [
                match.end()
                for match in re.finditer(IMAGE_RESET.pattern, prefix, re.I | re.M)
            ]
        stream.seek(state["offset"])
        new = stream.read(max_bytes)
    state.update(offset=state["offset"] + len(new), inode=stat.st_ino)
    now = datetime.now(timezone.utc)
    job_path = config_path.parent / "job/status.json"
    job = json.loads(job_path.read_text()) if job_path.is_file() else {}
    started = job.get("started_utc") or state["first_observed_utc"]
    elapsed = (
        now - datetime.fromisoformat(started.replace("Z", "+00:00"))
    ).total_seconds()
    event_text = state["event_tail"] + new.decode("ascii")
    lines = event_text.splitlines(keepends=True)
    state["event_tail"] = lines.pop() if lines and not lines[-1].endswith("\n") else ""
    events = []
    line_offset = state["offset"] - len(event_text.encode("ascii"))
    for line in lines:
        line_offset += len(line.encode("ascii"))
        cycle = re.search(r"Testing convergence\s+in cycle\s+(\d+)", line)
        if cycle:
            state["last_event_cycle"] = int(cycle[1])
        if IMAGE_RESET.fullmatch(line):
            state["image_reset_offsets"].append(line_offset)
        if RESET.search(line):
            state["reset_offsets"].append(line_offset)
            events.append(
                {
                    "line": line.strip(),
                    "source_end_byte": line_offset,
                    "after_cycle": state.get("last_event_cycle"),
                }
            )
    state["reset_count"] = len(state["reset_offsets"])
    state["bfgs_memory_reset_count"] = state["reset_count"]
    state["image_reset_count"] = len(state["image_reset_offsets"])
    if events:
        with (output / "resets.jsonl").open("a") as stream:
            stream.write(
                json.dumps(
                    {
                        "observed_utc": now.isoformat(),
                        "elapsed_s": elapsed,
                        "events": events,
                    }
                )
                + "\n"
            )
    payload = state["pending"].encode("ascii") + new
    matches = list(BLOCK.finditer(payload))
    for match in matches:
        record = parse_block(match.group(), cfg["regular_images"])
        state["report_count"] += 1
        record.update(
            report_count=state["report_count"],
            reset_count=sum(
                offset <= state["offset"] - len(payload) + match.end()
                for offset in state["reset_offsets"]
            ),
            bfgs_memory_reset_count=sum(
                offset <= state["offset"] - len(payload) + match.end()
                for offset in state["reset_offsets"]
            ),
            image_reset_count=sum(
                offset <= state["offset"] - len(payload) + match.end()
                for offset in state["image_reset_offsets"]
            ),
            elapsed_s=elapsed,
            observed_utc=now.isoformat(),
            source_end_byte=state["offset"] - len(payload) + match.end(),
            source_start_byte=state["offset"] - len(payload) + match.start(),
            native_output=str(native),
            elapsed_reference="launcher started_utc"
            if job.get("started_utc")
            else "monitor first observation",
        )
        record["chemical_gaps"] = chemical_gaps(config_path.parent, record["cycle"])
        atomic_json(output / "current_profile.json", record)
        write_csv(output / "current_profile.csv", record["images"])
        summary = {k: v for k, v in record.items() if not isinstance(v, (dict, list))}
        with (output / "history.jsonl").open("a") as stream:
            stream.write(json.dumps(summary, allow_nan=False) + "\n")
        for checkpoint in CHECKPOINTS:
            path = output / f"checkpoint_{checkpoint:04d}.json"
            if record["cycle"] >= checkpoint and not path.exists():
                atomic_json(path, {"requested_cycle": checkpoint, "actual": record})
    remainder = payload[matches[-1].end() :] if matches else payload
    marker = remainder.rfind(b"NEB Report\n")
    state["pending"] = (remainder[marker:] if marker >= 0 else remainder[-32:]).decode(
        "ascii"
    )
    atomic_json(state_path, state)
    profile = output / "current_profile.json"
    return {
        "status": "complete_report_available"
        if profile.exists()
        else "waiting_for_complete_report",
        "trial": config_path.parent.name,
        "terminal": job.get("ended_utc") is not None,
        "report_count": state["report_count"],
        "reset_count": state["reset_count"],
        "bfgs_memory_reset_count": state["bfgs_memory_reset_count"],
        "image_reset_count": state["image_reset_count"],
        "bytes_consumed": state["offset"],
        "native_size_at_read": stat.st_size,
    }
