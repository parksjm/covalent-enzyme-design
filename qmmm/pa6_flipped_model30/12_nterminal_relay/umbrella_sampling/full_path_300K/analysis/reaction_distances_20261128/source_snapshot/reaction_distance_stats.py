"""Descriptive reaction distances for independent, biased umbrella windows."""

from __future__ import annotations

import ast
import csv
import hashlib
from pathlib import Path
from typing import Any

import numpy as np

DISTANCES = ("r_og1_hg1", "r_nleave_hg1", "r_nalpha_hg1", "r_c_og1", "r_c_n")
LABELS = (
    "Thr233 OG1–HG1",
    "HG1–substrate amide N",
    "Thr233 Nα–HG1",
    "Thr233 OG1–carbonyl C",
    "Carbonyl C–amide N",
)
ROLES = {
    "carbon": 19916,
    "leaving_n": 19929,
    "og1": 13916,
    "hg1": 13917,
    "nalpha": 13905,
}


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a source file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_row(row: dict[str, str]) -> dict[str, Any]:
    """Extract distances and tracked hydrogen ownership without bond inference."""
    result: dict[str, Any] = {k: float(row[k]) for k in DISTANCES if k in row}
    # The complete contact table includes this distance regardless of ownership.
    contacts = ast.literal_eval(row["hydrogen_owner_distances_A"])
    result["r_nalpha_hg1"] = float(contacts[13917][13905])
    if "q_act_A" in row:
        residual = result["r_og1_hg1"] - result["r_nalpha_hg1"] - float(row["q_act_A"])
        if abs(residual) > 1e-6:
            raise ValueError("Direct distance / activation-coordinate mismatch")
    for coordinate, left, right in (
        ("q_del_A", "r_nalpha_hg1", "r_nleave_hg1"),
        ("q_acyl_A", "r_c_n", "r_c_og1"),
    ):
        if (
            coordinate in row
            and abs(result[left] - result[right] - float(row[coordinate])) > 1e-6
        ):
            raise ValueError(f"Direct distance / {coordinate} mismatch")
    if not all(np.isfinite(result[k]) and result[k] > 0 for k in DISTANCES):
        raise ValueError("Invalid distance")
    result["time_ps"] = float(row["time_ps"])
    for key in ("hg1_owner", "nalpha_h1_owner", "nalpha_h2_owner"):
        result[key] = int(row[key])
    return result


def load_sources(base: Path, manifest: dict[str, Any]) -> tuple[list[dict], dict]:
    """Read exactly the native blocks declared by the movie frame plan."""
    rows, provenance = [], {}
    horizon = float(manifest["sampled_horizon_per_window_ps"])
    for window in manifest["windows"]:
        blocks = list(
            dict.fromkeys(
                f["block"] for f in manifest["frame_plan"] if f["window"] == window
            )
        )
        if len(blocks) != manifest["production_blocks_per_window"]:
            raise ValueError(f"Manifest block count mismatch: {window}")
        window_rows = []
        for block in blocks:
            path = base / window / block / "raw_features.csv"
            provenance[str(path.resolve())] = sha256(path)
            with path.open(newline="") as handle:
                records = list(csv.DictReader(handle))
            if not records:
                raise ValueError(f"Empty native feature table: {path}")
            for frame, record in enumerate(records, 1):
                row = extract_row(record)
                row.update(window=window, block=block, frame=frame)
                window_rows.append(row)
        times = np.array([r["time_ps"] for r in window_rows])
        steps = np.diff(times)
        if (
            len(steps) == 0
            or np.any(steps <= 0)
            or not np.allclose(steps, np.median(steps), atol=1e-4, rtol=0)
        ):
            raise ValueError(f"Irregular native timestamps: {window}")
        origin = times[-1] - horizon
        if not np.isclose(times[0] - origin, np.median(steps), atol=1e-4):
            raise ValueError(f"Native times do not span declared horizon: {window}")
        for row in window_rows:
            row["production_time_ps"] = round(row["time_ps"] - origin, 4)
        rows.extend(window_rows)
    return rows, provenance


def retained_periods(rows: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    """Select one initial five-ps discard and first/last retained five ps."""
    retained = [r for r in rows if r["production_time_ps"] > 5]
    if not retained or retained[-1]["production_time_ps"] < 15:
        raise ValueError("Need at least ten retained ps after the five-ps discard")
    early = [r for r in retained if r["production_time_ps"] <= 10]
    late = [
        r
        for r in retained
        if r["production_time_ps"] > retained[-1]["production_time_ps"] - 5
    ]
    return retained, early, late


def campaign_metadata(
    movie: list[dict], rows: list[dict], windows: list[str], charts: dict[str, str]
) -> dict[str, Any]:
    """Derive sampling labels and phase boundaries from exact retained rows."""
    groups = []
    for window in windows:
        native = [r for r in rows if r["window"] == window]
        frames = [r["movie_frame0based"] for r in movie if r["window"] == window]
        if not frames or frames != list(range(frames[0], frames[-1] + 1)):
            raise ValueError(f"Movie window is absent or noncontiguous: {window}")
        retained, _, _ = retained_periods(native)
        groups.append(
            {
                "window": window,
                "chart": charts[window],
                "movie_start": frames[0],
                "movie_stop": frames[-1] + 1,
                "movie_count": len(frames),
                "retained_count": len(retained),
                "horizon_ps": native[-1]["production_time_ps"],
            }
        )
    boundaries = [
        i for i in range(1, len(groups)) if groups[i]["chart"] != groups[i - 1]["chart"]
    ]
    return {
        "windows": groups,
        "phase_boundary_window_positions": [i + 0.5 for i in boundaries],
        "phase_boundary_movie_frames": [
            groups[i]["movie_start"] - 0.5 for i in boundaries
        ],
        "native_count": len(rows),
        "retained_count": sum(g["retained_count"] for g in groups),
        "movie_count": len(movie),
    }


def distribution_diagnostics(
    values: np.ndarray, early_values: np.ndarray, late_values: np.ndarray
) -> dict[str, float | None]:
    """Measure bin-dependent early/late overlap and adjacent-frame correlation."""
    low = np.floor(values.min() / 0.05) * 0.05
    high = np.ceil(values.max() / 0.05) * 0.05 + 0.05
    edges = np.arange(low, high + 0.025, 0.05)
    early = np.histogram(early_values, bins=edges)[0] / len(early_values)
    late = np.histogram(late_values, bins=edges)[0] / len(late_values)
    # Correlation is undefined for constant trajectories; report JSON null.
    lag = (
        float(np.corrcoef(values[:-1], values[1:])[0, 1])
        if np.std(values) > 0
        else None
    )
    return {
        "early_late_histogram_intersection_0p05A": float(np.minimum(early, late).sum()),
        "lag1_correlation": lag,
    }


def summarize(rows: list[dict], windows: list[str]) -> list[dict]:
    """Summarize post-discard distributions and first/last retained five ps."""
    output = []
    for window in windows:
        selected, early, late = retained_periods(
            [r for r in rows if r["window"] == window]
        )
        for key in DISTANCES:
            values = np.array([r[key] for r in selected])
            result = {
                "window": window,
                "distance": key,
                "n": len(values),
                "mean_A": float(values.mean()),
                "std_A": float(values.std()),
                "min_A": float(values.min()),
                "max_A": float(values.max()),
            }
            for q in (5, 25, 50, 75, 95):
                result[f"p{q}_A"] = float(np.percentile(values, q))
            early_values = np.array([r[key] for r in early])
            late_values = np.array([r[key] for r in late])
            result["early_5ps_mean_A"] = float(early_values.mean())
            result["late_5ps_mean_A"] = float(late_values.mean())
            result["late_minus_early_A"] = float(
                late_values.mean() - early_values.mean()
            )
            result.update(distribution_diagnostics(values, early_values, late_values))
            output.append(result)
    return output


def ownership_summary(rows: list[dict], windows: list[str]) -> list[dict]:
    """Report early/late owner fractions for HG1 and both original N hydrogens."""
    output = []
    for window in windows:
        retained, early, late = retained_periods(
            [r for r in rows if r["window"] == window]
        )
        for period, selected in (
            ("early_5ps", early),
            ("late_5ps", late),
            ("retained_after_5ps", retained),
        ):
            for key in ("hg1_owner", "nalpha_h1_owner", "nalpha_h2_owner"):
                owners = [r[key] for r in selected]
                output.append(
                    {
                        "window": window,
                        "period": period,
                        "hydrogen": key,
                        "owner_fractions": {
                            str(k): owners.count(k) / len(owners)
                            for k in sorted(set(owners))
                        },
                    }
                )
    return output


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write a nonempty rectangular descriptive table."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def movie_rows(rows: list[dict], plan: list[dict]) -> list[dict]:
    """Join exact movie frame indices to native frame observables."""
    lookup = {(r["window"], r["block"], r["frame"]): r for r in rows}
    result = []
    for expected, frame in enumerate(plan):
        if frame["output_frame0based"] != expected:
            raise ValueError("Nonsequential movie frame plan")
        row = lookup[(frame["window"], frame["block"], frame["source_frame1based"])]
        if abs(row["time_ps"] - frame["window_time_ps"]) > 1e-4:
            raise ValueError("Movie/native timestamp mismatch")
        result.append({"movie_frame0based": expected, **row})
    return result
