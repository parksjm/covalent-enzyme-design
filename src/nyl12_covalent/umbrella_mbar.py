"""PyMBAR diagnostics for one-dimensional umbrella-sampling campaigns."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from pymbar import FES, MBAR, timeseries

GAS_CONSTANT_KCAL_MOL_K = 0.00198720425864083


def load_umbrella_samples(
    campaign: Path,
    stage: str,
    *,
    require_successful_process: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Load one explicitly selected stage from every umbrella window.

    Args:
        campaign: Campaign directory containing ``inputs.json`` and windows.
        stage: Exact stage-directory name; glob patterns are rejected.
        require_successful_process: Require an exit-status-zero ``process.json``.
            Set false only for an immutable snapshot of an active campaign.

    Returns:
        Campaign metadata and ordered window sample records.

    Raises:
        ValueError: If the stage is ambiguous, incomplete, or inconsistent.
    """
    if any(character in stage for character in "*?["):
        raise ValueError("MBAR smoke tests require one exact stage name")
    inputs_path = campaign / "inputs.json"
    inputs = json.loads(inputs_path.read_text())
    rows: list[dict[str, Any]] = []
    for window in inputs["windows"]:
        stage_directory = campaign / window["label"] / stage
        coordinate_path = stage_directory / "umbrella_coordinate.dat"
        if not coordinate_path.is_file():
            raise ValueError(f"Missing umbrella coordinates: {coordinate_path}")
        if require_successful_process:
            process_path = stage_directory / "process.json"
            if not process_path.is_file():
                raise ValueError(f"Missing process record: {process_path}")
            process = json.loads(process_path.read_text())
            if process.get("exit_status") != 0:
                raise ValueError(f"Umbrella stage is not complete: {stage_directory}")
        samples = _read_dumpave(coordinate_path)
        rows.append(
            {
                "window": str(window["label"]),
                "target_A": float(window["target_A"]),
                "samples_A": samples,
                "source": str(coordinate_path.resolve()),
                "sha256": hashlib.sha256(coordinate_path.read_bytes()).hexdigest(),
            }
        )
    if len(rows) < 2:
        raise ValueError("At least two complete windows are required")
    return inputs, rows


def load_umbrella_sample_blocks(
    sources: list[tuple[Path, str]],
    *,
    require_successful_process: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Load and concatenate matching production blocks across campaigns.

    Args:
        sources: Ordered ``(campaign, exact_stage_name)`` pairs.
        require_successful_process: Require every selected stage to have an
            exit-status-zero process record.

    Returns:
        Common campaign metadata and ordered multi-block window records.

    Raises:
        ValueError: If no sources are supplied or window metadata differ.
    """
    if not sources:
        raise ValueError("At least one umbrella source block is required")
    loaded = [
        load_umbrella_samples(
            campaign,
            stage,
            require_successful_process=require_successful_process,
        )
        for campaign, stage in sources
    ]
    reference_inputs, reference_windows = loaded[0]
    temperature = float(reference_inputs["temperature_K"])
    force_constant = float(reference_inputs["force_constant_kcal_mol_A2"])
    for inputs, windows in loaded[1:]:
        if not np.isclose(float(inputs["temperature_K"]), temperature):
            raise ValueError("Umbrella source temperatures differ")
        if not np.isclose(float(inputs["force_constant_kcal_mol_A2"]), force_constant):
            raise ValueError("Umbrella source force constants differ")
        reference_grid = [(row["window"], row["target_A"]) for row in reference_windows]
        current_grid = [(row["window"], row["target_A"]) for row in windows]
        if current_grid != reference_grid:
            raise ValueError("Umbrella source window grids differ")

    combined: list[dict[str, Any]] = []
    for window_index, reference in enumerate(reference_windows):
        blocks = [windows[window_index] for _, windows in loaded]
        samples = [np.asarray(block["samples_A"], dtype=float) for block in blocks]
        digests = [block["sha256"] for block in blocks]
        combined.append(
            {
                "window": reference["window"],
                "target_A": reference["target_A"],
                "samples_A": np.concatenate(samples),
                "sample_blocks_A": samples,
                "source": ";".join(block["source"] for block in blocks),
                "sha256": hashlib.sha256("".join(digests).encode()).hexdigest(),
                "source_blocks": [
                    {
                        "source": block["source"],
                        "sha256": block["sha256"],
                        "sample_count": int(block["samples_A"].size),
                    }
                    for block in blocks
                ],
            }
        )
    return reference_inputs, combined


def run_mbar_smoke(
    windows: list[dict[str, Any]],
    *,
    temperature_k: float,
    force_constant_kcal_mol_a2: float,
    bins: int,
) -> dict[str, Any]:
    """Run decorrelation, MBAR overlap, and histogram-FES diagnostics.

    Amber's flat-bottom restraint with ``r2 == r3`` uses
    ``k * (coordinate - target)**2`` outside its zero-width center.

    Args:
        windows: Ordered records returned by :func:`load_umbrella_samples`.
        temperature_k: Simulation temperature in kelvin.
        force_constant_kcal_mol_a2: Amber ``rk2 == rk3`` value.
        bins: Number of common histogram bins for the provisional FES.

    Returns:
        JSON-serializable smoke-test diagnostics and provisional PMF points.

    Raises:
        ValueError: If inputs cannot support the requested diagnostic.
    """
    if temperature_k <= 0.0 or force_constant_kcal_mol_a2 <= 0.0:
        raise ValueError("Temperature and force constant must be positive")
    if bins < 4:
        raise ValueError("At least four FES bins are required")
    targets = np.asarray([row["target_A"] for row in windows], dtype=float)
    if not np.all(np.diff(targets) > 0.0):
        raise ValueError("Umbrella targets must be strictly increasing")

    decorrelated: list[np.ndarray] = []
    window_summary: list[dict[str, Any]] = []
    for row in windows:
        samples = np.asarray(row["samples_A"], dtype=float)
        if samples.size < 4 or not np.isfinite(samples).all():
            raise ValueError(f"Insufficient finite samples in {row['window']}")
        raw_blocks = row.get("sample_blocks_A", [samples])
        selected_blocks: list[np.ndarray] = []
        block_summary: list[dict[str, Any]] = []
        for block_index, raw_block in enumerate(raw_blocks, start=1):
            block = np.asarray(raw_block, dtype=float)
            if block.size < 4 or not np.isfinite(block).all():
                raise ValueError(
                    f"Insufficient finite samples in {row['window']} block {block_index}"
                )
            block_inefficiency = float(timeseries.statistical_inefficiency(block))
            indices = np.asarray(
                timeseries.subsample_correlated_data(
                    block,
                    g=block_inefficiency,
                    conservative=True,
                ),
                dtype=int,
            )
            selected_block = block[indices]
            if selected_block.size < 2:
                raise ValueError(
                    f"Fewer than two decorrelated samples in {row['window']} "
                    f"block {block_index}"
                )
            selected_blocks.append(selected_block)
            block_summary.append(
                {
                    "block": block_index,
                    "raw_sample_count": int(block.size),
                    "statistical_inefficiency": block_inefficiency,
                    "decorrelated_sample_count": int(selected_block.size),
                }
            )
        selected = np.concatenate(selected_blocks)
        if selected.size < 2:
            raise ValueError(f"Fewer than two decorrelated samples in {row['window']}")
        decorrelated.append(selected)
        window_summary.append(
            {
                "window": row["window"],
                "target_A": row["target_A"],
                "raw_sample_count": int(samples.size),
                "statistical_inefficiency": float(samples.size / selected.size),
                "decorrelated_sample_count": int(selected.size),
                "mean_A": float(np.mean(samples)),
                "standard_deviation_A": float(np.std(samples, ddof=1)),
                "source": row["source"],
                "sha256": row["sha256"],
                "blocks": block_summary,
            }
        )

    counts = np.asarray([values.size for values in decorrelated], dtype=int)
    coordinate = np.concatenate(decorrelated)
    beta = 1.0 / (GAS_CONSTANT_KCAL_MOL_K * temperature_k)
    reduced_bias = (
        beta
        * force_constant_kcal_mol_a2
        * np.square(coordinate[None, :] - targets[:, None])
    )
    mbar = MBAR(reduced_bias, counts, verbose=False)
    overlap = mbar.compute_overlap()
    mbar_effective = np.asarray(mbar.compute_effective_sample_number(), dtype=float)

    edges = np.linspace(float(coordinate.min()), float(coordinate.max()), bins + 1)
    # ``np.digitize`` assigns a value exactly equal to the final edge to an
    # out-of-range bin.  Extend only that edge by one floating-point unit so
    # the observed maximum remains inside the requested histogram domain.
    edges[-1] = np.nextafter(edges[-1], np.inf)
    centers = 0.5 * (edges[:-1] + edges[1:])
    fes = FES(reduced_bias, counts, mbar_options={"verbose": False})
    fes.generate_fes(
        np.zeros(coordinate.size),
        coordinate,
        fes_type="histogram",
        histogram_parameters={"bin_edges": edges},
    )
    occupied_indices = sorted(
        int(label[0]) for label in fes.histogram_data["bin_label"]
    )
    occupied_centers = centers[occupied_indices]
    fes_result = fes.get_fes(
        occupied_centers,
        reference_point="from-lowest",
        uncertainty_method="analytical",
    )
    thermal_energy = GAS_CONSTANT_KCAL_MOL_K * temperature_k
    profile = [
        {
            "coordinate_A": float(position),
            "free_energy_kcal_mol": float(value * thermal_energy),
            "uncertainty_kcal_mol": float(error * thermal_energy),
        }
        for position, value, error in zip(
            occupied_centers,
            fes_result["f_i"],
            fes_result["df_i"],
            strict=True,
        )
    ]
    overlap_matrix = np.asarray(overlap["matrix"], dtype=float)
    adjacent = [
        float(min(overlap_matrix[index, index + 1], overlap_matrix[index + 1, index]))
        for index in range(len(windows) - 1)
    ]
    return {
        "interpretation": "pipeline smoke test; not a converged or promotable PMF",
        "temperature_K": temperature_k,
        "force_constant_kcal_mol_A2": force_constant_kcal_mol_a2,
        "amber_restraint_energy_convention": "k*(xi-target)^2",
        "window_count": len(windows),
        "raw_sample_count": int(sum(row["raw_sample_count"] for row in window_summary)),
        "decorrelated_sample_count": int(counts.sum()),
        "windows": window_summary,
        "mbar_effective_sample_number": mbar_effective.tolist(),
        "mbar_overlap_scalar": float(overlap["scalar"]),
        "mbar_overlap_eigenvalues": np.asarray(overlap["eigenvalues"]).tolist(),
        "mbar_overlap_matrix": overlap_matrix.tolist(),
        "minimum_bidirectional_adjacent_overlap": min(adjacent),
        "fes_bin_edges_A": edges.tolist(),
        "provisional_pmf": profile,
    }


def _read_dumpave(path: Path) -> np.ndarray:
    values: list[float] = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 2:
            raise ValueError(f"Malformed DUMPAVE line {line_number} in {path}")
        try:
            int(fields[0])
            values.append(float(fields[1]))
        except ValueError as error:
            raise ValueError(
                f"Non-numeric DUMPAVE line {line_number} in {path}"
            ) from error
    result = np.asarray(values, dtype=float)
    if result.size == 0 or not np.isfinite(result).all():
        raise ValueError(f"No finite umbrella samples in {path}")
    return result
