"""Immutable compact-block receipts and a runnable relay analysis request."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


def digest(path: Path) -> str:
    """Return the SHA256 identity of an existing file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _resolve(parent: Path, name: str) -> Path:
    return (parent / name).resolve(strict=True)


def load_block(
    receipt_path: Path,
    *,
    temperature: float,
    hamiltonian: str,
    diagnostic_columns: list[str],
) -> dict[str, Any]:
    """Validate immutable acceptance/data/process proofs without chemistry cuts.

    Args:
        receipt_path: Immutable accepted production-block receipt.
        temperature: Exact common physical temperature.
        hamiltonian: Common physical Hamiltonian source identity.
        diagnostic_columns: Numeric orthogonal diagnostics used for decorrelation.

    Returns:
        Analysis-ready arrays, all original rows, and immutable source identities.
    """
    proof = json.loads(receipt_path.read_text())
    if (
        proof["accepted"] is not True
        or proof["exit_status"] != 0
        or proof["frames_censored"] != 0
        or proof["coordinates_finite"] is not True
        or proof["energies_finite"] is not True
    ):
        raise ValueError(
            "Production acceptance/process/finite/censoring contract failed"
        )
    if (
        proof["temperature_K"] != temperature
        or proof["physical_hamiltonian_sha256"] != hamiltonian
    ):
        raise ValueError("Physical temperature or common Hamiltonian mismatch")
    csv_path = _resolve(receipt_path.parent, proof["csv_file"])
    process_path = _resolve(receipt_path.parent, proof["process_file"])
    expected = {
        str(csv_path): proof["csv_sha256"],
        str(process_path): proof["process_sha256"],
    }
    if not proof["source_sha256"]:
        raise ValueError("Immutable source proof mapping cannot be empty")
    for filename, sha in proof["source_sha256"].items():
        source = _resolve(receipt_path.parent, filename)
        if source.name in {"status.json", "campaign_status.json"}:
            raise ValueError("Mutable status metadata cannot establish source identity")
        expected[str(source)] = sha
    for filename, sha in expected.items():
        if digest(Path(filename)) != sha:
            raise ValueError(
                f"Immutable accepted-block source hash mismatch: {filename}"
            )
    process = json.loads(process_path.read_text())
    if process.get("exit_status") != 0:
        raise ValueError("Immutable process proof was not successful")
    with csv_path.open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("Accepted block has no observable rows")
    if "frames" in proof and len(rows) != proof["frames"]:
        raise ValueError("Accepted frame count differs from CSV")
    feature_columns = [
        name + "_A" if name + "_A" in rows[0] else name
        for name in ("q_acyl", "q_act", "q_del")
    ]
    columns = feature_columns + diagnostic_columns
    values = np.array([[float(row[name]) for name in columns] for row in rows])
    mask_text = [str(row["native_bias_evaluated"]).lower() for row in rows]
    if any(value not in {"true", "false", "1", "0"} for value in mask_text):
        raise ValueError(
            "Native parity mask must contain explicit boolean or zero/one values"
        )
    mask = np.array([value in {"true", "1"} for value in mask_text])
    native = np.array(
        [
            float(row["native_bias_kcal_mol"]) if evaluated else np.nan
            for row, evaluated in zip(rows, mask, strict=True)
        ]
    )
    parity = dict(
        native_bias_evaluated=mask,
        native_bias_source=proof["native_bias_source"],
        native_cv_quantization_A=proof["native_cv_quantization_A"],
        native_energy_quantization_kcal_mol=proof.get(
            "native_energy_quantization_kcal_mol", 0.0
        ),
        coordinate_serialization_tolerance_A=proof[
            "coordinate_serialization_tolerance_A"
        ],
    )
    if parity["native_bias_source"] != "amber_dumpave_independent_CVs":
        raise ValueError(
            "Require independently emitted Amber DUMPAVE CV parity evidence"
        )
    if parity["native_bias_source"] == "amber_dumpave_independent_CVs":
        parity["native_cv_values_A"] = np.array(
            [
                [
                    float(row[name]) if evaluated else np.nan
                    for name in ("native_cv_1_A", "native_cv_2_A")
                ]
                for row, evaluated in zip(rows, mask, strict=True)
            ]
        )
        if (
            int(mask.sum()) != proof["native_bias_evaluated_frames"]
            or int((~mask).sum()) != proof["native_bias_unevaluated_terminal_frames"]
        ):
            raise ValueError("Native parity mask/proof counts differ")
    time = np.array([float(row["time_ps"]) for row in rows])
    if not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
        raise ValueError("Per-block times must be finite and strictly increasing")
    if not np.isfinite(values).all() or not np.isfinite(native[mask]).all():
        raise ValueError("Nonfinite compact production observables or native bias")
    return dict(
        window=proof["window"],
        values=values,
        native_bias_kcal_mol=native,
        **parity,
        sampled_duration_ps=proof["sampled_duration_ps"],
        original_rows=rows,
        receipt_sha256=digest(receipt_path),
        source_sha256=expected,
        diagnostic_columns=diagnostic_columns,
    )


def analyze_request(request_path: Path) -> dict[str, Any]:
    """Execute a source-bound compact analysis request without launching MD.

    Args:
        request_path: JSON request with windows, path_3cv_A, temperature_K,
            physical_hamiltonian_sha256, s_bin_edges, es_search_bounds,
            neighbor_pairs, diagnostic_columns and accepted_block_receipts.

    Returns:
        Analysis and uncensored source/proton-owner records suitable for JSON.
    """
    from .relay_umbrella_analysis import analyze

    cfg = json.loads(request_path.read_text())
    hamiltonian = cfg["physical_hamiltonian_sha256"]
    if len(hamiltonian) != 64 or any(c not in "0123456789abcdef" for c in hamiltonian):
        raise ValueError("Hamiltonian identity must be an explicit SHA256 digest")
    receipts = [
        _resolve(request_path.parent, name) for name in cfg["accepted_block_receipts"]
    ]
    if len(set(receipts)) != len(receipts):
        raise ValueError("An accepted block receipt cannot be counted twice")
    blocks = [
        load_block(
            path,
            temperature=cfg["temperature_K"],
            hamiltonian=hamiltonian,
            diagnostic_columns=cfg["diagnostic_columns"],
        )
        for path in receipts
    ]
    source_sets = [frozenset(block["source_sha256"].items()) for block in blocks]
    if len(set(source_sets)) != len(source_sets):
        raise ValueError("The same accepted production data cannot be counted twice")
    result = analyze(
        cfg["windows"],
        blocks,
        path=np.array(cfg["path_3cv_A"]),
        temperature=cfg["temperature_K"],
        edges=np.array(cfg["s_bin_edges"]),
        es_bounds=tuple(cfg["es_search_bounds"]),
        acyl_bounds=tuple(cfg["acyl_search_bounds"]),
        neighbor_pairs=cfg["neighbor_pairs"],
    )
    if cfg.get("bootstrap"):
        from .relay_umbrella_bootstrap import bootstrap

        settings = cfg["bootstrap"]
        result["gates"]["bootstrap_95CI"] = bootstrap(
            cfg["windows"],
            blocks,
            analysis_options=dict(
                path=np.array(cfg["path_3cv_A"]),
                temperature=cfg["temperature_K"],
                edges=np.array(cfg["s_bin_edges"]),
                es_bounds=tuple(cfg["es_search_bounds"]),
                acyl_bounds=tuple(cfg["acyl_search_bounds"]),
                neighbor_pairs=cfg["neighbor_pairs"],
            ),
            replicates=settings["replicates"],
            chunk_frames=settings["chunk_frames"],
            seed=settings["seed"],
        )
    result.update(
        request_sha256=digest(request_path),
        path_3cv_sha256=hashlib.sha256(
            np.asarray(cfg["path_3cv_A"], dtype="<f8").tobytes()
        ).hexdigest(),
        physical_hamiltonian_sha256=hamiltonian,
        temperature_K=cfg["temperature_K"],
        diagnostic_columns=cfg["diagnostic_columns"],
        accepted_sources=[
            dict(
                receipt_sha256=block["receipt_sha256"],
                source_sha256=block["source_sha256"],
            )
            for block in blocks
        ],
    )
    from .relay_umbrella_analysis_reports import chemical_diagnostics

    result["chemical_diagnostics"] = chemical_diagnostics(
        blocks, cfg["diagnostic_columns"]
    )
    return result
