"""Explicit independent-state loading for supplemental relay umbrella analyses."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .relay_fixed_burnin import fixed_window_data
from .relay_pooled_parity import validate_parity
from .relay_umbrella_analysis_io import digest, load_block
from .relay_umbrella_cv import RelayBiasWindow

BASE = Path(
    "/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K"
)
FIXED = BASE / "analysis/fixed_burnin_20261126"
OLD = BASE / "additional_sampling_20261126/prepared.json"
NEW = BASE / "activation_support_20261127/prepared.json"


def checkpoint_receipts(root: Path, label: str, horizon: int) -> list[Path] | None:
    """Require the exact accepted production prefix without reading derived trees.

    Args:
        root: One cohort's physical source root.
        label: Declared native window name, preserving seed lineage.
        horizon: Inclusive production duration in5ps increments.

    Returns:
        Ordered immutable receipts orNone until all required blocks exist.
    """
    if horizon < 10 or horizon % 5:
        raise ValueError("Require10psorlatercomplete5pshorizon")
    receipts = [
        root / label / f"{i:02d}_production_5ps/analysis_acceptance.json"
        for i in range(3, horizon // 5 + 3)
    ]
    return receipts if all(p.is_file() for p in receipts) else None


def load_cohort(
    prepared_path: Path, horizon: int, reference: dict, group: str, expected: int
) -> dict | None:
    """Validate and retain each trajectory as its own MBAR/resampling state.

    Args:
        prepared_path: Immutable cohort preparation manifest.
        horizon: This cohort's exact common production horizon.
        reference: Main physical Hamiltonian/temperature/diagnostic contract.
        group: Explicit lineage namespace, not a bias-center identity.
        expected: Exact required number of native windows.

    Returns:
        Selected fixed-discard arrays and provenance, orNone if not complete.
    """
    prepared = json.loads(prepared_path.read_text())
    windows = prepared["windows"]
    if len(windows) != expected or len({w["label"] for w in windows}) != expected:
        raise ValueError("Cohortwindowinventorychanged")
    if (
        prepared["physical_hamiltonian_sha256"]
        != reference["physical_hamiltonian_sha256"]
        or prepared["temperature_K"] != reference["temperature_K"]
    ):
        raise ValueError("CannotpooldifferentphysicalHamiltonianortemperature")
    all_paths = []
    for window in windows:
        paths = checkpoint_receipts(prepared_path.parent, window["label"], horizon)
        if paths is None:
            return None
        all_paths.append(paths)
    states = []
    selected = []
    proofs = []
    samples = []
    for window, paths in zip(windows, all_paths, strict=True):
        bias = RelayBiasWindow(**window["bias"])
        blocks = []
        for path in paths:
            block = load_block(
                path,
                temperature=reference["temperature_K"],
                hamiltonian=reference["physical_hamiltonian_sha256"],
                diagnostic_columns=reference["diagnostic_columns"],
            )
            if block["window"] != window["label"]:
                raise ValueError("Receiptlineage/windowmismatch")
            parity = validate_parity(block, bias)
            blocks.append(block)
            proofs.append(
                dict(
                    path=str(path.resolve()),
                    receipt_sha256=block["receipt_sha256"],
                    source_sha256=block["source_sha256"],
                    parity=parity,
                    state_id=f"{group}/{window['label']}",
                )
            )
        values, counts, selection = fixed_window_data([window], blocks, 5, horizon)
        if counts.tolist() != [(horizon - 5) * 20]:
            raise ValueError("Independentproductiondiscard/countmismatch")
        native_rows = [
            (p, row)
            for p, b in zip(paths, blocks, strict=True)
            for row in b["original_rows"]
        ][100 : horizon * 20]
        samples.extend(
            [
                dict(
                    receipt=str(p.resolve()),
                    source_frame1based=int(row["frame"]),
                    native_time_ps=float(row["time_ps"]),
                )
                for p, row in native_rows
            ]
        )
        selected.append(values[:, :3])
        states.append(
            dict(
                state_id=f"{group}/{window['label']}",
                window=window["label"],
                group=group,
                bias=window["bias"],
                retained_frames=len(values),
                production_horizon_ps=horizon,
                discard_once_ps=5,
                source_lineage=window.get("source_lineage", {}),
                starting_restart_sha256=window.get(
                    "source_restart_sha256", window.get("seed_source_sha256")
                ),
                source_stage_seeds=[
                    stage["seed"] for stage in window["stages"] if "seed" in stage
                ],
                window_statistics=selection[0],
            )
        )
    return dict(
        q=np.concatenate(selected),
        counts=np.array([s["retained_frames"] for s in states]),
        states=states,
        proofs=proofs,
        samples=samples,
        prepared_path=str(prepared_path.resolve()),
        prepared_sha256=digest(prepared_path),
    )


def frozen_main(cache: Path) -> dict:
    """Cache the exact frozen45ps main baseline after full native parity checks.

    Args:
        cache: New pooled-analysis cache directory; no native artifacts are written.

    Returns:
        The57 independent frozen45ps states and selected physical observations.
    """
    cache.mkdir(parents=True, exist_ok=True)
    arrays = cache / "main57_frozen45ps.npz"
    metadata = cache / "main57_frozen45ps.json"
    report_path = FIXED / "fixed_45ps_analysis.json"
    report = json.loads(report_path.read_text())
    request = Path(report["checkpoint"]) / "analysis_request.json"
    config = json.loads(request.read_text())
    identity = {
        str(p): digest(p) for p in [report_path, request, BASE / "prepared.json"]
    }
    if arrays.exists():
        info = json.loads(metadata.read_text())
        if info["frozen_source_sha256"] != identity or info["arrays_sha256"] != digest(
            arrays
        ):
            raise ValueError("Frozenmainbaselineidentitychanged")
        data = np.load(arrays)
        return dict(info, q=data["q"], counts=data["counts"], config=config)
    result = load_cohort(BASE / "prepared.json", 45, config, "main57", 57)
    if result is None:
        raise ValueError("Requiredfrozen45psbaselineincomplete")
    if [
        dict(receipt_sha256=p["receipt_sha256"], source_sha256=p["source_sha256"])
        for p in result["proofs"]
    ] != report["accepted_sources"]:
        # The fixed report orders receipts lexicographically; scientific identities
        # must match regardless of explicit state-group traversal order.
        expected = {p["receipt_sha256"] for p in report["accepted_sources"]}
        if {p["receipt_sha256"] for p in result["proofs"]} != expected:
            raise ValueError("Frozen45psaccepted-sourceidentitymismatch")
    cached45 = np.load(
        FIXED / "physical_coordinate_diagnostic_20261127/matched_45ps_weights.npz"
    )
    if not np.array_equal(result["q"], cached45["q3cv_A"]) or not np.array_equal(
        result["counts"], cached45["counts"]
    ):
        raise ValueError("Frozenmainobservationsdifferfromvalidated45pscache")
    np.savez_compressed(arrays, q=result.pop("q"), counts=result.pop("counts"))
    result.update(frozen_source_sha256=identity, arrays_sha256=digest(arrays))
    metadata.write_text(json.dumps(result, indent=2) + "\n")
    data = np.load(arrays)
    return dict(result, q=data["q"], counts=data["counts"], config=config)


def pooled_data(cache: Path, horizon: int, include_new: bool) -> dict | None:
    """Pool61or94 independent states, preserving duplicate biases and seed starts.

    Args:
        cache: Derived frozen-baseline cache directory.
        horizon: Common supplemental production horizon,10through30ps.
        include_new: Include all33 new support states when true.

    Returns:
        Source-bound arrays/configuration orNone until exact common coverage.
    """
    if horizon not in (10, 15, 20, 25, 30):
        raise ValueError("Unauthorizedsupplementalhorizon")
    main = frozen_main(cache)
    groups = [main]
    for path, group, count in [(OLD, "old_targets4", 4)] + (
        [(NEW, "activation_support33", 33)] if include_new else []
    ):
        result = load_cohort(path, horizon, main["config"], group, count)
        if result is None:
            return None
        groups.append(result)
    states = [s for g in groups for s in g["states"]]
    if len({s["state_id"] for s in states}) != len(states):
        raise ValueError("Duplicateindependentstateidentity")
    proofs = [p for g in groups for p in g["proofs"]]
    csv_identities = [
        next(v for k, v in p["source_sha256"].items() if k.endswith(".csv"))
        for p in proofs
    ]
    if len(set(csv_identities)) != len(csv_identities):
        raise ValueError("Samephysicalblockduplicatedacrosspooledstates")
    return dict(
        q=np.concatenate([g["q"] for g in groups]),
        counts=np.concatenate([g["counts"] for g in groups]),
        states=states,
        proofs=proofs,
        samples=[s for g in groups for s in g["samples"]],
        config=main["config"],
        main_frozen_horizon_ps=45,
        supplemental_horizon_ps=horizon,
        bootstrap_eligible=horizon >= 15,
        prepared_sources=[
            dict(path=g["prepared_path"], sha256=g["prepared_sha256"]) for g in groups
        ],
    )
