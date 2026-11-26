"""Strict legacy-US adapter for the shared all-frame fixed-burn-in estimator."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import Rst7

from .neb_umbrella_bias import HarmonicWindow, amber_restraints
from .relay_umbrella_analysis import KB
from .relay_umbrella_analysis_io import digest
from .relay_umbrella_cv import RelayBiasWindow, cross_bias_energies, raw_features
from .relay_umbrella_projection import project_path

DIAGNOSTICS = [
    "r_c_n",
    "r_c_og1",
    "r_og1_hg1",
    "r_nleave_hg1",
    "cn_nleave_hpool",
    "cn_og1_hpool",
    "cn_progress",
    "cn_nalpha_hpool",
    "r_nalpha_h1",
    "r_nalpha_h2",
]


def direct_reference(prepared: dict[str, Any]) -> tuple[np.ndarray, dict[str, str]]:
    """Recover the original 24-image NEB geometry in the same three-CV metric.

    Args:
        prepared: Original direct campaign manifest, including source restart hashes.

    Returns:
        Ordered q_acyl/q_act/q_del points and verified reference-file identities.
    """
    images: dict[int, tuple[str, str]] = {}
    for window in prepared["windows"]:
        identity = (window["source_restart"], window["source_restart_sha256"])
        number = int(window["seed_image"])
        if number in images and images[number] != identity:
            raise ValueError("A NEB image has conflicting restart identities")
        images[number] = identity
    if sorted(images) != list(range(1, 25)):
        raise ValueError("Require exactly the original 24 ordered NEB images")
    points, hashes = [], {}
    for filename, sha in (images[index] for index in sorted(images)):
        path = Path(filename)
        if digest(path) != sha:
            raise ValueError("Original NEB reference restart changed: " + filename)
        restart = Rst7.open(str(path))
        features = raw_features(
            np.asarray(restart.coordinates).reshape(-1, 3), restart.box
        )
        points.append([features[name] for name in ("q_acyl", "q_act", "q_del")])
        hashes[filename] = sha
    return np.asarray(points), hashes


def load_direct_dataset(
    campaign: Path,
    augmented: Path,
    relay_request: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, str]]:
    """Validate augmented direct data and recorded cross-biases without censoring.

    Args:
        campaign: Existing direct-US campaign root.
        augmented: New sparse-coordinate augmentation root, window/stage ordered.
        relay_request: Frozen relay request defining the shared bins and basin ranges.

    Returns:
        Shared-statistics config, all continuous production blocks, and source hashes.
        Native biases remain the original q_acyl/q_H restraints, not relay biases.
    """
    prepared_path = campaign / "prepared.json"
    if digest(prepared_path) != (campaign / "prepared.sha256").read_text().strip():
        raise ValueError("Original direct prepared identity changed")
    prepared = json.loads(prepared_path.read_text())
    relay = json.loads(relay_request.read_text())
    if prepared["temperature_K"] != relay["temperature_K"]:
        raise ValueError("Physical temperatures differ")
    windows = []
    for original in prepared["windows"]:
        bias = original["bias"]
        windows.append(
            dict(
                label=original["label"],
                bias=dict(
                    feature_names=["q_acyl", "q_h"],
                    centers_A=[bias["q_acyl_center"], bias["q_h_center"]],
                    force_constants=[bias["k_acyl"], bias["k_h"]],
                    outer_extent_A=bias["outer_extent"],
                ),
            )
        )
    if len(windows) != 63 or len({w["label"] for w in windows}) != 63:
        raise ValueError("Expected 63 unique legacy umbrella states")
    bias_objects = [RelayBiasWindow(**row["bias"]) for row in windows]
    blocks, sources = (
        [],
        {
            str(prepared_path): digest(prepared_path),
            str(relay_request): digest(relay_request),
        },
    )
    topology = campaign / "inputs/system.prmtop"
    sources[str(topology)] = digest(topology)
    for index, original in enumerate(prepared["windows"]):
        for stage in original["stages"]:
            if stage["kind"] != "production":
                continue
            native = campaign / original["label"] / stage["label"]
            derived = augmented / original["label"] / stage["label"]
            receipt_path = derived / "augmentation.json"
            receipt = json.loads(receipt_path.read_text())
            csv_path = derived / "observables.csv"
            if digest(csv_path) != receipt["csv_sha256"]:
                raise ValueError("Augmented direct CSV changed")
            # The augmentation checked native trajectory/CSV identities on its
            # owning host. Recheck all compact originals on the analysis host.
            for filename in (
                "process.json",
                "acceptance.json",
                "observables/observables.json",
                "observables/observables.csv",
            ):
                path = native / filename
                expected = receipt["source_sha256"][str(path)]
                if digest(path) != expected:
                    raise ValueError(
                        "Compact source changed since augmentation: " + str(path)
                    )
            if (native / "umbrella.RST").read_text().strip() != amber_restraints(
                HarmonicWindow(**original["bias"]), prepared["roles_1based"]
            ).strip():
                raise ValueError("Recorded native restraint differs from declared bias")
            with csv_path.open() as stream:
                rows = list(csv.DictReader(stream))
            values = np.array(
                [
                    [
                        float(row[key])
                        for key in ["q_acyl", "q_act", "q_del", *DIAGNOSTICS]
                    ]
                    for row in rows
                ]
            )
            features = dict(
                zip(("q_acyl", "q_act", "q_del"), values[:, :3].T, strict=True)
            )
            features["q_h"] = values[:, 1] + values[:, 2]
            cross = cross_bias_energies(features, bias_objects)
            npz_path = native / "samples_and_cross_bias.npz"
            with np.load(npz_path, allow_pickle=False) as stored:
                if not np.allclose(
                    stored["cv_values"],
                    np.column_stack([values[:, 0], features["q_h"]]),
                    atol=1e-9,
                    rtol=0,
                ):
                    raise ValueError(
                        "Original direct CV archive differs from coordinates"
                    )
                if not np.allclose(stored["bias_kcal_mol"], cross, atol=1e-7, rtol=0):
                    raise ValueError("All-window direct cross-bias parity failed")
                if not np.allclose(
                    stored["u_kn_bias_only"],
                    cross / (KB * prepared["temperature_K"]),
                    atol=1e-7,
                    rtol=0,
                ):
                    raise ValueError(
                        "Original reduced-bias units or temperature differ"
                    )
            sources[str(receipt_path)] = digest(receipt_path)
            sources[str(npz_path)] = digest(npz_path)
            blocks.append(
                dict(
                    window=original["label"],
                    values=values,
                    original_rows=rows,
                    sampled_duration_ps=5.0,
                    augmentation_receipt=str(receipt_path),
                    native_bias_kcal_mol=cross[index],
                )
            )
    path, reference_sources = direct_reference(prepared)
    sources.update(reference_sources)
    # Preserve the relay histogram grid and extend only its outermost edges.
    # Terminal-tangent fluctuations must not be clipped or silently censored.
    base_edges = np.asarray(relay["s_bin_edges"], dtype=float)
    width = np.diff(base_edges)[0]
    if not np.allclose(np.diff(base_edges), width, rtol=0, atol=1e-12):
        raise ValueError("Shared histogram edges must be equally spaced")
    projection = project_path(
        np.concatenate([b["values"][:, :3] for b in blocks]), path
    )
    lower, upper = float(projection["s"].min()), float(projection["s"].max())
    left = min(0, int(np.floor((lower - base_edges[0]) / width)))
    right = max(len(base_edges) - 1, int(np.floor((upper - base_edges[0]) / width)) + 1)
    edges = base_edges[0] + np.arange(left, right + 1) * width
    config = dict(
        windows=windows,
        temperature_K=prepared["temperature_K"],
        path_3cv_A=path.tolist(),
        s_bin_edges=edges.tolist(),
        relay_base_bin_edges=base_edges.tolist(),
        bin_extension=dict(
            rule="Extend outer edges only on the unchanged grid to retain all frames",
            raw_s_range=[lower, upper],
            added_left_bins=-left,
            added_right_bins=right - len(base_edges) + 1,
            common_bin_width=width,
            frames_censored=0,
        ),
        es_search_bounds=relay["es_search_bounds"],
        acyl_search_bounds=relay["acyl_search_bounds"],
        fixed_reference_s=0.0125,
        diagnostic_columns=DIAGNOSTICS,
    )
    return config, blocks, sources
