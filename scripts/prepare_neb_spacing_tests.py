"""Prepare immutable shared spacing-test coordinates; leave engines untouched."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import shutil
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from nyl12_covalent.neb_spacing import (  # noqa:E402
    load_spacing_source,
    spacing_controls,
    predict_spacing_band,
)  # noqa:E402
from nyl12_covalent.native_sander_cli import source_sha256  # noqa:E402
from nyl12_covalent.relay_case import write_json  # noqa:E402
from nyl12_covalent.relay_neb_trace import proton_audit  # noqa:E402
from scripts.prepare_dimer_continuation import write_restart  # noqa:E402

LOGGER = logging.getLogger(__name__)


def prepare(config_path: Path, snapshot: Path, output: Path) -> dict[str, Any]:
    """Publish shared source band and both bounded13-control sampling variants.

    Args:
        config_path: Original NEB config; never modified.
        snapshot: Immutable accepted monitor update, selected explicitly.
        output: Shared new input root, conventionally04_neb/spacing_tests/inputs.

    Returns:
        Uniform/adaptive manifest paths. No trial config or engine is created.

    Raises:
        ValueError: If source or endpoint/proton gates fail.
        FileExistsError: If any requested source/variant directory already exists.
    """
    config_path, snapshot, output = (
        p.resolve() for p in (config_path, snapshot, output)
    )
    model, reference, source, evidence = load_spacing_source(config_path, snapshot)
    output.mkdir(parents=True, exist_ok=True)
    source_folder = output / f"source_cycle{evidence['cycle']}"
    source_folder.mkdir()
    shutil.copy2(snapshot / "band_snapshot.json", source_folder / "band_snapshot.json")
    shutil.copy2(
        snapshot / "source_nebpath.xyz.gz", source_folder / "source_nebpath.xyz.gz"
    )
    shutil.copy2(config_path, source_folder / "original_neb_config.json")
    shutil.copy2(
        Path(model.config["topology"]), source_folder / "source_topology.prmtop"
    )
    shutil.copy2(
        Path(model.config["source_restart"]),
        source_folder / "case_source_restart.ncrst",
    )
    full_reference = np.asarray(reference.coordinates).reshape(-1, 3)
    write_restart(
        source_folder / "native_frozen_reference.ncrst", full_reference, reference
    )
    active = np.flatnonzero(model.active)
    lookup = {int(atom): i for i, atom in enumerate(active)}
    qm = np.array([lookup[int(atom)] for atom in model.qm])

    def audit(band: np.ndarray) -> dict[str, Any]:
        """Audit all proton labels and local reaction angle/adjacent active gaps."""
        rows = []
        ids = {key: lookup[value - 1] for key, value in model.atom_map.items()}
        for i, points in enumerate(band):
            xyz = full_reference.copy()
            xyz[active] = points
            row = proton_audit(model, xyz)
            c, o, og, n = (
                points[ids[k]] for k in ("carbon", "oxygen", "thr_og", "leaving_n")
            )
            a, b = og - c, o - c
            row.update(
                image=i + 1,
                OG1_C_A=float(np.linalg.norm(og - c)),
                C_Nleave_A=float(np.linalg.norm(c - n)),
                C_O_A=float(np.linalg.norm(c - o)),
                attack_angle_degrees=float(
                    np.degrees(
                        np.arccos(
                            np.clip(
                                np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b),
                                -1,
                                1,
                            )
                        )
                    )
                ),
            )
            rows.append(row)
        return dict(
            images=rows,
            qm_adjacent_arcs_A=np.linalg.norm(
                np.diff(band[:, qm], axis=0).reshape(len(band) - 1, -1), axis=1
            ).tolist(),
            active_adjacent_arcs_A=np.linalg.norm(
                np.diff(band, axis=0).reshape(len(band) - 1, -1), axis=1
            ).tolist(),
            max_adjacent_atom_movement_A=float(
                np.linalg.norm(np.diff(band, axis=0), axis=2).max(initial=0)
            ),
        )

    source_audit = audit(source)
    source_paths = []
    for i, points in enumerate(source):
        xyz = full_reference.copy()
        xyz[active] = points
        path = source_folder / f"image_{i + 1:02d}.ncrst"
        write_restart(path, xyz, reference)
        source_paths.append(str(path))
    evidence.update(
        images=source_paths,
        audit=source_audit,
        source_code_sha256={
            str(p): source_sha256(p)
            for p in (Path(__file__), ROOT / "src/nyl12_covalent/neb_spacing.py")
        },
    )
    write_json(source_folder / "source_manifest.json", evidence)
    results = {}
    for name, density, counts in (("uniform", 1.0, (24, 48)), ("adaptive", 2.0, (24,))):
        folder = output / name
        folder.mkdir()
        controls, sampling = spacing_controls(source, qm, focus_weight=density)
        if not np.array_equal(controls[0], source[0]) or not np.array_equal(
            controls[-1], source[-1]
        ):
            raise ValueError("Control endpoint source identity changed")
        images = []
        for i, points in enumerate(controls):
            xyz = full_reference.copy()
            xyz[active] = points
            path = folder / f"control_{i:02d}.ncrst"
            write_restart(path, xyz, reference)
            images.append(str(path))
        predictions = {}
        for count in counts:
            predicted, metadata = predict_spacing_band(controls, qm, count)
            metadata["audit"] = audit(predicted)
            metadata["midpoint_audit"] = audit((predicted[:-1] + predicted[1:]) / 2)
            metadata["preflight_passed"] = (
                metadata["endpoint_within_1e_minus6_A"]
                and not metadata["dropped_controls_zero_based"]
            )
            predictions[str(count)] = metadata
        manifest = dict(
            images=images,
            source_cycle=evidence["cycle"],
            source_snapshot=str(snapshot),
            source_manifest=str(source_folder / "source_manifest.json"),
            sampling=sampling,
            endpoint_source_exact=True,
            endpoint_tolerance_A=1e-6,
            predictions=predictions,
            audits=audit(controls),
            image_sha256={path: source_sha256(Path(path)) for path in images},
            preflight_passed=all(
                row["preflight_passed"] for row in predictions.values()
            ),
            fixed_exterior_exact=True,
            real_QM_H_labels_preserved=True,
            method="same102QM54H8links/charge-1/all4662active; coordinate-only noengine",
            source_is_optimized_stationary_band=False,
        )
        write_json(folder / "manifest.json", manifest)
        results[name] = str(folder / "manifest.json")
        LOGGER.info(
            "Published%s13 controls; preflight=%s predictions=%s",
            name,
            manifest["preflight_passed"],
            {
                n: dict(slots=r["slots_1based"], drift=r["endpoint_max_drift_A"])
                for n, r in predictions.items()
            },
        )
    write_json(output / "seed_inputs.json", results)
    return results


def main() -> None:
    """Prepare explicitly selected retained source coordinates only."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    prepare(args.config, args.snapshot, args.output)


if __name__ == "__main__":
    main()
