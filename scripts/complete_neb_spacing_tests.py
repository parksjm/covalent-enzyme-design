"""Finish manifests from retained spacing inputs without overwriting failed outputs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
from pathlib import Path
import sys

import numpy as np
from parmed.amber import Rst7

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from nyl12_covalent.neb_spacing import predict_spacing_band, spacing_controls  # noqa:E402
from nyl12_covalent.neb_spacing_audit import spacing_image_audit  # noqa:E402
from nyl12_covalent.native_sander_cli import source_sha256  # noqa:E402
from nyl12_covalent.relay_case import load_case, write_json  # noqa:E402
from scripts.prepare_dimer_continuation import write_restart  # noqa:E402


def complete(output: Path) -> None:
    """Continue retained source/controls; midpoint diagnostics never gate actual images."""
    source_folder = output / "source_cycle268"
    source = json.loads((source_folder / "source_manifest.json").read_text())
    cfg = json.loads((source_folder / "original_neb_config.json").read_text())
    if (
        hashlib.sha256(
            gzip.decompress((source_folder / "source_nebpath.xyz.gz").read_bytes())
        ).hexdigest()
        != source["source_xyz_sha256"]
    ):
        raise ValueError("Retained XYZ source payload mismatch")
    model = load_case(Path(cfg["case_path"]))
    reference = Rst7.open(str(source_folder / "native_frozen_reference.ncrst"))
    frozen = np.asarray(reference.coordinates).reshape(-1, 3)
    active = np.flatnonzero(model.active)
    lookup = {int(atom): i for i, atom in enumerate(active)}
    qm = np.array([lookup[int(i)] for i in model.qm])
    frames = []
    for i in range(1, 25):
        xyz = np.asarray(
            Rst7.open(str(source_folder / f"image_{i:02d}.ncrst")).coordinates
        ).reshape(-1, 3)
        if not np.array_equal(xyz[~model.active], frozen[~model.active]):
            raise ValueError("Retained source frozen exterior mismatch")
        frames.append(xyz[active])
    frames = np.stack(frames)
    source_reference = frozen.copy()
    source_reference[active] = frames[0]

    def audit(band: np.ndarray) -> dict:
        """Collect strict actual-image flags and explicit reactive/midpoint diagnostics."""
        rows = []
        ids = {key: lookup[value - 1] for key, value in model.atom_map.items()}
        for i, points in enumerate(band):
            xyz = frozen.copy()
            xyz[active] = points
            row = spacing_image_audit(model, xyz, source_reference)
            c, o, og, n = (
                points[ids[k]] for k in ("carbon", "oxygen", "thr_og", "leaving_n")
            )
            a, b = og - c, o - c
            row.update(
                image=i + 1,
                OG1_C_A=float(np.linalg.norm(og - c)),
                C_Nleave_A=float(np.linalg.norm(c - n)),
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
            strict_H_gate_passed=all(r["strict_H_gate_passed"] for r in rows),
            qm_adjacent_arcs_A=np.linalg.norm(
                np.diff(band[:, qm], axis=0).reshape(len(band) - 1, -1), axis=1
            ).tolist(),
            active_adjacent_arcs_A=np.linalg.norm(
                np.diff(band, axis=0).reshape(len(band) - 1, -1), axis=1
            ).tolist(),
        )

    for name, density in (("uniform", 1.0), ("adaptive", 2.0)):
        folder = output / name
        folder.mkdir(exist_ok=True)
        computed, sampling = spacing_controls(frames, qm, focus_weight=density)
        images = []
        controls = []
        for i, points in enumerate(computed):
            path = folder / f"control_{i:02d}.ncrst"
            if path.exists():
                xyz = np.asarray(Rst7.open(str(path)).coordinates).reshape(-1, 3)
                if not np.array_equal(xyz[active], points):
                    raise ValueError(
                        "Retained controls differ from deterministic source sampling"
                    )
            else:
                xyz = frozen.copy()
                xyz[active] = points
                write_restart(path, xyz, reference)
            controls.append(xyz[active])
            images.append(str(path.resolve()))
        controls = np.stack(controls)
        control_audit = audit(controls)
        predictions = {}
        counts = (24, 48) if name == "uniform" else (24,)
        for count in counts:
            band, report = predict_spacing_band(controls, qm, count)
            report["audit"] = audit(band)
            report["midpoint_audit"] = audit((band[:-1] + band[1:]) / 2)
            report["preflight_passed"] = bool(
                report["endpoint_within_1e_minus6_A"]
                and not report["dropped_controls_zero_based"]
                and report["audit"]["strict_H_gate_passed"]
            )
            predictions[str(count)] = report
            manifest = dict(
                images=images,
                source_cycle=268,
                source_snapshot=source["snapshot"],
                source_manifest=str((source_folder / "source_manifest.json").resolve()),
                sampling=sampling,
                endpoint_source_exact=bool(
                    np.array_equal(controls[0], frames[0])
                    and np.array_equal(controls[-1], frames[-1])
                ),
                endpoint_tolerance_A=1e-6,
                predictions=dict(predictions),
                audits=control_audit,
                image_sha256={p: source_sha256(Path(p)) for p in images},
                preflight_passed=predictions["24"]["preflight_passed"],
                midpoint_diagnostic_not_actual_image_gate=True,
                fixed_exterior_exact=True,
                real_QM_H_labels_preserved=True,
                source_is_optimized_stationary_band=False,
            )
            target = folder / (
                "manifest.json" if count == 24 else "manifest_24_48.json"
            )
            write_json(target, manifest)
            logging.info(
                "READY%s%d gate=%s endpointdrift=%g manifest=%s",
                name,
                count,
                report["preflight_passed"],
                report["endpoint_max_drift_A"],
                target,
            )
            for kind in ("audit", "midpoint_audit"):
                for row in report[kind]["images"]:
                    if row["offending_hydrogens"]:
                        logging.info(
                            "%s%d %s%d offenders=%s",
                            name,
                            count,
                            kind,
                            row["image"],
                            row["offending_hydrogens"],
                        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    logging.basicConfig(level=logging.INFO)
    complete(parser.parse_args().inputs.resolve())
