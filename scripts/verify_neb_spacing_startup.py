"""Verify first append-only native history frames against frozen band predictions."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys
import time
from typing import Any

import numpy as np
from parmed.amber import Rst7
from parmed.periodic_table import Element

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from nyl12_covalent.neb_spacing import predict_spacing_band  # noqa:E402
from nyl12_covalent.native_sander_cli import source_sha256  # noqa:E402
from nyl12_covalent.relay_case import load_case, write_json  # noqa:E402

NATIVE = Path(
    "/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/libdlfind/src/dlf_neb.F90"
)
LOGGER = logging.getLogger(__name__)


class PendingFrameError(RuntimeError):
    """A first native history frame has not yet been written completely."""


def read_first_frame(path: Path, elements: list[str]) -> tuple[np.ndarray, bytes]:
    """Read only the first complete fixed-width XYZ frame and verify stable prefix.

    Args:
        path: Native append-only per-image neb_####.xyz history.
        elements: Expected active4662 element order from fixed original topology.

    Returns:
        Coordinates in angstrom and exact immutable first-frame bytes.

    Raises:
        PendingFrameError: If file is missing or first frame is incomplete.
        ValueError: If atom count, native format, elements or finiteness fail.
        RuntimeError: If inode or first-frame prefix changes during capture.
    """
    if not path.exists():
        raise PendingFrameError(str(path))
    before = path.stat()
    with path.open("rb") as stream:
        lines = [stream.readline() for _ in range(len(elements) + 2)]
    if any(not line.endswith(b"\n") for line in lines):
        raise PendingFrameError(str(path))
    if int(lines[0]) != len(elements) or lines[1].strip():
        raise ValueError("Native first-frame header/order mismatch")
    xyz = []
    for line, element in zip(lines[2:], elements, strict=True):
        record = line.rstrip(b"\r\n").decode("ascii")
        if len(record) != 47 or record[:2].strip().upper() != element.upper():
            raise ValueError("Native element order/fixed-width mapping mismatch")
        xyz.append([float(record[start : start + 15]) for start in (2, 17, 32)])
    coordinates = np.asarray(xyz, dtype=float)
    if not np.isfinite(coordinates).all():
        raise ValueError("Nonfinite first-frame coordinates")
    prefix = b"".join(lines)
    with path.open("rb") as stream:
        if stream.read(len(prefix)) != prefix:
            raise RuntimeError("Native first-frame prefix changed")
    after = path.stat()
    if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
        raise RuntimeError("Native history file was replaced during prefix capture")
    return coordinates, prefix


def verify_trial(
    trial: Path, model: Any, *, max_wait_seconds: float = 120.0
) -> dict[str, Any]:
    """Compare actual first native frames with the frozen source-control predictor.

    Args:
        trial: Existing running trial directory; scientific files are read-only.
        model: Common original validated102QM54H case.
        max_wait_seconds: Bounded wait for initial history completion.

    Returns:
        Immutable startup geometry/mapping evidence, separately from convergence.

    Raises:
        ValueError: If identities/config/frozen input exterior are inconsistent.
        PendingFrameError: If first-frame completeness remains pending at deadline.
        FileExistsError: If verification evidence already exists.
    """
    trial = trial.resolve()
    cfgpath = trial / "config.json"
    cfg = json.loads(cfgpath.read_text())
    if (trial / "startup_verification.json").exists():
        raise FileExistsError(trial / "startup_verification.json")
    active = np.flatnonzero(model.active)
    if (
        cfg["active_atoms"] != 4662
        or cfg["qm_charge"] != -1
        or cfg["qm_atom_ids"] != (model.qm + 1).tolist()
        or cfg["active_atom_ids"] != (active + 1).tolist()
    ):
        raise ValueError("Trial QM/active/charge identities differ from frozen case")
    lookup = {int(atom): i for i, atom in enumerate(active)}
    qm = np.array([lookup[int(atom)] for atom in model.qm])
    h = np.array([lookup[int(atom)] for atom in model.hydrogens])
    elements = [Element[model.topology.atoms[i].atomic_number] for i in active]
    controls = []
    input_hashes = {str(cfgpath): source_sha256(cfgpath)}
    reference = Rst7.open(cfg["reactant"])
    frozen = np.asarray(reference.coordinates).reshape(-1, 3)
    for name in cfg["images"]:
        path = Path(name)
        rst = Rst7.open(str(path))
        xyz = np.asarray(rst.coordinates).reshape(-1, 3)
        if not np.array_equal(
            xyz[~model.active], frozen[~model.active]
        ) or not np.array_equal(rst.box, reference.box):
            raise ValueError("Control frozen exterior/box changed")
        controls.append(xyz[active])
        input_hashes[str(path)] = source_sha256(path)
    controls = np.stack(controls)
    count = cfg["regular_images"]
    predicted, prediction = predict_spacing_band(controls, qm, count)
    deadline = time.monotonic() + max_wait_seconds
    while True:
        try:
            captured = [
                read_first_frame(
                    Path(cfg["calculation"]) / f"neb_{i:04d}.xyz", elements
                )
                for i in range(1, count + 1)
            ]
            break
        except PendingFrameError as error:
            if time.monotonic() >= deadline:
                raise
            LOGGER.info("Explicitstartup_pending %s; bounded retry", error)
            time.sleep(min(10.0, max(0.0, deadline - time.monotonic())))
    actual = np.stack([row[0] for row in captured])
    difference = actual - predicted
    component = float(np.abs(difference).max())
    norm = np.linalg.norm(difference, axis=2)
    endpoints = [
        float(np.linalg.norm(actual[0] - controls[0], axis=1).max()),
        float(np.linalg.norm(actual[-1] - controls[-1], axis=1).max()),
    ]
    passed = bool(component <= 1e-7 and max(endpoints) <= 1e-6)
    output = trial / "startup_verification_sources"
    output.mkdir()
    prefixes = []
    for i, (_, payload) in enumerate(captured, 1):
        target = output / f"neb_{i:04d}.first_frame.xyz.gz"
        with gzip.open(target, "wb") as stream:
            stream.write(payload)
        prefixes.append(
            dict(
                image=i,
                native_history=str(Path(cfg["calculation"]) / f"neb_{i:04d}.xyz"),
                first_frame_sha256=hashlib.sha256(payload).hexdigest(),
                retained_prefix=str(target),
            )
        )
    recipe = {}
    for source in (Path(__file__), ROOT / "src/nyl12_covalent/neb_spacing.py", NATIVE):
        copied = output / source.name
        shutil.copy2(source, copied)
        digest = source_sha256(source)
        if source_sha256(copied) != digest:
            raise ValueError("Startup recipe snapshot mismatch")
        recipe[str(source)] = digest
    mdin = Path(cfg["base_input"])
    input_hashes[str(mdin)] = source_sha256(mdin)
    np.savez_compressed(
        output / "predicted_and_actual_active.npz",
        predicted=predicted,
        actual=actual,
        active_atom_ids=active + 1,
        qm_positions=qm,
        hydrogen_positions=h,
    )
    report = dict(
        status="passed" if passed else "failed",
        passed=passed,
        trial=str(trial),
        regular_images=count,
        authoritative_frame="First complete append-only neb_####.xyz frame, native initialization before optimizer motion",
        native_write_order_evidence=dict(
            source=str(NATIVE),
            lines="377-405; dlf_neb_init writes initialXYZ before gradient/optimizer updates",
        ),
        first_frame_prefixes=prefixes,
        element_order_verified=True,
        active_atoms=4662,
        qm_atoms=102,
        real_H=54,
        qm_charge=-1,
        max_component_error_A=component,
        max_vector_error_A=float(norm.max()),
        qm_max_vector_error_A=float(norm[:, qm].max()),
        H_max_vector_error_A=float(norm[:, h].max()),
        print_precision_decimal_places=7,
        coordinate_component_tolerance_A=1e-7,
        endpoint_vector_errors_A=endpoints,
        endpoint_tolerance_A=1e-6,
        input_frozen_exterior_bitwise=True,
        periodic_box_bitwise=True,
        prediction=prediction,
        input_sha256=input_hashes,
        recipe_sha256=recipe,
        native_numbered_NetCDF_restarts_used=False,
        new_energy_or_force_evaluations=False,
        force_convergence_assessed=False,
        chemical_stationarity_assessed=False,
    )
    write_json(trial / "startup_verification.json", report)
    LOGGER.info(
        "startup%s passed=%s max_component=%g QMmax=%g Hmax=%g endpointmax=%g",
        trial.name,
        passed,
        component,
        report["qm_max_vector_error_A"],
        report["H_max_vector_error_A"],
        max(endpoints),
    )
    return report


def main() -> None:
    """Verify all explicitly selected trials with one common validated model."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trials", type=Path, nargs="+")
    parser.add_argument("--max-wait-seconds", type=float, default=120.0)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    first = json.loads((args.trials[0] / "config.json").read_text())
    model = load_case(Path(first["case_path"]))
    for trial in args.trials:
        verify_trial(trial, model, max_wait_seconds=args.max_wait_seconds)


if __name__ == "__main__":
    main()
