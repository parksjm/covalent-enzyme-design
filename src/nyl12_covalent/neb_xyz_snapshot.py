"""Reconstruct coherent NEB full-system restarts from native active-atom XYZ."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask
from parmed.periodic_table import Element

from nyl12_covalent.neb_progress import NebReport, parse_last_complete_neb_report
from nyl12_covalent.neb_snapshot import (
    _CYCLE,
    _SAVE,
    _require_finished_run,
    _require_safe_writer,
    _signature,
)

NATIVE_DEFECT = (
    "Amber dlfind_module.write_custom_rst_file passes x,x to write_nc_restart "
    "instead of coordinate/x(lcrd); numbered native NetCDF restarts are unusable"
)


def reconstruct_xyz_band(
    payload: bytes,
    topology: AmberParm,
    reference: Rst7,
    active_mask: str,
    expected_images: int,
    expected_active: int,
) -> tuple[list[np.ndarray], list[int]]:
    """Map native Angstrom XYZ frames into the fixed reference atom order.

    Args:
        payload: Immutable native nebpath.xyz bytes.
        topology: Original topology loaded with the fixed reference coordinates.
        reference: Fixed E-S restart providing the frozen exterior and cell.
        active_mask: Exact native DL-FIND mask evaluated on that reference.
        expected_images: Regular band frame count excluding climbing slot.
        expected_active: Expected movable atom count after mass > 0.1 filtering.

    Returns:
        Full-system coordinate arrays and ordered zero-based active indices.

    Raises:
        ValueError: If mapping, frame format, elements or anchor geometry fail.
    """
    if expected_images < 2 or expected_active < 1:
        raise ValueError("Positive active count and at least two images required")
    selected = AmberMask(topology, active_mask).Selection()
    active = [
        i for i, atom in enumerate(topology.atoms) if selected[i] and atom.mass > 0.1
    ]
    if len(active) != expected_active:
        raise ValueError("Fixed-reference active atom count mismatch")
    base = np.asarray(reference.coordinates).reshape(-1, 3)
    if len(base) != len(topology.atoms) or not np.isfinite(base).all():
        raise ValueError("Invalid reference atom count/coordinates")
    expected_elements = []
    for index in active:
        number = topology.atoms[index].atomic_number
        if number <= 0 or number >= len(Element):
            raise ValueError("Active atom has no unambiguous element")
        expected_elements.append(Element[number].upper())
    active_set = set(active)
    anchors = [
        (bond.atom1.idx, bond.atom2.idx)
        for bond in topology.bonds
        if {bond.atom1.name, bond.atom2.name}
        in ({"CA", "CB"}, {"N", "CA"}, {"CA", "C"})
        and (bond.atom1.idx in active_set or bond.atom2.idx in active_set)
    ]
    if not anchors:
        raise ValueError("No nonreactive backbone/CA-CB anchor bonds available")
    if not payload.endswith(b"\n"):
        raise ValueError("Incomplete native XYZ final record")
    lines = payload.decode("ascii").splitlines()
    cursor = 0
    frames = []
    for _ in range(expected_images):
        if cursor + expected_active + 2 > len(lines):
            raise ValueError("Incomplete native XYZ frame")
        if int(lines[cursor].strip()) != expected_active or lines[cursor + 1].strip():
            raise ValueError("XYZ atom count or native blank comment mismatch")
        coordinates = []
        for offset, element in enumerate(expected_elements):
            line = lines[cursor + 2 + offset]
            if len(line) != 47 or line[:2].strip().upper() != element:
                raise ValueError("XYZ element/order or fixed-width column mismatch")
            coordinates.append(
                [float(line[start : start + 15]) for start in (2, 17, 32)]
            )
        xyz = base.copy()
        xyz[active] = coordinates
        if not np.isfinite(xyz).all():
            raise ValueError("Nonfinite reconstructed XYZ coordinates")
        for first, second in anchors:
            distance = float(np.linalg.norm(xyz[first] - xyz[second]))
            baseline = float(np.linalg.norm(base[first] - base[second]))
            if not 1.0 < distance < 2.0 or abs(distance - baseline) > 0.35:
                raise ValueError("Nonreactive backbone/CA-CB anchor geometry failed")
        frames.append(xyz)
        cursor += expected_active + 2
    if any(line.strip() for line in lines[cursor:]):
        raise ValueError("Unexpected extra XYZ frames or trailing records")
    return frames, active


def snapshot_complete_band_from_xyz(
    mdout: Path,
    xyz_path: Path,
    destination: Path,
    config: dict,
    expected_images: int = 16,
    *,
    writer_pid: int | None = None,
    completed_job_status: Path | None = None,
) -> NebReport:
    """Snapshot a proven complete native XYZ band into validated full NetCDF restarts.

    Args:
        mdout: Native Amber output from the verified unbuffered writer.
        xyz_path: Native nebpath.xyz, overwritten after each complete report.
        destination: New immutable snapshot directory; must not already exist.
        config: Campaign topology/reactant/active_mask/active_atoms configuration.
        expected_images: Number of regular path images, excluding CI slot.
        writer_pid: Live serial Sander process with unbuffered output. Omit only
            for the terminal-run path validated by ``completed_job_status``.
        completed_job_status: Terminal job status required when no writer is live.

    Returns:
        Energy report associated with reconstructed full-system snapshots.

    Raises:
        ValueError: If report/write ordering, mapping or geometry is unproven.
        RuntimeError: If live inputs change during capture, or fixed inputs change during reconstruction.
        FileExistsError: If destination exists.
        OSError: If prerequisites cannot be read.
    """
    if writer_pid is None:
        if completed_job_status is None:
            raise ValueError("Offline capture requires terminal job status")
        writer = _require_finished_run(completed_job_status, mdout.parent)
    else:
        if completed_job_status is not None:
            raise ValueError("Supply either a live writer or terminal job status")
        writer = _require_safe_writer(writer_pid)
    reference_path = Path(config["reactant"])
    topology_path = Path(config["topology"])
    paths = [mdout, xyz_path, reference_path, topology_path]
    if completed_job_status is not None:
        paths.append(completed_job_status)
    before = [_signature(path) for path in paths]
    output, xyz_bytes = mdout.read_bytes(), xyz_path.read_bytes()
    terminal_status_sha256 = (
        hashlib.sha256(completed_job_status.read_bytes()).hexdigest()
        if completed_job_status is not None
        else None
    )
    report = parse_last_complete_neb_report(output, expected_images)
    tail = output[report.report_end_byte :]
    cycle = _CYCLE.search(tail)
    if b"NEB Report" in tail or report.cycle is None or cycle is None:
        raise ValueError("Latest XYZ band has no closed native save loop")
    saves = list(_SAVE.finditer(tail[: cycle.start()]))
    if [int(match[1]) for match in saves] != list(range(1, expected_images + 1)):
        raise ValueError("Incomplete ordered XYZ/native save loop")
    # Close the live-file consistency window before slow topology loading.
    # The captured byte strings are immutable; the writer may advance afterward.
    if [_signature(path) for path in paths] != before or mdout.read_bytes() != output:
        raise RuntimeError("XYZ snapshot inputs changed during byte capture")
    if writer_pid is None:
        current_writer = _require_finished_run(completed_job_status, mdout.parent)
    else:
        current_writer = _require_safe_writer(writer_pid)
    if xyz_path.read_bytes() != xyz_bytes or current_writer != writer:
        raise RuntimeError("XYZ data or writer changed during byte capture")
    if [_signature(path) for path in paths] != before:
        raise RuntimeError("XYZ inputs changed before byte capture closed")
    reference = Rst7.open(str(reference_path))
    topology = AmberParm(str(topology_path), xyz=str(reference_path))
    frames, active = reconstruct_xyz_band(
        xyz_bytes,
        topology,
        reference,
        config["active_mask"],
        expected_images,
        config["active_atoms"],
    )
    reference_hash = hashlib.sha256(reference_path.read_bytes()).hexdigest()
    topology_hash = hashlib.sha256(topology_path.read_bytes()).hexdigest()
    if [_signature(path) for path in paths[2:]] != before[2:]:
        raise RuntimeError("Fixed reference or topology changed while reconstructing")
    if writer_pid is None:
        current_writer = _require_finished_run(completed_job_status, mdout.parent)
        if current_writer != writer or [_signature(path) for path in paths] != before:
            raise RuntimeError(
                "Terminal job status or native files changed during reconstruction"
            )
    destination.mkdir(parents=True, exist_ok=False)
    restarts = []
    for image, coordinates in enumerate(frames, 1):
        target = destination / f"optimized.{image:04d}.ncrst"
        restart = Rst7(natom=len(coordinates))
        restart.coordinates = coordinates
        restart.box = reference.box
        restart.write(str(target), netcdf=True)
        check = Rst7.open(str(target)).coordinates.reshape(-1, 3)
        if not np.allclose(check, coordinates, atol=1e-10, rtol=0):
            raise RuntimeError("Reconstructed NetCDF coordinate roundtrip failed")
        restarts.append(
            {
                "image": image,
                "name": target.name,
                "source": str(xyz_path.resolve()),
                "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
            }
        )
    metadata = {
        "report": asdict(report),
        "source": str(mdout.resolve()),
        "snapshot_utc": datetime.now(timezone.utc).isoformat(),
        "writer_pid": writer_pid,
        "writer_proof": writer,
        "terminal_job_status": (
            {
                "path": str(completed_job_status.resolve()),
                "sha256": terminal_status_sha256,
            }
            if completed_job_status is not None
            else None
        ),
        "capture_policy": (
            "stable inputs and active unbuffered writer verified before reconstruction"
            if writer_pid is not None
            else "terminal job status and absence of a live writer verified before reconstruction"
        ),
        "coordinate_association": "closed native XYZ loop; stable snapshot; fixed-reference reconstruction",
        "coordinate_units": "Angstrom; native write_xyz converts Bohr before output",
        "native_restart_defect": NATIVE_DEFECT,
        "reference": {"path": str(reference_path), "sha256": reference_hash},
        "topology": {"path": str(topology_path), "sha256": topology_hash},
        "active_mask": config["active_mask"],
        "active_indices_0based": active,
        "mapped_elements": [
            Element[topology.atoms[index].atomic_number] for index in active
        ],
        "nonreactive_anchor_policy": "All topology-bonded CA-CB, N-CA, CA-C touching active shell: 1-2 Angstrom and <=0.35 Angstrom change from reference",
        "active_filter": "ascending topology indices; active mask and mass > 0.1",
        "natom": reference.natom,
        "restraints": False,
        "source_xyz_sha256": hashlib.sha256(xyz_bytes).hexdigest(),
        "scientific_status": "reconstructed NEB band; independent energy validation required",
        "restarts": restarts,
    }
    (destination / "source_nebpath.xyz").write_bytes(xyz_bytes)
    (destination / "source_output_snapshot.mdout").write_bytes(output)
    (destination / "band_snapshot.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    return report
