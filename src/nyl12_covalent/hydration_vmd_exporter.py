"""Build and atomically refresh the stripped hydration VMD bundle."""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from parmed.amber import AmberParm, NetCDFTraj

from .hydration_vmd_export_io import (
    ATOM_MAP,
    QM_RANGES,
    accepted_extensions,
    cpptraj_input,
    dcd_dimensions,
    sha,
    selected_indices,
    topology_cpptraj_input,
    trajectory_times,
    water_residues,
)
from .hydration_vmd_viewer import viewer_text
from .relay_progress_vmd import substrate_indices
from .vmd_chain_groups import display_chain_groups
from .vmd_movie_policy import mapped_cc_bonds

LOGGER = logging.getLogger(__name__)


def _build_bundle(
    source: Path, output: Path, cpptraj: Path, summary: Path
) -> dict[str, Any]:
    """Build and verify one fresh bundle in a new staging directory."""
    topology = source / "inputs/system.prmtop"
    restart = source / "inputs/endpoint.ncrst"
    heat, equil = source / "heat/trajectory.nc", source / "equil/trajectory.nc"
    extensions = accepted_extensions(source)
    parm = AmberParm(str(topology))
    waters = water_residues(summary, parm)
    kept = selected_indices(parm, waters)
    output.mkdir(parents=True)
    topology_control = output / "topology.cpptraj"
    topology_control.write_text(
        topology_cpptraj_input(topology, waters, output / "stripped.parm7"),
        encoding="utf-8",
    )
    control = output / "trajectory.cpptraj"
    control.write_text(
        cpptraj_input(
            topology,
            restart,
            heat,
            equil,
            waters,
            extensions,
            output / "hydration.dcd",
        ),
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(cpptraj.parent.parent / "lib")
    for name, input_path in (
        ("topology_cpptraj.log", topology_control),
        ("trajectory_cpptraj.log", control),
    ):
        with (output / name).open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                [str(cpptraj), "-i", str(input_path)],
                check=False,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                cwd=output,
            )
        if completed.returncode:
            tail = (output / name).read_text(encoding="utf-8", errors="replace")[-5000:]
            raise RuntimeError(
                f"cpptraj failed for {name} with exit {completed.returncode}:\n{tail}"
            )
    frame_times = [0.0]
    stage_times: dict[str, list[float]] = {"endpoint": [0.0]}
    for path in (heat, equil, *(segment / "trajectory.nc" for segment in extensions)):
        trajectory = NetCDFTraj.open_old(str(path))
        dimensions = trajectory._ncfile.variables["coordinates"].shape
        if dimensions[1:] != (len(parm.atoms), 3):
            raise ValueError(f"Unexpected trajectory atom dimensions: {path}")
        trajectory.close()
        times = trajectory_times(path)
        if len(times) != dimensions[0]:
            raise ValueError(
                f"Trajectory time and coordinate frame counts differ: {path}"
            )
        offset = 0.0
        if frame_times and times[0] <= frame_times[-1]:
            offset = frame_times[-1]
        shifted = times + offset
        if shifted[0] <= frame_times[-1]:
            raise ValueError(
                f"Trajectory segment does not continue chronologically: {path}"
            )
        frame_times.extend(float(value) for value in shifted)
        stage_times[path.parent.name] = [float(value) for value in shifted]
    nframes = len(frame_times)
    exported = AmberParm(str(output / "stripped.parm7"))
    if len(exported.atoms) != len(kept):
        raise ValueError(
            "cpptraj stripped topology atom count does not match selection"
        )
    actual_frames, actual_atoms = dcd_dimensions(output / "hydration.dcd")
    if actual_frames != nframes or actual_atoms != len(kept):
        raise ValueError(
            f"DCD has {(actual_frames, actual_atoms)}; expected {(nframes, len(kept))}"
        )
    mapping = {int(old): new for new, old in enumerate(kept)}
    qm_old = [idx for lo, hi in QM_RANGES for idx in range(lo - 1, hi)]
    qm_new = [mapping[idx] for idx in qm_old]
    water_atoms_old = []
    for residue_id in waters:
        water_atoms_old.extend(atom.idx for atom in parm.residues[residue_id - 1].atoms)
    water_atoms_new = [mapping[idx] for idx in water_atoms_old]
    atom_indices = {key: mapping[value - 1] for key, value in ATOM_MAP.items()}
    for key, value in {"Phe113N": 1638, "Asn186ND2": 2732}.items():
        atom_indices[key] = mapping[value - 1]
    candidate_o = next(
        atom.idx
        for atom in parm.residues[waters[0] - 1].atoms
        if atom.atomic_number == 8
    )
    atom_pairs = {
        "OG1_C199": (atom_indices["OG1"], atom_indices["C199"]),
        "C199_Nleave": (atom_indices["C199"], atom_indices["Nleave"]),
        "OG1_HG1": (atom_indices["OG1"], atom_indices["HG1"]),
        "Nalpha_HG1": (atom_indices["Nalpha"], atom_indices["HG1"]),
        "Nleave_HG1": (atom_indices["Nleave"], atom_indices["HG1"]),
        "O199_Phe113N": (atom_indices["O199"], atom_indices["Phe113N"]),
        "O199_Asn186ND2": (atom_indices["O199"], atom_indices["Asn186ND2"]),
        "candidate_O_C199": (mapping[candidate_o], atom_indices["C199"]),
        "candidate_O_Nalpha": (mapping[candidate_o], atom_indices["Nalpha"]),
    }
    bd = tuple(mapping[value - 1] for value in (13916, 19916, 19917))
    center_atoms = [atom_indices["OG1"], atom_indices["C199"], *water_atoms_new]
    (output / "view_hydration.vmd").write_text(
        viewer_text(
            nframes,
            len(kept),
            qm_new,
            water_atoms_new,
            atom_indices["HG1"],
            atom_pairs,
            bd,
            center_atoms,
            substrate_indices=[
                mapping[index]
                for index in sorted(
                    substrate_indices(
                        parm, ATOM_MAP["C199"] - 1, ATOM_MAP["Nleave"] - 1
                    )
                )
            ],
            cc_bonds=mapped_cc_bonds(parm, mapping),
            chain_groups=display_chain_groups(
                parm, {old + 1: new for old, new in mapping.items()}
            ),
        ),
        encoding="utf-8",
    )
    (output / "frame_manifest.json").write_text(
        json.dumps(
            {
                "frame_count": nframes,
                "frame_order": [
                    "inputs/endpoint.ncrst (0 ps)",
                    "heat/trajectory.nc",
                    "equil/trajectory.nc",
                    *(f"{segment.name}/trajectory.nc" for segment in extensions),
                ],
                "source_time_ps": {
                    "units": "picosecond",
                    "by_stage": stage_times,
                    "all_frames_chronological": True,
                },
                "frame_times_ps": frame_times,
                "atom_count": len(kept),
                "stripped_to_original_atom_indices_0based": kept.tolist(),
                "fixed_water_residues_1based": waters,
                "accepted_extension_segments": [segment.name for segment in extensions],
                "qm_atom_count": len(qm_new),
                "source_sha256": {
                    str(path.relative_to(source)): sha(path)
                    for path in (
                        topology,
                        restart,
                        heat,
                        equil,
                        *(segment / "trajectory.nc" for segment in extensions),
                    )
                },
                "water_summary_sha256": sha(summary),
                "scope": "hydration reconnaissance; no reaction or converged sampling",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    candidates = json.loads(summary.read_text(encoding="utf-8"))["waters"][:5]
    (output / "water_candidates.json").write_text(
        json.dumps({"schema_version": 1, "waters": candidates}, indent=2) + "\n",
        encoding="utf-8",
    )
    (output / "README.md").write_text(
        "# Deacylation hydration movie\n\n"
        "From this directory run `vmd -e view_hydration.vmd`. The bundle is self-contained. "
        "It displays protein NewCartoon hidden by default, the PA6 substrate, physical QM atoms, "
        "and five fixed candidate waters. Bulk solvent and ions are stripped. Frames are the "
        "initial endpoint followed by all heating and equilibration frames in chronological order. "
        "This is a hydration reconnaissance, not reaction sampling or a converged ensemble. "
        "The initial endpoint, heating, equilibration, and accepted extension segments are "
        "shown in chronological order.\n",
        encoding="utf-8",
    )
    payload = sorted(
        path
        for path in output.iterdir()
        if path.is_file() and path.name != "MANIFEST.sha256"
    )
    (output / "MANIFEST.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in payload), encoding="utf-8"
    )
    return {"frames": nframes, "atoms": len(kept), "waters": waters}


def _publish_bundle(staged: Path, output: Path) -> None:
    """Atomically publish a verified staged bundle and retain prior manifests."""
    if output.exists():
        history_path = output / "manifest_history.json"
        history = (
            json.loads(history_path.read_text(encoding="utf-8"))
            if history_path.exists()
            else []
        )
        history.append(
            {
                "previous_manifest_sha256": sha(output / "MANIFEST.sha256"),
                "previous_manifest": (output / "MANIFEST.sha256").read_text(
                    encoding="utf-8"
                ),
            }
        )
        (staged / "manifest_history.json").write_text(
            json.dumps(history, indent=2) + "\n", encoding="utf-8"
        )
        payload = sorted(
            path
            for path in staged.iterdir()
            if path.is_file() and path.name != "MANIFEST.sha256"
        )
        (staged / "MANIFEST.sha256").write_text(
            "".join(f"{sha(path)}  {path.name}\n" for path in payload),
            encoding="utf-8",
        )
        previous = output.with_name(f".{output.name}.previous")
        if previous.exists():
            raise FileExistsError(f"Unresolved prior bundle backup exists: {previous}")
        output.rename(previous)
        try:
            staged.rename(output)
        except Exception:
            previous.rename(output)
            raise
        clutter = output.parent / "clutter"
        clutter.mkdir(parents=True, exist_ok=True)
        retained = (
            clutter / f"{output.name}_previous_{sha(previous / 'MANIFEST.sha256')[:12]}"
        )
        if retained.exists():
            raise FileExistsError(f"Prior bundle archive already exists: {retained}")
        previous.rename(retained)
    else:
        staged.rename(output)


def export(source: Path, output: Path, cpptraj: Path, summary: Path) -> dict[str, Any]:
    """Refresh the portable VMD bundle atomically from current accepted segments."""
    source, output, cpptraj, summary = (
        path.resolve() for path in (source, output, cpptraj, summary)
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix=f".{output.name}.stage-", dir=output.parent))
    staged = temp / output.name
    try:
        result = _build_bundle(source, staged, cpptraj, summary)
        _publish_bundle(staged, output)
    except Exception as exc:
        failure_dir = source / "clutter" / f"vmd_export_failure_{temp.name.lstrip('.')}"
        failure_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(temp), str(failure_dir))
        raise RuntimeError(f"Failed export evidence retained at {failure_dir}") from exc
    temp.rmdir()
    return result
