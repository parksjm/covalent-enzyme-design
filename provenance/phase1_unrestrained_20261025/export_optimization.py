"""Publish bounded-memory phase1 optimization trajectories for VMD."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
RUN = (
    ROOT
    / "output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2/phase1_unrestrained"
)
PUBLISH = ROOT / "qmmm/pa6_asp272/phase1_unrestrained/vmd"
MASK = "!:1-1313,3328,10325,20216"
ATOM_COUNT = 19976
CPPTRAJ = Path("/data/jpc/code/ambertools26/venv/bin/cpptraj")
BOX_LABELS = ("x", "y", "z", "alpha", "beta", "gamma")


def sha256(path: Path) -> str:
    """Return the SHA256 digest of a file using bounded memory."""
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atom_signature(parm: Any) -> list[tuple[str, str]]:
    """Return ordered atom and residue names for topology comparison."""
    return [(atom.name, atom.residue.name) for atom in parm.atoms]


def stable_copy(source: Path, destination: Path) -> str:
    """Copy a live source only if its identity and size remain stable."""
    before = source.stat()
    digest = hashlib.sha256()
    with source.open("rb") as incoming, destination.open("wb") as outgoing:
        for block in iter(lambda: incoming.read(4 * 1024 * 1024), b""):
            outgoing.write(block)
            digest.update(block)
    after = source.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ino) != (
        after.st_size,
        after.st_mtime_ns,
        after.st_ino,
    ):
        destination.unlink()
        raise RuntimeError(f"Source changed while snapshotting: {source}")
    if destination.stat().st_size != after.st_size:
        raise RuntimeError(f"Snapshot size mismatch: {source}")
    return digest.hexdigest()


def _counts(path: Path) -> tuple[int, int]:
    """Use the shared bounded-memory DCD header validator."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from nyl12_covalent.path_inspection import dcd_counts

    frames, atoms = dcd_counts(path)
    if atoms != ATOM_COUNT:
        raise ValueError(f"Unexpected DCD atom count: {atoms}")
    return frames, atoms


def run_cpptraj(script: Path, log: Path) -> None:
    """Run cpptraj with the project AmberTools runtime environment."""
    env = os.environ.copy()
    env["AMBERHOME"] = "/data/jpc/code/ambertools26/venv"
    env["LD_LIBRARY_PATH"] = "/data/jpc/code/ambertools26/venv/lib"
    with log.open("w") as output:
        subprocess.run(
            [str(CPPTRAJ), "-i", str(script)],
            stdout=output,
            stderr=subprocess.STDOUT,
            check=True,
            env=env,
        )


def publish_once(
    preparation_path: Path,
    layout_path: Path,
    viewer_helper: Path,
    reactive_helper: Path,
) -> dict[str, Any]:
    """Publish one stable snapshot of initial, trajectory, and final frames."""
    preparation = json.loads(preparation_path.read_text())
    layout = json.loads(layout_path.read_text())
    topology = Path(preparation["topology"])
    starting = Path(preparation["start_coordinates"])
    trajectory = RUN / "optimization.nc"
    final = RUN / "optimized.ncrst"
    process_path = RUN / "process.json"
    process = json.loads(process_path.read_text()) if process_path.is_file() else {}
    final_included = final.is_file() and process.get("exit_status") == 0
    for required in (
        preparation_path,
        topology,
        starting,
        trajectory,
        viewer_helper,
        reactive_helper,
    ):
        if not required.is_file():
            raise FileNotFoundError(required)
    box = preparation["box"]
    if len(box) != 6:
        raise ValueError("Preparation manifest must provide six fixed box values")
    PUBLISH.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="phase1-vmd-", dir=RUN) as scratch_name:
        scratch = Path(scratch_name)
        trajectory_copy = scratch / "optimization.nc"
        source_hashes = {
            "preparation": sha256(preparation_path),
            "topology": sha256(topology),
            "starting_coordinates": sha256(starting),
            "trajectory": stable_copy(trajectory, trajectory_copy),
            "publisher": sha256(Path(__file__)),
            "viewer_helper": sha256(viewer_helper),
            "reactive_helper": sha256(reactive_helper),
        }
        sources = [starting, trajectory_copy]
        source_labels = ["starting_frame", "optimization_frames"]
        if final_included:
            final_copy = scratch / "optimized.ncrst"
            source_hashes["optimized_coordinates"] = stable_copy(final, final_copy)
            sources.append(final_copy)
            source_labels.append("exact_final_restart")
        staged_parm = scratch / "stripped.parm7"
        parm_script = scratch / "strip.cpptraj"
        parm_script.write_text(
            f"parm {topology}\nparmstrip {MASK}\nparmwrite out {staged_parm}\n"
        )
        run_cpptraj(parm_script, scratch / "strip.log")
        from parmed.amber import AmberParm

        stripped = AmberParm(str(staged_parm))
        if len(stripped.atoms) != ATOM_COUNT:
            raise ValueError(f"Stripped topology has {len(stripped.atoms)} atoms")
        source_to_display = []
        next_index = 0
        kept_residues = set(range(1, 1314)) | {3328, 10325, 20216}
        full = AmberParm(str(topology))
        for atom in full.atoms:
            if atom.residue.idx + 1 in kept_residues:
                source_to_display.append(atom.idx)
                next_index += 1
        if next_index != ATOM_COUNT:
            raise ValueError("Topology residue mask does not produce 19976 atoms")
        reference_path = ROOT / "qmmm/pa6_asp272/vmd/stripped.parm7.gz"
        if reference_path.is_file():
            reference = AmberParm(str(reference_path))
            if _atom_signature(stripped) != _atom_signature(reference):
                raise ValueError(
                    "Stripped atom order differs from canonical VMD topology"
                )
        staged_dcd = scratch / "optimization.dcd"
        commands = [f"parm {topology}"] + [f"trajin {path}" for path in sources]
        commands.extend(
            [
                "box "
                + " ".join(
                    f"{label} {value:.16g}" for label, value in zip(BOX_LABELS, box)
                ),
                f"strip {MASK}",
                f"trajout {staged_dcd} dcd namdcell",
                "run",
            ]
        )
        dcd_script = scratch / "publish.cpptraj"
        dcd_script.write_text("\n".join(commands) + "\n")
        run_cpptraj(dcd_script, scratch / "publish.log")
        shutil.copy2(scratch / "strip.log", PUBLISH / "strip.cpptraj.log")
        shutil.copy2(scratch / "publish.log", PUBLISH / "publish.cpptraj.log")
        frame_count, atom_count = _counts(staged_dcd)
        staged_top = scratch / "stripped.parm7.gz"
        with staged_parm.open("rb") as incoming, staged_top.open("wb") as raw:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw, mtime=0
            ) as compressed:
                shutil.copyfileobj(incoming, compressed)
        if (
            sha256(staged_parm)
            != hashlib.sha256(gzip.decompress(staged_top.read_bytes())).hexdigest()
        ):
            raise ValueError("Compressed stripped topology digest mismatch")
        shutil.copy2(viewer_helper, PUBLISH / viewer_helper.name)
        shutil.copy2(
            ROOT / "qmmm/pa6_asp272/vmd/qm_indices.tcl", PUBLISH / "qm_indices.tcl"
        )
        shutil.copy2(reactive_helper, PUBLISH / reactive_helper.name)
        viewer = PUBLISH / "view_optimization.vmd"
        viewer.write_text(
            "set ::nyl12_dir [file dirname [file normalize [info script]]]\n"
            "foreach name {qm_indices.tcl nyl12_es_vmd_common.tcl} {source [file join $::nyl12_dir $name]}\n"
            f"nyl12_load stripped.parm7 optimization.dcd {frame_count} $::nyl12_qm_indices $::nyl12_stripped_atom_count\n"
            'puts "Unrestrained phase1 optimization: first frame is the exact input restart; final frame is the exact endpoint when present."\n'
        )

        def unique_pairs(pairs: list[list[str]]) -> list[list[str]]:
            seen: set[frozenset[str]] = set()
            unique: list[list[str]] = []
            for pair in pairs:
                key = frozenset(pair)
                if key not in seen:
                    seen.add(key)
                    unique.append(pair)
            return unique

        layout["stage_measurements"] = {
            "unrestrained": unique_pairs(
                layout["stage_measurements"]["1"] + layout["stage_measurements"]["2"]
            )
        }
        layout["stage_dynamic_bonds"] = {
            "unrestrained": unique_pairs(
                layout["stage_dynamic_bonds"]["1"] + layout["stage_dynamic_bonds"]["2"]
            )
        }
        from nyl12_covalent.vmd_reactive_distances import append_reactive_labels

        append_reactive_labels(
            viewer,
            [{"stage": "unrestrained"}] * frame_count,
            layout,
        )
        for name, staged in (
            ("optimization.dcd", staged_dcd),
            ("stripped.parm7.gz", staged_top),
        ):
            os.replace(staged, PUBLISH / name)
        (PUBLISH / "layout_config.json").write_text(json.dumps(layout, indent=2) + "\n")
        manifest = {
            "source_hashes": source_hashes,
            "frames": frame_count,
            "atoms": atom_count,
            "frame_sources": source_labels,
            "starting_frame_included": True,
            "complete_optimization_frames_included": True,
            "exact_final_restart_included": final_included,
            "final_restart_condition": "Included only after process.json exit_status is 0",
            "optimizer_frames_are_not_physical_time": True,
            "strip_mask": MASK,
            "stripped_atom_count": atom_count,
            "retained_residues_1based": sorted(kept_residues),
            "stripped_to_full_atom_index_0based": source_to_display,
            "viewer": str(viewer),
            "strip_topology_sha256": sha256(staged_parm),
        }
        manifest_path = PUBLISH / "source_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        status = {
            "state": "published",
            "frames": frame_count,
            "atoms": atom_count,
            "manifest": str(manifest_path),
            "trajectory_sha256": sha256(PUBLISH / "optimization.dcd"),
        }
        (PUBLISH / "export_status.json").write_text(json.dumps(status, indent=2) + "\n")
        hashes = [
            f"{sha256(path)}  {path.name}"
            for path in sorted(PUBLISH.iterdir())
            if path.is_file() and path.name != "MANIFEST.sha256"
        ]
        (PUBLISH / "MANIFEST.sha256").write_text("\n".join(hashes) + "\n")
        return status


def main() -> None:
    """Run a one-shot export or watch the active optimization until its restart exists."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation", type=Path, default=RUN / "preparation.json")
    parser.add_argument(
        "--layout",
        type=Path,
        default=ROOT / "configs/qmmm_inspection_layout_20261025.json",
    )
    parser.add_argument(
        "--viewer-helper",
        type=Path,
        default=ROOT / "qmmm/pa6_asp272/vmd/nyl12_es_vmd_common.tcl",
    )
    parser.add_argument(
        "--reactive-helper",
        type=Path,
        default=ROOT / "scripts/nyl12_reactive_distances.tcl",
    )
    parser.add_argument(
        "--watch", action="store_true", help="Republish when source files change"
    )
    parser.add_argument("--poll-seconds", type=float, default=15.0)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    previous: tuple[int, int, int, int] | None = None
    while True:
        trajectory = RUN / "optimization.nc"
        final = RUN / "optimized.ncrst"
        process_path = RUN / "process.json"
        process = json.loads(process_path.read_text()) if process_path.exists() else {}
        signature = (
            trajectory.stat().st_size,
            trajectory.stat().st_mtime_ns,
            final.stat().st_mtime_ns if final.exists() else 0,
            process_path.stat().st_mtime_ns if process_path.exists() else 0,
        )
        if signature != previous:
            try:
                result = publish_once(
                    args.preparation,
                    args.layout,
                    args.viewer_helper,
                    args.reactive_helper,
                )
                logging.info(
                    "Published %d frames with %d atoms",
                    result["frames"],
                    result["atoms"],
                )
            except Exception as error:
                PUBLISH.mkdir(parents=True, exist_ok=True)
                (PUBLISH / "export_status.json").write_text(
                    json.dumps(
                        {
                            "state": "failed",
                            "error_type": type(error).__name__,
                            "error": str(error),
                            "trajectory_bytes": trajectory.stat().st_size
                            if trajectory.exists()
                            else None,
                            "updated_unix": time.time(),
                        },
                        indent=2,
                    )
                    + "\n"
                )
                logging.exception("Trajectory publication failed")
                if not args.watch:
                    raise
            previous = signature
        if not args.watch or process.get("exit_status") is not None:
            break
        time.sleep(args.poll_seconds)


if __name__ == "__main__":
    main()
