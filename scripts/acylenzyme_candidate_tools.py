"""Export coherent NEB bands and execute unbiased candidate trial releases."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

from scripts.export_direct_pt_progress import (
    STRIP_MASK,
    dcd_counts,
    run_cpptraj,
    vmd_script,
    write_checksum_manifest,
)
from scripts.run_acylenzyme_neb import (
    geometry,
    input_text,
    native_converged,
    run_engine,
    write_json,
)


def export_band(
    config: dict,
    snapshot: Path,
    output: Path,
    cycle: int,
    *,
    generation_name: str | None = None,
    description: str | None = None,
) -> Path:
    """Export the validated full-system band with dynamic reactive VMD bonds.

    Args:
        config: Campaign configuration with topology and Sander paths.
        snapshot: Immutable coherent band snapshot directory.
        output: Progress publication directory.
        cycle: Numeric native cycle or saved-history frame identifier.
        generation_name: Optional explicit generation directory name.
        description: Optional label embedded in the VMD playback script.

    Returns:
        Immutable export generation; latest points to this generation.
    """
    manifest = json.loads((snapshot / "band_snapshot.json").read_text())
    generation = output / (generation_name or f"cycle_{cycle:06d}")
    generation.mkdir(parents=True, exist_ok=False)
    retained = generation / "source_restarts"
    retained.mkdir()
    paths = []
    for row in manifest["restarts"]:
        target = retained / row["name"]
        shutil.copyfile(snapshot / row["name"], target)
        if hashlib.sha256(target.read_bytes()).hexdigest() != row["sha256"]:
            raise ValueError(
                "Retained export input hash does not match coherent snapshot"
            )
        paths.append(target)
    topology = Path(config["topology"])
    cpptraj = Path(config["sander"]).parent / "cpptraj"
    trajectory = generation / "trajectory.cpptraj"
    trajectory.write_text(
        "\n".join(
            [f"parm {topology}"]
            + [f"trajin {path.resolve()}" for path in paths]
            + [f"strip {STRIP_MASK}", "trajout reaction_path.dcd dcd", "run", "quit"]
        )
        + "\n"
    )
    top_input = generation / "topology.cpptraj"
    top_input.write_text(
        f"parm {topology}\nparmstrip {STRIP_MASK}\nparmwrite out stripped.parm7\nquit\n"
    )
    run_cpptraj(cpptraj, trajectory, generation, generation / "trajectory_export.log")
    run_cpptraj(cpptraj, top_input, generation, generation / "topology_export.log")
    frames, atoms = dcd_counts(generation / "reaction_path.dcd")
    if frames != len(paths) or atoms < 19929:
        raise ValueError("Exported NEB band has invalid frame/atom count")
    script = vmd_script(
        frames, description=description or f"provisional NEB cycle {cycle}"
    )
    script += "\n# Leaving amine connectivity changes dynamically.\n"
    script += "remove_reactive_bond $molid 19915 19928\n"
    script += "remove_reactive_bond $molid 13916 19928\n"
    for first, second, cutoff in [(19915, 19928, 1.90), (13916, 19928, 1.40)]:
        script += (
            f"mol representation DynamicBonds {cutoff:.2f} 0.12 12.0\n"
            "mol color Name\n"
            f'mol selection "index {first} {second}"\n'
            "mol material Opaque\nmol addrep $molid\n"
            f"label add Bonds $molid/{first} $molid/{second}\n"
        )
    (generation / "reaction_path.vmd").write_text(script)
    shutil.copyfile(snapshot / "band_snapshot.json", generation / "band_snapshot.json")
    write_json(
        generation / "export_provenance.json",
        {
            "coordinate_source": "validated coherent full-system band",
            "archived_coordinates": "immutable stripped reaction_path.dcd",
            "full_restarts_retained": "source_restarts/; original per-image hashes in band_snapshot.json",
            "cpptraj_input_paths_status": "retained full-system inputs; commands rerunnable on original host",
            "source_identity": "band_snapshot.json records original restart paths and hashes",
        },
    )
    write_checksum_manifest(generation)
    pending = output / ".latest.pending"
    pending.symlink_to(generation.name)
    pending.replace(output / "latest")
    return generation


def prepare_candidate(snapshot: Path, campaign: Path, image: int) -> Path:
    """Retain one matched restart and its original report/hash manifest.

    Args:
        snapshot: Validated coherent full band.
        campaign: Campaign root.
        image: Regular band image number.

    Returns:
        New per-image candidate directory; existing candidates are rejected.
    """
    metadata = json.loads((snapshot / "band_snapshot.json").read_text())
    row = next(row for row in metadata["restarts"] if row["image"] == image)
    data = (snapshot / row["name"]).read_bytes()
    if hashlib.sha256(data).hexdigest() != row["sha256"]:
        raise ValueError("Candidate restart does not match validated snapshot hash")
    candidate = campaign / "candidate_minima" / f"image_{image:02d}"
    candidate.mkdir(parents=True, exist_ok=False)
    (candidate / "candidate.ncrst").write_bytes(data)
    metadata.update(
        candidate_image=image,
        candidate_sha256=row["sha256"],
        scientific_status="provisional dip; unbiased trial pending",
    )
    write_json(candidate / "candidate.json", metadata)
    return candidate


def launch_candidate(config_path: Path, candidate: Path, journal: Path) -> None:
    """Launch one independent logged tmux release through the stable job wrapper.

    Args:
        config_path: Original campaign configuration.
        candidate: Retained restart and provenance directory.
        journal: Active append-only journal for the independent trial.
    """
    script = Path(__file__).with_name("monitor_acylenzyme_candidates.py").resolve()
    wrapper = script.with_name("run_nyl12_job.py")
    cfg = json.loads(config_path.read_text())
    image = json.loads((candidate / "candidate.json").read_text())["candidate_image"]
    command = [
        sys.executable,
        str(wrapper),
        "launch",
        "--session",
        f"acyl_trial_{candidate.parent.parent.name}_{image:02d}",
        "--output-dir",
        str(candidate / "job"),
        "--cwd",
        str(script.parents[1]),
    ]
    for path in [
        config_path,
        candidate / "candidate.ncrst",
        script,
        Path(__file__),
        script.with_name("run_acylenzyme_neb.py"),
        Path(cfg["topology"]),
        Path(cfg["reactant"]),
        Path(cfg["base_input"]),
    ]:
        command += ["--input", str(path.resolve())]
    command += [
        "--",
        sys.executable,
        str(script),
        "--config",
        str(config_path),
        "--release",
        str(candidate),
        "--journal",
        str(journal.resolve()),
    ]
    subprocess.run(command, check=True)


def release_candidate(config_path: Path, candidate: Path) -> None:
    """Run an unbiased trial and record convergence and geometry without gating chemistry.

    Args:
        config_path: Campaign configuration preserving method and reference.
        candidate: Immutable candidate input/provenance directory.
    """
    cfg = json.loads(config_path.read_text())
    metadata = json.loads((candidate / "candidate.json").read_text())
    start = candidate / "candidate.ncrst"
    if hashlib.sha256(start.read_bytes()).hexdigest() != metadata["candidate_sha256"]:
        raise ValueError("Candidate input hash changed")
    folder = candidate / "unrestrained_release"
    text = input_text(Path(cfg["base_input"]).read_text(), Path(cfg["reactant"]))
    result = run_engine(cfg, folder, start, text)
    write_json(
        candidate / "trial_acceptance.json",
        {
            "native_converged": native_converged(folder / "minimized.mdout"),
            "initial_geometry_A": geometry(start),
            "final_geometry_A": geometry(result),
            "reaction_restraints": False,
            "acylenzyme_gate_applied": False,
            "stationary_point_claim": False,
            "source_candidate_sha256": metadata["candidate_sha256"],
            "optimized_restart_sha256": hashlib.sha256(result.read_bytes()).hexdigest(),
        },
    )
