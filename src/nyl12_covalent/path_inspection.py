"""Publish immutable, progressively extended reaction-path inspection bundles."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import struct
import subprocess
from pathlib import Path
from typing import Any

from parmed.amber import AmberParm, Rst7

from nyl12_covalent.path_distances import distance_jumps, measure_distances


def digest(path: Path) -> str:
    """Return the SHA-256 digest of a source artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dcd_counts(path: Path) -> tuple[int, int]:
    """Read and validate a standard little-endian coordinate DCD header."""
    with path.open("rb") as handle:
        if handle.read(4) != struct.pack("<i", 84):
            raise ValueError("Unsupported DCD record format")
        data = handle.read(84)
        if len(data) != 84 or data[:4] != b"CORD":
            raise ValueError("Invalid coordinate DCD header")
        frames = struct.unpack("<i", data[4:8])[0]
        if handle.read(4) != struct.pack("<i", 84):
            raise ValueError("Invalid DCD header trailer")
        title_size = struct.unpack("<i", handle.read(4))[0]
        if title_size < 0:
            raise ValueError("Invalid DCD title")
        handle.seek(title_size + 4, 1)
        if handle.read(4) != struct.pack("<i", 4):
            raise ValueError("Invalid DCD atom-count record")
        atoms = struct.unpack("<i", handle.read(4))[0]
        if frames < 1 or atoms < 1:
            raise ValueError("Empty DCD")
        return frames, atoms


def collect_endpoints(run: Path, config: dict[str, Any]) -> dict[str, Any]:
    """Read finalized accepted windows and an explicitly rejected endpoint.

    Args:
        run: Existing path directory with status.json and inputs.json.
        config: Explicit atom identities, pairs, threshold and export paths.

    Returns:
        Auditable ordered endpoint data, hashes and discontinuity flags.
    """
    status = json.loads((run / "status.json").read_text())
    inputs = json.loads((run / "inputs.json").read_text())
    labels = list(status["completed_windows"])
    if status["state"] == "stopped" and config["include_rejected_endpoint"]:
        rejected = []
        for path in sorted(run.glob("*/acceptance.json")):
            if path.parent.name not in labels and path.parent.name not in config.get(
                "excluded_endpoint_directories", []
            ):
                rejected.append(path.parent.name)
        labels.extend(rejected)
    if not labels:
        raise ValueError("No finalized endpoints available")
    hashes = {str(run / "inputs.json"): digest(run / "inputs.json")}
    topology = Path(inputs["topology"])
    hashes[str(topology)] = digest(topology)
    rows = []
    parm = AmberParm(str(topology))
    for label in labels:
        folder = run / label
        acceptance = json.loads((folder / "acceptance.json").read_text())
        process = json.loads((folder / "process.json").read_text())
        if process["exit_status"] != 0:
            raise ValueError(f"Unfinished or failed engine: {label}")
        restart = folder / "optimized.ncrst"
        coordinates = Rst7.open(str(restart)).coordinates.reshape((len(parm.atoms), 3))
        for key, identity in config["atom_identities"].items():
            atom = parm.atoms[int(key) - 1]
            if (
                atom.name != identity["name"]
                or atom.residue.idx + 1 != identity["residue"]
            ):
                raise ValueError(f"Reaction atom identity mismatch: {key}")
        accepted = label in status["completed_windows"]
        if accepted and not acceptance["native_converged"]:
            raise ValueError(f"Completed window lacks convergence: {label}")
        rows.append(
            {
                "window": label,
                "stage": process.get("stage", "direct"),
                "accepted": accepted,
                "restart": str(restart),
                "distances_A": measure_distances(coordinates, config["pairs"]),
                "acceptance": acceptance,
            }
        )
        for path in (restart, folder / "acceptance.json", folder / "process.json"):
            hashes[str(path)] = digest(path)
    return {
        "rows": rows,
        "source_sha256": hashes,
        "run_state": status,
        "pairs": config["pairs"],
        "box": inputs["box"],
        "jump_threshold_A": config["jump_threshold_A"],
        "jumps": distance_jumps(rows, config["jump_threshold_A"]),
    }


def publish_inspection(run: Path, config_path: Path, plot_script: Path) -> Path:
    """Publish one verified endpoint movie and distance plot per new step.

    Args:
        run: Scientific run to inspect, never modified by this exporter.
        config_path: Explicit local executables, viewer assets and atom mapping.
        plot_script: Matplotlib worker run using the declared Python environment.

    Returns:
        Immutable generation directory; an atomic latest link identifies it.
    """
    config = json.loads(config_path.read_text())
    data = collect_endpoints(run, config)
    for path in (
        config_path,
        plot_script,
        Path(__file__),
        Path(__file__).with_name("path_distances.py"),
    ):
        data["source_sha256"][str(path)] = digest(path)
    output = Path(config["output"])
    output.mkdir(exist_ok=True)
    identity = hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
    latest = output / "latest"
    if latest.exists():
        previous = json.loads((latest / "reaction_distances.json").read_text())
        if previous["identity"] == identity:
            return latest.resolve()
    generation = output / f"update_{len(list(output.glob('update_*'))) + 1:04d}"
    generation.mkdir()
    data["identity"] = identity
    (generation / "reaction_distances.json").write_text(
        json.dumps(data, indent=2) + "\n"
    )
    names = list(config["pairs"])
    with (generation / "reaction_distances.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["window", "stage", "accepted", *names])
        for row in data["rows"]:
            writer.writerow(
                [
                    row["window"],
                    row["stage"],
                    row["accepted"],
                    *[row["distances_A"][name] for name in names],
                ]
            )
    assets = Path(config["viewer_assets"])
    topology = json.loads((run / "inputs.json").read_text())["topology"]
    commands = [f"parm {topology}"]
    commands += [f"trajin {row['restart']}" for row in data["rows"]]
    box = data["box"]
    commands += [
        f"box x {box[0]} y {box[1]} z {box[2]} alpha {box[3]} beta {box[4]} gamma {box[5]}",
        f"strip {config['strip_mask']}",
        "trajout reaction_path.dcd dcd namdcell",
        "run",
        "quit",
    ]
    cpp = generation / "export.cpptraj"
    cpp.write_text("\n".join(commands) + "\n")
    env = os.environ.copy()
    env.update(config["export_environment"])
    with (generation / "export.log").open("w") as handle:
        subprocess.run(
            [config["cpptraj"], "-i", str(cpp)],
            cwd=generation,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )
    expected = (len(data["rows"]), config["expected_stripped_atoms"])
    if dcd_counts(generation / "reaction_path.dcd") != expected:
        raise ValueError("Exported path atom/frame mismatch")
    for name in ("stripped.parm7.gz", "qm_indices.tcl", "nyl12_es_vmd_common.tcl"):
        shutil.copy2(assets / name, generation / name)
        data["source_sha256"][str(assets / name)] = digest(assets / name)
    loader = """set ::nyl12_dir [file dirname [file normalize [info script]]]
foreach name {qm_indices.tcl nyl12_es_vmd_common.tcl} {
    source [file join $::nyl12_dir $name]
}
nyl12_load stripped.parm7 reaction_path.dcd FRAMES $::nyl12_qm_indices $::nyl12_stripped_atom_count
puts "One endpoint per reaction window; inspect reaction_distances.json for acceptance."
""".replace("FRAMES", str(expected[0]))
    (generation / "reaction_path.vmd").write_text(loader)
    with (generation / "plot.log").open("w") as handle:
        subprocess.run(
            [
                config["plot_python"],
                str(plot_script),
                str(generation / "reaction_distances.json"),
            ],
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )
    (generation / "reaction_distances.json").write_text(
        json.dumps(data, indent=2) + "\n"
    )
    manifest = "".join(
        f"{digest(path)}  {path.name}\n"
        for path in sorted(generation.iterdir())
        if path.is_file()
    )
    (generation / "MANIFEST.sha256").write_text(manifest)
    if config.get("download_directory"):
        from nyl12_covalent.path_publish import publish_download

        publish_download(generation, Path(config["download_directory"]))
    temporary = output / ".latest.pending"
    temporary.symlink_to(generation.name)
    temporary.replace(latest)
    return generation
