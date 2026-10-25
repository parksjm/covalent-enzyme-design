"""Export fixed-water optimizer display coordinates using installed cpptraj."""

from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
from parmed.amber import AmberParm, NetCDFTraj


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nyl12_covalent.vmd_export_io import sha, validate_dcd, select_display_atoms  # noqa: E402


def export(spec_path: Path, output: Path) -> None:
    """Create a fresh curated package with verified full-solvent DCD coordinates.

    Args:
        spec_path: Explicit source paths, engine configuration and executable.
        output: Fresh output directory.

    Raises:
        ValueError: Source inventory, frame count or export verification fails.
    """
    spec = json.loads(spec_path.read_text())
    paths = {key: Path(value) for key, value in spec.items()}
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    original_hashes = {key: sha(path) for key, path in paths.items()}
    output.mkdir(parents=True, exist_ok=False)
    source = output / "source"
    source.mkdir()
    shutil.copy2(paths["helper"], source / paths["helper"].name)
    for path in (
        Path(__file__),
        Path(__file__).parents[1] / "src/nyl12_covalent/vmd_export_io.py",
        paths["viewer"],
        paths["full_viewer"],
        paths["common"],
        spec_path,
    ):
        shutil.copy2(path, source / path.name)
    topology = AmberParm(str(paths["topology"]))
    if len(topology.atoms) != 103251:
        raise ValueError("Expected original 103251-atom full-solvent topology")
    config = json.loads(paths["engine_config"].read_text())
    mask = config["qm_options"]["qmmask"]
    if not mask.startswith("@"):
        raise ValueError("Explicit numeric QM atom mask required")
    indices = [int(value) - 1 for value in mask[1:].split(",")]
    if (
        len(indices) != 98
        or len(set(indices)) != 98
        or not all(0 <= value < 103251 for value in indices)
    ):
        raise ValueError("Expected exactly 98 distinct physical QM atom indices")
    states, frames = [], []
    for key, kind in [
        ("initial", "original_initial"),
        ("first_trajectory", "first_pilot_saved_state"),
        ("continuation_trajectory", "continuation_saved_state"),
        ("final", "exact_final_endpoint"),
    ]:
        if "trajectory" in key:
            traj = NetCDFTraj.open_old(str(paths[key]))
            coordinates, boxes = np.array(traj.coordinates), np.array(traj.box)
            traj.close()
        else:
            parm = AmberParm(str(paths["topology"]), xyz=str(paths[key]))
            coordinates, boxes = np.array([parm.coordinates]), np.array([parm.box])
        expected_count = {
            "initial": 1,
            "first_trajectory": 25,
            "continuation_trajectory": 21,
            "final": 1,
        }[key]
        if coordinates.shape != (expected_count, 103251, 3) or boxes.shape != (
            expected_count,
            6,
        ):
            raise ValueError(f"Unexpected source frame/box inventory: {key}")
        for index, (xyz, box) in enumerate(zip(coordinates, boxes)):
            if not np.isfinite(xyz).all() or not np.isfinite(box).all():
                raise ValueError("Nonfinite source coordinates or box")
            raw_box = box.copy()
            derived = "trajectory" in key and np.all(raw_box == 0)
            if derived:
                box = states[0][1].copy()
            if states and not np.allclose(box, states[0][1], atol=1e-10, rtol=0):
                raise ValueError("Nonzero source box conflicts with declared fixed box")
            states.append((xyz, box))
            frames.append(
                {
                    "display_frame_index_0based": len(frames),
                    "source_kind": kind,
                    "source_path": str(paths[key]),
                    "source_state_index_1based": index + 1,
                    "periodic_box_metadata": {
                        "derived_from_fixed_reference": bool(derived),
                        "raw_box_values": raw_box.tolist(),
                        "box_reference_path": str(
                            paths["initial"] if derived else paths[key]
                        ),
                        "box_reference_sha256": original_hashes[
                            "initial" if derived else key
                        ],
                        "reason": "Explicit ntb1 fixed-box reconstruction for zero optimizer NetCDF cells; no coordinate transformation"
                        if derived
                        else "Direct source restart box",
                    },
                }
            )
    kept, water_ids = select_display_atoms(topology, states, indices)
    residue_ids = sorted({topology.atoms[int(i)].residue.idx + 1 for i in kept})
    ranges = []
    start = previous = residue_ids[0]
    for value in residue_ids[1:] + [None]:
        if value is not None and value == previous + 1:
            previous = value
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = value
    keep_mask = ":" + ",".join(ranges)
    stripped = output / "stripped.parm7"
    script = source / "stripped_topology.cpptraj"
    script.write_text(
        f"parm {paths['topology']}\nparmstrip !{keep_mask}\nparmwrite out {stripped}\n"
    )
    with (source / "stripped_topology.log").open("w") as log:
        subprocess.run(
            [str(paths["cpptraj"]), "-i", str(script)],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    remap = {int(old): new for new, old in enumerate(kept)}
    stripped_indices = [remap[index] for index in indices]
    for name, keys in [
        (
            "optimization",
            ["initial", "first_trajectory", "continuation_trajectory", "final"],
        ),
        ("initial_final", ["initial", "final"]),
    ]:
        commands = (
            [f"parm {paths['topology']}"]
            + [f"trajin {paths[key]}" for key in keys]
            + ([f"strip !{keep_mask}"] if name == "optimization" else [])
            + [
                "box "
                + " ".join(
                    f"{label} {value:.16g}"
                    for label, value in zip(
                        ("x", "y", "z", "alpha", "beta", "gamma"), states[0][1]
                    )
                )
            ]
            + [f"trajout {output / (name + '.dcd')} dcd namdcell", "run"]
        )
        script = source / (name + ".cpptraj")
        script.write_text("\n".join(commands) + "\n")
        with (source / (name + ".log")).open("w") as log:
            subprocess.run(
                [str(paths["cpptraj"]), "-i", str(script)],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
    validations = {
        "optimization": validate_dcd(
            output / "optimization.dcd", [(xyz[kept], box) for xyz, box in states]
        ),
        "initial_final": validate_dcd(
            output / "initial_final.dcd", [states[0], states[-1]]
        ),
    }
    with (
        paths["topology"].open("rb") as incoming,
        (output / "md_solvated.parm7.gz").open("wb") as raw,
    ):
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as outgoing:
            shutil.copyfileobj(incoming, outgoing)
    with gzip.open(output / "md_solvated.parm7.gz", "rb") as handle:
        if hashlib.sha256(handle.read()).hexdigest() != original_hashes["topology"]:
            raise ValueError("Compressed topology digest mismatch")
    with (
        stripped.open("rb") as incoming,
        (output / "stripped.parm7.gz").open("wb") as raw,
    ):
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as outgoing:
            shutil.copyfileobj(incoming, outgoing)
    for key in ("viewer", "full_viewer", "common"):
        shutil.copy2(paths[key], output / paths[key].name)
    (output / "qm_indices.tcl").write_text(
        "set ::nyl12_full_qm_indices {"
        + " ".join(map(str, indices))
        + "}\n"
        + "set ::nyl12_qm_indices {"
        + " ".join(map(str, stripped_indices))
        + "}\n"
        + f"set ::nyl12_stripped_atom_count {len(kept)}\n"
    )
    (output / "frame_manifest.json").write_text(
        json.dumps(
            {
                "frames": frames,
                "initial_final_source_display_frames": [0, 47],
                "optimizer_states_not_MD_times": True,
                "duplicates_preserved": True,
                "physical_QM_atoms": 98,
                "stripped_atom_count": len(kept),
                "retained_water_topology_residues_1based": water_ids,
                "water_criterion": "Fixed union: oxygen within4A Euclidean distance of any QM heavy atom in any48 source frames; no reimaging",
                "stripped_to_original_atom_indices_0based": kept.tolist(),
                "virtual_link_caps_inserted": 0,
                "validation": validations,
            },
            indent=2,
        )
        + "\n"
    )
    (output / "README.md").write_text(
        "# Nyl12 optimization VMD inspection\n\nIn VMD: source /absolute/path/nyl12_es_optimization.vmd\n"
        "Default48frames keep all protein and full LIG with all H, plus fixed whole waters whose oxygen comes within4A of any QM heavy atom in any source frame (Euclidean, no reimaging). Bulk water and ions are stripped. Retained original water IDs and exact atom remapping are in frame_manifest.json.\n"
        "Commands: initial, final, pocket, system. All98 physical QM atoms use Licorice. Eight virtual link caps are not inserted.\n"
        "Frame0=ORIGINAL initial; then25first-pilot saved states,21continuation saved states and exact107cycle endpoint. Duplicates remain. Optimizer frames are not MD times or a reaction path; last saved state can differ from exact endpoint.\n"
        "Optional full-solvent comparison: source /absolute/path/nyl12_es_initial_final.vmd (2frames,103251atoms, all waters/ions/H). Or use full after loading the default viewer.\n"
        "gzip must be on PATH; topology expands beside viewer. DCD float32 display coordinates and all periodic boxes are checked against original sources. Original fullprecision coordinates remain on workstation. Native VMD GUI unavailable here.\n"
    )
    shutil.copy2(paths["handoff"], output / "standalone_handoff_20261020.md")
    if any(sha(path) != original_hashes[key] for key, path in paths.items()):
        raise ValueError("Input source changed during export")
    (output / "source_manifest.json").write_text(
        json.dumps(
            {
                "sources": {
                    key: {
                        "path": str(path),
                        "sha256": original_hashes[key],
                        "bytes": path.stat().st_size,
                    }
                    for key, path in paths.items()
                },
                "exporter_sha256": sha(Path(__file__)),
                "legacy_viewer_sha256": original_hashes["legacy_viewer"],
            },
            indent=2,
        )
        + "\n"
    )
    payload = [
        "md_solvated.parm7.gz",
        "initial_final.dcd",
        "optimization.dcd",
        paths["viewer"].name,
        paths["full_viewer"].name,
        paths["common"].name,
        "stripped.parm7.gz",
        "qm_indices.tcl",
        "frame_manifest.json",
        "source_manifest.json",
        "README.md",
        "standalone_handoff_20261020.md",
    ]
    (output / "MANIFEST.sha256").write_text(
        "".join(f"{sha(output / name)}  {name}\n" for name in payload)
    )
    logging.info(
        "Verified stripped48-frame optimizer movie and full-solvent2-frame comparison"
    )


def main() -> None:
    """Parse explicit source specification and export destination."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.spec, args.output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
