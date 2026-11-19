"""Inspect retained PA6 ES seeds; export mapped pockets without optimization."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.active_site_geometry import measure_active_site  # noqa: E402
from nyl12_covalent.es_water_inventory import inventory_es_waters  # noqa: E402
from nyl12_covalent.es_seed_legacy_geometry import METRICS, _compare, inspect_legacy_phi  # noqa: E402, F401
from nyl12_covalent.es_seed_display import _align, _display, _viewer, _write_pdb  # noqa: E402, F401

LOGGER = logging.getLogger(__name__)


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atom(parm: Any, row: dict) -> Any:
    residue = parm.residues[row["residue"] - 1]
    atoms = [atom for atom in residue.atoms if atom.name == row["atom"]]
    if residue.name != row["resname"] or len(atoms) != 1:
        raise ValueError(f"Missing or ambiguous mapped atom: {row}")
    return atoms[0]


def _geometry(parm: Any, mapping: dict) -> dict:
    xyz = parm.coordinates.copy()
    keys = {
        "nucleophile": "Thr929.OG1",
        "oxyanion_donor_backbone": "Phe113.N",
        "oxyanion_donor_sidechain": "Asn186.ND2",
        "scissile_carbon": "substrate.C54",
        "scissile_oxygen": "substrate.O49",
        "leaving_nitrogen": "substrate.N59",
    }
    coordinates = {
        name: tuple(xyz[_atom(parm, mapping[role]).idx]) for role, name in keys.items()
    }
    coordinates["Thr929.CB"] = tuple(
        xyz[_atom(parm, {"residue": 929, "resname": "THR", "atom": "CB"}).idx]
    )
    coordinates["substrate.C62"] = tuple(
        xyz[_atom(parm, {"residue": 1313, "resname": "LIG", "atom": "C62"}).idx]
    )
    result = measure_active_site(coordinates, "pa6", "es")
    rounded = {
        key: tuple(float(f"{value:.3f}") for value in point)
        for key, point in coordinates.items()
    }
    result["NC_rounded_phi"] = measure_active_site(rounded, "pa6", "es")["metrics"][
        "herschlag_phi_projected"
    ]
    result["NC_rounding_policy"] = (
        "Diagnostic rounding of writtenNC to3 decimals; actual historicalPDB separately verified; rawNCphi remains primary"
    )
    return result


def _pocket(parm: Any, mapping: dict, waters: dict) -> set[int]:
    xyz = parm.coordinates.copy()
    centers = np.asarray(
        [
            xyz[_atom(parm, mapping[key]).idx]
            for key in [
                "thr_n_terminus",
                "nucleophile",
                "scissile_carbon",
                "scissile_oxygen",
            ]
        ]
    )
    residues = {1313} | {water["residue_index"] for water in waters["waters"]}
    for residue in parm.residues[:1312]:
        heavy = [atom.idx for atom in residue.atoms if atom.atomic_number > 1]
        if heavy and np.min(np.linalg.norm(xyz[heavy, None] - centers, axis=2)) <= 6:
            residues.add(residue.idx + 1)
    return residues


def inspect_seeds(reuse_audit: Path, manifest: Path, output: Path) -> dict:
    """Write verified seed geometry, local-water inventories and mapped PDBs.

    Args:
        reuse_audit: Completed source audit containing nine extracted frames.
        manifest: Explicit PA6 source and atom mappings.
        output: New directory; existing paths are rejected.

    Returns:
        Inspection report. Geometry mismatches are retained before raising.
    """
    audit, config = (
        json.loads(reuse_audit.read_text()),
        json.loads(manifest.read_text()),
    )
    frames = audit["selected_frames"]
    if len(frames) != 9 or not audit["acceptance"]["reuse_audit_complete"]:
        raise ValueError("Expected a completed reuse audit with nine selected frames")
    for source in [
        manifest,
        Path(config["geometry_csv"]),
        Path(config["system"]) / "md_solvated.prmtop",
    ]:
        matches = [
            row
            for row in audit["sources"]
            if Path(row["path"]).resolve() == source.resolve()
        ]
        if len(matches) != 1 or _hash(source) != matches[0]["sha256"]:
            raise ValueError(f"Audit source hash mismatch: {source}")
    output.mkdir(exist_ok=False)
    topology = Path(config["system"]) / "md_solvated.prmtop"
    mapping = config["atom_map"]
    water_roles = {key: mapping[key] for key in ["thr_n_terminus", "nucleophile"]}
    water_roles.update(
        {
            f"asp{number}_{atom.lower()}": {
                "residue": residue,
                "resname": "ASP",
                "atom": atom,
            }
            for number, residue in [(272, 968), (274, 970)]
            for atom in ["OD1", "OD2"]
        }
    )
    with Path(config["geometry_csv"]).open(newline="") as handle:
        csv_rows = list(csv.DictReader(handle))
    source_paths = [
        reuse_audit,
        manifest,
        topology,
        Path(config["geometry_csv"]),
        Path(__file__),
        ROOT / "src/nyl12_covalent/active_site_geometry.py",
        ROOT / "src/nyl12_covalent/es_water_inventory.py",
        ROOT / "src/nyl12_covalent/es_seed_display.py",
        ROOT / "src/nyl12_covalent/es_seed_legacy_geometry.py",
        Path(
            "/data/jpc/enzyme-substrate-modeling/output_nyl12_pa6_md_pilot_100ns_20260818/monitor_active_site_geometry.py"
        ),
    ]
    report = {
        "state": "ES seeds; not optimized TI/TS",
        "water_roles": water_roles,
        "frames": [],
    }
    reference_ca = None
    for index, frame in enumerate(frames):
        coordinates = Path(frame["solvated_coordinates"]["path"])
        if _hash(coordinates) != frame["solvated_coordinates"]["sha256"]:
            raise ValueError(f"Source frame hash mismatch: {coordinates}")
        source_paths.append(coordinates)
        rows = [
            row
            for row in csv_rows
            if row["model"] == frame["model"]
            and row["replica"] == frame["replica"]
            and int(row["frame"]) == frame["frame"]
        ]
        if len(rows) != 1:
            raise ValueError(
                "Expected one original geometry CSV row for selected frame"
            )
        parm = AmberParm(str(topology), xyz=str(coordinates))
        geometry = _geometry(parm, mapping)
        legacy_path = (
            Path(config["geometry_csv"]).parent
            / frame["model"]
            / frame["replica"]
            / f"herschlag_phi_atoms.pdb.{frame['frame']}"
        )
        legacy = inspect_legacy_phi(parm, legacy_path)
        geometry["historical_PDB"] = legacy
        source_paths.append(legacy_path)
        comparison = _compare(geometry["metrics"], rows[0], legacy["phi_degrees"])
        if not all(item["passed"] for item in comparison.values()):
            (output / "geometry_mismatch.json").write_text(
                json.dumps({"frame": frame, "comparison": comparison}, indent=2)
            )
            raise ValueError(
                f"Geometry mismatch: {frame['replica']} frame {frame['full_frame']}"
            )
        waters = inventory_es_waters(topology, coordinates, water_roles)
        reference_ca, alignment = _align(parm, reference_ca)
        name = f"{frame['replica']}_frame_{frame['full_frame']:04d}"
        pocket = _write_pdb(
            parm,
            _pocket(parm, mapping, waters),
            [water["residue_index"] for water in waters["waters"]],
            output / f"{name}.pdb",
        )
        if index == 0:
            context = _write_pdb(
                parm, set(range(1, 1314)), [], output / "fullprotein_substrate.pdb"
            )
            (output / "context_mapping.json").write_text(json.dumps(context, indent=2))
        (output / f"{name}_mapping.json").write_text(json.dumps(pocket, indent=2))
        report["frames"].append(
            {
                "name": name,
                "source": frame,
                "geometry": geometry,
                "comparison": comparison,
                "water_inventory": waters,
                "display_alignment": alignment,
            }
        )
        LOGGER.info("Verified %s: %d local waters", name, waters["local_water_count"])
    _viewer(output, [frame["name"] for frame in report["frames"]])
    report["sources"] = [
        {"path": str(path), "sha256": _hash(path)} for path in source_paths
    ]
    (output / "inspection.json").write_text(json.dumps(report, indent=2) + "\n")
    (output / "README.txt").write_text(
        "Nine retained classical ES seeds, not optimized TI/TS.\nIn PyMOL: run /absolute/path/view_es_seeds.py (portable from any working directory).\nAlternatively cd to this folder and @view_es_seeds.pml.\nEnable one pocket object at a time; context is first-frame reference only.\nProtein beta chains E-H retain precursor numbering233-328.\nWater contacts and explicit H angles are descriptive; no relay assigned.\nRawNCdistances/aBD match CSV within0.001 A/0.01degree; source-matched3-decimalPDBphi within1e-5degree.\n"
    )
    frozen = output / "source"
    frozen.mkdir()
    for path in source_paths[:2] + [
        path for path in source_paths if path.suffix == ".py"
    ]:
        shutil.copy2(path, frozen / path.name)
    for frame in frames:
        path = (
            Path(config["geometry_csv"]).parent
            / frame["model"]
            / frame["replica"]
            / f"herschlag_phi_atoms.pdb.{frame['frame']}"
        )
        shutil.copy2(path, frozen / f"{frame['replica']}_{path.name}")
    outputs = [
        {"path": str(path.relative_to(output)), "sha256": _hash(path)}
        for path in sorted(output.rglob("*"))
        if path.is_file()
    ]
    (output / "output_hashes.json").write_text(json.dumps(outputs, indent=2) + "\n")
    return report


def main() -> None:
    """Parse explicit input/output locations and perform read-only inspection."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-audit", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    inspect_seeds(args.reuse_audit, args.manifest, args.output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
