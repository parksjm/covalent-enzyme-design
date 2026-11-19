"""Run a bounded PROPKA screen on matched Tier A Nyl12 MD snapshots."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path
import subprocess

import MDAnalysis as mda

LOGGER = logging.getLogger(__name__)
CASES = ("neutral_replica_02", "ash272_replica_1", "ash274_replica_1")
CHAIN_LIMITS = (232, 464, 696, 928, 1024, 1120, 1216, 1312)
CHAIN_STARTS = (1, 233, 465, 697, 929, 1025, 1121, 1217)
CHAINS = "ABCDEFGH"


def sha256(path: Path) -> str:
    """Hash an input file in bounded memory.

    Args:
        path: Existing input file.

    Returns:
        Hexadecimal SHA-256 digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def chain_residue(resid: int) -> tuple[str, int]:
    """Map Amber residue IDs to eight mature Nyl12 chains.

    Args:
        resid: One-based Amber residue number.

    Returns:
        PDB chain letter and one-based chain-local residue number.
    """
    for chain, start, limit in zip(CHAINS, CHAIN_STARTS, CHAIN_LIMITS, strict=True):
        if start <= resid <= limit:
            return chain, resid - start + 1
    raise ValueError(f"Residue {resid} is outside the declared protein")


def representative(csv_path: Path) -> tuple[int, dict[str, str]]:
    """Choose the temporal median Tier A frame without energy cherry-picking.

    Args:
        csv_path: Matched per-frame geometry audit.

    Returns:
        One-based trajectory frame and its recorded geometry.
    """
    with csv_path.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["tier"] == "A"]
    if not rows:
        raise ValueError(f"No Tier A frames in {csv_path}")
    row = rows[len(rows) // 2]
    return int(row["frame"]), row


def write_protein_pdb(universe: mda.Universe, path: Path) -> None:
    """Write protein heavy atoms with explicit mature-chain boundaries.

    Args:
        universe: Universe positioned at the selected trajectory frame.
        path: New protein-only PDB file.
    """
    protein = universe.select_atoms("resid 1:1312")
    if len({atom.resid for atom in protein}) != 1312:
        raise ValueError("Expected exactly 1312 protein residues")
    previous_chain = ""
    serial = 0
    with path.open("w") as handle:
        for atom in protein:
            name = str(atom.name).strip()
            if name.startswith("H"):
                continue
            chain, local = chain_residue(int(atom.resid))
            if chain != previous_chain and previous_chain:
                handle.write("TER\n")
            previous_chain = chain
            serial += 1
            resname = "ASP" if atom.resname == "ASH" else str(atom.resname)
            if len(resname) != 3 or not name:
                raise ValueError(f"Unrecognized protein atom {atom}")
            x, y, z = atom.position
            handle.write(
                f"ATOM  {serial:5d} {name:>4s} {resname:3s} {chain}{local:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}"
                f"          {name[0]:>2s}\n"
            )
        handle.write("TER\nEND\n")


def main() -> None:
    """Select declared snapshots, write protein PDBs, and run PROPKA."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--propka", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    comparison = json.loads(args.comparison.read_text())
    manifest = json.loads(args.manifest.read_text())
    case_map = {case["name"]: case for case in manifest["cases"]}
    result_map = {item["name"]: item for item in comparison["replicas"]}
    records = []
    for name in CASES:
        case = case_map[name]
        frame, geometry = representative(Path(result_map[name]["frame_csv"]))
        source = [Path(case["topology"]), *map(Path, case["trajectories"])]
        universe = mda.Universe(str(source[0]), [str(p) for p in source[1:]])
        if len(universe.trajectory) != case["expected_frames"]:
            raise ValueError(f"Unexpected trajectory length for {name}")
        universe.trajectory[frame - 1]
        folder = args.output / name
        folder.mkdir()
        pdb = folder / "protein.pdb"
        write_protein_pdb(universe, pdb)
        command = [str(args.propka), "--pH", "8.0", pdb.name]
        completed = subprocess.run(
            command, cwd=folder, capture_output=True, text=True, check=False
        )
        (folder / "propka.log").write_text(completed.stdout + completed.stderr)
        if completed.returncode != 0:
            raise RuntimeError(f"PROPKA failed for {name}: {completed.returncode}")
        output = folder / "protein.pka"
        if not output.exists():
            raise FileNotFoundError(output)
        records.append(
            {
                "case": name,
                "state": case["state"],
                "frame": frame,
                "geometry": geometry,
                "source_sha256": {str(path): sha256(path) for path in source},
                "protein_pdb_sha256": sha256(pdb),
                "propka_output_sha256": sha256(output),
                "propka_command": command,
                "protein_only": True,
                "ash_heavy_residues_written_as_asp": True,
            }
        )
        LOGGER.info("%s: Tier A frame %d, PROPKA exit 0", name, frame)
    (args.output / "manifest.json").write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
