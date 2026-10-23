"""Export the stopped donor-cluster trajectory and refresh the laptop package."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import shutil
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts.gaussian_log_to_xyz import parse_log, trajectory_text  # noqa: E402
from scripts.prepare_nyl12_inspection_download import build_package  # noqa: E402


def main() -> None:
    """Retain all complete frames with explicit H and fixed input connectivity."""
    base = Path(__file__).resolve().parent
    pilot = base / "svp_pilot_01"
    out = base / "trajectory_snapshot_01"
    out.mkdir(exist_ok=False)
    raw = (pilot / "optimization.log").read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    (out / "optimization.log.txt").write_bytes(raw)
    parsed = parse_log(raw.decode(), orientation="input", allow_incomplete=True)
    mol = Chem.SDMolSupplier(str(pilot / "native_seed.sdf"), removeHs=False)[0]
    atoms = json.loads((pilot / "atom_map.json").read_text())["atoms"]
    origins = {row["origin_id"]: row["index"] for row in atoms}
    if mol is None or mol.GetNumAtoms() != 63:
        raise ValueError("Expected63atom graph")
    info = {
        "THR": ("THR", "E", 1),
        "NME": ("NME", "E", 2),
        "PA6": ("TIL", "I", 1),
        "WAT": ("WAT", "W", 1),
        "NMA": ("NMA", "B", 1),
        "ACM": ("ACM", "C", 1),
    }
    for atom, row in zip(mol.GetAtoms(), atoms, strict=True):
        origin = row["origin_id"]
        resname, chain, number = info[origin.split(":")[0]]
        name = origin.split(":")[-1]
        if origin == "THR:HG1":
            name = "H3"
        if len(name) > 4:
            name = f"H{atom.GetIdx() + 1}"
        atom.SetMonomerInfo(
            Chem.AtomPDBResidueInfo(
                atomName=f"{name:>4}",
                residueName=resname,
                residueNumber=number,
                chainId=chain,
                isHeteroAtom=True,
            )
        )
    pdb, metrics = (
        ["REMARK Fixed input connectivity; diagnostic failed TI1 pilot.\n"],
        [],
    )
    for count, frame in enumerate(parsed.frames, 1):
        if list(frame.atomic_numbers) != [a.GetAtomicNum() for a in mol.GetAtoms()]:
            raise ValueError("Element order mismatch")
        xyz = np.asarray(frame.coordinates)
        conf = Chem.Conformer(63)
        conf.Set3D(True)
        for i, point in enumerate(xyz):
            conf.SetAtomPosition(i, tuple(float(v) for v in point))
        mol.RemoveAllConformers()
        mol.AddConformer(conf)
        block = Chem.MolToPDBBlock(mol)
        pdb.append(f"MODEL     {count:4d}\n")
        pdb.extend(line + "\n" for line in block.splitlines() if line.strip() != "END")
        pdb.append("ENDMDL\n")

        def distance(a: str, b: str) -> float:
            return float(np.linalg.norm(xyz[origins[a]] - xyz[origins[b]]))

        metrics.append(
            {
                "frame": count,
                "log_line": frame.line,
                "OG1_C54_A": distance("THR:OG1", "PA6:C54"),
                "N_H3_A": distance("THR:N", "THR:HG1"),
                "OG1_H3_A": distance("THR:OG1", "THR:HG1"),
                "water_N_A": distance("WAT:O", "THR:N"),
                "water_OG1_A": distance("WAT:O", "THR:OG1"),
            }
        )
    pdb.append("END\n")
    (out / "optimization_trajectory.pdb").write_text("".join(pdb))
    (out / "optimization_trajectory.xyz").write_text(
        trajectory_text(parsed, "optimization.log", digest)
    )
    with (out / "trajectory_metrics.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
    (out / "view_optimization.pml").write_text(
        "load optimization_trajectory.pdb, ti1_progress, discrete=0\nhide everything, ti1_progress\nshow sticks, ti1_progress\nset valence, 0\nutil.cnc ti1_progress\nset movie_fps, 3\n"
        + f"mset 1 -{len(metrics)}\nzoom ti1_progress\n"
    )
    report = {
        "frames": len(metrics),
        "atoms_per_frame": 63,
        "explicit_H_per_frame": 37,
        "source_sha256": digest,
        "source_log": str(pilot / "optimization.log"),
        "last_frame": metrics[-1],
        "scientific_acceptance": False,
        "note": "Complete Input orientation tables including trial steps. Cycle-limited, unconverged and TI1 state-retention failure. PDB connectivity is fixed input graph; not evidence bonds survived.",
    }
    (out / "trajectory_report.json").write_text(json.dumps(report, indent=2) + "\n")
    shutil.copy2(
        pilot / "scientific_acceptance.json", out / "scientific_acceptance.json"
    )
    (out / "README.md").write_text(
        "# Failed TI1 pilot trajectory\n\nRun @view_optimization.pml then mplay from this directory in PyMOL. VMD can open optimization_trajectory.xyz. All37H including water H are explicit. These are complete geometry tables, not MD or necessarily accepted optimization steps. Fixed PDB bonds are display annotations. Pilot reached20cycle limit and failed TI1 state retention. See CSV and acceptance report. This is a snapshot, not an automatic live feed.\n"
    )
    reread = Chem.MolFromPDBBlock(
        "".join(pdb), removeHs=False, sanitize=False, proximityBonding=False
    )
    if (
        reread is None
        or reread.GetNumAtoms() != 63
        or reread.GetNumConformers() != len(metrics)
    ):
        raise ValueError("PDB roundtrip inventory failure")
    if sum(a.GetAtomicNum() == 1 for a in reread.GetAtoms()) != 37:
        raise ValueError("Missing explicit H in PDB")
    current = ROOT / "output_laptop_inspection_20261021"
    prior = json.loads((current / "provenance.json").read_text())
    files = {name: Path(row["source"]) for name, row in prior["files"].items()}
    for source in out.iterdir():
        if source.name != "optimization.log.txt":
            files["qm/ti1_oxyanion_cluster/optimization_progress/" + source.name] = (
                source
            )
    for name in ("status.json", "command.json"):
        snapshot = out / ("pilot_" + name)
        shutil.copy2(pilot / "job" / name, snapshot)
        files["qm/ti1_oxyanion_cluster/pilot/job/" + name] = snapshot
    staging = ROOT / "output_laptop_inspection_20261021_donor_movie_staging"
    backup = ROOT / "output_laptop_inspection_20261021_revision_23"
    if backup.exists() or staging.exists():
        raise FileExistsError("Fresh export paths required")
    build_package(files, staging)
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if hashlib.sha256((staging / name).read_bytes()).hexdigest() != expected:
            raise ValueError("Package hash failure")
    current.rename(backup)
    staging.rename(current)
    logging.info(
        "Verified%sframes,63atoms37H per frame; export refreshed. Last:%s",
        len(metrics),
        metrics[-1],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
