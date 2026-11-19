"""Prepare a bounded, atom-preserving direct acylenzyme endpoint hypothesis."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
import re
import shutil
import sys

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.acylenzyme_seed import build_acylenzyme_seed  # noqa: E402

LOGGER = logging.getLogger(__name__)
EXPECTED = "d57540b6bd5ede95c360174ea4d1bd92c80625355d3c2ccd13f321443b9674ad"


def prepare(project: Path) -> None:
    """Search fourteen explicit internal geometries and retain every report."""
    base = (
        project
        / "qmmm/pa6_flipped_model30/05_neutral_reaction_path/04_direct_acylenzyme"
    )
    source = base / "work/01_coupled_coordinate_drive/asynchronous_20/optimized.ncrst"
    topology = base / "inputs/system.prmtop"
    mdin = base / "inputs/fixed_shell_compact_residues.mdin"
    output = base / "acylenzyme_endpoint"
    if output.exists():
        raise FileExistsError(output)
    if hashlib.sha256(source.read_bytes()).hexdigest() != EXPECTED:
        raise ValueError("Window20 source hash mismatch")
    reference = Rst7.open(str(source))
    xyz = reference.coordinates.reshape(-1, 3)
    parm = AmberParm(str(topology), xyz=str(source))
    match = re.search(r"active='([^']+)'", mdin.read_text())
    if match is None:
        raise ValueError("Missing active mask")
    active = np.array(AmberMask(parm, match.group(1)).Selection(), dtype=bool)
    if active.sum() != 4165:
        raise ValueError("Original 4165 movable shell not preserved")
    output.mkdir()
    inputs = output / "inputs"
    inputs.mkdir()
    snapshots = output / "source_snapshot"
    snapshots.mkdir()
    files = [
        source,
        topology,
        mdin,
        Path(__file__),
        project / "src/nyl12_covalent/acylenzyme_seed.py",
    ]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    for p in files:
        shutil.copy2(p, (inputs if p in files[:3] else snapshots) / p.name)
    best_score = None
    best_seed = None
    best_report = None
    reports = []
    radii = {1: 1.2, 6: 1.7, 7: 1.55, 8: 1.52, 16: 1.8}
    for cn in (2.6, 2.85):
        for bend in (-30, -20, -10, 0, 10, 20, 30):
            seed, report = build_acylenzyme_seed(
                parm,
                xyz,
                ester_distance=1.4,
                leaving_distance=cn,
                orient_leaving=True,
                source_proton_state="leaving_n",
                leaving_bend_degrees=bend,
            )
            displacements = np.linalg.norm(seed - xyz, axis=1)
            if np.any(displacements[~active] > 1e-10):
                raise ValueError("Builder moved frozen atoms")
            ratios = []
            for contact in report["close_contacts"]:
                atoms = [parm.atoms[i - 1] for i in contact["atoms_1based"]]
                contact["identities"] = [
                    f"{a.residue.name}{a.residue.idx + 1}:{a.name}(Z={a.atomic_number})"
                    for a in atoms
                ]
                contact["heavy_pair"] = all(a.atomic_number > 1 for a in atoms)
                ratio = contact["distance_A"] / (
                    0.65 * sum(radii.get(a.atomic_number, 1.8) for a in atoms)
                )
                contact["cutoff_ratio"] = ratio
                ratios.append(ratio)
            report.update(
                {
                    "max_displacement_A": float(displacements.max()),
                    "max_displacement_atom_1based": int(displacements.argmax() + 1),
                    "frozen_max_displacement_A": float(displacements[~active].max()),
                    "leaving_N_original_HN_19937_A": float(
                        np.linalg.norm(seed[19928] - seed[19936])
                    ),
                    "OG1_HG1_A": float(np.linalg.norm(seed[13915] - seed[13916])),
                    "displacements_A": {
                        str(i + 1): float(displacements[i])
                        for i in np.flatnonzero(displacements > 1e-10)
                    },
                }
            )
            score = (min(ratios, default=1.0), -sum(1 - r for r in ratios))
            report["clash_score"] = score
            reports.append(report)
            if best_score is None or score > best_score:
                best_score, best_seed, best_report = score, seed.copy(), report
    if best_seed is None:
        raise ValueError("No seeds generated")
    restart = Rst7(natom=len(best_seed))
    restart.coordinates = best_seed
    restart.box = reference.box
    restart.write(str(inputs / "acylenzyme_seed.ncrst"), netcdf=True)
    parm.coordinates = best_seed
    parm.save(str(inputs / "acylenzyme_seed.pdb"), overwrite=False)
    protocol = {
        "status": "geometry review pending; no energy optimization launched",
        "source_sha256": hashes,
        "method_planned": "DFTB3/3OB periodic QM/MM; 98 QM + 8 links; charge -2",
        "movable_atoms": 4165,
        "outer_atoms_frozen": True,
        "reaction_restraints_planned": "none; no proton holding or restrained preoptimization",
        "atom_count": len(xyz),
        "topology_charge": float(sum(a.charge for a in parm.atoms)),
        "selection": "maximize worst contact/cutoff ratio then minimize summed overlap; all contacts treated equally",
        "construction": "Fixed Thr OG1/acyl tail; nearest two-sphere C placement; planar C=O1.23; detached rigid fragment translated to explicit C-N length and bend; HG1 open valence on leaving N; bounded orientation <=15 degrees",
        "classical_product_topology_valid": False,
        "best_seed": best_report,
        "all_candidate_reports": reports,
    }
    (output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    LOGGER.info(
        "Best C-N %.2f, bend %s, maxmove %.3f, contacts %d",
        best_report["distances_A"]["C_Nleaving"],
        best_report["leaving_bend_degrees"],
        best_report["max_displacement_A"],
        len(best_report["close_contacts"]),
    )
    LOGGER.info("Best report: %s", json.dumps(best_report))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    prepare(ROOT)
