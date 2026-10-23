"""Measure proton-transfer coordinates and candidate oxyanion donors; no optimization."""

from __future__ import annotations

import hashlib
import json
import logging
import sys
from pathlib import Path

import numpy as np
import parmed as pmd

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts.gaussian_log_to_xyz import parse_log  # noqa: E402
from nyl12_covalent.active_site_geometry import angle_degrees  # noqa: E402
from nyl12_covalent.local_repair_model import identity  # noqa: E402


def main() -> None:
    """Record distances in source/final snapshots and donor contacts in the source protein."""
    out = Path(__file__).resolve().parent
    pilot = ROOT / "output_nyl12_waterfree_reference_20261023/waterfree_pilot_01"
    data = (pilot / "optimization.log").read_bytes()
    (out / "stopped_optimization.log.txt").write_bytes(data)
    frames = parse_log(data.decode(), orientation="input", allow_incomplete=True).frames
    mapping = json.loads((pilot / "atom_map.json").read_text())["atom_index_by_origin"]

    def state(frame_index: int) -> dict:
        xyz = np.asarray(frames[frame_index].coordinates)

        def dist(a: str, b: str) -> float:
            return float(np.linalg.norm(xyz[mapping[a]] - xyz[mapping[b]]))

        return {
            "source_log_line": frames[frame_index].line,
            "OG1_C54_A": dist("THR:OG1", "PA6:C54"),
            "C54_O49_A": dist("PA6:C54", "PA6:O49"),
            "ammonium_H": {
                h: {"N_H_A": dist("THR:N", h), "OG1_H_A": dist("THR:OG1", h)}
                for h in ("THR:H1", "THR:H2", "THR:HG1")
            },
        }

    source = (
        ROOT
        / "output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/calculation"
    )
    model = pmd.load_file(
        str(source / "model.prmtop"), xyz=str(source / "prepared.rst7")
    )
    coords = np.asarray(model.coordinates)
    oxygen = next(a for a in model.residues[328].atoms if a.name == "O49")
    donors = []
    for atom in model.atoms:
        if atom.residue.idx >= 328 or atom.atomic_number not in (7, 8):
            continue
        hydrogens = [a for a in atom.bond_partners if a.atomic_number == 1]
        for hydrogen in hydrogens:
            distance = float(np.linalg.norm(coords[hydrogen.idx] - coords[oxygen.idx]))
            chain, resid = identity(atom.residue.idx)
            donors.append(
                {
                    "chain": chain,
                    "resid": resid,
                    "resname": atom.residue.name,
                    "topology_residue_1based": atom.residue.idx + 1,
                    "donor": atom.name,
                    "hydrogen": hydrogen.name,
                    "H_O49_A": distance,
                    "donor_O49_A": float(
                        np.linalg.norm(coords[atom.idx] - coords[oxygen.idx])
                    ),
                    "D_H_O49_angle_deg": angle_degrees(
                        coords[atom.idx], coords[hydrogen.idx], coords[oxygen.idx]
                    ),
                }
            )
    result = {
        "stopped_log_sha256": hashlib.sha256(data).hexdigest(),
        "frames": len(frames),
        "initial": state(0),
        "last_complete": state(-1),
        "candidate_donors_source_AF3_not_validated": sorted(
            donors, key=lambda x: x["H_O49_A"]
        )[:8],
        "interpretation": "Geometry evidence only; no replacement calculation or parameter acceptance",
    }
    (out / "report.json").write_text(json.dumps(result, indent=2) + "\n")
    logging.info("Wrote stopped-state/protein-donor coordinate review")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
