"""Report proximity and existing QM membership without modifying coordinates."""

from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from parmed.amber import AmberParm
from parmed.amber.mask import AmberMask

QM = "@1621-1655,2727-2734,13905-13931,14555-14560,14574-14579,19914-19917,19925-19931,19937-19941"
TARGETS = {
    "Thr_N": 13905,
    "Thr_OG1": 13916,
    "substrate_C": 19916,
    "substrate_O": 19917,
    "leaving_N": 19929,
}


def audit(topology: Path, coordinates: Path) -> dict:
    """Measure heavy/polar contacts and selected possible relay partners."""
    parm = AmberParm(str(topology), xyz=str(coordinates))
    xyz = parm.coordinates
    qm = np.asarray(AmberMask(parm, QM).Selection(), dtype=bool)
    result = {
        "sources": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (topology, coordinates)
        },
        "distance_convention": "Cartesian distances; no minimum-image wrapping",
        "qm_residues": [],
        "nearby_polar_residues": [],
        "selected_residues": [],
    }
    for residue in parm.residues:
        atoms = list(residue.atoms)
        ids = [a.idx for a in atoms]
        if qm[ids].any():
            result["qm_residues"].append(
                {
                    "amber_resid": residue.idx + 1,
                    "name": residue.name,
                    "qm_atoms": [a.name for a in atoms if qm[a.idx]],
                    "mm_atoms": [a.name for a in atoms if not qm[a.idx]],
                }
            )
        polar = [a for a in atoms if a.atomic_number in (7, 8, 16)]
        if not polar:
            continue
        best = min(
            (float(np.linalg.norm(xyz[a.idx] - xyz[t - 1])), a.name, label)
            for a in polar
            for label, t in TARGETS.items()
            if a.idx != t - 1
        )
        if best[0] < 4.0:
            result["nearby_polar_residues"].append(
                {
                    "amber_resid": residue.idx + 1,
                    "name": residue.name,
                    "distance_A": best[0],
                    "atom": best[1],
                    "target": best[2],
                    "has_qm_atoms": bool(qm[ids].any()),
                }
            )
        if residue.idx + 1 in (113, 156, 159, 186, 929, 930, 967, 968, 969, 970):
            details = []
            targets = {**TARGETS}
            for a in parm.residues[185].atoms:
                if a.name in ("OD1", "ND2"):
                    targets["Asn186_" + a.name] = a.idx + 1
            for a in polar:
                details.append(
                    {
                        "name": a.name,
                        "index": a.idx + 1,
                        "qm": bool(qm[a.idx]),
                        "distances_A": {
                            k: float(np.linalg.norm(xyz[a.idx] - xyz[i - 1]))
                            for k, i in targets.items()
                        },
                    }
                )
            result["selected_residues"].append(
                {
                    "amber_resid": residue.idx + 1,
                    "name": residue.name,
                    "polar_atoms": details,
                }
            )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = audit(args.topology, args.coordinates)
    data["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open("x") as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")
