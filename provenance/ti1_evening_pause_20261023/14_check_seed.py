"""Measure user-edited water geometry and report original full-pocket contacts."""

from __future__ import annotations
import hashlib
import json
import logging
from pathlib import Path
import numpy as np
import parmed as pmd
from rdkit import Chem
from nyl12_covalent.active_site_geometry import angle_degrees
from nyl12_covalent.gaussian_checks import validate_geometry

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = Path(__file__).resolve().parent


def main() -> None:
    """Preserve user coordinates and produce measurements without optimization."""
    source = ROOT / "output_nyl12_oxyanion_cluster_20261023/seed_01"
    path = BASE / "native_seed_jmp.xyz"
    lines = path.read_text().splitlines()
    mol = Chem.SDMolSupplier(str(source / "native_seed.sdf"), removeHs=False)[0]
    xyz = np.array([[float(v) for v in line.split()[1:]] for line in lines[2:]])
    if (
        int(lines[0]) != 63
        or xyz.shape != (63, 3)
        or [a.GetSymbol() for a in mol.GetAtoms()]
        != [line.split()[0] for line in lines[2:]]
    ):
        raise ValueError("XYZ inventory/order mismatch")
    ids = json.loads((source / "atom_map.json").read_text())["atom_index_by_origin"]

    def distance(a: str, b: str) -> float:
        return float(np.linalg.norm(xyz[ids[a]] - xyz[ids[b]]))

    n, og, ow = (xyz[ids[s]] for s in ("THR:N", "THR:OG1", "WAT:O"))
    angles = {
        h: angle_degrees(n, xyz[ids[h]], ow) for h in ("THR:H1", "THR:H2", "THR:HG1")
    }
    water_angles = {h: angle_degrees(ow, xyz[ids[h]], og) for h in ("WAT:H1", "WAT:H2")}
    context = (
        ROOT
        / "output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/calculation"
    )
    model = pmd.load_file(
        str(context / "model.prmtop"), xyz=str(context / "prepared.rst7")
    )
    contacts = []
    for atom in model.atoms:
        if atom.residue.idx == 329:
            continue
        contacts.append(
            {
                "residue_topology_1based": atom.residue.idx + 1,
                "resname": atom.residue.name,
                "atom": atom.name,
                "element_Z": atom.atomic_number,
                "water_O_distance_A": float(
                    np.linalg.norm(ow - model.coordinates[atom.idx])
                ),
            }
        )
    report = {
        "user_xyz_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "geometry": validate_geometry(
            mol, [a.GetAtomicNum() for a in mol.GetAtoms()], xyz
        ),
        "water_O_N_A": distance("WAT:O", "THR:N"),
        "water_O_OG1_A": distance("WAT:O", "THR:OG1"),
        "N_H_Owater_angles_deg": angles,
        "Owater_H_OG1_angles_deg": water_angles,
        "water_OH_lengths_A": [distance("WAT:O", h) for h in ("WAT:H1", "WAT:H2")],
        "HOH_deg": angle_degrees(xyz[ids["WAT:H1"]], ow, xyz[ids["WAT:H2"]]),
        "nearest_full_context_atoms": sorted(
            contacts, key=lambda row: row["water_O_distance_A"]
        )[:15],
        "nearest_full_substrate_heavy_atoms": sorted(
            [
                r
                for r in contacts
                if r["residue_topology_1based"] == 329 and r["element_Z"] > 1
            ],
            key=lambda row: row["water_O_distance_A"],
        )[:5],
        "full_context_note": "Original unrelaxed alpha/beta pocket; distances reported, not a full vdW/energy acceptance. User seed unchanged.",
    }
    (BASE / "seed_checks.json").write_text(json.dumps(report, indent=2) + "\n")
    conf = mol.GetConformer()
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, tuple(float(v) for v in p))
    with Chem.SDWriter(str(BASE / "native_seed_jmp.sdf")) as writer:
        writer.write(mol)
    logging.info(
        "Water N/OG1 distances %.3f/%.3f; NH angles %s; OH angles %s",
        report["water_O_N_A"],
        report["water_O_OG1_A"],
        angles,
        water_angles,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
