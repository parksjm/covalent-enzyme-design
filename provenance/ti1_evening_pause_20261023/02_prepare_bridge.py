"""Construct an explicit NH-water-OG1 bridge, then stage an SMD diagnostic."""

from __future__ import annotations
import hashlib
import importlib.util
import json
import logging
import shutil
from pathlib import Path
import numpy as np
from rdkit import Chem
from nyl12_covalent.active_site_geometry import angle_degrees
from nyl12_covalent.gaussian_checks import validate_geometry

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = Path(__file__).resolve().parent
OLD = ROOT / "output_nyl12_oxyanion_cluster_20261023"


def main() -> None:
    """Search deterministic bridge orientations without moving the TI1 core."""
    seed = OLD / "seed_01"
    mol = Chem.SDMolSupplier(str(seed / "native_seed.sdf"), removeHs=False)[0]
    mapping = json.loads((seed / "atom_map.json").read_text())
    ids = mapping["atom_index_by_origin"]
    xyz = np.array(
        [
            [float(v) for v in s.split()[1:]]
            for s in (seed / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    n, og = xyz[ids["THR:N"]], xyz[ids["THR:OG1"]]
    axis = (og - n) / np.linalg.norm(og - n)
    tangent = np.cross(axis, [1.0, 0.0, 0.0])
    tangent /= np.linalg.norm(tangent)
    other = np.cross(axis, tangent)
    water = [ids[s] for s in ("WAT:O", "WAT:H1", "WAT:H2")]
    nh = [ids[s] for s in ("THR:H1", "THR:H2", "THR:HG1")]
    best = None
    tested = 0
    for radius in (2.7, 2.8, 2.9, 3.0):
        height = np.sqrt(radius**2 - (np.linalg.norm(og - n) / 2) ** 2)
        for theta in np.linspace(0, 2 * np.pi, 180, endpoint=False):
            oxygen = (n + og) / 2 + height * (
                np.cos(theta) * tangent + np.sin(theta) * other
            )
            donor = max(nh, key=lambda h: angle_degrees(n, xyz[h], oxygen))
            angle = angle_degrees(n, xyz[donor], oxygen)
            if angle < 150:
                continue
            ohaxis = (og - oxygen) / np.linalg.norm(og - oxygen)
            cross = np.cross(ohaxis, axis)
            cross /= np.linalg.norm(cross)
            cross2 = np.cross(ohaxis, cross)
            for phi in np.linspace(0, 2 * np.pi, 36, endpoint=False):
                points = xyz.copy()
                points[water[0]] = oxygen
                points[water[1]] = oxygen + 0.96 * ohaxis
                direction = np.cos(np.deg2rad(104.5)) * ohaxis + np.sin(
                    np.deg2rad(104.5)
                ) * (np.cos(phi) * cross + np.sin(phi) * cross2)
                points[water[2]] = oxygen + 0.96 * direction
                tested += 1
                graph = validate_geometry(
                    mol, [a.GetAtomicNum() for a in mol.GetAtoms()], points
                )
                if not graph["accepted"]:
                    continue
                # Prefer linear NH donation, then fewer close nonbonded water contacts.
                remaining = [
                    i
                    for i in range(63)
                    if i not in water and i not in (ids["THR:N"], ids["THR:OG1"], donor)
                ]
                separation = np.linalg.norm(
                    points[water][:, None, :] - points[remaining][None, :, :], axis=2
                )
                score = (
                    angle - float(np.sum(np.maximum(0, 2.0 - separation) ** 2)) * 100
                )
                if best is None or score > best[0]:
                    best = (score, points, donor, angle, radius)
    if best is None:
        raise ValueError(
            "No clash-screened bridge satisfies NH angle>=150deg; noGaussian staged"
        )
    _, points, donor, angle, radius = best
    out = BASE / "seed_01"
    shutil.copytree(seed, out, ignore=shutil.ignore_patterns("source"))
    conf = mol.GetConformer()
    for i, p in enumerate(points):
        conf.SetAtomPosition(i, tuple(float(v) for v in p))
    with Chem.SDWriter(str(out / "native_seed.sdf")) as writer:
        writer.write(mol)
    (out / "native_seed.xyz").write_text(
        "63\nScripted NH-water-OG1 bridge; unoptimized SMD diagnostic\n"
        + "\n".join(
            a.GetSymbol() + " " + " ".join(f"{v:.10f}" for v in p)
            for a, p in zip(mol.GetAtoms(), points, strict=True)
        )
        + "\n"
    )
    report = json.loads((out / "preflight.json").read_text())
    report["water_bridge"] = {
        "NH_donor_origin": mapping["atoms"][donor]["origin_id"],
        "N_H_Owater_angle_deg": angle,
        "N_Owater_A": radius,
        "OG1_Owater_A": radius,
        "Owater_H_OG1_angle_deg": 180.0,
        "OH_A": 0.96,
        "HOH_deg": 104.5,
        "tested_orientations": tested,
        "nonwater_max_displacement_A": float(
            np.max(
                np.abs(
                    points[[i for i in range(63) if i not in water]]
                    - xyz[[i for i in range(63) if i not in water]]
                )
            )
        ),
    }
    report["coordinate_provenance_note"] = (
        "Nonwater unchanged from original seed; water coordinates replaced by deterministic bridge search. Original coordinate_provenance.json retained as historical source mapping only."
    )
    (out / "preflight.json").write_text(json.dumps(report, indent=2) + "\n")
    (out / "water_placement.json").write_text(
        json.dumps(report["water_bridge"], indent=2) + "\n"
    )
    for name in (
        "assess_pilot.py",
        "prepare_oxyanion_cluster.py",
        "donor_fragments.py",
    ):
        shutil.copy2(OLD / name, BASE / name)
    text = (
        (OLD / "prepare_svp_pilot.py")
        .read_text()
        .replace("oxyanion63_svp_20261023_pilot01", "bridge63_smd_20261023_pilot01")
        .replace(
            "#p wB97XD/def2SVP Opt=", "#p wB97XD/def2SVP SCRF=(SMD,Solvent=Water) Opt="
        )
        .replace(
            "wB97XD/def2SVP; no extra dispersion; gas phase",
            "wB97XD/def2SVP; SMD water; no extra dispersion",
        )
        .replace(
            "remaining protein electrostatics absent",
            "SMD water continuum; remaining specific protein electrostatics absent",
        )
    )
    (BASE / "prepare_smd_pilot.py").write_text(text)
    spec = importlib.util.spec_from_file_location(
        "smd_prep", BASE / "prepare_smd_pilot.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.main()
    snapshot = BASE / "svp_pilot_01/source/bridge_preparation.py"
    shutil.copy2(Path(__file__), snapshot)
    (BASE / "svp_pilot_01/source/bridge_seed_preflight.json").write_bytes(
        (out / "preflight.json").read_bytes()
    )
    (BASE / "bridge_source_hashes.json").write_text(
        json.dumps(
            {
                str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (
                    Path(__file__),
                    seed / "native_seed.xyz",
                    seed / "atom_map.json",
                    seed / "native_seed.sdf",
                )
            },
            indent=2,
        )
        + "\n"
    )
    logging.info("Bridge accepted; %s", report["water_bridge"])


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
