"""Map the unrelaxed native TI1 seed and physically preconditioned relay water."""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
from pathlib import Path

import numpy as np
import parmed as pmd
from rdkit import Chem
from rdkit.Chem.rdMolTransforms import GetDihedralDeg

from nyl12_covalent.gaussian_checks import validate_geometry
from nyl12_covalent.ti1_reference_graph import build_ti1_reference_graph

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = ROOT / "output_nyl12_native_reference_20261023"
SOURCE = (
    ROOT
    / "output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/calculation"
)
PA6 = ROOT / "provenance/xa_ma_capped_reference_20261022/pa6_l4_source.mol2"


def main() -> None:
    """Construct declared cap H, preserve all other named coordinates, and screen."""
    model = pmd.load_file(
        str(SOURCE / "model.prmtop"), xyz=str(SOURCE / "prepared.rst7")
    )
    wet = pmd.load_file(
        str(SOURCE / "model.prmtop"),
        xyz=str(SOURCE / "00_tip3p_water_only/original_order.rst7"),
    )
    if len(model.atoms) != 5054 or len(model.residues) != 330:
        raise ValueError("Unexpected source inventory")
    rows = {
        name: model.residues[idx]
        for name, idx in [("THR", 232), ("NEXT", 233), ("PA6", 328), ("WAT", 329)]
    }
    names = {
        key: {atom.name: atom.idx for atom in residue.atoms}
        for key, residue in rows.items()
    }
    if any(len(names[key]) != len(residue.atoms) for key, residue in rows.items()):
        raise ValueError("Source atom names are not unique")
    xyz_source = np.asarray(model.coordinates)
    xyz_water = np.asarray(wet.coordinates)
    water_indices = [atom.idx for atom in rows["WAT"].atoms]
    nonwater = [i for i in range(5054) if i not in water_indices]
    nonwater_shift = float(
        np.max(np.linalg.norm(xyz_water[nonwater] - xyz_source[nonwater], axis=1))
    )
    if nonwater_shift > 1e-6:
        raise ValueError("Water-only preconditioning changed nonwater coordinates")
    molecule, metadata = build_ti1_reference_graph(PA6, "R")
    points, provenance = [], []
    for atom in molecule.GetAtoms():
        origin = atom.GetProp("origin_id")
        group, name = origin.split(":")
        constructed = False
        if group == "NME":
            mapping = {
                "N": "N",
                "CH3": "CA",
                "H": "H",
                "HH31": "HA",
                "HH32": "C",
                "HH33": "CB",
            }
            source_name = mapping[name]
            source_index = names["NEXT"][source_name]
            point = xyz_source[source_index].copy()
            if name in {"HH32", "HH33"}:
                center = xyz_source[names["NEXT"]["CA"]]
                direction = point - center
                point = center + 1.09 * direction / np.linalg.norm(direction)
                constructed = True
        else:
            source_name = "H3" if origin == "THR:HG1" else name
            source_index = names[group][source_name]
            point = (xyz_water if group == "WAT" else xyz_source)[source_index].copy()
        points.append(point)
        provenance.append(
            {
                "index": atom.GetIdx(),
                "origin_id": origin,
                "source_atom_index_1based": source_index + 1,
                "source_atom_name": source_name,
                "constructed_cap_H": constructed,
                "source_coordinates": "water_only/original_order.rst7"
                if group == "WAT"
                else "prepared.rst7",
            }
        )
    xyz = np.asarray(points)
    conformer = Chem.Conformer(103)
    conformer.Set3D(True)
    for idx, point in enumerate(xyz):
        conformer.SetAtomPosition(idx, tuple(float(v) for v in point))
    molecule.AddConformer(conformer)
    Chem.AssignAtomChiralTagsFromStructure(molecule, replaceExistingTags=True)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    indices = metadata["atom_index_by_origin"]
    carbon = molecule.GetAtomWithIdx(indices["PA6:C54"])
    if not carbon.HasProp("_CIPCode"):
        raise ValueError("Source attack configuration cannot be assigned")
    cip = carbon.GetProp("_CIPCode")
    expected, metadata = build_ti1_reference_graph(PA6, cip)
    expected.AddConformer(conformer)
    checks = validate_geometry(
        expected, [a.GetAtomicNum() for a in expected.GetAtoms()], xyz
    )
    summary = {
        "accepted": checks["accepted"],
        "geometry": checks,
        "attacked_carbon_cip_from_geometry": cip,
        "formula": metadata["formula"],
        "formal_charge": 0,
        "explicit_H": sum(a.GetAtomicNum() == 1 for a in expected.GetAtoms()),
        "max_nonwater_preconditioning_displacement_A": nonwater_shift,
        "water_oxygen_physical_displacement_A": float(
            np.linalg.norm(xyz_water[names["WAT"]["O"]] - xyz_source[names["WAT"]["O"]])
        ),
        "water_preconditioning_converged": False,
        "scope": "Starting coordinates only; not optimized or parameter validated",
    }
    summary["reactive_distances_A"] = {
        label: float(np.linalg.norm(xyz[indices[left]] - xyz[indices[right]]))
        for label, left, right in [
            ("OG1_C54", "THR:OG1", "PA6:C54"),
            ("C54_O49", "PA6:C54", "PA6:O49"),
            ("C54_N59", "PA6:C54", "PA6:N59"),
            ("C54_C62", "PA6:C54", "PA6:C62"),
        ]
    }
    summary["psi_deg"] = GetDihedralDeg(
        conformer, *[indices[x] for x in ("THR:N", "THR:CA", "THR:C", "NME:N")]
    )
    summary["chi1_deg"] = GetDihedralDeg(
        conformer, *[indices[x] for x in ("THR:N", "THR:CA", "THR:CB", "THR:OG1")]
    )
    out = BASE / "native_seed_01"
    out.mkdir(exist_ok=False)
    with Chem.SDWriter(str(out / "native_seed.sdf")) as writer:
        writer.write(expected)
    lines = [
        "103",
        f"Unoptimized native TI1 starting seed; geometry-derived attack {cip}",
    ]
    lines.extend(
        f"{atom.GetSymbol():2s} " + " ".join(f"{v:.10f}" for v in point)
        for atom, point in zip(expected.GetAtoms(), xyz, strict=True)
    )
    (out / "native_seed.xyz").write_text("\n".join(lines) + "\n")
    for name, data in [
        ("atom_map.json", metadata),
        ("coordinate_provenance.json", provenance),
        ("preflight.json", summary),
    ]:
        (out / name).write_text(json.dumps(data, indent=2) + "\n")
    (out / "view_seed.pml").write_text(
        "load native_seed.sdf, native_ti1_seed\nset valence, 0\nshow sticks, native_ti1_seed\nshow spheres, native_ti1_seed and elem O\nset sphere_scale, 0.2\nutil.cnc native_ti1_seed\nzoom native_ti1_seed\n"
    )
    source_dir = out / "source"
    source_dir.mkdir()
    paths = [
        Path(__file__),
        PA6,
        SOURCE / "model.prmtop",
        SOURCE / "prepared.rst7",
        SOURCE / "00_tip3p_water_only/original_order.rst7",
        ROOT / "src/nyl12_covalent/ti1_reference_graph.py",
        ROOT / "src/nyl12_covalent/gaussian_checks.py",
    ]
    for path in paths:
        shutil.copy2(path, source_dir / path.name)
    hashes = {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths
    }
    (out / "source_manifest.json").write_text(json.dumps(hashes, indent=2) + "\n")
    logging.info(
        "Native seed: attack %s, accepted %s, %s", cip, checks["accepted"], summary
    )
    if not checks["accepted"]:
        raise SystemExit(3)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
