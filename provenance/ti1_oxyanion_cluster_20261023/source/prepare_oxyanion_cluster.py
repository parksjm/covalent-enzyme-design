"""Build the approved unoptimized TI1 + N-methylacetamide + acetamide cluster."""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
from pathlib import Path

import numpy as np
import parmed as pmd
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from donor_fragments import _donor

from nyl12_covalent.active_site_geometry import angle_degrees
from nyl12_covalent.gaussian_checks import validate_geometry

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = Path(__file__).resolve().parent
SEED = ROOT / "output_nyl12_truncated_reference_20261023/truncated_pilot_01"
NATIVE = ROOT / "output_nyl12_native_reference_20261023/native_seed_01"
SOURCE = ROOT / (
    "output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/"
    "attempt_05_tip3p_preconditioned/calculation"
)


def _json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n")


def _read_seed() -> tuple[Chem.Mol, np.ndarray, list[dict]]:
    molecule = Chem.SDMolSupplier(str(SEED / "native_seed.sdf"), removeHs=False)[0]
    if molecule is None or molecule.GetNumAtoms() != 42:
        raise ValueError("Expected intact42atom explicit-H seed")
    rows = json.loads((SEED / "atom_map.json").read_text())["atoms"]
    lines = (SEED / "native_seed.xyz").read_text().splitlines()
    xyz = np.asarray([[float(v) for v in s.split()[1:]] for s in lines[2:]])
    if xyz.shape != (42, 3) or not np.isfinite(xyz).all():
        raise ValueError("Invalid42atom seed XYZ")
    if [s.split()[0] for s in lines[2:]] != [r["element"] for r in rows]:
        raise ValueError("Seed XYZ element order differs from atom map")
    for atom, row in zip(molecule.GetAtoms(), rows, strict=True):
        atom.SetProp("origin_id", row["origin_id"])
    return molecule, xyz, rows


def _frame_check(model: pmd.Structure, xyz: np.ndarray) -> tuple[list[dict], dict]:
    original = json.loads((NATIVE / "coordinate_provenance.json").read_text())
    truncated = json.loads((SEED / "coordinate_provenance.json").read_text())
    wet_path = SOURCE / "00_tip3p_water_only/original_order.rst7"
    wet = pmd.load_file(str(SOURCE / "model.prmtop"), xyz=str(wet_path))
    records, shifts, wet_shifts = [], [], []
    for index, row in enumerate(truncated):
        record = dict(row, index=index, seed42_index_0based=index)
        if not row["constructed_cap_H"]:
            native_row = original[row["native_source_index_0based"]]
            record.update(native_row, index=index, seed42_index_0based=index)
            if native_row["constructed_cap_H"]:
                records.append(record)
                continue
            src_idx = native_row["source_atom_index_1based"] - 1
            water = row["origin_id"].startswith("WAT:")
            reference = wet.coordinates if water else model.coordinates
            shift = float(np.linalg.norm(xyz[index] - reference[src_idx]))
            (wet_shifts if water else shifts).append(shift)
            if shift > 1e-8:
                raise ValueError(
                    f"Source coordinate-frame mismatch: {row['origin_id']}"
                )
        records.append(record)
    return records, {
        "common_nonwater_atoms": len(shifts),
        "common_nonwater_max_displacement_A": max(shifts),
        "water_atoms": len(wet_shifts),
        "water_source_max_displacement_A": max(wet_shifts),
        "alignment_applied": False,
    }


def build_cluster() -> tuple[Chem.Mol, np.ndarray, list[dict], dict]:
    """Construct approved donor graphs and verify exact core/source frame agreement.

    Returns:
        Molecular graph, coordinates, atom provenance, and frame verification.
    """
    core, xyz, _ = _read_seed()
    model = pmd.load_file(
        str(SOURCE / "model.prmtop"), xyz=str(SOURCE / "prepared.rst7")
    )
    if (len(model.atoms), len(model.residues)) != (5054, 330):
        raise ValueError("Unexpected Amber source inventory")
    if [model.residues[i - 1].name for i in (112, 113, 186)] != ["ILE", "PHE", "ASN"]:
        raise ValueError("Source donor residue identities differ")
    records, frame = _frame_check(model, xyz)
    nma = [(112, x) for x in ("CA", "HA", "C", "O")]
    nma += [(113, x) for x in ("N", "H", "CA", "HA")]
    acm = [(186, x) for x in ("CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21", "HD22")]
    specifications = [
        (
            "NMA",
            nma,
            [
                (112, "CA", "N", "HCAP_N"),
                (112, "CA", "CB", "HCAP_CB"),
                (113, "CA", "C", "HCAP_C"),
                (113, "CA", "CB", "HCAP_CB"),
            ],
        ),
        ("ACM", acm, [(186, "CB", "CA", "HCAP_CA")]),
    ]
    points = list(xyz)
    for prefix, selected, cuts in specifications:
        donor, donor_points, donor_records = _donor(model, prefix, selected, cuts)
        core = Chem.CombineMols(core, donor)
        points.extend(donor_points)
        records.extend(donor_records)
    core.RemoveAllConformers()
    conformer = Chem.Conformer(core.GetNumAtoms())
    conformer.Set3D(True)
    for index, point in enumerate(points):
        conformer.SetAtomPosition(index, tuple(float(v) for v in point))
        records[index]["index"] = index
    core.AddConformer(conformer)
    Chem.SanitizeMol(core)
    return core, np.asarray(points), records, frame


def assess_cluster(molecule: Chem.Mol, xyz: np.ndarray) -> dict:
    """Screen graph inventory, coordinate-derived stereo and donor contacts.

    Args:
        molecule: Constructed explicit-H molecular graph.
        xyz: Coordinates in the original protein frame.

    Returns:
        Separate chemical inventory and conservative geometry screening results.
    """
    indices = {a.GetProp("origin_id"): a.GetIdx() for a in molecule.GetAtoms()}
    formula, charge = CalcMolFormula(molecule), Chem.GetFormalCharge(molecule)
    fragments = sorted(len(x) for x in Chem.GetMolFrags(molecule))
    if (molecule.GetNumAtoms(), formula, charge, fragments) != (
        63,
        "C15H37N5O6",
        0,
        [3, 9, 12, 39],
    ):
        raise ValueError(
            "63atom cluster graph/formula/charge/component inventory mismatch"
        )
    expected = {"THR:CA": "S", "THR:CB": "R", "PA6:C54": "R"}
    derived = Chem.Mol(molecule)
    Chem.AssignAtomChiralTagsFromStructure(derived, replaceExistingTags=True)
    Chem.AssignStereochemistry(derived, cleanIt=True, force=True)
    stereo = {
        name: derived.GetAtomWithIdx(indices[name]).GetProp("_CIPCode")
        for name in expected
    }
    if stereo != expected:
        raise ValueError(f"Coordinate-derived stereochemistry mismatch: {stereo}")
    geometry = validate_geometry(
        molecule, [a.GetAtomicNum() for a in molecule.GetAtoms()], xyz
    )
    contacts = []
    oxygen = xyz[indices["PA6:O49"]]
    for nitrogen, hydrogens in [
        ("NMA:113:N", ["NMA:113:H"]),
        ("ACM:186:ND2", ["ACM:186:HD21", "ACM:186:HD22"]),
    ]:
        n = xyz[indices[nitrogen]]
        for hydrogen in hydrogens:
            h = xyz[indices[hydrogen]]
            contacts.append(
                {
                    "donor": nitrogen,
                    "hydrogen": hydrogen,
                    "N_O_A": float(np.linalg.norm(n - oxygen)),
                    "H_O_A": float(np.linalg.norm(h - oxygen)),
                    "N_H_O_deg": angle_degrees(n, h, oxygen),
                }
            )
    return {
        "chemical_inventory_passed": True,
        "atom_count": 63,
        "formula": formula,
        "formal_charge": charge,
        "multiplicity": 1,
        "explicit_H": 37,
        "component_atom_counts": fragments,
        "coordinate_stereocenters": stereo,
        "geometry": geometry,
        "geometry_screen_passed": geometry["accepted"],
        "donor_contacts": contacts,
        "optimized": False,
        "parameter_validated": False,
    }


def main() -> None:
    """Write an independent seed and immutable hashed inputs, with no QM job."""
    molecule, xyz, records, frame = build_cluster()
    report = assess_cluster(molecule, xyz)
    report["coordinate_frame"] = frame
    report["retained_seed42_max_displacement_A"] = float(
        np.max(np.abs(xyz[:42] - _read_seed()[1]))
    )
    out = BASE / "seed_01"
    out.mkdir(exist_ok=False)
    block = Chem.MolToMolBlock(molecule, forceV3000=True).splitlines()
    in_atoms = False
    for i, line in enumerate(block):
        if line == "M  V30 BEGIN ATOM":
            in_atoms = True
        elif line == "M  V30 END ATOM":
            in_atoms = False
        elif in_atoms:
            fields = line.split()
            fields[4:7] = [f"{v:.10f}" for v in xyz[int(fields[2]) - 1]]
            block[i] = "M  V30 " + " ".join(fields[2:])
    (out / "native_seed.sdf").write_text("\n".join(block) + "\n$$$$\n")
    lines = ["63", "Unoptimized TI1 + Phe113 NMA + Asn186 acetamide; neutral singlet"]
    lines += [
        f"{a.GetSymbol():2s} " + " ".join(f"{v:.10f}" for v in p)
        for a, p in zip(molecule.GetAtoms(), xyz, strict=True)
    ]
    (out / "native_seed.xyz").write_text("\n".join(lines) + "\n")
    atoms = [
        dict(
            index=a.GetIdx(),
            origin_id=a.GetProp("origin_id"),
            element=a.GetSymbol(),
            formal_charge=a.GetFormalCharge(),
        )
        for a in molecule.GetAtoms()
    ]
    _json(
        out / "atom_map.json",
        {
            "formula": report["formula"],
            "formal_charge": 0,
            "multiplicity": 1,
            "atom_index_by_origin": {r["origin_id"]: r["index"] for r in atoms},
            "atoms": atoms,
        },
    )
    _json(out / "coordinate_provenance.json", records)
    _json(out / "preflight.json", report)
    (out / "view_seed.pml").write_text(
        "load native_seed.sdf, oxyanion_ti1\nset valence, 0\nshow sticks, oxyanion_ti1\nshow sticks, oxyanion_ti1 and hydro\nutil.cnc oxyanion_ti1\nzoom oxyanion_ti1\n"
    )
    snapshot = out / "source"
    snapshot.mkdir()
    paths = [
        Path(__file__),
        BASE / "donor_fragments.py",
        SEED / "native_seed.sdf",
        SEED / "native_seed.xyz",
        SEED / "atom_map.json",
        SEED / "coordinate_provenance.json",
        NATIVE / "coordinate_provenance.json",
        SOURCE / "model.prmtop",
        SOURCE / "prepared.rst7",
        SOURCE / "00_tip3p_water_only/original_order.rst7",
        ROOT / "src/nyl12_covalent/gaussian_checks.py",
        ROOT / "src/nyl12_covalent/active_site_geometry.py",
    ]
    manifest = []
    for i, path in enumerate(paths):
        target = snapshot / f"{i:02d}_{path.name}"
        shutil.copy2(path, target)
        manifest.append(
            {
                "source": str(path),
                "snapshot": str(target.relative_to(out)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    _json(out / "source_manifest.json", manifest)
    logging.info(
        "63atom seed constructed; geometry screen passed: %s",
        report["geometry_screen_passed"],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
