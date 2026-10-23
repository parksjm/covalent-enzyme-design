"""Prepare a source-mapped native-state TI1 pilot with ethyl substrate boundaries."""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from nyl12_covalent.gaussian_checks import validate_geometry

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = Path(__file__).resolve().parent
SOURCE = ROOT / "output_nyl12_native_reference_20261023/native_seed_01"
OUT = BASE / "truncated_pilot_01"
LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Cap two alkyl cuts, screen chemistry, and stage an independent Gaussian job."""
    molecule = Chem.SDMolSupplier(str(SOURCE / "native_seed.sdf"), removeHs=False)[0]
    if molecule is None or molecule.GetNumAtoms() != 103:
        raise ValueError("Expected original explicit-H 103-atom seed")
    metadata = json.loads((SOURCE / "atom_map.json").read_text())
    original = metadata["atom_index_by_origin"]
    xyz = np.asarray(
        [
            [float(v) for v in line.split()[1:]]
            for line in (SOURCE / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    for row in metadata["atoms"]:
        molecule.GetAtomWithIdx(row["index"]).SetProp("origin_id", row["origin_id"])
    selected_heavy = {"C54", "O49", "N59", "C62", "C70", "C67", "C75"}
    keep = {
        a.GetIdx()
        for a in molecule.GetAtoms()
        if not a.GetProp("origin_id").startswith("PA6:")
        or a.GetProp("origin_id").split(":")[1] in selected_heavy
    }
    for idx in list(keep):
        keep.update(
            a.GetIdx()
            for a in molecule.GetAtomWithIdx(idx).GetNeighbors()
            if a.GetAtomicNum() == 1
        )
    retained = sorted(keep)
    rw = Chem.RWMol(molecule)
    for idx in reversed(range(103)):
        if idx not in keep:
            rw.RemoveAtom(idx)
    positions = [xyz[idx].copy() for idx in retained]
    provenance = [
        {
            "origin_id": molecule.GetAtomWithIdx(idx).GetProp("origin_id"),
            "native_source_index_0based": idx,
            "constructed_cap_H": False,
        }
        for idx in retained
    ]
    current = {a.GetProp("origin_id"): a.GetIdx() for a in rw.GetAtoms()}
    cuts = [("PA6:C70", "PA6:C78", "PA6:HC70"), ("PA6:C75", "PA6:C79", "PA6:HC75")]
    for center, removed, cap_origin in cuts:
        center_idx = current[center]
        atom = rw.GetAtomWithIdx(center_idx)
        atom.SetNoImplicit(True)
        if len(atom.GetNeighbors()) != 3:
            raise ValueError(f"Unexpected cap boundary valence: {center}")
        hydrogen = Chem.Atom(1)
        hydrogen.SetProp("origin_id", cap_origin)
        hydrogen_idx = rw.AddAtom(hydrogen)
        rw.AddBond(center_idx, hydrogen_idx, Chem.BondType.SINGLE)
        direction = xyz[original[removed]] - xyz[original[center]]
        positions.append(
            xyz[original[center]] + 1.09 * direction / np.linalg.norm(direction)
        )
        provenance.append(
            {
                "origin_id": cap_origin,
                "native_source_index_0based": None,
                "constructed_cap_H": True,
                "cut": [center, removed],
                "CH_distance_A": 1.09,
            }
        )
    product = rw.GetMol()
    Chem.SanitizeMol(product)
    product.RemoveAllConformers()
    conformer = Chem.Conformer(product.GetNumAtoms())
    conformer.Set3D(True)
    for idx, point in enumerate(positions):
        conformer.SetAtomPosition(idx, tuple(float(v) for v in point))
    product.AddConformer(conformer)
    Chem.AssignStereochemistry(product, cleanIt=True, force=True)
    indices = {a.GetProp("origin_id"): a.GetIdx() for a in product.GetAtoms()}
    expected_stereo = {
        indices[name]: code
        for name, code in [("THR:CA", "S"), ("THR:CB", "R"), ("PA6:C54", "R")]
    }
    actual_stereo = {
        idx: product.GetAtomWithIdx(idx).GetProp("_CIPCode") for idx in expected_stereo
    }
    formula = CalcMolFormula(product)
    hydrogens = sum(a.GetAtomicNum() == 1 for a in product.GetAtoms())
    if (product.GetNumAtoms(), formula, Chem.GetFormalCharge(product), hydrogens) != (
        42,
        "C10H25N3O4",
        0,
        25,
    ):
        raise ValueError("Truncated atom/formula/charge/H inventory mismatch")
    if expected_stereo != actual_stereo:
        raise ValueError("Graph stereochemistry changed during truncation")
    geometry = validate_geometry(
        product, [a.GetAtomicNum() for a in product.GetAtoms()], np.asarray(positions)
    )
    OUT.mkdir(exist_ok=False)
    with Chem.SDWriter(str(OUT / "native_seed.sdf")) as writer:
        writer.write(product)
    lines = [
        "42",
        "Unoptimized native-state R TI1; two ethyl boundaries; explicit water H",
    ]
    lines += [
        f"{atom.GetSymbol():2s} " + " ".join(f"{v:.10f}" for v in point)
        for atom, point in zip(product.GetAtoms(), positions, strict=True)
    ]
    (OUT / "native_seed.xyz").write_text("\n".join(lines) + "\n")
    metadata = {
        "formula": formula,
        "formal_charge": 0,
        "attacked_carbon_cip": "R",
        "atom_index_by_origin": indices,
        "atoms": [
            {
                "index": a.GetIdx(),
                "origin_id": a.GetProp("origin_id"),
                "element": a.GetSymbol(),
                "formal_charge": a.GetFormalCharge(),
            }
            for a in product.GetAtoms()
        ],
    }
    preflight = {
        "accepted": geometry["accepted"],
        "geometry": geometry,
        "atom_count": 42,
        "explicit_H": hydrogens,
        "formula": formula,
        "formal_charge": 0,
        "retained_original_atoms": 40,
        "new_cap_H": 2,
        "retained_coordinate_max_displacement_A": float(
            np.max(np.linalg.norm(np.asarray(positions[:40]) - xyz[retained], axis=1))
        ),
        "cuts": cuts,
        "source": str(SOURCE),
        "forcefield_validated": False,
    }
    mobile_heavy = {"THR:OG1", "PA6:C54", "PA6:O49", "PA6:N59", "WAT:O"}
    mobile = [
        a.GetAtomicNum() == 1 or a.GetProp("origin_id") in mobile_heavy
        for a in product.GetAtoms()
    ]
    frozen = [i for i, allowed in enumerate(mobile) if not allowed]
    chi = [indices[name] + 1 for name in ("THR:N", "THR:CA", "THR:CB", "THR:OG1")]
    old_policy = (
        ROOT
        / "output_nyl12_native_reference_20261023/native_frozen_scaffold_pilot_01/pilot_policy.json"
    )
    policy = json.loads(old_policy.read_text())
    policy.update(
        {
            "scope": "42atom native-state frozen-scaffold pilot; no FF/RESP acceptance",
            "frozen_atom_indices_0based": frozen,
            "chi1_constraint_indices_1based": chi,
            "atom_count": 42,
            "nprocshared": 16,
            "memory": "32GB",
            "cpu_affinity": "0-15, sixteen distinct physical cores",
            "scratch": "/data/scratch/gaussian_scratch/jpc/truncated_ti1_20261023_pilot01",
        }
    )
    for name, record in [
        ("atom_map.json", metadata),
        ("preflight.json", preflight),
        ("coordinate_provenance.json", provenance),
        ("pilot_policy.json", policy),
    ]:
        (OUT / name).write_text(json.dumps(record, indent=2) + "\n")
    if not geometry["accepted"]:
        raise ValueError(
            "Truncated geometry failed preflight; no Gaussian input staged"
        )
    input_lines = [
        "%chk=optimization.chk",
        "%mem=32GB",
        "%nprocshared=16",
        "#p wB97XD/def2TZVP Opt=(ModRedundant,Tight,MaxCycles=150,MaxStep=5)",
        " SCF=(Tight,XQC) Integral=UltraFine NoSymm",
        "",
        "Truncated R TI1 native-state pilot; no parameter acceptance",
        "",
        "0 1",
    ]
    input_lines += [
        f"{atom.GetSymbol():2s} {0 if allowed else -1:2d} "
        + " ".join(f"{v:.10f}" for v in point)
        for atom, point, allowed in zip(
            product.GetAtoms(), positions, mobile, strict=True
        )
    ]
    input_lines += ["", "D " + " ".join(map(str, chi)) + " F", ""]
    (OUT / "input.gjf").write_text("\n".join(input_lines) + "\n")
    shutil.copy2(BASE / "assess_pilot.py", OUT / "assess_pilot.py")
    (OUT / "view_seed.pml").write_text(
        "load native_seed.sdf, truncated_ti1\nset valence, 0\nshow sticks, truncated_ti1\nutil.cnc truncated_ti1\nzoom truncated_ti1\n"
    )
    wrapper = (
        "#!/usr/bin/env bash\nset -eo pipefail\nexport g16root=/usr/local\n"
        "source /usr/local/g16/bsd/g16.profile\n"
        f'export GAUSS_SCRDIR={policy["scratch"]}\nmkdir -- "$GAUSS_SCRDIR"\n'
        f"cd -- {OUT}\n/usr/bin/taskset -c 0-15 /usr/local/g16/g16 < input.gjf > optimization.log 2>&1\n"
        f"PYTHONPATH={ROOT}/src {ROOT}/.venv/bin/python {OUT}/assess_pilot.py\n"
    )
    (OUT / "run_gaussian.sh").write_text(wrapper)
    snapshot = OUT / "source"
    snapshot.mkdir()
    sources = [
        Path(__file__),
        BASE / "assess_pilot.py",
        SOURCE / "native_seed.sdf",
        SOURCE / "native_seed.xyz",
        SOURCE / "atom_map.json",
        old_policy,
        ROOT / "src/nyl12_covalent/gaussian_checks.py",
        ROOT / "src/nyl12_covalent/active_site_geometry.py",
    ]
    for i, source in enumerate(sources):
        shutil.copy2(source, snapshot / f"{i}_{source.name}")
    hashes = {
        str(source): hashlib.sha256(source.read_bytes()).hexdigest()
        for source in sources
    }
    (OUT / "source_manifest.json").write_text(json.dumps(hashes, indent=2) + "\n")
    LOGGER.info(
        "Truncated seed passes:42atoms25H,net0,R;30mobile12frozen;16physicalcores"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
