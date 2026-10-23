"""Prepare a water-free TI1 candidate with explicit Cartesian GIC constraints."""

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
SOURCE = ROOT / "output_nyl12_truncated_reference_20261023/truncated_pilot_01"
OUT = BASE / "waterfree_pilot_01"


def main() -> None:
    """Remove only source water, retain coordinates/state, and stage the GIC input."""
    molecule = Chem.SDMolSupplier(str(SOURCE / "native_seed.sdf"), removeHs=False)[0]
    if molecule is None or molecule.GetNumAtoms() != 42:
        raise ValueError(
            "Expected original 42-atom starting graph, not failed endpoint"
        )
    metadata = json.loads((SOURCE / "atom_map.json").read_text())
    old_origins = metadata["atom_index_by_origin"]
    inverse = {index: origin for origin, index in old_origins.items()}
    xyz = np.asarray(
        [
            [float(v) for v in line.split()[1:]]
            for line in (SOURCE / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    remove = sorted(
        index for origin, index in old_origins.items() if origin.startswith("WAT:")
    )
    if (
        set(old_origins) & {"WAT:O", "WAT:H1", "WAT:H2"}
        != {"WAT:O", "WAT:H1", "WAT:H2"}
        or len(remove) != 3
    ):
        raise ValueError("Source water inventory differs from explicit H2O")
    keep = [i for i in range(42) if i not in remove]
    rw = Chem.RWMol(molecule)
    for idx in reversed(remove):
        rw.RemoveAtom(idx)
    product = rw.GetMol()
    Chem.SanitizeMol(product)
    product.RemoveAllConformers()
    conformer = Chem.Conformer(39)
    conformer.Set3D(True)
    for idx, old_idx in enumerate(keep):
        conformer.SetAtomPosition(idx, tuple(float(v) for v in xyz[old_idx]))
    product.AddConformer(conformer)
    Chem.AssignStereochemistry(product, cleanIt=True, force=True)
    indices = {inverse[old]: new for new, old in enumerate(keep)}
    formula = CalcMolFormula(product)
    if (product.GetNumAtoms(), formula, Chem.GetFormalCharge(product)) != (
        39,
        "C10H23N3O3",
        0,
    ):
        raise ValueError("Water-free inventory/formula/charge mismatch")
    geometry = validate_geometry(
        product, [a.GetAtomicNum() for a in product.GetAtoms()], xyz[keep]
    )
    expected = {
        indices[name]: code
        for name, code in [("THR:CA", "S"), ("THR:CB", "R"), ("PA6:C54", "R")]
    }
    if geometry["observed_stereocenters"] != expected or not geometry["accepted"]:
        raise ValueError("Water-free source geometry/stereochemistry failed")
    OUT.mkdir(exist_ok=False)
    with Chem.SDWriter(str(OUT / "native_seed.sdf")) as writer:
        writer.write(product)
    lines = [
        "39",
        "Water-free native-state R TI1 starting candidate; no optimized-state claim",
    ]
    lines += [
        f"{atom.GetSymbol():2s} " + " ".join(f"{v:.10f}" for v in point)
        for atom, point in zip(product.GetAtoms(), xyz[keep], strict=True)
    ]
    (OUT / "native_seed.xyz").write_text("\n".join(lines) + "\n")
    rows = [{**metadata["atoms"][old], "index": new} for new, old in enumerate(keep)]
    mapping = {
        "formula": formula,
        "formal_charge": 0,
        "attacked_carbon_cip": "R",
        "atom_index_by_origin": indices,
        "atoms": rows,
    }
    policy = json.loads((SOURCE / "pilot_policy.json").read_text())
    frozen_origins = [inverse[i] for i in policy["frozen_atom_indices_0based"]]
    frozen = [indices[origin] for origin in frozen_origins]
    chi = [indices[name] + 1 for name in ("THR:N", "THR:CA", "THR:CB", "THR:OG1")]
    policy.update(
        {
            "scope": "water-free39 native-state constrained geometry candidate; no FF/RESP acceptance",
            "atom_count": 39,
            "frozen_atom_indices_0based": frozen,
            "mobile_heavy_origins": ["THR:OG1", "PA6:C54", "PA6:O49", "PA6:N59"],
            "chi1_constraint_indices_1based": chi,
            "optimizer_coordinates": "ReadAllGIC: all Cartesian coordinates plus chi1",
            "maxstep": 2,
            "water_present": False,
            "scratch": "/data/scratch/gaussian_scratch/jpc/waterfree39_gic_20261023_attempt01",
        }
    )
    preflight = {
        "accepted": True,
        "geometry": geometry,
        "formula": formula,
        "atom_count": 39,
        "explicit_H": 23,
        "formal_charge": 0,
        "removed_origins": [inverse[i] for i in remove],
        "retained_coordinate_max_displacement_A": 0.0,
        "source": str(SOURCE / "native_seed.xyz"),
        "source_is_unoptimized_start_not_failed_endpoint": True,
        "forcefield_accepted": False,
    }
    provenance = [
        {
            "index": new,
            "origin_id": inverse[old],
            "source_42atom_index_0based": old,
            "coordinates_changed": False,
        }
        for new, old in enumerate(keep)
    ]
    for name, record in [
        ("atom_map.json", mapping),
        ("pilot_policy.json", policy),
        ("preflight.json", preflight),
        ("coordinate_provenance.json", provenance),
    ]:
        (OUT / name).write_text(json.dumps(record, indent=2) + "\n")
    header = [
        "%chk=optimization.chk",
        "%mem=32GB",
        "%nprocshared=16",
        "#p wB97XD/def2TZVP Opt=(ReadAllGIC,Tight,MaxCycles=150,MaxStep=2)",
        " SCF=(Tight,XQC) Integral=UltraFine NoSymm",
        "",
        "Water-free39 native-state TI1; identical physical freezes via GIC",
        "",
        "0 1",
    ]
    header += [
        f"{atom.GetSymbol():2s} " + " ".join(f"{v:.10f}" for v in point)
        for atom, point in zip(product.GetAtoms(), xyz[keep], strict=True)
    ]
    header += [""]
    for i in range(39):
        for axis in "XYZ":
            option = "(Frozen)" if i in frozen else ""
            header.append(f"C{axis}{i + 1}{option}={axis}({i + 1})")
    header += ["CHI(Frozen)=D(" + ",".join(map(str, chi)) + ")", ""]
    (OUT / "input.gjf").write_text("\n".join(header) + "\n")
    shutil.copy2(BASE / "assess_pilot.py", OUT / "assess_pilot.py")
    (OUT / "view_seed.pml").write_text(
        "load native_seed.sdf, waterfree_ti1\nset valence, 0\nshow sticks, waterfree_ti1\nutil.cnc waterfree_ti1\nzoom waterfree_ti1\n"
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
        SOURCE / "pilot_policy.json",
        ROOT / "src/nyl12_covalent/gaussian_checks.py",
    ]
    for i, source in enumerate(sources):
        shutil.copy2(source, snapshot / f"{i}_{source.name}")
    hashes = {
        str(source): hashlib.sha256(source.read_bytes()).hexdigest()
        for source in sources
    }
    (OUT / "source_manifest.json").write_text(json.dumps(hashes, indent=2) + "\n")
    logging.info("Prepared39atoms23H net0 R;12frozen,27mobile,chi1 retained via GIC")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
