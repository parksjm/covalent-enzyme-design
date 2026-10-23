"""Stage a frozen-scaffold native TI1 pilot; never fit parameters."""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
from pathlib import Path

import numpy as np
from rdkit import Chem

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = ROOT / "output_nyl12_native_reference_20261023"
LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Check the source seed and stage one explicit, nonautomatic pilot."""
    seed = BASE / "native_seed_01"
    molecules = list(Chem.SDMolSupplier(str(seed / "native_seed.sdf"), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("Expected one explicit-H native seed")
    molecule = molecules[0]
    preflight = json.loads((seed / "preflight.json").read_text())
    if not preflight["accepted"]:
        raise ValueError("Native seed preflight failed; do not stage Gaussian")
    if molecule.GetNumAtoms() != 103 or Chem.GetFormalCharge(molecule) != 0:
        raise ValueError("Native inventory/charge mismatch")
    metadata = json.loads((seed / "atom_map.json").read_text())
    origins = metadata["atom_index_by_origin"]
    xyz_lines = (seed / "native_seed.xyz").read_text().splitlines()[2:]
    xyz = np.asarray([[float(v) for v in line.split()[1:]] for line in xyz_lines])
    if xyz.shape != (103, 3) or not np.isfinite(xyz).all():
        raise ValueError("Malformed high-precision seed coordinates")
    mobile_origins = {"THR:OG1", "PA6:C54", "PA6:O49", "PA6:N59", "WAT:O"}
    mobile = [atom.GetAtomicNum() == 1 for atom in molecule.GetAtoms()]
    for origin in mobile_origins:
        mobile[origins[origin]] = True
    frozen = [i for i, allowed in enumerate(mobile) if not allowed]
    out = BASE / "native_frozen_scaffold_pilot_01"
    out.mkdir(exist_ok=False)
    for path in seed.iterdir():
        if path.is_file():
            shutil.copy2(path, out / path.name)
    lines = [
        "%chk=optimization.chk",
        "%mem=32GB",
        "%nprocshared=8",
        "#p wB97XD/def2TZVP Opt=(ModRedundant,Tight,MaxCycles=150,MaxStep=5)",
        " SCF=(Tight,XQC) Integral=UltraFine NoSymm",
        "",
        "Native TI1 frozen-scaffold pilot; not a charge-fit reference",
        "",
        "0 1",
    ]
    for atom, point, allowed in zip(molecule.GetAtoms(), xyz, mobile, strict=True):
        flag = 0 if allowed else -1
        lines.append(
            f"{atom.GetSymbol():2s} {flag:2d} " + " ".join(f"{v:.10f}" for v in point)
        )
    chi = [origins[name] + 1 for name in ("THR:N", "THR:CA", "THR:CB", "THR:OG1")]
    lines.extend(["", "D " + " ".join(map(str, chi)) + " F", ""])
    (out / "input.gjf").write_text("\n".join(lines) + "\n")
    scratch = "/data/scratch/gaussian_scratch/jpc/native_ti1_20261023_pilot01"
    wrapper = (
        "#!/usr/bin/env bash\nset -eo pipefail\n"
        "export g16root=/usr/local\nsource /usr/local/g16/bsd/g16.profile\n"
        f'export GAUSS_SCRDIR={scratch}\nmkdir -- "$GAUSS_SCRDIR"\n'
        f"cd -- {out}\n/usr/local/g16/g16 < input.gjf > optimization.log 2>&1\n"
        f"PYTHONPATH={ROOT}/src {ROOT}/.venv/bin/python {out}/assess_pilot.py\n"
    )
    (out / "run_gaussian.sh").write_text(wrapper)
    shutil.copy2(BASE / "assess_pilot.py", out / "assess_pilot.py")
    code = out / "source"
    code.mkdir()
    for name in ("gaussian_checks.py", "ti1_reference_graph.py"):
        shutil.copy2(ROOT / "src/nyl12_covalent" / name, code / name)
    shutil.copy2(Path(__file__), code / Path(__file__).name)
    policy = {
        "scope": "conditional local state-retention pilot; no FF/RESP/native acceptance",
        "method": "wB97XD/def2TZVP; no extra dispersion; gas phase",
        "charge": 0,
        "multiplicity": 1,
        "mobile_heavy_origins": sorted(mobile_origins),
        "all_hydrogens_mobile": True,
        "frozen_atom_indices_0based": frozen,
        "chi1_constraint_indices_1based": chi,
        "torsion_source": "AF3 seed, explicitly distinct from n12c torsion reference",
        "reactive_bond_restraints": False,
        "environment_limit": "Frozen scaffold retains conformation; missing protein electrostatics are not restored.",
        "prospective_checks": {
            "frozen_atom_tolerance_A": 0.0001,
            "chi1_tolerance_deg": 0.1,
            "reactive_bond_windows_A": {
                "OG1_C54": [1.30, 1.75],
                "C54_O49": [1.25, 1.65],
                "C54_N59": [1.25, 1.70],
                "C54_C62": [1.35, 1.75],
            },
            "oxyanion_H_minimum_A": 1.20,
        },
        "scratch": scratch,
    }
    (out / "pilot_policy.json").write_text(json.dumps(policy, indent=2) + "\n")
    files = [p for p in out.rglob("*") if p.is_file()]
    manifest = {
        str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in files
    }
    (out / "staging_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Staged %s; %d mobile atoms, %d frozen", out, sum(mobile), len(frozen))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
