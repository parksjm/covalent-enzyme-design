"""Stage a short def2-SVP state-retention pilot for the approved donor cluster."""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
from pathlib import Path

import numpy as np
from rdkit import Chem

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = Path(__file__).resolve().parent


def main() -> None:
    """Require seed acceptance and retain explicit boundaries and proton checks."""
    seed = BASE / "seed_01"
    report = json.loads((seed / "preflight.json").read_text())
    if not report["chemical_inventory_passed"] or not report["geometry_screen_passed"]:
        raise ValueError("Unacceptable donor seed; Gaussian must not launch")
    molecule = Chem.SDMolSupplier(str(seed / "native_seed.sdf"), removeHs=False)[0]
    if molecule is None or molecule.GetNumAtoms() != 63:
        raise ValueError("Expected 63-atom donor cluster")
    xyz = np.asarray(
        [
            [float(v) for v in line.split()[1:]]
            for line in (seed / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    mapping = json.loads((seed / "atom_map.json").read_text())
    origins = mapping["atom_index_by_origin"]
    if xyz.shape != (63, 3) or not np.isfinite(xyz).all():
        raise ValueError("Invalid seed coordinates")
    source_policy = (
        ROOT
        / "output_nyl12_truncated_reference_20261023/truncated_pilot_01/pilot_policy.json"
    )
    policy = json.loads(source_policy.read_text())
    frozen = list(policy["frozen_atom_indices_0based"])
    frozen += [
        row["index"]
        for row in mapping["atoms"]
        if row["origin_id"].startswith(("NMA:", "ACM:")) and row["element"] != "H"
    ]
    if len(frozen) != 21 or len(set(frozen)) != 21:
        raise ValueError("Expected 12 core and 9 donor frozen heavy atoms")
    chi = [origins[name] + 1 for name in ("THR:N", "THR:CA", "THR:CB", "THR:OG1")]
    scratch = "/data/scratch/gaussian_scratch/jpc/oxyanion63_svp_20261023_pilot01"
    policy.update(
        {
            "scope": "63-atom short constrained state-retention pilot; no FF/RESP or minimum acceptance",
            "method": "wB97XD/def2SVP; no extra dispersion; gas phase",
            "atom_count": 63,
            "frozen_atom_indices_0based": frozen,
            "frozen_donor_policy": "All 9 amide heavy atoms fixed at source protein positions; donor H mobile. Represents held protein donor geometry, not a free cluster minimum.",
            "chi1_constraint_indices_1based": chi,
            "optimizer_coordinates": "ReadAllGIC Cartesian plus frozen chi1",
            "max_cycles": 20,
            "maxstep": 2,
            "water_present": True,
            "scratch": scratch,
            "environment_limit": "Explicit relay water and two held amide donors; remaining protein electrostatics absent. Neither state retention nor convergence validates a native enzyme minimum.",
        }
    )
    policy["prospective_checks"]["ammonium_proton_limits"] = {
        "N_H_minimum_A": 0.85,
        "N_H_maximum_A": 1.25,
        "OG1_H_minimum_A": 1.20,
    }
    out = BASE / "svp_pilot_01"
    out.mkdir(exist_ok=False)
    for name in (
        "native_seed.sdf",
        "native_seed.xyz",
        "atom_map.json",
        "coordinate_provenance.json",
        "preflight.json",
    ):
        shutil.copy2(seed / name, out / name)
    shutil.copy2(BASE / "assess_pilot.py", out / "assess_pilot.py")
    (out / "pilot_policy.json").write_text(json.dumps(policy, indent=2) + "\n")
    lines = [
        "%chk=optimization.chk",
        "%mem=32GB",
        "%nprocshared=16",
        "#p wB97XD/def2SVP Opt=(ReadAllGIC,Tight,MaxCycles=20,MaxStep=2)",
        " SCF=(Tight,XQC) Integral=UltraFine NoSymm",
        "",
        "R TI1 plus relay water and oxyanion donors; short state-retention pilot",
        "",
        "0 1",
    ]
    lines += [
        f"{atom.GetSymbol():2s} " + " ".join(f"{v:.10f}" for v in point)
        for atom, point in zip(molecule.GetAtoms(), xyz, strict=True)
    ]
    lines += [""]
    for i in range(63):
        for axis in "XYZ":
            option = "(Frozen)" if i in frozen else ""
            lines.append(f"C{axis}{i + 1}{option}={axis}({i + 1})")
    lines += ["CHI(Frozen)=D(" + ",".join(map(str, chi)) + ")", ""]
    (out / "input.gjf").write_text("\n".join(lines) + "\n")
    wrapper = f"""#!/usr/bin/env bash
set -eo pipefail
export g16root=/usr/local
source /usr/local/g16/bsd/g16.profile
export GAUSS_SCRDIR={scratch}
mkdir -- "$GAUSS_SCRDIR"
cd -- {out}
engine_status=0
if /usr/bin/taskset -c 0-15 /usr/local/g16/g16 < input.gjf > optimization.log 2>&1; then
    engine_status=0
else
    engine_status=$?
fi
assessment_status=0
if PYTHONPATH={ROOT}/src {ROOT}/.venv/bin/python {out}/assess_pilot.py; then
    assessment_status=0
else
    assessment_status=$?
fi
if [[ $engine_status -ne 0 ]]; then
    exit "$engine_status"
fi
exit "$assessment_status"
"""
    (out / "run_gaussian.sh").write_text(wrapper)
    snapshot = out / "source"
    snapshot.mkdir()
    paths = [
        Path(__file__),
        BASE / "assess_pilot.py",
        BASE / "prepare_oxyanion_cluster.py",
        BASE / "donor_fragments.py",
        source_policy,
        ROOT / "src/nyl12_covalent/gaussian_checks.py",
    ]
    paths += [
        seed / name
        for name in (
            "native_seed.sdf",
            "native_seed.xyz",
            "atom_map.json",
            "preflight.json",
        )
    ]
    for i, path in enumerate(paths):
        shutil.copy2(path, snapshot / f"{i:02d}_{path.name}")
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (out / "source_manifest.json").write_text(json.dumps(hashes, indent=2) + "\n")
    logging.info(
        "Staged63atom SVP pilot:21 frozen heavy,42 mobile,chi1 retained; all reactive bonds/protons free"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
