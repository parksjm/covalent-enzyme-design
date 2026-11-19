#!/usr/bin/env python3
"""Build neutral-beta model30 with released modular PA6 and minimize only."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import re
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import parmed  # noqa: E402
from nyl12_covalent.raw_pa6_modular import convert_heavy_complex  # noqa: E402
from nyl12_covalent.modular_pa6_md_setup import (  # noqa: E402
    _box_volume,
    _ion_counts,
    _require_clean_leap,
    _run_tleap,
    _write_leap_input,
    _write_leap_pdb,
)
from nyl12_covalent.modular_pa6_topology import (  # noqa: E402
    sha256,
    validate_release,
    validate_topology,
)

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Prepare a checked solvated topology and run unrestrained GPU minimization."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "pdb", "mol2", "mapping", "release", "tleap", "engine"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    system = args.root / "00_system"
    minimization = args.root / "01_mm_minimization"
    system.mkdir()
    minimization.mkdir()
    release = validate_release(args.release)
    for source, name in (
        (args.pdb, "raw_model30.pdb"),
        (args.mol2, "whole_l4_source.mol2"),
        (args.mapping, "whole_l4_context_mapping.json"),
    ):
        shutil.copy2(source, system / name)
    conversion = convert_heavy_complex(
        args.pdb,
        args.mol2,
        args.mapping,
        args.release,
        system / "complex_heavy_modular.pdb",
    )
    (system / "conversion.json").write_text(json.dumps(conversion, indent=2) + "\n")
    _write_leap_pdb(
        system / "complex_heavy_modular.pdb", system / "input_complex_leap.pdb"
    )
    for source in args.release.iterdir():
        if source.is_file():
            shutil.copy2(source, system / source.name)
    _write_leap_input(system / "leap_count.in", final_ions=None)
    LOGGER.info("Building neutral-beta ff19SB/GAFF2/TIP3P model30")
    _require_clean_leap(
        _run_tleap(args.tleap, system, "leap_count.in", "leap_count.log")
    )
    dry = validate_topology(system / "md_dry.prmtop", args.release, solvated=False)
    volume = _box_volume(system / "solvated_noions.inpcrd")
    sodium, chloride, pairs = _ion_counts(dry["total_charge"], volume)
    _write_leap_input(system / "leap_final.in", final_ions=(sodium, chloride))
    _require_clean_leap(
        _run_tleap(args.tleap, system, "leap_final.in", "leap_final.log")
    )
    solvated = validate_topology(
        system / "md_solvated.prmtop", args.release, solvated=True
    )
    topology = parmed.load_file(str(system / "md_solvated.prmtop"))
    beta = {}
    for residue in (929, 1025, 1121, 1217):
        atoms = {atom.name: atom for atom in topology.residues[residue - 1].atoms}
        if "H3" in atoms or any(
            abs(atoms[name].charge - charge) > 1e-6
            for name, charge in (("N", -0.4844), ("H1", 0.1224), ("H2", 0.1224))
        ):
            raise ValueError(f"Neutral Thr N terminus validation failed at {residue}")
        beta[residue] = {name: atoms[name].charge for name in ("N", "H1", "H2")}
    acceptance = {
        "status": "topology_accepted_not_geometry_optimized",
        "release_id": release["release_id"],
        "dry": dry,
        "solvated": solvated,
        "salt_pairs_150mM": pairs,
        "sodium": sodium,
        "chloride": chloride,
        "box_volume_A3": volume,
        "beta_NH2": beta,
        "protonation": "Four neutral beta Thr N termini; other ff19SB LEaP standard defaults; raw HIS resolves to HIE; standard ASP/GLU anionic LYS/ARG cationic; alpha N termini charged; protein C termini charged.",
        "histidine_inventory": [
            {"residue": i + 1, "name": r.name}
            for i, r in enumerate(topology.residues[:1312])
            if r.name in ("HIE", "HID", "HIP", "HIS")
        ],
    }
    (system / "acceptance.json").write_text(json.dumps(acceptance, indent=2) + "\n")
    mdin = minimization / "min.in"
    mdin.write_text(
        "Model30 fully unrestrained MM minimization\n&cntrl\n imin=1, maxcyc=10000, ncyc=1000, drms=0.001,\n ntb=1, cut=10.0, ntr=0, nmropt=0, ntc=1, ntf=1,\n ntpr=100, ntwx=10, ntwr=1000, ioutfm=1, ntxo=2,\n/\n"
    )
    command = [
        str(args.engine),
        "-O",
        "-i",
        str(mdin),
        "-o",
        str(minimization / "min.out"),
        "-p",
        str(system / "md_solvated.prmtop"),
        "-c",
        str(system / "md_solvated.inpcrd"),
        "-r",
        str(minimization / "min.rst7"),
        "-x",
        str(minimization / "min.nc"),
        "-inf",
        str(minimization / "min.mdinfo"),
    ]
    (minimization / "command.json").write_text(
        json.dumps(
            {
                "command": command,
                "input_sha256": {
                    str(p): sha256(p)
                    for p in (
                        mdin,
                        system / "md_solvated.prmtop",
                        system / "md_solvated.inpcrd",
                        args.engine,
                    )
                },
                "constraints": "None: ntr=0 nmropt=0 ntc=1 ntf=1",
                "convergence": "drms=0.001 with maxcyc=10000; completion does not guarantee convergence",
            },
            indent=2,
        )
        + "\n"
    )
    LOGGER.info("Starting fully unrestrained model30 MM minimization")
    with (minimization / "engine.log").open("w") as handle:
        result = subprocess.run(
            command,
            cwd=minimization,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    (minimization / "process_acceptance.json").write_text(
        json.dumps({"returncode": result.returncode}, indent=2) + "\n"
    )
    if result.returncode:
        raise RuntimeError(
            f"pmemd.cuda exited {result.returncode}; inspect engine.log and min.out"
        )
    endpoint = parmed.load_file(
        str(system / "md_solvated.prmtop"), xyz=str(minimization / "min.rst7")
    )
    endpoint.save(str(minimization / "endpoint_solvated.pdb"), overwrite=False)
    endpoint.strip(":WAT,Na+,Cl-")
    endpoint.save(str(minimization / "endpoint_complex.pdb"), overwrite=False)
    output = (minimization / "min.out").read_text()
    records = re.findall(
        r"^\s*\d+\s+[-+0-9.E]+\s+([-+0-9.E]+)\s+[-+0-9.E]+\s+\S+\s+\d+",
        output,
        re.MULTILINE,
    )
    if not records:
        raise ValueError(
            "Cannot extract endpoint RMS gradient from minimization output"
        )
    gradient = float(records[-1])
    (minimization / "scientific_acceptance.json").write_text(
        json.dumps(
            {
                "final_rms_gradient": gradient,
                "drms_target": 0.001,
                "gradient_converged": gradient <= 0.001,
                "status": "MM_endpoint_requires_phi_tier_native_screen_before_further_work",
                "endpoint_pdb_sha256": sha256(minimization / "endpoint_complex.pdb"),
            },
            indent=2,
        )
        + "\n"
    )
    LOGGER.info(
        "MM endpoint exported; RMS gradient %.6g; no MD or QM/MM launched", gradient
    )


if __name__ == "__main__":
    main()
