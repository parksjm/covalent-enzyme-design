"""Extract, minimize, and geometry-audit the selected ASH274 E-S frame."""

from __future__ import annotations

import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone
from typing import Any

import numpy as np
import parmed

ROOT = Path(__file__).resolve().parents[5]
COMPARISON = ROOT / "qmmm/pa6_flipped_model30/07_ash274_comparison"
SYSTEM = COMPARISON / "00_system"
REPLICA = COMPARISON / "03_md/ash274_replica_1"
SELECTION = (
    ROOT
    / "qmmm/pa6_flipped_model30/09_propka_static_tier_a/screen_attempt1/manifest.json"
)
OUTPUT = Path(__file__).resolve().parent
TOPOLOGY = SYSTEM / "ash274_solvated.prmtop"
MD_TOPOLOGY = REPLICA / "system.prmtop"
TRAJ10 = REPLICA / "05_production_10ns.nc"
TRAJ90 = REPLICA / "06_production_90ns.nc"
CPPTRAJ = Path("/data/jpc/code/ambertools26/venv/bin/cpptraj")
ENGINE = Path("/data/jpc/code/pmemd26/bin/pmemd.cuda_SPFP")
RUN_SCRIPT = OUTPUT / "run.sh"
LOGGER = logging.getLogger(__name__)
EXPECTED_CHARGE = -0.004004946642281126
FRAME = 3492
LOCAL_FRAME = 2492
FRAME_TIME_PS = 35420.0
PROTEIN_SITES = (
    (929, 968, 970),
    (1025, 1064, 1066),
    (1121, 1160, 1162),
    (1217, 1256, 1258),
)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a retained input or output."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atom_index(structure: Any, residue_number: int, name: str) -> int:
    """Resolve one atom by one-based Amber residue number and exact atom name."""
    if residue_number < 1 or residue_number > len(structure.residues):
        raise ValueError(f"Residue outside topology: {residue_number}")
    atoms = [
        atom.idx
        for atom in structure.residues[residue_number - 1].atoms
        if atom.name == name
    ]
    if len(atoms) != 1:
        raise ValueError(
            f"Expected one atom {residue_number}:{name}, found {len(atoms)}"
        )
    return int(atoms[0])


def validate_microstate(structure: Any, acceptance: dict[str, Any]) -> None:
    """Require ASP272/ASH274, neutral Thr N termini, and the recorded system charge."""
    if acceptance.get("protonated_residue") != 274 or not acceptance.get("neutral_Thr"):
        raise ValueError(
            "ASH274 preparation acceptance does not match requested microstate"
        )
    if abs(sum(atom.charge for atom in structure.atoms) - EXPECTED_CHARGE) > 1e-5:
        raise ValueError("ASH274 solvated topology charge differs from accepted charge")
    if len(structure.atoms) != 104934 or len(structure.residues) != 29781:
        raise ValueError("Unexpected full-solvent ASH274 topology size")
    for thr, asp272, ash274 in PROTEIN_SITES:
        residues = [structure.residues[number - 1] for number in (thr, asp272, ash274)]
        if [residue.name for residue in residues] != ["THR", "ASP", "ASH"]:
            raise ValueError(
                f"Unexpected catalytic residue mapping at {thr}/{asp272}/{ash274}"
            )
        names = {atom.name: atom for atom in residues[0].atoms}
        if "H3" in names or not {"N", "H1", "H2", "OG1"}.issubset(names):
            raise ValueError(f"Neutral Thr N atom inventory mismatch at residue {thr}")
        for name, charge in (("N", -0.4844), ("H1", 0.1224), ("H2", 0.1224)):
            if abs(names[name].charge - charge) > 1e-6:
                raise ValueError(f"Neutral Thr charge mismatch at {thr}:{name}")
        ash_names = {atom.name: atom for atom in residues[2].atoms}
        if (
            "HD2" not in ash_names
            or ash_names["HD2"] not in ash_names["OD2"].bond_partners
        ):
            raise ValueError(
                f"ASH274 HD2-OD2 connectivity mismatch at residue {ash274}"
            )
        if "HD2" in {atom.name for atom in residues[1].atoms}:
            raise ValueError(f"ASP272 unexpectedly protonated at residue {asp272}")


def geometry(structure: Any) -> dict[str, float]:
    """Measure the catalytic attack, oxyanion contacts, and signed phi."""
    sys.path.insert(0, str(ROOT / "src"))
    from nyl12_covalent.active_site_geometry import (
        angle_degrees,
        projected_herschlag_phi,
    )

    indices = {
        "og1": atom_index(structure, 929, "OG1"),
        "cb": atom_index(structure, 929, "CB"),
        "phe_n": atom_index(structure, 113, "N"),
        "asn_nd2": atom_index(structure, 186, "ND2"),
        "carbon": atom_index(structure, 1314, "C"),
        "oxygen": atom_index(structure, 1314, "O"),
        "amide_n": atom_index(structure, 1315, "N"),
        "adjacent_c": atom_index(structure, 1314, "C5"),
    }
    xyz = np.asarray(structure.coordinates, dtype=float)
    point = {key: xyz[index] for key, index in indices.items()}

    def distance(left: str, right: str) -> float:
        """Return the distance between two mapped active-site atoms."""
        return float(np.linalg.norm(point[left] - point[right]))

    return {
        "d1": distance("og1", "carbon"),
        "d2": distance("phe_n", "oxygen"),
        "d3": distance("asn_nd2", "oxygen"),
        "attack": angle_degrees(point["og1"], point["carbon"], point["oxygen"]),
        "projected_phi": projected_herschlag_phi(
            tuple(point["og1"]),
            tuple(point["oxygen"]),
            tuple(point["carbon"]),
            tuple(point["amide_n"]),
            tuple(point["adjacent_c"]),
        ),
        "scissile_CO_A": distance("carbon", "oxygen"),
        "scissile_CN_A": distance("carbon", "amide_n"),
    }


def selected_record() -> dict[str, Any]:
    """Return the declared frame record after checking frame and source hashes."""
    records = json.loads(SELECTION.read_text())
    matches = [item for item in records if item["case"] == "ash274_replica_1"]
    if len(matches) != 1:
        raise ValueError("Expected one ASH274 replica 1 selection record")
    record = matches[0]
    if (
        record["state"] != "ASP272_ASH274"
        or int(record["frame"]) != FRAME
        or float(record["geometry"]["time_ps"]) != FRAME_TIME_PS
        or record["geometry"]["tier"] != "A"
    ):
        raise ValueError(
            "Saved selection does not identify requested Tier A frame 3492"
        )
    for source, expected in record["source_sha256"].items():
        path = Path(source)
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Selected-frame source hash mismatch: {path}")
    if sha256(MD_TOPOLOGY) != sha256(TOPOLOGY):
        raise ValueError("Production topology differs from accepted ASH274 topology")
    return record


def main() -> None:
    """Extract, minimize, and acceptance-audit the selected full-solvent frame."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if any(
        (OUTPUT / name).exists()
        for name in ("status.json", "solvated.rst7", "optimized.rst7")
    ):
        raise FileExistsError("This ASH274 selected-MM attempt already has outputs")
    selection = selected_record()
    acceptance = json.loads((SYSTEM / "preparation_acceptance.json").read_text())
    topology = parmed.load_file(str(TOPOLOGY))
    validate_microstate(topology, acceptance)
    inputs = [
        TOPOLOGY,
        SYSTEM / "preparation_acceptance.json",
        SYSTEM / "preparation_manifest.json",
        SELECTION,
        TRAJ10,
        TRAJ90,
        MD_TOPOLOGY,
        CPPTRAJ,
        ENGINE,
        RUN_SCRIPT,
        Path(__file__),
    ]
    extract_input = OUTPUT / "extract.cpptraj"
    extracted = OUTPUT / "solvated.rst7"
    extract_input.write_text(
        f"parm {TOPOLOGY.resolve()}\n"
        f"trajin {TRAJ90.resolve()} {LOCAL_FRAME} {LOCAL_FRAME} 1\n"
        "autoimage anchor :1313-1316\n"
        f"trajout {extracted.resolve()} restart\nrun\n"
    )
    environment = os.environ.copy()
    environment["AMBERHOME"] = "/data/jpc/code/ambertools26/venv"
    environment["LD_LIBRARY_PATH"] = "/data/jpc/code/ambertools26/venv/lib"
    environment["OMP_NUM_THREADS"] = "1"
    with (OUTPUT / "cpptraj.log").open("w") as log:
        subprocess.run(
            [str(CPPTRAJ), "-i", str(extract_input)],
            cwd=OUTPUT,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    starting = parmed.load_file(str(TOPOLOGY), xyz=str(extracted))
    if len(starting.atoms) != len(topology.atoms):
        raise ValueError("Extracted restart does not contain the full solvated system")
    before_geometry = geometry(starting)
    for key in ("d1", "d2", "d3", "attack", "projected_phi"):
        expected = float(selection["geometry"][key if key != "attack" else "attack"])
        if abs(before_geometry[key] - expected) > 0.02:
            raise ValueError(f"Extracted frame geometry mismatch for {key}")
    from nyl12_covalent.modular_es_snapshot_selection import classify_geometry

    if (
        classify_geometry(before_geometry) != "A"
        or before_geometry["projected_phi"] <= 0
    ):
        raise ValueError("Extracted frame fails positive-phi Tier A gate")
    mdin = OUTPUT / "minimize.mdin"
    mdin.write_text(
        "ASH274 ASP272/ASH274 neutral-N E-S replica_1_tier_A unrestrained MM minimization\n"
        "&cntrl\n imin=1, maxcyc=10000, ncyc=5000, ntmin=1,\n"
        " ntb=1, igb=0, cut=10.0, ntc=1, ntf=1,\n"
        " ntr=0, nmropt=0, drms=0.01,\n"
        " ntpr=25, ntwx=25, ntwr=250, ioutfm=1, ntxo=2,\n/\n"
    )
    command = [
        str(ENGINE),
        "-O",
        "-i",
        str(mdin),
        "-p",
        str(TOPOLOGY),
        "-c",
        str(extracted),
        "-o",
        str(OUTPUT / "minimized.mdout"),
        "-r",
        str(OUTPUT / "optimized.rst7"),
        "-x",
        str(OUTPUT / "optimization.nc"),
        "-inf",
        str(OUTPUT / "mdinfo"),
    ]
    process = {
        "command": command,
        "gpu_visibility": "CUDA_VISIBLE_DEVICES=0",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "status": "running",
    }
    (OUTPUT / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    run_environment = environment | {"CUDA_VISIBLE_DEVICES": "0"}
    with (OUTPUT / "engine.log").open("w") as log:
        result = subprocess.run(
            command,
            cwd=OUTPUT,
            env=run_environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    process.update(
        returncode=result.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    (OUTPUT / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    if result.returncode != 0:
        (OUTPUT / "status.json").write_text(
            json.dumps({"exit_status": result.returncode}, indent=2) + "\n"
        )
        raise RuntimeError(f"PMEMD minimization exited {result.returncode}")
    mdout = (OUTPUT / "minimized.mdout").read_text(errors="replace")
    if "FINAL RESULTS" not in mdout or any(
        term in mdout for term in ("NaN", "Infinity", "FATAL ERROR")
    ):
        raise ValueError(
            "PMEMD output lacks normal final results or reports a numerical failure"
        )
    final = parmed.load_file(str(TOPOLOGY), xyz=str(OUTPUT / "optimized.rst7"))
    validate_microstate(parmed.load_file(str(TOPOLOGY)), acceptance)
    after_geometry = geometry(final)
    tier = classify_geometry(after_geometry)
    if tier not in {"A", "B"} or after_geometry["projected_phi"] <= 0:
        raise ValueError("Minimized frame failed existing Tier A/B positive-phi gate")
    source_hashes = {str(path.resolve()): sha256(path) for path in inputs}
    output_hashes = {
        name: sha256(OUTPUT / name)
        for name in (
            "extract.cpptraj",
            "cpptraj.log",
            "solvated.rst7",
            "minimize.mdin",
            "minimized.mdout",
            "optimized.rst7",
            "optimization.nc",
            "engine.log",
        )
    }
    report = {
        "case": "replica_1_tier_A",
        "microstate": {
            "ASP272": "ASP",
            "ASP274": "ASH(HD2 bonded to OD2)",
            "Thr233_N": "neutral NH2",
        },
        "source_frame": {
            "global_frame": FRAME,
            "local_frame_90ns": LOCAL_FRAME,
            "time_ps": FRAME_TIME_PS,
            "tier": "A",
            "geometry": selection["geometry"],
        },
        "atom_mapping": {
            "thr233_OG1": "Amber residue 929 OG1",
            "oxyanion_C": "Amber residue 1314 C",
            "oxyanion_O": "Amber residue 1314 O",
            "scissile_N": "Amber residue 1315 N",
            "adjacent_C": "Amber residue 1314 C5",
            "oxyanion_partner_Phe_N": "Amber residue 113 N",
            "oxyanion_partner_Asn_ND2": "Amber residue 186 ND2",
        },
        "system": {
            "atom_count": len(topology.atoms),
            "residue_count": len(topology.residues),
            "total_charge_e": sum(atom.charge for atom in topology.atoms),
            "full_solvent_retained": True,
        },
        "method": "Amber PMEMD26 CUDA SPFP; unrestrained periodic MM minimization; maxcyc=10000, ncyc=5000, drms=0.01",
        "geometry_before_mm": before_geometry,
        "geometry_after_mm": after_geometry,
        "geometric_acceptance": {
            "accepted": True,
            "tier_after_mm": tier,
            "positive_projected_phi": True,
            "criterion": "existing modular Tier A/B geometry classifier and positive projected Herschlag phi",
        },
        "force_convergence": "not asserted beyond Amber drms criterion; preparation minimization only",
        "source_sha256": source_hashes,
        "output_sha256": output_hashes,
        "status": "MM minimization complete; no QM/MM launched",
    }
    (OUTPUT / "preparation.json").write_text(json.dumps(report, indent=2) + "\n")
    (OUTPUT / "status.json").write_text(
        json.dumps(
            {
                "exit_status": 0,
                "geometric_acceptance": report["geometric_acceptance"],
                "qm_mm_launched": False,
            },
            indent=2,
        )
        + "\n"
    )
    LOGGER.info(
        "ASH274 MM minimization complete; post-MM tier %s; positive phi retained", tier
    )


if __name__ == "__main__":
    main()
