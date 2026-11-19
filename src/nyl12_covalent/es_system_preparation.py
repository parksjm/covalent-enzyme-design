"""Prepare a source-pinned solvated Nyl12 PA6 E-S system with AmberTools."""

from __future__ import annotations

import hashlib
import json
import math
import os
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Any

from parmed.amber import AmberParm


BETA_N_TERMINI = (929, 1025, 1121, 1217)
CHAIN_ENDS = (232, 464, 696, 928, 1024, 1120, 1216, 1312)
N_TERMINI = (1, 233, 465, 697, 929, 1025, 1121, 1217)
LIGAND_RESIDUE = 1313
SALT_MOLAR = 0.150


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def off_atom_names(path: Path) -> set[str]:
    """Read the atom-name inventory from a LIG Amber OFF library."""
    names: set[str] = set()
    active = False
    for line in path.read_text().splitlines():
        if line.startswith("!entry.LIG.unit.atoms table"):
            active = True
            continue
        if line.startswith("!"):
            active = False
        elif active and line.strip():
            names.add(shlex.split(line)[0])
    if not names:
        raise ValueError(f"No LIG atoms found in {path}")
    return names


def write_leap_pdb(source: Path, destination: Path, ligand_names: set[str]) -> None:
    """Write sequential Amber numbering, explicit termini, and one PA6 ligand."""
    output = [f"REMARK source {source.resolve()}"]
    observed_ligand: list[str] = []
    terminators = 0
    for line in source.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = int(line[22:26])
        if residue > LIGAND_RESIDUE:
            continue
        atom = line[12:16].strip()
        if residue in N_TERMINI and atom == "H":
            atom = "H1"
        if residue == LIGAND_RESIDUE:
            observed_ligand.append(atom)
        output.append(f"{line[:12]}{atom:>4}{line[16:]}")
        if residue in CHAIN_ENDS and atom == "OXT":
            output.append("TER")
            terminators += 1
    if terminators != len(CHAIN_ENDS):
        raise ValueError(f"Expected {len(CHAIN_ENDS)} protein chain terminators")
    if len(observed_ligand) != len(set(observed_ligand)):
        raise ValueError("Duplicate PA6 ligand atom names")
    if set(observed_ligand) != ligand_names:
        raise ValueError("PA6 ligand atoms do not match the OFF library")
    output.extend(("TER", "END"))
    destination.write_text("\n".join(output) + "\n")


def leap_input(params: Path, sodium: int | None, chloride: int | None) -> str:
    """Return ff19SB/GAFF1/TIP3P LEaP input with neutral beta termini."""
    neutralization = "\n".join(
        f"remove complex complex.{residue}.H3\n"
        f"set complex.{residue}.N charge -0.4844\n"
        f"set complex.{residue}.H1 charge 0.1224\n"
        f"set complex.{residue}.H2 charge 0.1224"
        for residue in BETA_N_TERMINI
    )
    stem = "solvated_noions"
    ions = ""
    if sodium is not None and chloride is not None:
        stem = "md_solvated"
        ions = f"addIonsRand complex Na+ {sodium} Cl- {chloride} 6.0\n"
    return f"""source leaprc.protein.ff19SB
source leaprc.gaff
source leaprc.water.tip3p
loadamberparams {params / "lig.frcmod"}
loadoff {params / "lig.lib"}
complex = loadpdb input_complex_min_leap.pdb
{neutralization}
check complex
charge complex
saveamberparm complex md_dry.prmtop md_dry.inpcrd
savepdb complex md_dry.pdb
solvateOct complex TIP3PBOX 12.0
{ions}check complex
charge complex
saveamberparm complex {stem}.prmtop {stem}.inpcrd
savepdb complex {stem}.pdb
quit
"""


def _run_tleap(tleap: Path, input_name: str, output: Path, log_name: str) -> None:
    environment = os.environ.copy()
    amberhome = tleap.resolve().parent.parent
    environment["AMBERHOME"] = str(amberhome)
    environment["LD_LIBRARY_PATH"] = (
        str(amberhome / "lib") + ":" + environment.get("LD_LIBRARY_PATH", "")
    )
    with (output / log_name).open("w") as log:
        subprocess.run(
            [str(tleap.resolve()), "-f", input_name],
            cwd=output,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            text=True,
        )


def _charge(log: Path) -> float:
    values = [
        float(line.rsplit(":", 1)[1])
        for line in log.read_text().splitlines()
        if "Total unperturbed charge:" in line
    ]
    if not values:
        raise ValueError(f"No charge record in {log}")
    return values[-1]


def _volume(coordinates: Path) -> float:
    values = [
        float(value) for value in coordinates.read_text().splitlines()[-1].split()
    ]
    if len(values) != 6:
        raise ValueError("Expected six periodic box values")
    a, b, c, alpha, beta, gamma = values
    ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
    factor = math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)
    return a * b * c * factor


def _ion_counts(charge: float, volume: float) -> tuple[int, int, int]:
    pairs = round(SALT_MOLAR * 6.02214076e23 * 1.0e-27 * volume)
    integer_charge = round(charge)
    if abs(charge - integer_charge) > 0.01:
        raise ValueError(f"Nonintegral dry-system charge: {charge}")
    return (
        pairs + max(0, -integer_charge),
        pairs + max(0, integer_charge),
        pairs,
    )


def prepare_es_system(
    source_pdb: Path,
    ligand_library: Path,
    ligand_frcmod: Path,
    tleap: Path,
    output: Path,
    source_label: str = "Model 23",
    source_snapshot_name: str = "source_model_23.pdb",
) -> dict[str, Any]:
    """Build and audit a fresh Nyl12 ff19SB/GAFF1/TIP3P system."""
    required = (source_pdb, ligand_library, ligand_frcmod, tleap)
    if any(not path.is_file() for path in required):
        raise FileNotFoundError("A required Model 23 preparation input is missing")
    if not source_label.strip():
        raise ValueError("source_label must be nonempty")
    snapshot = Path(source_snapshot_name)
    if snapshot.name != source_snapshot_name or snapshot.suffix.lower() != ".pdb":
        raise ValueError("source_snapshot_name must be a PDB basename")
    output.mkdir(parents=True, exist_ok=False)
    params = output / "parameters"
    params.mkdir()
    shutil.copy2(ligand_library, params / "lig.lib")
    shutil.copy2(ligand_frcmod, params / "lig.frcmod")
    shutil.copy2(source_pdb, output / snapshot)
    write_leap_pdb(
        source_pdb,
        output / "input_complex_min_leap.pdb",
        off_atom_names(ligand_library),
    )
    (output / "leap_count.in").write_text(leap_input(params.resolve(), None, None))
    _run_tleap(tleap, "leap_count.in", output, "leap_count.log")
    charge = _charge(output / "leap_count.log")
    volume = _volume(output / "solvated_noions.inpcrd")
    sodium, chloride, pairs = _ion_counts(charge, volume)
    (output / "leap_solvate.in").write_text(
        leap_input(params.resolve(), sodium, chloride)
    )
    _run_tleap(tleap, "leap_solvate.in", output, "leap_solvate.log")
    topology = AmberParm(
        str(output / "md_solvated.prmtop"),
        xyz=str(output / "md_solvated.inpcrd"),
    )
    for residue in BETA_N_TERMINI:
        names = {atom.name for atom in topology.residues[residue - 1].atoms}
        if not {"N", "H1", "H2", "OG1", "HG1"} <= names or "H3" in names:
            raise ValueError(f"Incorrect neutral Thr inventory at residue {residue}")
    ligand = {atom.name for atom in topology.residues[LIGAND_RESIDUE - 1].atoms}
    if ligand != off_atom_names(params / "lig.lib"):
        raise ValueError("Solvated topology changed the PA6 atom inventory")
    sources = [*required, Path(__file__).resolve()]
    report = {
        "source_label": source_label,
        "source_snapshot": str((output / snapshot).resolve()),
        "force_field": "ff19SB + GAFF1 + TIP3P",
        "charge_method": "historical PA6 library; not RESP",
        "atom_count": len(topology.atoms),
        "residue_count": len(topology.residues),
        "dry_charge": charge,
        "box_volume_A3": volume,
        "salt_pairs_150_mM": pairs,
        "sodium_ions": sodium,
        "chloride_ions": chloride,
        "total_partial_charge": float(sum(atom.charge for atom in topology.atoms)),
        "sources": [
            {"path": str(path.resolve()), "sha256": sha256(path)} for path in sources
        ],
    }
    (output / "preparation.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
