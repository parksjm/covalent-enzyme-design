"""Construct and validate Nyl12 systems containing released modular PA6 L4."""

from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from nyl12_covalent.modular_pa6_md_protocol import (
    cpptraj_analysis,
    mdin_files,
    replica_runner,
)
from nyl12_covalent.modular_pa6_topology import (
    sha256,
    validate_release,
    validate_topology,
)


BETA_NTERM_RESIDUES = (929, 1025, 1121, 1217)
CHAIN_END_RESIDUES = (232, 464, 696, 928, 1024, 1120, 1216, 1312, 1316)
NTERMINUS_RESIDUES = (1, 233, 465, 697, 929, 1025, 1121, 1217)
AVOGADRO = 6.02214076e23


@dataclass(frozen=True)
class SystemSpec:
    """One accepted Nyl12/PA6 starting complex."""

    label: str
    pdb_path: Path
    conversion_manifest: Path


def prepare_campaign(
    systems: list[SystemSpec],
    release_dir: Path,
    output_dir: Path,
    tleap: Path,
    pmemd_cuda: Path,
    *,
    replicas: int = 3,
) -> dict[str, Any]:
    """Build, validate, and write a portable PA6 L4 MD campaign.

    Args:
        systems: Accepted modular PDB inputs and conversion manifests.
        release_dir: Pinned PA6 L4 release directory.
        output_dir: New immutable campaign directory.
        tleap: AmberTools tleap executable.
        pmemd_cuda: Pinned Amber pmemd.cuda executable path on run hosts.
        replicas: Number of independent 100 ns replicas per system.

    Returns:
        Complete campaign manifest.
    """
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    if not systems or replicas < 1:
        raise ValueError("At least one system and one replica are required")
    release = validate_release(release_dir)
    output_dir.mkdir(parents=True)
    rows = []
    replica_paths: list[str] = []
    for system_index, spec in enumerate(systems, start=1):
        rows.append(
            _prepare_system(
                spec, release_dir, output_dir / "systems" / spec.label, tleap
            )
        )
        for replica in range(1, replicas + 1):
            relative = Path("systems") / spec.label / f"replica_{replica:02d}"
            _write_replica(
                output_dir / relative,
                810_000 + system_index * 100 + replica * 17,
                pmemd_cuda,
            )
            replica_paths.append(str(relative))
    _write_queues(output_dir, replica_paths)
    manifest = {
        "schema_version": 1,
        "status": "prepared_and_topology_validated_not_yet_simulated",
        "release_id": release["release_id"],
        "release_manifest_sha256": sha256(release_dir / "manifest.json"),
        "force_field": "ff19SB + GAFF2 modular PA6 RESP/improper release + TIP3P",
        "pb_radii": "mbondi3",
        "salt_molar": 0.150,
        "water_buffer_angstrom": 12.0,
        "production_ns_per_replica": 100,
        "replicas_per_system": replicas,
        "systems": rows,
        "replica_directories": replica_paths,
        "pmemd_cuda_path": str(pmemd_cuda),
        "tleap_sha256": sha256(tleap),
        "input_sha256": {
            str(spec.pdb_path.resolve()): sha256(spec.pdb_path) for spec in systems
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def _prepare_system(
    spec: SystemSpec, release_dir: Path, system_dir: Path, tleap: Path
) -> dict[str, Any]:
    conversion = json.loads(spec.conversion_manifest.read_text(encoding="utf-8"))
    if conversion.get("status") != "accepted_legacy_whole_l4_to_modular_pdb_conversion":
        raise ValueError(f"Unaccepted conversion manifest: {spec.conversion_manifest}")
    if conversion.get("output_pdb_sha256") != sha256(spec.pdb_path):
        raise ValueError(f"Converted PDB hash mismatch: {spec.pdb_path}")
    system_dir.mkdir(parents=True)
    _write_leap_pdb(spec.pdb_path, system_dir / "input_complex_leap.pdb")
    for path in release_dir.iterdir():
        if path.is_file() and path.name not in {"README.md", "UPSTREAM.json"}:
            shutil.copy2(path, system_dir / path.name)
    _write_leap_input(system_dir / "leap_count.in", final_ions=None)
    count_log = _run_tleap(tleap, system_dir, "leap_count.in", "leap_count.log")
    _require_clean_leap(count_log)
    dry_report = validate_topology(
        system_dir / "md_dry.prmtop", release_dir, solvated=False
    )
    volume = _box_volume(system_dir / "solvated_noions.inpcrd")
    sodium, chloride, pairs = _ion_counts(dry_report["total_charge"], volume)
    _write_leap_input(system_dir / "leap_final.in", final_ions=(sodium, chloride))
    final_log = _run_tleap(tleap, system_dir, "leap_final.in", "leap_final.log")
    _require_clean_leap(final_log)
    final_report = validate_topology(
        system_dir / "md_solvated.prmtop", release_dir, solvated=True
    )
    report = {
        "label": spec.label,
        "source_pdb_sha256": sha256(spec.pdb_path),
        "dry_topology": dry_report,
        "solvated_topology": final_report,
        "box_volume_angstrom3": volume,
        "salt_pairs_150mM": pairs,
        "sodium_ions": sodium,
        "chloride_ions": chloride,
        "checks": {
            "accepted_conversion_hash": True,
            "leap_logs_have_zero_errors": True,
            "modular_sequence_and_atoms": True,
            "modular_types_and_charges_match_release": True,
            "three_inter_residue_pa6_bonds": True,
            "mbondi3_radii": True,
            "solvated_system_neutral": abs(final_report["total_charge"]) < 1.0e-4,
        },
    }
    if not all(report["checks"].values()):
        raise ValueError(f"Topology acceptance failed for {spec.label}")
    (system_dir / "acceptance.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


def _write_leap_pdb(source: Path, destination: Path) -> None:
    output: list[str] = []
    for line in source.read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = int(line[22:26])
        atom = line[12:16].strip()
        if residue in NTERMINUS_RESIDUES and atom == "H":
            line = f"{line[:12]} H1 {line[16:]}"
        output.append(line)
        if residue in CHAIN_END_RESIDUES and atom == "OXT":
            output.append("TER")
    output.append("END")
    destination.write_text("\n".join(output) + "\n", encoding="ascii")


def _leap_header() -> str:
    edits = "\n".join(
        f"remove complex complex.{residue}.H3\n"
        f"set complex.{residue}.N charge -0.4844\n"
        f"set complex.{residue}.H1 charge 0.1224\n"
        f"set complex.{residue}.H2 charge 0.1224"
        for residue in BETA_NTERM_RESIDUES
    )
    return f"""source leaprc.protein.ff19SB
source leaprc.water.tip3p
source leaprc.pa6_shared_p6i
complex = loadpdb input_complex_leap.pdb
{edits}
check complex
charge complex
"""


def _write_leap_input(path: Path, final_ions: tuple[int, int] | None) -> None:
    text = _leap_header()
    text += "saveamberparm complex md_dry.prmtop md_dry.inpcrd\n"
    text += "savepdb complex md_dry.pdb\nsolvateOct complex TIP3PBOX 12.0\n"
    if final_ions is None:
        text += "saveamberparm complex solvated_noions.prmtop solvated_noions.inpcrd\n"
    else:
        sodium, chloride = final_ions
        text += f"addIonsRand complex Na+ {sodium} Cl- {chloride} 6.0\n"
        text += "check complex\ncharge complex\n"
        text += "saveamberparm complex md_solvated.prmtop md_solvated.inpcrd\n"
        text += "savepdb complex md_solvated.pdb\n"
    path.write_text(text + "quit\n", encoding="ascii")


def _run_tleap(tleap: Path, cwd: Path, input_name: str, log_name: str) -> str:
    result = subprocess.run(
        [str(tleap), "-f", input_name],
        cwd=cwd,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (cwd / log_name).write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"tleap failed ({result.returncode}); see {cwd / log_name}")
    return result.stdout


def _require_clean_leap(log: str) -> None:
    errors = [int(value) for value in re.findall(r"Errors\s*=\s*(\d+)", log)]
    if "FATAL:" in log or any(value != 0 for value in errors):
        raise ValueError("LEaP reported a fatal error or nonzero error count")


def _box_volume(path: Path) -> float:
    values = [float(value) for value in path.read_text().splitlines()[-1].split()]
    if len(values) != 6:
        raise ValueError(f"Expected six periodic box values in {path}")
    a_len, b_len, c_len, alpha, beta, gamma = values
    ca, cb, cg = (math.cos(math.radians(value)) for value in (alpha, beta, gamma))
    return (
        a_len
        * b_len
        * c_len
        * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)
    )


def _ion_counts(charge: float, volume_angstrom3: float) -> tuple[int, int, int]:
    pairs = round(0.150 * AVOGADRO * 1.0e-27 * volume_angstrom3)
    integral_charge = round(charge)
    if abs(charge - integral_charge) > 1.0e-4:
        raise ValueError(f"Dry system charge is not integral: {charge:.8f}")
    return pairs + max(0, -integral_charge), pairs + max(0, integral_charge), pairs


def _write_replica(path: Path, seed: int, pmemd_cuda: Path) -> None:
    path.mkdir(parents=True)
    for name, text in mdin_files(seed).items():
        (path / name).write_text(text, encoding="ascii")
    (path / "run_replica.sh").write_text(replica_runner(pmemd_cuda), encoding="ascii")
    (path / "analyze_replica.cpptraj").write_text(cpptraj_analysis(), encoding="ascii")


def _write_queues(root: Path, replicas: list[str]) -> None:
    queues = root / "queues"
    queues.mkdir()
    assignments = {
        "violet_gpu0": [],
        "violet_gpu1": [],
        "azalea_gpu0": [],
        "azalea_gpu1": [],
    }
    for index, replica in enumerate(replicas):
        list(assignments.values())[index % 4].append(replica)
    for name, entries in assignments.items():
        (queues / f"{name}.txt").write_text("\n".join(entries) + "\n", encoding="utf-8")
    (root / "host_assignment.json").write_text(
        json.dumps(assignments, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
