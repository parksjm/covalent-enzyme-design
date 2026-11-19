"""Verify staged ASH272/ASH274 inputs and run two independent 100 ns replicas."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import logging
from math import cos, pi, sqrt
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone

import parmed as pmd


def digest(path: Path) -> str:
    """Return SHA256 of a retained file."""
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def record(path: Path, data: object) -> None:
    """Write a machine-readable process or acceptance record."""
    path.write_text(json.dumps(data, indent=2) + "\n")


def execute(
    command: list[str], folder: Path, name: str, gpu: int | None = None
) -> None:
    """Execute a logged process with an explicit optional GPU assignment."""
    data = {
        "command": command,
        "cwd": str(folder),
        "gpu": gpu,
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    record(folder / f"{name}.process.json", data)
    environment = os.environ.copy()
    if gpu is not None:
        environment["CUDA_VISIBLE_DEVICES"] = str(gpu)
    with (folder / f"{name}.engine.log").open("w") as handle:
        process = subprocess.run(
            command,
            cwd=folder,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    data.update(
        returncode=process.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    record(folder / f"{name}.process.json", data)
    if process.returncode:
        raise RuntimeError(f"{name} exited {process.returncode}")


def verify(
    original: pmd.Structure, corrected: pmd.Structure, target_residue: int = 272
) -> dict[str, object]:
    """Require four added protons, +4 charge, neutral Thr and unchanged substrate."""
    if target_residue not in {272, 274}:
        raise ValueError("Target residue must be Asp272 or Asp274")
    old_charge = sum(atom.charge for atom in original.atoms)
    charge = sum(atom.charge for atom in corrected.atoms)
    if (
        abs(charge - old_charge - 4) > 1e-5
        or len(corrected.atoms) != len(original.atoms) + 4
    ):
        raise ValueError("Expected four added atoms and exactly +4 e dry-charge change")
    sites = tuple(
        site + (2 if target_residue == 274 else 0) for site in (968, 1064, 1160, 1256)
    )
    for site in sites:
        residue = corrected.residues[site - 1]
        atoms = {atom.name: atom for atom in residue.atoms}
        if residue.name != "ASH" or "HD2" not in atoms:
            raise ValueError(f"ASH/HD2 absent at {site}")
        if atoms["OD2"] not in atoms["HD2"].bond_partners:
            raise ValueError(f"HD2 not bonded to OD2 at {site}")
        partner = site + (2 if target_residue == 272 else -2)
        if corrected.residues[partner - 1].name != "ASP":
            raise ValueError(f"Partner Asp changed at {partner}")
    for site in (929, 1025, 1121, 1217):
        atoms = {atom.name: atom for atom in corrected.residues[site - 1].atoms}
        if "H3" in atoms:
            raise ValueError(f"Thr233 has H3 at {site}")
        for name, expected in (("N", -0.4844), ("H1", 0.1224), ("H2", 0.1224)):
            if name not in atoms or abs(atoms[name].charge - expected) > 1e-6:
                raise ValueError(f"Neutral Thr233 charge mismatch at {site}:{name}")

    def substrate(structure: pmd.Structure) -> list[tuple[object, ...]]:
        return [
            (
                residue.idx,
                residue.name,
                atom.name,
                atom.type,
                atom.charge,
                atom.mass,
                atom.epsilon,
                atom.rmin,
            )
            for residue in structure.residues
            if residue.name.startswith("P6")
            for atom in residue.atoms
        ]

    before = substrate(original)
    if not before or before != substrate(corrected):
        raise ValueError(
            "Substrate mapping, RESP charges or nonbonded parameters changed"
        )
    return {
        "original_dry_charge_e": old_charge,
        "dry_charge_e": charge,
        "charge_change_e": charge - old_charge,
        "four_HD2": True,
        "protonated_residue": target_residue,
        "neutral_Thr": True,
        "substrate_parameters_unchanged": True,
    }


def prepare_system(args: argparse.Namespace) -> tuple[Path, Path]:
    """Finish LEaP with recounted ions and perform the matched MM minimization."""
    folder = args.system_dir
    manifest = json.loads((folder / "preparation_manifest.json").read_text())
    for filename, expected in manifest["prepared_sha256"].items():
        if digest(folder / filename) != expected:
            raise ValueError(f"Staged input changed: {filename}")
    for filename, expected in manifest["source_sha256"].items():
        if digest(args.source_system_dir / filename) != expected:
            raise ValueError(f"Original release changed: {filename}")
    if args.reuse_count:
        if args.count_log is None:
            raise ValueError("--reuse-count requires its retained --count-log")
        log = args.count_log
    else:
        if (folder / "count.prmtop").exists() or (folder / "count.inpcrd").exists():
            raise FileExistsError("Count outputs exist; use explicit --reuse-count")
        execute([str(args.tleap), "-f", "leap_count.in"], folder, "leap_count")
        log = folder / "leap_count.engine.log"
    if re.search(r"FATAL|Exiting LEaP: Errors =\s*[1-9]", log.read_text()):
        raise ValueError("LEaP count log reports errors")
    prefix = (folder / "leap_count.in").read_text().split("solvateOct", 1)[0]
    saves = re.findall(r"saveamberparm\s+complex\s+(\S+)\s+(\S+)", prefix, re.I)
    if len(saves) != 1:
        raise ValueError("Expected one dry-topology save")
    dry_name = saves[0][0]
    old = pmd.load_file(str(args.source_system_dir / dry_name))
    dry = pmd.load_file(str(folder / dry_name))
    target_residue = int(manifest["protonated_residue"])
    acceptance = verify(old, dry, target_residue)
    count = pmd.load_file(
        str(folder / "count.prmtop"), xyz=str(folder / "count.inpcrd")
    )
    a, b, c, alpha, beta, gamma = map(float, count.box)
    ca, cb, cg = (cos(angle * pi / 180) for angle in (alpha, beta, gamma))
    volume = a * b * c * sqrt(1 + 2 * ca * cb * cg - ca * ca - cb * cb - cg * cg)
    pairs = round(0.150 * volume * 6.02214076e-4)
    charge = float(acceptance["dry_charge_e"])
    sodium, chloride = pairs + max(0, -round(charge)), pairs + max(0, round(charge))
    residual = charge - round(charge)
    final = folder / f"leap_ash{target_residue}_final.in"
    if final.exists():
        raise FileExistsError(final)
    final.write_text(
        prefix
        + "solvateOct complex TIP3PBOX 12.0\n"
        + f"addIonsRand complex Na+ {sodium} Cl- {chloride} 6.0\n"
        + f"saveAmberParm complex ash{target_residue}_solvated.prmtop "
        f"ash{target_residue}_solvated.inpcrd\n"
        + f"savePdb complex ash{target_residue}_solvated.pdb\nquit\n"
    )
    execute([str(args.tleap), "-f", str(final)], folder, "leap_final")
    if re.search(
        r"FATAL|Exiting LEaP: Errors =\s*[1-9]",
        (folder / "leap_final.engine.log").read_text(),
    ):
        raise ValueError("Final LEaP log reports errors")
    topology = folder / f"ash{target_residue}_solvated.prmtop"
    solvated = pmd.load_file(str(topology))
    if abs(sum(atom.charge for atom in solvated.atoms) - residual) > 1e-5:
        raise ValueError("Final charge differs from inherited release residual")
    acceptance.update(
        box_volume_angstrom3=volume,
        salt_molar=0.150,
        salt_pairs=pairs,
        sodium=sodium,
        chloride=chloride,
        inherited_charge_residual_e=residual,
        count_log_sha256=digest(log),
        count_topology_sha256=digest(folder / "count.prmtop"),
        count_coordinates_sha256=digest(folder / "count.inpcrd"),
    )
    record(folder / "preparation_acceptance.json", acceptance)
    mdin = folder / "min_unrestrained.in"
    mdin.write_text(
        f"Matched ASH{target_residue} MM minimization\n&cntrl\n"
        " imin=1,maxcyc=10000,ncyc=1000,drms=.001,ntb=1,cut=10,\n"
        " ntr=0,nmropt=0,ntc=1,ntf=1,ntpr=100,ntwx=10,ntwr=1000,\n"
        " ioutfm=1,ntxo=2,\n/\n"
    )
    restart = folder / "min_unrestrained.rst7"
    execute(
        [
            str(args.engine),
            "-O",
            "-i",
            str(mdin),
            "-o",
            str(folder / "min_unrestrained.out"),
            "-p",
            str(topology),
            "-c",
            str(folder / f"ash{target_residue}_solvated.inpcrd"),
            "-r",
            str(restart),
            "-x",
            str(folder / "min_unrestrained.nc"),
        ],
        folder,
        "min_unrestrained",
        0,
    )
    if re.search(
        r"\b(?:NaN|Infinity)\b|FATAL|ERROR",
        (folder / "min_unrestrained.out").read_text(),
        re.I,
    ):
        raise ValueError("Minimization reports numerical/engine error")
    pmd.load_file(str(topology), xyz=str(restart))
    record(
        folder / "minimization_acceptance.json",
        {
            "process_complete": True,
            "restart_readable": True,
            "force_convergence": "not asserted; MD preparation only",
            "restart_sha256": digest(restart),
            "topology_sha256": digest(topology),
        },
    )
    return topology, restart


def main() -> None:
    """Run preparation followed by independent GPU zero/one 10+90 ns replicas."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("system-dir", "source-system-dir", "engine", "tleap", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--reuse-count", action="store_true")
    parser.add_argument("--count-log", type=Path)
    args = parser.parse_args()
    for name in ("system_dir", "source_system_dir", "engine", "tleap", "python"):
        setattr(args, name, getattr(args, name).resolve(strict=True))
    if args.count_log is not None:
        args.count_log = args.count_log.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    logging.basicConfig(level=logging.INFO)
    root = Path(__file__).resolve().parents[1]
    snapshot = args.output / "source_snapshot"
    files = [
        "scripts/run_ash272_overnight.py",
        "scripts/prepare_flipped_model30_ash272.py",
        "scripts/run_model30_md_replica.py",
        "scripts/continue_model30_md.py",
        "src/nyl12_covalent/pilot_md_protocol.py",
    ]
    for relative in files:
        target = snapshot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / relative, target)
    record(
        args.output / "source_hashes.json",
        {name: digest(snapshot / name) for name in files},
    )
    record(args.output / "status.json", {"exit_status": None})
    try:
        topology, coordinates = prepare_system(args)

        def replica(gpu: int) -> None:
            target_residue = json.loads(
                (args.system_dir / "preparation_manifest.json").read_text()
            )["protonated_residue"]
            label = f"ash{target_residue}_replica_{gpu + 1}"
            folder = args.output / label
            execute(
                [
                    str(args.python),
                    str(snapshot / "scripts/run_model30_md_replica.py"),
                    "--topology",
                    str(topology),
                    "--coordinates",
                    str(coordinates),
                    "--engine",
                    str(args.engine),
                    "--output",
                    str(folder),
                    "--seed",
                    str(732101 + gpu * 101),
                    "--label",
                    label,
                ],
                args.output,
                f"{label}_10ns",
                gpu,
            )
            pilot_job = folder / "pilot_completion"
            pilot_job.mkdir()
            record(pilot_job / "status.json", {"exit_status": 0})
            execute(
                [
                    str(args.python),
                    str(snapshot / "scripts/continue_model30_md.py"),
                    "--replica-dir",
                    str(folder),
                    "--source-job-dir",
                    str(pilot_job),
                    "--engine",
                    str(args.engine),
                    "--seed",
                    str(842101 + gpu * 101),
                ],
                args.output,
                f"{label}_90ns",
                gpu,
            )

        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(replica, gpu) for gpu in (0, 1)]
            for future in futures:
                future.result()
    except BaseException:
        record(args.output / "status.json", {"exit_status": 1})
        raise
    record(
        args.output / "status.json",
        {
            "exit_status": 0,
            "scientific_acceptance": "Sampling and chemical retention require analysis",
        },
    )


if __name__ == "__main__":
    main()
