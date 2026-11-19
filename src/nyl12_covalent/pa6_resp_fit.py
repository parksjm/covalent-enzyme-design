"""Run and audit a shared four-conformer RESP fit for PA6 L4."""

from __future__ import annotations

import hashlib
import json
import logging
import re
import subprocess
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from nyl12_covalent.chemistry import load_pa6
from nyl12_covalent.gaussian_esp import read_espgen
from nyl12_covalent.ti1_reference_graph import validate_pa6_zwitterion

LOGGER = logging.getLogger(__name__)
ATOM_COUNT = 79
CONFORMER_COUNT = 4


def fit_pa6_multiconformer_resp(
    source_mol2: Path,
    esp_paths: Sequence[Path],
    amber_bin: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Fit shared two-stage RESP charges and create candidate GAFF2 files.

    The four ESP datasets must use the canonical 79-atom PA6 L4 order. The
    generated MOL2 and frcmod are candidates pending bonded-parameter and MD
    validation; this function does not certify a production force field.

    Args:
        source_mol2: Canonical explicit-hydrogen PA6 L4 zwitterion MOL2.
        esp_paths: Four accepted Amber ``espgen`` datasets in declared order.
        amber_bin: AmberTools directory containing antechamber, respgen, resp,
            and parmchk2.
        output_dir: New, nonexisting directory for all fitting artifacts.

    Returns:
        Serializable fitting manifest and scientific checks.

    Raises:
        FileExistsError: If ``output_dir`` already exists.
        FileNotFoundError: If an input or required AmberTools executable is absent.
        RuntimeError: If a tool fails or a fitting invariant is violated.
        ValueError: If the molecular graph or ESP input is malformed.
    """
    if output_dir.exists():
        raise FileExistsError(f"RESP output directory already exists: {output_dir}")
    if len(esp_paths) != CONFORMER_COUNT or len(set(esp_paths)) != CONFORMER_COUNT:
        raise ValueError("Exactly four distinct ESP datasets are required")
    molecule = load_pa6(source_mol2).mol
    validate_pa6_zwitterion(molecule)
    names = [atom.GetProp("name") for atom in molecule.GetAtoms()]
    if len(names) != ATOM_COUNT or len(set(names)) != ATOM_COUNT:
        raise ValueError("Canonical PA6 atom names must be unique and complete")

    tools = {
        name: amber_bin / name
        for name in ("antechamber", "respgen", "resp", "parmchk2")
    }
    missing = [
        str(path)
        for path in [source_mol2, *esp_paths, *tools.values()]
        if not path.is_file()
    ]
    if missing:
        raise FileNotFoundError(f"Missing RESP prerequisites: {missing}")
    for esp_path in esp_paths:
        read_espgen(esp_path, ATOM_COUNT)

    output_dir.mkdir(parents=True)
    ac_path = output_dir / "pa6_l4_zwitterion.ac"
    stage1_input = output_dir / "resp_stage1.in"
    stage2_input = output_dir / "resp_stage2.in"
    combined_esp = output_dir / "four_conformers.esp"
    combined_esp.write_bytes(b"".join(path.read_bytes() for path in esp_paths))

    _run(
        "antechamber_ac",
        [
            tools["antechamber"],
            "-i",
            source_mol2,
            "-fi",
            "mol2",
            "-o",
            ac_path,
            "-fo",
            "ac",
            "-at",
            "gaff2",
            "-nc",
            "0",
            "-pf",
            "y",
        ],
        output_dir,
    )
    for stage, destination in (("resp1", stage1_input), ("resp2", stage2_input)):
        _run(
            f"respgen_{stage}",
            [
                tools["respgen"],
                "-i",
                ac_path,
                "-o",
                destination,
                "-f",
                stage,
                "-e",
                "2",
                "-n",
                str(CONFORMER_COUNT),
            ],
            output_dir,
        )
        _validate_resp_input(destination)

    stage1_qout = output_dir / "resp_stage1.qout"
    stage2_qout = output_dir / "resp_stage2.qout"
    _run_resp(tools["resp"], output_dir, 1, stage1_input, combined_esp, stage1_qout)
    _run_resp(
        tools["resp"],
        output_dir,
        2,
        stage2_input,
        combined_esp,
        stage2_qout,
        stage1_qout,
    )
    charges = parse_shared_resp_charges(stage2_qout, ATOM_COUNT, CONFORMER_COUNT)
    charge_path = output_dir / "pa6_l4_resp_charges.txt"
    charge_path.write_text("\n".join(f"{charge:.8f}" for charge in charges) + "\n")

    typed_mol2 = output_dir / "pa6_l4_zwitterion_resp_gaff2.mol2"
    frcmod = output_dir / "pa6_l4_zwitterion_resp_gaff2.frcmod"
    _run(
        "antechamber_apply_resp",
        [
            tools["antechamber"],
            "-i",
            ac_path,
            "-fi",
            "ac",
            "-o",
            typed_mol2,
            "-fo",
            "mol2",
            "-at",
            "gaff2",
            "-c",
            "rc",
            "-cf",
            charge_path,
            "-nc",
            "0",
            "-pf",
            "y",
        ],
        output_dir,
    )
    _run(
        "parmchk2_gaff2",
        [tools["parmchk2"], "-i", typed_mol2, "-f", "mol2", "-o", frcmod, "-s", "2"],
        output_dir,
    )
    stage2_text = (output_dir / "resp_stage2.out").read_text()
    rrms = _extract_rrms(stage2_text)
    manifest = {
        "schema_version": 1,
        "scope": "PA6 L4 zwitterion shared four-conformer two-stage RESP fit",
        "conformer_order": [path.parent.name for path in esp_paths],
        "atom_count": ATOM_COUNT,
        "atom_names_in_order": names,
        "net_charge": float(np.sum(charges)),
        "stage2_esp_relative_rms": rrms,
        "checks": {
            "four_valid_esp_datasets": True,
            "cross_conformer_charges_identical": True,
            "net_charge_zero": bool(abs(np.sum(charges)) <= 1.0e-5),
            "gaff2_candidate_written": typed_mol2.is_file(),
            "parmchk2_candidate_written": frcmod.is_file(),
        },
        "forcefield_md_certified": False,
        "next_required_review": [
            "inspect parmchk2 missing/analogous parameters",
            "validate nylon bonded terms and torsion scans",
            "isolated and solvated minimization/MD validation",
        ],
        "input_sha256": {
            str(path): _sha256(path) for path in [source_mol2, *esp_paths]
        },
        "tool_sha256": {name: _sha256(path) for name, path in tools.items()},
    }
    if not all(manifest["checks"].values()):
        raise RuntimeError(f"RESP/GAFF2 candidate checks failed: {manifest['checks']}")
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Prepared non-certified PA6 RESP/GAFF2 candidate at %s", output_dir)
    return manifest


def parse_shared_resp_charges(
    qout_path: Path, atom_count: int, conformer_count: int
) -> np.ndarray:
    """Return one charge vector after verifying all conformer copies agree.

    Args:
        qout_path: RESP qout containing all fitted molecular copies.
        atom_count: Atoms in each conformer.
        conformer_count: Number of conformer copies.

    Returns:
        One shared atomic-charge vector.

    Raises:
        ValueError: If count, finiteness, equality, or net charge is invalid.
    """
    values = np.asarray(
        [float(item) for item in qout_path.read_text().split()], dtype=float
    )
    expected = atom_count * conformer_count
    if values.size != expected or not np.isfinite(values).all():
        raise ValueError(
            f"Expected {expected} finite RESP charges; found {values.size}"
        )
    matrix = values.reshape(conformer_count, atom_count)
    if not np.allclose(matrix, matrix[0], atol=5.0e-6, rtol=0.0):
        raise ValueError("RESP conformer charge vectors are not identical")
    if abs(float(np.sum(matrix[0]))) > 1.0e-5:
        raise ValueError("Shared RESP charge vector is not net neutral")
    return matrix[0]


def _validate_resp_input(path: Path) -> None:
    text = path.read_text()
    if not re.search(r"\bnmol\s*=\s*4\s*,", text) or text.count("    0   79") != 4:
        raise RuntimeError(f"RESP input lacks four neutral 79-atom blocks: {path}")
    for conformer in range(2, CONFORMER_COUNT + 1):
        if f"    1    1    {conformer}    1" not in text:
            raise RuntimeError(f"RESP input lacks cross-conformer constraints: {path}")


def _run_resp(
    executable: Path,
    output_dir: Path,
    stage: int,
    input_path: Path,
    esp_path: Path,
    qout_path: Path,
    qin_path: Path | None = None,
) -> None:
    command: list[Path | str] = [
        executable,
        "-O",
        "-i",
        input_path.name,
        "-o",
        f"resp_stage{stage}.out",
        "-p",
        f"resp_stage{stage}.pch",
        "-t",
        qout_path.name,
        "-e",
        esp_path.name,
    ]
    if qin_path is not None:
        command.extend(["-q", qin_path.name])
    _run(f"resp_stage{stage}", command, output_dir, cwd=output_dir)


def _run(
    tag: str,
    command: Sequence[Path | str],
    output_dir: Path,
    *,
    cwd: Path | None = None,
) -> None:
    argv = [str(item) for item in command]
    (output_dir / f"{tag}.command.json").write_text(
        json.dumps({"argv": argv}, indent=2) + "\n"
    )
    with (
        (output_dir / f"{tag}.stdout.log").open("x") as stdout,
        (output_dir / f"{tag}.stderr.log").open("x") as stderr,
    ):
        try:
            subprocess.run(argv, check=True, stdout=stdout, stderr=stderr, cwd=cwd)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"{tag} failed with status {exc.returncode}; see {tag}.stderr.log"
            ) from exc


def _extract_rrms(text: str) -> float:
    matches = re.findall(r"ESP relative RMS .*?\s([-+\d.]+)\s*$", text, re.MULTILINE)
    if len(matches) != 1:
        raise RuntimeError("Could not identify one RESP relative RMS statistic")
    return float(matches[0])


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
