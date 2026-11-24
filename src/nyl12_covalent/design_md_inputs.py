"""Derive portable MD inputs from an accepted four-ligand relaxation."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import parmed

from .design_model_relaxation import _final_results
from .design_relaxation_topology import validate_four_ligand_topology

LOGGER = logging.getLogger(__name__)


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _inventory(root: Path) -> dict[str, str]:
    return {
        str(p.relative_to(root)): _digest(p)
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def _write(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def _artifact(path: Path, root: Path) -> dict[str, str]:
    return {"path": str(path.relative_to(root)), "sha256": _digest(path)}


def _mapping(top: Any) -> list[dict[str, Any]]:
    residues = top.residues
    if len(residues) < 1328:
        raise ValueError("Topology lacks the required protein and four ligands")

    def normalize(name: str) -> str:
        return "HIS" if name in {"HID", "HIE", "HIP"} else name

    for start, length in ((0, 232), (928, 96)):
        copies = [
            [
                normalize(r.name)
                for r in residues[start + i * length : start + (i + 1) * length]
            ]
            for i in range(4)
        ]
        if any(copy != copies[0] for copy in copies[1:]):
            raise ValueError("Protein sequence differs across symmetry copies")
    rows = []
    for i, site in enumerate(("AE", "BF", "CG", "DH")):
        for sequence, base in (
            (113, 113),
            (156, 156),
            (186, 186),
            (233, 929),
            (272, 968),
            (274, 970),
        ):
            alpha = sequence < 233
            index = base + i * (232 if alpha else 96)
            name = residues[index - 1].name
            expected = {
                113: "PHE",
                156: "LYS",
                186: "ASN",
                233: "THR",
                272: "ASP",
                274: "ASP",
            }[sequence]
            if name != expected:
                raise ValueError(
                    f"Unexpected protected residue {site}/{sequence}: {name}"
                )
            rows.append(
                dict(
                    site=site,
                    chain=site[0 if alpha else 1],
                    sequence_residue=sequence,
                    structure_residue=sequence if alpha else sequence - 232,
                    amber_residue=index,
                    resname=name,
                )
            )
    return rows


def _verify_source(
    source: Path, release: Path
) -> tuple[dict[str, Any], Any, list[dict[str, Any]]]:
    acceptance = json.loads((source / "acceptance.json").read_text())
    process = json.loads((source / "process_acceptance.json").read_text())
    method = acceptance["method"]
    if (
        acceptance.get("status") != "minimization_complete"
        or acceptance.get("force_converged") is not True
        or acceptance.get("process_exit_code") != 0
        or process.get("status") != "completed"
        or not process.get("events")
        or any(e.get("returncode") != 0 for e in process["events"])
    ):
        raise ValueError("Native relaxation acceptance or process completion failed")
    if (
        method.get("protein_force_field") != "ff19SB"
        or method.get("water_model") != "TIP3P"
        or method.get("ligand_force_field") != "released modular PA6 RESP/GAFF2"
        or method.get("protonation")
        != "Neutral beta Thr; ASP272/ASP274 anionic at all sites."
        or method.get("restraints") != "none; ntr=0, nmropt=0, ntc=1, ntf=1"
        or method.get("drms_target_kcal_mol_angstrom") != 0.05
    ):
        raise ValueError("Native method or protonation differs from required protocol")
    original_release = source / "inputs" / "parameters"
    if (
        not release.is_dir()
        or _inventory(release) != _inventory(original_release)
        or not _inventory(release)
    ):
        raise ValueError("Parameter release differs from the preserved native release")
    top_path, rst = (
        source / "system" / "complex.prmtop",
        source / "minimization" / "min.rst7",
    )
    commands = [
        json.loads(line)
        for line in (source / "command_log.jsonl").read_text().splitlines()
    ]
    candidates = [
        c
        for c in commands
        if c.get("event") == "planned"
        and "-p" in c.get("argv", [])
        and "-r" in c.get("argv", [])
    ]
    if len(candidates) != 1:
        raise ValueError("Expected exactly one native minimization command")
    command = candidates[0]
    argv = command["argv"]
    paths = [Path(argv[argv.index(flag) + 1]) for flag in ("-p", "-r", "-i")]
    if paths[0].resolve() != top_path.resolve() or paths[1].resolve() != rst.resolve():
        raise ValueError("Native command topology or restart path mismatch")
    for path in (paths[0], paths[2]):
        if command["input_sha256"].get(str(path)) != _digest(path):
            raise ValueError("Native minimization input hash mismatch")
    mdin = paths[2].read_text().replace(" ", "").lower()
    if "ntr=0" not in mdin or "nmropt=0" not in mdin:
        raise ValueError("Native minimization input is not unrestrained")
    if not any(
        e.get("argv") == argv and e.get("returncode") == 0 for e in process["events"]
    ):
        raise ValueError(
            "Native minimization command has no matching successful process event"
        )
    _, rms = _final_results(
        (source / "minimization" / "min.out").read_text(errors="replace")
    )
    if rms > 0.05 or not math.isclose(
        rms, acceptance["final_rms_gradient"], rel_tol=1e-10, abs_tol=1e-12
    ):
        raise ValueError("Final RMS gradient does not satisfy the recorded acceptance")
    audit = validate_four_ligand_topology(top_path, release, solvated=True)
    if audit["release_id"] != acceptance["parameter_release_id"]:
        raise ValueError("Topology release identity differs from acceptance")
    top = parmed.load_file(str(top_path), xyz=str(rst))
    xyz, box = np.asarray(top.coordinates), np.asarray(top.box)
    if (
        xyz.shape != (len(top.atoms), 3)
        or box.shape != (6,)
        or not np.isfinite(xyz).all()
        or not np.isfinite(box).all()
        or np.any(box <= 0)
    ):
        raise ValueError("Restart coordinates or periodic box are invalid")
    export = acceptance["heavy_export"]
    stripped_top = source / "minimization" / "endpoint_complex.prmtop"
    stripped_rst = source / "minimization" / "endpoint_complex.rst7"
    if (
        _digest(stripped_top) != export["topology_sha256"]
        or _digest(stripped_rst) != export["restart_sha256"]
    ):
        raise ValueError("Stripped endpoint hash mismatch")
    stripped = parmed.load_file(str(stripped_top), xyz=str(stripped_rst))
    kept = [a for a in top.atoms if a.residue.name not in {"WAT", "Na+", "Cl-"}]

    def identities(atoms: Any) -> list[tuple[str, str, float]]:
        return [(a.name, a.residue.name, float(a.charge)) for a in atoms]

    if identities(kept) != identities(stripped.atoms) or not np.allclose(
        xyz[[a.idx for a in kept]], np.asarray(stripped.coordinates), rtol=0, atol=1e-6
    ):
        raise ValueError("Solvated restart differs from hashed stripped endpoint")
    return acceptance, top, _mapping(top)


def prepare_design_dataset(
    source_relaxation: Path, output: Path, dataset_id: str, release: Path
) -> dict[str, Any]:
    """Derive a self-contained portable dataset from a validated native relaxation.

    Args:
        source_relaxation: Native relaxation directory and its immutable receipts.
        output: New dataset root; existing paths are rejected.
        dataset_id: Explicit portable manifest identifier.
        release: Exact parameter release used by the relaxation.

    Returns:
        Portable datasets manifest, also written to ``datasets.json``.

    Raises:
        ValueError: If provenance, scientific acceptance or coordinates fail.
        FileExistsError: If the output already exists.
    """
    if output.exists():
        raise FileExistsError(output)
    if (
        not dataset_id
        or dataset_id in {".", ".."}
        or Path(dataset_id).name != dataset_id
    ):
        raise ValueError("dataset_id must be a single nonempty path component")
    source = source_relaxation.resolve(strict=True)
    acceptance, _, rows = _verify_source(source, release)
    output.mkdir(parents=True)
    folder = output / dataset_id
    folder.mkdir()
    evidence = folder / "evidence"
    evidence.mkdir()
    for relative in (
        "acceptance.json",
        "process_acceptance.json",
        "command_log.jsonl",
        "minimization/min.in",
        "minimization/min.out",
    ):
        shutil.copy2(source / relative, evidence / Path(relative).name)
    shutil.copytree(release, output / "parameters")
    top, rst = folder / "topology.prmtop", folder / "start.rst7"
    shutil.copy2(source / "system" / "complex.prmtop", top)
    shutil.copy2(source / "minimization" / "min.rst7", rst)
    provenance = {
        "source_relaxation": str(source),
        "upstream_receipts": _inventory(evidence),
        "parameter_files": _inventory(output / "parameters"),
        "derived": True,
    }
    prep, mini = (
        folder / "preparation_acceptance.json",
        folder / "minimization_acceptance.json",
    )
    _write(
        prep,
        {
            "derived": True,
            "source_provenance": provenance,
            "parameter_release_id": acceptance["parameter_release_id"],
        },
    )
    _write(
        mini,
        {
            "derived": True,
            "process_complete": True,
            "restart_readable": True,
            "coordinates_sha256": _digest(rst),
            "topology_sha256": _digest(top),
            "source_provenance": provenance,
        },
    )
    entry = {
        "microstate": "ASP272/ASP274",
        "terminal_protonation": "NH2",
        "restart_stage": "unrestrained_minimized",
        "protected_residues": rows,
    }
    for name, path in (
        ("topology", top),
        ("restart", rst),
        ("minimizedstart", rst),
        ("preparation_acceptance", prep),
        ("minimization_acceptance", mini),
    ):
        entry[name] = _artifact(path, output)
    manifest = {"schema_version": 1, "datasets": {dataset_id: entry}}
    _write(output / "datasets.json", manifest)
    LOGGER.info("Derived portable MD dataset %s from %s", dataset_id, source)
    return manifest
