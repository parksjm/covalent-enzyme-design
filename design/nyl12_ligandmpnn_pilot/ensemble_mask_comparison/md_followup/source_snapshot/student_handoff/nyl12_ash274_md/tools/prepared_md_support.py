"""Explicit input, symmetry and GPU gates for the portable prepared-MD runner."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
from typing import Any

from prepared_terminal_state import terminal_contract, validate_terminal

PROTECTED = {113: "PHE", 156: "LYS", 186: "ASN", 233: "THR"}
SITES = ("AE", "BF", "CG", "DH")


def digest(path: Path) -> str:
    """Return SHA-256 without loading full topology/trajectory files in memory."""
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def artifact(record: dict[str, str], base: Path) -> Path:
    """Resolve and verify a manifest-pinned input; never substitute another file."""
    path = Path(record["path"])
    path = (
        (base / path).resolve(strict=True)
        if not path.is_absolute()
        else path.resolve(strict=True)
    )
    if digest(path) != record["sha256"]:
        raise ValueError(f"Artifact SHA-256 mismatch: {path}")
    return path


def replicas(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Validate independent positive stage seeds, replica names and explicit GPU IDs."""
    rows = config["replicas"]
    if not rows:
        raise ValueError("At least one explicitly configured replica is required")
    names, seeds = set(), set()
    for row in rows:
        name, seed, gpu = row["name"], row["seed"], str(row["gpu"])
        if not name or any(
            c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-"
            for c in name
        ):
            raise ValueError("Replica name must be a portable alphanumeric identifier")
        if (
            name in names
            or type(seed) is not int
            or seed <= 0
            or seed + 50000 > 2147483647
        ):
            raise ValueError("Replica names and positive32bitseeds must be distinct")
        stage_seeds = {seed + 10000 * i for i in range(6)}
        if seeds & stage_seeds:
            raise ValueError(
                "Replica stage seeds overlap; choose independent seed ranges"
            )
        if not gpu.isdecimal():
            raise ValueError(
                "GPU must be an explicitly selected numeric NVIDIA device ID"
            )
        names.add(name)
        seeds.update(stage_seeds)
    return rows


def validate_protected(
    topology: Any,
    rows: list[dict[str, Any]],
    terminal_protonation: str = "NH2",
    expected_thr_charge: float = -0.001,
) -> None:
    """Require six protected identities across all four explicitly mapped sites.

    Args:
        topology: Loaded original ParmEd topology.
        rows: Site/chain/full-sequence/Amber-index/resname correspondence records.
        terminal_protonation: Legacy-neutral NH2 or explicitly declared NH3+.
        expected_thr_charge: Whole Thr charge; legacy-neutral recipe is -0.001e.

    Raises:
        ValueError: If mapping, catalytic identity or neutral-Thr inventory differs.
    """
    expected = {(site, residue) for site in SITES for residue in (*PROTECTED, 272, 274)}
    if (
        len(rows) != 24
        or {(r["site"], r["sequence_residue"]) for r in rows} != expected
    ):
        raise ValueError("Require exactly24mappedprotectedresidues acrossAE/BF/CG/DH")
    amber_ids = [r["amber_residue"] for r in rows]
    if len(set(amber_ids)) != 24:
        raise ValueError("Protected Amber residue indices must be distinct")
    for row in rows:
        site, sequence = row["site"], row["sequence_residue"]
        if row["chain"] != (site[0] if sequence < 233 else site[1]):
            raise ValueError(
                "Protected source chain differs from explicit symmetry-site mapping"
            )
        index = row["amber_residue"]
        if type(index) is not int or not 1 <= index <= len(topology.residues):
            raise ValueError("Protected Amber residue index is invalid")
        residue = topology.residues[index - 1]
        allowed = {PROTECTED[sequence]} if sequence in PROTECTED else {"ASP", "ASH"}
        if residue.name != row["resname"] or residue.name not in allowed:
            raise ValueError(
                f"Protected catalytic identity mismatch: {site}/{sequence}"
            )
        names = {a.name for a in residue.atoms}
        if sequence == 233:
            validate_terminal(residue, terminal_protonation, expected_thr_charge)
        if sequence in (272, 274) and (("HD2" in names) != (residue.name == "ASH")):
            raise ValueError("Declared Asp proton inventory differs from topology")


def load_dataset(
    manifest_path: Path, dataset_id: str
) -> tuple[dict[str, Any], Path, Path, dict[str, str]]:
    """Validate a prepared dataset, matching minimized receipt and all four active sites.

    Args:
        manifest_path: Dataset index with relative artifact paths and explicit hashes.
        dataset_id: Explicit dataset identifier, including microstate identity.

    Returns:
        Dataset record, topology, minimized restart and verified artifact identities.

    Raises:
        ValueError: If accepted preparation, coordinate matching or mapping fails.
    """
    import numpy as np
    import parmed

    document = json.loads(manifest_path.read_text())
    dataset = document["datasets"][dataset_id]
    if dataset["restart_stage"] != "unrestrained_minimized":
        raise ValueError(
            "Run the supplied authoritative MMminimization first; rawsolvatedstarts are not accepted"
        )
    if not dataset["microstate"]:
        raise ValueError("Explicit chemical microstate is required")
    files = {
        name: artifact(dataset[name], manifest_path.parent)
        for name in (
            "topology",
            "restart",
            "preparation_acceptance",
            "minimization_acceptance",
        )
    }
    receipt = json.loads(files["minimization_acceptance"].read_text())
    if (
        receipt.get("process_complete") is not True
        or receipt.get("restart_readable") is not True
    ):
        raise ValueError("Accepted unrestrained-minimization receipt is required")
    if receipt["coordinates_sha256"] != digest(files["restart"]) or receipt[
        "topology_sha256"
    ] != digest(files["topology"]):
        raise ValueError(
            "Minimization receipt does not certify these exact topology/start coordinates"
        )
    top = parmed.load_file(str(files["topology"]), xyz=str(files["restart"]))
    xyz = np.asarray(top.coordinates)
    box = np.asarray(top.box)
    if (
        xyz.shape != (len(top.atoms), 3)
        or not np.isfinite(xyz).all()
        or box.shape != (6,)
        or not np.isfinite(box).all()
        or np.any(box <= 0)
    ):
        raise ValueError(
            "Prepared topology/restart must have finite matching coordinates and periodic box"
        )
    states = {
        "ASP272/ASH274": ("ASP", "ASH"),
        "ASP272/ASP274": ("ASP", "ASP"),
        "ASH272/ASP274": ("ASH", "ASP"),
    }
    if dataset["microstate"] not in states:
        raise ValueError(
            "Microstate requires explicit reviewed preparation; supportedWTstates are ASP/ASH,ASP/ASP,ASH/ASP"
        )
    expected = dict(zip((272, 274), states[dataset["microstate"]], strict=True))
    if any(
        r["resname"] != expected[r["sequence_residue"]]
        for r in dataset["protected_residues"]
        if r["sequence_residue"] in expected
    ):
        raise ValueError("Dataset microstate differs from protected Asp/ASH mapping")
    terminal, charge, basis = terminal_contract(dataset)
    validate_protected(top, dataset["protected_residues"], terminal, charge)
    dataset = {
        **dataset,
        "validated_terminal_contract": {
            "state": terminal,
            "whole_thr_charge_e": charge,
            "basis": basis,
            "copy_count": 4,
            "four_thr_charge_sum_e": 4 * charge,
        },
    }
    return (
        dataset,
        files["topology"],
        files["restart"],
        {str(p): digest(p) for p in (manifest_path, *files.values())},
    )


def engines(config: dict[str, Any], base: Path) -> tuple[Path, Path]:
    """Require explicit AMBERHOME and a pmemd.cuda executable; no CPU fallback."""
    resolved = []
    for key in ("amber_home", "cuda_engine"):
        raw = os.path.expandvars(os.path.expanduser(config[key]))
        if "$" in raw:
            raise ValueError(f"Unresolved explicit environment variable in {key}")
        path = Path(raw)
        resolved.append(
            (base / path).resolve(strict=True)
            if not path.is_absolute()
            else path.resolve(strict=True)
        )
    home, engine = resolved
    if not engine.name.startswith("pmemd.cuda") or not os.access(engine, os.X_OK):
        raise ValueError(
            "An explicit executable pmemd.cuda GPU engine is required; CPU fallback is forbidden"
        )
    if not (home / "lib").is_dir():
        raise FileNotFoundError(home / "lib")
    return home, engine


def require_gpu(gpu: str) -> dict[str, str]:
    """Require NVIDIA driver visibility for the chosen device before launch/execution."""
    result = subprocess.run(
        ["nvidia-smi", "--query-gpu=index,name,uuid", "--format=csv,noheader"],
        check=True,
        capture_output=True,
        text=True,
    )
    lines = result.stdout.splitlines()
    matches = [line for line in lines if line.split(",", 1)[0].strip() == gpu]
    if len(matches) != 1:
        raise ValueError(f"GPU{gpu} unavailable; no CPU fallback permitted")
    return {
        "gpu_index": gpu,
        "device": matches[0],
        "check": "nvidia-smi driver visibility; pmemd.cuda failure remains fatal",
    }
