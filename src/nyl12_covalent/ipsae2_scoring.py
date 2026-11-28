"""Isolated upstream ipSAE2 execution on verified, unmodified native models."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REVISION = "fe560f4c2461eeb91aebded61ac8c1e83b264f60"


def file_sha256(path: str | Path) -> str:
    """Return the SHA256 of a file without loading it all into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def runtime_digest(runtime: dict[str, Any]) -> str:
    """Identify the exact configured executable and pinned source hashes."""
    return hashlib.sha256(json.dumps(runtime, sort_keys=True).encode()).hexdigest()


def verify_runtime(runtime: dict[str, Any]) -> None:
    """Reject unpinned sources or a changed isolated executable."""
    if runtime["revision"] != REVISION or runtime["version"] != "0.1.11":
        raise ValueError("ipSAE2 revision/version differs from reviewed release")
    if not runtime["command"] or not runtime["source_sha256"]:
        raise ValueError("Explicit isolated command and source hashes required")
    for name, expected in runtime["source_sha256"].items():
        if file_sha256(name) != expected:
            raise ValueError(f"Runtime hash mismatch: {name}")
    executable = str(Path(runtime["command"][0]).resolve())
    if executable not in {str(Path(p).resolve()) for p in runtime["source_sha256"]}:
        raise ValueError("Isolated executable must be included in source_sha256")


def require_qualification(receipt: dict[str, Any], runtime: dict[str, Any]) -> None:
    """Require accepted WT/design technical qualification for this runtime."""
    if receipt.get("status") != "accepted":
        raise ValueError("Accepted qualification receipt required for production")
    if receipt.get("runtime_sha256") != runtime_digest(runtime):
        raise ValueError("Qualification does not identify this runtime")
    if set(receipt.get("qualified_roles", [])) != {"WT", "design"}:
        raise ValueError("Qualification must include WT and eligible design")


def _number(value: str) -> float | None:
    if value in {"NA", "", "nan", "None"}:
        return None
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite upstream number: {value}")
    return number


def parse_scores(
    path: Path, protein_chains: list[str], ligand_chains: list[str]
) -> list[dict[str, Any]]:
    """Preserve raw rows and explicitly extract protein-aligned directions.

    Missing upstream scores stay null. The TSV cannot distinguish all causes
    of NA; no inferred physical cause or hard-zero replacement is assigned.
    """
    rows = []
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream, delimiter="\t"):
            left, right = raw["chn1"].split("+"), raw["chn2"].split("+")
            row: dict[str, Any] = {"raw": raw}
            for key in ("ipsae", "ipsae12", "ipsae21"):
                row[key] = _number(raw[key])
            row["undefined_reason"] = (
                "upstream_NA_cause_not_distinguished_in_TSV"
                if any(row[key] is None for key in ("ipsae", "ipsae12", "ipsae21"))
                else None
            )
            contacts = [
                _number(raw[k]) for k in ("n_scored_contact12", "n_scored_contact21")
            ]
            row["no_scored_contact"] = (
                all(value == 0 for value in contacts)
                if all(value is not None for value in contacts)
                else None
            )
            row["physical_contact_status"] = "unknown_from_TSV"
            direction = None
            if set(left) <= set(protein_chains) and set(right) <= set(ligand_chains):
                direction = "12"
            elif set(right) <= set(protein_chains) and set(left) <= set(ligand_chains):
                direction = "21"
            row["protein_aligned_direction"] = direction
            row["protein_aligned_ipsae"] = (
                row[f"ipsae{direction}"] if direction else None
            )
            row["protein_interface"] = set(left) <= set(protein_chains) and set(
                right
            ) <= set(protein_chains)
            rows.append(row)
    if not rows:
        raise ValueError(f"No score rows: {path}")
    return rows


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def score_model_raw(
    model_manifest: dict[str, Any],
    mapping_receipt: dict[str, Any],
    runtime: dict[str, Any],
    output_dir: Path,
) -> dict[str, Any]:
    """Score native coordinates for explicit qualification or a gated caller.

    This low-level API intentionally does not assert global qualification so
    an explicit WT/design qualification runner can generate that evidence.
    Production callers must first call :func:`require_qualification`.
    """
    verify_runtime(runtime)
    if mapping_receipt.get("accepted") is not True:
        raise ValueError("Accepted native mapping receipt required")
    inputs = model_manifest["artifacts"]
    for key in ("cif", "pae", "confidence"):
        artifact = inputs[key]
        if file_sha256(artifact["path"]) != artifact["sha256"]:
            raise ValueError(f"Native artifact hash mismatch: {key}")
    if any(
        mapping_receipt["artifacts"][key]["sha256"] != inputs[key]["sha256"]
        for key in ("cif", "pae", "confidence")
    ):
        raise ValueError("Mapping receipt does not match native triplet")
    proteins = mapping_receipt["protein_chains"]
    ligands = mapping_receipt["ligand_chains"]
    if len(proteins) != 8 or len(ligands) != 4 or len(set(proteins + ligands)) != 12:
        raise ValueError("Expected eight disjoint protein chains and four ligands")
    groups = mapping_receipt.get("protomer_groups")
    if groups is not None and (
        len(groups) != 4
        or any(len(group) != 2 for group in groups)
        or sorted(chain for group in groups for chain in group) != sorted(proteins)
        or mapping_receipt.get("physical_pairs_verified") is not True
    ):
        raise ValueError(
            "Protomer grouping requires verified physical alpha/beta pairs"
        )
    output_dir.mkdir(parents=True, exist_ok=False)
    _write_json(output_dir / "model_manifest.json", model_manifest)
    _write_json(output_dir / "mapping_receipt.json", mapping_receipt)
    _write_json(output_dir / "runtime.json", runtime)
    runs: dict[str, list[list[str]]] = {
        "chain_pairs": [],
        "whole_protein": [proteins],
    }
    if groups:
        runs["protomers"] = groups
    result: dict[str, Any] = {
        "status": "running",
        "case_id": model_manifest["case_id"],
        "runtime_sha256": runtime_digest(runtime),
        "runs": {},
        "interpretation": "Supplementary confidence; no catalytic or PA6 threshold claim",
        "ligand_profiles": {ligand: [] for ligand in ligands},
        "protomer_grouping": "verified"
        if groups
        else "not_requested_no_verified_pairs",
    }
    for name, group_list in runs.items():
        tsv = output_dir / f"{name}.tsv"
        command = [
            *runtime["command"],
            *[inputs[key]["path"] for key in ("cif", "pae", "confidence")],
            "--all-pairs",
            "--output-file",
            str(tsv.resolve()),
        ]
        for group in group_list:
            command.extend(["--group", "+".join(group)])
        record = {
            "command": command,
            "started_utc": datetime.now(timezone.utc).isoformat(),
        }
        _write_json(output_dir / f"{name}_process.json", record)
        with (output_dir / f"{name}.log").open("w") as log:
            process = subprocess.run(
                command, stdout=log, stderr=subprocess.STDOUT, check=False
            )
        record.update(
            returncode=process.returncode,
            finished_utc=datetime.now(timezone.utc).isoformat(),
        )
        _write_json(output_dir / f"{name}_process.json", record)
        if process.returncode:
            result.update(status="failed", failed_run=name)
            _write_json(output_dir / "scores.json", result)
            raise RuntimeError(
                f"ipSAE2 {name} exited {process.returncode}; see {output_dir}"
            )
        rows = parse_scores(tsv, proteins, ligands)
        result["runs"][name] = {"tsv_sha256": file_sha256(tsv), "rows": rows}
        for row in rows:
            for ligand in ligands:
                if (
                    ligand in (row["raw"]["chn1"], row["raw"]["chn2"])
                    and row["protein_aligned_direction"]
                ):
                    result["ligand_profiles"][ligand].append({"run": name, **row})
    if any(not rows for rows in result["ligand_profiles"].values()):
        raise ValueError("Missing per-ligand score profile; all four must be retained")
    result["status"] = "completed"
    _write_json(output_dir / "scores.json", result)
    return result
