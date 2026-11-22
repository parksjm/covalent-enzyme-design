"""Validation helpers for frozen generation inputs and mutations."""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
from typing import Any

ARMS = ("focused_19_positions", "broader_29_positions")
SEEDS = [2026112201, 2026112202, 2026112203, 2026112204, 2026112205]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _msa_hash(source: Path, native: str, sequence: str) -> str:
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields, rows = reader.fieldnames, list(reader)
    if not fields or not {"key", "sequence"}.issubset(fields):
        raise ValueError(f"MSA {source} must have key and sequence columns")
    if sum(row["key"] == "0" and row["sequence"] == native for row in rows) != 1:
        raise ValueError(f"MSA {source} must have one native query row")
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        changed = dict(row)
        if row["key"] == "0":
            changed["sequence"] = sequence
        writer.writerow(changed)
    return hashlib.sha256(buffer.getvalue().encode()).hexdigest()


def _yaml_identity(path: Path) -> dict[str, Any]:
    data = _read_json(path)
    for entry in data["sequences"]:
        protein = entry.get("protein")
        if protein:
            msa = Path(protein["msa"])
            if not msa.is_file():
                raise FileNotFoundError(msa)
            protein["msa"] = {"sha256": _sha(msa)}
    for template in data.get("templates", []):
        cif = Path(template["cif"])
        if not cif.is_file():
            raise FileNotFoundError(cif)
        template["cif"] = {"sha256": _sha(cif)}
    return data


def _manifest_masks(
    path: Path,
) -> tuple[
    dict[str, dict[str, set[int]]],
    dict[str, set[int]],
    dict[str, Any],
    dict[str, Path],
]:
    manifest = _read_json(path)
    if manifest.get("expected_jobs") != 8 or len(manifest.get("jobs", [])) != 8:
        raise ValueError("Generation manifest must contain eight frozen jobs")
    # The frozen source snapshot keeps this preflight beside the manifest.
    protocol_root = path.resolve().parent.parent
    protocol = _read_json(protocol_root / "protocol.json")
    if tuple(protocol["arms"]) != ARMS or len(protocol["parents_proposed"]) != 2:
        raise ValueError("Frozen mask protocol has unexpected arms or parent count")
    parent_proposal = _read_json(
        protocol_root / "preflight/parent_and_mask_proposal.json"
    )
    approved_parent_hashes = {row["sha256"] for row in parent_proposal["parents"]}
    frozen = _read_json(protocol_root / "preflight/final_masks.json")
    masks: dict[str, dict[str, set[int]]] = {}
    protected: dict[str, set[int]] | None = None
    seen: set[tuple[str, str, float]] = set()
    parent_paths: dict[str, Path] = {}
    for job in manifest["jobs"]:
        argv = job["argv"]
        try:
            temperature = float(argv[argv.index("--temperature") + 1])
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"Generation job {job.get('job_id')} lacks temperature"
            ) from exc
        arm, parent = job["arm"], str(job["parent"])
        if arm not in ARMS:
            raise ValueError(f"Unknown design arm: {arm}")
        try:
            pdb_path = Path(argv[argv.index("--pdb_path") + 1]).resolve()
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"Generation job {job.get('job_id')} lacks parent PDB"
            ) from exc
        if _sha(pdb_path) not in approved_parent_hashes:
            raise ValueError(
                f"Generation job {job.get('job_id')} uses an unapproved parent PDB"
            )
        if parent in parent_paths and parent_paths[parent] != pdb_path:
            raise ValueError(f"Generation parent {parent} maps to multiple PDBs")
        parent_paths[parent] = pdb_path
        key = (arm, parent, temperature)
        if key in seen:
            raise ValueError(f"Duplicate generation stratum: {key}")
        seen.add(key)
        audit = job.get("audit_manifest", {})
        raw_mask, config = audit.get("design_mask", {}), audit.get("config", {})
        current = {
            "alpha": set(raw_mask.get("alpha_positions", [])),
            "beta": set(raw_mask.get("beta_positions", [])),
        }
        current_protected = {
            "alpha": set(config.get("protected_alpha_positions", [])),
            "beta": set(config.get("protected_beta_positions", [])),
        }
        if not all(current.values()) or not all(current_protected.values()):
            raise ValueError(
                f"Generation job {job.get('job_id')} has incomplete masks/protections"
            )
        expected_mask_key = (
            "focused_positions" if arm == ARMS[0] else "expanded_positions"
        )
        expected = {
            family: set(frozen[expected_mask_key][family])
            for family in ("alpha", "beta")
        }
        hard = {
            family: set(frozen["hard_protected_positions"][family])
            for family in ("alpha", "beta")
        }
        if current != expected:
            raise ValueError(
                f"Generation arm {arm} differs from frozen {expected_mask_key}"
            )
        if current_protected != hard:
            raise ValueError(
                "Generation hard protections differ from the audited frozen protections"
            )
        if arm in masks and masks[arm] != current:
            raise ValueError(f"Mask varies across source jobs in arm {arm}")
        if protected is not None and protected != current_protected:
            raise ValueError("Hard protections vary across generation jobs")
        masks[arm], protected = current, current_protected
    expected_strata = {
        (arm, parent, float(temp))
        for arm in ARMS
        for parent in parent_paths
        for temp in protocol["temperatures"]
    }
    if (
        seen != expected_strata
        or len(parent_paths) != 2
        or {_sha(parent_path) for parent_path in parent_paths.values()}
        != approved_parent_hashes
        or protected is None
    ):
        raise ValueError(
            "Generation manifest does not cover the full two-by-two-by-two design"
        )
    return masks, protected, manifest, parent_paths


def _candidate_mutations(
    alpha: str, beta: str, native_alpha: str, native_beta: str
) -> list[dict[str, Any]]:
    result = []
    for family, sequence, native, offset in (
        ("alpha", alpha, native_alpha, 0),
        ("beta", beta, native_beta, 232),
    ):
        for pos, (before, after) in enumerate(zip(native, sequence, strict=True), 1):
            if before != after:
                result.append(
                    {
                        "family": family,
                        "physical_position": pos,
                        "merged_position": pos + offset,
                        "from": before,
                        "to": after,
                    }
                )
    return result


def _yaml_identity_data(
    data: dict[str, Any], alpha_hash: str, beta_hash: str
) -> dict[str, Any]:
    """Normalize a prepared YAML to scientific content hashes, not paths."""
    value = json.loads(json.dumps(data))
    proteins = [entry["protein"] for entry in value["sequences"] if "protein" in entry]
    if len(proteins) != 2:
        raise ValueError("Four-ligand template must contain two proteins")
    for protein, digest in zip(proteins, (alpha_hash, beta_hash), strict=True):
        protein["msa"] = {"sha256": digest}
    for template in value.get("templates", []):
        cif = Path(template["cif"])
        if not cif.is_file():
            raise FileNotFoundError(cif)
        template["cif"] = {"sha256": _sha(cif)}
    return value
