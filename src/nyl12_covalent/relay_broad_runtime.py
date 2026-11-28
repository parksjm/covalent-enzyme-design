"""Strict exports, native execution receipts, and scores for the broad pilot."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np

from .relay_ensemble_sequences import summarize_score


def write_new(path: Path, data: Any) -> None:
    """Write a new immutable JSON artifact, refusing replacement."""
    with path.open("x") as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write("\n")


def digest(path: Path) -> str:
    """Return the exact file digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def journal(root: Path, phase: str, action: str, reason: str) -> None:
    """Append a required action event, failing closed if logging fails."""
    subprocess.run(
        [
            sys.executable,
            str(root / "scripts/record_agent_action.py"),
            "--log",
            str(root / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "broad_controller",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            reason,
            "--rollback",
            "Retain immutable artifacts; stop only verified controller PID",
        ],
        check=True,
    )


def native_run(
    root: Path, argv: list[str], output: Path, environment: dict[str, str]
) -> int:
    """Retain one exact command, log and independent process receipt."""
    journal(root, "intent", "Execute " + str(output), json.dumps(argv))
    write_new(
        output.with_suffix(".command.json"), {"argv": argv, "environment": environment}
    )
    started = datetime.now(timezone.utc).isoformat()
    with output.with_suffix(".log").open("x") as handle:
        result = subprocess.run(
            argv,
            env={**os.environ, **environment},
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    write_new(
        output.with_suffix(".process.json"),
        {
            "exit_code": result.returncode,
            "started_utc": started,
            "ended_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    journal(root, "outcome", "Execute " + str(output), f"exit={result.returncode}")
    return result.returncode


def verify_parent(parent: dict[str, Any]) -> None:
    """Require exact hashes and independently accepted chemistry receipts."""
    for field in (
        "pdb_path",
        "canonical_mapping",
        "pa6_graph_receipt",
        "scientific_acceptance",
    ):
        if digest(Path(parent[field])) != parent[field + "_sha256"]:
            raise ValueError(f"Accepted parent hash mismatch: {field}")
    for field in ("pa6_graph_receipt", "scientific_acceptance"):
        if json.loads(Path(parent[field]).read_text())["accepted"] is not True:
            raise ValueError(f"Parent rejected by {field}")
    if parent["sequence_scope"] != "native_WT":
        raise ValueError("Only native-sequence RFD3 parents are permitted")


def export_local(source: Path, destination: Path) -> None:
    """Export beta-local numbering, preserving every atom and coordinate."""
    lines = []
    for line in source.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[21] in "EFGH":
            number = int(line[22:26])
            if number > 232:
                number -= 232
            if not 1 <= number <= 96:
                raise ValueError("Invalid beta residue numbering")
            line = line[:22] + f"{number:4d}" + line[26:]
        lines.append(line)
    with destination.open("x") as handle:
        handle.write("\n".join(lines) + "\n")


def parser_audit(pdb: Path, engine: Path) -> dict[str, Any]:
    """Run the actual native CPU parser and check full catalytic context."""
    sys.path.insert(0, str(engine.parent))
    import data_utils

    protein = data_utils.parse_PDB(
        str(pdb),
        device="cpu",
        parse_all_atoms=True,
        parse_atoms_with_zero_occupancy=True,
    )[0]
    names = [
        f"{c}{p}"
        for c, p in zip(
            protein["chain_letters"], protein["R_idx"].tolist(), strict=True
        )
    ]
    expected = [
        f"{c}{p}" for c in "ABCDEFGH" for p in range(1, 233 if c in "ABCD" else 97)
    ]
    if names != expected or not bool(protein["mask"].all()):
        raise ValueError("Native parser did not preserve complete canonical chains")
    if int(protein["Y_m"].sum()) != 33:
        raise ValueError("Native parser must retain exactly 33 ligand heavy atoms")
    fixed = [
        ("A113", "N", 0),
        ("A156", "NZ", 35),
        ("A186", "ND2", 15),
        ("E1", "N", 0),
        ("E1", "OG1", 9),
        ("E40", "OD1", 16),
        ("E40", "OD2", 17),
        ("E42", "OD1", 16),
        ("E42", "OD2", 17),
    ]
    for residue, atom, index in fixed:
        if not protein["xyz_37_m"][names.index(residue), index]:
            raise ValueError(f"Native parser lost fixed sidechain {residue}:{atom}")
    return {
        "accepted": True,
        "protein_rows": len(names),
        "ligand_heavy_atoms": 33,
        "fixed_context": [f"{r}:{a}" for r, a, _ in fixed],
        "pdb_sha256": digest(pdb),
        "zero_occupancy_parsing": True,
    }


def score_summary(
    stats: dict[str, Any],
    sequence: str,
    global_mask: dict[str, list[int]],
    local_mask: dict[str, list[int]],
) -> dict[str, float]:
    """Report global mutable and distinct local pocket likelihood diagnostics."""
    broad = summarize_score(stats, sequence, global_mask)
    local = summarize_score(stats, sequence, local_mask)
    return {**broad, **{f"local_{key}": value for key, value in local.items()}}


def rank_broad(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Preselect by balanced local-pocket diagnostic, explicitly not activity."""
    candidates = [item for item in candidates if item["mutation_count"] > 0]
    if not candidates:
        raise ValueError("No non-WT broad candidate was generated")
    for item in candidates:
        ts = [s["local_AE_design_nll"] for s in item["scores"] if s["state"] == "TS"]
        if not ts or not np.isfinite(ts).all():
            raise ValueError("Missing finite local-pocket TS scores")
        item["local_ts_mean_nll"] = float(np.mean(ts))
        item["local_ts_worst_nll"] = float(np.max(ts))
    return sorted(
        candidates,
        key=lambda x: (
            x["local_ts_mean_nll"],
            x["local_ts_worst_nll"],
            x["mutation_count"],
            x["candidate_id"],
        ),
    )
