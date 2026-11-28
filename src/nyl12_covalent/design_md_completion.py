"""Verify matched MD replica completion and selection eligibility."""

from __future__ import annotations
import hashlib
import json
from pathlib import Path

from nyl12_covalent.native_thr_policy import geometry_eligible


def sha(path: Path) -> str:
    """Return a streaming source SHA-256 digest."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_completion(
    replica_dir: Path, topology_hash: str, start_hash: str
) -> list[Path]:
    """Validate completed native production stages and their input identities.

    Args:
        replica_dir: Exact selected replica directory.
        topology_hash: Expected dataset topology SHA-256.
        start_hash: Expected minimized source restart SHA-256.

    Returns:
        The two ordered production trajectory paths.

    Raises:
        ValueError: If process, scientific stage or source identity gates fail.
    """
    receipt = json.loads((replica_dir / "100ns_process_complete.json").read_text())
    if receipt.get("completed") is not True:
        raise ValueError("Replica completion receipt is not accepted")
    dynamics = replica_dir / "dynamics"
    first = json.loads((dynamics / "05_production_10ns.process.json").read_text())
    continuation = json.loads((dynamics / "extension_100ns.json").read_text())
    if first.get("returncode") != 0 or continuation.get("returncode") != 0:
        raise ValueError("Production native process failed")
    if (
        first.get("reference_sha256") != start_hash
        or sha(dynamics / "system.prmtop") != topology_hash
    ):
        raise ValueError("Production topology or original start differs from dataset")
    for record, names in (
        (first, ("system.prmtop", "04_equil_npt.rst7", "05_production_10ns.in")),
        (
            continuation,
            ("system.prmtop", "05_production_10ns.rst7", "06_production_90ns.in"),
        ),
    ):
        recorded = record["input_sha256"]
        for name in names:
            matches = [
                value for path, value in recorded.items() if Path(path).name == name
            ]
            if len(matches) != 1 or sha(dynamics / name) != matches[0]:
                raise ValueError(f"Production input hash mismatch: {name}")
    for stem, steps, final in (
        ("05_production_10ns", 5000000, 10500.0),
        ("06_production_90ns", 45000000, 100500.0),
    ):
        gate = json.loads((dynamics / f"{stem}.acceptance.json").read_text())
        if (
            gate.get("completed_steps") != steps
            or gate.get("final_time_ps") != final
            or gate.get("temperature_gate_pass") is not True
        ):
            raise ValueError(f"Production scientific stage failed: {stem}")
        if sha(dynamics / f"{stem}.rst7") != gate["restart_sha256"]:
            raise ValueError(f"Production restart hash mismatch: {stem}")
    if continuation.get("total_production_ns") != 100:
        raise ValueError("Continuation duration is not 100 ns")
    return [
        dynamics / f"{stem}.nc" for stem in ("05_production_10ns", "06_production_90ns")
    ]


def add_selection_eligibility(row: dict) -> dict:
    """Add Tier A/B eligibility while preserving native-Thr diagnostic fields.

    Args:
        row: Saved original geometry and native-gate record.

    Returns:
        Geometry record with Tier A/B-only selection eligibility.
    """
    return {**row, "selection_eligible": geometry_eligible(row["tier"])}


def replica_path(follow: Path, dataset: str, replica: str) -> Path:
    """Resolve a supported dataset and exact replica directory.

    Args:
        follow: Matched MD follow-up directory.
        dataset: WT, broader_mask07, focused_mask03, or V111C.
        replica: Replica number 01 or 02.

    Returns:
        The explicit replica path.
    """
    if dataset not in (
        "WT",
        "broader_mask07",
        "focused_mask03",
        "V111C",
    ) or replica not in ("01", "02"):
        raise ValueError("Unsupported dataset or replica")
    return follow / "md" / dataset / f"replica_{replica}"
