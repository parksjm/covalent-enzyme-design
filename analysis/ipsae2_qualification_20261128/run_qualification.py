"""Qualify pinned ipSAE2 on native WT/design data and independent formulas."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import subprocess
from typing import Any

import numpy as np
from ipsae.readers import load_model

from nyl12_covalent.ipsae2_scoring import (
    file_sha256,
    runtime_digest,
    score_model_raw,
)

LOGGER = logging.getLogger(__name__)
TOLERANCE = 5.1e-7  # Upstream TSV reports six decimal places.


def write_new(path: Path, value: Any) -> None:
    """Write immutable finite qualification evidence."""
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def direct_score(
    pae: np.ndarray, tokens: Any, left: list[str], right: list[str]
) -> float | None:
    """Evaluate the documented directional formula without calling scores.py.

    Args:
        pae: Original matrix, row alignment frame and column scored token.
        tokens: Read-only upstream token identities and frame-validity flags.
        left: Native protein/ligand alignment chains.
        right: Native scored chains, all of one molecular type here.

    Returns:
        Maximum frame score, or None when no valid frame exists.
    """
    frames = np.flatnonzero(
        tokens["chain"].isin(left).to_numpy() & tokens["frame_valid"].to_numpy()
    )
    scored = np.flatnonzero(tokens["chain"].isin(right).to_numpy())
    types = set(tokens.iloc[scored]["type"])
    if len(types) != 1 or not len(scored):
        raise ValueError("Qualification requires nonempty single-type scored groups")
    if not len(frames):
        return None
    sub = np.asarray(pae[np.ix_(frames, scored)], dtype=float)
    if types == {"prot"}:
        mask = sub <= 10.0
        counts = mask.sum(axis=1)
        scale = np.maximum(1.0, 1.24 * np.cbrt(np.maximum(counts, 1) - 15.0) - 1.8)
    elif types == {"lig"}:
        mask = np.ones(sub.shape, dtype=bool)
        counts = mask.sum(axis=1)
        scale = np.full(len(frames), 2.0)
    else:
        raise ValueError(f"Unqualified scored molecular types: {types}")
    kernel = 1.0 / (1.0 + np.square(sub / scale[:, None]))
    means = np.divide(
        np.where(mask, kernel, 0).sum(axis=1),
        counts,
        out=np.zeros(len(frames)),
        where=counts > 0,
    )
    return float(np.max(means))


def verify_numerics(model: dict[str, Any], result: dict[str, Any]) -> dict[str, Any]:
    """Check raw-axis identity and every directional TSV row numerically."""
    paths = [model["artifacts"][key]["path"] for key in ("cif", "pae", "confidence")]
    data = load_model(paths, verbose=False)
    with np.load(paths[1], allow_pickle=False) as loaded:
        raw = loaded["pae"]
    np.testing.assert_array_equal(data.pae, raw)
    if not np.isfinite(raw).all() or (raw < 0).any():
        raise ValueError("Invalid raw PAE")
    errors, comparisons, grouped_overrides = [], 0, 0
    expected_rows = {"chain_pairs": 66, "whole_protein": 70, "protomers": 88}
    for run, output in result["runs"].items():
        if len(output["rows"]) != expected_rows[run]:
            raise ValueError(f"Unexpected all-pair count for {run}")
        for row in output["rows"]:
            source = row["raw"]
            left, right = source["chn1"].split("+"), source["chn2"].split("+")
            for key, a, b in (("ipsae12", left, right), ("ipsae21", right, left)):
                expected = direct_score(raw, data.tokens, a, b)
                observed = row[key]
                if expected is None or observed is None:
                    if expected is not observed:
                        raise ValueError("Undefined-direction policy mismatch")
                else:
                    error = abs(expected - observed)
                    if error > TOLERANCE:
                        raise ValueError(f"Directional numerical mismatch {error}")
                    errors.append(error)
                comparisons += 1
            if run != "chain_pairs" and row["protein_aligned_direction"]:
                if row["ipsae"] != row["protein_aligned_ipsae"]:
                    grouped_overrides += 1
    for ligand, rows in result["ligand_profiles"].items():
        ordinary = [
            r["protein_aligned_ipsae"] for r in rows if r["run"] == "chain_pairs"
        ]
        whole = [
            r["protein_aligned_ipsae"]
            for r in rows
            if r["run"] == "whole_protein" and "+" in r["raw"]["chn1"]
        ]
        if len(ordinary) != 8 or len(whole) != 1:
            raise ValueError(f"Missing native per-ligand coverage: {ligand}")
        if abs(max(ordinary) - whole[0]) > TOLERANCE:
            raise ValueError("Whole-protein grouping differs from frame union")
    return {
        "directional_values_checked": comparisons,
        "max_absolute_error": max(errors),
        "TSV_rounding_tolerance": TOLERANCE,
        "raw_PAE_no_transposition": True,
        "grouped_scalar_requires_direction_override_count": grouped_overrides,
        "all_four_ligands_retained": True,
    }


def qualify(base: Path, output: Path) -> None:
    """Run scoring, independent arithmetic, replay and source preservation checks."""
    output.mkdir()
    runtime = json.loads((base / "runtime.json").read_text())
    upstream = json.loads(
        (base / "environment/qualification_acceptance.json").read_text()
    )
    pairs = {
        m["case_id"]: m
        for m in json.loads((base / "physical_pairing.json").read_text())["models"]
    }
    models = json.loads((base / "qualification_models.json").read_text())["models"]
    records = []
    for source in models:
        case = source["case"]
        mapping_path = base / f"{case}_native_mapping.json"
        mapping = json.loads(mapping_path.read_text())
        mapping.update(
            protomer_groups=pairs[case]["protomer_groups"],
            physical_pairs_verified=True,
            physical_pairing_evidence_sha256=file_sha256(
                base / "physical_pairing.json"
            ),
        )
        model = {"case_id": case, "artifacts": mapping["artifacts"]}
        before = {k: file_sha256(v["path"]) for k, v in mapping["artifacts"].items()}
        scored = output / case
        result = score_model_raw(model, mapping, runtime, scored)
        checks = verify_numerics(model, result)
        process = json.loads((scored / "chain_pairs_process.json").read_text())
        repeat = scored / "chain_pairs_repeat.tsv"
        argv = list(process["command"])
        argv[argv.index("--output-file") + 1] = str(repeat)
        with (scored / "repeat.log").open("x") as log:
            subprocess.run(argv, check=True, stdout=log, stderr=subprocess.STDOUT)
        if repeat.read_bytes() != (scored / "chain_pairs.tsv").read_bytes():
            raise ValueError("Repeated scoring is not byte-identical")
        after = {k: file_sha256(v["path"]) for k, v in mapping["artifacts"].items()}
        if after != before:
            raise ValueError("Qualification modified a native artifact")
        record = {
            "case_id": case,
            "mapping_receipt_sha256": file_sha256(mapping_path),
            "checks": checks,
            "repeat_byte_identical": True,
            "native_artifacts_unchanged": True,
        }
        records.append(record)
        write_new(scored / "numerical_qualification.json", record)
        LOGGER.info("Qualified %s: %s", case, checks)
    acceptance = {
        "status": "accepted",
        "runtime_sha256": runtime_digest(runtime),
        "qualified_roles": ["WT", "design"],
        "upstream_test_receipt_sha256": file_sha256(
            base / "environment/qualification_acceptance.json"
        ),
        "upstream_test_receipt": upstream,
        "models": records,
        "runner_sha256": file_sha256(Path(__file__)),
        "scope": "Technical input mapping and scoring only; not catalytic validation or a PA6 confidence cutoff",
    }
    write_new(output / "qualification_acceptance.json", acceptance)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    qualify(args.base.resolve(), args.output.resolve())
