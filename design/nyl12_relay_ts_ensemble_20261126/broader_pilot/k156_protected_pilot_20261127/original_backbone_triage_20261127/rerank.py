"""Rerank existing MD-backbone candidates without RFD3 context scores."""
from __future__ import annotations
import copy
import csv
import hashlib
import json
import logging
from pathlib import Path
import sys
from nyl12_covalent.relay_broad_runtime import rank_broad


def main(campaign: Path) -> None:
    """Retain source hashes, original-parent scores and diagnostic ranking."""
    manifest_path = campaign / "manifest.json"
    ranked_path = campaign / "generation_run01/ranked_candidates.json"
    manifest = json.loads(manifest_path.read_text())
    originals = set(manifest["request"]["original_parent_ids"])
    candidates = []
    for source in json.loads(ranked_path.read_text()):
        if source["parent"] not in originals:
            continue
        if len(source["lineage"]) != 1:
            raise ValueError("Mixed generation lineages need explicit review")
        row = copy.deepcopy(source)
        row["scores"] = [s for s in row["scores"] if s["state"] == "ES" or s["template"] in originals]
        if {s["template"] for s in row["scores"] if s["state"] == "TS"} != originals:
            raise ValueError("Missing original TS cross-score")
        candidates.append(row)
    if len(candidates) != 12:
        raise ValueError("Expected exactly12 original-backbone candidates")
    ranked = rank_broad(candidates)
    selected = []
    for parent in manifest["request"]["original_parent_ids"]:
        selected.append(next(r["candidate_id"] for r in ranked if r["parent"] == parent))
    payload = {
        "status": "diagnostic_sequence_triage_not_activity_prediction",
        "source_sha256": {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in [manifest_path, ranked_path, Path(__file__)]},
        "candidate_count": 12,
        "rank_basis": "original-parent local TS mean NLL, worst NLL, mutation count, identity",
        "two_for_structural_review": selected,
        "selection_caveat": "One peroriginalgenerationparent; final selection requires sequence diversity,fold,maturation,assembly andmatchedWTgeometryreview",
        "candidates": ranked,
    }
    out = Path(__file__).parent
    with (out / "ranked_original_candidates.json").open("x") as handle:
        json.dump(payload, handle, indent=2, allow_nan=False)
        handle.write("\n")
    fields = ["candidate_id", "parent", "mutation_count", "local_ts_mean_nll", "local_ts_worst_nll"]
    with (out / "ranking.csv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: row[k] for k in fields} for row in ranked)
    logging.info("Retained %d original-backbone candidates; structural-review pair %s", len(ranked), selected)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main(Path(sys.argv[1]).resolve())
