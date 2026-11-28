"""Review assisted relay hypotheses independently of direct-route control exports."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nyl12_covalent.leaving_group_relay_assisted import (  # noqa: E402
    audit_shortlist,
    read_assisted,
    review_source_hashes,
    review_status,
    select_assisted,
)
from nyl12_covalent.leaving_group_relay_export import (  # noqa: E402
    digest,
    evaluate_native,
    export_frame,
)
from nyl12_covalent.leaving_group_relay_selection import refinement_frames  # noqa: E402


def main() -> None:
    """Retain raw classifications, audit ≤20 frames and export labeled hypotheses."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--screen-dir", required=True, action="append", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    config = json.loads(args.config.read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    rows, provenance, totals = read_assisted(
        args.screen_dir, config, args.output / "all_classified_paths.jsonl"
    )
    checked = evaluate_native(audit_shortlist(rows), config)
    selected = select_assisted(checked)
    cases = {case["name"]: case for case in config["cases"]}
    exports = [
        export_frame(
            row, cases[row["case"]], args.output, f"{row['selection_role']}_{i}"
        )
        for i, row in enumerate(selected, 1)
    ]
    frames = refinement_frames(checked, totals)
    names = {row["case"] for row in frames}
    refined_config = {
        **config,
        "stride": 1,
        "cases": [case for case in config["cases"] if case["name"] in names],
    }
    status = review_status(selected)
    report = {
        "status": status,
        "exports": exports,
        "native_frames_evaluated": len(checked),
        "distinct_assisted_frames": len(rows),
        "config_sha256": digest(args.config),
        "source_sha256": review_source_hashes(Path(__file__)),
        "provenance": provenance,
        "mechanism_validated": False,
        "direct_controls_are_new_relays": False,
    }
    for name, value in (
        ("ranked_assisted_candidates.json", rows),
        ("native_audit.json", checked),
        ("selection.json", report),
        ("refinement_frames.json", frames),
        ("refinement_cases.json", refined_config),
        ("config.json", config),
    ):
        (args.output / name).write_text(json.dumps(value, indent=2) + "\n")
    lines = [
        f"Status: {status}",
        "",
        "Assisted hypotheses require at least one water or internal residue.",
        "Direct OG1→leaving-N routes are counted separately as controls.",
        "Only strict all-leg assisted paths with canonical Tier A/B are jointly qualified.",
        "Orientation nearpasses and strict non-Tier paths remain diagnostic; native failures remain explicit.",
        "Conditional Nalpha activation and leaving-amide acceptance do not establish a complete mechanism.",
        "",
        f"Path counts: {json.dumps(provenance['path_counts'], sort_keys=True)}",
        "",
        "|Case|Frame|Export role|Source class|Violation|Tier|Native pass|",
        "|---|---:|---|---|---:|---|---|",
    ]
    lines.extend(
        f"|{r['case']}|{r['frame']}|{r['selection_role']}|{r['source_class']}|{r['violation']}|{r['tier']}|{r['native_pass']}|"
        for r in selected
    )
    (args.output / "summary.md").write_text("\n".join(lines) + "\n")
    logging.info("%s; %d native frames, %d exports", status, len(checked), len(exports))


if __name__ == "__main__":
    main()
