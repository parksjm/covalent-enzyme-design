"""Select bounded observed relay candidates and export unchanged source frames."""

from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nyl12_covalent.leaving_group_relay_export import (  # noqa: E402
    digest,
    evaluate_native,
    export_frame,
)
from nyl12_covalent.leaving_group_relay_selection import (  # noqa: E402
    final_selection,
    read_candidates,
    refinement_frames,
    shortlist,
)


def main() -> None:
    """Produce ranked diagnostics, refinement frames and up to four exact exports."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--screen-dir", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    config = json.loads(args.config.read_text())
    rows = read_candidates(args.screen_dir, config)
    checked = evaluate_native(shortlist(rows), config)
    selected = final_selection(checked)
    totals = {}
    source_hashes = {}
    for directory in args.screen_dir:
        for summary in directory.rglob("summary.json"):
            data = json.loads(summary.read_text())
            for case in data["cases"]:
                totals[case["case"]] = case["total_frames"]
                source_hashes.update(case["input_sha256"])
    if not totals:
        raise ValueError("No screen summary with total_frames found")
    args.output.mkdir(parents=True, exist_ok=False)
    frames = refinement_frames(checked, totals)
    (args.output / "refinement_frames.json").write_text(
        json.dumps(frames, indent=2) + "\n"
    )
    (args.output / "ranked_candidates.json").write_text(
        json.dumps(checked, indent=2) + "\n"
    )
    with (args.output / "ranked_candidates.csv").open("x", newline="") as handle:
        if checked:
            writer = csv.DictWriter(handle, fieldnames=list(checked[0]))
            writer.writeheader()
            writer.writerows(
                {
                    k: json.dumps(v) if isinstance(v, (list, dict)) else v
                    for k, v in row.items()
                }
                for row in checked
            )
    cases = {case["name"]: case for case in config["cases"]}
    exports = []
    for i, row in enumerate(selected, 1):
        for source in (cases[row["case"]]["topology"], row["trajectory"]):
            if digest(Path(source)) != source_hashes[source]:
                raise ValueError(f"Source differs from screened input: {source}")
        exports.append(
            export_frame(
                row, cases[row["case"]], args.output, f"{row['selection_role']}_{i}"
            )
        )
    status = "qualified_strict_tier_ab" if selected else "no_qualified_candidate"
    manifest = {
        "status": status,
        "exports": exports,
        "distinct_screen_frames": len(rows),
        "native_frames_evaluated": len(checked),
        "config_sha256": digest(args.config),
        "source_input_sha256": source_hashes,
    }
    (args.output / "selection.json").write_text(json.dumps(manifest, indent=2) + "\n")
    lines = [
        f"Status: {status}",
        "",
        f"Native geometry evaluated for {len(checked)} distinct frames; {len(selected)} selected.",
        "Selection requires strict all-leg hydrogen bonds and canonical Tier A/B; native failures remain labeled.",
        "Independent cases are preferred; frames within one case are separated by at least 1000 ps.",
        "NH2 Nalpha-start paths remain conditional activation hypotheses; activation-leg failures are labeled activation_unresolved. Full OG1-start paths rank first within a geometry tier. Signed phi orientations are retained. Leaving-amide N acceptance remains conditional; this screen does not establish a mechanism.",
        "",
        "|Case|Frame|Role|Tier|Phi orientation|Native pass|",
        "|---|---:|---|---|---|---|",
    ]
    lines.extend(
        f"|{r['case']}|{r['frame']}|{r['selection_role']}|{r['tier']}|{r['phi_orientation']}|{r['native_pass']}|"
        for r in selected
    )
    if not selected:
        lines.append(
            "No strict Tier A/B candidate; ranked rows are diagnostic nearpasses only."
        )
    (args.output / "summary.md").write_text("\n".join(lines) + "\n")
    logging.info(
        "%s: %d exports, %d refinement frames", status, len(exports), len(frames)
    )


if __name__ == "__main__":
    main()
