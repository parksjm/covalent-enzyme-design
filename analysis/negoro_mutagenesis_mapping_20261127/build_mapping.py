"""Reproduce the preserved NylC-to-Nyl12 mapping without realigning sequences."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path
import platform
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
ALIGNMENT = Path(
    "analysis/nyl12_active_site_conservation_20261110/nyl12_nylc_numbering_alignment.fasta"
)
REFERENCE = Path(
    "analysis/nyl12_active_site_conservation_20261110/public_nylc_3axg.fasta"
)
INVENTORY = Path(
    "docs/figures/nyl12_sequence_protection_20261127/canonical_position_inventory.json"
)
ANCHORS = {
    134: "F101",
    146: "F113",
    189: "K156",
    219: "N186",
    266: "N232",
    267: "T233",
    306: "D272",
    307: "G273",
    308: "D274",
}


def sha256(path: Path) -> str:
    """Return an exact-byte SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_mapping(root: Path) -> dict[str, Any]:
    """Count nongap residues in both preserved alignment rows.

    Args:
        root: Repository root containing the three immutable source inputs.

    Returns:
        The original mapping payload with all 355 reference positions.

    Raises:
        ValueError: If source sequences, lengths, or anchor mappings disagree.
    """
    alignment, reference, inventory = [
        root / path for path in (ALIGNMENT, REFERENCE, INVENTORY)
    ]
    lines = alignment.read_text().splitlines()
    query, target = lines[1], lines[3]
    rows = json.loads(inventory.read_text())
    query_sequence = "".join(row["wt"] for row in rows)
    target_sequence = "".join(
        line for line in reference.read_text().splitlines() if not line.startswith(">")
    )
    if not (
        len(query) == len(target)
        and query.replace("-", "") == query_sequence
        and target.replace("-", "") == target_sequence
        and len(query_sequence) == 328
        and len(target_sequence) == 355
    ):
        raise ValueError(
            "Preserved alignment does not match the complete source sequences"
        )
    query_position = target_position = 0
    output = []
    for column, (query_aa, target_aa) in enumerate(zip(query, target), 1):
        if query_aa != "-":
            query_position += 1
        if target_aa != "-":
            target_position += 1
        if target_aa == "-":
            continue
        record = rows[query_position - 1] if query_aa != "-" else None
        output.append(
            dict(
                nylc_position=target_position,
                nylc_aa=target_aa,
                nylc_label=f"{target_aa}{target_position}",
                nyl12_position=query_position if record else None,
                nyl12_aa=query_aa if record else None,
                nyl12_label=record["wt_label"] if record else None,
                protected=record["protected"] if record else None,
                protection_reasons=record["reasons"] if record else None,
                alignment_column=column,
                identity=query_aa == target_aa,
                assignment="homolog_alignment" if record else "query_gap",
            )
        )
    if len(output) != 355 or any(
        output[position - 1]["nyl12_label"] != label
        for position, label in ANCHORS.items()
    ):
        raise ValueError("Mapping length or source-specific anchor validation failed")
    return dict(
        method="Count nongap residues independently in both rows of existing global sequence alignment; match2 mismatch-1 gap-2; no new alignment or constant-offset transformation",
        source_reference="3AXG NylCA/NylC 355-residue reference; not mutant7YU0/7YU1",
        limitations="Homolog sequence correspondence; uncertain local alignments and structural equivalence require review. Experimental variants should be mapped to this sequence only after verifying source construct numbering and identity.",
        source_hashes={
            str(path.relative_to(root)): sha256(path)
            for path in (alignment, reference, inventory)
        },
        query_matches_current_census_WT=True,
        query_length=328,
        reference_length=355,
        rows=output,
    )


def main() -> None:
    """Write a new mapping directory and verify against preserved evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--compare",
        type=Path,
        default=Path(__file__).with_name("nylc_to_nyl12_protection_mapping.json"),
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    root = args.root.resolve()
    payload = build_mapping(root)
    if payload != json.loads(args.compare.read_text()):
        raise ValueError("Reproduced payload differs from preserved evidence")
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output / "nylc_to_nyl12_protection_mapping.json"
    output.write_text(json.dumps(payload, indent=2) + "\n")
    with output.with_suffix(".csv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(payload["rows"][0]))
        writer.writeheader()
        writer.writerows(payload["rows"])
    provenance = dict(
        python=platform.python_version(),
        implementation=platform.python_implementation(),
        script_sha256=sha256(Path(__file__)),
        source_hashes=payload["source_hashes"],
        comparison=str(args.compare.resolve()),
        comparison_sha256=sha256(args.compare),
        exact_payload_match=True,
        rows=355,
        anchors=ANCHORS,
        alignment_recomputed=False,
    )
    (args.output / "reproduction.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    logging.info(
        "Exact355-row payload reproduced; nine anchors verified; original files unchanged"
    )


if __name__ == "__main__":
    main()
