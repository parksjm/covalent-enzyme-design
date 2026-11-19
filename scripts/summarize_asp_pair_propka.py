"""Extract Nyl12 catalytic pKa values from a recorded PROPKA screen."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

TARGETS = {("ASP", 40): "Asp272", ("ASP", 42): "Asp274", ("N+", 1): "Thr233_N"}
PATTERN = re.compile(r"^\s*(ASP|N\+)\s+(\d+)\s+([E-H])\s+(-?\d+\.\d+)\s+")


def sha256(path: Path) -> str:
    """Return the hexadecimal SHA-256 of a small source file.

    Args:
        path: Existing source file.

    Returns:
        SHA-256 hex digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_summary(path: Path) -> dict[tuple[str, str], float]:
    """Read only the final PROPKA pKa summary table.

    Args:
        path: PROPKA ``.pka`` output.

    Returns:
        Predicted pKa keyed by beta chain and catalytic group.
    """
    text = path.read_text()
    if text.count("SUMMARY OF THIS PREDICTION") != 1:
        raise ValueError(f"Missing or duplicated summary: {path}")
    lines = text.split("SUMMARY OF THIS PREDICTION", 1)[1].splitlines()
    values: dict[tuple[str, str], float] = {}
    for line in lines:
        match = PATTERN.match(line)
        if match is None:
            continue
        label = TARGETS.get((match.group(1), int(match.group(2))))
        if label is None:
            continue
        key = (match.group(3), label)
        if key in values:
            raise ValueError(f"Duplicate catalytic pKa {key} in {path}")
        values[key] = float(match.group(4))
    if len(values) != 12:
        raise ValueError(f"Expected 12 catalytic pKa values, found {len(values)}")
    return values


def main() -> None:
    """Write CSV and JSON summaries next to the immutable PROPKA outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen", type=Path, required=True)
    args = parser.parse_args()
    manifest = args.screen / "manifest.json"
    runs = json.loads(manifest.read_text())
    rows = []
    source_hashes = {}
    for run in runs:
        path = args.screen / run["case"] / "protein.pka"
        source_hashes[str(path)] = sha256(path)
        for (chain, group), pka in sorted(parse_summary(path).items()):
            rows.append(
                {
                    "case": run["case"],
                    "state": run["state"],
                    "frame": run["frame"],
                    "chain": chain,
                    "group": group,
                    "pka": pka,
                }
            )
    csv_path = args.screen / "catalytic_pka.csv"
    json_path = args.screen / "catalytic_pka.json"
    if csv_path.exists() or json_path.exists():
        raise FileExistsError("Refusing to overwrite existing pKa summaries")
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(
        json.dumps(
            {
                "rows": rows,
                "source_sha256": source_hashes,
                "screen_manifest_sha256": sha256(manifest),
                "interpretation": "Protein-only, static PROPKA predictions; not microstate populations.",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
