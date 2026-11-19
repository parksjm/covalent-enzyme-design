"""Copy compact pKa evidence and filter only accepted precursor table identifiers."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path("/data/jpc/covalent-enzyme-design")
OUT = Path(__file__).resolve().parent


def digest(path: Path) -> str:
    """Return immutable source/destination SHA256."""
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main() -> None:
    """Copy requested records and retain exact row-filter/source-hash provenance."""
    records = []
    sources = [
        (
            "analysis/negoro_2023_protonation/crystal_propka_20261116/attempt3/results/asp306_asp308_pka.csv",
            "nylc_mature_corrected_attempt3.csv",
        ),
        (
            "qmmm/pa6_flipped_model30/09_propka_static_tier_a/screen_attempt1/catalytic_pka.csv",
            "nyl12_catalytic_pka_36_values.csv",
        ),
        (
            "analysis/negoro_2023_protonation/terminal_pka_20261116/before_after.json",
            "terminal_before_after.json",
        ),
        (
            "analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/scientific_acceptance.json",
            "mcce_nylc_scientific_acceptance.json",
        ),
        (
            "analysis/nyl12_mcce_step1_preflight_20261116/neutral_replica_02_frame284_attempt2/scientific_acceptance.json",
            "mcce_nyl12_preflight_scientific_acceptance.json",
        ),
    ]
    for relative, name in sources:
        source, target = ROOT / relative, OUT / name
        if target.exists():
            raise FileExistsError(target)
        if name == "nyl12_catalytic_pka_36_values.csv":
            with source.open(newline="") as stream:
                if len(list(csv.DictReader(stream))) != 36:
                    raise ValueError("Expected full36-value Nyl12 pKa table")
        shutil.copy2(source, target)
        if digest(source) != digest(target):
            raise ValueError("Copied pKa evidence differs from source")
        records.append(
            dict(
                source=str(source),
                source_sha256=digest(source),
                destination=name,
                destination_sha256=digest(target),
                bytes=target.stat().st_size,
                operation="byte_identical_copy",
            )
        )
    precursor = (
        ROOT
        / "analysis/negoro_2023_protonation/crystal_propka_20261116/attempt2/results/asp306_asp308_pka.csv"
    )
    with precursor.open(newline="") as stream:
        reader = csv.DictReader(stream)
        fields, rows = reader.fieldnames, list(reader)
    selected = [row for row in rows if row["structure"] in {"7YU0", "7YU1"}]
    if {row["structure"] for row in selected} != {"7YU0", "7YU1"}:
        raise ValueError("Both requested precursor structures must be present")
    target = OUT / "nylc_precursor_7YU0_7YU1_attempt2.csv"
    with target.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(selected)
    records.append(
        dict(
            source=str(precursor),
            source_sha256=digest(precursor),
            destination=target.name,
            destination_sha256=digest(target),
            bytes=target.stat().st_size,
            operation="filter_rows_only_no_recalculation",
            column="structure",
            included_identifiers=["7YU0", "7YU1"],
            source_rows=len(rows),
            output_rows=len(selected),
            excluded_identifiers=sorted(
                {row["structure"] for row in rows} - {"7YU0", "7YU1"}
            ),
        )
    )
    claims = {}
    for name in (
        "mcce_nylc_scientific_acceptance.json",
        "mcce_nyl12_preflight_scientific_acceptance.json",
    ):
        receipt = json.loads((OUT / name).read_text())
        claims[name] = dict(
            source_status=receipt.get("state", receipt.get("status")),
            accepted_pKa_result=False,
            scope="Retained negative/preflight evidence only",
        )
    mapping = dict(
        records=records,
        extraction_script_sha256=digest(Path(__file__)),
        obsolete_5XYG_attempt2_rows_included=False,
        recalculated_pKa_values=False,
        source_data_modified=False,
        jobs_launched=False,
        mcce_receipt_scope=claims,
    )
    with (OUT / "source_mapping.json").open("x") as stream:
        json.dump(mapping, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
