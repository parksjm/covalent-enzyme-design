"""Preserve shipped Amber DFTB3 inputs and compare printed regression energies."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
AMBERHOME = Path("/data/jpc/code/ambertools26/venv")
SHIPPED = Path("/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/test/qmmm_DFTB/dftb3")
CASES = {
    "water": ("water", "Run.h2o.go", "h2o", "h2o.go.out.save"),
    "nma_water": ("nma-h2o-box", "Run.nma-spcfw-15.md", "nma-spcfw-15", "nma-spcfw-15.md.out.save"),
}


def record(path: Path) -> dict[str, str | int]:
    """Hash an unchanged source artifact."""
    before = path.stat()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Source changed during hashing: {path}")
    return {"path": str(path), "bytes": after.st_size, "sha256": digest}


def prepare() -> None:
    """Copy immutable test inputs and extract exactly one mdin heredoc."""
    source_records = []
    for name, (folder, run, stem, reference) in CASES.items():
        source = ROOT / "source" / name
        attempt = ROOT / "attempts" / name
        source.mkdir(parents=True, exist_ok=False)
        attempt.mkdir(parents=True, exist_ok=False)
        for filename in (run, stem + ".prmtop", stem + ".inpcrd", reference):
            original = SHIPPED / folder / filename
            source_records.append(record(original))
            shutil.copy2(original, source / filename)
        text = (source / run).read_text()
        marker = "cat > mdin <<EOF\n"
        if text.count(marker) != 1:
            raise ValueError(f"Expected exactly one mdin heredoc: {run}")
        _, _, tail = text.partition(marker)
        mdin, separator, _ = tail.partition("\nEOF\n")
        if not separator:
            raise ValueError(f"Missing mdin heredoc terminator: {run}")
        (source / "mdin").write_text(mdin + "\n")
        for filename in ("mdin", stem + ".prmtop", stem + ".inpcrd"):
            shutil.copy2(source / filename, attempt / filename)
    source_records.append(record(AMBERHOME / "bin/sander"))
    source_records.extend(record(AMBERHOME / "dat/slko/3ob-3-1" / f"{a}-{b}.skf")
                          for a in "CHNO" for b in "CHNO")
    source_records.append(record(Path(__file__)))
    (ROOT / "source_manifest.json").write_text(json.dumps(source_records, indent=2) + "\n")


def series(text: str, field: str) -> list[float]:
    """Read one named printed energy series, including summary records."""
    number = r"([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)"
    if field == "MIN_ENERGY":
        matches = re.findall(r"NSTEP\s+ENERGY[^\n]*\n\s*\d+\s+" + number, text)
    else:
        matches = re.findall(r"\b" + field + r"\s*=\s*" + number, text)
    values = [float(value.replace("D", "E").replace("d", "e")) for value in matches]
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError(f"Missing or nonfinite {field} energy series")
    return values


def compare() -> None:
    """Keep process completion distinct from regression energy agreement."""
    results = {}
    for name, (_, _, _, reference) in CASES.items():
        actual_path = ROOT / "attempts" / name / "engine.mdout"
        actual = actual_path.read_text()
        expected = (ROOT / "source" / name / reference).read_text()
        status = json.loads((ROOT / "jobs" / name / "status.json").read_text())
        fields = ("MIN_ENERGY", "DFTBESCF") if name == "water" else ("Etot", "EKtot", "EPtot", "DFTBESCF")
        comparisons = {}
        for field in fields:
            observed, saved = series(actual, field), series(expected, field)
            same_length = len(observed) == len(saved)
            differences = [a - b for a, b in zip(observed, saved, strict=False)]
            max_difference = max(abs(value) for value in differences)
            comparisons[field] = {
                "observed": observed, "reference": saved,
                "difference_kcal_per_mol": differences,
                "length_matches": same_length,
                "max_absolute_difference_kcal_per_mol": max_difference,
                "absolute_tolerance_kcal_per_mol": 0.0002,
                "matches": same_length and max_difference <= 0.0002 + 1e-10,
            }
        results[name] = {
            "process_status": status,
            "regression_energy_matches": all(row["matches"] for row in comparisons.values()),
            "energy_series": comparisons,
            "parameter_and_version_lines": [line for line in actual.splitlines()
                if any(token in line for token in ("Amber", "AMBER", ".skf", "Hubbard", "DFTB3", "dU/dq"))],
            "output_record": record(actual_path),
        }
    result = {
        "engine_regressions": results,
        "scope": "Shipped DFTB3 engine regressions only; no enzyme or mechanism validation.",
        "water_tolerance_note": "0.0002 kcal/mol comparison threshold chosen explicitly for printed minimization energies; NMA uses shipped Run tolerance.",
    }
    (ROOT / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")


def main() -> None:
    """Run explicit prepare or compare action."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "compare"))
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    else:
        compare()


if __name__ == "__main__":
    main()
