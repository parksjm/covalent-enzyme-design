"""Verify insertion handling, conservation denominators and residue mapping."""

from pathlib import Path
import json
import sys

import pytest

from scripts.analyze_nyl12_site_conservation import (
    aligned_sequence,
    analyze_msa,
    global_alignment,
    main,
    site_counts,
)


def test_lowercase_insertions_preserve_match_columns() -> None:
    assert aligned_sequence("ATaaD-GtD") == "ATD-GD"
    for invalid in ("", "abcd", "AC?D", "AC.D"):
        with pytest.raises(ValueError):
            aligned_sequence(invalid)


def test_gap_and_ambiguous_denominators() -> None:
    result = site_counts(["D", "D", "N", "-", "X"], 1, "D")
    assert result["rows"] == 5
    assert result["nongap"] == 4
    assert result["query_fraction_nongap"] == 0.5
    assert result["query_fraction_standard"] == pytest.approx(2 / 3)
    assert result["gap_fraction"] == 0.2
    assert site_counts(["-"], 1, "D")["query_fraction_nongap"] is None
    with pytest.raises(ValueError):
        site_counts([], 1, "D")


def test_deduplication_and_coverage_filter(tmp_path: Path) -> None:
    query = "T" * 39 + "DGD" + "A" * 54
    alternative = "S" + query[1:]
    gapped = "-" * 30 + query[30:]
    path = tmp_path / "msa.csv"
    path.write_text(
        f"key,sequence\n0,{query}\n-1,{query}\n-1,{alternative}\n-1,{alternative[:10]}aa{alternative[10:]}\n-1,{gapped}\n"
    )
    result = analyze_msa(path, "beta", 0.8)
    assert result["raw_rows"] == 5
    assert result["unique_rows_including_query"] == 3
    assert result["accepted_unique_nonquery"] == 1
    assert result["rejected_coverage"] == 1
    assert (
        result["sites"]["272"]["filtered_unique_nonquery"]["query_fraction_nongap"] == 1
    )
    path.write_text("key,sequence\n0,AAA\n-1,AA\n")
    with pytest.raises(ValueError, match="lengths"):
        analyze_msa(path, "beta", 0.8)


def test_global_alignment_maps_insertions_and_deletions() -> None:
    a, b, mapping = global_alignment("ACDEFG", "ACQDEFG")
    assert (a, b) == ("AC-DEFG", "ACQDEFG")
    assert mapping[3] == 4
    a, b, mapping = global_alignment("ACQDEFG", "ACDEFG")
    assert (a, b) == ("ACQDEFG", "AC-DEFG")
    assert 3 not in mapping


def test_cli_retains_hashed_sources_and_refuses_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    alpha_query = "A" * 112 + "F" + "A" * 42 + "K" + "A" * 29 + "N" + "A" * 46
    beta_query = "T" * 39 + "DGD" + "A" * 54
    alpha, beta, fasta = (
        tmp_path / name for name in ("alpha.csv", "beta.csv", "ref.fasta")
    )
    for path, query in ((alpha, alpha_query), (beta, beta_query)):
        path.write_text(f"key,sequence\n0,{query}\n-1,S{query[1:]}\n")
    fasta.write_text(f">public_test_reference\n{alpha_query}{beta_query}\n")
    output = tmp_path / "audit"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "analyze",
            "--alpha-msa",
            str(alpha),
            "--beta-msa",
            str(beta),
            "--nylc-fasta",
            str(fasta),
            "--nylc-origin-url",
            "https://example.org/public.fasta",
            "--output",
            str(output),
        ],
    )
    main()
    report = json.loads((output / "conservation.json").read_text())
    assert len(report["reference"]["sha256"]) == 64
    assert report["correspondence"]["272"]["NylC_position"] == 272
    assert (output / "public_nylc_3axg.fasta").read_text() == fasta.read_text()
    assert (
        "Nyl12_modeled_alpha_plus_beta"
        in (output / "nyl12_nylc_numbering_alignment.fasta").read_text()
    )
    with pytest.raises(FileExistsError):
        main()
