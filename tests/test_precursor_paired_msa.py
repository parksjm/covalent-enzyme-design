"""Biological row provenance and cleavage-boundary regressions for paired MSAs."""

import csv
import hashlib
from pathlib import Path

import pytest

from nyl12_covalent.precursor_paired_msa import (
    read_precursor_a3m,
    split_precursor_records,
    write_paired_precursor_msas,
)

QUERY = ("query", "ACDEFG")


def test_boundary_insertions_follow_beta_and_preserve_homolog_origin() -> None:
    records = [QUERY, ("same_label", "VtTKmnLMN"), ("same_label", "Y-QP-R")]
    result = split_precursor_records(records, alpha_query="ACD", beta_query="EFG")
    assert result["alpha_rows"] == [
        {"sequence": "ACD", "key": -1},
        {"sequence": "VtTK", "key": 1},
        {"sequence": "Y-Q", "key": 2},
    ]
    assert result["beta_rows"] == [
        {"sequence": "EFG", "key": -1},
        {"sequence": "mnLMN", "key": 1},
        {"sequence": "P-R", "key": 2},
    ]
    assert result["row_origins"] == [
        {"key": 1, "source_row": 1, "header": "same_label"},
        {"key": 2, "source_row": 2, "header": "same_label"},
    ]
    assert result["boundary_after_query_residue1"] == 3
    assert result["query_length"] == 6


def test_native_segment_duplicates_removed_from_both_sides_without_rekeying() -> None:
    records = [
        QUERY,
        ("a", "VTKLMN"),
        ("b", "V-TKRS"),
        ("c", "VT-KRS"),
        ("d", "HIKLMN"),
        ("e", "VTKQRS"),
        ("f", "---QRS"),
        ("g", "HIQ---"),
        ("h", "ACDLMN"),
        ("i", "QRSQ-W"),
    ]
    result = split_precursor_records(records, alpha_query="ACD", beta_query="EFG")
    assert [r["key"] for r in result["alpha_rows"]] == [-1, 1, 2, 9]
    assert [r["key"] for r in result["beta_rows"]] == [-1, 1, 2, 9]
    assert [r["source_row"] for r in result["excluded_rows"]] == list(range(3, 9))
    assert result["paired_homolog_row_count"] == 3


@pytest.mark.parametrize(
    ("records", "alpha", "beta"),
    [
        ([], "ACD", "EFG"),
        ([QUERY], "", "ACDEFG"),
        ([QUERY], "ACDEFG", ""),
        ([QUERY], "aCD", "EFG"),
        ([QUERY], "AC-", "EFG"),
        ([("query", "AAAAAA")], "ACD", "EFG"),
    ],
)
def test_exact_mature_queries_required(
    records: list[tuple[str, str]], alpha: str, beta: str
) -> None:
    with pytest.raises(ValueError, match="exactly equal"):
        split_precursor_records(records, alpha_query=alpha, beta_query=beta)


@pytest.mark.parametrize(
    ("row", "message"),
    [
        (("", "VTKLMN"), "requires a header"),
        (("empty", ""), "requires a header"),
        (("short", "VTKLM"), "query columns"),
        (("long", "VTKLMNP"), "query columns"),
        (("bad", "VTK?MN"), "Unsupported"),
    ],
)
def test_bad_homolog_rows_fail(row: tuple[str, str], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        split_precursor_records([QUERY, row], alpha_query="ACD", beta_query="EFG")


def test_query_only_and_duplicate_only_sources_fail() -> None:
    for records in [[QUERY], [QUERY, ("copy", QUERY[1])]]:
        with pytest.raises(ValueError, match="No independently retained"):
            split_precursor_records(records, alpha_query="ACD", beta_query="EFG")


def test_reader_handles_multiline_and_validated_single_precursor_metadata(
    tmp_path: Path,
) -> None:
    path = tmp_path / "source.a3m"
    path.write_text("#6\t1\n\n>query\nACD\nEFG\n>origin accession\nVTK\nmnLMN\n")
    assert read_precursor_a3m(path) == [QUERY, ("origin accession", "VTKmnLMN")]


@pytest.mark.parametrize(
    "text",
    [
        "",
        "ACDEFG\n>query\nACDEFG\n",
        ">\nACDEFG\n",
        "#3,3\t1,1\n>query\nACDEFG\n",
        "#6\t2\n>query\nACDEFG\n",
        "#7\t1\n>query\nACDEFG\n",
        "#6\t1\n#6\t1\n>query\nACDEFG\n",
        ">query\nACDEFG\n#6\t1\n",
        "#bad\n>query\nACDEFG\n",
    ],
)
def test_reader_rejects_ambiguous_or_invalid_sources(tmp_path: Path, text: str) -> None:
    path = tmp_path / "bad.a3m"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_precursor_a3m(path)


def test_writer_preserves_source_and_hashes_real_native_csv_rows(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.a3m"
    original = b">query\nACDEFG\n>same\nVTKmnLMN\n>same\nY-QP-R\n"
    source.write_bytes(original)
    output = tmp_path / "paired"
    receipt = write_paired_precursor_msas(
        source, output, alpha_query="ACD", beta_query="EFG"
    )
    assert source.read_bytes() == original
    assert receipt["source_sha256"] == hashlib.sha256(original).hexdigest()
    for family, expected in [
        ("alpha", ["ACD", "VTK", "Y-Q"]),
        ("beta", ["EFG", "mnLMN", "P-R"]),
    ]:
        path = output / f"{family}_paired.csv"
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        assert [r["sequence"] for r in rows] == expected
        assert [int(r["key"]) for r in rows] == [-1, 1, 2]
        assert (
            receipt["output_sha256"][path.name]
            == hashlib.sha256(path.read_bytes()).hexdigest()
        )
    assert (output / "precursor_pairing_receipt.json").is_file()
    with pytest.raises(FileExistsError):
        write_paired_precursor_msas(source, output, alpha_query="ACD", beta_query="EFG")


def test_invalid_alignment_does_not_create_output_directory(tmp_path: Path) -> None:
    source = tmp_path / "bad.a3m"
    source.write_text(">query\nACDEFG\n")
    output = tmp_path / "should_not_exist"
    with pytest.raises(ValueError):
        write_paired_precursor_msas(source, output, alpha_query="ACD", beta_query="EFG")
    assert not output.exists()
