"""Regressions for sparse crystal templates across masked sequence gaps."""

from copy import deepcopy
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from nyl12_covalent.boltz_template_alignment import (
    bounded_template_rows,
    make_bounded_template_featurizer,
)


def _fixture() -> tuple[Any, list[Any], dict[int, int]]:
    query_chains, source_chains, query, source, records = [], [], [], [], []
    authors = {}
    for number, chain in enumerate("ABCDEFGH"):
        q_asym, t_asym = number + 20, number + 40
        query_chains.append({"name": chain, "asym_id": q_asym})
        source_chains.append({"name": chain, "asym_id": t_asym})
        length = 232 if chain in "ABCD" else 96
        for residue in range(length):
            query.append(
                {"asym_id": q_asym, "res_idx": residue, "token_idx": len(query)}
            )
        retained = (
            list(range(8, 95)) + list(range(104, 223))
            if chain in "ABCD"
            else list(range(1, 94))
        )
        for compact, author in enumerate(retained):
            token = {
                "asym_id": t_asym,
                "res_idx": compact,
                "token_idx": len(source),
                "source_coordinate": (float(author), float(number), 1.0),
            }
            authors[len(source)] = author
            source.append(token)
        segments = (
            [(7, 94, 0, 87), (103, 222, 87, 206)]
            if chain in "ABCD"
            else [(0, 93, 0, 93)]
        )
        for qs, qe, ts, te in segments:
            records.append(
                SimpleNamespace(
                    name="crystal",
                    force=False,
                    query_chain=chain,
                    template_chain=chain,
                    query_st=qs,
                    query_en=qe,
                    template_st=ts,
                    template_en=te,
                )
            )
    query_chains.append({"name": "I", "asym_id": 99})
    query.append({"asym_id": 99, "res_idx": 0, "token_idx": len(query)})
    data = SimpleNamespace(
        structure=SimpleNamespace(chains=query_chains),
        templates={"crystal": SimpleNamespace(chains=source_chains)},
        template_tokens={"crystal": source},
        tokens=query,
        record=SimpleNamespace(templates=records),
    )
    return data, records, authors


def test_sparse_segments_restore_all_author_indices_without_mutation() -> None:
    data, records, authors = _fixture()
    original = deepcopy(data)
    rows = bounded_template_rows(data, "crystal", records, 3)
    assert len(rows) == 1196
    assert len({r["q_idx"] for r in rows}) == 1196
    for row in rows:
        query = data.tokens[row["q_idx"]]
        assert query["res_idx"] + 1 == authors[row["token"]["token_idx"]]
        assert row["pdb_id"] == 3
        if query["asym_id"] < 24:
            assert query["res_idx"] + 1 not in range(95, 104)
            assert query["res_idx"] + 1 not in range(223, 233)
    assert not any(data.tokens[r["q_idx"]]["asym_id"] == 99 for r in rows)
    assert data == original


def test_original_constant_offset_leaks_loop_and_record_trim_cannot_fix_it() -> None:
    data, records, _ = _fixture()
    second = records[1]
    original_start = second.template_st - 1 - (second.query_st - 1)
    trimmed_start = second.template_st - second.query_st
    assert original_start == trimmed_start == -16
    leaked = {
        t["res_idx"] - original_start + 1
        for t in data.template_tokens["crystal"]
        if t["asym_id"] == 40
    }
    assert set(range(95, 104)) <= leaked
    bounded = bounded_template_rows(data, "crystal", records[:2], 0)
    actual = {data.tokens[r["q_idx"]]["res_idx"] + 1 for r in bounded}
    assert actual == set(range(8, 95)) | set(range(104, 223))


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("name", "wrong", "matching"),
        ("force", True, "nonforced"),
        ("query_chain", "I", "protein-only"),
        ("template_chain", "I", "protein-only"),
        ("query_st", -1, "Invalid"),
        ("query_en", 7, "Invalid"),
        ("template_st", -1, "Invalid"),
        ("template_en", 0, "Invalid"),
        ("query_en", 93, "equal-length"),
    ],
)
def test_invalid_records_fail(field: str, value: Any, message: str) -> None:
    data, records, _ = _fixture()
    record = deepcopy(records[0])
    setattr(record, field, value)
    with pytest.raises(ValueError, match=message):
        bounded_template_rows(data, "crystal", [record], 0)


def test_overlap_and_missing_source_or_query_fail() -> None:
    data, records, _ = _fixture()
    with pytest.raises(ValueError, match="twice"):
        bounded_template_rows(data, "crystal", [records[0], records[0]], 0)
    data.tokens = [q for q in data.tokens if q["token_idx"] != 7]
    with pytest.raises(ValueError, match="correspondence is absent"):
        bounded_template_rows(data, "crystal", records, 0)
    data, records, _ = _fixture()
    data.template_tokens["crystal"].pop(0)
    with pytest.raises(ValueError, match="lacks observed"):
        bounded_template_rows(data, "crystal", records, 0)


def test_empty_records_fail() -> None:
    data, _, _ = _fixture()
    with pytest.raises(ValueError, match="No explicit"):
        bounded_template_rows(data, "crystal", [], 0)


def test_wrapper_preserves_native_calculator_coordinates_and_grouping() -> None:
    data, _, _ = _fixture()
    calls = []

    def native(value: Any, rows: list[dict[str, Any]], maximum: int) -> dict:
        calls.append((value, rows, maximum))
        return {
            "coordinates": np.array([r["token"]["source_coordinate"] for r in rows])
        }

    torch = SimpleNamespace(tensor=np.asarray, stack=np.stack, float32=np.float32)
    featurize = make_bounded_template_featurizer(native, torch)
    result = featurize(data, 1400)
    assert len(calls) == 1
    assert calls[0][0] is data and calls[0][2] == 1400
    assert result["coordinates"].shape == (1, 1196, 3)
    assert np.array_equal(
        result["coordinates"][0],
        np.array([t["source_coordinate"] for t in data.template_tokens["crystal"]]),
    )
    assert result["template_force"].tolist() == [False]
    assert np.isinf(result["template_force_threshold"]).all()
    data.record.templates = []
    with pytest.raises(ValueError, match="required"):
        featurize(data, 1400)
