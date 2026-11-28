"""Identity, precision and graph regressions for the explicit export adapter."""

import numpy as np
import pytest

from nyl12_covalent.rfd3_export_adapter import (
    adapt_coordinates,
    bf16_coordinates,
    reassociate_typed_graph,
)


def fixture() -> tuple:
    """Return a complete terminal inventory and one deliberately mobile atom."""
    keys = [
        (c, 232 if c in "ABCD" else 96, n)
        for c in "ABCDEFGH"
        for n in ("C", "O", "OXT")
    ]
    keys += [("A", 100, "CA")]
    source = np.random.default_rng(8).normal(size=(len(keys), 3)) + 60
    native_keys = [key for key in keys if key[2] != "OXT"]
    center = np.array([60.125, 60.25, 60.5], dtype=np.float32)
    order = [keys.index(key) for key in native_keys]
    native = bf16_coordinates(source.astype(np.float32)[order] - center)
    native[-1] += 0.123
    return keys, source, native_keys, native, {"A100"}, center


def test_restores_only_fixed_and_preserves_mobile() -> None:
    args = fixture()
    raw = args[3].copy()
    result, evidence = adapt_coordinates(*args)
    np.testing.assert_array_equal(result[:-1], args[1][:-1])
    np.testing.assert_array_equal(result[-1], raw[-1].astype(float) + args[-1])
    np.testing.assert_array_equal(args[3], raw)
    assert evidence["accepted"]
    assert evidence["fixed_components_compared"] == 48
    assert len(evidence["restored_source_OXT"]) == 8


def test_invalid_inventory_precision_masks_and_coordinates_fail() -> None:
    for change in (
        "duplicate",
        "missing",
        "shape",
        "nan",
        "precision",
        "mobile_oxt",
        "allfixed",
    ):
        args = list(fixture())
        if change == "duplicate":
            args[2][-1] = args[2][0]
        elif change == "missing":
            args[2] = args[2][1:]
        elif change == "shape":
            args[3] = args[3][1:]
        elif change == "nan":
            args[3][0, 0] = np.nan
        elif change == "precision":
            args[3][0, 0] += 0.01
        elif change == "mobile_oxt":
            args[4].add("A232")
        else:
            args[4] = set()
        with pytest.raises(ValueError):
            adapt_coordinates(*args)
    with pytest.raises(ValueError):
        bf16_coordinates(np.array([np.nan]))


def test_typed_graph_is_reindexed_without_inference() -> None:
    keys = [("I", i, "C") for i in (1, 2, 3)]
    bonds = np.array([[0, 1, 1], [1, 2, 2]])
    np.testing.assert_array_equal(
        reassociate_typed_graph(keys, bonds, keys[::-1]), [[2, 1, 1], [1, 0, 2]]
    )
    for target, graph in (
        (keys[:-1], bonds),
        (keys + keys[:1], bonds),
        (keys, np.array([])),
        (keys, np.array([[0, 4, 1]])),
        (keys, np.array([[0, 1, 0]])),
    ):
        with pytest.raises(ValueError):
            reassociate_typed_graph(keys, graph, target)
