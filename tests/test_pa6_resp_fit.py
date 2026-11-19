"""Tests for strict shared multi-conformer RESP charge handling."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent.pa6_resp_fit import parse_shared_resp_charges


def _write_qout(path: Path, values: np.ndarray) -> None:
    path.write_text(
        "\n".join(" ".join(f"{value:.8f}" for value in row) for row in values)
    )


def test_shared_resp_charges_accept_identical_neutral_copies(tmp_path: Path) -> None:
    row = np.asarray([0.2, -0.3, 0.1])
    qout = tmp_path / "charges.qout"
    _write_qout(qout, np.tile(row, (4, 1)))

    observed = parse_shared_resp_charges(qout, atom_count=3, conformer_count=4)

    assert np.allclose(observed, row)


@pytest.mark.parametrize("failure", ["count", "difference", "charge", "nonfinite"])
def test_shared_resp_charges_reject_invalid_fit(tmp_path: Path, failure: str) -> None:
    matrix = np.tile(np.asarray([0.2, -0.3, 0.1]), (4, 1))
    if failure == "count":
        matrix = matrix.ravel()[:-1].reshape(1, -1)
    elif failure == "difference":
        matrix[2, 1] += 0.01
    elif failure == "charge":
        matrix[:, 0] += 0.01
    else:
        matrix[1, 0] = np.nan
    qout = tmp_path / "charges.qout"
    _write_qout(qout, matrix)

    with pytest.raises(ValueError):
        parse_shared_resp_charges(qout, atom_count=3, conformer_count=4)
