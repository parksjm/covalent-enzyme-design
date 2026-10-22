"""Tests for optional Gaussian formatted-checkpoint gradients."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.gaussian_fchk import parse_gaussian_fchk


def _record(label: str, kind: str, value: str) -> str:
    """Format a scalar formatted-checkpoint record."""
    return f"{label:<40}   {kind}   {value}"


def _fchk(gradient: str | None = None) -> str:
    """Build a minimal valid two-atom formatted-checkpoint fixture."""
    lines = [
        "minimal fixture",
        "generated for parser test",
        _record("Number of atoms", "I", "2"),
        _record("Charge", "I", "0"),
        _record("Multiplicity", "I", "1"),
        _record("Atomic numbers", "I", "N= 2"),
        " 6 8",
        _record("Current cartesian coordinates", "R", "N= 6"),
        " 0.0 1.0 2.0 3.0 4.0 5.0",
    ]
    if gradient is not None:
        lines.extend(
            [
                _record("Cartesian Gradient", "R", f"N= {len(gradient.split())}"),
                gradient,
            ]
        )
    lines.extend(
        [
            _record("Cartesian Force Constants", "R", "N= 21"),
            " " + " ".join(str(value) for value in range(21)),
        ]
    )
    return "\n".join(lines) + "\n"


def test_gradient_parses_d_exponents_with_shape_and_units() -> None:
    result = parse_gaussian_fchk(_fchk("1.0D-3 -2.0D-3 3.0D-3 4.0D-3 5.0D-3 -6.0D-3"))

    assert result.gradient_eh_per_bohr is not None
    assert result.gradient_eh_per_bohr.shape == (2, 3)
    np.testing.assert_allclose(
        result.gradient_eh_per_bohr,
        [[0.001, -0.002, 0.003], [0.004, 0.005, -0.006]],
    )


def test_absent_gradient_is_none() -> None:
    result = parse_gaussian_fchk(_fchk())

    assert result.gradient_eh_per_bohr is None


@pytest.mark.parametrize(
    ("kind", "values", "count", "message"),
    [
        ("R", "1 2 3 4 5", 5, "expected 6"),
        ("I", "1 2 3 4 5 6", 6, "real-valued"),
        ("R", "1 2 3 4 5 NaN", 6, "malformed real"),
        ("R", "1 2 3 4 5 1e999", 6, "non-finite"),
    ],
)
def test_malformed_gradient_is_rejected(
    kind: str, values: str, count: int, message: str
) -> None:
    fixture = _fchk()
    gradient_record = (
        _record("Cartesian Gradient", kind, f"N= {count}") + "\n" + values + "\n"
    )
    fixture = fixture.replace(
        _record("Cartesian Force Constants", "R", "N= 21"),
        gradient_record + _record("Cartesian Force Constants", "R", "N= 21"),
    )

    with pytest.raises(ValueError, match=message):
        parse_gaussian_fchk(fixture)


def test_duplicate_gradient_is_rejected() -> None:
    record = _record("Cartesian Gradient", "R", "N= 6") + "\n1 2 3 4 5 6\n"
    fixture = _fchk("1 2 3 4 5 6")
    fixture = fixture.replace(record, record + record)

    with pytest.raises(ValueError, match="duplicate"):
        parse_gaussian_fchk(fixture)
