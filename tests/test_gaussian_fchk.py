"""Tests for strict Gaussian formatted-checkpoint parsing."""

import numpy as np
import pytest

from nyl12_covalent.gaussian_fchk import parse_gaussian_fchk


def _record(label: str, kind: str, value: str) -> str:
    return f"{label:<40}   {kind}   {value}"


def _array(label: str, kind: str, values: list[str], width: int = 5) -> str:
    rows = [
        " ".join(values[start : start + width])
        for start in range(0, len(values), width)
    ]
    return "\n".join([_record(label, kind, f"N= {len(values)}"), *rows])


def _character_array(label: str, values: list[str], width: int = 5) -> str:
    rows = [
        "".join(f"{value:<12}" for value in values[start : start + width])
        for start in range(0, len(values), width)
    ]
    return "\n".join([_record(label, "C", f"N= {len(values)}"), *rows])


def _fixture() -> str:
    # Lower-triangle row-major values for a 6x6 matrix (two atoms).
    packed = [str(value) for value in range(1, 22)]
    return "\n".join(
        [
            "Title Card Required",
            "SP RHF/6-31G(d) FChk",
            _record("Number of atoms", "I", "2"),
            _record("Charge", "I", "-1"),
            _record("Multiplicity", "I", "1"),
            _array("Atomic numbers", "I", ["6", "8"]),
            _array(
                "Current cartesian coordinates",
                "R",
                ["0.0", "1.0D+00", "-2.0", "3.0", "4.0", "5.0"],
                width=4,
            ),
            _array("Cartesian Force Constants", "R", packed, width=7),
        ]
    )


def test_parses_native_units_and_packed_lower_triangle() -> None:
    result = parse_gaussian_fchk(
        _fixture(),
        expected_atomic_numbers=(6, 8),
        expected_charge=-1,
        expected_multiplicity=1,
    )

    assert result.atomic_numbers.tolist() == [6, 8]
    assert result.coordinates_bohr.shape == (2, 3)
    assert result.coordinates_bohr[0].tolist() == [0.0, 1.0, -2.0]
    assert result.charge == -1
    assert result.multiplicity == 1
    assert result.hessian_eh_per_bohr2.shape == (6, 6)
    assert result.hessian_eh_per_bohr2[0, 0] == 1.0
    assert result.hessian_eh_per_bohr2[1, 0] == 2.0
    assert result.hessian_eh_per_bohr2[1, 1] == 3.0
    assert result.hessian_eh_per_bohr2[5, 5] == 21.0
    np.testing.assert_array_equal(
        result.hessian_eh_per_bohr2, result.hessian_eh_per_bohr2.T
    )


@pytest.mark.parametrize(
    ("expected", "message"),
    [
        ((6, 7), "atomic numbers"),
        ((6, 8), "charge"),
    ],
)
def test_rejects_identity_mismatch(expected: tuple[int, ...], message: str) -> None:
    kwargs = (
        {"expected_atomic_numbers": expected}
        if message == "atomic numbers"
        else {"expected_charge": 0}
    )
    with pytest.raises(ValueError, match=message):
        parse_gaussian_fchk(_fixture(), **kwargs)


def test_rejects_truncated_array() -> None:
    text = _fixture().replace(
        _record("Current cartesian coordinates", "R", "N= 6"),
        _record("Current cartesian coordinates", "R", "N= 7"),
    )
    with pytest.raises(ValueError, match="truncated array"):
        parse_gaussian_fchk(text)


def test_rejects_duplicate_required_record() -> None:
    charge_record = _record("Charge", "I", "-1")
    text = _fixture().replace(charge_record, f"{charge_record}\n{charge_record}")
    with pytest.raises(ValueError, match="duplicate"):
        parse_gaussian_fchk(text)


def test_rejects_incorrect_hessian_triangle_length() -> None:
    full_packed = _array(
        "Cartesian Force Constants", "R", [str(i) for i in range(1, 22)], 7
    )
    short_packed = _array(
        "Cartesian Force Constants", "R", [str(i) for i in range(1, 21)], 7
    )
    text = _fixture().replace(full_packed, short_packed)
    with pytest.raises(ValueError, match="3N triangle"):
        parse_gaussian_fchk(text)


def test_rejects_nonfinite_real() -> None:
    with pytest.raises(ValueError, match="malformed real"):
        parse_gaussian_fchk(_fixture().replace("1.0D+00", "NaN"))


def test_skips_fixed_width_character_arrays_and_unrelated_duplicate_fields() -> None:
    text = "\n".join(
        [
            _fixture(),
            _character_array(
                "Full Title",
                ["neutral beta", "-N charge-", "template;", "not enzyme", "refinement"],
            ),
            _character_array("Route", ["#p hf/6-31g", " opt test"]),
            _character_array("Atom Types", ["", "", "", "", "", ""], width=5),
            _record("Force Field", "I", "0"),
            _record("Force Field", "I", "1"),
            _character_array("Gaussian Version", ["ES64L-G16", "RevC.02"]),
        ]
    )

    result = parse_gaussian_fchk(text)

    assert result.atomic_numbers.tolist() == [6, 8]
    assert result.hessian_eh_per_bohr2.shape == (6, 6)
