import pytest

from nyl12_covalent.es_dftb_unrestrained import compact_mask, dftb_input


def test_compact_mask_merges_only_consecutive_indices() -> None:
    assert compact_mask([1, 2, 3, 6, 8, 9]) == "@1-3,6,8-9"


def test_compact_mask_rejects_duplicates() -> None:
    with pytest.raises(ValueError, match="unique"):
        compact_mask([1, 1])


def test_dftb_input_is_fully_unrestrained_and_saves_trajectory() -> None:
    text = dftb_input(
        "@1-3",
        active_mask="(:1-1313)|(:929,1313<:8.0)",
        case_label="Template Model 7",
    )
    assert text.startswith("Template Model 7")
    assert "ntr=0" in text
    assert "nmropt=0" in text
    assert "active='(:1-1313)|(:929,1313<:8.0)'" in text
    assert "ntwx=25" in text
    assert "restraint_wt" not in text
    assert "DISANG" not in text


def test_dftb_input_rejects_empty_active_mask() -> None:
    with pytest.raises(ValueError, match="Amber atom mask"):
        dftb_input("@1-3", active_mask="")
