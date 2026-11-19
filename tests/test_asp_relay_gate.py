"""Tests for the explicit Asp272 relay-preorganization gate."""

import pytest

from nyl12_covalent.asp_relay_gate import relay_contact_gate


def test_relay_gate_accepts_predeclared_contact_window() -> None:
    result = relay_contact_gate({"OG1_Asp272_OD1_A": 2.98, "HG1_Asp272_OD1_A": 2.18})
    assert result["accepted"] is True
    assert "native Thr233" in result["interpretation"]


@pytest.mark.parametrize(
    ("og1", "hg1"),
    [(2.49, 1.9), (3.21, 1.9), (2.9, 1.49), (2.9, 2.31)],
)
def test_relay_gate_rejects_out_of_range_contact(og1: float, hg1: float) -> None:
    result = relay_contact_gate({"OG1_Asp272_OD1_A": og1, "HG1_Asp272_OD1_A": hg1})
    assert result["accepted"] is False
