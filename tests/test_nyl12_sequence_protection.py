"""Hard maturation/catalytic protections survive stale or permissive masks."""

import pytest

from nyl12_covalent.ligandmpnn_pilot_audit import _audit_ties_and_changes
from nyl12_covalent.nyl12_sequence_protection import (
    MINIMUM_PROTECTED_ALPHA,
    MINIMUM_PROTECTED_BETA_LOCAL,
    MINIMUM_PROTECTED_PRECURSOR,
)
from nyl12_covalent.relay_broad_pilot import CORE_PROTECTED
from nyl12_covalent.relay_ensemble_sequences import audit_sequence


def _native() -> dict[str, str]:
    return {**dict.fromkeys("ABCD", "A" * 232), **dict.fromkeys("EFGH", "A" * 96)}


def test_hard_set_covers_full_tail_and_recovered_legacy_rules() -> None:
    required = {113, 156, 186, *range(221, 236), *range(272, 277)}
    assert MINIMUM_PROTECTED_PRECURSOR == required
    assert MINIMUM_PROTECTED_ALPHA <= CORE_PROTECTED["alpha"]
    assert MINIMUM_PROTECTED_BETA_LOCAL <= CORE_PROTECTED["beta"]
    assert {1, 2, 3, 40, 41, 42, 43, 44} == MINIMUM_PROTECTED_BETA_LOCAL


@pytest.mark.parametrize("position", sorted(MINIMUM_PROTECTED_PRECURSOR))
def test_every_required_position_rejected_even_with_permissive_mask(
    position: int,
) -> None:
    native = _native()
    family = "alpha" if position <= 232 else "beta"
    local = position if position <= 232 else position - 232
    changed = dict(native)
    for chain in "ABCD" if family == "alpha" else "EFGH":
        value = changed[chain]
        changed[chain] = value[: local - 1] + "C" + value[local:]
    mask = {"alpha": list(range(1, 233)), "beta": list(range(1, 97))}
    with pytest.raises(ValueError, match="protected"):
        audit_sequence(":".join(changed.values()), native, mask)
    with pytest.raises(ValueError, match="protected"):
        _audit_ties_and_changes(
            changed,
            native,
            set(mask["alpha"]),
            set(mask["beta"]),
            set(),
            set(),
            "stale_protection_config",
            1,
        )


def test_other_mutations_remain_allowed_and_input_sets_are_not_mutated() -> None:
    native = _native()
    changed = dict(native)
    for chain in "ABCD":
        changed[chain] = "AC" + native[chain][2:]
    alpha, beta = set(), set()
    changes, _, _ = _audit_ties_and_changes(
        changed, native, {2}, set(), alpha, beta, "valid", 1
    )
    assert changes == [
        {
            "family": "alpha",
            "physical_position": 2,
            "merged_position": 2,
            "from": "A",
            "to": "C",
        }
    ]
    assert alpha == beta == set()


def test_k156_remains_protected_even_when_in_a_stale_design_mask() -> None:
    native = _native()
    changed = dict(native)
    for chain in "ABCD":
        value = native[chain]
        native[chain] = value[:155] + "K" + value[156:]
        changed[chain] = value[:155] + "E" + value[156:]
    mask = {"alpha": [156], "beta": []}
    assert 156 in MINIMUM_PROTECTED_PRECURSOR
    assert 156 in CORE_PROTECTED["alpha"]
    with pytest.raises(ValueError, match="protected"):
        audit_sequence(":".join(changed.values()), native, mask)
    with pytest.raises(ValueError, match="protected"):
        _audit_ties_and_changes(
            changed, native, {156}, set(), set(), set(), "K156E_rejected", 1
        )
