"""Native potential, independent coordinates and strict input checks."""

from pathlib import Path

import pytest

from nyl12_covalent.neb_umbrella_bias import (
    HarmonicWindow,
    amber_restraints,
    bias_energy,
)
from nyl12_covalent.neb_umbrella_inputs import (
    STAGE_STEPS,
    active_group_cards,
    audited_qmmm_block,
    build_mdin,
)

OLD = """old audited input
&cntrl
 imin=1,
/
&qmmm
 qm_theory='DFTB3', qmcharge=-2, dftb_telec=0.0,
 qmmask='@1-98', scfconv=1.d-10, qmmm_int=1,
/
&dlfind
 active=':1-20,25', neb='FROZEN',
/
"""
ROLES = {"carbon": 11, "leaving_n": 22, "thr_og": 33, "proton": 44}


@pytest.mark.parametrize(
    "delta,expected", [(-3, 16), (-2, 8), (-1, 2), (0, 0), (1, 2), (2, 8), (3, 16)]
)
def test_native_harmonic_and_linear_tails(delta: float, expected: float) -> None:
    window = HarmonicWindow(4, 5, 2, 3, outer_extent=2)
    assert bias_energy((4 + delta, 5), window) == pytest.approx(expected)


def test_independent_two_coordinate_energy() -> None:
    w = HarmonicWindow(1, -1, 2, 7)
    assert bias_energy((2, 1), w) == 30  # 2*1²+7*2²; no half factor.
    text = amber_restraints(w, ROLES)
    assert text.count("&rst") == 2
    assert "iat=11,22,11,33,0,0,0,0" in text
    assert "iat=33,44,22,44,0,0,0,0" in text
    assert text.count("rstwt=1,-1,0,0") == 2
    assert "rk2=7, rk3=7" in text


@pytest.mark.parametrize(
    "kwargs",
    [
        {"k_h": 0},
        {"k_acyl": -1},
        {"outer_extent": 0},
        {"q_h_center": float("nan")},
        {"q_acyl_center": float("inf")},
    ],
)
def test_invalid_window(kwargs: dict) -> None:
    values = dict(q_acyl_center=0, q_h_center=0, k_acyl=1, k_h=1)
    values.update(kwargs)
    with pytest.raises(ValueError):
        HarmonicWindow(**values)


@pytest.mark.parametrize("values", [(0,), (0, float("nan")), (0, 0, 0)])
def test_invalid_samples(values: tuple) -> None:
    with pytest.raises(ValueError):
        bias_energy(values, HarmonicWindow(0, 0, 1, 1))


@pytest.mark.parametrize("replacement", [0, -1, True, 22])
def test_invalid_role(replacement: int) -> None:
    roles = {**ROLES, "carbon": replacement}
    with pytest.raises(ValueError):
        amber_restraints(HarmonicWindow(0, 0, 1, 1), roles)


@pytest.mark.parametrize("stage", STAGE_STEPS)
def test_native_stage_preserves_method_and_boundary(stage: str) -> None:
    text = build_mdin(OLD, stage, 731111, ":1-20,25", Path("two_cv.RST"))
    assert "&dlfind" not in text
    assert "ibelly=1, bellymask=''" in text
    assert "RES 1 20\nRES 25 25\nEND\nEND" in text
    assert text.index("qmshake=0") < text.index("Exact old-band movable selection")
    assert text.index("RES 25 25") < text.index("&wt")
    assert "ntr=0" in text and "ntc=1, ntf=1" in text
    assert "qmshake=0" in text and "ig=731111" in text
    assert f"nstlim={STAGE_STEPS[stage]}" in text and "dt=0.0005000" in text
    assert ("irest=0, ntx=1" in text) == (stage == "heat")
    assert ("type='TEMP0'" in text) == (stage == "heat")
    assert "qmmask='@1-98', scfconv=1.d-10, qmmm_int=1," in text


@pytest.mark.parametrize(
    "key,replacement",
    [
        ("qmcharge=-2", "qmcharge=-1"),
        ("qm_theory='DFTB3'", "qm_theory='PM6'"),
        ("dftb_telec=0.0", "dftb_telec=300"),
        ("qmmask='@1-98'", "qmmask='@1-99'"),
        ("qmmask='@1-98'", "qmmask='@1-97,97'"),
        ("qmmask='@1-98'", "qmmask=':1'"),
        ("qmmask='@1-98'", "qmmask='@98-1'"),
        ("dftb_telec=0.0", "qmshake=1, dftb_telec=0.0"),
        ("dftb_telec=0.0", "dftb_telec=0.0, dftb_telec=0.0"),
    ],
)
def test_wrong_hamiltonian_rejected(key: str, replacement: str) -> None:
    with pytest.raises(ValueError):
        audited_qmmm_block(OLD.replace(key, replacement))


def test_existing_qmshake_and_single_atoms() -> None:
    text = OLD.replace("@1-98", "@1-97,98").replace(
        "qmmm_int=1,", "qmmm_int=1, qmshake=0,"
    )
    block = audited_qmmm_block(text)
    assert block.count("qmshake") == 1


@pytest.mark.parametrize(
    "changes",
    [
        {"stage": "bad"},
        {"seed": -1},
        {"seed": 2147483648},
        {"seed": True},
        {"active_mask": ":1-21"},
        {"restraint_path": Path("bad path")},
    ],
)
def test_bad_stage_configuration(changes: dict) -> None:
    args = dict(
        old_mdin=OLD,
        stage="heat",
        seed=1,
        active_mask=":1-20,25",
        restraint_path=Path("ok.RST"),
    )
    args.update(changes)
    with pytest.raises(ValueError):
        build_mdin(**args)


def test_missing_or_duplicate_block_fails() -> None:
    with pytest.raises(ValueError):
        audited_qmmm_block("")
    with pytest.raises(ValueError):
        audited_qmmm_block(OLD + OLD)


@pytest.mark.parametrize("mask", [":0", ":4-2", ":1-2-3", ":bad"])
def test_invalid_group_mask(mask: str) -> None:
    with pytest.raises(ValueError):
        active_group_cards(mask)


def test_atom_group_encoding() -> None:
    assert "ATOM 1 2\nATOM 7 7" in active_group_cards("@1-2,7")
