"""Scientific regression checks for direct acyl-enzyme endpoint identity."""

from nyl12_covalent.direct_acyl_drive import (
    candidate_geometry,
    chemistry_retained,
    composite,
    direct_restraint,
)


def _geometry() -> dict[str, float]:
    return {
        "retained_QM_H_parent_min_A": 1.0,
        "retained_QM_H_parent_max_A": 1.1,
        "N_H1_A": 1.02,
        "N_H2_A": 1.02,
        "OG1_H1_A": 2.0,
        "OG1_H2_A": 2.0,
        "substrate_N_HN_A": 1.02,
        "N_HG1_A": 2.0,
        "substrate_N_HG1_A": 1.02,
        "OG1_HG1_A": 1.8,
        "OG1_carbonyl_C_A": 1.5,
        "carbonyl_C_O_A": 1.25,
        "carbonyl_C_N_A": 2.5,
        "proton_difference_A": 0.78,
    }


def test_acyl_candidate_requires_cleavage_and_carbonyl() -> None:
    values = _geometry()
    assert candidate_geometry(values)
    values.update(carbonyl_C_N_A=1.4, carbonyl_C_O_A=1.4)
    assert chemistry_retained(values)
    assert not candidate_geometry(values)


def test_proton_destination_and_retained_hydrogens() -> None:
    values = _geometry()
    values["N_HG1_A"] = 1.1
    assert not chemistry_retained(values)
    values["N_HG1_A"] = 2.0
    values["substrate_N_HN_A"] = 1.6
    assert not chemistry_retained(values)
    values["substrate_N_HN_A"] = 1.02
    values["retained_QM_H_parent_max_A"] = 1.6
    assert not chemistry_retained(values)


def test_shared_proton_and_restraint_mapping() -> None:
    values = _geometry()
    values.update(substrate_N_HG1_A=1.5, OG1_HG1_A=1.5)
    assert chemistry_retained(values)
    assert not candidate_geometry(values)
    assert abs(composite(values) - 1.78) < 1e-10
    text = direct_restraint(0.3)
    assert "iat=13916,13917,13917,19929,13916,19916,19916,19929" in text
    assert "rstwt=1,-1,-1,1" in text
    assert "13905" not in text


def test_optional_proton_guard_preserves_default_coordinate() -> None:
    original = direct_restraint(0.3)
    guarded = direct_restraint(0.3, substrate_proton_guard=True)
    assert guarded.startswith(original)
    assert guarded.count("&rst") == 2
    assert "iat=19929,13917, r1=0.8, r2=0.95, r3=1.25, r4=1.5" in guarded
