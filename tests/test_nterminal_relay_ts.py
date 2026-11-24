"""Verify conditional N-terminal relay TS preparation without running science."""

import sys
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import nterminal_relay_ts as ts


SOURCE = """Title
&cntrl
 ntr=0, nmropt=0,
/
&qmmm
 dftb_telec=0.0, scfconv=1.0d-12, qmmask='@1-2',
/
&dlfind
 wtmask='@1-2', active='@1-2',
/
"""


def rows(energies):
    return [
        {"label": str(i), "accepted": True, "physical_kcal_mol": e, "restart": str(i)}
        for i, e in enumerate(energies)
    ]


def test_peak_is_highest_driven_and_requires_strict_neighbors():
    assert ts._peak(rows([0, 3, 1, 4, 2]))[1]["label"] == "3"
    assert ts._peak(rows([0, 3, 1, 4, 5])) is None
    assert ts._peak(rows([0, 2, 2, 0])) is None
    with pytest.raises(ValueError):
        ts._peak(rows([0, float("nan"), 1]))


def test_mode_preserves_exterior_and_norm():
    active = np.array([True, True, False])
    ref = np.zeros((3, 3))
    before = ref.copy()
    peak = ref.copy()
    after = ref.copy()
    after[0, 0] = 2
    mode = ts._mode(before, peak, after, active, ref)
    assert np.linalg.norm(mode[active] - peak[active]) == pytest.approx(0.01)
    assert np.array_equal(mode[~active], ref[~active])
    after[-1, 0] = 1
    with pytest.raises(ValueError, match="frozen"):
        ts._mode(before, peak, after, active, ref)
    with pytest.raises(ValueError, match="Degenerate"):
        ts._mode(ref, ref, ref, active, ref)


def test_dimer_input_preserves_control_and_masks(tmp_path):
    text = ts._dimer_input(SOURCE, tmp_path / "ref", tmp_path / "mode")
    assert "dimer=.true." in text and "tole=1.d-7" in text and "scfconv=1.0d-12" in text
    assert "wtmask='@1-2', active='@1-2'" in text
    for bad in (
        SOURCE.replace("nmropt=0", "nmropt=1"),
        SOURCE.replace("wtmask='@1-2'", "wtmask='@3'"),
        SOURCE.replace("&dlfind", "&wt type='END' /\n&dlfind"),
    ):
        with pytest.raises(ValueError):
            ts._dimer_input(bad, tmp_path / "ref", tmp_path / "mode")


def test_no_peak_records_no_launch(tmp_path, monkeypatch):
    monkeypatch.setattr(ts, "journal", lambda *args: None)
    result = ts.run_ts(tmp_path, tmp_path, rows([0, 1, 2]))
    assert result["status"] == "not_launched"
    assert (tmp_path / "03_ts_search/decision.json").is_file()


@pytest.mark.parametrize(
    "scenario", ["success", "native_fail", "force_fail", "exterior"]
)
def test_native_dispatch_and_independent_force(tmp_path, monkeypatch, scenario):
    monkeypatch.setattr(ts, "journal", lambda *args: None)
    model = tmp_path / "model/inputs"
    model.mkdir(parents=True)
    reference = model / "optimized.ncrst"
    reference.write_text("ref")
    run_rows = rows([0, 2, 1])
    for i, row in enumerate(run_rows):
        path = tmp_path / f"{i}.ncrst"
        path.write_text(str(i))
        row["restart"] = str(path)
    candidate = tmp_path / "candidate.ncrst"
    candidate.write_text("candidate")
    active = np.array([True, True, False])

    class State:
        def __init__(self, root, restart):
            self.active = active
            self.xyz = np.zeros((3, 3))
            if str(restart).endswith("2.ncrst"):
                self.xyz[0, 0] = 1
            if str(restart).endswith("candidate.ncrst") and scenario == "exterior":
                self.xyz[2, 0] = 1
            self.restart = object()
            self.mdin_source = SOURCE

        def gradient(self, xyz, evidence, *, precision):
            assert precision
            return -1.0, np.full((2, 3), 1e-3 if scenario == "force_fail" else 1e-5)

    monkeypatch.setattr(ts, "NativeState", State)
    monkeypatch.setattr(
        ts, "write_restart", lambda path, xyz, ref: path.write_text("mode")
    )

    def run_native(campaign, folder, start, mdin, restraints):
        assert (
            restraints is None
            and folder.name == "native_dimer"
            and "dimer=.true." in mdin
        )
        return {
            "accepted": scenario != "native_fail",
            "native_converged": True,
            "restart": str(candidate),
        }

    monkeypatch.setitem(
        sys.modules,
        "nyl12_covalent.nterminal_relay_native",
        SimpleNamespace(run_native=run_native),
    )
    if scenario == "success":
        result = ts.run_ts(tmp_path, tmp_path, run_rows)
        assert result["status"] == "awaiting_curvature_verification"
        assert result["maximum_gradient_Eh_per_bohr"] < 1e-4
    else:
        with pytest.raises(ValueError):
            ts.run_ts(tmp_path, tmp_path, run_rows)
