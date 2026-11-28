"""Validate independent pooling, native parity and paired-block uncertainty gates."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent import relay_pooled_analysis as analysis
from nyl12_covalent import relay_pooled_estimator as estimator
from nyl12_covalent import relay_pooled_sources as sources
from nyl12_covalent.relay_pooled_parity import validate_parity
from nyl12_covalent.relay_umbrella_analysis_io import load_block
from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow


def test_real_native_parity_and_failures() -> None:
    """Both independently emitted CVs and biasenergies must pass, includingmask."""
    cfg = json.loads(
        (
            sources.FIXED
            / "automatic_progression_20261126/checkpoint_45ps/analysis_request.json"
        ).read_text()
    )
    path = Path(cfg["accepted_block_receipts"][0])
    block = load_block(
        path,
        temperature=300,
        hamiltonian=cfg["physical_hamiltonian_sha256"],
        diagnostic_columns=cfg["diagnostic_columns"],
    )
    bias = RelayBiasWindow(**cfg["windows"][0]["bias"])
    proof = validate_parity(block, bias)
    assert proof["evaluated_frames"] == 99
    bad = dict(block, native_bias_kcal_mol=block["native_bias_kcal_mol"] + 10)
    with pytest.raises(ValueError, match="biasparity"):
        validate_parity(bad, bias)
    bad = dict(block, native_cv_values_A=block["native_cv_values_A"] + 1)
    with pytest.raises(ValueError, match="CVparity"):
        validate_parity(bad, bias)
    bad = dict(block, native_bias_evaluated=np.zeros(100, bool))
    with pytest.raises(ValueError, match="terminalframe"):
        validate_parity(bad, bias)


def test_exact_source_prefix_no_derived_discovery(tmp_path: Path) -> None:
    """Ignore equilandderivedreceipts; requireeverycompleteproductionblock."""
    assert sources.checkpoint_receipts(tmp_path, "w", 10) is None
    for stage in [3, 4]:
        p = tmp_path / "w" / f"{stage:02d}_production_5ps/analysis_acceptance.json"
        p.parent.mkdir(parents=True)
        p.write_text("{}")
    assert len(sources.checkpoint_receipts(tmp_path, "w", 10)) == 2
    with pytest.raises(ValueError):
        sources.checkpoint_receipts(tmp_path, "w", 5)


def test_actual61_lineages_counts_and_discard() -> None:
    """Oppositestartsstaydistinct;old4haveonlyONEretained5pschunkat10ps."""
    out = sources.BASE / "analysis/pooled_support_20261127"
    data = sources.pooled_data(out / "baseline_cache", 10, False)
    assert data is not None and len(data["states"]) == 61
    assert data["counts"].tolist() == [800] * 57 + [100] * 4
    assert data["q"].shape == (46000, 3) and not data["bootstrap_eligible"]
    assert len({s["state_id"] for s in data["states"]}) == 61
    assert all(s["discard_once_ps"] == 5 for s in data["states"])
    with pytest.raises(ValueError):
        sources.pooled_data(out / "baseline_cache", 12, False)


def test_independent_resampling_and_no_one_chunk_fallback() -> None:
    """Resamplewithinstateoffsets;duplicatebiasesneverjustifyconcatenatingstarts."""
    counts = np.array([800, 200, 200])
    ids = estimator.resampled_indices(counts, 100, 17)
    assert len(ids) == 1200
    assert (ids[:800] < 800).all()
    assert ((ids[800:1000] >= 800) & (ids[800:1000] < 1000)).all()
    assert (ids[1000:] >= 1000).all()
    with pytest.raises(ValueError):
        estimator.resampled_indices(np.array([800, 100]), 100, 3)
    with pytest.raises(ValueError):
        estimator.paired_bootstrap(np.empty((2, 900)), np.array([800, 100]), {}, 100)
    with pytest.raises(ValueError):
        estimator.paired_bootstrap(np.empty((2, 1200)), counts, {}, 100, workers=8)


def test_histograms_use_actual_global_weights() -> None:
    """LocalhistogramsretainKishESSanddominant-framefractionswithoutKDE."""
    lw = np.log(np.array([0.8, 0.1, 0.1]))
    defs = dict(
        x=dict(
            values=np.array([0.2, 0.3, 0.8]),
            mask=np.ones(3, bool),
            edges=np.array([0, 0.5, 1]),
            reference=0,
        )
    )
    out = estimator.evaluate(lw, defs, np.arange(3))["x"]
    assert out["counts"] == [2, 1]
    assert out["maximum_frame_weight_fraction"][0] == pytest.approx(8 / 9)
    assert out["pmf"][0] == pytest.approx(0)
    assert out["pmf"][1] > 0


def test_paired_bootstrap_same_draws_all_coordinates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Allprojectionintervalsreceiveidenticalresampledglobalfit,notindependentdraws."""
    counts = np.array([200, 200])
    reduced = np.zeros((2, 400))
    p = dict(
        values=np.r_[np.zeros(100), np.ones(100), np.zeros(100), np.ones(100)],
        mask=np.ones(400, bool),
        edges=np.array([-0.5, 0.5, 1.5]),
        reference=0,
    )
    monkeypatch.setattr(
        estimator,
        "fit_logweights",
        lambda u, n: (np.full(u.shape[1], -np.log(u.shape[1])), dict(connected=True)),
    )
    result = estimator.paired_bootstrap(
        reduced, counts, dict(first=p, second=p), 100, workers=1
    )
    assert result["first"] == result["second"]
    assert len(result["first"]["replicates_pmf"]) == 200


def test_point_only_gate_skips_all_bootstraps(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Common10pswritesexplicitunresolvedCIandcannotstartshorter/IIDfallback."""
    report = dict(
        states=61, frames=3, counts=[3], f_k=[0], projections={}, connected=True
    )
    monkeypatch.setattr(
        analysis, "point", lambda data: (report, {}, np.zeros((1, 3)), np.zeros(3))
    )
    monkeypatch.setattr(analysis, "plot_report", lambda *args: None)
    monkeypatch.setattr(
        analysis,
        "paired_bootstrap",
        lambda *args, **kwargs: pytest.fail("bootstrapmustnotrunat10ps"),
    )
    data = dict(
        q=np.zeros((3, 3)),
        counts=np.array([3]),
        states=[],
        supplemental_horizon_ps=10,
        bootstrap_eligible=False,
    )
    out = tmp_path / "point"
    analysis.analyze(data, out)
    acceptance = json.loads((out / "acceptance.json").read_text())
    assert (
        acceptance["bootstrap_refits"] == 0 and "UNRESOLVED" in acceptance["CI_status"]
    )


def test_real_cached_point_reprojection_and_plot_no_refit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reproduceaccepted61pointplotsfromcachedweightswithoutanotherMBARfit."""
    root = sources.BASE / "analysis/pooled_support_20261127"
    data = json.loads((root / "first61_inputs.json").read_text())
    arrays = np.load(root / "first61_inputs.npz")
    data.update(q=arrays["q"], counts=arrays["counts"])
    accepted = json.loads((root / "first61_point10ps/point.json").read_text())
    weights = np.load(root / "first61_point10ps/global_weights.npz")
    monkeypatch.setattr(
        estimator,
        "fit_logweights",
        lambda u, n: (
            weights["logweights"],
            dict(
                connected=accepted["connected"],
                overlap=np.array(accepted["overlap"]),
                f_k=np.array(accepted["f_k"]),
            ),
        ),
    )
    result = analysis.analyze(data, tmp_path / "reprojection")
    assert result["counts"] == accepted["counts"]
    np.testing.assert_allclose(
        np.array(result["projections"]["hard"]["pmf"], float),
        np.array(accepted["projections"]["hard"]["pmf"], float),
        equal_nan=True,
    )
    assert (tmp_path / "reprojection/pooled_point.png").is_file()
    assert (tmp_path / "reprojection/global_weights.npz").is_file()


def test_baseline_cache_build_and_corruption_detection(tmp_path: Path) -> None:
    """Hashvalidatedcachebuildpreservesall57lineagesanddetectssubsequentchanges."""
    data = sources.frozen_main(tmp_path / "cache")
    assert len(data["states"]) == 57 and data["counts"].tolist() == [800] * 57
    assert data["q"].shape == (45600, 3)
    (tmp_path / "cache/main57_frozen45ps.npz").write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="identitychanged"):
        sources.frozen_main(tmp_path / "cache")


def test_eligible_branch_produces_both_paired_CIs_and_halves(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """At15ps, bothblockwidthsandindependent-statehalvesareexplicitlyrequested."""
    report = dict(
        states=94,
        frames=400,
        counts=[200, 200],
        f_k=[0, 0],
        projections={},
        connected=True,
    )
    calls = []

    def fake_point(data: dict) -> tuple:
        calls.append(("point", data["counts"].tolist()))
        return dict(report), {}, np.zeros((2, len(data["q"]))), np.zeros(len(data["q"]))

    def fake_boot(
        u: np.ndarray, n: np.ndarray, p: dict, chunk: int, workers: int
    ) -> dict:
        calls.append(("bootstrap", chunk, workers))
        return {}

    monkeypatch.setattr(analysis, "point", fake_point)
    monkeypatch.setattr(analysis, "plot_report", lambda *args: None)
    monkeypatch.setattr(analysis, "paired_bootstrap", fake_boot)
    data = dict(
        q=np.zeros((400, 3)),
        counts=np.array([200, 200]),
        states=[],
        supplemental_horizon_ps=15,
        bootstrap_eligible=True,
    )
    analysis.analyze(data, tmp_path / "eligible")
    assert ("bootstrap", 100, 4) in calls and ("bootstrap", 50, 4) in calls
    assert calls.count(("point", [100, 100])) == 2
    assert (tmp_path / "eligible/bootstrap_5ps.json").exists()
    assert (tmp_path / "eligible/bootstrap_2.5ps.json").exists()
    assert (
        json.loads((tmp_path / "eligible/acceptance.json").read_text())[
            "bootstrap_refits"
        ]
        == 400
    )


def test_watch_10_then_15_to30_nonprimary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Automaticisolatedmonitorrunsallfivecohortcheckpointsat4CPUwithoutprimarywrites."""
    from scripts import watch_relay_pooled_support as watcher

    old = tmp_path / "old/prepared.json"
    new = tmp_path / "new/prepared.json"
    extension = old.parent / "extension_20261127/manifest.json"
    for p in [old, new, extension]:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("{}")
    monkeypatch.setattr(watcher, "OLD", old)
    monkeypatch.setattr(watcher, "NEW", new)
    events = []
    monkeypatch.setattr(watcher, "journal", lambda *args: None)
    monkeypatch.setattr(
        watcher,
        "pooled_data",
        lambda c, h, include: dict(
            supplemental_horizon_ps=h, bootstrap_eligible=h >= 15
        ),
    )

    def fake_analysis(data: dict, out: Path, workers: int) -> None:
        events.append(
            (data["supplemental_horizon_ps"], data["bootstrap_eligible"], workers)
        )
        out.mkdir()
        (out / "acceptance.json").write_text(
            json.dumps(dict(states=94, horizon_ps=data["supplemental_horizon_ps"]))
        )

    monkeypatch.setattr(watcher, "analyze", fake_analysis)
    watcher.watch(tmp_path, 1, 60)
    assert events == [
        (10, False, 4),
        (15, True, 4),
        (20, True, 4),
        (25, True, 4),
        (30, True, 4),
    ]
    assert (
        json.loads((tmp_path / "monitor_state.json").read_text())["primary_replacement"]
        is False
    )
    watcher.watch(tmp_path, 1, 60)
    assert len(events) == 5
    with pytest.raises(ValueError):
        watcher.watch(tmp_path, 100, 60)
