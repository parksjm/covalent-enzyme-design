"""Synthetic statistical and source-integrity checks for the relay adapter."""

import json
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent.relay_umbrella_analysis import (
    KB,
    analyze,
    path_histogram,
    select_uncorrelated,
)
from nyl12_covalent.relay_umbrella_analysis_io import (
    analyze_request,
    digest,
    load_block,
)
from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from nyl12_covalent.relay_umbrella_projection import project_path


def fixtures(separation: float = 0.6) -> tuple[list, list, dict]:
    """Create exact equilibrium samples under two different bias chart families."""
    rng = np.random.default_rng(728)
    windows, blocks = [], []
    physical_k, restraint_k, temperature = 1.0, 1.0, 300.0
    for i, (pair, center) in enumerate(
        [
            (("q_acyl", "q_act"), -separation),
            (("q_acyl", "q_del"), separation),
        ]
    ):
        bias = RelayBiasWindow(pair, (center, 0.0), (restraint_k, restraint_k))
        windows.append(
            dict(
                label=f"w{i}",
                bias=dict(
                    feature_names=pair,
                    centers_A=bias.centers_A,
                    force_constants=bias.force_constants,
                    outer_extent_A=10.0,
                ),
            )
        )
        means = np.array([restraint_k * center / (physical_k + restraint_k), 0.0, 0.0])
        force = np.array([physical_k + restraint_k, physical_k, physical_k])
        force[1 if i == 0 else 2] += restraint_k
        samples = rng.normal(
            means, np.sqrt(KB * temperature / (2 * force)), size=(12000, 3)
        )
        raw = dict(zip(("q_acyl", "q_act", "q_del"), samples.T, strict=True))
        blocks.append(
            dict(
                window=f"w{i}",
                values=samples,
                native_bias_kcal_mol=cross_bias_energies(raw, [bias])[0],
                sampled_duration_ps=5.0,
            )
        )
    kwargs = dict(
        path=np.array([[-2.0, 0, 0], [2.0, 0, 0]]),
        temperature=temperature,
        edges=np.linspace(-0.1, 1.1, 17),
        es_bounds=(0.4, 0.6),
        acyl_bounds=(0.8, 1.0),
        neighbor_pairs=[[0, 1]],
    )
    return windows, blocks, kwargs


def test_projection_extends_endpoints_and_reports_ambiguity() -> None:
    result = project_path(
        np.array([[-1, 0, 0], [3, 1, 0]]), np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    )
    np.testing.assert_allclose(result["s"], [-0.5, 1.5])
    np.testing.assert_allclose(result["transverse_A"], [0, 1])
    with pytest.raises(ValueError, match="differ"):
        project_path(np.zeros((1, 3)), np.zeros((2, 3)))
    with pytest.raises(ValueError, match="shape"):
        project_path(np.zeros((1, 2)), np.zeros((2, 3)))
    with pytest.raises(ValueError, match="Finite"):
        project_path(np.full((1, 3), np.nan), np.ones((2, 3)))


def test_shared_multichart_recovers_known_unbiased_harmonic() -> None:
    windows, blocks, kwargs = fixtures()
    result = analyze(windows, blocks, **kwargs)
    assert result["status"] == "preliminary_connected_not_converged"
    assert result["bias_free_energies_are_not_PMF"]
    pmf = result["pmf"]
    centers = np.array(pmf["centers"])
    measured = np.array([np.nan if v is None else v for v in pmf["pmf_kcal_mol"]])
    x = -2 + 4 * centers
    use = np.abs(x) < 0.9
    expected = x**2
    difference = measured[use] - expected[use]
    assert np.max(np.abs(difference - np.mean(difference))) < 0.13
    assert result["converged"] is False


def test_missing_and_disconnected_cannot_produce_global_pmf() -> None:
    windows, blocks, kwargs = fixtures()
    assert analyze(windows, [], **kwargs)["status"] == "readiness_no_accepted_data"
    assert analyze(windows, blocks[:1], **kwargs)["pmf"] is None
    windows, blocks, kwargs = fixtures(20.0)
    # Same physical U0, but virtually no bridge-state overlap.
    result = analyze(windows, blocks, **kwargs)
    assert result["status"] == "disconnected_no_global_PMF"
    assert result["pmf"] is None


def test_native_parity_and_invalid_inputs_fail() -> None:
    windows, blocks, kwargs = fixtures()
    blocks[0]["native_bias_kcal_mol"][0] += 1
    with pytest.raises(ValueError, match="parity"):
        analyze(windows, blocks, **kwargs)
    with pytest.raises(ValueError, match="Positive temperature"):
        analyze(windows, [], **(kwargs | {"temperature": 0}))
    with pytest.raises(ValueError, match="unique"):
        analyze([windows[0], windows[0]], [], **kwargs)
    with pytest.raises(ValueError, match="four finite"):
        select_uncorrelated(np.zeros((3, 4)))
    indices, summary = select_uncorrelated(np.ones((10, 4)))
    assert len(indices) == 10 and summary["discarded_frames"] == 0


def test_histogram_does_not_interpolate_gaps_or_invent_reference() -> None:
    data = np.array([0.1, 0.9])
    result = path_histogram(
        data, np.array([0.5, 0.5]), np.linspace(0, 1, 6), 300, (0, 0.2), (0.8, 1.0)
    )
    assert result["minimum_to_barrier_kcal_mol"] is None
    assert None in result["pmf_kcal_mol"]
    assert (
        path_histogram(
            data, np.ones(2) / 2, np.linspace(0, 1, 6), 300, (0.3, 0.5), (0.8, 1.0)
        )["status"]
        == "no_ES_reference_support"
    )
    with pytest.raises(ValueError, match="edges"):
        path_histogram(
            data, np.ones(2) / 2, np.array([0.2, 0.5, 0.8]), 300, (0, 0.2), (0.8, 1.0)
        )
    with pytest.raises(ValueError, match="interval"):
        path_histogram(
            data, np.ones(2) / 2, np.linspace(0, 1, 6), 300, (0.2, 0.1), (0.8, 1.0)
        )


def receipt(tmp_path: Path) -> tuple[Path, dict]:
    """Write a tiny accepted compact block with immutable source proofs."""
    csv_path = tmp_path / "raw_features.csv"
    csv_path.write_text(
        "time_ps,q_acyl,q_act,q_del,native_bias_kcal_mol,cn,owner,native_bias_evaluated,native_cv_1_A,native_cv_2_A\n"
        + "\n".join(f"{i},0,0,0,0,1,Nalpha,1,0,0" for i in range(5))
    )
    process = tmp_path / "process.json"
    process.write_text('{"exit_status": 0}')
    source = tmp_path / "input.json"
    source.write_text('{"immutable":true}')
    proof = dict(
        accepted=True,
        exit_status=0,
        frames_censored=0,
        coordinates_finite=True,
        energies_finite=True,
        temperature_K=300,
        physical_hamiltonian_sha256="a" * 64,
        csv_file=csv_path.name,
        csv_sha256=digest(csv_path),
        process_file=process.name,
        process_sha256=digest(process),
        source_sha256={source.name: digest(source)},
        window="w0",
        sampled_duration_ps=5,
        native_bias_source="amber_dumpave_independent_CVs",
        native_cv_quantization_A=0.0005,
        coordinate_serialization_tolerance_A=0.0,
        native_bias_evaluated_frames=5,
        native_bias_unevaluated_terminal_frames=0,
        frames=5,
    )
    path = tmp_path / "accepted.json"
    path.write_text(json.dumps(proof))
    return path, proof


def test_receipt_preserves_owners_and_rejects_mutation(tmp_path: Path) -> None:
    path, proof = receipt(tmp_path)
    result = load_block(
        path, temperature=300, hamiltonian="a" * 64, diagnostic_columns=["cn"]
    )
    assert result["original_rows"][0]["owner"] == "Nalpha"
    assert result["values"].shape == (5, 4)
    with pytest.raises(ValueError, match="temperature"):
        load_block(path, temperature=301, hamiltonian="a" * 64, diagnostic_columns=[])
    (tmp_path / "input.json").write_text("changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        load_block(path, temperature=300, hamiltonian="a" * 64, diagnostic_columns=[])
    proof["accepted"] = False
    path.write_text(json.dumps(proof))
    with pytest.raises(ValueError, match="contract"):
        load_block(path, temperature=300, hamiltonian="a" * 64, diagnostic_columns=[])


def test_request_empty_readiness_and_duplicate_rejection(tmp_path: Path) -> None:
    windows, _, kwargs = fixtures()
    cfg = dict(
        windows=windows,
        path_3cv_A=kwargs["path"].tolist(),
        temperature_K=300,
        physical_hamiltonian_sha256="a" * 64,
        s_bin_edges=kwargs["edges"].tolist(),
        es_search_bounds=[0.4, 0.6],
        acyl_search_bounds=[0.8, 1.0],
        neighbor_pairs=[[0, 1]],
        diagnostic_columns=[],
        accepted_block_receipts=[],
    )
    request = tmp_path / "request.json"
    request.write_text(json.dumps(cfg))
    assert analyze_request(request)["status"] == "readiness_no_accepted_data"
    path, _ = receipt(tmp_path)
    cfg["accepted_block_receipts"] = [path.name, path.name]
    request.write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match="counted twice"):
        analyze_request(request)


def test_readiness_validates_reference_and_neighbor_inputs() -> None:
    windows, _, kwargs = fixtures()
    with pytest.raises(ValueError, match="Consecutive"):
        analyze(windows, [], **(kwargs | dict(path=np.zeros((2, 3)))))
    with pytest.raises(ValueError, match="histogram edges"):
        analyze(windows, [], **(kwargs | dict(edges=np.array([0, 0, 1]))))
    with pytest.raises(ValueError, match="search intervals"):
        analyze(windows, [], **(kwargs | dict(es_bounds=(1, 0))))
    with pytest.raises(ValueError, match="Neighbor pairs"):
        analyze(windows, [], **(kwargs | dict(neighbor_pairs=[[0, 2]])))


def test_native_dump_rounding_and_terminal_mask_retain_every_sample() -> None:
    windows, blocks, kwargs = fixtures()
    for index, block in enumerate(blocks):
        names = windows[index]["bias"]["feature_names"]
        columns = [("q_acyl", "q_act", "q_del").index(name) for name in names]
        actual = block["values"][:, columns]
        rounded = np.round(actual, 3)
        features = {name: rounded[:, j] for j, name in enumerate(names)}
        block["native_bias_kcal_mol"] = cross_bias_energies(
            features, [RelayBiasWindow(**windows[index]["bias"])]
        )[0]
        block["native_cv_values_A"] = rounded
        block["native_cv_quantization_A"] = 0.0005
        block["native_bias_evaluated"] = np.ones(len(actual), bool)
        block["native_bias_evaluated"][-1] = False
        block["native_bias_kcal_mol"][-1] = np.nan
        block["native_cv_values_A"][-1] = np.nan
    result = analyze(windows, blocks, **kwargs)
    assert result["status"] == "preliminary_connected_not_converged"
    assert all(item["raw_frames"] == 12000 for item in result["blocks"])
    assert all(
        item["native_bias_unevaluated_terminal_frames"] == 1
        for item in result["blocks"]
    )
    blocks[0]["native_bias_evaluated"][0] = False
    with pytest.raises(ValueError, match="terminal"):
        analyze(windows, blocks, **kwargs)
