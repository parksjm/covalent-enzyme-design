"""Check physical distance extraction and independent window statistics."""

import json
from pathlib import Path

import pytest

from nyl12_covalent.reaction_distance_stats import (
    DISTANCES,
    campaign_metadata,
    extract_row,
    load_sources,
    movie_rows,
    ownership_summary,
    sha256,
    summarize,
    write_csv,
)


def record() -> dict[str, str]:
    """Return one explicit contact record."""
    return {
        "r_og1_hg1": "1.8",
        "r_nleave_hg1": "3",
        "r_c_og1": "2.4",
        "r_c_n": "1.4",
        "time_ps": "6.05",
        "hg1_owner": "13905",
        "nalpha_h1_owner": "13905",
        "nalpha_h2_owner": "13916",
        "hydrogen_owner_distances_A": "{13917: {13905: 1.1}}",
    }


def test_extract_exact_contact_and_reject_invalid() -> None:
    """Use the direct distance irrespective of the current proton owner."""
    row = record()
    assert extract_row(row)["r_nalpha_hg1"] == 1.1
    row["r_c_n"] = "nan"
    with pytest.raises(ValueError, match="Invalid distance"):
        extract_row(row)


def test_discard_statistics_ownership_and_join(tmp_path: Path) -> None:
    """Discard only the initial five ps and preserve exact movie indices."""
    rows = []
    for i in range(900):
        rows.append(
            {
                **extract_row(record()),
                "window": "A",
                "block": "03",
                "frame": i + 1,
                "production_time_ps": (i + 1) * 0.05,
                "r_c_n": float(i),
            }
        )
    stats = summarize(rows, ["A"])
    cn = next(r for r in stats if r["distance"] == "r_c_n")
    assert cn["n"] == 800
    assert cn["early_5ps_mean_A"] == 149.5
    assert cn["late_5ps_mean_A"] == 849.5
    assert cn["late_minus_early_A"] == 700
    assert ownership_summary(rows, ["A"])[0]["owner_fractions"] == {"13905": 1}
    plan = [
        {
            "output_frame0based": 0,
            "window": "A",
            "block": "03",
            "source_frame1based": 1,
            "window_time_ps": 6.05,
        }
    ]
    assert movie_rows(rows, plan)[0]["movie_frame0based"] == 0
    plan[0]["window_time_ps"] = 99
    with pytest.raises(ValueError, match="timestamp"):
        movie_rows(rows, plan)
    plan[0]["output_frame0based"] = 2
    with pytest.raises(ValueError, match="Nonsequential"):
        movie_rows(rows, plan)
    with pytest.raises(ValueError, match="ten retained"):
        summarize(rows[:200], ["A"])
    path = tmp_path / "stats.csv"
    write_csv(path, stats)
    assert "late_minus_early_A" in path.read_text()
    assert len(sha256(path)) == 64
    assert len(stats) == len(DISTANCES)


@pytest.mark.parametrize("blocks", [9, 12])
def test_load_sources_exact_declared_blocks(tmp_path: Path, blocks: int) -> None:
    """Load all 900 samples with explicit production times and source hashes."""
    import csv

    manifest = {
        "windows": ["A"],
        "sampled_horizon_per_window_ps": blocks * 5,
        "production_blocks_per_window": blocks,
        "frame_plan": [],
    }
    for block in range(3, 3 + blocks):
        path = tmp_path / "A" / f"{block:02d}_production_5ps/raw_features.csv"
        path.parent.mkdir(parents=True)
        with path.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(record()))
            writer.writeheader()
            for frame in range(1, 101):
                row = record()
                row["time_ps"] = str(6 + (block - 3) * 5 + frame * 0.05)
                writer.writerow(row)
        manifest["frame_plan"].append({"window": "A", "block": path.parent.name})
    rows, sources = load_sources(tmp_path, manifest)
    assert len(rows) == blocks * 100 and len(sources) == blocks
    assert rows[-1]["production_time_ps"] == blocks * 5
    assert rows[0]["production_time_ps"] == 0.05
    (tmp_path / "A/03_production_5ps/raw_features.csv").write_text("a\n")
    with pytest.raises(ValueError, match="Empty native"):
        load_sources(tmp_path, manifest)
    json.dumps(sources)


def test_coordinate_sign_consistency() -> None:
    """Require the physically defined signed coordinate to match direct contacts."""
    row = record()
    row["q_act_A"] = "0.7"
    assert extract_row(row)["r_nalpha_hg1"] == 1.1
    row["q_act_A"] = "-0.7"
    with pytest.raises(ValueError, match="activation-coordinate"):
        extract_row(row)


@pytest.mark.parametrize("coordinate,value", [("q_del_A", -1.9), ("q_acyl_A", -1.0)])
def test_other_coordinate_consistency(coordinate: str, value: float) -> None:
    """Check delivery and acyl coordinates against separately stored distances."""
    row = record()
    row[coordinate] = str(value)
    extract_row(row)
    row[coordinate] = str(value + 1)
    with pytest.raises(ValueError, match="Direct distance"):
        extract_row(row)


def test_sixty_ps_non90_movie_metadata() -> None:
    """Derive60ps horizons, unequal movie counts and exact phase boundary."""
    rows = [
        {**extract_row(record()), "window": w, "production_time_ps": i * 0.05}
        for w in ["A", "D"]
        for i in range(1, 1201)
    ]
    movie = [
        {"window": "A" if i < 120 else "D", "movie_frame0based": i} for i in range(220)
    ]
    meta = campaign_metadata(
        movie, rows, ["A", "D"], {"A": "activation", "D": "delivery"}
    )
    assert meta["retained_count"] == 2200
    assert meta["movie_count"] == 220
    assert meta["phase_boundary_movie_frames"] == [119.5]
    assert meta["phase_boundary_window_positions"] == [1.5]
    assert [r["movie_count"] for r in meta["windows"]] == [120, 100]
    assert [r["horizon_ps"] for r in meta["windows"]] == [60, 60]
    assert summarize(rows, ["A"])[0]["n"] == 1100
    with pytest.raises(ValueError, match="noncontiguous"):
        campaign_metadata(movie[1:], rows, ["A", "missing"], {"A": "activation"})
