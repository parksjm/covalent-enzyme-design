"""Check mapping and periodic imaging against retained four-site endpoints."""

from pathlib import Path
import json
import numpy as np
import parmed
import pytest
from nyl12_covalent.design_md_geometry import (
    frame_rows,
    image_near,
    index_mapping,
    native_reference,
)


def test_triclinic_periodic_image() -> None:
    box = np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    points = np.array([[9.8, 0.0, 0.0], [-9.8, 0.0, 0.0]])
    result = image_near(points, np.zeros(3), box)
    assert result[:, 0] == pytest.approx([-0.2, 0.2])
    with pytest.raises(ValueError):
        image_near(points, np.zeros(3), np.zeros(6))


@pytest.mark.parametrize(
    "name,expected_tiers,expected_native",
    [
        ("WT", ["A"] * 4, [False] * 4),
        ("broader_mask07", ["A", "B", "A", "B"], [True, False, True, False]),
    ],
)
def test_original_mappings_match_minimized_sites(
    name: str, expected_tiers: list[str], expected_native: list[bool]
) -> None:
    root = Path(__file__).resolve().parents[1]
    p = root / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup"
    entry = json.loads((p / "systems" / name / "datasets.json").read_text())[
        "datasets"
    ][name]
    top = (p / "systems" / name / entry["topology"]["path"]).resolve()
    rst = (p / "systems" / name / entry["restart"]["path"]).resolve()
    t = parmed.load_file(str(top), xyz=str(rst))
    source = Path(
        json.loads((p / "selection.json").read_text())["source_selection"][name]
    )
    m = index_mapping(t, source)
    rows = frame_rows(
        np.asarray(t.coordinates),
        np.asarray(t.box),
        m,
        native_reference(root / "cif/n12c.cif"),
    )
    assert [r["tier"] for r in rows] == expected_tiers
    assert [r["native_pass"] for r in rows] == expected_native
    assert all(r["herschlag_phi"] < 0 for r in rows)
