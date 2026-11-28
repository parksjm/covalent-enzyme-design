"""Validate editable scientific slide construction and source failures."""

from __future__ import annotations

import hashlib
import importlib
import json
from pathlib import Path

import pytest
from PIL import Image

from scripts.make_presentation_contact_sheet import make_contact_sheet

presentation_runtime = pytest.importorskip("pptx")
presentation_builder = importlib.import_module("scripts.build_nyl12_progress_deck")
Presentation = presentation_runtime.Presentation
build_deck = presentation_builder.build_deck
main = presentation_builder.main
validate_source = presentation_builder.validate_source


def _source(tmp_path: Path) -> Path:
    figure = tmp_path / "figure.png"
    Image.new("RGB", (200, 80), "white").save(figure)
    tall = tmp_path / "tall.png"
    Image.new("RGB", (80, 200), "white").save(tall)
    slides = [
        {"kind": "title", "title": "Research", "bullets": ["Scope", "Limits"]},
        {"kind": "text", "title": "Evidence", "bullets": ["Measured", "Modeled"]},
        {
            "kind": "image",
            "title": "Profile",
            "image": "figure.png",
            "takeaway": "Potential energy",
            "caption": "Not a PMF",
        },
        {
            "kind": "split",
            "title": "Audit",
            "image": "tall.png",
            "bullets": ["Finite sampling"],
        },
        {
            "kind": "table",
            "title": "Rates",
            "columns": ["T", "Rate"],
            "rows": [[300, "1e-4"], [348.15, "1e-2"]],
            "widths": [6, 6.2],
        },
        {
            "kind": "cards",
            "title": "Decisions",
            "cards": [{"title": "Pass", "text": "Geometry"}] * 3,
        },
        {
            "kind": "snapshots",
            "title": "Snapshot placeholders",
            "cards": [{"title": "E-S", "text": "Manual snapshot required"}] * 4,
        },
    ]
    for slide in slides:
        slide["notes"] = "Scientific caveat retained."
        slide["sources"] = ["/absolute/primary/source.json"]
    path = tmp_path / "slides.json"
    path.write_text(json.dumps({"title": "Test", "slides": slides}))
    return path


def test_build_retains_notes_figures_tables_and_receipt(tmp_path: Path) -> None:
    source = _source(tmp_path)
    output = tmp_path / "report.pptx"
    receipt = build_deck(source, output)
    result = Presentation(output)
    assert len(result.slides) == receipt["slide_count"] == 7
    assert receipt["pptx_sha256"] == hashlib.sha256(output.read_bytes()).hexdigest()
    assert len(receipt["assets_sha256"]) == 2
    for slide in result.slides:
        assert "Scientific caveat retained" in slide.notes_slide.notes_text_frame.text
        assert (
            "/absolute/primary/source.json" in slide.notes_slide.notes_text_frame.text
        )
        for shape in slide.shapes:
            assert shape.left >= 0 and shape.top >= 0
            assert shape.left + shape.width <= result.slide_width
            assert shape.top + shape.height <= result.slide_height
    tables = [shape.table for shape in result.slides[4].shapes if shape.has_table]
    assert tables[0].cell(2, 0).text == "348.15"
    assert output.with_suffix(".build.json").is_file()


@pytest.mark.parametrize(
    "spec",
    [
        {"slides": []},
        {"slides": [{"kind": "unknown", "title": "No"}]},
        {"slides": [{"kind": "text"}]},
        {
            "slides": [
                {"kind": "table", "title": "Ragged", "columns": ["A"], "rows": [[1, 2]]}
            ]
        },
        {"slides": [{"kind": "cards", "title": "Empty", "cards": []}]},
    ],
)
def test_invalid_scientific_slide_source_fails(spec: dict, tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        validate_source(spec, tmp_path)


def test_missing_figure_fails_before_output(tmp_path: Path) -> None:
    spec = {"slides": [{"kind": "image", "title": "Missing", "image": "no.png"}]}
    with pytest.raises(FileNotFoundError):
        validate_source(spec, tmp_path)


def test_command_line_build(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = _source(tmp_path)
    output = tmp_path / "cli.pptx"
    monkeypatch.setattr(
        "sys.argv", ["builder", "--source", str(source), "--output", str(output)]
    )
    main()
    assert len(Presentation(output).slides) == 7


def test_cover_date_and_scientific_callout_layouts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from datetime import datetime
    from zoneinfo import ZoneInfo

    monkeypatch.setattr(
        presentation_builder,
        "build_date",
        lambda: datetime(2026, 11, 25, tzinfo=ZoneInfo("America/New_York")),
    )
    figure = tmp_path / "figure.png"
    Image.new("RGB", (160, 90), "white").save(figure)
    source = tmp_path / "layouts.json"
    source.write_text(
        json.dumps(
            {
                "title": "Layouts",
                "slides": [
                    {
                        "kind": "title",
                        "title": "Mechanism-guided redesign",
                        "subtitle": "QM/MM pathways and sequence design",
                    },
                    {
                        "kind": "process",
                        "title": "Modeling workflow",
                        "lanes": [{"title": "Path", "steps": ["MD", "QM/MM", "PMF"]}],
                    },
                    {
                        "kind": "annotated_image",
                        "title": "Direct-seeded data",
                        "image": "figure.png",
                        "callouts": [
                            "Mixed proton-transfer channels",
                            "Sampling not converged",
                        ],
                    },
                    {
                        "kind": "placeholder",
                        "title": "Free-energy comparison",
                        "status": "Analysis pending",
                        "message": "Common-reference comparison is not yet available.",
                    },
                ],
            }
        )
    )
    output = tmp_path / "layouts.pptx"
    receipt = build_deck(source, output)
    deck = Presentation(output)
    cover_text = " ".join(
        shape.text for shape in deck.slides[0].shapes if shape.has_text_frame
    )
    assert "Mechanism-guided redesign" in cover_text
    assert "25 November 2026" in cover_text
    assert receipt["build_date"] == "2026-11-25"
    callout_text = " ".join(
        shape.text for shape in deck.slides[2].shapes if shape.has_text_frame
    )
    assert "Mixed proton-transfer channels" in callout_text
    placeholder_text = " ".join(
        shape.text for shape in deck.slides[3].shapes if shape.has_text_frame
    )
    assert "Common-reference comparison is not yet available" in placeholder_text


def test_contact_sheet_orders_every_numbered_slide(tmp_path: Path) -> None:
    for number in (2, 1, 3):
        Image.new(
            "RGB", (160, 90), f"#{number}{number}{number}{number}{number}{number}"
        ).save(tmp_path / f"slide-{number:02}.png")
    output = tmp_path / "contact_sheet.png"
    count = make_contact_sheet(tmp_path, output)
    assert count == 3
    with Image.open(output) as sheet:
        assert sheet.size == (4 * (480 + 12) + 12, 1 * (270 + 12) + 12)
