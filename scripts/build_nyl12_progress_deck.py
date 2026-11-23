"""Build an editable scientific progress deck from reviewed JSON slide sources.

Run in an isolated report environment; python-pptx is not a production dependency.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
from typing import Any

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.util import Inches, Pt

LOGGER = logging.getLogger(__name__)
NAVY = "17324D"
TEAL = "008C95"
INK = "223548"
MUTED = "526577"
PALE = "EFF5F8"
WHITE = "FFFFFF"
FONT = "DejaVu Sans"


def validate_source(source: dict[str, Any], base: Path) -> None:
    """Reject incomplete slides and unreadable figure inputs before building.

    Args:
        source: Reviewed slide specification.
        base: Directory containing the specification and relative assets.
    """
    slides = source.get("slides", [])
    if not slides:
        raise ValueError("The deck must contain slides")
    supported = {"title", "text", "image", "split", "table", "cards", "snapshots"}
    for index, spec in enumerate(slides, 1):
        if not spec.get("title") or spec.get("kind") not in supported:
            raise ValueError(f"Invalid slide {index}")
        if spec["kind"] in {"image", "split"}:
            path = base / spec["image"]
            with Image.open(path) as figure:
                figure.verify()
        if spec["kind"] == "table":
            rows = [spec["columns"], *spec["rows"]]
            if not rows[0] or any(len(row) != len(rows[0]) for row in rows):
                raise ValueError(f"Ragged table on slide {index}")
        if spec["kind"] in {"cards", "snapshots"}:
            if not 1 <= len(spec.get("cards", [])) <= 4:
                raise ValueError(f"Invalid card count on slide {index}")


def _box(
    slide: Any,
    text: str,
    x: float,
    y: float,
    width: float,
    height: float,
    size: float = 22,
    color: str = INK,
    bold: bool = False,
) -> Any:
    shape = slide.shapes.add_textbox(
        Inches(x), Inches(y), Inches(width), Inches(height)
    )
    frame = shape.text_frame
    frame.clear()
    frame.word_wrap = True
    for index, line in enumerate(text.split("\n")):
        paragraph = frame.paragraphs[0] if index == 0 else frame.add_paragraph()
        paragraph.text = line
        paragraph.font.name = FONT
        paragraph.font.size = Pt(size)
        paragraph.font.bold = bold
        paragraph.font.color.rgb = RGBColor.from_string(color)
        paragraph.space_after = Pt(12)
    return shape


def _rect(
    slide: Any, x: float, y: float, width: float, height: float, color: str
) -> None:
    shape = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(x),
        Inches(y),
        Inches(width),
        Inches(height),
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = RGBColor.from_string(color)
    shape.line.fill.background()


def _picture(slide: Any, path: Path, area: tuple[float, ...]) -> None:
    x, y, width, height = area
    with Image.open(path) as figure:
        ratio = figure.width / figure.height
    actual_w, actual_h = width, width / ratio
    if actual_h > height:
        actual_h, actual_w = height, height * ratio
    slide.shapes.add_picture(
        str(path),
        Inches(x + (width - actual_w) / 2),
        Inches(y + (height - actual_h) / 2),
        width=Inches(actual_w),
        height=Inches(actual_h),
    )


def _table(slide: Any, spec: dict[str, Any]) -> None:
    rows = [spec["columns"], *spec["rows"]]
    shape = slide.shapes.add_table(
        len(rows),
        len(rows[0]),
        Inches(0.55),
        Inches(1.65),
        Inches(12.2),
        Inches(spec.get("table_height", 3.8)),
    )
    table = shape.table
    if "widths" in spec:
        for column, width in zip(table.columns, spec["widths"], strict=True):
            column.width = Inches(width)
    for i, row in enumerate(rows):
        for j, value in enumerate(row):
            cell = table.cell(i, j)
            cell.text = str(value)
            cell.margin_left = Inches(0.12)
            cell.margin_top = Inches(0.1)
            cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor.from_string(
                NAVY if i == 0 else WHITE if i % 2 == 0 else PALE
            )
            for paragraph in cell.text_frame.paragraphs:
                paragraph.font.name = FONT
                paragraph.font.size = Pt(spec.get("font_size", 19))
                paragraph.font.bold = i == 0
                paragraph.font.color.rgb = RGBColor.from_string(
                    WHITE if i == 0 else INK
                )


def _cards(slide: Any, spec: dict[str, Any]) -> None:
    cards = spec["cards"]
    columns = 2 if len(cards) in {2, 4} else len(cards)
    width = (12.2 - 0.3 * (columns - 1)) / columns
    height = 2.15 if len(cards) == 4 else 4.2
    for index, card in enumerate(cards):
        x = 0.55 + (index % columns) * (width + 0.3)
        y = 1.72 + (index // columns) * (height + 0.2)
        _rect(slide, x, y, width, height, PALE)
        _box(
            slide,
            card["title"],
            x + 0.2,
            y + 0.15,
            width - 0.4,
            0.5,
            size=20,
            color=TEAL,
            bold=True,
        )
        _box(
            slide,
            card["text"],
            x + 0.2,
            y + 0.72,
            width - 0.4,
            height - 0.82,
            size=spec.get("font_size", 19),
        )


def build_deck(source_path: Path, destination: Path) -> dict[str, Any]:
    """Create slides, full source notes, and a hashed build receipt.

    Args:
        source_path: JSON specification with reviewed claims and figure paths.
        destination: Output PowerPoint file.

    Returns:
        Build receipt containing slide titles, hashes and source provenance.
    """
    source = json.loads(source_path.read_text())
    validate_source(source, source_path.parent)
    presentation = Presentation()
    presentation.slide_width = Inches(13.333)
    presentation.slide_height = Inches(7.5)
    assets = {}
    for index, spec in enumerate(source["slides"], 1):
        slide = presentation.slides.add_slide(presentation.slide_layouts[6])
        _rect(slide, 0, 0, 13.333, 0.12, TEAL)
        _box(
            slide,
            spec["title"],
            0.55,
            0.35,
            12.2,
            0.8,
            size=spec.get("title_size", 29),
            color=NAVY,
            bold=True,
        )
        if spec.get("takeaway"):
            _box(
                slide,
                spec["takeaway"],
                0.55,
                1.15,
                12.2,
                0.5,
                size=17,
                color=TEAL,
                bold=True,
            )
        kind = spec["kind"]
        if kind in {"title", "text"}:
            text = "\n".join(spec.get("bullets", []))
            _box(slide, text, 0.65, 1.95, 11.9, 4.55, size=spec.get("font_size", 24))
        elif kind in {"image", "split"}:
            path = source_path.parent / spec["image"]
            area = (0.55, 1.75, 12.2, 4.75)
            if kind == "split":
                _box(
                    slide,
                    "\n".join(spec["bullets"]),
                    0.55,
                    1.9,
                    4.55,
                    4.55,
                    size=spec.get("font_size", 20),
                )
                area = (5.35, 1.75, 7.4, 4.75)
            _picture(slide, path, area)
            assets[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        elif kind == "table":
            _table(slide, spec)
            _box(
                slide,
                spec.get("caption", ""),
                0.65,
                6.0,
                12.0,
                0.65,
                size=16,
                color=MUTED,
            )
        else:
            _cards(slide, spec)
        footer = spec.get("footer", "Nyl12–PA6 | Research progress | 23 Nov 2026")
        _box(slide, footer, 0.55, 7.05, 11.6, 0.25, size=10, color=MUTED)
        _box(slide, str(index), 12.35, 7.05, 0.4, 0.25, size=10, color=MUTED)
        if kind in {"image", "split"} and spec.get("caption"):
            _box(slide, spec["caption"], 0.65, 6.55, 12.0, 0.4, size=13, color=MUTED)
        notes = [spec.get("notes", ""), "Sources:", *spec.get("sources", [])]
        slide.notes_slide.notes_text_frame.text = "\n".join(notes)
    destination.parent.mkdir(parents=True, exist_ok=True)
    presentation.save(destination)
    receipt = {
        "title": source["title"],
        "slide_count": len(source["slides"]),
        "slide_titles": [slide["title"] for slide in source["slides"]],
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "builder_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "pptx_sha256": hashlib.sha256(destination.read_bytes()).hexdigest(),
        "assets_sha256": assets,
        "scientific_status": source.get("scientific_status"),
    }
    destination.with_suffix(".build.json").write_text(
        json.dumps(receipt, indent=2) + "\n"
    )
    LOGGER.info("Built %s (%d slides)", destination, receipt["slide_count"])
    return receipt


def main() -> None:
    """Parse report inputs and build the PowerPoint."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    build_deck(args.source, args.output)


if __name__ == "__main__":
    main()
