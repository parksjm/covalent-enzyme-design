"""Editable PowerPoint shape helpers for scientific slide content."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from PIL import Image
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.util import Inches, Pt

from scripts.presentation_layouts import FONT, INK, MUTED, NAVY, PALE, TEAL, WHITE


def box(
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
    """Add a wrapped text box with consistent typography."""
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


def rect(
    slide: Any, x: float, y: float, width: float, height: float, color: str
) -> None:
    """Add a filled rectangle without an outline."""
    shape = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(width), Inches(height)
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = RGBColor.from_string(color)
    shape.line.fill.background()


def picture(slide: Any, path: Path, area: tuple[float, ...]) -> None:
    """Place an image inside an area while preserving its aspect ratio."""
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


def table(slide: Any, spec: dict[str, Any]) -> None:
    """Draw a styled table from slide-source rows and column widths."""
    rows = [spec["columns"], *spec["rows"]]
    shape = slide.shapes.add_table(
        len(rows),
        len(rows[0]),
        Inches(0.55),
        Inches(1.65),
        Inches(12.2),
        Inches(spec.get("table_height", 3.8)),
    )
    target = shape.table
    if "widths" in spec:
        for column, width in zip(target.columns, spec["widths"], strict=True):
            column.width = Inches(width)
    for row_index, row in enumerate(rows):
        for column_index, value in enumerate(row):
            cell = target.cell(row_index, column_index)
            cell.text = str(value)
            cell.margin_left = Inches(0.12)
            cell.margin_top = Inches(0.1)
            cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor.from_string(
                NAVY if row_index == 0 else WHITE if row_index % 2 == 0 else PALE
            )
            for paragraph in cell.text_frame.paragraphs:
                paragraph.font.name = FONT
                paragraph.font.size = Pt(spec.get("font_size", 19))
                paragraph.font.bold = row_index == 0
                paragraph.font.color.rgb = RGBColor.from_string(
                    WHITE if row_index == 0 else INK
                )


def cards(slide: Any, spec: dict[str, Any]) -> None:
    """Draw up to four text cards in a balanced grid."""
    items = spec["cards"]
    columns = 2 if len(items) in {2, 4} else len(items)
    width = (12.2 - 0.3 * (columns - 1)) / columns
    height = 2.15 if len(items) == 4 else 4.2
    for index, item in enumerate(items):
        x = 0.55 + (index % columns) * (width + 0.3)
        y = 1.72 + (index // columns) * (height + 0.2)
        rect(slide, x, y, width, height, PALE)
        box(
            slide,
            item["title"],
            x + 0.2,
            y + 0.15,
            width - 0.4,
            0.5,
            size=24,
            color=TEAL,
            bold=True,
        )
        box(
            slide,
            item["text"],
            x + 0.2,
            y + 0.72,
            width - 0.4,
            height - 0.82,
            size=spec.get("font_size", 25),
        )


def process(slide: Any, spec: dict[str, Any]) -> None:
    """Draw editable process lanes with explicit labels and ordered steps."""
    lanes = spec["lanes"]
    top = 2.0 if spec.get("takeaway") else 1.65
    lane_gap = 0.18
    bottom = 6.05 if spec.get("caption") else 6.55
    lane_height = min(
        1.8,
        (bottom - top - lane_gap * (len(lanes) - 1)) / len(lanes),
    )
    for lane_index, lane in enumerate(lanes):
        y = top + lane_index * (lane_height + lane_gap)
        box(slide, lane["title"], 0.55, y, 1.9, 0.7, size=20, color=TEAL, bold=True)
        steps = lane["steps"]
        gap = 0.12
        width = (10.15 - gap * (len(steps) - 1)) / len(steps)
        for index, step in enumerate(steps):
            x = 2.6 + index * (width + gap)
            rect(slide, x, y, width, lane_height, PALE)
            box(
                slide,
                step,
                x + 0.08,
                y + 0.08,
                width - 0.16,
                lane_height - 0.16,
                size=spec.get("font_size", 19),
            )
            if index < len(steps) - 1:
                box(
                    slide,
                    "→",
                    x + width,
                    y + (lane_height - 0.4) / 2,
                    gap,
                    0.4,
                    size=22,
                    color=TEAL,
                    bold=True,
                )
    if spec.get("caption"):
        box(slide, spec["caption"], 0.65, 6.18, 12.0, 0.55, size=16, color=MUTED)


def placeholder(slide: Any, spec: dict[str, Any]) -> None:
    """Show an explicit placeholder for an unfinished scientific output."""
    rect(slide, 0.8, 1.8, 11.7, 4.45, PALE)
    box(
        slide,
        spec["status"].upper(),
        1.2,
        2.15,
        10.9,
        0.65,
        size=24,
        color=TEAL,
        bold=True,
    )
    box(
        slide,
        spec.get("message", "No result is available yet."),
        1.2,
        3.0,
        10.9,
        1.25,
        size=24,
        color=NAVY,
        bold=True,
    )
    if spec.get("detail"):
        box(slide, spec["detail"], 1.2, 4.45, 10.9, 1.25, size=22, color=MUTED)


def annotated_image(slide: Any, path: Path, labels: list[str]) -> None:
    """Place a figure beside editable, readable scientific callouts."""
    picture(slide, path, (0.55, 1.75, 7.35, 4.3))
    if len(labels) == 1:
        start, spacing, height = 2.7, 0.0, 1.6
    elif len(labels) == 2:
        start, spacing, height = 2.05, 1.9, 1.6
    else:
        start, spacing, height = 1.8, 1.5, 1.35
    for index, label in enumerate(labels):
        y = start + index * spacing
        rect(slide, 8.2, y, 4.5, height, PALE)
        box(
            slide,
            label,
            8.45,
            y + 0.12,
            4.0,
            height - 0.24,
            size=24,
            color=NAVY,
            bold=True,
        )
