"""Shared visual layouts and date formatting for scientific presentations."""

from __future__ import annotations

from datetime import datetime
from zoneinfo import ZoneInfo

from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.util import Inches, Pt

NAVY = "17324D"
TEAL = "008C95"
INK = "223548"
MUTED = "526577"
PALE = "EFF5F8"
WHITE = "FFFFFF"
FONT = "DejaVu Sans"


def build_date() -> datetime:
    """Return the local build date in the project's America/New_York timezone."""
    return datetime.now(ZoneInfo("America/New_York"))


def add_cover(
    slide: object,
    title: str,
    subtitle: str,
    date_label: str,
    cover_image: str | None,
) -> None:
    """Draw a distinct title slide with date and optional scientific image."""
    panel = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, Inches(0), Inches(0), Inches(13.333), Inches(7.5)
    )
    panel.fill.solid()
    panel.fill.fore_color.rgb = RGBColor.from_string(NAVY)
    panel.line.fill.background()
    text_width = 7.2 if cover_image else 11.3
    _cover_text(slide, title, 0.8, 1.75, text_width, 2.0, 34, WHITE, True)
    _cover_text(slide, subtitle, 0.82, 4.0, text_width - 0.4, 1.05, 22, "DCE7ED")
    _cover_text(slide, date_label, 0.82, 6.55, text_width, 0.35, 16, "DCE7ED")
    if cover_image:
        from PIL import Image

        with Image.open(cover_image) as image:
            ratio = image.width / image.height
        width, height = 5.1, 5.1 / ratio
        if height > 5.4:
            width, height = 5.4 * ratio, 5.4
        slide.shapes.add_picture(
            cover_image,
            Inches(8.0 + (5.0 - width) / 2),
            Inches(1.05 + (5.4 - height) / 2),
            width=Inches(width),
            height=Inches(height),
        )


def add_footer(slide: object, label: str, slide_number: int, date_label: str) -> None:
    """Draw a consistent footer with project, current date, and slide number."""
    _cover_text(
        slide, f"{label}  |  {date_label}", 0.55, 7.05, 11.6, 0.25, 10, "526577"
    )
    _cover_text(slide, str(slide_number), 12.35, 7.05, 0.4, 0.25, 10, "526577")


def _cover_text(
    slide: object,
    text: str,
    x: float,
    y: float,
    width: float,
    height: float,
    size: float,
    color: str,
    bold: bool = False,
) -> None:
    """Add one formatted text box used by cover and footer layouts."""
    shape = slide.shapes.add_textbox(
        Inches(x), Inches(y), Inches(width), Inches(height)
    )
    frame = shape.text_frame
    frame.clear()
    frame.word_wrap = True
    paragraph = frame.paragraphs[0]
    paragraph.text = text
    paragraph.font.name = FONT
    paragraph.font.size = Pt(size)
    paragraph.font.bold = bold
    paragraph.font.color.rgb = RGBColor.from_string(color)
