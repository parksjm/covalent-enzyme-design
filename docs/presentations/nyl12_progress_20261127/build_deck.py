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
from pptx.util import Inches

from scripts.presentation_layouts import NAVY, add_cover, add_footer, build_date
from scripts.presentation_shapes import (
    box as _box,
    cards as _cards,
    annotated_image as _annotated_image,
    picture as _picture,
    placeholder as _placeholder,
    process as _process,
    rect as _rect,
    table as _table,
)

LOGGER = logging.getLogger(__name__)
TEAL = "008C95"
INK = "223548"
MUTED = "526577"


def validate_source(source: dict[str, Any], base: Path) -> None:
    """Reject incomplete slides and unreadable figure inputs before building.

    Args:
        source: Reviewed slide specification.
        base: Directory containing the specification and relative assets.
    """
    slides = source.get("slides", [])
    if not slides:
        raise ValueError("The deck must contain slides")
    supported = {
        "title",
        "text",
        "image",
        "annotated_image",
        "split",
        "table",
        "cards",
        "snapshots",
        "process",
        "placeholder",
        "two_images",
        "four_images",
    }
    for index, spec in enumerate(slides, 1):
        if not spec.get("title") or spec.get("kind") not in supported:
            raise ValueError(f"Invalid slide {index}")
        if spec["kind"] in {"image", "annotated_image", "split"} or (
            spec["kind"] == "title" and spec.get("image")
        ):
            path = base / spec["image"]
            with Image.open(path) as figure:
                figure.verify()
        if spec["kind"] == "annotated_image" and not spec.get("callouts"):
            raise ValueError(f"Annotated figure on slide {index} needs callouts")
        if spec["kind"] == "table":
            rows = [spec["columns"], *spec["rows"]]
            if not rows[0] or any(len(row) != len(rows[0]) for row in rows):
                raise ValueError(f"Ragged table on slide {index}")
        if spec["kind"] in {"cards", "snapshots"}:
            if not 1 <= len(spec.get("cards", [])) <= 4:
                raise ValueError(f"Invalid card count on slide {index}")
        if spec["kind"] == "process":
            lanes = spec.get("lanes", [])
            if not lanes or any(not lane.get("steps") for lane in lanes):
                raise ValueError(f"Invalid process lanes on slide {index}")
        if spec["kind"] == "placeholder" and not spec.get("status"):
            raise ValueError(f"Placeholder slide {index} needs an explicit status")


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
    today = build_date()
    date_label = today.strftime("%d %B %Y").lstrip("0")
    presentation = Presentation()
    presentation.slide_width = Inches(13.333)
    presentation.slide_height = Inches(7.5)
    assets = {}
    for index, spec in enumerate(source["slides"], 1):
        slide = presentation.slides.add_slide(presentation.slide_layouts[6])
        if spec["kind"] == "title":
            cover_image = (
                str(source_path.parent / spec["image"]) if spec.get("image") else None
            )
            if cover_image:
                assets[cover_image] = hashlib.sha256(
                    Path(cover_image).read_bytes()
                ).hexdigest()
            add_cover(
                slide, spec["title"], spec.get("subtitle", ""), date_label, cover_image
            )
            notes = [spec.get("notes", ""), "Sources:", *spec.get("sources", [])]
            slide.notes_slide.notes_text_frame.text = "\n".join(notes)
            continue
        _rect(slide, 0, 0, 13.333, 0.12, TEAL)
        _box(
            slide,
            spec["title"],
            0.55,
            0.28,
            12.2,
            1.05,
            size=spec.get("title_size", 34),
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
                size=24,
                color=TEAL,
                bold=True,
            )
        kind = spec["kind"]
        if kind in {"title", "text"}:
            text = "\n".join(spec.get("bullets", []))
            _box(slide, text, 0.65, 1.95, 11.9, 4.55, size=spec.get("font_size", 28))
        elif kind in {"image", "split", "annotated_image"}:
            path = source_path.parent / spec["image"]
            area = (0.55, 1.75, 12.2, 4.3)
            if kind == "annotated_image":
                _annotated_image(slide, path, spec["callouts"])
                area = None
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
                area = (5.35, 1.75, 7.4, 4.3)
            if area is not None:
                _picture(slide, path, area)
            assets[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        elif kind in {"two_images", "four_images"}:
            positions = (
                [(0.55, 1.85, 6.0, 3.9), (6.8, 1.85, 6.0, 3.9)]
                if kind == "two_images"
                else [(0.55 + i * 3.13, 1.95, 2.95, 4.1) for i in range(4)]
            )
            for relative, label, area in zip(
                spec["images"], spec["labels"], positions, strict=True
            ):
                path = source_path.parent / relative
                if kind == "four_images":
                    left, top, right, bottom = spec["matched_crop_box"]
                    ratio = (right - left) / (bottom - top)
                    width, height = (
                        min(area[2], area[3] * ratio),
                        min(area[3], area[2] / ratio),
                    )
                    picture = slide.shapes.add_picture(
                        str(path),
                        Inches(area[0] + (area[2] - width) / 2),
                        Inches(area[1]),
                        width=Inches(width),
                        height=Inches(height),
                    )
                    picture.crop_left, picture.crop_right = (
                        left / 1600,
                        1 - right / 1600,
                    )
                    picture.crop_top, picture.crop_bottom = (
                        top / 1400,
                        1 - bottom / 1400,
                    )
                else:
                    _picture(slide, path, area)
                _box(
                    slide,
                    label,
                    area[0],
                    area[1] - 0.28,
                    area[2],
                    0.4,
                    size=24,
                    color=TEAL,
                    bold=True,
                )
                assets[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
            _box(
                slide,
                spec.get("caption", ""),
                0.65,
                6.35,
                12,
                0.5,
                size=20,
                color=MUTED,
            )
        elif kind == "table":
            _table(slide, spec)
            _box(
                slide,
                spec.get("caption", ""),
                0.65,
                6.0,
                12.0,
                0.65,
                size=spec.get("caption_size", 20),
                color=MUTED,
            )
        elif kind in {"cards", "snapshots"}:
            _cards(slide, spec)
        elif kind == "process":
            _process(slide, spec)
        elif kind == "placeholder":
            _placeholder(slide, spec)
        add_footer(
            slide,
            spec.get("footer", "Nyl12–PA6 | Research progress"),
            index,
            date_label,
        )
        if kind in {"image", "split", "annotated_image"} and spec.get("caption"):
            _box(
                slide,
                spec["caption"],
                0.65,
                6.18,
                12.0,
                0.62,
                size=spec.get("caption_size", 20),
                color=MUTED,
            )
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
        "build_date": today.date().isoformat(),
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
