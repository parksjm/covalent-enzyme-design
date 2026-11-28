"""Create a compact contact sheet from numbered presentation slide renders."""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

from PIL import Image, ImageOps

LOGGER = logging.getLogger(__name__)


def make_contact_sheet(source_dir: Path, destination: Path) -> int:
    """Arrange numbered slide PNGs in a 4-column review sheet.

    Args:
        source_dir: Directory containing slide-NN.png files.
        destination: Contact sheet PNG output path.

    Returns:
        Number of slides placed on the sheet.
    """
    slides = sorted(
        source_dir.glob("slide-*.png"),
        key=lambda path: int(path.stem.rsplit("-", maxsplit=1)[1]),
    )
    if not slides:
        raise ValueError(f"No numbered slide renders found in {source_dir}")
    columns = 4
    rows = math.ceil(len(slides) / columns)
    cell_width, cell_height, margin = 480, 270, 12
    canvas = Image.new(
        "RGB",
        (
            columns * (cell_width + margin) + margin,
            rows * (cell_height + margin) + margin,
        ),
        "#E8EEF2",
    )
    for index, path in enumerate(slides):
        with Image.open(path) as source:
            thumbnail = ImageOps.fit(
                source.convert("RGB"),
                (cell_width, cell_height),
                method=Image.Resampling.LANCZOS,
            )
        x = margin + (index % columns) * (cell_width + margin)
        y = margin + (index // columns) * (cell_height + margin)
        canvas.paste(thumbnail, (x, y))
    destination.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(destination, format="PNG", optimize=True)
    return len(slides)


def main() -> None:
    """Parse paths and write a slide contact sheet."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    count = make_contact_sheet(args.source_dir, args.output)
    LOGGER.info("Created contact sheet for %d slides: %s", count, args.output)


if __name__ == "__main__":
    main()
