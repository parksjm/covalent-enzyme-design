"""Write lightweight SVG and PNG relaxed-scan energy profiles."""

from __future__ import annotations

import shutil
import subprocess
import math
from pathlib import Path


def write_profile_plot(
    output_dir: Path,
    label: str,
    angles: dict[str, list[float]],
    relative_energies: dict[str, list[float]],
) -> None:
    """Write a forward/reverse SVG and PNG plot using ImageMagick.

    Args:
        output_dir: Destination directory.
        label: Scan name displayed in the plot title.
        angles: Optimized dihedral angles, in degrees, for each direction.
        relative_energies: Energies in kcal/mol for each direction.

    Raises:
        RuntimeError: If ImageMagick is unavailable or conversion fails.
        ValueError: If the profiles are empty or their arrays differ in length.
    """
    if not angles or set(angles) != set(relative_energies):
        raise ValueError("Angle and energy directions must be identical and nonempty")
    for direction in angles:
        if not angles[direction] or len(angles[direction]) != len(
            relative_energies[direction]
        ):
            raise ValueError(f"Mismatched or empty profile: {direction}")
    unwrapped_angles = {
        direction: unwrap_angles_degrees(values) for direction, values in angles.items()
    }
    converter = shutil.which("convert")
    if converter is None:
        raise RuntimeError("ImageMagick convert is required to create scan PNG plots")
    width, height = 900, 560
    left, right, top, bottom = 90, 30, 55, 70
    x_values = [value for series in unwrapped_angles.values() for value in series]
    y_values = [value for series in relative_energies.values() for value in series]
    min_x, max_x = min(x_values), max(x_values)
    min_y, max_y = min(0.0, min(y_values)), max(y_values)
    max_y = max_y if max_y > min_y else min_y + 1.0
    plot_width, plot_height = width - left - right, height - top - bottom

    def point_xy(angle: float, energy: float) -> tuple[float, float]:
        x = left + (angle - min_x) * plot_width / (max_x - min_x or 1.0)
        y = top + (max_y - energy) * plot_height / (max_y - min_y)
        return x, y

    colors = {"forward": "#1769aa", "reverse": "#d95f02"}
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="30" text-anchor="middle" font-size="20" font-family="sans-serif">{_xml_escape(label)} relaxed scan</text>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height - bottom}" stroke="black"/>',
        f'<line x1="{left}" y1="{height - bottom}" x2="{width - right}" y2="{height - bottom}" stroke="black"/>',
        f'<text x="{width / 2}" y="{height - 18}" text-anchor="middle" font-size="15" font-family="sans-serif">Continuous unwrapped dihedral (degrees)</text>',
        f'<text transform="translate(20 {height / 2}) rotate(-90)" text-anchor="middle" font-size="15" font-family="sans-serif">Relative electronic energy (kcal/mol)</text>',
    ]
    for direction in angles:
        coords = [
            point_xy(angle, energy)
            for angle, energy in zip(
                unwrapped_angles[direction], relative_energies[direction], strict=True
            )
        ]
        points_text = " ".join(f"{x:.2f},{y:.2f}" for x, y in coords)
        lines.append(
            f'<polyline points="{points_text}" fill="none" stroke="{colors[direction]}" stroke-width="2.5"/>'
        )
        lines.extend(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.5" fill="{colors[direction]}"/>'
            for x, y in coords
        )
    lines.extend(
        [
            f'<text x="{width - right - 130}" y="{top + 8}" fill="{colors["forward"]}" font-size="14" font-family="sans-serif">Forward</text>',
            f'<text x="{width - right - 130}" y="{top + 28}" fill="{colors["reverse"]}" font-size="14" font-family="sans-serif">Reverse</text>',
            "</svg>",
        ]
    )
    svg_path = output_dir / "energy_profile.svg"
    png_path = output_dir / "energy_profile.png"
    svg_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    result = subprocess.run(
        [converter, str(svg_path), str(png_path)],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not png_path.is_file():
        raise RuntimeError(
            f"Could not convert scan plot to PNG: {result.stderr.strip()}"
        )


def unwrap_angles_degrees(angles: list[float]) -> list[float]:
    """Unwrap a scan's signed angles continuously in their existing order.

    Args:
        angles: Ordered signed dihedral values in degrees.

    Returns:
        Continuous values adjusted by multiples of 360 degrees.

    Raises:
        ValueError: If any angle is nonfinite or outside ``[-180, 180]``.
    """
    if any(
        not math.isfinite(angle) or not -180.0 <= angle <= 180.0 for angle in angles
    ):
        raise ValueError("Raw dihedral angles must be finite and in [-180, 180]")
    if not angles:
        return []
    result = [float(angles[0])]
    for previous_raw, current_raw in zip(angles, angles[1:], strict=False):
        delta = (current_raw - previous_raw + 180.0) % 360.0 - 180.0
        result.append(result[-1] + delta)
    return result


def _xml_escape(value: str) -> str:
    """Escape text inserted into SVG content."""
    return value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
