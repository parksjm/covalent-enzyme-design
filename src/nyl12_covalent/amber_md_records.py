"""Separate instantaneous Amber MD records from average/fluctuation summaries."""

from __future__ import annotations
import re

_STATISTICS = re.compile(
    r"^\s*(?:A\s*V\s*E\s*R\s*A\s*G\s*E\s*S|"
    r"R\s*M\s*S\s+F\s*L\s*U\s*C\s*T\s*U\s*A\s*T\s*I\s*O\s*N\s*S)\b",
    re.IGNORECASE | re.MULTILINE,
)


def instantaneous_md_text(text: str) -> str:
    """Return native instantaneous records before any statistics-only section.

    Args:
        text: Complete native Amber MD output.

    Returns:
        Prefix excluding averages and RMS fluctuations, regardless of spacing.
    """
    return _STATISTICS.split(text, maxsplit=1)[0]
