"""Run one source-bound, bias-free plateau release with the frozen native code."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import sys
from typing import Any

LOGGER = logging.getLogger(__name__)


def verified_record(plan: dict[str, Any], image: int) -> dict[str, Any]:
    """Validate the shared sources and requested nonduplicate starting image.

    Args:
        plan: Retained preparation manifest.
        image: Zero-based image in the original combined PES plot.

    Returns:
        The verified starting-image record.

    Raises:
        ValueError: If scope, identity, reuse status or a content hash is invalid.
        FileExistsError: If native outputs already exist.
    """
    matches = [row for row in plan["images"] if row["image"] == image]
    if len(matches) != 1 or not 11 <= image <= 21:
        raise ValueError("Require one requested plateau image in 11–21")
    row = matches[0]
    if row["reuse_existing_release"]:
        raise ValueError("This image reuses an accepted release; do not rerun it")
    campaign = Path(plan["campaign"]).resolve(strict=True)
    output = Path(row["native_directory"]).resolve()
    if not output.is_relative_to(campaign) or output == campaign:
        raise ValueError("Native output must remain inside the original campaign")
    if output.exists():
        raise FileExistsError(f"Preserve existing native attempt: {output}")
    hashes = {**plan["source_sha256"], row["start_restart"]: row["start_sha256"]}
    for filename, expected in hashes.items():
        observed = hashlib.sha256(Path(filename).read_bytes()).hexdigest()
        if observed != expected:
            raise ValueError(f"Source/input identity changed: {filename}")
    return row


def run_release(plan_path: Path, image: int) -> dict[str, Any]:
    """Execute the unchanged native all-mobile-region minimizer for one image.

    Args:
        plan_path: Source-hashed preparation manifest.
        image: Original combined-plot image index.

    Returns:
        Native scientific acceptance, independently checked by the frozen code.

    Raises:
        RuntimeError: If native convergence or independent force acceptance fails.
    """
    plan = json.loads(plan_path.read_text())
    row = verified_record(plan, image)
    campaign = Path(plan["campaign"])
    source = campaign / "source_snapshot"
    sys.path.insert(0, str(source))
    sys.path.insert(0, str(source / "src"))
    from nyl12_covalent.nterminal_relay_native import run_native
    from nyl12_covalent.relay_release import minimization_input
    from nyl12_covalent.stationary_native import NativeState

    state = NativeState(campaign / "model")
    result = run_native(
        campaign,
        Path(row["native_directory"]),
        Path(row["start_restart"]),
        minimization_input(state.mdin_source, maxcycles=10000),
        restraints=None,
    )
    if result.get("accepted") is not True:
        raise RuntimeError("Release failed native/independent force acceptance")
    LOGGER.info(
        "Image %s release accepted; inspect chemical retention separately", image
    )
    return result


def main() -> None:
    """Run one explicitly prepared plateau image without reaction restraints."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--image", type=int, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    run_release(args.plan, args.image)


if __name__ == "__main__":
    main()
