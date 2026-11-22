"""Refresh the viewer while retaining the original bundle for rollback."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import logging
import shutil
from pathlib import Path

from parmed.amber import AmberParm

LOGGER = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Return the SHA256 of a local file."""
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main() -> None:
    """Validate the existing bundle and regenerate only its VMD viewer."""
    repair = Path(__file__).resolve().parent
    bundle = repair.parent
    checksum_path = bundle / "MANIFEST.sha256"
    entries = [line.split(maxsplit=1) for line in checksum_path.read_text().splitlines()]
    for expected, name in entries:
        if digest(bundle / name) != expected:
            raise ValueError(f"Bundle checksum mismatch: {name}")
    metadata = json.loads((bundle / "manifest.json").read_text())
    config = json.loads((bundle.parent / "prepared.json").read_text())
    topology_path = Path(config["topology"])
    if digest(topology_path) != metadata["source_sha256"][str(topology_path)]:
        raise ValueError("Original topology hash changed")
    spec = importlib.util.spec_from_file_location("repaired_movie_vmd", repair / "neb_us_movie_vmd.py")
    if spec is None or spec.loader is None:
        raise ImportError("Cannot load retained viewer renderer")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    script = module.render_movie_vmd(
        AmberParm(str(topology_path)),
        metadata["original_atom_ids_1based"],
        metadata["qm_atom_ids_1based"],
        metadata["roles_1based"],
        metadata["frames"],
    )
    original = repair / "original"
    original.mkdir(exist_ok=False)
    for name in ("view_umbrella.vmd", "MANIFEST.sha256"):
        shutil.copy2(bundle / name, original / name)
    viewer = bundle / "view_umbrella.vmd"
    viewer.write_text(script)
    receipt = {
        "change": "Explicit frame zero, valid two-argument VDW, fail-fast loader",
        "original_viewer_sha256": digest(original / viewer.name),
        "viewer_sha256": digest(viewer),
        "frames": metadata["frames"],
        "atoms": metadata["atoms"],
        "unchanged_data_sha256": {
            name: digest(bundle / name) for name in ("path.parm7", "path.dcd")
        },
        "GUI_VMD_test_performed": False,
        "rollback": "Restore view_umbrella.vmd and MANIFEST.sha256 from viewer_repair/original",
    }
    (repair / "repair_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    names = [name for _, name in entries] + [
        "viewer_repair/repair_receipt.json", "viewer_repair/refresh_viewer.py",
        "viewer_repair/neb_us_movie_vmd.py",
    ]
    checksum_path.write_text("".join(f"{digest(bundle / name)}  {name}\n" for name in names))
    LOGGER.info("Refreshed viewer: %s (%s frames); trajectory unchanged", viewer, metadata["frames"])


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
