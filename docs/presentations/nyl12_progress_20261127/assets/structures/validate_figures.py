"""Audit molecular PNG delivery and supplement the retained provenance."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from PIL import Image, ImageChops


def main() -> None:
    """Validate white opaque renders and save exact artifact hashes."""
    here = Path(__file__).resolve().parent
    acceptance_path = here / "figure_acceptance.json"
    acceptance = json.loads(acceptance_path.read_text())
    for name in [
        "af3_before",
        "af3_after",
        "boltz_before",
        "boltz_after",
        "ts_like",
        "assembly_loops",
        "loops_pocket",
    ]:
        path = here / f"{name}.png"
        im = Image.open(path).convert("RGBA")
        corners = [
            im.getpixel(x)
            for x in [
                (0, 0),
                (im.width - 1, 0),
                (0, im.height - 1),
                (im.width - 1, im.height - 1),
            ]
        ]
        assert corners[0] == (255, 255, 255, 255), (path, corners)
        if name != "loops_pocket":
            assert all(x == (255, 255, 255, 255) for x in corners), (path, corners)
        assert im.getchannel("A").getextrema() == (255, 255), path
        box = ImageChops.difference(
            im.convert("RGB"), Image.new("RGB", im.size, "white")
        ).getbbox()
        assert box is not None
        if name != "loops_pocket":
            assert box[0] > 0 and box[1] > 0
            assert box[2] < im.width and box[3] < im.height, (path, box)
        acceptance["figures"][name] = {
            "size": im.size,
            "foreground_bbox": box,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "exposed_background_white_and_image_opaque": True,
            "no_foreground_clipping": name != "loops_pocket",
            "crop_scope": "surrounding tetramer ribbons only; loops and ligand intact"
            if name == "loops_pocket"
            else "none",
        }
    acceptance["checks"] = {
        "ruff": "passed for current render scripts",
        "pytest": "not run; figure-only changes",
        "coverage": "not measured; rendering scripts",
    }
    acceptance_path.write_text(json.dumps(acceptance, indent=2) + "\n")
    provenance_path = here / "provenance.json"
    provenance = json.loads(provenance_path.read_text())
    provenance["figures"]["assembly_loops"]["loops"] = {
        "alpha_D_95_103": "deep blue",
        "beta_F_289_300": "orange-red",
    }
    provenance["figures"]["assembly_loops"]["occupied_site"] = "A/E ligand I"
    provenance["figures"]["loops_pocket"] = {
        "source": "ts_like",
        "parent_session": "assembly_loops.pse",
        "camera_method": "zoom actual D95-103/F289-300 loops and L4 with4A buffer",
        "context": "full tetramer retained, cartoon transparency0.55",
        "loops": "opaque separate tube objects, radius0.5A",
    }
    provenance["scripts"] = {
        name: hashlib.sha256((here / name).read_bytes()).hexdigest()
        for name in ["render_structures.py", "loops_pocket.py", "validate_figures.py"]
    }
    provenance["runtime"] = (
        "/data/jpc/code/miniforge3/envs/nyl12-pymol/bin/python; PyMOL3.1; two CPU threads; no GPUs"
    )
    provenance["reproduce"] = [
        "Run render_structures.py with documented PyMOL interpreter",
        "Run loops_pocket.py with same interpreter",
        "Run validate_figures.py with existing Pillow-enabled Python",
    ]
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
