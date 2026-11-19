"""Export four ASH272 QM/MM optimization trajectories for VMD."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "qmmm/pa6_flipped_model30/06_ash272_comparison"
CASES = ("replica_1_tier_A", "replica_1_tier_B", "replica_2_tier_A", "replica_2_tier_B")
LOGGER = logging.getLogger(__name__)


def sha(path: Path) -> str:
    """Return SHA-256 of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def indices(mask: str) -> list[int]:
    """Expand one-based Amber atom mask to zero-based VMD indices."""
    if not mask.startswith("@"):
        raise ValueError("Expected atom-index mask")
    result = []
    for item in mask[1:].split(","):
        bounds = list(map(int, item.split("-")))
        result.extend(range(bounds[0] - 1, bounds[-1]))
    if len(result) != len(set(result)):
        raise ValueError("Duplicate QM index")
    return result


def main() -> None:
    """Export compact DCD/topology bundles from the completed cases."""
    logging.basicConfig(level=logging.INFO)
    topology = BASE / "00_system/ash272_solvated.prmtop"
    for case in CASES:
        source = BASE / "06_qmmm_optimizations" / case
        calculation = source / ("attempt_01/calculation" if "tier_B" in case else "")
        prepared = json.loads((source / "preparation.json").read_text())
        qm = indices(prepared["qm_region"]["remapped_mask"])
        if len(qm) != 99 or prepared["qm_region"]["charge"] != -1:
            raise ValueError(f"Wrong QM inventory: {case}")
        initial = Path(prepared["starting_mm_endpoint"])
        trajectory = calculation / "optimization.nc"
        final = calculation / "optimized.ncrst"
        output = source / "vmd_optimization"
        output.mkdir(exist_ok=False)
        cpp_input = output / "export.cpptraj"
        cpp_input.write_text(
            f"parm {topology}\ntrajin {initial}\ntrajin {trajectory}\n"
            f"trajin {final}\nstrip :WAT,Na+,Cl-,NA,CL\n"
            f"parmwrite out {output / 'stripped.parm7'}\n"
            f"trajout {output / 'optimization.dcd'} dcd\nrun\n"
        )
        with (output / "cpptraj.log").open("w") as log:
            subprocess.run(
                ["/data/jpc/code/ambertools26/venv/bin/cpptraj", "-i", str(cpp_input)],
                cwd=ROOT,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        selected = " ".join(map(str, qm))
        (
            output / "view_optimization.vmd"
        ).write_text(f'''set here [file dirname [file normalize [info script]]]
mol new [file join $here stripped.parm7] type parm7 waitfor all
mol addfile [file join $here optimization.dcd] type dcd waitfor all
set m [molinfo top]
mol delrep 0 $m
mol representation NewCartoon 0.3 10.0 4.1 0
mol color ColorID 7
mol selection "protein"
mol addrep $m
mol showrep $m 0 0
mol representation Licorice 0.20 12.0 12.0
mol color Element
mol selection "index {selected}"
mol addrep $m
mol representation Licorice 0.18 12.0 12.0
mol color Element
mol selection "residue 1312 to 1315 and not element C and not (index {selected})"
mol addrep $m
mol representation Licorice 0.18 12.0 12.0
mol color ColorID 3
mol selection "residue 1312 to 1315 and element C and not (index {selected})"
mol addrep $m
color Display Background white
animate goto 0
puts "{case}: QM charge -1; initial, optimizer states, final restart."
''')
        paths = (topology, initial, trajectory, final, source / "preparation.json")
        (output / "source_hashes.json").write_text(
            json.dumps({str(path): sha(path) for path in paths}, indent=2) + "\n"
        )
        (output / "MANIFEST.sha256").write_text(
            "".join(
                f"{sha(path)}  {path.name}\n"
                for path in sorted(output.iterdir())
                if path.is_file()
            )
        )
        LOGGER.info("Exported %s", case)


if __name__ == "__main__":
    main()
