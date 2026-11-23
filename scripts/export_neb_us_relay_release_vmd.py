"""Export mapped, stripped native relay optimizer snapshots for VMD."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
from netCDF4 import Dataset  # noqa: E402
from parmed.amber import AmberParm, Rst7  # noqa: E402
from parmed.amber.mask import AmberMask  # noqa: E402
from nyl12_covalent.relay_release import ROLES  # noqa: E402
from nyl12_covalent.stationary_native import digest, write_json  # noqa: E402
from nyl12_covalent.vmd_export_io import validate_dcd  # noqa: E402

STRIP_MASK = ":WAT,Na+,Cl-,NA,CL"


def export(seed: Path, output: Path) -> None:
    """Snapshot available saved iterations and verify their mapped DCD export.

    Args:
        seed: Prepared/running release seed directory.
        output: Fresh immutable destination for this visualization checkpoint.
    """
    output.mkdir(parents=True, exist_ok=False)
    source = seed / "optimization.nc"
    snapshot = output / "source_snapshot"
    snapshot.mkdir()
    original_source = source
    with (
        Dataset(original_source) as incoming,
        Dataset(snapshot / source.name, "w", format="NETCDF3_64BIT_OFFSET") as outgoing,
    ):
        original_count = len(incoming.dimensions["frame"])
        if original_count < 1:
            raise ValueError("Native trajectory has no saved optimizer frame")
        original_selected = list(
            range(0, original_count, max(1, (original_count + 46) // 47))
        )
        outgoing.setncatts(incoming.__dict__)
        for key, dimension in incoming.dimensions.items():
            outgoing.createDimension(
                key, None if dimension.isunlimited() else len(dimension)
            )
        for key, variable in incoming.variables.items():
            target = outgoing.createVariable(key, variable.dtype, variable.dimensions)
            target.setncatts(variable.__dict__)
            if "frame" in variable.dimensions:
                if variable.dimensions[0] != "frame":
                    raise ValueError("Unsupported frame dimension order")
                for new, old in enumerate(original_selected):
                    target[new] = variable[old]
            else:
                target[:] = variable[:]
    shutil.copy2(Path(__file__), snapshot / Path(__file__).name)
    source = snapshot / source.name
    topology = seed / "inputs/system.prmtop"
    initial = seed / "inputs/optimized.ncrst"
    parm = AmberParm(str(topology))
    kept = np.flatnonzero(
        ~np.asarray(AmberMask(parm, STRIP_MASK).Selection(), dtype=bool)
    )
    remap = {int(original): new for new, original in enumerate(kept)}
    restart = Rst7.open(str(initial))
    xyz = np.asarray(restart.coordinates).reshape((-1, 3))
    states = [(xyz[kept], np.asarray(restart.box))]
    with Dataset(source) as nc:
        count = len(nc.dimensions["frame"])
        if count < 1:
            raise ValueError("Native trajectory has no saved optimizer frame")
        stride = max(1, (count + 46) // 47)
        selected = list(range(0, count, stride))
        for i in selected:
            states.append(
                (
                    np.asarray(nc["coordinates"][i, kept, :]),
                    np.r_[nc["cell_lengths"][i], nc["cell_angles"][i]],
                )
            )
    script = output / "export.cpptraj"
    script.write_text(
        f"parm {topology}\ntrajin {initial}\n"
        f"trajin {source} 1 {count} {stride}\nstrip {STRIP_MASK}\n"
        f"trajout {output / 'optimization.dcd'} dcd namdcell\nrun\n"
        "quit\n"
    )
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = "/data/jpc/code/ambertools26/venv/lib"
    with (output / "cpptraj.log").open("w") as log:
        subprocess.run(
            ["/data/jpc/code/ambertools26/venv/bin/cpptraj", "-i", str(script)],
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            check=True,
        )
    topology_script = output / "topology.cpptraj"
    topology_script.write_text(
        f"parm {topology}\nparmstrip {STRIP_MASK}\n"
        f"parmwrite out {output / 'stripped.parm7'}\nquit\n"
    )
    with (output / "topology_cpptraj.log").open("w") as log:
        subprocess.run(
            [
                "/data/jpc/code/ambertools26/venv/bin/cpptraj",
                "-i",
                str(topology_script),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            check=True,
        )
    validation = validate_dcd(output / "optimization.dcd", states)
    mapped = {key: remap[value - 1] for key, value in ROLES.items()}
    lines = [
        "set here [file dirname [file normalize [info script]]]",
        "mol new [file join $here stripped.parm7] type parm7 waitfor all",
        "mol addfile [file join $here optimization.dcd] type dcd waitfor all",
        "mol delrep 0 top",
        "mol representation NewCartoon; mol selection protein; mol color Chain; mol addrep top",
        "mol representation Licorice; mol selection {resname LIG or resid 929}; mol color Element; mol addrep top",
    ]
    lines += [
        "proc relay_remove_bond {a b} {",
        "  foreach {first second} [list $a $b $b $a] {",
        '    set selected [atomselect top "index $first"]',
        "    set neighbors [lindex [$selected getbonds] 0]",
        "    set position [lsearch -exact $neighbors $second]",
        "    if {$position >= 0} {$selected setbonds [list [lreplace $neighbors $position $position]]}",
        "    $selected delete",
        "  }",
        "}",
    ]
    for a, b in [("HG", "OG"), ("H1", "N"), ("H2", "N"), ("HN", "NL"), ("C", "NL")]:
        lines.append(f"relay_remove_bond {mapped[a]} {mapped[b]}")
    alternatives = [
        (h, a, 1.2) for h in ("H1", "H2", "HG", "HN") for a in ("N", "OG", "NL", "O")
    ]
    alternatives += [("C", "OG", 1.85), ("C", "NL", 1.85)]
    for a, b, cutoff in alternatives:
        lines.append(
            f"mol representation DynamicBonds {cutoff} 0.12 12; mol selection {{index {mapped[a]} {mapped[b]}}}; mol color Element; mol addrep top"
        )
    for a, b in [
        ("C", "OG"),
        ("C", "NL"),
        ("C", "O"),
        ("HG", "N"),
        ("H1", "N"),
        ("H2", "N"),
    ]:
        lines.append(
            f"label add Bonds [molinfo top]/{mapped[a]} [molinfo top]/{mapped[b]}"
        )
    lines.append(
        'puts "Frame0: exact sampled seed. Later frames: optimizer iterations, not MD time or reaction path."'
    )
    (output / "view_relay.vmd").write_text("\n".join(lines) + "\n")
    write_json(
        output / "manifest.json",
        {
            "seed": str(seed),
            "frame_count": len(states),
            "source_optimizer_frames_0based": original_selected,
            "original_source_frame_count_at_snapshot": original_count,
            "stripped_to_original_atom_indices_0based": kept.tolist(),
            "role_indices_stripped_0based": mapped,
            "validation": validation,
            "topology_sha256": digest(topology),
            "source_snapshot_sha256": digest(source),
            "seed_coordinate_sha256": digest(initial),
            "strip_mask": STRIP_MASK,
            "interpretation": "Optimizer states; not MD time, a timed reaction path, or evidence of convergence",
        },
    )
    (output / "README.md").write_text(
        "Source view_relay.vmd in VMD. Frame0 is the exact trajectory seed; later frames are saved native optimizer iterations. Numerical/scientific acceptance is separate.\n"
    )
    (output / "MANIFEST.sha256").write_text(
        "".join(
            f"{digest(p)}  {p.name}\n" for p in sorted(output.iterdir()) if p.is_file()
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    export(args.seed.resolve(), args.output.resolve())
