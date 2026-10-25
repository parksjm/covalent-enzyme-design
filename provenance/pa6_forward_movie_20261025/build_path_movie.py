"""Assemble actual forward optimizer frames for interactive VMD inspection."""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
import struct
import subprocess
from pathlib import Path

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = ROOT / "output_nyl12_pa6_dftb_pilot_20261024"
SOURCE = BASE / "coordinate_driving_forward_reverse_attempt3"
OUTPUT = BASE / "overnight_review_20261025/forward_movie"
VIEWER = ROOT / "output_vmd_inspection_20261024"


def header(path: Path) -> tuple[int, int]:
    """Read standard little-endian DCD frame and atom counts."""
    with path.open("rb") as handle:
        record = struct.unpack("<i", handle.read(4))[0]
        assert record == 84, (path, record)
        data = handle.read(record)
        assert data[:4] == b"CORD"
        frames = struct.unpack("<i", data[4:8])[0]
        assert struct.unpack("<i", handle.read(4))[0] == record
        title_size = struct.unpack("<i", handle.read(4))[0]
        handle.seek(title_size + 4, 1)
        assert struct.unpack("<i", handle.read(4))[0] == 4
        atoms = struct.unpack("<i", handle.read(4))[0]
    return frames, atoms


def main() -> None:
    """Concatenate all retained frames and preserve source/frame provenance."""
    OUTPUT.mkdir()
    hashes = {}
    windows = []
    commands = [f"parm {VIEWER / 'stripped.parm7'}"]
    offset = 0
    for index in range(18):
        folder = SOURCE / f"forward_{index:02d}"
        paths = [
            folder / name
            for name in (
                "optimization_stripped.dcd",
                "display_frames.json",
                "process.json",
                "acceptance.json",
            )
        ]
        for path in paths:
            hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        trajectory, display_path, process_path, acceptance_path = paths
        display = json.loads(display_path.read_text())
        process = json.loads(process_path.read_text())
        acceptance = json.loads(acceptance_path.read_text())
        frames, atoms = header(trajectory)
        assert atoms == 19976
        assert display["exact_endpoint_added"]
        assert frames == display["saved_optimizer_frames"] + 1
        windows.append(
            {
                "window": index,
                "name": folder.name,
                "frame_count": frames,
                "vmd_first_frame": offset,
                "vmd_last_frame": offset + frames - 1,
                "exact_endpoint_frame": offset + frames - 1,
                "saved_optimizer_frames": display["saved_optimizer_frames"],
                "target_attack_A": process["target_attack_A"],
                "target_proton_A": process["target_proton_A"],
                "achieved_distances_A": acceptance["distances_A"],
                "RMS": acceptance["RMS"],
                "source": str(trajectory),
            }
        )
        offset += frames
        commands.append(f"trajin {trajectory} namdcell")
    commands.extend(["trajout forward_driven_path.dcd dcd namdcell", "run", "quit"])
    script = OUTPUT / "concatenate.cpptraj"
    script.write_text("\n".join(commands) + "\n")
    command = ["/data/jpc/code/ambertools26/venv/bin/cpptraj", "-i", str(script)]
    with (OUTPUT / "cpptraj.log").open("w") as handle:
        subprocess.run(
            command, cwd=OUTPUT, stdout=handle, stderr=subprocess.STDOUT, check=True
        )
    assert header(OUTPUT / "forward_driven_path.dcd") == (offset, 19976)
    for name in ("stripped.parm7.gz", "qm_indices.tcl", "nyl12_es_vmd_common.tcl"):
        path = VIEWER / name
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        shutil.copy2(path, OUTPUT / name)
    hashes[str(VIEWER / "stripped.parm7")] = hashlib.sha256(
        (VIEWER / "stripped.parm7").read_bytes()
    ).hexdigest()
    hashes[str(Path(__file__).resolve())] = hashlib.sha256(
        Path(__file__).read_bytes()
    ).hexdigest()
    shutil.copy2(__file__, OUTPUT / "build_path_movie.snapshot.py")
    loader = """# Actual forward optimizer states; FAILED TI1 construction, reverse NOT RUN.
set ::nyl12_dir [file dirname [file normalize [info script]]]
foreach name {qm_indices.tcl nyl12_es_vmd_common.tcl} {
    source [file join $::nyl12_dir $name]
}
nyl12_load stripped.parm7 forward_driven_path.dcd FRAMES $::nyl12_qm_indices $::nyl12_stripped_atom_count
puts "Forward-only failed TI1 construction: FRAMES optimizer frames, no reverse, no interpolation."
puts "Use VMD animation controls; initial, final, pocket, system are available."
""".replace("FRAMES", str(offset))
    (OUTPUT / "forward_driven_path.vmd").write_text(loader)
    (OUTPUT / "README.txt").write_text(
        "Open with: vmd -e forward_driven_path.vmd\n"
        "FAILED TI1 construction: forward scan only; reverse NOT RUN.\n"
        "These are optimizer states, not reaction time or a validated reaction path.\n"
        "All actual source frames retained in order, with no alignment/interpolation.\n"
        "ntwx=10 optimizer snapshots plus each exact endpoint; early convergence "
        "can yield endpoint-only windows.\n"
        "Frame numbering in window_to_frame_manifest.json is zero-based (VMD).\n"
        "Endpoint attack 2.60843 A vs target 1.55 A; the candidate TI1 gate failed.\n"
        "98 physical QM atoms shown as Licorice; full stripped selection has 19976 atoms.\n"
    )
    manifest = {
        "scientific_status": "FAILED TI1 construction; reverse NOT RUN",
        "frames": offset,
        "atoms": 19976,
        "frame_index_base": 0,
        "transforms": "none; no interpolation or alignment",
        "windows": windows,
        "source_sha256": hashes,
        "command": command,
        "process_acceptance": True,
        "scientific_acceptance": False,
    }
    (OUTPUT / "window_to_frame_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    logging.info("Verified %d frames and 19976 atoms; retained all 18 windows", offset)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
