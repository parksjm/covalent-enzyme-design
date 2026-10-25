"""Export a bounded-size snapshot of completed Asp272 optimizer windows."""
from __future__ import annotations

import hashlib
import json
import logging
import math
import runpy
import shutil
import subprocess
from pathlib import Path

ROOT = Path('/data/jpc/covalent-enzyme-design')
BASE = ROOT / 'output_nyl12_pa6_dftb_pilot_20261024'
RUN = BASE / 'asp272_relay_attempt2'
OUT = RUN / 'vmd_progress_snapshot_20261025'
VIEW = ROOT / 'output_vmd_inspection_20261024'
HEADER_SOURCE = ROOT / 'provenance/pa6_forward_movie_20261025/build_path_movie.py'
header = runpy.run_path(str(HEADER_SOURCE))['header']


def main() -> None:
    """Verify source sizes, concatenate sampled evaluations, and publish viewer."""
    OUT.mkdir()
    status = json.loads((RUN / 'status.json').read_text())
    labels = status['completed_windows']
    if not labels:
        raise ValueError('No completed windows')
    sources = [RUN / label / 'optimization_stripped.dcd' for label in labels]
    counts = [header(path)[0] for path in sources]
    stride = max(1, math.ceil(sum(counts) / 280))
    commands = [f"parm {VIEW / 'stripped.parm7'}"]
    records, hashes = [], {}
    offset = 0
    for label, path, count in zip(labels, sources, counts, strict=True):
        if header(path)[1] != 19976:
            raise ValueError(f'Wrong atom count: {path}')
        folder = RUN / label
        display = json.loads((folder / 'display_frames.json').read_text())
        accept = json.loads((folder / 'acceptance.json').read_text())
        if count != display['saved_evaluation_frames'] + 1:
            raise ValueError(f'Frame mismatch: {label}')
        if not accept['native_converged']:
            raise ValueError(f'Unconverged: {label}')
        frames = list(range(1, count + 1, stride))
        commands.append(f'trajin {path} 1 {count} {stride} namdcell')
        if frames[-1] != count:
            frames.append(count)
            commands.append(f'trajin {path} {count} {count} namdcell')
        records.append({'window': label, 'first_vmd_frame': offset,
                        'endpoint_vmd_frame': offset + len(frames) - 1,
                        'source_frames_one_based': frames,
                        'distances_A': accept['distances_A']})
        offset += len(frames)
        for source in (path, folder / 'display_frames.json', folder / 'acceptance.json'):
            hashes[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
    commands += ['trajout asp272_relay_progress.dcd dcd namdcell', 'run', 'quit']
    cpp = OUT / 'concatenate.cpptraj'
    cpp.write_text('\n'.join(commands) + '\n')
    with (OUT / 'cpptraj.log').open('w') as handle:
        subprocess.run(['/data/jpc/code/ambertools26/venv/bin/cpptraj', '-i', str(cpp)],
                       cwd=OUT, stdout=handle, stderr=subprocess.STDOUT, check=True)
    trajectory = OUT / 'asp272_relay_progress.dcd'
    if header(trajectory) != (offset, 19976):
        raise ValueError('Output atom/frame mismatch')
    if trajectory.stat().st_size > 100_000_000:
        raise ValueError('Movie exceeds laptop transfer limit')
    loader = '''# Actual optimizer evaluations, including trial steps; not physical time.
set ::nyl12_dir [file dirname [file normalize [info script]]]
foreach name {qm_indices.tcl nyl12_es_vmd_common.tcl} {
    source [file join $::nyl12_dir $name]
}
nyl12_load stripped.parm7 asp272_relay_progress.dcd FRAMES $::nyl12_qm_indices $::nyl12_stripped_atom_count
puts "Asp272 relay: completed windows only; optimization evaluations, not MD."
'''.replace('FRAMES', str(offset))
    (OUT / 'asp272_relay_progress.vmd').write_text(loader)
    manifest = {'frames': offset, 'atoms': 19976, 'snapshot_status': status,
                'evaluation_stride': stride, 'windows': records,
                'source_sha256': hashes, 'interpolation_or_alignment': False,
                'note': 'Includes optimizer trial evaluations; not physical time or a validated path.'}
    for source in (Path(__file__), HEADER_SOURCE, VIEW / 'stripped.parm7',
                   VIEW / 'qm_indices.tcl', VIEW / 'nyl12_es_vmd_common.tcl'):
        hashes[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
    (OUT / 'asp272_relay_frames.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (OUT / 'asp272_relay_README.md').write_text(
        f'# Asp272 relay progress\n\n{len(labels)} completed windows, {offset} frames. '
        f'Saved evaluation stride {stride}, with every exact endpoint included. '
        'All hydrogens, protein, substrate, and three fixed waters retained; '
        'bulk water and ions removed. QM in Licorice; non-QM ligand carbon orange; '
        'protein cartoon available but off. Trial evaluations may move backward. '
        'No physical time, interpolation, or alignment.\n\n'
        'In VMD Tk Console: `source asp272_relay_progress.vmd`\n')
    for name in ('asp272_relay_progress.dcd', 'asp272_relay_progress.vmd',
                 'asp272_relay_frames.json', 'asp272_relay_README.md'):
        shutil.copy2(OUT / name, VIEW / name)
    # Extend the existing checksum manifest without dropping older entries.
    old = (VIEW / 'MANIFEST.sha256').read_text()
    entries = []
    for name in ('asp272_relay_progress.dcd', 'asp272_relay_progress.vmd',
                 'asp272_relay_frames.json', 'asp272_relay_README.md'):
        entries.append(hashlib.sha256((VIEW / name).read_bytes()).hexdigest() + '  ' + name)
    (VIEW / 'MANIFEST.sha256').write_text(old.rstrip() + '\n' + '\n'.join(entries) + '\n')
    logging.info('Verified %d windows, %d frames, %d atoms, %.1f MB',
                 len(labels), offset, 19976, trajectory.stat().st_size / 1e6)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
