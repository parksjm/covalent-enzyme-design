"""Summarize verified ES seeds without mechanism assignment or re-selection."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

BASE = Path(__file__).resolve().parent
OUTPUT = BASE / 'seed_inspection_attempt4'
SOURCE = OUTPUT / 'inspection.json'
report = json.loads(SOURCE.read_text())
rows = []
for frame in report['frames']:
    metrics = frame['geometry']['metrics']
    waters = frame['water_inventory']
    atoms = waters['selected_atoms']
    nitrogen = atoms['thr_n_terminus']['coordinates_A']
    row = {'seed': frame['name'], 'physical_time_ns': frame['source']['physical_time_ns'],
           'd1_A': metrics['d1_thr929_og1_scissile_c'],
           'd2_A': metrics['d2_phe113_n_scissile_o'],
           'd3_A': metrics['d3_asn186_nd2_scissile_o'],
           'aBD_degrees': metrics['a_BD_og1_c_o'],
           'raw_NC_phi_degrees': metrics['herschlag_phi_projected'],
           'historical_PDB_phi_degrees': frame['geometry']['historical_PDB']['phi_degrees'],
           'local_water_count': waters['local_water_count'],
           'geometric_only_candidate_count': waters['geometric_only_candidate_count']}
    for residue in [272, 274]:
        row[f'nearest_Asp{residue}_O_ThrN_A'] = min(
            math.dist(nitrogen, atoms[f'asp{residue}_{oxygen}']['coordinates_A'])
            for oxygen in ['od1', 'od2'])
    rows.append(row)
with (OUTPUT / 'seed_summary.csv').open('x', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
notes = ['ES seeds; not optimized TI/TS. Water counts describe these selected snapshots only.',
         'Raw NC phi remains primary; historical PDB phi is independently source-precision verified.',
         'Asp distances are descriptive; residue correspondence does not assign catalysis.', '']
notes += ['| ' + ' | '.join(rows[0]) + ' |', '| ' + ' | '.join(['---'] * len(rows[0])) + ' |']
for row in rows:
    notes.append('| ' + ' | '.join(f'{value:.4f}' if isinstance(value, float) else str(value) for value in row.values()) + ' |')
(OUTPUT / 'seed_summary.md').write_text('\n'.join(notes) + '\n')
provenance = [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
              for path in [Path(__file__), SOURCE, OUTPUT / 'seed_summary.csv', OUTPUT / 'seed_summary.md']]
(OUTPUT / 'summary_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
