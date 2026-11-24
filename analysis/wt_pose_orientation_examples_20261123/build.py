"""Extract and align two WT active sites; preserve source and measurement provenance."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from nyl12_covalent.catalytic_alignment_geometry import fit_coordinates
from nyl12_covalent.active_site_geometry import projected_herschlag_phi

ROOT = Path('/data/jpc/covalent-enzyme-design')
OUT = Path(__file__).resolve().parent
CSV = ROOT / 'analysis/native_thr_gate_calibration_20261123/inputs/occupied_sites.csv'
BASE = ROOT / 'design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation/models/WT'
CASES = [('flipped', '04', 'AEI', '2026112204', 'A'), ('positive', '02', 'BFJ', '2026112202', 'B')]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    lines = [line for line in path.read_text().splitlines() if line[:6].strip() in {'ATOM', 'HETATM'}]
    result = {}
    for line in lines:
        key = (line[21], int(line[22:26]), line[12:16].strip())
        if key in result:
            raise ValueError(f'duplicate atom {key}')
        result[key] = (line, np.array([float(line[i:i+8]) for i in (30, 38, 46)]))
    return result

def point(atoms, chain, residue, name):
    return atoms[(chain, residue, name)][1]

rows = list(csv.DictReader(CSV.open()))
data = []
for label, seed, chains, number, tier in CASES:
    path = BASE / f'seed_{seed}/relaxation/relaxed.pdb'
    atoms = read(path)
    row = next(r for r in rows if r['variant'] == 'WT' and r['seed'] == number and r['stage'] == 'relaxed' and r['site'] == chains[:2])
    if row['ligand'] != f'{chains[2]}:1' or row['tier'] != tier:
        raise ValueError('Unexpected source mapping or tier')
    n = point(atoms, chains[1], 1, 'N')
    ligand = lambda atom: point(atoms, chains[2], 1, atom)
    distance = float(np.linalg.norm(n - ligand('N59')))
    phi = projected_herschlag_phi(point(atoms, chains[1], 1, 'OG1'), ligand('O49'), ligand('C54'), ligand('N59'), ligand('C62'))
    if abs(distance - float(row['nterm_to_leaving_n'])) > 1e-8 or abs(phi - float(row['herschlag_phi'])) > 1e-8:
        raise ValueError('Independent source geometry differs from CSV')
    fit = np.array([point(atoms, chains[0] if residue < 233 else chains[1], residue if residue < 233 else residue - 232, atom) for residue in (113, 186, 272, 274) for atom in ('N', 'CA', 'C', 'O')])
    data.append((label, seed, chains, path, atoms, row, distance, phi, fit))
rotation, center, target, rmsd = fit_coordinates(data[1][-1], data[0][-1])
manifest = {'alignment': {'method': 'Kabsch; positive onto flipped; no Thr233 or ligand atoms in fit', 'residues_canonical': [113, 186, 272, 274], 'atoms_per_residue': ['N', 'CA', 'C', 'O'], 'matched_atoms': 16, 'rmsd_angstrom': rmsd}, 'source_csv': str(CSV), 'source_csv_sha256': sha(CSV), 'examples': []}
for label, seed, chains, path, atoms, row, distance, phi, fit in data:
    selected = []
    for (chain, residue, name), (line, xyz) in atoms.items():
        if ((chain == chains[0] and residue in (113, 186)) or (chain == chains[1] and residue in (1, 40, 42)) or (chain == chains[2] and residue == 1)):
            if label == 'positive': xyz = (xyz - center) @ rotation + target
            selected.append(line[:30] + ''.join(f'{value:8.3f}' for value in xyz) + line[54:])
    out = OUT / f'{label}.pdb'
    out.write_text('\n'.join(selected) + '\nEND\n')
    manifest['examples'].append({'label': label, 'seed': row['seed'], 'site': row['site'], 'source_pdb': str(path), 'source_pdb_sha256': sha(path), 'extracted_pdb_sha256': sha(out), 'nalpha_atom': f'{chains[1]}:1:N', 'leaving_atom': f'{chains[2]}:1:N59', 'nalpha_to_leaving_n_angstrom': distance, 'projected_phi_degrees': phi, 'tier': tier, 'native_pass': row['native_pass'], 'display_coordinates': 'source' if label == 'flipped' else 'rigid Kabsch transform, PDB coordinates rounded to 0.001 A'})
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
