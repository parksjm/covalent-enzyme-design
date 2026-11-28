"""Render source-bound molecular figures using the existing PyMOL environment."""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pymol2

ROOT = Path('/data/jpc/covalent-enzyme-design')
HERE = Path(__file__).resolve().parent
BASE = ROOT / 'design/nyl12_relay_ts_ensemble_20261126'
WT = BASE / 'broader_pilot/wt_protein_template_controls_20261127/amber_relaxation_20261127'


def atoms(path: Path) -> dict[tuple[str, str, str], np.ndarray]:
    """Return native atom coordinates keyed by canonical identity."""
    return {(s[21], s[22:26].strip(), s[12:16].strip()): np.array(
        [float(s[30:38]), float(s[38:46]), float(s[46:54])])
        for s in path.read_text().splitlines() if s.startswith(('ATOM  ', 'HETATM'))}


def frame(data: dict, alpha: str, beta: str) -> tuple[np.ndarray, np.ndarray]:
    """Construct a camera with Phe left, Asn right and Thr below."""
    phe, asn, thr = [data[key] for key in [(alpha, '113', 'CA'),
                      (alpha, '186', 'CA'), (beta, '233', 'CA')]]
    x = asn - phe
    x /= np.linalg.norm(x)
    y = (phe + asn) / 2 - thr
    y -= x * np.dot(y, x)
    y /= np.linalg.norm(y)
    z = np.cross(x, y)
    return np.stack([x, y, z], axis=1), (phe + asn + thr) / 3


def main() -> None:
    """Produce molecular PNGs, sessions and a source/camera provenance record."""
    logging.basicConfig(level=logging.INFO)
    sources = {}
    for model in ['af3', 'boltz']:
        for state in ['before', 'after']:
            case = f'{model}_template_seed20261126'
            sources[f'{model}_{state}'] = WT / case / 'postprocess' / f'{case}_{state}' / 'canonical_precursor.pdb'
    sources['ts_like'] = BASE / 'templates/medoid_01_delivery_full_010_frame_068/complex_heavy.pdb'
    ref = atoms(sources['boltz_after'])
    basis, center = frame(ref, 'B', 'F')
    provenance = {'sources': {}, 'figures': {}, 'native_files_modified': False,
                  'hydrogens_added': False, 'threads': 2,
                  'connectivity': 'PyMOL PDB connectivity; display valence=0, no force-field changes',
                  'scope': 'WT local panels display full heavy-only L4 and Phe113, Asn186, Thr233; no hidden ligand atoms.'}
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        for key, path in sources.items():
            cmd.reinitialize()
            cmd.bg_color('white')
            for setting, value in {'max_threads': 2, 'valence': 0, 'orthoscopic': 1,
                'ray_opaque_background': 1, 'antialias': 2, 'depth_cue': 0,
                'ray_shadows': 0, 'ray_trace_mode': 1, 'stick_radius': .18,
                'sphere_scale': .23, 'specular': .15}.items():
                cmd.set(setting, value)
            cmd.load(str(path), 'mol')
            raw = atoms(path)
            if key == 'ts_like':
                rotation, origin = frame(raw, 'A', 'E')
                transform = lambda xyz: (xyz - origin) @ rotation
                alpha, beta, ligand = 'A', 'E', 'I'
                alignment = {'method': 'own catalytic reference frame; no coordinate optimization'}
            else:
                shared = sorted(k for k in ref.keys() & raw.keys() if k[2] == 'CA'
                    and k[0] in 'BF' and int(k[1]) not in {113, 186, 233}
                    and not 95 <= int(k[1]) <= 103 and not 289 <= int(k[1]) <= 300)
                p, q = np.array([raw[k] for k in shared]), np.array([ref[k] for k in shared])
                pc, qc = p.mean(axis=0), q.mean(axis=0)
                u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
                correction = np.eye(3)
                correction[2, 2] = np.linalg.det(u @ vt)
                fit = u @ correction @ vt
                transform = lambda xyz: ((xyz - pc) @ fit + qc - center) @ basis
                alpha, beta, ligand = 'B', 'F', 'J'
                alignment = {'method': 'Kabsch on canonical B/F CA excluding catalytic residues and highlighted loops',
                    'atoms': len(shared), 'rmsd_A': float(np.sqrt(np.mean(np.sum(((p-pc) @ fit + qc-q)**2, axis=1))))}
            native = cmd.get_coords('mol')
            cmd.load_coords(transform(native), 'mol')
            local = f'(chain {alpha} and resi 113+186) or (chain {beta} and resi 233) or chain {ligand}'
            cmd.hide('everything')
            cmd.show('sticks', 'mol and (' + local + ')')
            cmd.set_color('protein_c', [0.0, .72, .80])
            cmd.set_color('substrate_c', [1.0, .48, .06])
            cmd.set_color('light_h', [.72, .72, .72])
            cmd.color('protein_c', 'mol and elem C')
            cmd.color('substrate_c', f'mol and chain {ligand} and elem C')
            for element, color in [('N','blue'), ('O','red'), ('H','light_h'), ('S','yellow')]:
                cmd.color(color, 'mol and elem ' + element)
            if key != 'ts_like':
                for name, first, second in [
                    ('attack', f'chain {beta} and resi 233 and name OG1', f'chain {ligand} and name C55'),
                    ('oxyanion', f'chain {alpha} and resi 186 and name ND2', f'chain {ligand} and name O50')]:
                    cmd.distance(name, 'mol and '+first, 'mol and '+second)
                    cmd.hide('labels', name)
                    cmd.color('gray50', name)
                    cmd.set('dash_radius', .035, name)
            cmd.reset()
            # One shared reference camera and scale for all four WT panels.
            if key != 'ts_like':
                anchors = transform(np.array([raw[k] for k in [(alpha,'113','CA'),(alpha,'186','CA'),(beta,'233','CA')]]))
                if key == 'af3_before':
                    # Use fixed world-space center/extent, independent of ligand pose.
                    cmd.pseudoatom('camera_center', pos=[0,0,0])
                    cmd.zoom('camera_center', buffer=11.5)
                    camera = cmd.get_view()
                    cmd.delete('camera_center')
                cmd.set_view(camera)
            else:
                cmd.zoom('visible', buffer=2.3)
            cmd.png(str(HERE / f'{key}.png'), 1600, 1400, dpi=300, ray=1, quiet=1)
            cmd.save(str(HERE / f'{key}.pse'))
            provenance['sources'][key] = {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
            provenance['figures'][key] = {'canonical_site': alpha+beta, 'ligand_chain': ligand,
                'camera': list(cmd.get_view()), 'alignment': alignment,
                'display_coordinates': {f'{ch}:{resi}:{atom}': transform(raw[(ch,resi,atom)]).tolist()
                    for ch,resi,atom in [(alpha,'113','CA'),(alpha,'186','CA'),(beta,'233','CA')]}}
            logging.info('Rendered %s', key)
            if key == 'ts_like':
                cmd.hide('everything')
                cmd.show('cartoon', 'mol and chain A+B+C+D+E+F+G+H')
                for chains, color in [('A+E','green'),('B+F','cyan'),('C+G','magenta'),('D+H','yellow')]:
                    cmd.color(color, 'mol and chain '+chains)
                cmd.set_color('loop_alpha', [.15,.25,.85])
                cmd.set_color('loop_beta', [.9,.25,.1])
                cmd.color('loop_alpha','mol and chain A and resi 95-103')
                cmd.color('loop_beta','mol and chain E and resi 289-300')
                cmd.set('cartoon_loop_radius', .45, 'mol and chain A and resi 95-103')
                cmd.set('cartoon_loop_radius', .45, 'mol and chain E and resi 289-300')
                cmd.show('sticks', 'mol and (chain I or (chain E and resi 233))')
                cmd.color('substrate_c','mol and chain I and elem C')
                cmd.color('red','mol and (chain I or (chain E and resi 233)) and elem O')
                cmd.color('blue','mol and (chain I or (chain E and resi 233)) and elem N')
                cmd.zoom('mol', buffer=3)
                cmd.png(str(HERE/'assembly_loops.png'),1800,1600,dpi=300,ray=1,quiet=1)
                cmd.save(str(HERE/'assembly_loops.pse'))
                provenance['figures']['assembly_loops']={'source':'ts_like','camera':list(cmd.get_view()),
                    'loops':{'alpha_A_95_103':'deep blue','beta_E_289_300':'orange-red'},
                    'canonical_chains':'A/E green, B/F cyan, C/G magenta, D/H yellow; ligand I orange C',
                    'note':'Assembly context from TS-like template, not reference crystal; label ligand as modeled L4/PA6 orientation.'}
    (HERE/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
