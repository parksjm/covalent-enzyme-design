"""Audit saved endpoint geometry without modifying coordinates."""
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import parmed as pmd
from parmed import unit as u
from parmed.geometry import box_lengths_and_angles_to_vectors

ROOT = Path('/data/jpc/covalent-enzyme-design')
OUT = Path(__file__).parent
BASE = ROOT / 'output_nyl12_pa6_dftb_pilot_20261024'
CONFIG = ROOT / 'configs/pa6_es_dftb_pilot_20261024.json'
TOP = Path(json.loads(CONFIG.read_text())['topology'])
STATES = {
    'baseline_ES': BASE / 'analytic_relaxation_continuation_attempt1/optimized.ncrst',
    'failed_forward': BASE / 'coordinate_driving_forward_reverse_attempt3/forward_17/optimized.ncrst',
    'strict_ES': BASE / 'es_tight_convergence_attempt1/optimized.ncrst',
}
LABELS = {'ThrN':13905,'ThrH1':13906,'ThrH2':13907,'ThrOG':13916,'ThrHG':13917,
          'Asp272OD1':14559,'Asp272OD2':14560,'Asp274OD1':14578,'Asp274OD2':14579}
structure = pmd.load_file(str(TOP))
waters = [a for a in structure.atoms if a.residue.name in ('WAT','HOH') and a.atomic_number == 8]
results = {}
for label, path in STATES.items():
    restart = pmd.amber.Rst7.open(str(path))
    xyz = np.asarray(restart.coordinates).reshape(-1,3)
    box = restart.box
    if box is None:
        raise ValueError('Full solvated restart has no periodic box')
    vectors = np.array(box_lengths_and_angles_to_vectors(*box).value_in_unit(u.angstrom), dtype=float)
    shifts = np.array(list(itertools.product((-1,0,1), repeat=3))) @ vectors
    inverse = np.linalg.inv(vectors)

    def delta(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """Return nearest lattice-image displacement using 27 candidates."""
        raw = b-a
        central = raw - np.rint(raw @ inverse) @ vectors
        candidates = central + shifts
        return candidates[np.argmin(np.linalg.norm(candidates,axis=1))]

    def point(index: int) -> np.ndarray:
        return xyz[index-1]

    def distance(a: int, b: int) -> float:
        return float(np.linalg.norm(delta(point(a),point(b))))

    def angle(d: int, h: int, a: int) -> float:
        dh = delta(point(h),point(d))
        ah = delta(point(h),point(a))
        return float(np.degrees(np.arccos(np.clip(np.dot(dh,ah)/np.linalg.norm(dh)/np.linalg.norm(ah),-1,1))))

    contacts = []
    for acceptor, idx in LABELS.items():
        if 'Asp' not in acceptor:
            continue
        for hydrogen in ('ThrH1','ThrH2'):
            contacts.append({'acceptor':acceptor,'donor':'ThrN','hydrogen':hydrogen,
                             'D_A':distance(13905,idx),'H_A':distance(LABELS[hydrogen],idx),
                             'D_H_A_angle':angle(13905,LABELS[hydrogen],idx)})
    for acceptor in ('ThrN','Asp272OD1','Asp272OD2','Asp274OD1','Asp274OD2'):
        idx=LABELS[acceptor]
        contacts.append({'acceptor':acceptor,'donor':'ThrOG','hydrogen':'ThrHG',
                         'D_A':distance(13916,idx),'H_A':distance(13917,idx),
                         'D_H_A_angle':angle(13916,13917,idx)})
    nearby = {}
    for target in ('ThrN','ThrOG','Asp272OD1','Asp272OD2','Asp274OD1','Asp274OD2'):
        ranked=sorted(waters,key=lambda a:distance(LABELS[target],a.idx+1))[:3]
        nearby[target]=[{'amber_residue':a.residue.idx+1,'oxygen_atom':a.idx+1,
                         'distance':distance(LABELS[target],a.idx+1)} for a in ranked]
    results[label]={'box':list(map(float,box)),'contacts':contacts,'nearest_water_oxygens':nearby,
                    'covalent_distances':{'OG_HG':distance(13916,13917),'N_H1':distance(13905,13906),'N_H2':distance(13905,13907)}}
inputs=[CONFIG,TOP,*STATES.values(),Path(__file__)]
result={'units':'Angstrom and degrees; donor-H-acceptor angle; Amber one-based atom/residue IDs',
        'method':'Nearest periodic image via fractional centering and 27 triclinic images; waters selected by topology WAT/HOH oxygen',
        'hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'states':results,
        'scientific_acceptance':'Geometry diagnostic only; no endpoint convergence or mechanism inferred'}
(OUT/'geometry.json').write_text(json.dumps(result,indent=2)+'\n')
