"""Independently audit positional and driving restraint energies."""
import json
import re
from pathlib import Path
import numpy as np
from parmed.amber import AmberParm
ROOT=Path('/data/jpc/covalent-enzyme-design')
PILOT=ROOT/'output_nyl12_pa6_dftb_pilot_20261024'
RUN=PILOT/'coordinate_driving_forward_reverse_attempt3'
inputs=json.loads((RUN/'inputs.json').read_text())
ref=AmberParm(inputs['topology'],xyz=inputs['reference'])
variant=json.loads((ROOT/'configs/pa6_es_dftb_pilot_20261024.json').read_text())['variants'][0]
qm={atom['index']-1 for group in variant['groups'] for atom in group['atoms']}
mm=[i for i in range(len(ref.atoms)) if i not in qm]
rows=[]
for folder in sorted(RUN.glob('forward_*')):
    final=AmberParm(inputs['topology'],xyz=str(folder/'optimized.ncrst'))
    process=json.loads((folder/'process.json').read_text())
    acceptance=json.loads((folder/'acceptance.json').read_text())
    d=acceptance['distances_A']
    drive=25*((d['attack']-process['target_attack_A'])**2+(d['OG_H']-d['N_H']-process['target_proton_A'])**2)
    position=10*float(np.sum((final.coordinates[mm]-ref.coordinates[mm])**2))
    text=(folder/'minimized.mdout').read_text()
    printed=float(re.findall(r'RESTRAINT\s*=\s*([-+\d.Ee]+)',text)[-1])
    row={'window':folder.name,'driving_kcal_mol':drive,'position_kcal_mol':position,'printed_RESTRAINT':printed,
         'printed_minus_computed_kcal_mol':printed-position-drive,'RMS':acceptance['RMS']}
    if abs(row['printed_minus_computed_kcal_mol'])>0.0002:
        raise ValueError(row)
    rows.append(row)
(PILOT/'overnight_review_20261025/restraint_audit.json').write_text(json.dumps({'all18_verified':len(rows)==18,'rows':rows},indent=2)+'\n')
