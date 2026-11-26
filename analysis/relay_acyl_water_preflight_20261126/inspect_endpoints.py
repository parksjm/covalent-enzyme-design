"""Inspect fixed, source-verified acyl-side structures without changing them."""
from __future__ import annotations
import json
import logging
from pathlib import Path
import numpy as np
from parmed.amber import Rst7
from nyl12_covalent.deacylation_water_geometry import _topology, angle, sha256
from nyl12_covalent.design_md_geometry import image_near
from nyl12_covalent.relay_umbrella_cv import raw_features


def main() -> None:
    """Measure all waters in one endpoint and six prespecified US checkpoints."""
    root = Path('/data/jpc/covalent-enzyme-design')
    base = root/'qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K'
    prepared_path = base/'prepared.json'
    prepared = json.loads(prepared_path.read_text())
    topology, waters = _topology(base/'inputs/system.prmtop')
    seeds = [(base/'inputs/seeds/02_acylation_forward_endpoint.ncrst', None)]
    for i in range(21, 27):
        stage = base/f'delivery_full_{i:03d}/07_production_5ps'
        receipt = json.loads((stage/'acceptance.json').read_text())
        path = stage/'restart.ncrst'
        if not receipt['accepted'] or sha256(path) != receipt['checkpoint_sha256']:
            raise ValueError('Unaccepted or changed endpoint: '+str(path))
        seeds.append((path, stage/'acceptance.json'))
    records = []
    water_ids = [w['o'] for w in waters]
    for path, receipt in seeds:
        if receipt is None:
            expected = next(w['seed_source_sha256'] for w in prepared['windows'] if w['label']=='delivery_full_026')
            if sha256(path) != expected:
                raise ValueError('Acyl seed hash mismatch')
        restart = Rst7.open(str(path))
        xyz = np.asarray(restart.coordinates).reshape(-1,3)
        carbon = xyz[19915]
        oxy = image_near(xyz[[19916]], carbon, restart.box)[0]
        position = image_near(xyz[water_ids], carbon, restart.box)
        distances = np.linalg.norm(position-carbon, axis=1)
        angles = angle(position,carbon,np.broadcast_to(oxy,position.shape))
        near = np.flatnonzero((distances<=3.5)&(angles>=95)&(angles<=125))
        j = int(np.argmin(distances))
        records.append(dict(source=str(path),source_sha256=sha256(path),acceptance=str(receipt) if receipt else None,
            C_OG1_A=float(np.linalg.norm(xyz[19915]-xyz[13915])),
            C_Nleave_A=float(np.linalg.norm(xyz[19915]-xyz[19928])),
            features=raw_features(xyz,restart.box),
            nearest_water_residue=waters[j]['residue'],nearest_water_OC_A=float(distances[j]),nearest_angle_deg=float(angles[j]),
            candidate_waters=[dict(residue=waters[k]['residue'],OC_A=float(distances[k]),angle_deg=float(angles[k])) for k in near]))
    output = Path(__file__).with_name('report.json')
    if output.exists():
        raise FileExistsError(output)
    output.write_text(json.dumps(dict(scope='Prespecified relay acyl endpoint and delivery021..026 accepted25ps checkpoints; no simulations or ensemble claim',water_count=len(waters),prepared_sha256=sha256(prepared_path),script_sha256=sha256(Path(__file__)),records=records),indent=2)+'\n')
    for item in records:
        logging.info('%s: C-OG1 %.3f C-N %.3f closestwater %.3fA angle %.1fdeg candidates%d',Path(item['source']).parent.name,item['C_OG1_A'],item['C_Nleave_A'],item['nearest_water_OC_A'],item['nearest_angle_deg'],len(item['candidate_waters']))


if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)
    main()
