"""Record the external Nyl12 artifacts retained during the 2026-10-20 import."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path('/data/jpc/enzyme-substrate-modeling')
OLD = ROOT / 'output_nyl12_covalent_qmmm_20261019'
NEW = ROOT / 'output_nyl12_ti1_refinement_20261020'
paths: list[Path] = []
for substrate in ('pa6_l2', 'poly_xama_maxa'):
    model_name = f'nyl12_{substrate}_ti1_af3_seed20261019'
    paths.append(OLD / f'{substrate}_af3_ti1/prediction' / model_name / f'{model_name}_model.cif')
    paths.append(OLD / f'{substrate}_af3_ti1/canonical/audit.json')
    paths.extend(OLD / f'{substrate}_ti1_af3_neutral_beta' / name for name in (
        'ti1.prmtop', 'ti1.inpcrd', 'ti1_audit.json', 'ti1_provenance.json'
    ))
    paths.extend(NEW / substrate / 'water_search' / name for name in (
        'water_search.json', 'contact_diagnostic_v2.json'
    ))
paths.extend(NEW / 'jobs' / job / 'report.json' for job in (
    'qm_region_audit_pa6_r1', 'qm_region_audit_poly_r1'
))
paths.extend(NEW / 'gaussian_preflight_violet_20261020' / name for name in (
    'preflight.json', 'water.gjf', 'water.log'
))
paths.extend(NEW / 'jobs' / job / 'status.json' for job in (
    'g16_preflight', 'repository_export', 'repository_import'
))
records = []
for path in paths:
    data = path.read_bytes()
    records.append({'path': str(path), 'size_bytes': len(data),
                    'sha256': hashlib.sha256(data).hexdigest()})
output = Path(__file__).with_name('external_artifacts_20261020.json')
with output.open('x') as stream:
    json.dump({'purpose': 'external artifact locations and hashes; not copied into Git',
               'files': records}, stream, indent=2)
    stream.write('\n')
