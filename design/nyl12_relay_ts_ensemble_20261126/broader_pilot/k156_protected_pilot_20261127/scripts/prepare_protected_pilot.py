from pathlib import Path
import json,shutil,hashlib,sys
base=Path('design/nyl12_relay_ts_ensemble_20261126/broader_pilot').resolve();old=base/'protection_corrected_pilot_20261127';new=base/'k156_protected_pilot_20261127';new.mkdir();snapshot=new/'source_snapshot';snapshot.mkdir();shutil.copytree(Path('src').resolve(),snapshot/'src',ignore=shutil.ignore_patterns('__pycache__','*.pyc'));(snapshot/'scripts').mkdir()
for name in ['run_relay_broad_pilot.py','score_relay_ensemble_sequence.py']:shutil.copy2(Path('scripts')/name,snapshot/'scripts'/name)
m=json.loads((old/'effective_manifest.json').read_text());proof=Path('docs/figures/nyl12_sequence_protection_20261127/canonical_protection_restored_K156_20261127.json').resolve();p=json.loads(proof.read_text());mask=p['mutable'];assert p['counts']['fixed_total']==101 and p['counts']['mutable_total']==227 and 156 not in mask['alpha']
sys.path.insert(0,str(snapshot/'src'))
from nyl12_covalent.relay_broad_pilot import symmetry_manifest
flags=symmetry_manifest(mask)
for key,directory in [('generation_engine','native_generation'),('score_engine','native_scoring')]:
 engine=Path(m['native'][key]);dest=snapshot/directory/'ligandmpnn';shutil.copytree(engine.parent,dest,ignore=shutil.ignore_patterns('__pycache__','*.pyc'));m['native'][key]=str(dest/engine.name)
shutil.copytree(Path(m['native']['support']),snapshot/'runtime_support',ignore=shutil.ignore_patterns('__pycache__','*.pyc'));m['native']['support']=str(snapshot/'runtime_support')
m['native']['runner']=str(snapshot/'scripts/run_relay_broad_pilot.py');m['native']['score_wrapper']=str(snapshot/'scripts/score_relay_ensemble_sequence.py');m['native']['environment']['PYTHONPATH']=str(snapshot/'runtime_support');m['native']['environment']['PYTHONDONTWRITEBYTECODE']='1';m['native']['environment']['NUMEXPR_NUM_THREADS']='2'
for key in ['approved_for_launch','root_mask_decision','parser_audit']:m['native'].pop(key,None)
m['native'].update(mask=mask,**flags)
r=m['request'];r.update(campaign='nyl12_relay_k156_protected_pilot_20261127',status='prepared',protection_receipt=str(proof),protection_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),known_blockers=[]);r['generation'].update(mask=mask,**flags);r['execution_gate']['Negoro_audit_source']=str(Path('analysis/negoro_mutagenesis_mapping_20261127/nylc_critical_residue_evidence.json').resolve());r['execution_gate']['K156_fixed']=True
r['downstream']='Generation+125matchedscoresonly; no newRFD/Boltz/Amber authorized bythispilotmanifest'
(new/'manifest_draft.json').write_text(json.dumps(m,indent=2)+'\n')
