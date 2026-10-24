"""Freeze the reference-origin clarification and refresh the small-file export."""
from __future__ import annotations

import fcntl
import hashlib
import json
import logging
import shutil
import sys
from pathlib import Path

ROOT = Path('/data/jpc/covalent-enzyme-design')
sys.path.insert(0, str(ROOT))
from scripts.prepare_nyl12_inspection_download import build_package


def digest(path: Path) -> str:
    """Return a retained small artifact's SHA256."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    """Preserve the original checkpoint and verify the amended package."""
    archive = ROOT / 'provenance/pa6_reference_origin_update_20261024'
    current = ROOT / 'output_laptop_inspection_20261021'
    staging = ROOT / 'output_laptop_inspection_20261021_reference_origin_staging'
    backup = ROOT / 'output_laptop_inspection_20261021_revision_31'
    for path in (archive, staging, backup):
        if path.exists():
            raise FileExistsError(path)
    names = (
        'docs/pa6_es_reuse_audit_20261024.md',
        'docs/standalone_handoff_20261020.md',
        'manuscript/METHODS_DRAFT.md',
        'manuscript/SOURCE_LEDGER.md',
    )
    prior = json.loads((current / 'provenance.json').read_text())
    files = {name: Path(row['source']) for name, row in prior['files'].items()}
    before = {name: digest(current / name) for name in files}
    archive.mkdir()
    for name in names:
        target = archive / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, target)
        if digest(ROOT / name) != digest(target):
            raise ValueError(name)
        files[name] = target
    shutil.copy2(Path(__file__), archive / Path(__file__).name)
    journal = archive / 'agent_actions_at_checkpoint.jsonl'
    with (ROOT / 'output_agent_actions/20261021.jsonl').open('rb') as source:
        fcntl.flock(source, fcntl.LOCK_SH)
        try:
            with journal.open('xb') as target:
                shutil.copyfileobj(source, target)
        finally:
            fcntl.flock(source, fcntl.LOCK_UN)
    journal_name = 'provenance/agent_actions_at_export.jsonl'
    files[journal_name] = journal
    report = build_package(files, staging)
    for name, expected in before.items():
        if digest(current / name) != expected:
            raise ValueError(f'Current package changed: {name}')
        if name not in {*names, journal_name} and digest(staging / name) != expected:
            raise ValueError(f'Unrelated artifact changed: {name}')
    for line in (staging / 'MANIFEST.sha256').read_text().splitlines():
        expected, name = line.split('  ', 1)
        if digest(staging / name) != expected:
            raise ValueError(f'Package manifest mismatch: {name}')
    size = sum(path.stat().st_size for path in staging.rglob('*') if path.is_file())
    if size > 25_000_000:
        raise ValueError('Export exceeds 25 MB')
    current.rename(backup)
    staging.rename(current)
    verification = {
        'verified': True, 'payload_files': len(report['files']),
        'complete_package_bytes': size, 'retained_backup': str(backup),
        'original_checkpoint_preserved': 'provenance/pa6_es_reuse_checkpoint_20261024',
        'prior_payload_files_checked': len(before),
        'clarification': 'User is coauthor; source manuscript in revision; citation pending publication.',
    }
    (archive / 'package_verification.json').write_text(json.dumps(verification, indent=2) + '\n')
    manifest = {str(path.relative_to(archive)): digest(path)
                for path in sorted(archive.rglob('*')) if path.is_file()}
    (archive / 'source_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    for name, expected in manifest.items():
        if digest(archive / name) != expected:
            raise ValueError(f'Archive manifest mismatch: {name}')
    logging.info('Verified reference-origin amendment and %s-byte package', size)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    main()
