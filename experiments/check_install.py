#!/usr/bin/env python3
"""Checks that this folder is an unmodified copy of the package: compares every file with MANIFEST.txt and lists missing, changed and extra files
(extra .java files are the usual cause of compile errors when an older repository was copied into the same folder).
    python experiments/check_install.py"""
import hashlib, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
SKIP = ('experiments/work', 'experiments/work_quick', '__pycache__', 'MANIFEST.txt', 'RESULTS/')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
manifest = {}
for line in (ROOT / 'MANIFEST.txt').read_text().splitlines():
    h, rel = line.split('  ', 1)
    manifest[rel] = h
present = {p.relative_to(ROOT).as_posix(): p for p in ROOT.rglob('*') if p.is_file() and not any(s in p.relative_to(ROOT).as_posix() for s in SKIP)}
missing = sorted(set(manifest) - set(present)); extra = sorted(set(present) - set(manifest))
changed = sorted(r for r in manifest if r in present and sha(present[r]) != manifest[r])
print(f'{len(manifest)} files expected, {len(present)} found')
for name, lst in (('MISSING', missing), ('CHANGED', changed), ('EXTRA', extra)):
    print(f'{name}: {len(lst)}')
    for r in lst[:40]: print('   ', r)
print('OK: identical to the package' if not (missing or changed or extra) else 'NOT identical: re-extract the zip into an empty folder')
sys.exit(0 if not (missing or changed or extra) else 1)
