"""Restore byte-identical MAT files split only to satisfy GitHub file limits."""
import hashlib
import json
from pathlib import Path

root=Path(__file__).resolve().parents[1]
for entry in json.loads((root/'Data/chunks.json').read_text()):
    target=root/entry['path']
    if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()==entry['sha256']:
        continue
    temporary=target.with_suffix(target.suffix+'.restoring')
    digest=hashlib.sha256()
    with temporary.open('wb') as output:
        for name in entry['parts']:
            data=(root/name).read_bytes()
            output.write(data)
            digest.update(data)
    if digest.hexdigest()!=entry['sha256']:
        temporary.unlink()
        raise RuntimeError('Checksum mismatch: '+entry['path'])
    temporary.replace(target)
    print('Restored',entry['path'])
