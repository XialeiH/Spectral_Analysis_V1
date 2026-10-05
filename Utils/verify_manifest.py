"""Verify the published code and input files using only the Python standard library."""
import hashlib
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    records = json.loads((root / 'MANIFEST.json').read_text())
    errors = []
    for record in records:
        path = root / record['path']
        if not path.is_file():
            errors.append(f"Missing: {record['path']}")
            continue
        digest = hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(block)
        if digest.hexdigest() != record['sha256']:
            errors.append(f"Changed: {record['path']}")
    if errors:
        raise SystemExit('\n'.join(errors))
    print(f'PASS: {len(records)} published files match their SHA-256 checksums.')


if __name__ == '__main__':
    main()
