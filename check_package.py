#!/usr/bin/env python3
"""Standard-library verification of packaged file and canonical-data hashes."""
from __future__ import annotations
import gzip
import hashlib
import json
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / 'SHA256SUMS.json').read_text(encoding='utf-8'))
    for row in manifest['files']:
        path = root / row['path']
        data = path.read_bytes()
        if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
            raise SystemExit(f'File integrity failure: {row["path"]}')
    summary = json.loads((root / 'data' / 'summary.json').read_text(encoding='utf-8'))
    for section in ('polynomials', 'kernel'):
        record = summary[section]
        stream = gzip.decompress((root / 'data' / record['file']).read_bytes())
        if hashlib.sha256(stream).hexdigest() != record['uncompressed_sha256']:
            raise SystemExit(f'Canonical data hash mismatch: {section}')
        rows = [json.loads(line) for line in stream.splitlines()]
        if len(rows) != record['rows']:
            raise SystemExit(f'Row count mismatch: {section}')
    print(f'PASS: {len(manifest["files"])} packaged file hashes; '
          f'{summary["polynomials"]["rows"]} polynomial records; '
          f'{summary["kernel"]["rows"]} kernel records.')


if __name__ == '__main__':
    main()
