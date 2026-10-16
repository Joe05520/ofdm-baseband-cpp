#!/usr/bin/env python3
"""Expand losslessly compressed C++ trial CSVs without replacing modified data."""
from pathlib import Path
import gzip
root=Path(__file__).resolve().parents[1]
for src in sorted((root/'results').rglob('*.csv.gz')):
    dst=src.with_suffix('')
    data=gzip.decompress(src.read_bytes())
    if dst.exists() and dst.read_bytes()!=data:
        raise SystemExit(f'Refusing to replace modified evidence: {dst.relative_to(root)}')
    dst.write_bytes(data)
    print(dst.relative_to(root))
