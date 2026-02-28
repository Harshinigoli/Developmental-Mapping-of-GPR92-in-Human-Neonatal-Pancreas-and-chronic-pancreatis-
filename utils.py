from __future__ import annotations
from pathlib import Path
import re

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path

def sanitize_filename(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', name).strip('_')
