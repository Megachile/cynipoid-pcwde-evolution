"""Resolve pinned CAZy inputs strictly below the staged paper data root."""
from pathlib import Path


def resolve_input(data_root: Path, relative: str) -> Path:
    root = data_root.resolve(strict=True)
    rel = Path(relative)
    if not relative or rel.is_absolute() or '..' in rel.parts:
        raise ValueError(f'unsafe CAZy input path: {relative}')
    path = root / rel
    real = path.resolve(strict=True)
    if not real.is_relative_to(root) or not real.is_file():
        raise ValueError(f'CAZy input outside staged data root: {relative}')
    return real
