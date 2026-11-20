#!/usr/bin/env python3
"""Deterministic deposit archive and fail-closed checksum/member verification."""
import gzip
import hashlib
import os
import tarfile
from pathlib import Path

EPOCH = 1786924800
NAME = "pcwde_supporting_data"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_archive(outroot: Path) -> Path:
    root = outroot / NAME
    archive = outroot / f"{NAME}.tar.gz"
    receipt = outroot / "ARCHIVE.sha256"
    if not root.is_dir() or archive.exists() or receipt.exists():
        raise SystemExit("unsafe or existing archive output")
    paths = [root, *sorted(root.rglob("*"))]
    if any(path.is_symlink() or not (path.is_file() or path.is_dir()) for path in paths):
        raise SystemExit("archive source contains nonregular member")
    for path in paths:
        os.utime(path, (EPOCH, EPOCH), follow_symlinks=False)
    with archive.open("xb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=EPOCH) as compressed, tarfile.open(fileobj=compressed, mode="w") as tar:
        for path in paths:
            info = tar.gettarinfo(str(path), arcname=path.relative_to(outroot).as_posix())
            info.mode = 0o755 if path.is_dir() else 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            info.mtime = EPOCH
            if path.is_file():
                with path.open("rb") as handle:
                    tar.addfile(info, handle)
            else:
                tar.addfile(info)
    receipt.write_text(f"{sha256(archive)}  {archive.name}\n", encoding="utf-8")
    validate_archive(outroot)
    return archive


def validate_archive(outroot: Path) -> int:
    root = outroot / NAME
    archive = outroot / f"{NAME}.tar.gz"
    receipt = outroot / "ARCHIVE.sha256"
    if not root.is_dir() or archive.is_symlink() or receipt.is_symlink() or not archive.is_file() or not receipt.is_file():
        raise SystemExit("archive or checksum sidecar missing or unsafe")
    if receipt.read_text(encoding="utf-8") != f"{sha256(archive)}  {archive.name}\n":
        raise SystemExit("archive checksum sidecar mismatch")
    expected = {path.relative_to(outroot).as_posix(): path for path in [root, *root.rglob("*")]}
    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        names = [member.name for member in members]
        if len(names) != len(set(names)) or set(names) != set(expected):
            raise SystemExit("archive member inventory mismatch")
        for member in members:
            path = Path(member.name)
            if path.is_absolute() or ".." in path.parts or path.parts[0] != NAME:
                raise SystemExit("archive member path escape")
            if not (member.isfile() or member.isdir()):
                raise SystemExit("archive contains nonregular member")
            original = expected[member.name]
            if member.isfile() != original.is_file() or member.mode != (0o644 if member.isfile() else 0o755) or member.uid != 0 or member.gid != 0 or member.mtime != EPOCH:
                raise SystemExit("archive member metadata mismatch")
            if member.isfile():
                digest = hashlib.sha256()
                with tar.extractfile(member) as handle:
                    for block in iter(lambda: handle.read(1024 * 1024), b""):
                        digest.update(block)
                if digest.hexdigest() != sha256(original):
                    raise SystemExit("archive member checksum mismatch")
    return len(members)
