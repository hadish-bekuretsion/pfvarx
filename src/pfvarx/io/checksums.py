"""File provenance: SHA-256, size, mtime and role recording."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class FileProvenance:
    path: str
    role: str
    sha256: str
    size_bytes: int
    mtime_epoch: float
    exists: bool = True

    def as_dict(self) -> dict:
        return asdict(self)


def sha256_file(path: str | Path, chunk_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(path: str | Path, role: str, compute_hash: bool = True) -> FileProvenance:
    p = Path(path)
    if not p.exists():
        return FileProvenance(str(p), role, "", 0, 0.0, exists=False)
    st = p.stat()
    return FileProvenance(
        path=str(p),
        role=role,
        sha256=sha256_file(p) if compute_hash else "",
        size_bytes=st.st_size,
        mtime_epoch=st.st_mtime,
    )


def write_provenance(records: list[FileProvenance], out_path: str | Path) -> Path:
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = [r.as_dict() for r in records]
    tmp = out.with_suffix(out.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2))
    os.replace(tmp, out)
    return out
