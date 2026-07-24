"""Persistent atom-order identity for ONIOM export/import round trips."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Sequence


MARKER_NAME = "MLMM_REF_PDB_ORDER_V1_SHA256"
_MARKER_RE = re.compile(
    rf"(?<![A-Za-z0-9_]){MARKER_NAME}=([0-9a-f]{{64}})(?![A-Za-z0-9_])"
)
_DIGEST_DOMAIN = b"mlmm-ref-pdb-order/v1\0"


def _pdb_identity_rows(path: Path | str) -> list[list[str]]:
    rows: list[list[str]] = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        padded = line.ljust(80)
        rows.append(
            [
                padded[0:6].strip(),
                padded[6:11].strip(),
                padded[12:16],
                padded[16:17],
                padded[17:20].strip(),
                padded[21:22],
                padded[22:26].strip(),
                padded[26:27],
                padded[72:76].strip(),
                padded[76:78].strip().title(),
            ]
        )
    if not rows:
        raise ValueError(f"Reference PDB has no ATOM/HETATM rows: {path}")
    return rows


def pdb_order_digest(path: Path | str) -> str:
    """Hash ordered PDB atom identities while excluding mutable geometry fields."""
    payload = json.dumps(
        _pdb_identity_rows(path),
        ensure_ascii=True,
        separators=(",", ":"),
    ).encode("ascii")
    return hashlib.sha256(_DIGEST_DOMAIN + payload).hexdigest()


def format_order_marker(digest: str) -> str:
    """Return the versioned marker after validating its digest."""
    value = str(digest)
    if re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError("ONIOM reference-order digest must be 64 lowercase hex characters.")
    return f"{MARKER_NAME}={value}"


def extract_embedded_order_digest(path: Path | str) -> str | None:
    """Read one valid embedded marker, rejecting malformed or duplicate markers."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    if MARKER_NAME not in text:
        return None
    matches = _MARKER_RE.findall(text)
    if len(matches) != 1:
        detail = "malformed" if not matches else "duplicated or conflicting"
        raise ValueError(f"ONIOM reference-order marker is {detail}.")
    return matches[0]


def elements_are_unique(elements: Sequence[str]) -> bool:
    """Return whether ordered element symbols alone identify every atom."""
    normalized = [str(element).strip().title() for element in elements]
    return len(normalized) == len(set(normalized))
