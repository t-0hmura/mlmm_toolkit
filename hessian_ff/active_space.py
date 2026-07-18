"""M54: ordered active-atom validation + parsed-topology content identity.

Two small, pure pieces shared by hessian_ff's public entry points:

* :func:`validate_active_atoms` — one ordered validator below both
  ``build_analytical_hessian`` and ``workflows.torch_hessian`` so an invalid
  active list (empty / duplicate / negative / ``>= natom`` / float / bool /
  numeric-string) fails with one typed ``ValueError`` *before* any native
  extension is probed or a Hessian block is allocated.  The caller's valid
  order is preserved exactly (no sorting, no dedup-drop).
* :func:`topology_identity` — a content identity (resolved path, SHA-256,
  size, parser/schema version) for a parsed prmtop.  The digest is
  authoritative; a same-path byte replacement produces a different digest, so a
  runtime keyed on it never reuses a stale parsed system.

This is a private mlmm_toolkit v0.3.3 implementation detail; it is not part of
any public PEScape contract and it imports nothing from pdb2reaction.
"""

from __future__ import annotations

import hashlib
import operator
from pathlib import Path
from typing import Any, Dict, Sequence, Union

PathLike = Union[str, Path]

# Bump when the parsed-topology schema changes so old runtime generations are
# not reused across an incompatible parser change.
TOPOLOGY_SCHEMA_VERSION = "hessian_ff.topology/v1"


def validate_active_atoms(natom: int, active_atoms: Sequence[Any]) -> list[int]:
    """Return the ordered, validated active-atom list or raise ``ValueError``.

    Accepts only true integer / index-protocol values.  Rejects bool, float
    (no silent truncation), numeric strings, empty lists, duplicates,
    negatives, and any index ``>= natom``.  The caller's valid order is
    preserved exactly.
    """

    n = int(natom)
    seq = list(active_atoms)
    if len(seq) == 0:
        raise ValueError("active atom list is empty")

    out: list[int] = []
    seen: set[int] = set()
    for value in seq:
        # bool is an int subclass; reject it explicitly so True/False can never
        # be silently read as 1/0.
        if isinstance(value, bool):
            raise ValueError(
                f"active atom index must be an integer, not bool: {value!r}"
            )
        # float would silently truncate; a numeric string is not an index.
        try:
            ia = operator.index(value)
        except TypeError as exc:
            raise ValueError(
                f"active atom index must be an integer, got {value!r}"
            ) from exc
        if ia < 0:
            raise ValueError(f"active atom index is negative: {ia}")
        if ia >= n:
            raise ValueError(f"active atom index out of range: {ia} (natom={n})")
        if ia in seen:
            raise ValueError(f"duplicate active atom index: {ia}")
        seen.add(ia)
        out.append(ia)
    return out


def _file_sha256(path: PathLike) -> str:
    hasher = hashlib.sha256()
    with open(Path(path), "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def topology_identity(prmtop: PathLike) -> Dict[str, Any]:
    """Content identity of a prmtop: resolved path, SHA-256, size, schema.

    The SHA-256 is authoritative; ``size`` is diagnostic / fast-path only.
    """

    resolved = Path(prmtop).resolve()
    return {
        "path": str(resolved),
        "sha256": _file_sha256(resolved),
        "size": int(resolved.stat().st_size),
        "schema_version": TOPOLOGY_SCHEMA_VERSION,
    }
