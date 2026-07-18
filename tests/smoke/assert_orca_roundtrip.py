#!/usr/bin/env python3
"""Validate geometry and layer identity across ORCA export/import."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def read_pdb(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    coords = []
    layers = []
    identities = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            layers.append(float(line[60:66]))
            # Atom/residue identity is the record fields, not trailing padding: the
            # restored PDB writes standard 80-column records while a ragged fixture
            # line may stop short, so compare identity whitespace-insensitively.
            identities.append((line[:30] + line[66:]).rstrip())
    if not coords:
        raise SystemExit(f"no PDB atoms in {path}")
    return np.asarray(coords), np.asarray(layers), identities


def layer_class(values: np.ndarray) -> np.ndarray:
    refs = np.asarray([0.0, 10.0, 20.0])
    return np.argmin(np.abs(values[:, None] - refs[None, :]), axis=1)


source = Path(sys.argv[1])
restored = Path(sys.argv[2])
src_xyz, src_layers, src_identities = read_pdb(source)
dst_xyz, dst_layers, dst_identities = read_pdb(restored)
if src_xyz.shape != dst_xyz.shape:
    raise SystemExit(f"atom-count mismatch: {src_xyz.shape} vs {dst_xyz.shape}")
if not np.allclose(src_xyz, dst_xyz, atol=1.0e-5, rtol=0.0):
    raise SystemExit("ORCA round-trip changed Cartesian coordinates")
if not np.array_equal(layer_class(src_layers), layer_class(dst_layers)):
    raise SystemExit("ORCA round-trip changed ML/movable/frozen layer identity")
if set(layer_class(src_layers).tolist()) != {0, 1, 2}:
    raise SystemExit("ORCA round-trip fixture does not exercise all three layer classes")
if src_identities != dst_identities:
    raise SystemExit("ORCA round-trip changed PDB atom/residue identity")
