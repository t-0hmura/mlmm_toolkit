"""Ordinal PDB parsing and fail-closed ML/MM atom identity resolution."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Sequence

from mlmm.domain.add_elem_info import guess_element
from mlmm.io.structure_formats import _hy36decode


@dataclass(frozen=True)
class PDBOrdinalAtom:
    """One PDB atom with distinct file-order and deposited serial indices."""

    idx: int
    serial: int | str
    record_name: str
    atom_name: str
    altloc: str
    resname: str
    chain_id: str
    resseq: str
    icode: str
    elem: str
    coord: tuple[float, float, float]

    @property
    def id(self) -> str:
        """Legacy ML-region identity used by the existing selection contract."""

        return f"{self.atom_name} {self.resname} {self.resseq}"

    @property
    def structural_key(self) -> tuple[str, str, str, str, str, str]:
        """Identity retained by prepared PDB files, excluding serial/coordinates."""

        return (
            self.chain_id,
            self.resname,
            self.resseq,
            self.icode,
            self.atom_name,
            self.altloc,
        )


@dataclass(frozen=True)
class ResolvedMLMMAtoms:
    """One canonical mapping shared by ML/MM calculators and DFT workflows."""

    full_atoms: tuple[PDBOrdinalAtom, ...]
    model_indices: tuple[int, ...]
    link_pairs: tuple[tuple[int, int], ...]
    link_element_pairs: tuple[tuple[str, str], ...]


# CHEMISTRY-RULE:9 parm7 1-based atom indexing (NOT PDB serial).
# Use the 1-based ATOM/HETATM file position as `idx`. parm7 atoms (parmed) are
# renumbered 1..N sequentially after tleap, so PDB serial fields cannot be used
# as parm7 indices because serial gaps break `real.atoms[idx-1]` lookups.
def parse_pdb_ordinal_atoms(path: Path | str) -> list[PDBOrdinalAtom]:
    """Parse PDB atoms with one-based file ordinals independent of serial gaps."""

    source = Path(path)
    atoms: list[PDBOrdinalAtom] = []
    with source.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            serial_field = line[6:11]
            try:
                serial: int | str = _hy36decode(5, serial_field)
            except ValueError:
                # Serial is diagnostic metadata, never the parm7 index.  Keep
                # an unusual deposited token without blocking ordinal mapping.
                serial = serial_field.strip()
            try:
                coord = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            except ValueError as exc:
                raise ValueError(
                    f"Invalid PDB coordinates at {source}:{line_number}."
                ) from exc
            if not all(math.isfinite(value) for value in coord):
                raise ValueError(
                    f"Non-finite PDB coordinates at {source}:{line_number}."
                )

            record_name = line[:6].strip().upper()
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            elem = line[76:78].strip().title()
            if not elem:
                elem = str(
                    guess_element(atom_name, resname, record_name == "HETATM") or ""
                ).title()
            if not elem:
                raise ValueError(
                    f"Cannot determine PDB element at {source}:{line_number}."
                )
            atoms.append(
                PDBOrdinalAtom(
                    idx=len(atoms) + 1,
                    serial=serial,
                    record_name=record_name,
                    atom_name=atom_name,
                    altloc=line[16:17].strip(),
                    resname=resname,
                    chain_id=line[21:22].strip(),
                    resseq=line[22:26].strip(),
                    icode=line[26:27].strip(),
                    elem=elem,
                    coord=coord,
                )
            )
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in the input PDB: {source}")
    return atoms


def _selector_metadata(
    atoms: Sequence[PDBOrdinalAtom],
) -> list[dict[str, str | int]]:
    metadata: list[dict[str, str | int]] = []
    for atom in atoms:
        try:
            resseq: str | int = int(atom.resseq)
        except ValueError:
            resseq = atom.resseq
        metadata.append(
            {
                "name": atom.atom_name,
                "resname": atom.resname,
                "resseq": resseq,
                "chain": atom.chain_id,
                "icode": atom.icode,
                "element": atom.elem,
            }
        )
    return metadata


def _resolve_model_indices(
    full_atoms: Sequence[PDBOrdinalAtom],
    model_atoms: Sequence[PDBOrdinalAtom],
) -> tuple[int, ...]:
    for source_name, atoms in (("input PDB", full_atoms), ("model PDB", model_atoms)):
        nonblank = [atom for atom in atoms if atom.altloc]
        if nonblank:
            first = nonblank[0]
            raise ValueError(
                f"{source_name} contains nonblank altLoc '{first.altloc}' at "
                f"{first.chain_id or '-'}:{first.resname}:{first.resseq}"
                f"{first.icode}:{first.atom_name}; normalize alternate conformers "
                "before resolving the ML region."
            )

    by_key: dict[tuple[str, str, str, str, str, str], list[int]] = {}
    for atom in full_atoms:
        by_key.setdefault(atom.structural_key, []).append(atom.idx)

    resolved: list[int] = []
    seen_model: set[tuple[str, str, str, str, str, str]] = set()
    seen_full: set[int] = set()
    for model_atom in model_atoms:
        key = model_atom.structural_key
        label = (
            f"{model_atom.chain_id or '-'}:{model_atom.resname}:"
            f"{model_atom.resseq}{model_atom.icode}:{model_atom.atom_name}"
        )
        if key in seen_model:
            raise ValueError(f"model_pdb contains duplicate atom identity {label}.")
        seen_model.add(key)
        matches = by_key.get(key, [])
        if len(matches) != 1:
            if not matches:
                raise ValueError(
                    f"model_pdb atom {label} is absent from the input PDB."
                )
            raise ValueError(
                f"model_pdb atom {label} matches {len(matches)} input atoms; "
                "retain unique chain and insertion-code identifiers."
            )
        idx = matches[0]
        full_atom = full_atoms[idx - 1]
        if full_atom.elem != model_atom.elem:
            raise ValueError(
                f"model_pdb atom {label} has element {model_atom.elem}, but the "
                f"matching input atom has element {full_atom.elem}."
            )
        if idx in seen_full:
            raise ValueError(
                f"Input atom ordinal {idx} is selected more than once by model_pdb."
            )
        seen_full.add(idx)
        resolved.append(idx)

    if any(right <= left for left, right in zip(resolved, resolved[1:])):
        raise ValueError(
            "model_pdb atom order differs from the input PDB; preserve full-system "
            "file order when creating the ML-region PDB."
        )
    return tuple(resolved)


def _resolve_manual_links(
    full_atoms: Sequence[PDBOrdinalAtom],
    model_indices: Sequence[int],
    manual_links: Sequence[Sequence[str]],
) -> tuple[tuple[int, int], ...]:
    # Imported lazily to keep structure-format imports acyclic.
    from mlmm.core.utils import resolve_atom_spec_index

    metadata = _selector_metadata(full_atoms)
    model_set = set(model_indices)
    pairs: list[tuple[int, int]] = []
    used_ml: set[int] = set()
    used_mm: set[int] = set()
    for pair_number, raw_pair in enumerate(manual_links, start=1):
        if len(raw_pair) != 2:
            raise ValueError(
                f"link_mlmm pair {pair_number} must contain exactly one ML and one "
                "MM atom selector."
            )
        ml_spec, mm_spec = (str(raw_pair[0]), str(raw_pair[1]))
        ml_idx = resolve_atom_spec_index(ml_spec, metadata) + 1
        mm_idx = resolve_atom_spec_index(mm_spec, metadata) + 1
        if ml_idx not in model_set:
            raise ValueError(
                f"link_mlmm pair {pair_number}: ML selector '{ml_spec}' is outside "
                "model_pdb."
            )
        if mm_idx in model_set:
            raise ValueError(
                f"link_mlmm pair {pair_number}: MM selector '{mm_spec}' is inside "
                "model_pdb."
            )
        if ml_idx in used_ml or mm_idx in used_mm:
            raise ValueError(
                "link_mlmm reuses an ML or MM endpoint; each boundary atom must "
                "appear in exactly one manual pair."
            )
        used_ml.add(ml_idx)
        used_mm.add(mm_idx)
        pairs.append((ml_idx, mm_idx))
    return tuple(pairs)


def _detect_link_pairs(
    full_atoms: Sequence[PDBOrdinalAtom],
    model_indices: Sequence[int],
) -> tuple[tuple[int, int], ...]:
    threshold = 1.7
    model_set = set(model_indices)
    coords = {atom.idx: atom.coord for atom in full_atoms}
    elems = {atom.idx: atom.elem for atom in full_atoms}
    pairs: list[tuple[int, int]] = []
    for ml_idx in model_indices:
        qcoord = coords[ml_idx]
        for atom in full_atoms:
            mm_idx = atom.idx
            if mm_idx in model_set:
                continue
            delta = (
                qcoord[0] - atom.coord[0],
                qcoord[1] - atom.coord[1],
                qcoord[2] - atom.coord[2],
            )
            distance = math.sqrt(sum(value * value for value in delta))
            if distance >= threshold:
                continue
            if (elems[ml_idx], elems[mm_idx]) in {
                ("C", "C"),
                ("C", "N"),
                ("N", "C"),
            }:
                pairs.append((ml_idx, mm_idx))

    ml_endpoints = [ml for ml, _ in pairs]
    mm_endpoints = [mm for _, mm in pairs]
    if len(set(ml_endpoints)) != len(ml_endpoints) or len(
        set(mm_endpoints)
    ) != len(mm_endpoints):
        raise ValueError(
            "Automatic link detection found a boundary atom in multiple pairs; "
            "specify link_mlmm manually."
        )
    return tuple(pairs)


def resolve_mlmm_atoms(
    input_pdb: Path | str,
    model_pdb: Path | str,
    manual_links: Sequence[Sequence[str]] | None = None,
) -> ResolvedMLMMAtoms:
    """Resolve the model region and boundary links once, without fuzzy fallback."""

    full_atoms = tuple(parse_pdb_ordinal_atoms(input_pdb))
    model_atoms = tuple(parse_pdb_ordinal_atoms(model_pdb))
    model_indices = _resolve_model_indices(full_atoms, model_atoms)
    link_pairs = (
        _resolve_manual_links(full_atoms, model_indices, manual_links)
        if manual_links is not None
        else _detect_link_pairs(full_atoms, model_indices)
    )
    elem_by_idx = {atom.idx: atom.elem for atom in full_atoms}
    link_element_pairs = tuple(
        (elem_by_idx[ml_idx], elem_by_idx[mm_idx])
        for ml_idx, mm_idx in link_pairs
    )
    return ResolvedMLMMAtoms(
        full_atoms=full_atoms,
        model_indices=model_indices,
        link_pairs=link_pairs,
        link_element_pairs=link_element_pairs,
    )


__all__ = [
    "PDBOrdinalAtom",
    "ResolvedMLMMAtoms",
    "parse_pdb_ordinal_atoms",
    "resolve_mlmm_atoms",
]
