"""Ordinal PDB atom parsing for topology-to-parm7 index mapping."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

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


__all__ = ["PDBOrdinalAtom", "parse_pdb_ordinal_atoms"]
