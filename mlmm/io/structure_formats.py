"""Format bridge for mmCIF and PDB-size-limit-safe structure handling.

The computational core intentionally continues to consume PDB.  This module
normalizes mmCIF files, and PDB files that exceed fixed-column limits, into a
temporary standards-compliant PDB while retaining the original atom-site
metadata.  Coordinates written by a workflow can then be rendered back to
mmCIF without losing multi-character chain IDs or large residue numbers.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import logging
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any, Iterable, Optional, Sequence

import numpy as np

from mlmm.io.altloc import (
    choose_altloc_label,
    occupancy_choice_key,
    parsed_occupancy,
)


logger = logging.getLogger(__name__)

CIF_SUFFIXES = frozenset({".cif", ".mmcif"})
_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
_MAX_RESIDUES = len(_CHAIN_IDS) * 9999


@dataclass(frozen=True)
class AtomSiteRecord:
    """Metadata for one atom in input-file order."""

    group_pdb: str
    element: str
    atom_name: str
    altloc: str
    resname: str
    chain_id: str
    resseq: str
    icode: str
    occupancy: float
    bfactor: float
    segid: str = ""
    occupancy_known: bool = True
    formal_charge: str = "."
    label_atom_id: str = ""
    label_comp_id: str = ""
    label_asym_id: str = ""
    label_entity_id: str = "1"
    label_seq_id: str = "."
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    @property
    def residue_key(self) -> tuple[str, str, str, str]:
        return (self.chain_id, self.resseq, self.icode, self.segid)

    def with_coordinates(self, xyz: Sequence[float]) -> "AtomSiteRecord":
        return replace(self, x=float(xyz[0]), y=float(xyz[1]), z=float(xyz[2]))


@dataclass(frozen=True)
class CoordinateTemplate:
    """Original atom-site metadata associated with an internal PDB."""

    records: tuple[AtomSiteRecord, ...]
    source_path: Path
    source_format: str
    reason: str

    @property
    def natoms(self) -> int:
        return len(self.records)


def validate_coordinate_template_symbols(
    frame_symbols: Sequence[Sequence[str]],
    template: CoordinateTemplate,
) -> None:
    """Validate each ordered frame element list against retained metadata."""

    expected = tuple(str(record.element).strip().title() for record in template.records)
    for frame_index, symbols in enumerate(frame_symbols, start=1):
        actual = tuple(str(symbol).strip().title() for symbol in symbols)
        if len(actual) != len(expected):
            raise ValueError(
                f"Atom count mismatch between XYZ frame {frame_index} and the retained "
                f"coordinate template: got {len(actual)}, expected {len(expected)}."
            )
        if actual != expected:
            mismatch = next(
                index
                for index, (actual_symbol, expected_symbol) in enumerate(
                    zip(actual, expected), start=1
                )
                if actual_symbol != expected_symbol
            )
            raise ValueError(
                "Ordered elements differ from the retained coordinate template "
                f"in XYZ frame {frame_index} at atom {mismatch}: got "
                f"{actual[mismatch - 1]}, expected {expected[mismatch - 1]}."
            )


_TEMPLATE_REGISTRY: dict[Path, CoordinateTemplate] = {}


def _registry_key(path: Path | str) -> Path:
    return Path(path).expanduser().resolve(strict=False)


def register_coordinate_template(
    pdb_path: Path | str,
    template: CoordinateTemplate,
) -> None:
    """Associate an internal PDB path with its original atom-site metadata."""

    _TEMPLATE_REGISTRY[_registry_key(pdb_path)] = template


def coordinate_template_for(path: Path | str) -> Optional[CoordinateTemplate]:
    """Return the registered template for *path*, if any."""

    return _TEMPLATE_REGISTRY.get(_registry_key(path))


def unregister_coordinate_template(path: Path | str) -> None:
    _TEMPLATE_REGISTRY.pop(_registry_key(path), None)


def is_cif_path(path: Path | str) -> bool:
    return Path(path).suffix.lower() in CIF_SUFFIXES


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def _column(data: dict[str, Any], key: str, nrows: int, default: str) -> list[str]:
    values = _as_list(data.get(key))
    if not values:
        return [default] * nrows
    if len(values) == 1 and nrows > 1:
        return values * nrows
    if len(values) != nrows:
        raise ValueError(
            f"mmCIF column {key!r} has {len(values)} rows; expected {nrows}."
        )
    return values


def _present(value: str) -> bool:
    return bool(value) and value not in {".", "?"}


def _first_present(*values: str, default: str = "") -> str:
    for value in values:
        if _present(value):
            return value
    return default


def _float_or(value: str, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _coherent_altloc_records(
    records: Sequence[AtomSiteRecord],
) -> tuple[list[AtomSiteRecord], int]:
    """Choose one labelled conformer per residue and retain shared atoms."""

    grouped: dict[tuple[str, str, str, str], list[AtomSiteRecord]] = {}
    order: list[tuple[str, str, str, str]] = []
    for record in records:
        key = record.residue_key
        if key not in grouped:
            grouped[key] = []
            order.append(key)
        grouped[key].append(record)

    selected: list[AtomSiteRecord] = []
    removed = 0
    for key in order:
        residue_records = grouped[key]
        labels: list[str] = []
        for record in residue_records:
            label = record.altloc.strip()
            if label and label not in labels:
                labels.append(label)
        # Small-molecule PDB files commonly reuse generic atom names such as
        # ``C`` or ``H`` within one residue.  Duplicate names alone do not
        # denote alternate conformers: only an explicit altloc label does.
        # Preserve the original atom order and multiplicity when no label is
        # present anywhere in the residue.
        if not labels:
            selected.extend(residue_records)
            continue
        chosen = choose_altloc_label(
            (
                record.altloc,
                record.occupancy if record.occupancy_known else None,
                record_index,
            )
            for record_index, record in enumerate(residue_records)
        )

        chosen_by_name: dict[str, list[AtomSiteRecord]] = {}
        for record in residue_records:
            if record.altloc.strip() == chosen:
                chosen_by_name.setdefault(record.atom_name.strip(), []).append(record)

        kept = 0
        emitted_chosen: set[str] = set()
        for record in residue_records:
            label = record.altloc.strip()
            name = record.atom_name.strip()
            if not label:
                # A labelled version of the same atom supersedes the shared
                # record. Otherwise preserve every blank record, including
                # repeated generic atom names in small-molecule PDB files.
                if name not in chosen_by_name:
                    selected.append(record)
                    kept += 1
                continue
            if label != chosen or name in emitted_chosen:
                continue
            candidates = chosen_by_name[name]
            winner = max(
                enumerate(candidates),
                key=lambda item: occupancy_choice_key(
                    item[1].occupancy if item[1].occupancy_known else None,
                    item[0],
                ),
            )[1]
            selected.append(replace(winner, altloc=""))
            emitted_chosen.add(name)
            kept += 1
        removed += len(residue_records) - kept
    return selected, removed


def read_mmcif_atom_sites(path: Path | str) -> list[AtomSiteRecord]:
    """Read the first coordinate model from an mmCIF atom_site loop."""

    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    path = Path(path)
    data = MMCIF2Dict(str(path))
    groups = _as_list(data.get("_atom_site.group_PDB"))
    if not groups:
        raise ValueError(f"No _atom_site rows found in mmCIF file: {path}")
    nrows = len(groups)

    element = _column(data, "_atom_site.type_symbol", nrows, "")
    label_atom = _column(data, "_atom_site.label_atom_id", nrows, "")
    auth_atom = _column(data, "_atom_site.auth_atom_id", nrows, "")
    altloc = _column(data, "_atom_site.label_alt_id", nrows, ".")
    label_comp = _column(data, "_atom_site.label_comp_id", nrows, "UNK")
    auth_comp = _column(data, "_atom_site.auth_comp_id", nrows, "")
    label_asym = _column(data, "_atom_site.label_asym_id", nrows, "")
    auth_asym = _column(data, "_atom_site.auth_asym_id", nrows, "")
    label_entity = _column(data, "_atom_site.label_entity_id", nrows, "1")
    label_seq = _column(data, "_atom_site.label_seq_id", nrows, ".")
    auth_seq = _column(data, "_atom_site.auth_seq_id", nrows, "")
    icodes = _column(data, "_atom_site.pdbx_PDB_ins_code", nrows, "?")
    xs = _column(data, "_atom_site.Cartn_x", nrows, "nan")
    ys = _column(data, "_atom_site.Cartn_y", nrows, "nan")
    zs = _column(data, "_atom_site.Cartn_z", nrows, "nan")
    occupancies = _column(data, "_atom_site.occupancy", nrows, "1.0")
    bfactors = _column(data, "_atom_site.B_iso_or_equiv", nrows, "0.0")
    charges = _column(data, "_atom_site.pdbx_formal_charge", nrows, ".")
    models = _column(data, "_atom_site.pdbx_PDB_model_num", nrows, "1")

    first_model = models[0]
    if any(model != first_model for model in models):
        logger.warning(
            "mmCIF input %s contains multiple coordinate models; using model %s only.",
            path,
            first_model,
        )

    records: list[AtomSiteRecord] = []
    for idx in range(nrows):
        if models[idx] != first_model:
            continue
        xyz = (_float_or(xs[idx], np.nan), _float_or(ys[idx], np.nan), _float_or(zs[idx], np.nan))
        if not np.all(np.isfinite(xyz)):
            raise ValueError(f"Non-finite coordinates in mmCIF atom_site row {idx + 1}: {path}")
        atom_name = _first_present(auth_atom[idx], label_atom[idx], default=element[idx])
        resname = _first_present(auth_comp[idx], label_comp[idx], default="UNK")
        chain = _first_present(auth_asym[idx], label_asym[idx], default="_")
        resseq = _first_present(auth_seq[idx], label_seq[idx], default=str(idx + 1))
        occupancy_value = parsed_occupancy(occupancies[idx])
        record = AtomSiteRecord(
            group_pdb=(groups[idx].upper() if groups[idx] else "HETATM"),
            element=element[idx].strip().title(),
            atom_name=atom_name.strip(),
            altloc=(altloc[idx] if _present(altloc[idx]) else "").strip(),
            resname=resname.strip(),
            chain_id=chain.strip(),
            resseq=resseq.strip(),
            icode=(icodes[idx] if _present(icodes[idx]) else "").strip(),
            occupancy=1.0 if occupancy_value is None else occupancy_value,
            bfactor=_float_or(bfactors[idx], 0.0),
            occupancy_known=occupancy_value is not None,
            formal_charge=charges[idx] if _present(charges[idx]) else ".",
            label_atom_id=label_atom[idx].strip() or atom_name.strip(),
            label_comp_id=label_comp[idx].strip() or resname.strip(),
            label_asym_id=label_asym[idx].strip() or chain.strip(),
            label_entity_id=label_entity[idx] if _present(label_entity[idx]) else "1",
            label_seq_id=label_seq[idx] if _present(label_seq[idx]) else ".",
            x=xyz[0],
            y=xyz[1],
            z=xyz[2],
        )
        if not record.element:
            raise ValueError(f"Missing type_symbol in mmCIF atom_site row {idx + 1}: {path}")
        records.append(record)

    records, removed = _coherent_altloc_records(records)
    if removed:
        logger.warning(
            "mmCIF input %s contained alternate locations; retained one coherent "
            "conformer per residue and removed %d atom-site row(s).",
            path,
            removed,
        )
    return records


_BASE36_UPPER = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
_BASE36_LOWER = "0123456789abcdefghijklmnopqrstuvwxyz"


def _base36_value(text: str, digits: str) -> int:
    value = 0
    for char in text:
        try:
            digit = digits.index(char)
        except ValueError as exc:
            raise ValueError(f"Invalid base-36 digit {char!r} in {text!r}.") from exc
        value = value * 36 + digit
    return value


def _hy36decode(width: int, text: str) -> int:
    """Decode decimal or upper/lower hybrid-36 fields used by large PDBs."""

    stripped = text.strip()
    if not stripped:
        return 0
    try:
        return int(stripped)
    except ValueError:
        pass
    if len(stripped) != width:
        raise ValueError(f"Invalid hybrid-36 value {text!r} for width {width}.")
    if stripped[0] in _BASE36_UPPER[10:]:
        return (
            _base36_value(stripped, _BASE36_UPPER)
            - 10 * (36 ** (width - 1))
            + 10**width
        )
    if stripped[0] in _BASE36_LOWER[10:]:
        return (
            _base36_value(stripped, _BASE36_LOWER)
            + 16 * (36 ** (width - 1))
            + 10**width
        )
    raise ValueError(f"Invalid hybrid-36 value {text!r} for width {width}.")


def _guess_pdb_element(atom_name: str, resname: str, is_hetatm: bool) -> str:
    from mlmm.domain.add_elem_info import guess_element

    value = guess_element(atom_name, resname, is_hetatm)
    return (value or "").title()


def _formal_charge_from_pdb(field: str) -> str:
    value = field.strip()
    if not value:
        return "."
    match = re.fullmatch(r"(?P<magnitude>\d)(?P<sign>[+-])", value)
    if match is None:
        return value
    magnitude = match.group("magnitude")
    return f"-{magnitude}" if match.group("sign") == "-" else magnitude


def _formal_charge_to_pdb(value: str) -> str:
    text = str(value or "").strip()
    if text in {"", ".", "?", "0", "+0", "-0"}:
        return "  "
    if re.fullmatch(r"\d[+-]", text):
        return text
    try:
        charge = int(text)
    except ValueError:
        return "  "
    if not 1 <= abs(charge) <= 9:
        return "  "
    return f"{abs(charge)}{'+' if charge > 0 else '-'}"


def read_pdb_atom_sites(
    path: Path | str,
    *,
    warn_altloc: bool = True,
) -> tuple[list[AtomSiteRecord], bool]:
    """Read PDB atom records, including common decimal-overflow variants."""

    path = Path(path)
    records: list[AtomSiteRecord] = []
    nonstandard = False
    in_first_model = True
    saw_model = False
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("MODEL"):
                if saw_model:
                    in_first_model = False
                else:
                    saw_model = True
                    in_first_model = True
                continue
            if line.startswith("ENDMDL"):
                if saw_model:
                    in_first_model = False
                continue
            if not in_first_model or not line.startswith(("ATOM  ", "HETATM")):
                continue

            serial_shift = 0
            while 11 + serial_shift < len(line) and line[11 + serial_shift].isdigit():
                serial_shift += 1
            if serial_shift:
                nonstandard = True
            offset = serial_shift
            record_name = line[:6].strip().upper()
            atom_name = line[12 + offset : 16 + offset].strip()
            altloc = line[16 + offset : 17 + offset].strip()
            resname = line[17 + offset : 20 + offset].strip()
            chain = line[21 + offset : 22 + offset].strip()

            residue_shift = 0
            resseq_field = line[22 + offset : 26 + offset]
            if 26 + offset < len(line) and line[26 + offset].isdigit():
                overflow_end = 26 + offset
                while overflow_end < len(line) and line[overflow_end].isdigit():
                    overflow_end += 1
                overflow_value = line[22 + offset : overflow_end].strip()
                if re.fullmatch(r"[-+]?\d{5,}", overflow_value):
                    resseq_field = overflow_value
                    residue_shift = overflow_end - (26 + offset)
                    nonstandard = True
            try:
                resseq_value = _hy36decode(4, resseq_field) if not residue_shift else int(resseq_field)
            except ValueError as exc:
                raise ValueError(
                    f"Cannot parse residue number at {path}:{line_number}: {resseq_field!r}"
                ) from exc
            if re.search(r"[A-Za-z]", resseq_field):
                nonstandard = True
            icode = line[26 + offset + residue_shift : 27 + offset + residue_shift].strip()
            coord_offset = offset + residue_shift
            try:
                x = float(line[30 + coord_offset : 38 + coord_offset])
                y = float(line[38 + coord_offset : 46 + coord_offset])
                z = float(line[46 + coord_offset : 54 + coord_offset])
            except ValueError as exc:
                raise ValueError(
                    f"Cannot parse coordinates at {path}:{line_number}."
                ) from exc
            occupancy_value = parsed_occupancy(
                line[54 + coord_offset : 60 + coord_offset].strip()
            )
            occupancy = 1.0 if occupancy_value is None else occupancy_value
            bfactor = _float_or(line[60 + coord_offset : 66 + coord_offset].strip(), 0.0)
            segid = line[72 + coord_offset : 76 + coord_offset].strip()
            element = line[76 + coord_offset : 78 + coord_offset].strip().title()
            formal_charge = _formal_charge_from_pdb(
                line[78 + coord_offset : 80 + coord_offset]
            )
            if not element:
                element = _guess_pdb_element(atom_name, resname, record_name == "HETATM")
            if not element:
                raise ValueError(
                    f"Cannot determine element at {path}:{line_number} ({atom_name!r})."
                )
            serial_field = line[6 : 11 + serial_shift]
            try:
                _hy36decode(5, serial_field)
            except ValueError:
                nonstandard = True
            if re.search(r"[A-Za-z]", serial_field):
                nonstandard = True
            records.append(
                AtomSiteRecord(
                    group_pdb=record_name,
                    element=element,
                    atom_name=atom_name,
                    altloc=altloc,
                    resname=resname,
                    chain_id=chain,
                    resseq=str(resseq_value),
                    icode=icode,
                    occupancy=occupancy,
                    bfactor=bfactor,
                    segid=segid,
                    occupancy_known=occupancy_value is not None,
                    formal_charge=formal_charge,
                    label_atom_id=atom_name,
                    label_comp_id=resname,
                    label_asym_id=chain or "_",
                    label_seq_id=str(resseq_value),
                    x=x,
                    y=y,
                    z=z,
                )
            )
    if not records:
        raise ValueError(f"No ATOM/HETATM records found in PDB file: {path}")
    records, removed = _coherent_altloc_records(records)
    if removed:
        nonstandard = True
        if warn_altloc:
            logger.warning(
                "PDB input %s contained alternate locations; retained one coherent "
                "conformer per residue and removed %d atom record(s) during normalization.",
                path,
                removed,
            )
    return records, nonstandard


def pdb_requires_normalization(path: Path | str) -> bool:
    """Return whether a PDB needs safe internal reindexing."""

    try:
        records, nonstandard = read_pdb_atom_sites(path, warn_altloc=False)
    except ValueError as exc:
        # Preserve the existing validation boundary for empty placeholder
        # files (notably CLI dry-run/help tests). Real geometry loading will
        # still report the missing atoms at the normal workflow stage.
        if "No ATOM/HETATM records" in str(exc):
            return False
        raise
    nres = len(dict.fromkeys(record.residue_key for record in records))
    return nonstandard or nres >= 10_000 or len(records) >= 99_999


def _pdb_atom_field(atom_name: str, element: str) -> str:
    name = atom_name.strip()[:4]
    if len(name) == 4:
        return name
    if len(element.strip()) == 1 and name and name[0].upper() == element[0].upper():
        return f" {name:<3}"
    return f"{name:<4}"


def _base36_name(value: int) -> str:
    """Return a four-character atom identifier for the internal PDB bridge."""

    if not 0 <= value < 36**4:
        raise ValueError("A single residue exceeds the internal PDB atom-name capacity.")
    digits = _BASE36_UPPER
    chars = []
    for _ in range(4):
        value, remainder = divmod(value, 36)
        chars.append(digits[remainder])
    return "".join(reversed(chars))


def _internal_atom_names(records: Sequence[AtomSiteRecord]) -> list[str]:
    """Return PDB-safe names that are unique within every internal residue.

    mmCIF permits atom identifiers that become non-unique after conversion to
    PDB's four-column atom-name field.  Bio.PDB indexes atoms by that field and
    otherwise drops later duplicates.  Preserve the first usable name and give
    only colliding records a bridge-local surrogate; original identifiers stay
    in the coordinate template used for selectors and output restoration.
    """

    candidates = [
        (record.atom_name.strip() or record.element.strip() or "X")[:4]
        for record in records
    ]
    reserved: dict[tuple[str, str, str, str], set[str]] = {}
    for record, candidate in zip(records, candidates):
        reserved.setdefault(record.residue_key, set()).add(candidate)

    used: dict[tuple[str, str, str, str], set[str]] = {}
    counters: dict[tuple[str, str, str, str], int] = {}
    names: list[str] = []
    for record, candidate in zip(records, candidates):
        key = record.residue_key
        residue_used = used.setdefault(key, set())
        if candidate not in residue_used:
            name = candidate
        else:
            counter = counters.get(key, 0)
            while True:
                name = _base36_name(counter)
                counter += 1
                if name not in reserved[key] and name not in residue_used:
                    break
            counters[key] = counter
        residue_used.add(name)
        names.append(name)
    return names


def _internal_atom_serial(atom_index: int) -> int:
    """Return a PDB-safe serial for a zero-based internal atom index."""

    if atom_index < 0:
        raise ValueError("Internal atom index must be non-negative.")
    return atom_index % 99_999 + 1


def write_internal_pdb(records: Sequence[AtomSiteRecord], path: Path | str) -> None:
    """Write records using safe one-character chains and four-digit residue IDs."""

    path = Path(path)
    residue_order = list(dict.fromkeys(record.residue_key for record in records))
    if len(residue_order) > _MAX_RESIDUES:
        raise ValueError(
            f"Structure has {len(residue_order)} residues; internal PDB bridge supports "
            f"at most {_MAX_RESIDUES:,}."
        )
    residue_map: dict[tuple[str, str, str, str], tuple[str, int]] = {}
    for index, key in enumerate(residue_order):
        residue_map[key] = (_CHAIN_IDS[index // 9999], index % 9999 + 1)

    internal_atom_names = _internal_atom_names(records)
    lines: list[str] = []
    previous_chain: Optional[str] = None
    for index, (record, internal_atom_name) in enumerate(
        zip(records, internal_atom_names)
    ):
        chain, resseq = residue_map[record.residue_key]
        if previous_chain is not None and chain != previous_chain:
            lines.append("TER\n")
        previous_chain = chain
        serial = _internal_atom_serial(index)
        record_name = "ATOM  " if record.group_pdb.upper() == "ATOM" else "HETATM"
        atom_field = _pdb_atom_field(internal_atom_name, record.element)
        resname = (record.resname or "UNK")[:3]
        if not (-999.999 <= record.x <= 9999.999 and -999.999 <= record.y <= 9999.999 and -999.999 <= record.z <= 9999.999):
            raise ValueError(
                "Coordinates exceed the fixed-column PDB range required by the internal bridge. "
                "Translate the structure closer to the origin before running mlmm."
            )
        lines.append(
            f"{record_name}{serial:5d} {atom_field}{' ':1}{resname:>3} {chain:1}{resseq:4d}{' ':1}   "
            f"{record.x:8.3f}{record.y:8.3f}{record.z:8.3f}"
            f"{record.occupancy:6.2f}{record.bfactor:6.2f}          "
            f"{record.element:>2s}{_formal_charge_to_pdb(record.formal_charge)}\n"
        )
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def normalize_structure_to_pdb(
    path: Path | str,
) -> tuple[Path, CoordinateTemplate, Path]:
    """Create and register a temporary internal PDB for CIF or oversized PDB input.

    Returns ``(internal_pdb, template, temporary_directory)``.
    """

    source = Path(path).expanduser().resolve()
    if is_cif_path(source):
        records = read_mmcif_atom_sites(source)
        source_format = "mmcif"
        reason = "mmCIF input"
    elif source.suffix.lower() == ".pdb":
        records, _ = read_pdb_atom_sites(source)
        source_format = "pdb"
        reason = "PDB fixed-column size limit"
    else:
        raise ValueError(f"Cannot normalize unsupported structure format: {source}")

    tmp_dir = Path(tempfile.mkdtemp(prefix="mlmm_structure_"))
    internal = tmp_dir / f"{source.stem}.pdb"
    try:
        write_internal_pdb(records, internal)
    except BaseException:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise
    template = CoordinateTemplate(
        records=tuple(records),
        source_path=source,
        source_format=source_format,
        reason=reason,
    )
    register_coordinate_template(internal, template)
    return internal, template, tmp_dir


def cleanup_normalized_structure(path: Path | str, tmp_dir: Optional[Path]) -> None:
    unregister_coordinate_template(path)
    if tmp_dir is not None:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def attach_template_metadata(structure: Any, template: CoordinateTemplate) -> None:
    """Attach original atom/residue identifiers to a Biopython Structure."""

    atoms = list(structure.get_atoms())
    if len(atoms) != template.natoms:
        raise ValueError(
            f"Internal PDB atom count ({len(atoms)}) does not match retained template "
            f"({template.natoms}) from {template.source_path}."
        )
    setattr(structure, "_mlmm_coordinate_template", template)
    for atom, record in zip(atoms, template.records):
        atom.xtra["mlmm_atom_site"] = record
        residue = atom.get_parent()
        residue.xtra.setdefault("mlmm_auth_chain", record.chain_id)
        residue.xtra.setdefault("mlmm_auth_resseq", record.resseq)
        residue.xtra.setdefault("mlmm_auth_icode", record.icode)
        residue.xtra.setdefault("mlmm_auth_resname", record.resname)


def residue_auth_identity(residue: Any) -> tuple[str, str, str, str]:
    """Return original chain, residue number, insertion code, and residue name."""

    chain = str(residue.xtra.get("mlmm_auth_chain", residue.get_parent().id or ""))
    resseq = str(residue.xtra.get("mlmm_auth_resseq", residue.id[1]))
    icode = str(residue.xtra.get("mlmm_auth_icode", residue.id[2] or "")).strip()
    resname = str(residue.xtra.get("mlmm_auth_resname", residue.get_resname()))
    return chain, resseq, icode, resname


def atom_site_from_biopython_atom(atom: Any) -> Optional[AtomSiteRecord]:
    value = atom.xtra.get("mlmm_atom_site")
    if not isinstance(value, AtomSiteRecord):
        return None
    return value.with_coordinates(atom.get_coord())


def template_from_selected_structure(
    structure: Any,
    selected_ids: Iterable[tuple],
    skip_map: dict[tuple, set[str]],
    *,
    link_coordinates: Sequence[Sequence[float]] = (),
) -> Optional[CoordinateTemplate]:
    """Build an output template while preserving selected atoms' original IDs."""

    selected = set(selected_ids)
    records: list[AtomSiteRecord] = []
    source_template = getattr(structure, "_mlmm_coordinate_template", None)
    if not isinstance(source_template, CoordinateTemplate):
        return None
    for atom in structure.get_atoms():
        residue = atom.get_parent()
        fid = residue.get_full_id()
        if fid not in selected or atom.get_name() in skip_map.get(fid, set()):
            continue
        record = atom_site_from_biopython_atom(atom)
        if record is None:
            raise ValueError("Selected atom is missing retained mmCIF/PDB metadata.")
        records.append(record)
    for index, xyz in enumerate(link_coordinates, start=1):
        records.append(
            AtomSiteRecord(
                group_pdb="HETATM",
                element="H",
                atom_name="HL",
                altloc="",
                resname="LKH",
                chain_id="_LINK",
                resseq=str(index),
                icode="",
                occupancy=1.0,
                bfactor=0.0,
                label_atom_id="HL",
                label_comp_id="LKH",
                label_asym_id="_LINK",
                label_entity_id=".",
                label_seq_id=".",
                x=float(xyz[0]),
                y=float(xyz[1]),
                z=float(xyz[2]),
            )
        )
    return CoordinateTemplate(
        records=tuple(records),
        source_path=source_template.source_path,
        source_format=source_template.source_format,
        reason=source_template.reason,
    )


def _cif_quote(value: Any) -> str:
    text = str(value)
    if not text:
        return "."
    if text in {".", "?"}:
        return text
    if re.search(r"\s", text) or text[0] in {"_", "#", "$", ";"} or "'" in text or '"' in text:
        if "'" not in text:
            return f"'{text}'"
        if '"' not in text:
            return f'"{text}"'
        raise ValueError(f"Cannot safely quote mmCIF token: {text!r}")
    return text


_ATOM_SITE_COLUMNS = (
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "pdbx_formal_charge",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_PDB_model_num",
)


def _pdb_coordinate_field(value: float) -> str:
    """Format one finite PDB coordinate without overflowing eight columns."""

    coordinate = float(value)
    if not np.isfinite(coordinate):
        raise ValueError("PDB coordinate frames must contain only finite values.")
    field = f"{coordinate:8.3f}"
    if len(field) != 8:
        raise ValueError(
            "Coordinates exceed the fixed-column PDB range required by the "
            "internal bridge. Translate the structure closer to the origin."
        )
    return field


def render_pdb_coordinate_frames(
    ref_pdb_path: Path | str,
    frame_symbols: Sequence[Sequence[str]],
    frames: Sequence[np.ndarray],
) -> str:
    """Render fully validated coordinate frames on one fixed PDB topology."""

    reference = Path(ref_pdb_path)
    ref_lines = [
        line
        for line in reference.read_text(encoding="utf-8").splitlines(keepends=True)
        if not (line.startswith(("MODEL", "ENDMDL")) or line.strip() == "END")
    ]
    atom_line_indices = [
        index
        for index, line in enumerate(ref_lines)
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    if not atom_line_indices:
        raise ValueError(f"No ATOM/HETATM records in reference PDB: {reference}")

    reference_records, _ = read_pdb_atom_sites(reference, warn_altloc=False)
    if len(reference_records) != len(atom_line_indices):
        raise ValueError(
            "Reference PDB topology contains multiple models or unresolved alternate "
            "locations; normalize it before coordinate overlay."
        )
    expected_symbols = tuple(record.element.title() for record in reference_records)
    if len(frame_symbols) != len(frames):
        raise ValueError("frame_symbols must contain one ordered element list per frame.")
    if not frames:
        raise ValueError("No coordinate frames were provided for PDB rendering.")

    validated: list[np.ndarray] = []
    for frame_index, (symbols, positions) in enumerate(
        zip(frame_symbols, frames), start=1
    ):
        actual_symbols = tuple(str(symbol).strip().title() for symbol in symbols)
        if len(actual_symbols) != len(expected_symbols):
            raise ValueError(
                f"Atom count mismatch in XYZ frame {frame_index}: got "
                f"{len(actual_symbols)}, expected {len(expected_symbols)}."
            )
        if actual_symbols != expected_symbols:
            mismatch = next(
                index
                for index, (actual, expected) in enumerate(
                    zip(actual_symbols, expected_symbols), start=1
                )
                if actual != expected
            )
            raise ValueError(
                f"Ordered elements differ from the reference topology in XYZ frame "
                f"{frame_index} at atom {mismatch}: got "
                f"{actual_symbols[mismatch - 1]}, expected {expected_symbols[mismatch - 1]}."
            )
        array = np.asarray(positions, dtype=float)
        if array.shape != (len(expected_symbols), 3):
            raise ValueError(
                f"Coordinate frame {frame_index} has shape {array.shape}; expected "
                f"({len(expected_symbols)}, 3)."
            )
        if not np.all(np.isfinite(array)):
            raise ValueError(f"XYZ frame {frame_index} contains non-finite coordinates.")
        for value in array.flat:
            _pdb_coordinate_field(float(value))
        validated.append(array)

    atom_line_set = set(atom_line_indices)
    multi_frame = len(validated) > 1
    output: list[str] = []
    for model_number, positions in enumerate(validated, start=1):
        if multi_frame:
            output.append(f"MODEL     {model_number:>4d}\n")
        atom_index = 0
        for line_index, line in enumerate(ref_lines):
            if line_index not in atom_line_set:
                output.append(line)
                continue
            x, y, z = positions[atom_index]
            output.append(
                line[:30]
                + _pdb_coordinate_field(x)
                + _pdb_coordinate_field(y)
                + _pdb_coordinate_field(z)
                + line[54:]
            )
            atom_index += 1
        if output and not output[-1].endswith("\n"):
            output.append("\n")
        if multi_frame:
            output.append("ENDMDL\n")
    return "".join(output)


def render_mmcif_frames(
    frames: Sequence[np.ndarray],
    template: CoordinateTemplate,
    *,
    occupancy_frames: Optional[Sequence[np.ndarray]] = None,
    bfactor_frames: Optional[Sequence[np.ndarray]] = None,
) -> str:
    """Render coordinate frames with retained atom-site metadata."""

    lines = ["data_mlmm\n", "#\n", "loop_\n"]
    lines.extend(f"_atom_site.{column}\n" for column in _ATOM_SITE_COLUMNS)
    if occupancy_frames is not None and len(occupancy_frames) != len(frames):
        raise ValueError("occupancy_frames must contain one array per coordinate frame.")
    if bfactor_frames is not None and len(bfactor_frames) != len(frames):
        raise ValueError("bfactor_frames must contain one array per coordinate frame.")
    atom_id = 1
    for model_number, positions in enumerate(frames, start=1):
        array = np.asarray(positions, dtype=float)
        if array.shape != (template.natoms, 3):
            raise ValueError(
                f"Coordinate frame has shape {array.shape}; expected ({template.natoms}, 3)."
            )
        if not np.all(np.isfinite(array)):
            raise ValueError(
                f"Coordinate frame {model_number} contains non-finite values."
            )
        occupancies = None
        if occupancy_frames is not None:
            occupancies = np.asarray(occupancy_frames[model_number - 1], dtype=float)
            if occupancies.shape != (template.natoms,):
                raise ValueError(
                    f"Occupancy frame has shape {occupancies.shape}; expected ({template.natoms},)."
                )
        bfactors = None
        if bfactor_frames is not None:
            bfactors = np.asarray(bfactor_frames[model_number - 1], dtype=float)
            if bfactors.shape != (template.natoms,):
                raise ValueError(
                    f"B-factor frame has shape {bfactors.shape}; expected ({template.natoms},)."
                )
        for atom_index, (record, xyz) in enumerate(zip(template.records, array)):
            occupancy = record.occupancy if occupancies is None else float(occupancies[atom_index])
            bfactor = record.bfactor if bfactors is None else float(bfactors[atom_index])
            values = (
                record.group_pdb,
                atom_id,
                record.element,
                record.label_atom_id or record.atom_name,
                record.altloc or ".",
                record.label_comp_id or record.resname,
                record.label_asym_id or record.chain_id or "_",
                record.label_entity_id or ".",
                record.label_seq_id or ".",
                record.icode or "?",
                f"{float(xyz[0]):.6f}",
                f"{float(xyz[1]):.6f}",
                f"{float(xyz[2]):.6f}",
                f"{occupancy:.2f}",
                f"{bfactor:.2f}",
                record.formal_charge or ".",
                record.resseq,
                record.resname,
                record.chain_id or "_",
                record.atom_name,
                model_number,
            )
            lines.append(" ".join(_cif_quote(value) for value in values) + "\n")
            atom_id += 1
    lines.append("#\n")
    return "".join(lines)


def write_mmcif_frames(
    frames: Sequence[np.ndarray],
    template: CoordinateTemplate,
    out_path: Path | str,
    *,
    occupancy_frames: Optional[Sequence[np.ndarray]] = None,
    bfactor_frames: Optional[Sequence[np.ndarray]] = None,
) -> None:
    """Write one or more coordinate frames with retained atom-site metadata."""

    content = render_mmcif_frames(
        frames,
        template,
        occupancy_frames=occupancy_frames,
        bfactor_frames=bfactor_frames,
    )
    Path(out_path).write_text(content, encoding="utf-8")


def _pdb_frame_data(
    path: Path | str,
) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    path = Path(path)
    frames: list[list[tuple[float, float, float]]] = []
    occupancy_frames: list[list[float]] = []
    bfactor_frames: list[list[float]] = []
    current: list[tuple[float, float, float]] = []
    current_occupancies: list[float] = []
    current_bfactors: list[float] = []
    saw_model = False

    def finish_frame() -> None:
        nonlocal current, current_occupancies, current_bfactors
        if not current:
            return
        frames.append(current)
        occupancy_frames.append(current_occupancies)
        bfactor_frames.append(current_bfactors)
        current = []
        current_occupancies = []
        current_bfactors = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                finish_frame()
                saw_model = True
                continue
            if line.startswith("ENDMDL"):
                finish_frame()
                continue
            if line.startswith(("ATOM  ", "HETATM")):
                current.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                current_occupancies.append(_float_or(line[54:60].strip(), 1.0))
                current_bfactors.append(_float_or(line[60:66].strip(), 0.0))
    if current or not saw_model:
        finish_frame()
    if not frames:
        raise ValueError(f"No coordinate frames found in PDB file: {path}")
    return (
        [np.asarray(frame, dtype=float) for frame in frames],
        [np.asarray(frame, dtype=float) for frame in occupancy_frames],
        [np.asarray(frame, dtype=float) for frame in bfactor_frames],
    )


def write_pdb_as_mmcif(
    pdb_path: Path | str,
    template: CoordinateTemplate,
    out_path: Path | str,
) -> None:
    frames, occupancies, bfactors = _pdb_frame_data(pdb_path)
    write_mmcif_frames(
        frames,
        template,
        out_path,
        occupancy_frames=occupancies,
        bfactor_frames=bfactors,
    )


def write_xyz_as_mmcif(
    xyz_path: Path | str,
    template: CoordinateTemplate,
    out_path: Path | str,
) -> None:
    from ase.io import read as ase_read

    atoms_frames = ase_read(str(xyz_path), index=":", format="xyz")
    if not atoms_frames:
        raise ValueError(f"No coordinate frames found in XYZ/TRJ file: {xyz_path}")
    write_mmcif_frames(
        [np.asarray(atoms.get_positions(), dtype=float) for atoms in atoms_frames],
        template,
        out_path,
    )


def register_output_template_and_write_cif(
    pdb_path: Path | str,
    template: Optional[CoordinateTemplate],
) -> Optional[Path]:
    """Register an internal PDB output and write its public mmCIF companion."""

    pdb_path = Path(pdb_path)
    out_cif = pdb_path.with_suffix(".cif")
    if template is None:
        unregister_coordinate_template(pdb_path)
        out_cif.unlink(missing_ok=True)
        return None
    try:
        write_pdb_as_mmcif(pdb_path, template, out_cif)
    except BaseException:
        unregister_coordinate_template(pdb_path)
        raise
    register_coordinate_template(pdb_path, template)
    return out_cif
