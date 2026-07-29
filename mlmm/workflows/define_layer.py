"""Define 3-layer ML/MM system based on distance from the ML region.

Example:
    mlmm define-layer -i system.pdb --model-pdb ml_region.pdb -o labeled.pdb

For detailed documentation, see: docs/define_layer.md
"""

from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import sys

import click
from mlmm.core.output import emit
import numpy as np

from mlmm.core.defaults import (
    BFACTOR_ML,
    BFACTOR_MOVABLE_MM,
    BFACTOR_FROZEN,
)
from mlmm.core.utils import prepare_input_structure
from mlmm.io.structure_formats import (
    coordinate_template_for,
    register_output_template_and_write_cif,
)

# Default radii for layer definition
# NOTE: reserved in 3-layer mode (no-op).
DEFAULT_RADIUS_PARTIAL_HESSIAN = 0.0  # Å
DEFAULT_RADIUS_FREEZE = 8.0  # Å


def _parse_pdb_atoms(pdb_path: Path) -> List[Dict[str, Any]]:
    """
    Parse ATOM/HETATM records from a PDB file.

    Returns a list of dictionaries with atom information:
    - idx: 0-based atom index
    - atom_name: atom name
    - res_name: residue name
    - chain_id: chain identifier
    - res_seq: residue sequence number
    - icode: insertion code
    - coord: (x, y, z) coordinates as numpy array
    - element: element symbol
    - res_key: (chain_id, res_seq, icode) tuple for residue identification
    """
    atoms = []
    with open(pdb_path, "r") as f:
        atom_idx = 0
        for line in f:
            # Reset atom_idx on MODEL records to match write_layered_pdb semantics
            # (single-model PDBs are still handled because the very first MODEL
            # leaves atom_idx==0). Stop after the first ENDMDL so multi-model
            # inputs are interpreted as model-1 only.
            if line.startswith("MODEL"):
                atom_idx = 0
                continue
            if line.startswith("ENDMDL"):
                if atoms:
                    break
                continue
            if not line.startswith(("ATOM  ", "HETATM")):
                continue

            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22].strip()
            res_seq_str = line[22:26].strip()
            try:
                res_seq = int(res_seq_str)
            except ValueError:
                res_seq = 0
            icode = line[26:27].strip()

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                x, y, z = 0.0, 0.0, 0.0

            element = line[76:78].strip() if len(line) >= 78 else ""

            atoms.append({
                "idx": atom_idx,
                "atom_name": atom_name,
                "res_name": res_name,
                "chain_id": chain_id,
                "res_seq": res_seq,
                "icode": icode,
                "coord": np.array([x, y, z]),
                "element": element,
                "res_key": (chain_id, res_seq, icode),
            })
            atom_idx += 1

    return atoms


def _get_ml_indices_from_model_pdb(
    input_atoms: List[Dict[str, Any]],
    model_pdb_path: Path,
    input_pdb_path: Optional[Path] = None,
) -> List[int]:
    """
    Get ML region atom indices by matching atoms from model_pdb to input atoms.

    Matching is done by (chain_id, res_seq, icode, res_name, atom_name) tuple
    so that multi-chain inputs with overlapping res_seq do not collide. If the
    model PDB has a blank chain_id field, we fall back to (res_name, res_seq,
    icode, atom_name) for that entry (single-chain model_pdb describing a
    multi-chain input).
    """
    # Prefer retained author identifiers when either side came through the
    # mmCIF/large-PDB bridge. Internal one-character chain IDs and 4-column
    # residue numbers are implementation details and may differ for a subset.
    input_template = (
        coordinate_template_for(input_pdb_path)
        if input_pdb_path is not None
        else None
    )
    model_template = coordinate_template_for(model_pdb_path)
    if input_template is not None or model_template is not None:
        def _identity(record):
            if isinstance(record, dict):
                return (
                    str(record["chain_id"]),
                    str(record["res_seq"]),
                    str(record["icode"]),
                    str(record["res_name"]),
                    str(record["atom_name"]),
                )
            return (
                str(record.chain_id),
                str(record.resseq),
                str(record.icode),
                str(record.resname),
                str(record.atom_name),
            )

        input_records = (
            input_template.records
            if input_template is not None
            else input_atoms
        )
        model_records = (
            model_template.records
            if model_template is not None
            else _parse_pdb_atoms(model_pdb_path)
        )
        input_chain_counts = Counter(_identity(record) for record in input_records)
        input_blank_counts = Counter(
            _identity(record)[1:] for record in input_records
        )
        unmatched = []
        for record in model_records:
            key = _identity(record)
            available = input_chain_counts if key[0] else input_blank_counts
            lookup = key if key[0] else key[1:]
            if available[lookup] <= 0:
                unmatched.append(key)
            else:
                available[lookup] -= 1
        model_chain_keys = set()
        model_blank_keys = set()
        for record in model_records:
            key = _identity(record)
            if key[0]:
                model_chain_keys.add(key)
            else:
                model_blank_keys.add(key[1:])
        matched = []
        for index, record in enumerate(input_records):
            key = _identity(record)
            if key in model_chain_keys or key[1:] in model_blank_keys:
                matched.append(index)
        if matched:
            if unmatched:
                raise ValueError(
                    "Model PDB atom identity is absent from the full input: "
                    f"{unmatched[0]!r}."
                )
            return matched
        # An extracted model PDB can deliberately use the internal identifiers
        # of a normalized full structure.  If author-identifier matching found
        # nothing, retain the established internal-PDB fallback below.

    # Parse model PDB
    model_atoms = _parse_pdb_atoms(model_pdb_path)

    # Create sets of ML atom identifiers (chain-aware + chain-blank fallback)
    ml_ids_chain = set()
    ml_ids_blank = set()
    for atom in model_atoms:
        chain = atom["chain_id"]
        key = (chain, atom["res_seq"], atom["icode"], atom["res_name"], atom["atom_name"])
        if chain:
            ml_ids_chain.add(key)
        else:
            ml_ids_blank.add(
                (atom["res_seq"], atom["icode"], atom["res_name"], atom["atom_name"])
            )

    input_chain_counts = Counter(
        (
            atom["chain_id"],
            atom["res_seq"],
            atom["icode"],
            atom["res_name"],
            atom["atom_name"],
        )
        for atom in input_atoms
    )
    input_blank_counts = Counter(
        (
            atom["res_seq"],
            atom["icode"],
            atom["res_name"],
            atom["atom_name"],
        )
        for atom in input_atoms
    )
    unmatched = []
    for atom in model_atoms:
        chain = atom["chain_id"]
        key = (
            chain,
            atom["res_seq"],
            atom["icode"],
            atom["res_name"],
            atom["atom_name"],
        )
        available = input_chain_counts if chain else input_blank_counts
        lookup = key if chain else key[1:]
        if available[lookup] <= 0:
            unmatched.append(key)
        else:
            available[lookup] -= 1
    if unmatched:
        raise ValueError(
            "Model PDB atom identity is absent from the full input: "
            f"{unmatched[0]!r}."
        )

    # Find matching atoms in input
    ml_indices = []
    for atom in input_atoms:
        key_chain = (
            atom["chain_id"], atom["res_seq"], atom["icode"],
            atom["res_name"], atom["atom_name"],
        )
        key_blank = (
            atom["res_seq"], atom["icode"], atom["res_name"], atom["atom_name"],
        )
        if key_chain in ml_ids_chain or key_blank in ml_ids_blank:
            ml_indices.append(atom["idx"])

    return sorted(ml_indices)


def _parse_indices_string(
    indices_str: str,
    one_based: bool,
) -> List[int]:
    """Parse comma-separated indices string into a list of 0-based integers."""
    indices = []
    for token in indices_str.split(","):
        token = token.strip()
        if not token:
            continue

        # Handle range notation (e.g., "1-5")
        if "-" in token and not token.startswith("-"):
            parts = token.split("-")
            if len(parts) == 2:
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                    for i in range(start, end + 1):
                        idx = i - 1 if one_based else i
                        indices.append(idx)
                    continue
                except ValueError:
                    pass

        try:
            idx = int(token)
            if one_based:
                idx -= 1
            indices.append(idx)
        except ValueError:
            raise click.BadParameter(f"Invalid index: {token}")

    return sorted(set(indices))


def compute_layer_indices_by_residue(
    atoms: List[Dict[str, Any]],
    ml_indices: List[int],
    radius_partial_hessian: float,
    radius_freeze: float,
) -> Dict[str, List[int]]:
    """
    Compute 3-layer indices based on distance from ML region.

    For residues WITHOUT ML atoms:
        The minimum distance from any atom in the residue to any ML atom is computed.
        The entire residue is assigned to the layer based on this minimum distance.

    For residues WITH ML atoms:
        ML atoms stay in Layer 1. Non-ML atoms (e.g., backbone) are classified
        individually based on their distance to the nearest ML atom.

    Parameters
    ----------
    atoms : List[Dict[str, Any]]
        List of atom dictionaries from _parse_pdb_atoms
    ml_indices : List[int]
        0-based indices of ML region atoms
    radius_partial_hessian : float
        Deprecated in 3-layer mode. Ignored by the 3-layer assignment.
    radius_freeze : float
        Distance cutoff (Å) for movable MM atoms

    Returns
    -------
    Dict[str, List[int]]
        Dictionary with keys: ml_indices, hess_mm_indices (empty in 3-layer mode),
        movable_mm_indices, frozen_indices
    """
    ml_set = set(ml_indices)
    n_atoms = len(atoms)

    # Get ML region coordinates
    ml_coords = np.array([atoms[i]["coord"] for i in ml_indices if i < n_atoms])
    if ml_coords.size == 0:
        raise ValueError("No valid ML atoms found in the input structure")

    # Group atoms by residue
    residue_atoms: Dict[Tuple, List[int]] = {}
    for atom in atoms:
        res_key = atom["res_key"]
        if res_key not in residue_atoms:
            residue_atoms[res_key] = []
        residue_atoms[res_key].append(atom["idx"])

    # Classify atoms based on distance to ML region
    # For residues WITHOUT ML atoms: classify entire residue by minimum distance
    # For residues WITH ML atoms: classify non-ML atoms individually by distance
    movable_mm_indices: List[int] = []
    frozen_indices: List[int] = []

    for atom_indices in residue_atoms.values():
        has_ml_atoms = any(idx in ml_set for idx in atom_indices)

        if has_ml_atoms:
            # This residue contains ML atoms
            # Classify non-ML atoms individually by their distance to ML region
            for idx in atom_indices:
                if idx in ml_set:
                    # Already an ML atom, skip
                    continue
                # Compute distance from this atom to nearest ML atom
                coord = atoms[idx]["coord"]
                dists = np.linalg.norm(ml_coords - coord, axis=1)
                min_dist = np.min(dists)

                if min_dist <= radius_freeze:
                    movable_mm_indices.append(idx)
                else:
                    frozen_indices.append(idx)
        else:
            # No ML atoms in this residue
            # Compute minimum distance from any atom in this residue to any ML atom
            residue_coords = np.array([atoms[i]["coord"] for i in atom_indices])
            min_dist = float("inf")
            for coord in residue_coords:
                dists = np.linalg.norm(ml_coords - coord, axis=1)
                min_dist = min(min_dist, np.min(dists))

            # Assign entire residue to a layer
            if min_dist <= radius_freeze:
                movable_mm_indices.extend(atom_indices)
            else:
                frozen_indices.extend(atom_indices)

    return {
        "ml_indices": sorted(ml_indices),
        # Compatibility key for existing downstream codepaths.
        "hess_mm_indices": [],
        "movable_mm_indices": sorted(movable_mm_indices),
        "frozen_indices": sorted(frozen_indices),
    }


def write_layered_pdb(
    input_pdb_path: Path,
    output_pdb_path: Path,
    layer_indices: Dict[str, List[int]],
) -> None:
    """
    Write a PDB file with B-factors set according to layer assignments.

    Parameters
    ----------
    input_pdb_path : Path
        Input PDB file
    output_pdb_path : Path
        Output PDB file
    layer_indices : Dict[str, List[int]]
        Dictionary with ml_indices, movable_mm_indices, frozen_indices
    """
    ml_set = set(layer_indices["ml_indices"])
    movable_mm_set = set(layer_indices["movable_mm_indices"])
    frozen_set = set(layer_indices["frozen_indices"])

    lines_out = []
    atom_idx = 0

    with open(input_pdb_path, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                atom_idx = 0
                lines_out.append(line)
                continue

            if line.startswith(("ATOM  ", "HETATM")):
                # Determine B-factor for this atom
                if atom_idx in ml_set:
                    bfac = BFACTOR_ML
                elif atom_idx in movable_mm_set:
                    bfac = BFACTOR_MOVABLE_MM
                elif atom_idx in frozen_set:
                    bfac = BFACTOR_FROZEN
                else:
                    # Default: treat as movable MM
                    bfac = BFACTOR_MOVABLE_MM

                # Replace B-factor (columns 61-66)
                if len(line) >= 66:
                    new_line = line[:60] + f"{bfac:6.2f}" + line[66:]
                else:
                    padded = line.rstrip("\n").ljust(66)
                    new_line = padded[:60] + f"{bfac:6.2f}" + "\n"
                lines_out.append(new_line)
                atom_idx += 1
            else:
                lines_out.append(line)

    with open(output_pdb_path, "w") as f:
        f.writelines(lines_out)


def _define_layers_pdb(
    input_pdb: Path,
    output_pdb: Path,
    model_pdb: Optional[Path] = None,
    model_indices: Optional[List[int]] = None,
    radius_partial_hessian: float = DEFAULT_RADIUS_PARTIAL_HESSIAN,
    radius_freeze: float = DEFAULT_RADIUS_FREEZE,
) -> Dict[str, List[int]]:
    """
    Main function to define 3-layer ML/MM system and write output PDB.

    Parameters
    ----------
    input_pdb : Path
        Input PDB file with full system
    output_pdb : Path
        Output PDB file with B-factors set to layer values
    model_pdb : Optional[Path]
        PDB file defining ML region atoms
    model_indices : Optional[List[int]]
        Explicit 0-based indices of ML region atoms (takes precedence over model_pdb)
    radius_partial_hessian : float
        Deprecated in 3-layer mode (ignored).
    radius_freeze : float
        Distance cutoff (Å) for movable MM atoms (default: 8.0)

    Returns
    -------
    Dict[str, List[int]]
        Layer indices dictionary
    """
    # Parse input PDB
    atoms = _parse_pdb_atoms(input_pdb)
    if not atoms:
        raise ValueError(f"No atoms found in {input_pdb}")

    # Get ML indices
    if model_indices is not None:
        ml_indices = sorted(set(model_indices))
    elif model_pdb is not None:
        ml_indices = _get_ml_indices_from_model_pdb(
            atoms, model_pdb, input_pdb_path=input_pdb
        )
    else:
        raise ValueError("Either model_pdb or model_indices must be provided")

    if not ml_indices:
        raise ValueError("No ML region atoms found")

    # Validate indices
    n_atoms = len(atoms)
    invalid = [i for i in ml_indices if i < 0 or i >= n_atoms]
    if invalid:
        raise ValueError(f"Invalid ML indices (out of range 0-{n_atoms-1}): {invalid}")

    # Compute layer indices
    layer_indices = compute_layer_indices_by_residue(
        atoms,
        ml_indices,
        radius_partial_hessian,
        radius_freeze,
    )

    # Write output PDB
    write_layered_pdb(input_pdb, output_pdb, layer_indices)

    return layer_indices


def define_layers(
    input_pdb: Path,
    output_pdb: Path,
    model_pdb: Optional[Path] = None,
    model_indices: Optional[List[int]] = None,
    radius_partial_hessian: float = DEFAULT_RADIUS_PARTIAL_HESSIAN,
    radius_freeze: float = DEFAULT_RADIUS_FREEZE,
) -> Dict[str, List[int]]:
    """Define layers through the PDB/mmCIF normalization bridge."""
    prepared_input = prepare_input_structure(Path(input_pdb))
    prepared_model = (
        prepare_input_structure(Path(model_pdb)) if model_pdb is not None else None
    )
    output_pdb = _effective_output_pdb(output_pdb)
    try:
        layer_indices = _define_layers_pdb(
            input_pdb=prepared_input.source_path,
            output_pdb=output_pdb,
            model_pdb=(
                prepared_model.source_path if prepared_model is not None else None
            ),
            model_indices=model_indices,
            radius_partial_hessian=radius_partial_hessian,
            radius_freeze=radius_freeze,
        )
        register_output_template_and_write_cif(
            output_pdb, prepared_input.structure_template
        )
        return layer_indices
    finally:
        prepared_input.cleanup()
        if prepared_model is not None:
            prepared_model.cleanup()


def _effective_output_pdb(output_path: Path) -> Path:
    """Return the actual PDB path used by the computational bridge."""

    output_path = Path(output_path)
    if output_path.suffix.lower() in {".cif", ".mmcif"}:
        return output_path.with_suffix(".pdb")
    return output_path



@click.command(
    help="Define 3-layer ML/MM system based on distance from ML region.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input PDB or mmCIF file containing the full system.",
)
@click.option(
    "--model-pdb",
    "model_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="PDB or mmCIF file defining atoms in the ML region.",
)
@click.option(
    "--model-indices",
    "model_indices_str",
    type=str,
    default=None,
    help="Comma-separated atom indices for ML region (e.g., '1,2,3,4' or '1-10,15,20-25'). "
         "Takes precedence over --model-pdb.",
)
@click.option(
    "--radius-partial-hessian",
    "radius_partial_hessian",
    type=float,
    default=DEFAULT_RADIUS_PARTIAL_HESSIAN,
    show_default=True,
    help="Deprecated in 3-layer mode (ignored).",
)
@click.option(
    "--radius-freeze",
    "radius_freeze",
    type=float,
    default=DEFAULT_RADIUS_FREEZE,
    show_default=True,
    help="Distance cutoff (Å) from ML region for movable MM atoms. "
         "Atoms beyond this distance are frozen.",
)
@click.option(
    "-o", "--output",
    "output_pdb",
    type=click.Path(path_type=Path, dir_okay=False),
    default=None,
    help="Output PDB file with B-factors set to layer values. "
         "Defaults to '<input>_layered.pdb'.",
)
@click.option(
    "--one-based/--zero-based",
    "one_based",
    default=True,
    show_default=True,
    help="Interpret --model-indices as 1-based (default) or 0-based.",
)
def cli(
    input_pdb: Path,
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    radius_partial_hessian: float,
    radius_freeze: float,
    output_pdb: Optional[Path],
    one_based: bool,
) -> None:
    """Define 3-layer ML/MM system based on distance from ML region."""

    # Validate input
    if model_pdb is None and model_indices_str is None:
        click.echo(
            "ERROR: Either --model-pdb or --model-indices must be provided.; "
            'recover: pass --model-pdb /path/to/ml_region.pdb OR --model-indices "1,5,8-12"',
            err=True,
        )
        sys.exit(1)

    # Validate radii
    if radius_partial_hessian < 0:
        click.echo(
            "ERROR: --radius-partial-hessian must be non-negative.; "
            "recover: use a value >= 0 (e.g., --radius-partial-hessian 5.0). Pass 0 to disable.",
            err=True,
        )
        sys.exit(1)
    if radius_freeze < 0:
        click.echo(
            "ERROR: --radius-freeze must be non-negative.; "
            "recover: use a value >= 0 (e.g., --radius-freeze 10.0). Pass 0 to disable.",
            err=True,
        )
        sys.exit(1)

    # Parse model indices if provided
    model_indices = None
    if model_indices_str is not None:
        try:
            model_indices = _parse_indices_string(model_indices_str, one_based)
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            sys.exit(1)

    # Default output path
    if output_pdb is None:
        output_pdb = input_pdb.parent / f"{input_pdb.stem}_layered.pdb"
    output_pdb = _effective_output_pdb(output_pdb)

    # Echo configuration
    emit("\n====== Define Layer Configuration ======\n", narrative=True)
    click.echo(f"Input PDB: {input_pdb}")
    if model_pdb is not None:
        click.echo(f"Model PDB: {model_pdb}")
    if model_indices is not None:
        click.echo(f"Model indices: {len(model_indices)} atoms")
    click.echo(f"Radius (partial Hessian): {radius_partial_hessian} Å")
    click.echo(f"Radius (freeze): {radius_freeze} Å")
    if abs(radius_partial_hessian) > 1.0e-12:
        click.echo(
            "[warn] --radius-partial-hessian is ignored in 3-layer mode.",
            err=True,
        )
    click.echo(f"Output PDB: {output_pdb}")
    click.echo()

    try:
        layer_indices = define_layers(
            input_pdb=input_pdb,
            output_pdb=output_pdb,
            model_pdb=model_pdb,
            model_indices=model_indices,
            radius_partial_hessian=radius_partial_hessian,
            radius_freeze=radius_freeze,
        )

        # Print summary
        emit("\n====== Layer Summary ======\n", narrative=True)
        click.echo(f"Layer 1 (ML, B={BFACTOR_ML:.0f}):         {len(layer_indices['ml_indices']):6d} atoms")
        click.echo(f"Layer 2 (Movable MM, B={BFACTOR_MOVABLE_MM:.0f}): {len(layer_indices['movable_mm_indices']):6d} atoms")
        click.echo(f"Layer 3 (Frozen MM, B={BFACTOR_FROZEN:.0f}):     {len(layer_indices['frozen_indices']):6d} atoms")
        click.echo()
        total = (
            len(layer_indices['ml_indices']) +
            len(layer_indices['movable_mm_indices']) +
            len(layer_indices['frozen_indices'])
        )
        click.echo(f"Total atoms: {total}")
        click.echo()
        click.echo(f"[output] Wrote '{output_pdb}'")
        output_cif = output_pdb.with_suffix(".cif")
        if output_cif.exists():
            click.echo(f"[output] Wrote '{output_cif}'")

    except ValueError as e:
        click.echo(f"ERROR: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"ERROR: Unexpected error: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    cli()
