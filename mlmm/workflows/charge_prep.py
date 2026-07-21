"""Product-local charge/spin preparation service (above ``core``).

These helpers derive a system's total charge from ``--ligand-charge`` metadata
(via ``extract``'s charge-summary engine) and resolve the final charge/spin for
a run.  They sit ABOVE ``core`` precisely because they consume a workflow-level
service (``mlmm.workflows.extract.compute_charge_summary``); keeping them out of
``mlmm.core.utils`` is what breaks the historical ``core.utils <-> extract``
import cycle (M39).  ``core`` must never import a workflow, so this preparation
step lives here and workflow subcommands import it from this module.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Tuple

import click

if TYPE_CHECKING:  # annotation only — no runtime import edge into ``core``
    from mlmm.core.utils import PreparedInputStructure

__all__ = [
    "resolve_charge_spin_or_raise",
]


def _round_charge_with_note(q: float, prefix: str = "") -> int:
    """Round a float charge to the nearest integer, with a note if not exact."""
    if not math.isfinite(q):
        raise click.BadParameter(f"Computed total charge is non-finite: {q!r}")
    q_int = int(round(q))
    if abs(float(q) - q_int) > 1e-6:
        click.echo(
            f"{prefix} NOTE: total charge = {q:+g} → rounded to integer {q_int:+d}."
        )
    return q_int


def _derive_charge_from_ligand_charge(
    pdb_path: Path,
    ligand_charge: Optional[str],
    *,
    prefix: str = "",
) -> Optional[int]:
    """Derive total system charge from a PDB file using ``--ligand-charge`` metadata.

    Returns ``None`` when *ligand_charge* is ``None`` or derivation fails.
    """
    if ligand_charge is None:
        return None
    try:
        from Bio import PDB as BioPDB
        # Lazy: keep the heavy ``extract`` module out of the import chain until a
        # ligand-charge derivation actually runs (workflow -> workflow edge).
        from mlmm.workflows.extract import compute_charge_summary, log_charge_summary

        parser = BioPDB.PDBParser(QUIET=True)
        complex_struct = parser.get_structure("complex", str(pdb_path))

        # Use only ML-region residues (B-factor ≈ 0) when layered PDB is available.
        # A residue is included if ANY of its atoms has B-factor < 1.0 (ML layer).
        ml_residue_ids = set()
        all_residue_ids = set()
        for res in complex_struct.get_residues():
            fid = res.get_full_id()
            all_residue_ids.add(fid)
            for atom in res.get_atoms():
                if atom.get_bfactor() < 1.0:
                    ml_residue_ids.add(fid)
                    break
        # Fall back to all residues if no B-factor layering is present
        # (i.e. every residue has B=0 means unlayered PDB).
        selected_ids = ml_residue_ids or all_residue_ids
        summary = compute_charge_summary(
            complex_struct, selected_ids, set(), ligand_charge
        )
        log_charge_summary(prefix, summary)
        q_total = float(summary.get("total_charge", 0.0))
        click.echo(
            f"{prefix} Charge summary (--ligand-charge):"
        )
        click.echo(
            f"  Protein: {summary.get('protein_charge', 0.0):+g},  "
            f"Ligand: {summary.get('ligand_total_charge', 0.0):+g},  "
            f"Ions: {summary.get('ion_total_charge', 0.0):+g},  "
            f"Total: {q_total:+g}"
        )
        return _round_charge_with_note(q_total, prefix)
    except Exception as e:
        click.echo(
            f"{prefix} NOTE: failed to derive charge from --ligand-charge: {e}",
            err=True,
        )
        return None


def resolve_charge_spin_or_raise(
    prepared: "PreparedInputStructure",
    charge: Optional[int],
    spin: Optional[int],
    *,
    spin_default: int = 1,
    charge_default: Optional[int] = None,
    ligand_charge: Optional[str] = None,
    prefix: str = "",
) -> Tuple[int, int]:
    """Resolve charge/spin from inputs.

    Priority: explicit ``-q/--charge`` > ``--ligand-charge`` derivation >
    ``charge_default``.  Raises :class:`click.ClickException` when charge
    cannot be resolved.
    """
    if charge is None and ligand_charge is not None:
        charge = _derive_charge_from_ligand_charge(
            prepared.source_path, ligand_charge, prefix=prefix,
        )
    if charge is None:
        if charge_default is None:
            raise click.ClickException(
                "Total charge is unresolved. Provide -q/--charge or --ligand-charge."
            )
        charge = charge_default
    if spin is None:
        spin = spin_default
    return int(charge), int(spin)
