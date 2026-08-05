"""Product-local charge/spin preparation service (above ``core``).

These helpers derive the selected ML model's net charge from ``--ligand-charge`` metadata
(via ``extract``'s charge-summary engine) and resolve the final charge/spin for
a run.  They sit ABOVE ``core`` precisely because they consume a workflow-level
service (``mlmm.workflows.extract.compute_charge_summary``); keeping them out of
``mlmm.core.utils`` keeps the ``core.utils <-> extract`` dependency acyclic.
``core`` must never import a workflow, so this preparation
step lives here and workflow subcommands import it from this module.
"""

from __future__ import annotations

import math
from numbers import Integral
from pathlib import Path
from typing import TYPE_CHECKING, Mapping, Optional, Tuple

import click

from mlmm.core.defaults import BFACTOR_TOLERANCE

if TYPE_CHECKING:  # annotation only — no runtime import edge into ``core``
    from mlmm.core.utils import PreparedInputStructure

__all__ = [
    "configured_model_charge_spin",
    "infer_present_terminal_cap_ids",
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


def infer_present_terminal_cap_ids(
    structure,
    selected_ids: set[tuple],
) -> tuple[set[tuple], set[tuple]]:
    """Return selected amino acids with explicit NH3+ or COO- terminal atoms."""

    from mlmm.workflows.extract import AMINO_ACIDS

    keep_ncap_ids: set[tuple] = set()
    keep_ccap_ids: set[tuple] = set()
    for res in structure.get_residues():
        fid = res.get_full_id()
        if fid not in selected_ids or res.get_resname().upper() not in AMINO_ACIDS:
            continue
        atom_names = {atom.get_name().strip().upper() for atom in res.get_atoms()}
        if {"H1", "H2", "H3"} <= atom_names:
            keep_ncap_ids.add(fid)
        if "OXT" in atom_names:
            keep_ccap_ids.add(fid)
    return keep_ncap_ids, keep_ccap_ids


def _derive_charge_from_ligand_charge(
    pdb_path: Path,
    ligand_charge: Optional[str],
    *,
    select_bfactor_layer: bool,
    prefix: str = "",
) -> Optional[int]:
    """Derive the selected ML model's net charge from ``--ligand-charge`` metadata.

    Returns ``None`` when *ligand_charge* is ``None`` or derivation fails.
    """
    if ligand_charge is None:
        return None
    from Bio import PDB as BioPDB
    # Lazy: keep the heavy ``extract`` module out of the import chain until a
    # ligand-charge derivation actually runs (workflow -> workflow edge).
    from mlmm.workflows.extract import (
        compute_charge_summary,
        log_charge_summary,
    )

    parser = BioPDB.PDBParser(QUIET=True)
    complex_struct = parser.get_structure("complex", str(pdb_path))
    residues = list(complex_struct.get_residues())
    all_residue_ids = {res.get_full_id() for res in residues}
    if select_bfactor_layer:
        for res in residues:
            selected = [
                abs(atom.get_bfactor()) <= BFACTOR_TOLERANCE
                for atom in res.get_atoms()
            ]
            if any(selected) and not all(selected):
                raise click.ClickException(
                    f"{prefix} B-factor ML selection splits residue "
                    f"{res.get_resname()} {res.id[1]}; provide -q/--charge "
                    "for this atom-level model selection."
                )
        # Use the same predicate the layer parser uses (BFACTOR_TOLERANCE), so the derived
        # ML-region charge covers exactly the residues that were layered as ML. A hard-coded
        # 0.5 window dropped residues whose ML atoms sit inside the documented +-1.0 tolerance,
        # silently under-counting the charge.
        selected_ids = {
            res.get_full_id()
            for res in residues
            if any(abs(atom.get_bfactor()) <= BFACTOR_TOLERANCE for atom in res.get_atoms())
        }
    else:
        selected_ids = set(all_residue_ids)
    if not selected_ids:
        raise click.ClickException(
            f"{prefix} No residues belong to the selected ML model."
        )

    keep_ncap_ids, keep_ccap_ids = infer_present_terminal_cap_ids(
        complex_struct,
        selected_ids,
    )

    summary = compute_charge_summary(
        complex_struct,
        selected_ids,
        set(),
        ligand_charge,
        keep_ncap_ids=keep_ncap_ids,
        keep_ccap_ids=keep_ccap_ids,
    )
    log_charge_summary(prefix, summary)
    q_total = float(summary.get("total_charge", 0.0))
    click.echo(f"{prefix} Charge summary for the selected ML model:")
    click.echo(
        f"  Protein: {summary.get('protein_charge', 0.0):+g},  "
        f"Ligand: {summary.get('ligand_total_charge', 0.0):+g},  "
        f"Ions: {summary.get('ion_total_charge', 0.0):+g},  "
        f"Total: {q_total:+g}"
    )
    return _round_charge_with_note(q_total, prefix)


def configured_model_charge_spin(
    yaml_cfg: Optional[Mapping[str, object]],
) -> Tuple[Optional[int], Optional[int]]:
    """Read canonical/legacy YAML electronic state without package defaults."""

    calc_yaml: Mapping[str, object] = {}
    if isinstance(yaml_cfg, Mapping):
        from mlmm.cli.decorators import canonicalize_calculator_section

        normalized = canonicalize_calculator_section(dict(yaml_cfg))
        raw_calc = normalized.get("calc")
        if isinstance(raw_calc, Mapping):
            calc_yaml = raw_calc

    def _configured_int(key: str) -> Optional[int]:
        raw = calc_yaml.get(key)
        section = "calc"
        if raw is None:
            return None
        if isinstance(raw, bool) or not isinstance(raw, Integral):
            raise click.BadParameter(
                f"{section}.{key} must be an "
                f"integer, got {raw!r}."
            )
        return int(raw)

    charge = _configured_int("model_charge")
    spin = _configured_int("model_mult")
    if spin is not None and spin < 1:
        raise click.BadParameter(
            f"ML-region multiplicity must be an integer >= 1, got {spin}."
        )
    return charge, spin


def _configured_model_pdb(
    yaml_cfg: Optional[Mapping[str, object]],
) -> Optional[Path]:
    """Return the canonical/legacy YAML model-PDB source, if configured."""

    if not isinstance(yaml_cfg, Mapping):
        return None
    from mlmm.cli.decorators import canonicalize_calculator_section

    normalized = canonicalize_calculator_section(dict(yaml_cfg))
    calc_yaml = normalized.get("calc")
    canonical = (
        calc_yaml.get("model_pdb")
        if isinstance(calc_yaml, Mapping)
        else None
    )
    raw = canonical
    return None if raw in (None, "") else Path(str(raw))


def resolve_charge_spin_or_raise(
    prepared: "PreparedInputStructure",
    charge: Optional[int],
    spin: Optional[int],
    *,
    spin_default: int = 1,
    charge_default: Optional[int] = None,
    ligand_charge: Optional[str] = None,
    model_pdb: Optional[Path] = None,
    model_indices_spec: Optional[str] = None,
    detect_layer: bool = True,
    yaml_cfg: Optional[Mapping[str, object]] = None,
    prefix: str = "",
) -> Tuple[int, int]:
    """Resolve charge/spin from inputs.

    Priority: explicit ``-q/-m`` > ``--ligand-charge`` derivation > configured
    ``calc.model_charge/model_mult`` (legacy ``mlmm`` section) > defaults.
    Raises :class:`click.ClickException` when charge cannot be resolved.
    """
    if charge is None and ligand_charge is not None:
        # The charge must be derived from the same effective ML-region source
        # that the calculator will use.  Every workflow passes the merged YAML
        # here, so resolve a YAML-only model_pdb once at this shared boundary.
        effective_model_pdb = (
            Path(model_pdb)
            if model_pdb is not None
            else _configured_model_pdb(yaml_cfg)
        )
        if model_indices_spec:
            raise click.ClickException(
                f"{prefix} Automatic charge derivation is unavailable with "
                "--model-indices; provide -q/--charge for the exact atom selection."
            )
        if effective_model_pdb is None and not detect_layer:
            raise click.ClickException(
                f"{prefix} Automatic charge derivation requires --detect-layer "
                "or --model-pdb; otherwise provide -q/--charge."
            )
        charge_source = (
            effective_model_pdb
            if effective_model_pdb is not None
            else prepared.source_path
        )
        try:
            charge = _derive_charge_from_ligand_charge(
                charge_source,
                ligand_charge,
                select_bfactor_layer=effective_model_pdb is None,
                prefix=prefix,
            )
        except click.ClickException:
            raise
        except Exception as exc:
            raise click.ClickException(
                f"{prefix} Failed to derive the selected ML-model charge: {exc}"
            ) from exc
    configured_charge, configured_spin = configured_model_charge_spin(yaml_cfg)
    if charge is None:
        charge = configured_charge
    if charge is None:
        if charge_default is None:
            raise click.ClickException(
                "ML-region charge is unresolved. Provide -q/--charge or "
                "--ligand-charge."
            )
        charge = charge_default
    if spin is None:
        spin = configured_spin
    if spin is None:
        spin = spin_default
    if int(spin) < 1:
        raise click.BadParameter(
            f"ML-region multiplicity must be an integer >= 1, got {spin}."
        )
    return int(charge), int(spin)
