# mlmm/extract.py

"""
extract — Automated binding‑pocket (active‑site) extractor
====================================================================

Usage (CLI)
-----------
    mlmm extract -i INPUT.pdb [INPUT2.pdb ...] -c <substrate_spec> \
        [-o OUTPUT.pdb ...] [-r <Å>] [--radius-het2het <Å>] \
        [--include-h2o/--no-include-h2o] [--exclude-backbone/--no-exclude-backbone] \
        [--add-linkh/--no-add-linkh] [--selected-resn "CHAIN:RES" ...] \
        [-l, --ligand-charge <number|"RES:Q,...">] [--verbose/--no-verbose]

Examples
--------
    # Minimal (ID-based substrate) with explicit total ligand charge
    mlmm extract -i complex.pdb -c A:123 -o pocket.pdb -l -3

    # Substrate provided as a PDB; per-resname charge mapping (others remain 0)
    mlmm extract -i complex.pdb -c substrate.pdb -o pocket.pdb \
        -l "GPP:-3,MMT:-1"

    # Name-based substrate selection including all matches (WARNING is logged)
    mlmm extract -i complex.pdb -c "GPP,MMT" -o pocket.pdb -l -4

    # Multi-structure to single multi-MODEL output with hetero-hetero proximity enabled
    mlmm extract -i complex1.pdb complex2.pdb -c A:123 \
        -o pocket_multi.pdb --radius-het2het 2.6 -l -3 --verbose

Description
-----------
Extracts an active‑site pocket around specified substrate residues from a protein–substrate complex,
applies biochemically aware truncation (backbone/side‑chain capping with safeguards), and can append
link hydrogens for cut bonds. Supports single structures and ensembles (multi‑MODEL or per‑file outputs).
Typical use cases include QM/MM, ML/MM, and cluster QM models.

Residue inclusion
-----
- Always include the substrate residues.
- Standard cutoff (``--radius``, default 2.6 Å):
  - If ``--no-exclude-backbone`` (default): include any residue if **any atom** is within the cutoff.
  - If ``--exclude-backbone``: for **amino‑acid residues**, the qualifying atom
    must be **non‑backbone** (not in {N, H*, CA, HA*, C, O, OXT}); non‑amino‑acid residues qualify by any atom.
- Independent hetero–hetero proximity (``--radius-het2het``):
  add residues if a **substrate hetero atom (non‑C/H)** is within the cutoff of a **protein hetero atom**.
  With ``--exclude-backbone``, amino‑acid neighbors must be **non‑backbone** atoms.
- Waters are included by default (``--include-h2o``; disable with ``--no-include-h2o``).
- ``--selected-resn`` force‑includes residues (chain and insertion codes supported).
- When ``--no-exclude-backbone`` and a selected residue’s **backbone atom** contacts the substrate
  (within either cutoff), include its peptide‑adjacent N‑side and C‑side neighbors (C–N ≤ 1.9 Å). For true termini,
  keep the respective terminal cap atoms (N/H* or C/O/OXT).
- **Disulfide safeguard:** if a selected CYS/CYX forms an SG–SG contact ≤ 2.5 Å, include both partners.
- **Proline safeguard:** if a selected **PRO** is not N‑terminal (peptide‑adjacent), include the immediately
  preceding amino acid. For that neighbor, **CA is always kept**, and when backbone exclusion is on,
  keep **C** and **O/OXT** (to preserve the bond into PRO–N).

Truncation (capping)
-----
- **Isolated residues** → keep **pure side‑chain** (remove N, CA, C, O, OXT and N/CA H*).
  - **PRO/HYP** retain N, CA, HA, H* to keep the ring.
- **Continuous peptide stretches** keep internal backbone; only **terminal caps** are removed
  (N‑cap: N/H*; C‑cap: C/O/OXT). TER‑aware segmentation prevents crossing chain breaks.
- With ``--exclude-backbone``, delete main‑chain atoms on all **non‑substrate amino acids**,
  except for the specific PRO/HYP retention and PRO‑adjacency preservation above.
- **Non‑amino‑acid residues**: atoms named like protein backbone ({"N","CA","HA","H","H1","H2","H3"})
  are **never deleted** by capping logic.

Link hydrogens (--add-linkh)
-----
- Adds **carbon‑only** link H at **1.09 Å** along the cut‑bond vector.
  Normal residues: checks **CB–CA**, **CA–N**, **CA–C**; **PRO/HYP**: **CA–C** only.
- If any are added, append a **TER** then a contiguous **HETATM** block with atoms named **``HL``**
  in residue **``LKH``** (chain ``L``); serials continue from the main block.
- In multi‑structure mode, link‑H **targets and ordering** are enforced to be **identical across models**;
  coordinates remain model‑specific. The flag is honored in both single‑ and multi‑structure modes.

Charge summary
-----
- **AMINO_ACIDS** dictionary supplies nominal integer charges for amino‑acid residues (common variants included).
- **ION** supplies charges for common ions (e.g., ZN, MG, FE2). Waters are 0.
- **Unknown residues** (not in AMINO_ACIDS/ION/WATER) are **0** unless ``--ligand-charge`` is given.
  - ``--ligand-charge <number>``: total charge distributed across **unknown substrate** residues
    (or across all unknowns if no unknown substrate).
  - ``--ligand-charge "RES1:Q1,RES2:Q2"``: set **per‑resname** charges; any other unknown residues remain 0.
- In multi‑structure mode the **charge summary is computed on the first input PDB** only.

Multi-structure ensembles
-----
- Accepts multiple input PDBs (same **atom count**; ordering is **assumed identical** and is
  **spot‑checked** at the beginning and end of the atom list).
- Each structure is selected independently; the **union** of selected residues is applied to all.
  Disulfides, PRO‑adjacency, and (if enabled) backbone‑contact neighbor augmentation are also unioned.
- Outputs:
  - Provide **one** output path → **multi‑MODEL** PDB (one MODEL per input).
  - Provide **N** output paths where **N == number of inputs** → **N** single‑model PDBs.
  - If ``-o`` is omitted with multiple inputs → per‑file outputs ``pocket_{original_filename}.pdb``.
- **Diagnostics:** atom counts (**raw** vs **after truncation**) are logged **per model**.

Substrate specification
-----
``-c/--center`` accepts:
- a **PDB path** (exact coordinate match on the first input; IDs propagated to others),
- a list of **residue IDs**: ``"123,124"``, ``"A:123,B:456"``, ``"123A"``, ``"A:123A"`` (insertion codes OK),
- or a list of **residue names** (case‑insensitive), e.g., ``"GPP,MMT"``.
  If multiple residues share the same name, **all** matches are used and a **WARNING** is logged.

Outputs (& Directory Layout)
----------------------------
<output>/ (default: pocket.pdb for single input; pocket_<source>.pdb per input when -o is omitted)
  ├─ pocket.pdb                         # default single-input pocket
  ├─ pocket_<original_filename>.pdb     # default per-input pocket when multiple inputs and -o omitted
  └─ <user_paths>.pdb                   # custom outputs; one path = multi-MODEL, N paths = per-structure

Link hydrogens, logs, and programmatic use
  ├─ Link-H block (when added) follows a TER as contiguous HETATM records named HL in residue LKH (chain L).
  ├─ INFO logs summarize residue selection, raw/kept atom counts, and the charge summary.
  └─ ``extract(..., api=True)`` / ``extract_api(...)`` returns ``{"outputs": [...], "counts": [...], "charge_summary": {...}}``.


Notes
-----
- **Defaults / behavior:**
  - ``--radius`` default: **2.6 Å**. If given **0**, internally nudged to **0.001 Å**.
  - ``--radius-het2het`` default: **0 Å** (off). Internally treated as **0.001 Å** if ``0`` is given.
  - ``--include-h2o`` default: **true**.
  - ``--exclude-backbone`` default: **false**.
  - ``--add-linkh`` default: **false**.
  - ``--ligand-charge`` default: **None** (unknown residues counted as 0 unless set).
  - Output default: single input → ``pocket.pdb``; multiple inputs → ``pocket_{original_filename}.pdb``.
- **Geometry thresholds and tolerances:**
  - Peptide adjacency: **C(prev)–N(next) ≤ 1.9 Å** (distance‑based; practical TER awareness).
  - Disulfide detection: **SG–SG ≤ 2.5 Å**.
  - Link‑H distance: **1.09 Å** (C–H) along the cut‑bond vector.
  - Exact match tolerance for substrate PDB: **1e‑3 Å** per atom.
- **Safeguards and special cases:**
  - **PRO/HYP** retain N, CA, HA, H* in isolated truncations; PRO’s **N‑side neighbor** is auto‑included
    when peptide‑adjacent; **CA** on that neighbor is always kept, and with backbone exclusion
    **C** and **O/OXT** are preserved to maintain the peptide bond into PRO–N.
  - **Non‑amino‑acid residues** never lose atoms named like backbone (``N, CA, HA, H, H1, H2, H3``).
  - **Waters** (HOH/WAT/H2O/DOD/TIP/TIP3/SOL) are always neutral (charge 0) and included by default.
- **Dependencies:** Python ≥ **3.10** (PEP 604 unions), Biopython ≥ **1.80**, NumPy.
- **Python API (for reference):**
  - ``extract(args: argparse.Namespace | None = None, api=False)`` — main entry (CLI or programmatic).
  - ``extract_api(...)`` — convenience wrapper that returns ``{'outputs','counts','charge_summary'}`` when used programmatically.
"""

from __future__ import annotations

import argparse
import logging
import io as _io
import os
import re
import sys
from typing import Dict, List, Set, Tuple, Iterable, Any, Optional, Sequence

import click
import numpy as np
from Bio import PDB
from Bio.PDB import NeighborSearch

# Public API
__all__ = ["extract", "extract_api"]


LOGGER = logging.getLogger(__name__)


class _ClickEchoHandler(logging.Handler):
    """Route logger records to click.echo with stderr for warnings/errors."""

    def emit(self, record: logging.LogRecord) -> None:
        try:
            click.echo(self.format(record), err=record.levelno >= logging.WARNING)
        except Exception:
            self.handleError(record)


def _configure_extract_logger(verbose: bool) -> None:
    """Configure module-local logger for CLI output without touching global logging."""
    level = logging.INFO if verbose else logging.WARNING
    LOGGER.setLevel(level)
    LOGGER.propagate = False

    for handler in list(LOGGER.handlers):
        LOGGER.removeHandler(handler)

    handler = _ClickEchoHandler()
    handler.setLevel(level)
    handler.setFormatter(
        logging.Formatter("%(message)s" if verbose else "%(levelname)s: %(message)s")
    )
    LOGGER.addHandler(handler)


# ---------------------------------------------------------------------
#   Constants
# ---------------------------------------------------------------------
BACKBONE_ATOMS: Set[str] = {
    "N", "C", "O", "CA", "OXT",
    "H", "H1", "H2", "H3", "HN", "HA", "HA2", "HA3",
}
# When --exclude-backbone true, remove the full main-chain set:
BACKBONE_ALL: Set[str] = BACKBONE_ATOMS

# Unified amino-acid dictionary: resname -> nominal integer charge
# (membership checks throughout the code use dictionary keys)
AMINO_ACIDS: Dict[str, int] = {
    # --- Standard 20 (L) ---
    "ALA":  0, "ARG": +1, "ASN":  0, "ASP": -1, "CYS":  0,
    "GLU": -1, "GLN":  0, "GLY":  0, "HIS":  0, "ILE":  0,
    "LEU":  0, "LYS": +1, "MET":  0, "PHE":  0, "PRO":  0,
    "SER":  0, "THR":  0, "TRP":  0, "TYR":  0, "VAL":  0,

    # --- Canonical extras ---
    "SEC":  0,   # selenocysteine
    "PYL": +1,   # pyrrolysine

    # --- Protonation / tautomers (Amber/CHARMM style) ---
    "HIP": +1,   # fully protonated His
    "HID":  0,   # Nδ-protonated His
    "HIE":  0,   # Nε-protonated His
    "ASH":  0,   # neutral Asp
    "GLH":  0,   # neutral Glu
    "LYN":  0,   # neutral Lys
    "ARN":  0,   # neutral Arg
    "TYM": -1,   # deprotonated Tyr (phenolate)

    # --- Phosphorylated residues ---
    "SEP": -2, "TPO": -2, "PTR": -2,
    "S1P": -1, "T1P": -1, "Y1P": -1,   # monoanionic phospho-Ser/Thr/Tyr

    # --- Phosphorylated histidines (phosaa19SB) ---
    "H1D":  0,  # ND1-phospho-His, neutral
    "H2D": -1,  # ND1-phospho-His, anionic
    "H1E":  0,  # NE2-phospho-His, neutral
    "H2E": -1,  # NE2-phospho-His, anionic
    
    # --- Cys family ---
    "CYX":  0,   # disulfide Cys
    "CSO":  0,   # Cys sulfenic acid
    "CSD": -1,   # Cys sulfinic acid
    "CSX":  0,   # generic Cys derivative
    "OCS": -1,   # cysteic acid
    "CYM": -1,   # deprotonated Cys

    # --- Lys variants / carboxylation ---
    "MLY": +1, "LLP": +1, "DLY": +1,
    "KCX": -1,   # Lysine Nz-Carboxylic Acid

    # --- D isomers (19 residues) ---
    "DAL":  0, "DAR": +1, "DSG": 0, "DAS": -1, "DCY": 0,
    "DGN":  0, "DGL": -1, "DHI": 0, "DIL":  0, "DLE": 0,
    "DLY": +1, "MED":  0, "DPN": 0, "DPR":  0, "DSN": 0,
    "DTH":  0, "DTR":  0, "DTY": 0, "DVA":  0,

    # --- Carboxylation / cyclization / others ---
    "CGU": -2,   # gamma-carboxy-glutamate
    "CGA": -1,   # carboxymethylated glutamate
    "PCA":  0,   # pyroglutamate
    "MSE":  0,   # selenomethionine
    "OMT":  0,   # methionine sulfone

    # --- Other modified residues possibly encountered ---
    "ASA": 0, "CIR": 0, "FOR": 0, "MVA": 0, "IIL": 0, "AIB": 0, "HTN": 0,
    "SAR": 0, "NMC": 0, "PFF": 0, "NFA": 0, "ALY": 0, "AZF": 0, "CNX": 0, "CYF": 0,

    # --- Hydroxyproline ---
    "HYP": 0,

    # --- All C-terminus ---
    "CALA": -1, "CARG":  0, "CASN": -1, "CASP": -2, "CCYS": -1,
    "CCYX": -1, "CGLN": -1, "CGLU": -2, "CGLY": -1, "CHID": -1,
    "CHIE": -1, "CHIP":  0, "CHYP": -1, "CILE": -1, "CLEU": -1,
    "CLYS":  0, "CMET": -1, "CPHE": -1, "CPRO": -1, "CSER": -1,
    "CTHR": -1, "CTRP": -1, "CTYR": -1, "CVAL": -1, "NHE": 0,
    "NME": 0,
    "CTER": -1,  # generic C-terminus
    
    # --- All N-terminus ---
    "NALA": +1, "NARG": +2, "NASN": +1, "NASP":  0, "NCYS": +1,
    "NCYX": +1, "NGLN": +1, "NGLU":  0, "NGLY": +1, "NHID": +1,
    "NHIE": +1, "NHIP": +2, "NILE": +1, "NLEU": +1, "NLYS": +2,
    "NMET": +1, "NPHE": +1, "NPRO": +1, "NSER": +1, "NTHR": +1,
    "NTRP": +1, "NTYR": +1, "NVAL": +1, "ACE": 0,
    "NTER": +1,  # generic N-terminus
}

# Common ions (by residue name) and their formal charges
ION: Dict[str, int] = {
    # +1
    "LI": +1, "NA": +1, "K": +1, "RB": +1, "CS": +1, "TL": +1, "AG": +1, "CU1": +1,
    "Ag": +1, "K+": +1, "Na+": +1, "NH4": +1, "H3O+": +1, "HE+": +1, "HZ+": +1, "Tl": +1,

    # +2
    "MG": +2, "CA": +2, "SR": +2, "BA": +2, "MN": +2, "FE2": +2, "CO": +2, "NI": +2,
    "CU": +2, "ZN": +2, "CD": +2, "HG": +2, "PB": +2, "Be": +2, "PD": +2, "PT": +2,
    "Sn": +2, "Ra": +2, "YB2": +2, "V2+": +2,

    # +3
    "FE": +3, "AU3": +3, "AL": +3, "GA": +3, "IN": +3,
    "CE": +3, "Ce": +3, "CR": +3, "Cr": +3, "Dy": +3, "EU": +3, "EU3": +3, "Er": +3,
    "GD3": +3, "LA": +3, "LU": +3, "Nd": +3, "PR": +3, "SM": +3, "Sm": +3, "TB": +3,
    "Tm": +3, "Y": +3, "Pu": +3,

    # +4
    "U4+": +4, "Th": +4, "Hf": +4, "Zr": +4,

    # -1
    "F": -1, "CL": -1, "BR": -1, "I": -1, "Cl-": -1, "IOD": -1,
}

DISULFIDE_CUTOFF = 2.5   # Å Sγ–Sγ (SG–SG)
EXACT_EPS = 1e-3         # Å tolerance for exact match
WATER_RES = {"HOH","WAT","H2O","DOD","TIP","TIP3","SOL"}

# Type for cross-structure residue identity (chain, hetflag, resseq, icode, resname)
ResidueKey = Tuple[str, str, int, str, str]

# ---------------------------------------------------------------------
#   Helpers
# ---------------------------------------------------------------------


def _extract_short_help() -> str:
    return "\n".join(
        [
            "Usage: mlmm extract [OPTIONS]",
            "",
            "Extract a binding pocket around substrate residues.",
            "",
            "Core options:",
            "  -i, --input PATH [PATH ...]      Input complex PDB file(s).",
            "  -c, --center SPEC                Substrate selector (PDB / residue IDs / residue names).",
            "  -o, --output PATH [PATH ...]     Output pocket PDB path(s).",
            "  -r, --radius FLOAT               Pocket inclusion radius in angstrom.",
            "  -l, --ligand-charge VALUE        Total or mapped ligand charge.",
            "  --help-advanced                  Show full extract options and exit.",
            "",
            "Use '--help-advanced' to see all extractor options.",
        ]
    )


# ── Native Click command ─────────────────────────────────────────────

_EXTRACT_ALL_FLAGS = (
    "-i", "--input",
    "-c", "--center",
    "-o", "--output",
    "-r", "--radius",
    "--radius-het2het",
    "--include-h2o", "--no-include-h2o",
    "--exclude-backbone", "--no-exclude-backbone",
    "--add-linkh", "--no-add-linkh",
    "--selected-resn",
    "-l",
    "--ligand-charge",
    "-v", "--verbose", "--no-verbose",
    "-h", "--help", "--help-advanced",
)


def _gather_extract_variadic(
    ctx_args: List[str],
    flag_names: Sequence[str],
) -> List[str]:
    """Collect variadic positional values after *flag_names* from Click extra args."""
    names_set = set(flag_names)
    stop_set = set(_EXTRACT_ALL_FLAGS)
    vals: List[str] = []
    i = 0
    while i < len(ctx_args):
        tok = ctx_args[i]
        if tok in names_set:
            j = i + 1
            while j < len(ctx_args) and ctx_args[j] not in stop_set:
                vals.append(ctx_args[j])
                j += 1
            i = j
        else:
            i += 1
    return vals


@click.command(
    name="extract",
    help=(
        "Extract a binding pocket around substrate residues (from a PDB or "
        "residue IDs/names), with biochemically aware truncation and optional "
        "link-H; supports multi-structure input and multi-MODEL output."
    ),
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
@click.option(
    "-i", "--input", "complex_pdb",
    type=str, multiple=True, required=True,
    help="Protein-substrate complex PDB(s). If multiple, they must have identical atom counts and ordering.",
)
@click.option(
    "-c", "--center", "substrate_pdb",
    type=str, required=True,
    help=(
        "Substrate specification: a PDB path, a comma/space-separated residue-ID list "
        "like '123,124' or 'A:123,B:456' (insertion codes supported), "
        "or a residue-name list like 'GPP,SAM'."
    ),
)
@click.option(
    "-o", "--output", "output_pdb",
    type=str, multiple=True, default=(),
    help=(
        "Output PDB path(s). One path for multi-MODEL PDB, or N paths for per-file output. "
        "If omitted: single input -> pocket.pdb; multiple inputs -> pocket_{filename}.pdb."
    ),
)
@click.option(
    "-r", "--radius",
    type=float, default=2.6, show_default=True,
    help="Cutoff (angstrom) around substrate atoms for pocket inclusion.",
)
@click.option(
    "--radius-het2het",
    type=float, default=0, show_default=True,
    help="Cutoff (angstrom) for substrate-protein hetero-atom proximity (non-C/H). 0 disables.",
)
@click.option(
    "--include-h2o/--no-include-h2o",
    "include_h2o",
    default=True, show_default=True,
    help="Include waters (HOH/WAT/H2O/DOD/TIP/TIP3/SOL).",
)
@click.option(
    "--exclude-backbone/--no-exclude-backbone",
    default=False, show_default=True,
    help="Delete main-chain atoms from non-substrate amino acids.",
)
@click.option(
    "--add-linkh/--no-add-linkh",
    "add_linkh",
    default=False, show_default=True,
    help="Add carbon-only link-H at 1.09 angstrom along cut-bond directions.",
)
@click.option(
    "--selected-resn",
    type=str, default="",
    help="Comma/space-separated residue IDs to force-include.",
)
@click.option(
    "--modified-residue",
    type=str, default="",
    help=(
        "Comma-separated residue names (with optional charge) to treat as amino acids "
        "for backbone truncation and charge assignment. "
        "Examples: 'HD1,HD2,HD3' (charge defaults to 0) or 'HD1:0,SEP:-2'."
    ),
)
@click.option(
    "-l",
    "--ligand-charge",
    type=str, default=None,
    help="Total charge number or per-resname mapping like 'GPP:-3,SAM:1'.",
)
@click.option(
    "-v", "--verbose/--no-verbose",
    default=True, show_default=True,
    help="Enable INFO-level logging.",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json next to the output PDB.",
)
@click.pass_context
def cli(
    ctx: click.Context,
    complex_pdb: Sequence[str],
    substrate_pdb: str,
    output_pdb: Sequence[str],
    radius: float,
    radius_het2het: float,
    include_h2o: bool,
    exclude_backbone: bool,
    add_linkh: bool,
    selected_resn: str,
    modified_residue: str,
    ligand_charge: Optional[str],
    verbose: bool,
    out_json: bool,
) -> None:
    # Recover variadic values after -i / -o from extra args (supports
    # space-separated syntax: ``-i a.pdb b.pdb`` in addition to ``-i a.pdb -i b.pdb``).
    extra_inputs = _gather_extract_variadic(ctx.args, ("-i", "--input"))
    input_list = list(complex_pdb) + extra_inputs if extra_inputs else list(complex_pdb)

    extra_outputs = _gather_extract_variadic(ctx.args, ("-o", "--output"))
    output_list: Optional[List[str]]
    if output_pdb or extra_outputs:
        output_list = list(output_pdb) + extra_outputs if extra_outputs else list(output_pdb)
    else:
        output_list = None

    ns = argparse.Namespace(
        complex_pdb=input_list,
        substrate_pdb=substrate_pdb,
        output_pdb=output_list,
        radius=radius,
        radius_het2het=radius_het2het,
        include_h2o=include_h2o,
        exclude_backbone=exclude_backbone,
        add_linkh=add_linkh,
        selected_resn=selected_resn,
        modified_residue=modified_residue,
        ligand_charge=ligand_charge,
        verbose=verbose,
    )
    result = extract(ns, api=out_json)

    if out_json and result is not None:
        from pathlib import Path as _Path
        from .utils import write_result_json

        first_output = (output_list or ["model.pdb"])[0]
        out_dir = _Path(first_output).resolve().parent

        counts = result.get("counts", [{}])
        first_counts = counts[0] if counts else {}
        charge_summary = result.get("charge_summary", {})

        result_data = {
            "status": "ok",
            "n_atoms_raw": first_counts.get("raw_atoms"),
            "n_atoms_extracted": first_counts.get("kept_atoms"),
            "total_charge": charge_summary.get("total_charge"),
            "protein_charge": charge_summary.get("protein_charge"),
            "ligand_total_charge": charge_summary.get("ligand_total_charge"),
            "ion_total_charge": charge_summary.get("ion_total_charge"),
            "ion_charges": charge_summary.get("ion_charges"),
            "unknown_residue_charges": charge_summary.get("unknown_residue_charges"),
            "n_link_hydrogens": result.get("n_link_hydrogens", 0),
            "files": {_Path(o).name: str(o) for o in result.get("outputs", [])},
            "center": substrate_pdb,
            "radius": radius,
            "exclude_backbone": exclude_backbone,
            "include_h2o": include_h2o,
            "ligand_charge_input": ligand_charge,
            "input_files": [str(p) for p in input_list],
        }
        write_result_json(out_dir, result_data, command="extract")


def _build_arg_parser(*, prog: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Extract a binding pocket around substrate residues (from a PDB or residue IDs/names), "
            "with biochemically aware truncation and optional link‑H; supports multi‑structure input "
            "and multi‑MODEL output. Also logs pocket charge summary."
        )
    )

    p.add_argument(
        "-i", "--input", dest="complex_pdb", required=True, nargs="+",
        metavar="complex.pdb",
        help="Protein–substrate complex PDB(s). If multiple, they must have identical atom counts and ordering."
    )
    p.add_argument(
        "-c", "--center", dest="substrate_pdb", required=True,
        metavar="substrate.pdb | '123,124' | 'A:123,B:456' | 'GPP,MMT'",
        help=("Substrate specification: either a PDB containing exactly the substrate residue(s), "
              "a comma/space‑separated residue‑ID list like '123,124' or 'A:123,B:456' "
              "(insertion codes supported: '123A' / 'A:123A'), "
              "or a comma/space‑separated **residue‑name** list like 'GPP,MMT'. "
              "When residue names are used and multiple residues share a name, all are used and a WARNING is logged.")
    )
    p.add_argument(
        "-o", "--output", dest="output_pdb", required=False, nargs="+",
        metavar="pocket.pdb", default=None,
        help=("Output PDB path(s). Provide one path to write a single multi‑MODEL PDB, "
              "or provide N paths where N == number of inputs to write N single‑model PDBs (one per input, in order). "
              "If omitted: single input → pocket.pdb; multiple inputs → pocket_{original_filename}.pdb.")
    )
    p.add_argument(
        "-r", "--radius", type=float, default=2.6,
        help=("Cutoff (Å) around substrate atoms. With --exclude-backbone, an **amino-acid** "
              "neighbor must have a **non-backbone** atom within this distance; otherwise **any atom** suffices. "
              "(default: 2.6)")
    )
    p.add_argument(
        "--radius-het2het", type=float, default=0,
        help=("Cutoff (Å) for substrate–protein hetero‑atom proximity (non‑C/H on both sides); "
              "applied independently of --radius. 0 conceptually disables this rule, "
              "but is internally treated as 0.001 Å. (default: 0)")
    )
    p.add_argument(
        "--include-h2o", "--include-H2O",
        dest="include_h2o",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include waters (HOH/WAT/H2O/DOD/TIP/TIP3/SOL). (default: True)"
    )
    p.add_argument(
        "--exclude-backbone",
        dest="exclude_backbone",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Delete main‑chain atoms (N, H*, CA, HA*, C, O, OXT) from non‑substrate amino acids; PRO/HYP keep N, CA, HA, H*. (default: False)"
    )
    p.add_argument(
        "--add-linkh", "--add-linkH",
        dest="add_linkh",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Add carbon‑only link‑H at 1.09 Å along cut‑bond directions; appended after a TER as HL/LKH HETATM records. (default: False)"
    )
    p.add_argument(
        "--selected-resn", dest="selected_resn", required=False, default="",
        help=("Comma/space‑separated residue IDs to force‑include (e.g., '123,124', 'A:123,B:456'; "
              "insertion codes allowed: '123A' / 'A:123A').")
    )
    p.add_argument(
        "-l", "--ligand-charge", type=str, default=None,
        help=("Either a single **number** giving the **total** charge to distribute across unknown residues "
              "(preferring unknown substrate), or a comma/space‑separated **per‑resname** list like "
              "'GPP:-3,MMT:-1'. In mapping mode, any other unknown residues remain 0.")
    )
    p.add_argument(
        "-v", "--verbose",
        dest="verbose",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=("Enable INFO-level logging."
              " default: True.")
    )
    return p


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """
    Parse CLI arguments.

    Parameters
    ----------
    argv : Optional[Sequence[str]]
        Command-line arguments to parse. If None, uses sys.argv.

    Returns
    -------
    argparse.Namespace
        Parameters for running the pocket extraction.
    """
    argv_list = list(argv) if argv is not None else None
    if argv_list is not None:
        wants_adv = "--help-advanced" in argv_list
        wants_help = ("--help" in argv_list) or ("-h" in argv_list)
        if wants_help and not wants_adv:
            click.echo(_extract_short_help())
            raise SystemExit(0)
        if wants_adv:
            argv_list = [a for a in argv_list if a != "--help-advanced"]
            if ("--help" not in argv_list) and ("-h" not in argv_list):
                argv_list.append("--help")

    p = _build_arg_parser(prog="mlmm extract")
    return p.parse_args(args=argv_list)


def load_structure(path: str, name: str) -> PDB.Structure.Structure:
    """
    Load a PDB file into a Biopython Structure object.
    """
    parser = PDB.PDBParser(QUIET=True)
    return parser.get_structure(name, path)


# ---------------------------------------------------------------------
#   Formatting helpers (for logging / API)
# ---------------------------------------------------------------------

def _fmt_res_id(res: PDB.Residue.Residue) -> str:
    """
    Return a compact residue tag like 'A:123A SER' or '123 SER'.
    """
    chain = res.get_parent().id or ""
    het, resseq, icode = res.id
    icode_txt = "" if icode == " " else icode
    chain_txt = f"{chain}:" if chain else ""
    return f"{chain_txt}{resseq}{icode_txt} {res.get_resname()}"


def _fmt_fid(structure, fid: Tuple) -> str:
    """
    Format a full-id into a human-friendly residue tag.
    """
    res: PDB.Residue.Residue = structure[fid[1]][fid[2]].child_dict[fid[3]]
    return _fmt_res_id(res)


# ---------------------------------------------------------------------
#   Substrate matching
# ---------------------------------------------------------------------

def is_exact_match(lig_atoms: Dict[str, PDB.Vector.Vector],
                   cand: PDB.Residue.Residue) -> bool:
    """
    Return True if candidate residue matches ligand atom names and positions within EXACT_EPS.
    """
    for name, vec in lig_atoms.items():
        if name not in cand:
            return False
        if (vec - cand[name].get_vector()).norm() > EXACT_EPS:
            return False
    return True


def find_substrate_residues(complex_struct, substrate_struct) -> List[PDB.Residue.Residue]:
    """
    Find substrate residues in the complex by **exact coordinate match** to a substrate PDB.
    """
    substrate_res_list = list(substrate_struct.get_residues())
    matched: List[PDB.Residue.Residue] = []
    for lig in substrate_res_list:
        lig_name = lig.get_resname()
        lig_atoms = {a.get_name(): a.get_vector() for a in lig}
        candidates = [r for r in complex_struct.get_residues()
                      if r.get_resname() == lig_name and len(r) == len(lig_atoms)]
        for cand in candidates:
            if is_exact_match(lig_atoms, cand):
                matched.append(cand)
                break
        else:
            chain_id = lig.get_full_id()[2] if len(lig.get_full_id()) > 2 else ""
            resseq = lig.id[1]
            icode = lig.id[2] if len(lig.id) > 2 else " "
            icode_str = "" if icode == " " else icode
            raise ValueError(
                f"Exact match not found for substrate residue {lig_name} chain {chain_id} {resseq}{icode_str}"
            )
    return matched


# ---------- Residue‑ID–based substrate selection ----------

_RES_TOKEN_RE = re.compile(r"""
    ^\s*
    (?:(?P<chain>[^:\s,]+)\s*:\s*)?   # optional chain like A or A_long
    (?P<resseq>\d+)                   # residue sequence number
    (?P<icode>[A-Za-z]?)              # optional insertion code (single letter)
    \s*$
""", re.VERBOSE)

def _parse_res_tokens(spec: str) -> List[Tuple[str | None, int, str | None]]:
    """
    Parse a residue specification string into (chain, resseq, icode) tuples.
    """
    if not spec or not spec.strip():
        raise ValueError("Empty -c/--center specification.")
    tokens = [t.strip() for t in re.split(r"[,\s]+", spec) if t.strip()]
    parsed: List[Tuple[str | None, int, str | None]] = []
    for tok in tokens:
        m = _RES_TOKEN_RE.match(tok)
        if not m:
            raise ValueError(
                f"Invalid residue specifier '{tok}'. Use '123', '123A', 'A:123', or 'A:123A'."
            )
        chain = m.group("chain")
        resseq = int(m.group("resseq"))
        icode = m.group("icode") or None
        parsed.append((chain, resseq, icode))
    return parsed


def find_substrate_by_idspec(complex_struct, spec: str) -> List[PDB.Residue.Residue]:
    """
    Resolve a comma/space-separated residue list into residues within the complex.

    Matching rules
    --------------
    * Chain may be omitted (matches all chains).
    * Insertion code may be omitted (matches any insertion code for that resseq).

    Returns
    -------
    list[Bio.PDB.Residue.Residue]
    """
    targets = _parse_res_tokens(spec)
    found: List[PDB.Residue.Residue] = []
    seen: Set[Tuple] = set()

    for chain_req, resseq_req, icode_req in targets:
        matches: List[PDB.Residue.Residue] = []
        for model in complex_struct:
            for chain in model:
                if chain_req is not None and chain.id != chain_req:
                    continue
                for res in chain.get_residues():
                    _, resseq, icode = res.id
                    if resseq != resseq_req:
                        continue
                    if icode_req is not None and icode != icode_req:
                        continue
                    fid = res.get_full_id()
                    if fid not in seen:
                        seen.add(fid)
                        matches.append(res)
        if not matches:
            chain_txt = f"{chain_req}:" if chain_req is not None else ""
            icode_txt = icode_req or ""
            raise ValueError(f"Residue '{chain_txt}{resseq_req}{icode_txt}' not found in complex.")
        found.extend(matches)

    return found

# ---------- Residue-name-based substrate selection ----------

def find_substrate_by_resname(complex_struct, spec: str) -> List[PDB.Residue.Residue]:
    """
    Resolve a comma/space-separated residue-name list (e.g., 'GPP,MMT') into residues in the complex.

    Behavior
    --------
    * Case-insensitive match against residue `resname`.
    * If multiple residues share the same name, **all** are included and a **WARNING** is logged.
    """
    if not spec or not spec.strip():
        raise ValueError("Empty -c/--center specification.")
    tokens = [t.strip().upper() for t in re.split(r"[,\s]+", spec) if t.strip()]
    found: List[PDB.Residue.Residue] = []
    seen_fids: Set[Tuple] = set()
    for rn in tokens:
        matches = [r for r in complex_struct.get_residues() if r.get_resname().upper() == rn]
        if not matches:
            raise ValueError(f"Residue name '{rn}' not found in complex.")
        if len(matches) > 1:
            try:
                sample = ", ".join(_fmt_res_id(r) for r in matches[:5])
            except Exception:
                sample = "(list omitted)"
            LOGGER.warning("[extract] Multiple residues with resname '%s' found (%d). Using all: %s",
                            rn, len(matches), sample)
        for r in matches:
            fid = r.get_full_id()
            if fid not in seen_fids:
                seen_fids.add(fid)
                found.append(r)
    return found


def resolve_substrate_residues(complex_struct, center_spec: str) -> List[PDB.Residue.Residue]:
    """
    Determine substrate residues from a PDB path, residue-ID list, or residue-name list.
    """
    if os.path.exists(center_spec):
        substrate_struct = load_structure(center_spec, "substrate")
        return find_substrate_residues(complex_struct, substrate_struct)
    # If it parses as ID-spec, treat as IDs (and propagate any not-found errors).
    try:
        _parse_res_tokens(center_spec)
        return find_substrate_by_idspec(complex_struct, center_spec)
    except ValueError:
        # Otherwise, interpret as residue-name list (e.g., 'GPP,MMT').
        return find_substrate_by_resname(complex_struct, center_spec)


# ---------------------------------------------------------------------
#   Polypeptide adjacency (C–N) helper
# ---------------------------------------------------------------------

def are_peptide_adjacent(prev_res: PDB.Residue.Residue,
                         next_res: PDB.Residue.Residue,
                         max_cn_dist: float = 1.9) -> bool:
    """
    Return True if prev_res—next_res are peptide-bond adjacent based on C(prev)–N(next) distance.

    Notes
    -----
    Distance‑based criterion; in practice this avoids crossing TER boundaries because missing
    atoms or long inter‑residue distances will fail the check.
    """
    if prev_res.get_resname() not in AMINO_ACIDS or next_res.get_resname() not in AMINO_ACIDS:
        return False
    if ("C" not in prev_res) or ("N" not in next_res):
        return False
    try:
        d = (prev_res["C"].get_vector() - next_res["N"].get_vector()).norm()
    except Exception:
        return False
    return (d == d) and (d <= max_cn_dist)  # d==d to filter NaN


# ---------------------------------------------------------------------
#   Residue selection around the substrate
# ---------------------------------------------------------------------

def select_residues(complex_struct,
                    substrate_res_list: List[PDB.Residue.Residue],
                    r_as: float,
                    r_het: float,
                    include_h2o: bool,
                    exclude_backbone: bool) -> Tuple[Set[Tuple], Set[Tuple]]:
    """
    Select pocket residues around the substrate.

    Selection rule
    --------------
    * Always include the substrate residues themselves.
    * Standard cutoff (`r_as`):
        - If `exclude_backbone` is **False**: include a residue if **any** atom is within `r_as`.
        - If `exclude_backbone` is **True**: for **amino acids**, require a **non‑backbone** atom
          to be within `r_as`; non‑amino‑acid residues are included if **any** atom is within `r_as`.
    * Hetero‑hetero cutoff (`r_het`):
        - Neighbor atom must be hetero (element not in {C,H}).
        - When `exclude_backbone` is **True** and the neighbor is an amino acid, that atom must
          also be **non‑backbone**.

    Returns
    -------
    (selected_ids, backbone_contact_ids)
      selected_ids : set of residue full-ids to output
      backbone_contact_ids : subset with any **backbone atom** within r_as or r_het of a substrate atom.
                             (Waters ignored; only relevant when exclude_backbone == False)
    """
    substrate_atoms = [a for lig in substrate_res_list for a in lig]
    substrate_het = [a for a in substrate_atoms if a.element not in ("C", "H")]
    ns = NeighborSearch(list(complex_struct.get_atoms()))

    selected_ids: Set[Tuple] = {res.get_full_id() for res in substrate_res_list}
    backbone_contact_ids: Set[Tuple] = set()

    def is_amino_backbone_atom(atom: PDB.Atom.Atom) -> bool:
        res = atom.get_parent()
        return (res.get_resname() in AMINO_ACIDS) and (atom.get_name() in BACKBONE_ATOMS)

    def maybe_add(atom, via_backbone: bool):
        res = atom.get_parent()
        if not include_h2o and res.get_resname() in WATER_RES:
            return
        fid = res.get_full_id()
        selected_ids.add(fid)
        if via_backbone and res.get_resname() in AMINO_ACIDS:
            backbone_contact_ids.add(fid)

    # standard radius: any atom within r_as (with backbone filter when exclude_backbone==True)
    for atom in substrate_atoms:
        for neigh in ns.search(atom.get_coord(), r_as):
            if exclude_backbone and is_amino_backbone_atom(neigh):
                continue  # require non-backbone atom for amino-acid residues
            via_backbone_neigh = (neigh.get_name() in BACKBONE_ATOMS)
            maybe_add(neigh, via_backbone_neigh)

    # hetero-hetero radius: both sides non-C/H (and non-backbone filter for amino acids when exclude_backbone==True)
    for atom in substrate_het:
        for neigh in ns.search(atom.get_coord(), r_het):
            if neigh.element in ("C", "H"):
                continue
            if exclude_backbone and is_amino_backbone_atom(neigh):
                continue
            via_backbone_neigh = (neigh.get_name() in BACKBONE_ATOMS)
            maybe_add(neigh, via_backbone_neigh)

    return selected_ids, backbone_contact_ids


# ---------------------------------------------------------------------
#   Disulfide augmentation
# ---------------------------------------------------------------------

def augment_disulfides(structure, selected_ids: Set[Tuple],
                       cutoff: float = DISULFIDE_CUTOFF):
    """
    Include Cys–Cys disulfide partners if either residue is selected (SG–SG ≤ cutoff).
    """
    sg_atoms = [r["SG"] for r in structure.get_residues()
                if r.get_resname() in {"CYS", "CYX"} and "SG" in r]

    if not sg_atoms:
        return

    ns = NeighborSearch(sg_atoms)
    for at in sg_atoms:
        for other in ns.search(at.get_coord(), cutoff):
            if other is at:
                continue
            f1 = at.get_parent().get_full_id()
            f2 = other.get_parent().get_full_id()
            if f1 in selected_ids or f2 in selected_ids:
                selected_ids.update((f1, f2))


# ---------------------------------------------------------------------
#   Proline augmentation (N-side neighbor inclusion; TER-aware)
# ---------------------------------------------------------------------

def augment_proline_prev_neighbor(structure, selected_ids: Set[Tuple]):
    """
    Ensure that if a selected PRO is not at the N-terminus, the immediately
    preceding (N-side) amino-acid residue is also selected.

    Notes
    -----
    Uses peptide adjacency (C–N ≤ 1.9 Å) to avoid crossing TER boundaries.
    """
    added = 0
    for fid in list(selected_ids):
        model_id, chain_id, res_id = fid[1], fid[2], fid[3]
        res: PDB.Residue.Residue = structure[model_id][chain_id].child_dict[res_id]
        if res.get_resname() != "PRO":
            continue
        chain = structure[model_id][chain_id]
        residues: List[PDB.Residue.Residue] = list(chain.get_residues())
        try:
            idx = next(i for i, r in enumerate(residues) if r.get_full_id() == fid)
        except StopIteration:
            continue
        prev_res = None
        for j in range(idx - 1, -1, -1):
            rj = residues[j]
            if rj.get_resname() in AMINO_ACIDS:
                prev_res = rj
                break
        if prev_res is None:
            continue
        if not are_peptide_adjacent(prev_res, res):
            continue
        prev_fid = prev_res.get_full_id()
        if prev_fid not in selected_ids:
            selected_ids.add(prev_fid)
            added += 1
    if added:
        LOGGER.info("[extract] Added %d N-side neighbor residues for PRO (TER-aware).", added)


# ---------------------------------------------------------------------
#   Backbone-contact neighbor augmentation (exclude_backbone == False; TER-aware)
# ---------------------------------------------------------------------

def augment_backbone_contact_neighbors(structure,
                                       selected_ids: Set[Tuple],
                                       backbone_contact_ids: Set[Tuple],
                                       substrate_ids: Set[Tuple]) -> Tuple[Set[Tuple], Set[Tuple]]:
    """
    If a non-substrate residue had **any backbone atom** within selection radii,
    include its immediate N- and C-side amino-acid neighbors **only if peptide-bond adjacent**.

    If a side has no peptide-adjacent neighbor (true terminus; e.g., separated by TER),
    mark the residue to **keep** the respective terminal atoms (N/H* for N-terminus; C/O/OXT for C-terminus).

    Returns
    -------
    keep_ncap_ids, keep_ccap_ids : sets of full-ids whose terminal caps must be preserved
    """
    keep_ncap_ids: Set[Tuple] = set()
    keep_ccap_ids: Set[Tuple] = set()
    added = 0
    termini_kept_n = 0
    termini_kept_c = 0

    for fid in list(backbone_contact_ids):
        if fid in substrate_ids:
            continue  # do not augment around substrate residues
        model_id, chain_id, res_id = fid[1], fid[2], fid[3]
        chain = structure[model_id][chain_id]
        residues: List[PDB.Residue.Residue] = list(chain.get_residues())
        try:
            idx = next(i for i, r in enumerate(residues) if r.get_full_id() == fid)
        except StopIteration:
            continue

        cur_res = residues[idx]

        # previous amino-acid — require peptide adjacency
        prev_res = None
        for j in range(idx - 1, -1, -1):
            rj = residues[j]
            if rj.get_resname() in AMINO_ACIDS:
                prev_res = rj
                break
        if prev_res is not None and are_peptide_adjacent(prev_res, cur_res):
            prev_fid = prev_res.get_full_id()
            if prev_fid not in selected_ids:
                selected_ids.add(prev_fid)
                added += 1
        else:
            keep_ncap_ids.add(fid)
            termini_kept_n += 1

        # next amino-acid — require peptide adjacency
        next_res = None
        for j in range(idx + 1, len(residues)):
            rj = residues[j]
            if rj.get_resname() in AMINO_ACIDS:
                next_res = rj
                break
        if next_res is not None and are_peptide_adjacent(cur_res, next_res):
            next_fid = next_res.get_full_id()
            if next_fid not in selected_ids:
                selected_ids.add(next_fid)
                added += 1
        else:
            keep_ccap_ids.add(fid)
            termini_kept_c += 1

    if added or termini_kept_n or termini_kept_c:
        LOGGER.info("[extract] Backbone-contact context (TER-aware): added %d neighbors; kept N-cap on %d, C-cap on %d residues.",
                     added, termini_kept_n, termini_kept_c)
    return keep_ncap_ids, keep_ccap_ids


# ---------------------------------------------------------------------
#   Backbone trimming / skip-map generation
# ---------------------------------------------------------------------

def mark_atoms_to_skip(structure, selected_ids: Set[Tuple], substrate_ids: Set[Tuple],
                       exclude_backbone: bool,
                       keep_ncap_ids: Set[Tuple] | None = None,
                       keep_ccap_ids: Set[Tuple] | None = None) -> Dict[Tuple, Set[str]]:
    """
    Decide which atoms to delete (truncation). Never delete substrate atoms.

    Returns
    -------
    dict[full-id -> set(atom_names_to_delete)]
    """
    keep_ncap_ids = keep_ncap_ids or set()
    keep_ccap_ids = keep_ccap_ids or set()

    # start with the original truncation logic (except for substrate residues)
    chain_map: Dict[Tuple[str, str], List[Tuple]] = {}
    for fid in selected_ids:
        if fid in substrate_ids:
            continue  # never delete atoms from substrate residues
        res = structure[fid[1]][fid[2]].child_dict[fid[3]]
        if res.get_resname() in WATER_RES:
            continue
        chain_map.setdefault((fid[1], fid[2]), []).append(fid)

    skip: Dict[Tuple, Set[str]] = {}

    # --- TER-aware segmentation: split by peptide adjacency in file order ---
    for (model, chain), fids in chain_map.items():
        chain_obj = structure[model][chain]
        residues_all: List[PDB.Residue.Residue] = list(chain_obj.get_residues())
        index_map: Dict[Tuple, int] = {r.get_full_id(): i for i, r in enumerate(residues_all)}

        # sort by file order
        fids.sort(key=lambda x: index_map.get(x, 10**9))

        # build segments by peptide-bond adjacency
        segs: List[List[Tuple]] = []
        cur_seg: List[Tuple] = []
        for k, fid in enumerate(fids):
            if not cur_seg:
                cur_seg = [fid]
                continue
            prev_fid = cur_seg[-1]
            prev_res = chain_obj.child_dict[prev_fid[3]]
            cur_res = chain_obj.child_dict[fid[3]]
            if are_peptide_adjacent(prev_res, cur_res):
                cur_seg.append(fid)
            else:
                segs.append(cur_seg)
                cur_seg = [fid]
        if cur_seg:
            segs.append(cur_seg)

        # apply cap deletions on these TER-aware segments
        for seg in segs:
            n_id, c_id = seg[0], seg[-1]
            single = len(seg) == 1

            def add(fid_local, names):
                skip.setdefault(fid_local, set()).update(names)

            n_res = chain_obj.child_dict[n_id[3]]
            c_res = chain_obj.child_dict[c_id[3]]

            # N-terminal cap deletion (only for amino acids; skip if PRO/HYP or explicitly kept)
            if (n_res.get_resname() in AMINO_ACIDS) and (n_res.get_resname() not in {"PRO", "HYP"}) and (n_id not in keep_ncap_ids):
                add(n_id, {"N", "H", "H1", "H2", "H3", "HN"})
            # C-terminal cap deletion (only for amino acids; skip if explicitly kept)
            if (c_res.get_resname() in AMINO_ACIDS) and (c_id not in keep_ccap_ids):
                add(c_id, {"C", "O", "OXT"})

            # Isolated stretch – remove CA/HA* (only for amino acids; except PRO/HYP)
            if single and (n_res.get_resname() in AMINO_ACIDS) and (n_res.get_resname() not in {"PRO", "HYP"}):
                add(n_id, {"CA", "HA", "HA2", "HA3"})

    # ---------------------------------------------------------------------
    #   Optional: remove *all* backbone atoms from every non-substrate residue
    #             PRO/HYP keep N, CA, and HA* to preserve the ring.
    # ---------------------------------------------------------------------
    if exclude_backbone:
        for fid in selected_ids:
            if fid in substrate_ids:
                continue
            res = structure[fid[1]][fid[2]].child_dict[fid[3]]
            if res.get_resname() in WATER_RES:
                continue
            if res.get_resname() in AMINO_ACIDS:
                if res.get_resname() in {"PRO", "HYP"}:
                    to_remove = BACKBONE_ALL - {"N", "CA", "HA", "H", "H1", "H2", "H3"}
                else:
                    to_remove = BACKBONE_ALL
                skip.setdefault(fid, set()).update(to_remove)

        # Preserve peptide carbonyl on the N-side neighbor of PRO
        for fid in selected_ids:
            res = structure[fid[1]][fid[2]].child_dict[fid[3]]
            if res.get_resname() != "PRO":
                continue
            chain = structure[fid[1]][fid[2]]
            residues: List[PDB.Residue.Residue] = list(chain.get_residues())
            try:
                idx = next(i for i, r in enumerate(residues) if r.get_full_id() == fid)
            except StopIteration:
                continue
            prev_res = None
            for j in range(idx - 1, -1, -1):
                rj = residues[j]
                if rj.get_resname() in AMINO_ACIDS:
                    prev_res = rj
                    break
            if prev_res is None:
                continue
            if not are_peptide_adjacent(prev_res, res):
                continue
            prev_fid = prev_res.get_full_id()
            if prev_fid in selected_ids:
                sk = skip.setdefault(prev_fid, set())
                for nm in ("C", "O", "OXT"):
                    if nm in sk:
                        sk.remove(nm)

    # Always keep CA on the N-side neighbor of PRO (independent of --exclude-backbone)
    for fid in selected_ids:
        res = structure[fid[1]][fid[2]].child_dict[fid[3]]
        if res.get_resname() != "PRO":
            continue
        chain = structure[fid[1]][fid[2]]
        residues: List[PDB.Residue.Residue] = list(chain.get_residues())
        try:
            idx = next(i for i, r in enumerate(residues) if r.get_full_id() == fid)
        except StopIteration:
            continue
        prev_res = None
        for j in range(idx - 1, -1, -1):
            rj = residues[j]
            if rj.get_resname() in AMINO_ACIDS:
                prev_res = rj
                break
        if prev_res is None:
            continue
        if not are_peptide_adjacent(prev_res, res):
            continue
        prev_fid = prev_res.get_full_id()
        if prev_fid in selected_ids:
            sk = skip.setdefault(prev_fid, set())
            if "CA" in sk:
                sk.remove("CA")

    return skip


def _atom_present_in_output(res: PDB.Residue.Residue, name: str, skip_set: Set[str]) -> bool:
    """
    True if the atom exists originally AND is not marked for deletion.
    """
    return (name in res) and (name not in skip_set)

def _atom_removed_by_truncation(res: PDB.Residue.Residue, name: str, skip_set: Set[str]) -> bool:
    """
    True if the atom exists originally AND is marked for deletion.
    """
    return (name in res) and (name in skip_set)

def compute_linkH_atoms(structure,
                        selected_ids: Set[Tuple],
                        skip_map: Dict[Tuple, Set[str]]) -> List[Tuple[float, float, float]]:
    """
    Identify severed bonds created by truncation and compute link‑H coordinates.

    Rules
    -----
    * Normal residues: place H along **CB→CA**, **CA→N**, **CA→C** if partner was removed.
    * PRO/HYP: place H along **CA→C** only.
    * Parent atom **must be Carbon**; H is placed along (parent → removed_partner) at **1.09 Å**.

    Returns
    -------
    list of (x, y, z) coordinates for link‑H atoms
    """
    link_coords: List[Tuple[float, float, float]] = []

    for fid in selected_ids:
        model_id, chain_id, res_id = fid[1], fid[2], fid[3]
        res: PDB.Residue.Residue = structure[model_id][chain_id].child_dict[res_id]
        if res.get_resname() in WATER_RES:
            continue
        skip_set = skip_map.get(fid, set())
        resname = res.get_resname()

        def _add_if_cut(parent_name: str, partner_name: str):
            if not _atom_present_in_output(res, parent_name, skip_set):
                return
            if not _atom_removed_by_truncation(res, partner_name, skip_set):
                return
            parent = res[parent_name]
            partner = res[partner_name]
            parent_elem = (parent.element or parent.get_name()[0]).upper()
            if parent_elem != "C":
                return
            v = np.array(partner.get_coord(), dtype=float) - np.array(parent.get_coord(), dtype=float)
            norm = np.linalg.norm(v)
            if not np.isfinite(norm) or norm < 1e-6:
                return
            v /= norm
            dist = 1.09  # C–H
            h = np.array(parent.get_coord(), dtype=float) + v * dist
            link_coords.append((float(h[0]), float(h[1]), float(h[2])))

        if resname in {"PRO", "HYP"}:
            _add_if_cut("CA", "C")
        else:
            _add_if_cut("CB", "CA")
            _add_if_cut("CA", "N")
            _add_if_cut("CA", "C")

    return link_coords


def _max_serial_from_pdb_text(pdb_text: str) -> int:
    """
    Find the maximum atom serial number in PDB text.
    """
    max_serial = 0
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                serial = int(line[6:11])
                if serial > max_serial:
                    max_serial = serial
            except Exception:
                continue
    return max_serial


def _format_linkH_block(link_coords: List[Tuple[float, float, float]],
                        start_serial: int,
                        chain_id: str = "L") -> str:
    """
    Format a contiguous HETATM block for link‑H atoms.

    Conventions
    -----------
    * Atom name: HL
    * Residue name: LKH
    * Chain: chain_id (default 'L')
    * Residue numbers: 1..N (one pseudo‑residue per H)
    """
    lines: List[str] = []
    serial = start_serial
    resseq = 1
    for (x, y, z) in link_coords:
        serial += 1
        line = (
            f"HETATM{serial:5d} "
            f"{'HL':>4s} "
            f"{'LKH':>3s} "
            f"{chain_id}"
            f"{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}"
            f"          {'H':>2s}"
        )
        lines.append(line)
        resseq += 1
    return ("\n".join(lines) + ("\n" if lines else ""))


# ---------------------------------------------------------------------
#   Charge calculation & logging
# ---------------------------------------------------------------------

def _sorted_fids_by_file_order(structure, fids: Iterable[Tuple]) -> List[Tuple]:
    """
    Sort full-ids by file order using a residue index map.
    """
    order: Dict[Tuple, int] = {}
    idx = 0
    for model in structure:
        for chain in model:
            for res in chain.get_residues():
                order[res.get_full_id()] = idx
                idx += 1
    return sorted(set(fids), key=lambda fid: order.get(fid, 10**12))

def _residue_key_from_res(res: PDB.Residue.Residue) -> ResidueKey:
    """
    Build a cross-structure residue key from a residue.
    """
    chain_id = res.get_parent().id
    hetflag, resseq, icode = res.id
    icode_str = icode if icode != " " else ""
    return (chain_id, hetflag, int(resseq), icode_str, res.get_resname())

def _residue_key_from_fid(structure, fid: Tuple) -> ResidueKey:
    """
    Build a cross-structure residue key from a full-id.
    """
    res = structure[fid[1]][fid[2]].child_dict[fid[3]]
    return _residue_key_from_res(res)

# ---- helper for parsing --ligand-charge (number or 'RES:Q' mapping) ----
def _parse_ligand_charge_option(ligand_charge: float | str | Dict[str, float] | None
                                ) -> Tuple[Optional[float], Optional[Dict[str, float]]]:
    """
    Returns
    -------
    (total_charge, mapping)
      total_charge : float | None
      mapping      : dict[RESNAME -> float] | None
    """
    if ligand_charge is None:
        return None, None
    if isinstance(ligand_charge, (int, float)):
        return float(ligand_charge), None
    if isinstance(ligand_charge, dict):
        mapping = {str(k).upper(): float(v) for k, v in ligand_charge.items()}
        return None, mapping
    if isinstance(ligand_charge, str):
        s = ligand_charge.strip()
        if not s:
            return None, None
        # try numeric
        try:
            return float(s), None
        except ValueError:
            pass
        # mapping: tokens "RES:Q"
        tokens = [t for t in re.split(r"[,\s]+", s) if t]
        mapping: Dict[str, float] = {}
        for tok in tokens:
            if ":" not in tok:
                raise ValueError(f"Invalid --ligand-charge token '{tok}'. Use 'RES:Q' (e.g., GPP:-3) or a number (e.g., -3).")
            res, qtxt = tok.split(":", 1)
            resname = res.strip().upper()
            if not resname:
                raise ValueError(f"Invalid --ligand-charge token '{tok}': empty residue name.")
            try:
                qval = float(qtxt.strip())
            except ValueError:
                raise ValueError(f"Invalid --ligand-charge token '{tok}': '{qtxt}' is not a number.")
            mapping[resname] = qval
        if not mapping:
            raise ValueError("Empty --ligand-charge mapping.")
        return None, mapping
    raise TypeError(f"Unsupported type for ligand_charge: {type(ligand_charge)!r}")

def compute_charge_summary(structure,
                           selected_ids: Set[Tuple],
                           substrate_ids: Set[Tuple],
                           ligand_charge: float | str | Dict[str, float] | None = None) -> Dict[str, Any]:
    """
    Compute pocket charge summary.

    Args
    ----
    structure : Bio.PDB.Structure.Structure
        The (first) structure to evaluate.
    selected_ids : set[tuple]
        Residues included in the pocket.
    substrate_ids : set[tuple]
        Residues designated as substrate.
    ligand_charge : float | str | dict[str,float] | None
        - float: total charge to assign across **unknown residues** (preferring unknown substrate).
        - str  : numeric string (total) or mapping like "GPP:-3,MMT:-1" (per‑resname).
        - dict : mapping {RESNAME: charge}. In mapping mode, other unknown residues remain 0.

    Returns
    -------
    dict with keys:
      - total_charge : float
      - protein_charge : float
      - ligand_total_charge : float
      - ion_total_charge : float
      - ion_charges : list[(str tag, float)]
      - unknown_residue_charges : dict[str -> float]  # for concise per‑resname log
    """
    per_map: Dict[ResidueKey, float] = {}
    aa_charge = 0.0
    total = 0.0

    fids_in_order = _sorted_fids_by_file_order(structure, selected_ids)

    # First pass: dictionary/ion/water charges; collect unknowns and ions
    unknown_fids: List[Tuple] = []
    unknown_substrate_fids: List[Tuple] = []
    ion_entries: List[Tuple[str, float]] = []

    for fid in fids_in_order:
        res = structure[fid[1]][fid[2]].child_dict[fid[3]]
        rn = res.get_resname().upper()
        key = _residue_key_from_res(res)
        if rn in WATER_RES:
            q = 0.0
        elif rn in AMINO_ACIDS:
            q = float(AMINO_ACIDS[rn])
            aa_charge += q
        elif rn in ION:
            q = float(ION[rn])
            ion_entries.append((_fmt_fid(structure, fid), q))
        else:
            q = 0.0
            unknown_fids.append(fid)
            if fid in substrate_ids:
                unknown_substrate_fids.append(fid)
        per_map[key] = q
        total += q

    # Apply --ligand-charge if provided
    total_spec, mapping_spec = _parse_ligand_charge_option(ligand_charge)

    if total_spec is not None:
        # Distribute total across unknown substrate if present, else across all unknowns
        targets = unknown_substrate_fids if unknown_substrate_fids else unknown_fids
        if targets:
            per_res_val = float(total_spec) / float(len(targets))
            for fid in targets:
                key = _residue_key_from_fid(structure, fid)
                per_map[key] = per_res_val
            # recompute totals
            total = sum(per_map.values())
            aa_charge = sum(q for k, q in per_map.items() if k[4] in AMINO_ACIDS)
    elif mapping_spec is not None:
        # Per‑resname mapping. Unspecified unknown residues remain 0.
        for fid in unknown_fids:
            res = structure[fid[1]][fid[2]].child_dict[fid[3]]
            rn = res.get_resname().upper()
            if rn in mapping_spec:
                key = _residue_key_from_fid(structure, fid)
                per_map[key] = float(mapping_spec[rn])
        # recompute totals
        total = sum(per_map.values())
        aa_charge = sum(q for k, q in per_map.items() if k[4] in AMINO_ACIDS)

    # Net ligand and ion charges
    unknown_keys = {_residue_key_from_fid(structure, fid) for fid in unknown_fids}
    ligand_total = sum(per_map[k] for k in unknown_keys)
    ion_total = sum(q for _, q in ion_entries)

    # Build per‑resname mapping for unknown residues (after applying any overrides)
    unknown_residue_charges: Dict[str, float] = {}
    for fid in unknown_fids:
        res = structure[fid[1]][fid[2]].child_dict[fid[3]]
        rn = res.get_resname().upper()
        key = _residue_key_from_fid(structure, fid)
        unknown_residue_charges[rn] = float(per_map[key])

    return {
        "total_charge": float(total),
        "protein_charge": float(aa_charge),
        "ligand_total_charge": float(ligand_total),
        "ion_total_charge": float(ion_total),
        "ion_charges": [(tag, float(q)) for tag, q in ion_entries],
        "unknown_residue_charges": unknown_residue_charges,
    }

def log_charge_summary(prefix: str,
                       summary: Dict[str, Any]):
    """
    Emit concise charge summary logs.
    """
    total = summary["total_charge"]
    protein = summary["protein_charge"]
    ligand = summary.get("ligand_total_charge", 0.0)
    ion_list: List[Tuple[str, float]] = summary.get("ion_charges", [])
    ion_total = summary.get("ion_total_charge", sum(q for _, q in ion_list))
    unk_map: Dict[str, float] = summary.get("unknown_residue_charges", {}) or {}

    if unk_map:
        items = ", ".join(f"{res}: {q:g}" for res, q in sorted(unk_map.items()))
        LOGGER.info("%s Per-resname ligand charges: %s", prefix, items)
    else:
        LOGGER.info("%s Per-resname ligand charges: (none)", prefix)

    LOGGER.info("%s Net protein charge: %+g", prefix, protein)
    LOGGER.info("%s Net ligand charge: %+g", prefix, ligand)
    if ion_list:
        LOGGER.info("%s Ion charges (each):", prefix)
        for tag, q in ion_list:
            LOGGER.info("  %s  ->  %+g", tag, q)
        LOGGER.info("%s Net ion charge: %+g", prefix, ion_total)
    else:
        LOGGER.info("%s Ion charges: (none)", prefix)
    LOGGER.info("%s Total pocket charge: %+g", prefix, total)


# =========================== Cross-structure helpers ===========================
#   Multi-model driver utilities
# ==============================================================================

def _build_key_maps(structure) -> Tuple[Dict[ResidueKey, Tuple], Dict[Tuple, ResidueKey]]:
    """
    Create maps between ResidueKey and full-id for a structure.
    """
    key2fid: Dict[ResidueKey, Tuple] = {}
    fid2key: Dict[Tuple, ResidueKey] = {}
    for model in structure:
        for chain in model:
            for res in chain.get_residues():
                key = _residue_key_from_res(res)
                fid = res.get_full_id()
                key2fid[key] = fid
                fid2key[fid] = key
    return key2fid, fid2key

def _keys_to_fids(structure, keys: Iterable[ResidueKey]) -> Set[Tuple]:
    """
    Translate a set of ResidueKeys into full-ids for this structure.
    """
    key2fid, _ = _build_key_maps(structure)
    fids: Set[Tuple] = set()
    missing: List[ResidueKey] = []
    for k in keys:
        fid = key2fid.get(k)
        if fid is None:
            missing.append(k)
        else:
            fids.add(fid)
    if missing:
        raise ValueError(f"Some residues not found in structure: {missing[:5]}{' ...' if len(missing)>5 else ''}")
    return fids

def _fids_to_keys(structure, fids: Iterable[Tuple]) -> Set[ResidueKey]:
    """
    Translate a set of full-ids into ResidueKeys.
    """
    return {_residue_key_from_fid(structure, fid) for fid in fids}

def _substrate_residues_for_structs(structs: List[PDB.Structure.Structure],
                                    center_spec: str) -> List[List[PDB.Residue.Residue]]:
    """
    Resolve substrate residues per structure.

    Behavior
    --------
    * If `center_spec` is a PDB path: exact‑match on the first structure only,
      then propagate to others by a residue‑ID list derived from the first match.
    * If `center_spec` is an ID list: apply to all structures.
    * If `center_spec` is a residue‑name list: apply to all structures; names may match multiple residues
      (all included; WARNING logged per structure).
    """
    if os.path.exists(center_spec):
        sub_first = resolve_substrate_residues(structs[0], center_spec)
        tokens = []
        for res in sub_first:
            chain = res.get_parent().id
            chain_txt = (chain or "").strip()
            het, num, icode = res.id
            icode_txt = "" if icode == " " else icode
            if chain_txt:
                tokens.append(f"{chain}:{num}{icode_txt}")
            else:
                tokens.append(f"{num}{icode_txt}")
        idspec = ",".join(tokens)
        out: List[List[PDB.Residue.Residue]] = []
        for si, st in enumerate(structs):
            out.append(find_substrate_by_idspec(st, idspec))
        return out
    else:
        # Distinguish ID-spec vs resname list by attempting to parse as IDs first.
        try:
            _parse_res_tokens(center_spec)
            return [find_substrate_by_idspec(st, center_spec) for st in structs]
        except ValueError:
            return [find_substrate_by_resname(st, center_spec) for st in structs]

def _disulfide_partner_keys(structure, candidate_keys: Set[ResidueKey],
                            cutoff: float = DISULFIDE_CUTOFF) -> Set[ResidueKey]:
    """
    Return ResidueKeys of disulfide partners to include for any selected CYS/CYX.
    """
    key2fid, _ = _build_key_maps(structure)
    sg_atoms: List[PDB.Atom.Atom] = []
    res_of_atom: Dict[PDB.Atom.Atom, ResidueKey] = {}
    for res in structure.get_residues():
        if res.get_resname() in {"CYS", "CYX"} and "SG" in res:
            at = res["SG"]
            sg_atoms.append(at)
            res_of_atom[at] = _residue_key_from_res(res)
    add: Set[ResidueKey] = set()
    if not sg_atoms:
        return add
    ns = NeighborSearch(sg_atoms)
    for at in sg_atoms:
        for other in ns.search(at.get_coord(), cutoff):
            if other is at:
                continue
            k1 = res_of_atom[at]
            k2 = res_of_atom[other]
            if (k1 in candidate_keys) or (k2 in candidate_keys):
                add.add(k1); add.add(k2)
    return add

def _assert_atom_ordering_identical(structs: List[PDB.Structure.Structure]):
    """
    Light consistency check across inputs:
    - Enforce identical atom counts.
    - Spot‑check ordering at the beginning and end of the atom list; if mismatched there (and overall lists differ),
      raise an error.
    """
    def signature(st: PDB.Structure.Structure) -> List[str]:
        sig: List[str] = []
        for model in st:
            for chain in model:
                for res in chain.get_residues():
                    het, resseq, icode = res.id
                    icode_txt = icode if icode != " " else ""
                    base = f"{chain.id}|{het}|{resseq}{icode_txt}|{res.get_resname()}"
                    for atom in res:
                        sig.append(base + f"|{atom.get_name()}")
        return sig
    sig0 = signature(structs[0])
    for i in range(1, len(structs)):
        sigi = signature(structs[i])
        if len(sigi) != len(sig0):
            raise ValueError(f"[multi] Atom count mismatch between input #1 and input #{i+1}: {len(sig0)} vs {len(sigi)}")
        check_pairs = [(0, min(10, len(sig0))),
                       (max(0, len(sig0)-10), len(sig0))]
        mismatch = False
        for a, b in check_pairs:
            if sig0[a:b] != sigi[a:b]:
                mismatch = True
                break
        if mismatch and sig0 != sigi:
            raise ValueError(f"[multi] Atom order mismatch between input #1 and input #{i+1}.")


def _strip_trailing_END(text: str) -> str:
    """
    Remove trailing 'END' lines and ensure a final newline.
    """
    lines = [ln for ln in text.splitlines() if ln.strip() != "END"]
    out = "\n".join(lines)
    if not out.endswith("\n"):
        out += "\n"
    return out


def _compute_linkH_defs(structure,
                        selected_ids: Set[Tuple],
                        skip_map: Dict[Tuple, Set[str]]) -> List[Tuple[Tuple[ResidueKey, str], Tuple[float, float, float]]]:
    """
    Deterministic list of link‑H definitions and coordinates.

    Returns
    -------
    list of ((ResidueKey, cut_type), (x, y, z)), where cut_type ∈ {"CB-CA","CA-N","CA-C"}.
    Ordering is by residue file order, then by cut_type in the sequence above.
    """
    out: List[Tuple[Tuple[ResidueKey, str], Tuple[float, float, float]]] = []
    for fid in _sorted_fids_by_file_order(structure, selected_ids):
        res: PDB.Residue.Residue = structure[fid[1]][fid[2]].child_dict[fid[3]]
        if res.get_resname() in WATER_RES:
            continue
        skip_set = skip_map.get(fid, set())
        key = _residue_key_from_res(res)

        def _maybe(parent_name: str, partner_name: str, cut_type: str):
            if not _atom_present_in_output(res, parent_name, skip_set):
                return
            if not _atom_removed_by_truncation(res, partner_name, skip_set):
                return
            parent = res[parent_name]
            partner = res[partner_name]
            parent_elem = (parent.element or parent.get_name()[0]).upper()
            if parent_elem != "C":
                return
            v = np.array(partner.get_coord(), dtype=float) - np.array(parent.get_coord(), dtype=float)
            norm = np.linalg.norm(v)
            if not np.isfinite(norm) or norm < 1e-6:
                return
            v /= norm
            dist = 1.09
            h = np.array(parent.get_coord(), dtype=float) + v * dist
            out.append(((key, cut_type), (float(h[0]), float(h[1]), float(h[2]))))

        if res.get_resname() in {"PRO", "HYP"}:
            _maybe("CA", "C", "CA-C")
        else:
            _maybe("CB", "CA", "CB-CA")
            _maybe("CA", "N",  "CA-N")
            _maybe("CA", "C",  "CA-C")

    return out


def extract_multi(args: argparse.Namespace, api=False) -> Dict[str, Any]:
    """
    Multi‑structure driver.

    Args
    ----
    args : argparse.Namespace
        Parsed CLI arguments (or equivalent) controlling selection, truncation, outputs.

    Returns
    -------
    dict
        {
          'outputs': List[str],
          'counts': List[{'raw_atoms': int, 'kept_atoms': int}],  # per model
          'charge_summary': {...},  # computed on model #1
        }
    """
    paths: List[str] = args.complex_pdb
    names = [f"complex{i+1}" for i in range(len(paths))]
    structs: List[PDB.Structure.Structure] = [load_structure(p, n) for p, n in zip(paths, names)]

    LOGGER.info("[extract:multi] Loaded %d structures.", len(structs))
    _assert_atom_ordering_identical(structs)

    # Substrates per structure (PDB-path -> first only, then propagate by IDs)
    subs_per_struct: List[List[PDB.Residue.Residue]] = _substrate_residues_for_structs(structs, args.substrate_pdb)

    # 1) Per-structure selection and backbone-contact → OR unify as keys
    union_sel_keys: Set[ResidueKey] = set()
    union_bb_contact_keys: Set[ResidueKey] = set()

    for st, subs in zip(structs, subs_per_struct):
        selected_ids, bb_contact_ids = select_residues(st, subs, args.radius, args.radius_het2het, args.include_h2o, args.exclude_backbone)
        union_sel_keys |= _fids_to_keys(st, selected_ids)
        union_bb_contact_keys |= _fids_to_keys(st, bb_contact_ids)

    LOGGER.info("[extract:multi] Initial union selection: %d residues; backbone-contact: %d residues.",
                 len(union_sel_keys), len(union_bb_contact_keys))

    # 1a) Force-include residues via --selected-resn (OR across structures)
    if getattr(args, "selected_resn", ""):
        forced_union: Set[ResidueKey] = set()
        for st in structs:
            forced_res = find_substrate_by_idspec(st, args.selected_resn)
            forced_union |= {_residue_key_from_res(r) for r in forced_res}
        if forced_union:
            LOGGER.info("[extract:multi] Force-include (--selected-resn): +%d residues.", len(forced_union))
            union_sel_keys |= forced_union

    # 2) Disulfide partners (OR across structures)
    dis_keys_union: Set[ResidueKey] = set()
    for st in structs:
        dis_keys_union |= _disulfide_partner_keys(st, union_sel_keys, DISULFIDE_CUTOFF)
    if dis_keys_union:
        LOGGER.info("[extract:multi] Disulfide partner addition (union): +%d residues.", len(dis_keys_union))
    union_sel_keys |= dis_keys_union

    # 3) Backbone-contact neighbor augmentation (if exclude_backbone == False)
    keep_ncap_union: Set[ResidueKey] = set()
    keep_ccap_union: Set[ResidueKey] = set()
    if not args.exclude_backbone and union_bb_contact_keys:
        added_neighbor_union: Set[ResidueKey] = set()
        for st, subs in zip(structs, subs_per_struct):
            sel_ids = _keys_to_fids(st, union_sel_keys)
            bb_ids = _keys_to_fids(st, union_bb_contact_keys & _fids_to_keys(st, sel_ids))
            sub_ids = {r.get_full_id() for r in subs}
            # single call performs neighbor augmentation and returns cap-preservation flags
            kn_fids, kc_fids = augment_backbone_contact_neighbors(st, sel_ids, bb_ids, sub_ids)
            after_keys = _fids_to_keys(st, sel_ids)
            added_neighbor_union |= (after_keys - union_sel_keys)
            keep_ncap_union |= _fids_to_keys(st, kn_fids)
            keep_ccap_union |= _fids_to_keys(st, kc_fids)
        if added_neighbor_union:
            LOGGER.info("[extract:multi] Backbone-contact neighbor addition (union): +%d residues.",
                         len(added_neighbor_union))
        union_sel_keys |= added_neighbor_union

    # 4) PRO N-side neighbor augmentation (OR across structures)
    pro_prev_add_union: Set[ResidueKey] = set()
    for st in structs:
        sel_ids = _keys_to_fids(st, union_sel_keys)
        augment_proline_prev_neighbor(st, sel_ids)
        added = _fids_to_keys(st, sel_ids) - union_sel_keys
        pro_prev_add_union |= added
    if pro_prev_add_union:
        LOGGER.info("[extract:multi] PRO N-side neighbor addition (union): +%d residues.",
                     len(pro_prev_add_union))
    union_sel_keys |= pro_prev_add_union

    # ==== Build skip maps per structure (using unified selection and cap-keep flags) ====
    selected_ids_per_struct: List[Set[Tuple]] = []
    skip_maps_per_struct: List[Dict[Tuple, Set[str]]] = []
    substrate_idsets_per_struct: List[Set[Tuple]] = []

    for st, subs in zip(structs, subs_per_struct):
        sel_fids = _keys_to_fids(st, union_sel_keys)
        selected_ids_per_struct.append(sel_fids)
        sub_ids = {r.get_full_id() for r in subs}
        substrate_idsets_per_struct.append(sub_ids)
        kn_fids = _keys_to_fids(st, keep_ncap_union) if (not args.exclude_backbone) else None
        kc_fids = _keys_to_fids(st, keep_ccap_union) if (not args.exclude_backbone) else None
        skip_map = mark_atoms_to_skip(st, sel_fids, sub_ids, args.exclude_backbone, kn_fids, kc_fids)
        skip_maps_per_struct.append(skip_map)

    # ==== Compute link‑H definitions for each model and ensure identical targets/order ====
    linkdefs_per_struct: List[List[Tuple[Tuple[ResidueKey, str], Tuple[float, float, float]]]] = []
    for st, sel_fids, skip_map in zip(structs, selected_ids_per_struct, skip_maps_per_struct):
        linkdefs = _compute_linkH_defs(st, sel_fids, skip_map)
        linkdefs_per_struct.append(linkdefs)
    ref_targets = [ld[0] for ld in linkdefs_per_struct[0]]
    for i in range(1, len(linkdefs_per_struct)):
        targets_i = [ld[0] for ld in linkdefs_per_struct[i]]
        if targets_i != ref_targets:
            raise RuntimeError(
                f"[multi] link-H targets/order differ between model #1 and model #{i+1}. "
                f"Ensure inputs and options produce identical truncation across models."
            )
    LOGGER.info("[extract:multi] link-H targets common across models: %d.", len(ref_targets))

    # ==== Write outputs ====
    per_file_outputs = (len(args.output_pdb) == len(paths))
    if not per_file_outputs and len(args.output_pdb) != 1:
        raise ValueError("[extract:multi] Provide either a single output path for a multi‑MODEL PDB "
                         "or exactly N output paths where N == number of inputs for per‑structure outputs.")

    io = PDB.PDBIO()
    model_texts: List[str] = []
    model_counts: List[Dict[str, int]] = []

    for m, (st, sel_fids, skip_map) in enumerate(zip(structs, selected_ids_per_struct, skip_maps_per_struct), start=1):
        io.set_structure(st)
        buf = _io.StringIO()
        io.save(buf, AS_Select(sel_fids, skip_map))
        main_text = _strip_trailing_END(buf.getvalue())

        # Atom-count diagnostics
        raw_atoms = sum(len(st[f[1]][f[2]].child_dict[f[3]]) for f in sel_fids)
        kept_atoms = sum(
            1 for fid in sel_fids
            for a in st[fid[1]][fid[2]].child_dict[fid[3]]
            if a.get_name() not in skip_map.get(fid, set())
        )
        LOGGER.info("[extract:multi] Raw atoms (model %d): %d", m, raw_atoms)
        LOGGER.info("[extract:multi] Atoms after truncation (model %d): %d", m, kept_atoms)
        model_counts.append({"raw_atoms": raw_atoms, "kept_atoms": kept_atoms})

        # Append TER + link‑H block (honor --add-linkh)
        link_coords = [coord for (_, coord) in linkdefs_per_struct[m-1]]
        if args.add_linkh and link_coords:
            if not main_text.endswith("\n"):
                main_text += "\n"
            parts = [main_text]
            last_line = main_text.splitlines()[-1].strip() if main_text.strip() else ""
            if last_line != "TER":
                parts.append("TER\n")
            start_serial = _max_serial_from_pdb_text(main_text)
            parts.append(_format_linkH_block(link_coords, start_serial))
            main_text = "".join(parts)

        model_texts.append(main_text)

    outputs: List[str] = []
    if per_file_outputs:
        for idx, text in enumerate(model_texts):
            content = text
            if not content.endswith("\n"):
                content += "\n"
            content += "END\n"
            out_path = args.output_pdb[idx]
            with open(out_path, "w") as fh:
                fh.write(content)
            outputs.append(out_path)
            LOGGER.info("[extract:multi] Single‑model pocket saved to %s", out_path)
    else:
        buf_models: List[str] = []
        for m, text in enumerate(model_texts, start=1):
            model_block = []
            model_block.append(f"MODEL     {m}\n")
            model_block.append(text)
            model_block.append("ENDMDL\n")
            buf_models.append("".join(model_block))
        out_path = args.output_pdb[0]
        with open(out_path, "w") as fh:
            for blk in buf_models:
                fh.write(blk)
            fh.write("END\n")
        outputs.append(out_path)
        LOGGER.info("[extract:multi] Multi‑MODEL pocket saved to %s", out_path)

    # ==== Charge summary (first model only) ====
    charge_summary = compute_charge_summary(
        structs[0],
        selected_ids_per_struct[0],
        substrate_idsets_per_struct[0],
        getattr(args, "ligand_charge", None)
    )
    log_charge_summary("[extract:multi]", charge_summary)

    if api==True:
        return {
            "outputs": outputs,
            "counts": model_counts,
            "charge_summary": charge_summary,
        }
    else:
        return


# ---------------------------------------------------------------------
#   PDB writer helper
# ---------------------------------------------------------------------
class AS_Select(PDB.Select):
    """
    Biopython Select subclass that filters residues/atoms according to skip map.
    """
    def __init__(self, selected_ids: Set[Tuple], skip_map: Dict[Tuple, Set[str]]):
        self.ids = selected_ids
        self.skip = skip_map

    def accept_residue(self, residue):
        return residue.get_full_id() in self.ids

    def accept_atom(self, atom):
        fid = atom.get_parent().get_full_id()
        return atom.get_name() not in self.skip.get(fid, set())


# ---------------------------------------------------------------------
#   Main driver (single or multi) — CLI or API
# ---------------------------------------------------------------------

def extract(args: argparse.Namespace | None = None, api=False) -> Dict[str, Any]:
    """
    Run from CLI (args=None → parse_args()) or as an API with a pre-built Namespace.

    Args
    ----
    args : argparse.Namespace | None
        If None, parse CLI args. Otherwise, use the provided Namespace.
    api : bool
        If True, return a structured result dictionary; if False (CLI), return None.

    Returns
    -------
    dict | None
        When api=True, returns { 'outputs', 'counts', 'charge_summary' }. Otherwise, None.
    """
    if args is None:
        args = parse_args()

    _configure_extract_logger(bool(args.verbose))

    # Augment AMINO_ACIDS with user-specified modified residues
    # Save original state so repeated API calls don't accumulate mutations.
    _amino_acids_snapshot = dict(AMINO_ACIDS)
    try:
        return _extract_body(args, api)
    finally:
        AMINO_ACIDS.clear()
        AMINO_ACIDS.update(_amino_acids_snapshot)


def _extract_body(args, api):
    """Inner body of extract(), separated for try/finally AMINO_ACIDS restoration."""
    _mod_res = getattr(args, 'modified_residue', '') or ''
    if _mod_res:
        for token in _mod_res.replace(' ', ',').split(','):
            token = token.strip()
            if not token:
                continue
            if ':' in token:
                name, charge_str = token.split(':', 1)
                AMINO_ACIDS[name.strip().upper()] = int(float(charge_str.strip()))
            else:
                AMINO_ACIDS[token.upper()] = 0
        LOGGER.info("[extract] Modified residues added to amino acid list: %s", _mod_res)

    if args.radius == 0.0:
        args.radius = 0.001
    if args.radius_het2het == 0.0:
        args.radius_het2het = 0.001

    # Log extract options
    LOGGER.info("[extract] Options: radius=%.2f, radius_het2het=%.2f, "
                "include_h2o=%s, exclude_backbone=%s, add_linkh=%s, "
                "selected_resn='%s'",
                args.radius, args.radius_het2het,
                args.include_h2o, args.exclude_backbone,
                args.add_linkh,
                getattr(args, 'selected_resn', ''))

    # default output names
    if args.output_pdb is None:
        if len(args.complex_pdb) > 1:
            # multiple inputs → per-file outputs: pocket_{original_filename}.pdb
            args.output_pdb = [
                f"pocket_{os.path.splitext(os.path.basename(p))[0]}.pdb"
                for p in args.complex_pdb
            ]
        else:
            args.output_pdb = ['pocket.pdb']

    # Single-structure path
    if len(args.complex_pdb) == 1:
        complex_struct = load_structure(args.complex_pdb[0], "complex")

        # Resolve substrate residues from PDB path or residue-ID/name list
        substrate_residues = resolve_substrate_residues(complex_struct, args.substrate_pdb)
        substrate_ids = {r.get_full_id() for r in substrate_residues}
        LOGGER.info("[extract] Substrate residues matched: resseq %s",
                     [r.id[1] for r in substrate_residues])

        selected_ids, backbone_contact_ids = select_residues(
            complex_struct, substrate_residues,
            args.radius, args.radius_het2het,
            args.include_h2o,
            args.exclude_backbone
        )

        # Force-include residues via --selected-resn
        if getattr(args, "selected_resn", ""):
            forced_res = find_substrate_by_idspec(complex_struct, args.selected_resn)
            add_n = 0
            for r in forced_res:
                fid = r.get_full_id()
                if fid not in selected_ids:
                    selected_ids.add(fid)
                    add_n += 1
            if add_n:
                LOGGER.info("[extract] Force-include (--selected-resn): +%d residues.", add_n)

        augment_disulfides(complex_struct, selected_ids)

        # Backbone-contact context (if enabled)
        keep_ncap_ids: Set[Tuple] = set()
        keep_ccap_ids: Set[Tuple] = set()
        if not args.exclude_backbone and backbone_contact_ids:
            kn, kc = augment_backbone_contact_neighbors(
                complex_struct, selected_ids, backbone_contact_ids, substrate_ids
            )
            keep_ncap_ids.update(kn)
            keep_ccap_ids.update(kc)

        # Ensure PRO's N-side neighbor is included (TER-aware)
        augment_proline_prev_neighbor(complex_struct, selected_ids)

        # Atom counts
        raw = sum(len(complex_struct[f[1]][f[2]].child_dict[f[3]]) for f in selected_ids)
        LOGGER.info("[extract] Raw atoms: %d", raw)

        skip_map = mark_atoms_to_skip(
            complex_struct, selected_ids, substrate_ids,
            args.exclude_backbone,
            keep_ncap_ids if not args.exclude_backbone else None,
            keep_ccap_ids if not args.exclude_backbone else None
        )

        kept_atoms = sum(
            1 for fid in selected_ids
            for a in complex_struct[fid[1]][fid[2]].child_dict[fid[3]]
            if a.get_name() not in skip_map.get(fid, set())
        )
        LOGGER.info("[extract] Atoms after truncation: %d", kept_atoms)

        # Warn about non-standard residues that may be amino acids
        _bb_full = {"N", "CA", "C", "O"}
        for fid in selected_ids:
            res = complex_struct[fid[1]][fid[2]].child_dict[fid[3]]
            resname = res.get_resname()
            if resname in AMINO_ACIDS or resname in WATER_RES:
                continue
            if fid in substrate_ids:
                continue
            res_atoms = {a.get_name() for a in res}
            if _bb_full.issubset(res_atoms):
                _resseq = res.get_id()[1]
                LOGGER.warning(
                    "[extract] WARNING: Residue %s %d may be an amino acid "
                    "(has N, CA, C, O) but is not recognized as a standard residue name. "
                    "Backbone truncation was not applied. "
                    "Consider preparing the pocket model manually.",
                    resname, _resseq,
                )

        # Save structure (and optionally append link‑H block)
        io = PDB.PDBIO()
        io.set_structure(complex_struct)

        buf = _io.StringIO()
        io.save(buf, AS_Select(selected_ids, skip_map))
        main_pdb_text = buf.getvalue()

        output_path = args.output_pdb[0]
        outputs: List[str] = []

        if args.add_linkh:
            link_coords = compute_linkH_atoms(complex_struct, selected_ids, skip_map)
            LOGGER.info("[extract] Link-H to add: %d", len(link_coords))

            lines = [ln for ln in main_pdb_text.splitlines() if ln.strip() != "END"]
            if lines and lines[-1].strip() == "TER":
                pass
            main_no_end = "\n".join(lines)
            if not main_no_end.endswith("\n"):
                main_no_end += "\n"

            final_parts = [main_no_end]
            if link_coords:
                final_parts.append("TER\n")
                start_serial = _max_serial_from_pdb_text(main_no_end)
                final_parts.append(_format_linkH_block(link_coords, start_serial))
            final_parts.append("END\n")

            with open(output_path, "w") as fh:
                fh.write("".join(final_parts))
            LOGGER.info("[extract] Binding-Pocket (Active Site) + link-H saved to %s", output_path)
            outputs.append(output_path)
        else:
            with open(output_path, "w") as fh:
                fh.write(main_pdb_text)
            LOGGER.info("[extract] Binding-Pocket (Active Site) saved to %s", output_path)
            outputs.append(output_path)

        # Charge summary (single model)
        charge_summary = compute_charge_summary(
            complex_struct, selected_ids, substrate_ids, getattr(args, "ligand_charge", None)
        )
        log_charge_summary("[extract]", charge_summary)
        
        if api:
            return {
                "outputs": outputs,
                "counts": [{"raw_atoms": raw, "kept_atoms": kept_atoms}],
                "charge_summary": charge_summary,
            }
        else:
            return

    # Multi-structure path
    return extract_multi(args, api=api)


def extract_api(complex_pdb: List[str],
                   center: str,
                   output: Optional[List[str]] = None,
                   radius: float = 2.6,
                   radius_het2het: float = 0.0,
                   include_h2o: bool = True,
                   exclude_backbone: bool = False,
                   add_linkh: bool = False,
                   selected_resn: str = "",
                   modified_residue: str = "",
                   ligand_charge: Optional[float | str | Dict[str, float]] = None,
                   verbose: bool = False) -> Dict[str, Any]:
    """
    Convenience API for programmatic use.

    Args
    ----
    complex_pdb : list[str]
        Input PDB path(s). len==1 → single, len>1 → multi.
    center : str
        Substrate spec: a PDB path, a residue‑ID list 'A:123,456' (insertion codes OK),
        or a residue‑name list 'GPP,MMT'.
    output : list[str] | None
        Output path(s): one path for multi‑MODEL PDB, or N paths for per‑file outputs.
        If None, defaults to ['pocket.pdb'].
    radius : float
        Atom–atom cutoff (Å) for inclusion around substrate atoms.
    radius_het2het : float
        Independent hetero‑hetero cutoff (Å) for non‑C/H pairs.
    include_h2o : bool
        Include waters in the selection.
    exclude_backbone : bool
        Remove backbone atoms on non‑substrate amino acids (with safeguards).
    add_linkh : bool
        Add link‑H atoms for cut bonds (carbon‑only) and append as HL/LKH HETATM records.
    selected_resn : str
        Additional residues to force‑include (comma/space separated).
    modified_residue : str
        Comma‑separated residue names (with optional charge) to treat as amino acids
        for backbone truncation and charge assignment. E.g. 'HD1,HD2' or 'HD1:0,SEP:-2'.
    ligand_charge : float | str | dict[str,float] | None
        Either a total charge (float/str) for unknown residues (prefer unknown substrate),
        or a mapping like {'GPP': -3, 'MMT': -1}. In mapping mode, other unknown residues remain 0.
    verbose : bool
        Enable INFO logging.

    Returns
    -------
    dict
        Same structure as `extract(..., api=True)`.
    """
    if not output:
        output = ['pocket.pdb']
    ns = argparse.Namespace(
        complex_pdb=complex_pdb,
        substrate_pdb=center,
        output_pdb=output,
        radius=radius,
        radius_het2het=radius_het2het,
        include_h2o=include_h2o,
        exclude_backbone=exclude_backbone,
        add_linkh=add_linkh,
        selected_resn=selected_resn,
        modified_residue=modified_residue,
        ligand_charge=ligand_charge,
        verbose=verbose,
    )
    return extract(ns, api=True)
