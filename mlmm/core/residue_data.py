"""Residue / ion / water dictionaries shared by L2 (workflows) and L3 (domain).

Foundation-tier data — no MLIP / numpy / torch / Bio dependencies, and no
Click / workflow imports.  ``mlmm.workflows.extract`` (charge inference),
``mlmm.domain.add_elem_info`` (element inference), and ``mlmm.workflows.mm_parm``
all consume these tables, so this is the single source of truth: editing one
table can no longer change element inference without changing charge inference.

The public tables are exported **read-only** (``MappingProxyType`` / ``frozenset``)
so a consumer cannot mutate the canonical catalog.  Callers that must add
entries (e.g. ``extract``'s ``--modified-residue``) copy into a request-local
``dict`` at their single mutation boundary; the canonical table stays intact.

For the amino-acid keys, ``rn = res.get_resname().upper()`` is used at the
charge-inference lookup, so uppercase keys are what the protein/terminus path
reads.  The ``ION`` table additionally carries Amber-style mixed-case aliases
for the raw-name path (``mm_parm``); their uppercase folds equal the uppercase
keys, so the upper-casing consumers see an identical table.
"""

from __future__ import annotations

from types import MappingProxyType
from typing import Dict, Mapping, Tuple

_AMINO_ACIDS: Dict[str, int] = {
    # --- Standard 20 (L) ---
    "ALA":  0, "ARG": +1, "ASN":  0, "ASP": -1, "CYS":  0,
    "GLU": -1, "GLN":  0, "GLY":  0, "HIS":  0, "ILE":  0,
    "LEU":  0, "LYS": +1, "MET":  0, "PHE":  0, "PRO":  0,
    "SER":  0, "THR":  0, "TRP":  0, "TYR":  0, "VAL":  0,

    # --- Canonical extras ---
    "SEC":  0,   # selenocysteine
    "PYL":  0,   # pyrrolysine

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
    "MLY": +1, "LLP":  0,
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

_ION: Dict[str, int] = {
    # +1
    "LI": +1, "NA": +1, "K": +1, "RB": +1, "CS": +1, "TL": +1, "AG": +1, "CU1": +1,
    "K+": +1, "NA+": +1, "NH4": +1, "H3O+": +1, "H3O": +1, "HE+": +1, "HZ+": +1,

    # +2
    "MG": +2, "CA": +2, "SR": +2, "BA": +2, "MN": +2, "FE2": +2, "CO": +2, "NI": +2,
    "CU": +2, "ZN": +2, "CD": +2, "HG": +2, "PB": +2, "BE": +2, "PD": +2, "PT": +2,
    "SN": +2, "RA": +2, "YB2": +2, "V2+": +2,

    # +3
    "FE": +3, "AU3": +3, "AL": +3, "GA": +3, "IN": +3,
    "CE": +3, "CR": +3, "DY": +3, "EU": +3, "EU3": +3, "ER": +3,
    "GD3": +3, "LA": +3, "LU": +3, "ND": +3, "PR": +3, "SM": +3, "TB": +3,
    "TM": +3, "Y": +3, "PU": +3,

    # +4
    "U4+": +4, "TH": +4, "HF": +4, "ZR": +4,

    # -1
    "F": -1, "CL": -1, "BR": -1, "I": -1, "CL-": -1, "IOD": -1,

    # Amber-style mixed-case aliases: some consumers (e.g. mm_parm) look up the
    # raw PDB residue name without upper-casing, and Amber writes ions like "Na+"
    # "Cl-", or "Be". Same charges as the upper-case keys above; kept so a single
    # canonical ION table serves both the upper-casing (extract) and raw-name paths.
    "Na+": +1, "Ag": +1, "Tl": +1,
    "Be": +2, "Ra": +2, "Sn": +2,
    "Ce": +3, "Cr": +3, "Dy": +3, "Er": +3, "Nd": +3, "Pu": +3, "Sm": +3, "Tm": +3,
    "Hf": +4, "Th": +4, "Zr": +4,
    "Cl-": -1,
}

# Amber terminal residue names whose AMINO_ACIDS value already BAKES IN the
# ionized-terminus formal charge (e.g. CGLU=-2, NLYS=+2). For these the
# kept-cap correction in compute_charge_summary must NOT be applied again.
C_TERMINAL_RESNAMES: frozenset = frozenset({
    "CALA", "CARG", "CASN", "CASP", "CCYS", "CCYX", "CGLN", "CGLU", "CGLY",
    "CHID", "CHIE", "CHIP", "CHYP", "CILE", "CLEU", "CLYS", "CMET", "CPHE",
    "CPRO", "CSER", "CTHR", "CTRP", "CTYR", "CVAL", "CTER",
})
N_TERMINAL_RESNAMES: frozenset = frozenset({
    "NALA", "NARG", "NASN", "NASP", "NCYS", "NCYX", "NGLN", "NGLU", "NGLY",
    "NHID", "NHIE", "NHIP", "NILE", "NLEU", "NLYS", "NMET", "NPHE", "NPRO",
    "NSER", "NTHR", "NTRP", "NTYR", "NVAL", "NTER",
})

WATER_RES: frozenset = frozenset({"HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SOL"})

DISULFIDE_CUTOFF = 2.5   # Å Sγ–Sγ (SG–SG)
EXACT_EPS = 1e-3         # Å tolerance for exact match

# Cross-structure residue identity: (chain, hetflag, resseq, icode, resname).
# The runtime key is built with ``str(resseq)`` (all fields are strings).
ResidueKey = Tuple[str, str, str, str, str]

# Public read-only views of the mutable source tables above.  Mutating either
# proxy raises ``TypeError``; consumers that must extend the catalog copy into a
# request-local ``dict`` first (see ``extract`` ``--modified-residue``).
AMINO_ACIDS: Mapping[str, int] = MappingProxyType(_AMINO_ACIDS)
ION: Mapping[str, int] = MappingProxyType(_ION)
