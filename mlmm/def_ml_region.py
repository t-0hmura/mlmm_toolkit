"""
def_ml_region.py — Automated ML‐region extractor for protein–substrate complexes
============================================================================

This script builds a Machine-Learning (ML) region around one or more substrate
residues in a protein–substrate complex, ready for ML/MM calculations.

Key features
------------
* **Substrate identification**  
  Every residue contained in <substrate_pdb> is matched *exactly* (all atomic
  coordinates within ±1 × 10⁻³ Å) to residues in <complex_pdb>.  
  The script aborts if any substrate residue cannot be matched unambiguously.

* **Residue selection around the substrate**  
  1. Any residue that has *any* atom within **--radius_ml** Å of a substrate atom
     is selected.  
  2. Additionally, hetero-atoms (non-C/H) of any residue that fall within
     **--radius_het** Å of a substrate hetero-atom pull in their parent residues.  
  3. Waters are included only when `--include_H2O true` is given.  
  4. Sequence neighbors (±1 in residue index) that were pulled in *only* via
     chain continuity are removed again unless they also satisfy the cut-off
     criteria above.

* **Side-chain/backbone truncation rules**  
  * **Isolated residues** (segment length = 1) are trimmed to *pure side-chain*:
    backbone atoms {N, CA, C, O, OXT, H/H1-3, HN, HA/HA2/HA3} are deleted  
    (exception – PRO never loses N, CA, HA* so the five-membered ring stays
    intact).  
  * **Continuous peptide stretches** retain their internal backbone; only
    terminal caps are removed (N-terminal: N/H*, C-terminal: C/O/OXT).  
  * **Cys–Cys disulfides**: if either cysteine of a C–Sγ···Sγ–C pair (≤ 2.5 Å)
    is in the ML region, its partner cysteine is added automatically.
  * **`--exclude_backbone true`**: for *all* amino-acid residues **except the
    substrate**, the entire main-chain scaffold (N, H*, CA, HA*, C, O) is removed,
    leaving only side-chain atoms.

* **Output**  
  A new PDB (<output_pdb>) that contains only the selected residues, already
  truncated according to the rules above and suitable for downstream ML/MM
  workflows.

Command-line usage
------------------
    def_ml_region -r complex.pdb -c substrate.pdb -o ml_region.pdb \
        [--radius_ml 2.6] [--radius_het 3.6] [--include_H2O true|false] \
        [--exclude_backbone true|false]

Example
    def_ml_region -r complex.pdb -c ligand.pdb -o ml_region.pdb --radius_ml 6.0 --radius_het 6.0 --include_H2O false --exclude_backbone true

Arguments
~~~~~~~~~
-r complex_pdb      Protein–substrate complex PDB file (all atoms).
-c substrate_pdb    PDB containing *exactly* the substrate residue(s) of interest.
-o output_pdb       Path where the trimmed ML-region PDB will be written.

Options
~~~~~~~
--radius_ml R       Distance (Å) from any substrate atom for residue inclusion
                    (default: 2.6 Å).
--radius_het R      Distance (Å) from substrate hetero-atoms to other hetero-atoms
                    for additional residue inclusion (default: 3.6 Å).
--include_H2O bool  Include crystallographic/explicit waters (default: false).
--exclude_backbone bool Remove all main-chain atoms (N, H, CA, HA, C, O) from every
                    amino-acid residue except the substrate (default: true).
                    For PRO, the N and CA atoms are retained.
                    If false, only terminal N and C caps are removed.
Dependencies
------------
* Python ≥ 3.9
* Biopython ≥ 1.80 (for `Bio.PDB`)
* NumPy (indirectly via Biopython)
"""

from __future__ import annotations
import argparse
from typing import Dict, List, Set, Tuple

from Bio import PDB
from Bio.PDB import NeighborSearch

# ---------------------------------------------------------------------
#   Constants
# ---------------------------------------------------------------------
BACKBONE_ATOMS: Set[str] = {
    "N", "C", "O", "CA", "OXT",
    "H", "H1", "H2", "H3", "HN", "HA", "HA2", "HA3",
}
# When --exclude_backbone true, remove the full main-chain set:
BACKBONE_ALL: Set[str] = BACKBONE_ATOMS

AMINO_ACIDS: Set[str] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL",

    # Non-standard amino acids (e.g., modified residues)
    "SEP", "TPO", "HIP", "PTR", "ASA", "CIR", "FOR", "MVA", "IIL", "AIB",
    "HTN", "CSO", "CSD", "CSX", "OCS", "SEC", "KCX", "MLY", "LLP", "PYL",
    "DLY", "HYP", "DPR", "CGA", "DGL", "PCA", "SAR", "NMC", "PFF", "NFA",
    "MSE", "OMT", "ASH", "HID", "HIE", "CYX"

}

DISULFIDE_CUTOFF = 2.5   # Å Sγ–Sγ
EXACT_EPS = 1e-3          # Å tolerance for exact match
WATER_NAMES = {"HOH", "WAT", "TIP3"}

# ---------------------------------------------------------------------
#   Helpers
# ---------------------------------------------------------------------

def str2bool(v: str) -> bool:
    if isinstance(v, bool):
        return v
    return v.lower() in {"true", "1", "yes", "y"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Extract an ML region around one or more substrate residues.")

    p.add_argument("-r", "--real",   dest="complex_pdb",   required=True,
                   metavar="complex.pdb",
                   help="Protein–substrate complex PDB file")
    p.add_argument("-c", "--center", dest="substrate_pdb", required=True,
                   metavar="ligand.pdb",
                   help="Substrate/ligand PDB file (one or more residues)")
    p.add_argument("-o", "--out",  dest="output_pdb",    required=True,
                   metavar="ml_region.pdb",
                   help="Output PDB file for the extracted ML region")

    p.add_argument("--radius_ml", type=float, default=2.6,
                   help="Residue cut-off in Å from any substrate atom "
                        "(default: 2.6)")
    p.add_argument("--radius_het", type=float, default=3.6,
                   help="Hetero-hetero cut-off in Å "
                        "(default: 3.6)")
    p.add_argument("--include_H2O", type=str2bool, default=False,
                   help="Include water molecules (default: false)")
    p.add_argument("--exclude_backbone", type=str2bool, default=True,
                   help="Delete all main-chain atoms from non-substrate residues "
                        "(default: true)")
    return p.parse_args()



def load_structure(path: str, name: str) -> PDB.Structure.Structure:
    parser = PDB.PDBParser(QUIET=True)
    return parser.get_structure(name, path)


# ---------------------------------------------------------------------
#   Substrate matching
# ---------------------------------------------------------------------

def is_exact_match(lig_atoms: Dict[str, PDB.Vector.Vector],
                   cand: PDB.Residue.Residue) -> bool:
    for name, vec in lig_atoms.items():
        if name not in cand:
            return False
        if (vec - cand[name].get_vector()).norm() > EXACT_EPS:
            return False
    return True


def find_substrate_residues(complex_struct, substrate_struct) -> List[PDB.Residue.Residue]:
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
            raise ValueError(
                f"Exact match not found for substrate residue {lig_name} {lig.id[1]}")
    return matched


# ---------------------------------------------------------------------
#   Residue selection around the substrate
# ---------------------------------------------------------------------

def select_residues(complex_struct, substrate_res_list: List[PDB.Residue.Residue],
                    r_ml: float, r_het: float, include_h2o: bool) -> Set[Tuple]:
    substrate_atoms = [a for lig in substrate_res_list for a in lig]
    substrate_het = [a for a in substrate_atoms if a.element not in ("C", "H")]
    ns = NeighborSearch(list(complex_struct.get_atoms()))

    selected_ids: Set[Tuple] = {res.get_full_id() for res in substrate_res_list}
    within_cutoff: Set[Tuple] = set(selected_ids)

    def maybe_add(atom):
        res = atom.get_parent()
        if not include_h2o and res.get_resname() in WATER_NAMES:
            return
        fid = res.get_full_id()
        selected_ids.add(fid)
        within_cutoff.add(fid)

    for atom in substrate_atoms:
        for neigh in ns.search(atom.get_coord(), r_ml):
            maybe_add(neigh)
    for atom in substrate_het:
        for neigh in ns.search(atom.get_coord(), r_het):
            if neigh.element in ("C", "H"):
                continue
            maybe_add(neigh)

    # remove sequence neighbors that were pulled in only via chain continuity
    for lig in substrate_res_list:
        model, chain, resinfo = lig.get_full_id()[1:4]
        seq = resinfo[1]
        chain_ent = complex_struct[model][chain]
        for off in (-1, 1):
            nbr = chain_ent.child_dict.get((' ', seq + off, ' '))
            if nbr:
                nbr_id = nbr.get_full_id()
                if nbr_id in selected_ids and nbr_id not in within_cutoff:
                    selected_ids.remove(nbr_id)

    return selected_ids


# ---------------------------------------------------------------------
#   Disulfide augmentation
# ---------------------------------------------------------------------

def augment_disulfides(structure, selected_ids: Set[Tuple],
                       cutoff: float = DISULFIDE_CUTOFF):
    sg_atoms = [r["SG"] for r in structure.get_residues()
                if r.get_resname() == "CYS" and "SG" in r]

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
#   Backbone trimming / skip-map generation
# ---------------------------------------------------------------------

def continuous_segments(sorted_ids: List[Tuple]):
    segs, cur, prev = [], [], None
    for fid in sorted_ids:
        idx = fid[3][1]
        if prev is None or idx == prev + 1:
            cur.append(fid)
        else:
            segs.append(cur)
            cur = [fid]
        prev = idx
    if cur:
        segs.append(cur)
    return segs


def mark_atoms_to_skip(structure, selected_ids: Set[Tuple], substrate_ids: Set[Tuple],
                       exclude_backbone: bool) -> Dict[Tuple, Set[str]]:
    """Return a mapping full-id → atoms to delete."""

    # start with the original truncation logic (except for substrate residues)
    chain_map: Dict[Tuple[str, str], List[Tuple]] = {}
    for fid in selected_ids:
        if fid in substrate_ids:
            continue  # never delete atoms from substrate residues
        res = structure[fid[1]][fid[2]].child_dict[fid[3]]
        if res.get_resname() in WATER_NAMES:
            continue
        chain_map.setdefault((fid[1], fid[2]), []).append(fid)

    skip: Dict[Tuple, Set[str]] = {}

    for (model, chain), fids in chain_map.items():
        fids.sort(key=lambda x: x[3][1])
        for seg in continuous_segments(fids):
            n_id, c_id = seg[0], seg[-1]
            single = len(seg) == 1

            def add(fid, names):
                skip.setdefault(fid, set()).update(names)

            n_res = structure[model][chain].child_dict[n_id[3]]
            c_res = structure[model][chain].child_dict[c_id[3]]

            # N-terminal cap deletion
            if n_res.get_resname() != "PRO":
                add(n_id, {"N", "H", "H1", "H2", "H3", "HN"})
            # C-terminal cap deletion
            add(c_id, {"C", "O", "OXT"})

            # Isolated stretch – remove CA/HA* (except PRO)
            if single and n_res.get_resname() != "PRO":
                add(n_id, {"CA", "HA", "HA2", "HA3"})

    # ---------------------------------------------------------------------
    #   Optional: remove *all* backbone atoms from every non-substrate residue
    # ---------------------------------------------------------------------
    if exclude_backbone:
        for fid in selected_ids:
            if fid in substrate_ids:
                continue
            res = structure[fid[1]][fid[2]].child_dict[fid[3]]
            if res.get_resname() in WATER_NAMES:
                continue
            if res.get_resname() in AMINO_ACIDS:
                skip.setdefault(fid, set()).update(BACKBONE_ALL)

    return skip


# ---------------------------------------------------------------------
#   PDB writer helper
# ---------------------------------------------------------------------
class MLSelect(PDB.Select):
    def __init__(self, selected_ids: Set[Tuple], skip_map: Dict[Tuple, Set[str]]):
        self.ids = selected_ids
        self.skip = skip_map

    def accept_residue(self, residue):
        return residue.get_full_id() in self.ids

    def accept_atom(self, atom):
        fid = atom.get_parent().get_full_id()
        return atom.get_name() not in self.skip.get(fid, set())


# ---------------------------------------------------------------------
#   Main driver
# ---------------------------------------------------------------------

def def_ml_region():
    args = parse_args()

    complex_struct = load_structure(args.complex_pdb, "complex")
    substrate_struct = load_structure(args.substrate_pdb, "substrate")

    substrate_residues = find_substrate_residues(complex_struct, substrate_struct)
    substrate_ids = {r.get_full_id() for r in substrate_residues}
    print("[def_ml_region] Substrate residues matched:", [r.id[1] for r in substrate_residues])

    selected_ids = select_residues(complex_struct, substrate_residues,
                                   args.radius_ml, args.radius_het,
                                   args.include_H2O)
    # count raw atoms
    raw = sum(len(complex_struct[f[1]][f[2]].child_dict[f[3]]) for f in selected_ids)
    print(f"[def_ml_region] Raw atoms: {raw}")

    augment_disulfides(complex_struct, selected_ids)

    skip_map = mark_atoms_to_skip(complex_struct, selected_ids, substrate_ids,
                                 args.exclude_backbone)

    kept_atoms = sum(
        1 for fid in selected_ids
        for a in complex_struct[fid[1]][fid[2]].child_dict[fid[3]]
        if a.get_name() not in skip_map.get(fid, set())
    )
    print(f"[def_ml_region] Atoms after truncation: {kept_atoms}")

    io = PDB.PDBIO()
    io.set_structure(complex_struct)
    io.save(args.output_pdb, MLSelect(selected_ids, skip_map))
    print(f"[def_ml_region] ML region saved to {args.output_pdb}")


if __name__ == "__main__":
    def_ml_region()
