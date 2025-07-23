"""
energy_summary.py
------------------------------------------------------------
(1) CLI tool  : Generate a Gibbs free-energy profile for
                Reactant / Transition State / Product.
(2) API funcs : analyze_single_structure()
                analyze_three_structures()

Default settings
----------------
Temperature      : 298.15 K

Examples – command line
-----------------------
$ energy_summary reaction.yml

yaml config example:
geom:
  reac: ./reac.xyz
  ts:   ./ts.xyz
  prod: ./prod.xyz
  freeze_atoms: [0, 1, 2, 5, 6, 7] # Except from Vibration Analysis
calc:
  real_pdb:   ./parm/complex.pdb
  real_parm7: ./parm/complex.parm7
  real_rst7:  ./parm/complex.rst7
  model_pdb:  ./parm/ml_region.pdb
  model_charge: -1
  model_mult: 1
  backend: aimnet2
  ml_device: auto
  ml_cuda_idx: 0
  mm_device: cpu
  mm_cuda_idx: 0
  mm_threads: 16
  vib_run: true

Examples – Jupyter / script
---------------------------
>>> from mlmm import analyze_single_structure, analyze_three_structures, mlmm
>>> calc = mlmm(**mlmm_kwargs)

# --- single structure ------------------------------------------------
>>> data = analyze_single_structure(
...     "int.xyz", calc, temp=298.15, vib_cutoff=100,
...     freeze_atoms=[0, 3, 9])
>>> print(data["G_kcal"])

# --- three structures + table + plots + JSON ------------------------
>>> results = analyze_three_structures(
...     "reac.xyz", "ts.xyz", "prod.xyz", calc,
...     temp=298.15, vib_cutoff=100, out_stem="my_reaction",
...     freeze_atoms=[0, 1, 2])
>>> print(results["ts"]["G_kcal"])
"""

# ---------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import argparse, json, sys, time
from pathlib import Path
from typing import Sequence, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import torch, yaml
from ase.io import read
from ase.data import atomic_masses              # >>> TR-PROJ >>>
from pysisyphus.constants import (
    ANG2BOHR, AU2KCALPERMOL, KB, JOULE2EV, PLANCK, C, NA
)

from mlmm import calc_freq_from_hessian, mlmm

# ---------------------------------------------------------------------
#   Constants
# ---------------------------------------------------------------------
KB_eV   = KB * JOULE2EV
H2KCAL  = AU2KCALPERMOL
EV2KCAL = (1.0 / JOULE2EV) * NA / 4184
CM2EV   = PLANCK * C * 1e-2 * JOULE2EV 

# ---------------------------------------------------------------------
#   Vibrational corrections
# ---------------------------------------------------------------------
def vib_corrections(hnu_eV: np.ndarray, temp: float):
    """Return ZPE, (U-TS), and G_vib in eV for positive hν list."""
    if hnu_eV.size == 0:
        return 0.0, 0.0, 0.0
    beta = hnu_eV / (KB_eV * temp)
    zpe  = 0.5 * hnu_eV
    u_v  = hnu_eV / np.expm1(beta)
    s_v  = KB_eV * (beta / np.expm1(beta) - np.log1p(-np.exp(-beta)))
    g_v  = zpe + u_v - temp * s_v
    return zpe.sum(), (u_v - temp * s_v).sum(), g_v.sum()

# ---------------------------------------------------------------------
#   Calculator factory
# ---------------------------------------------------------------------
def build_calc(cfg: dict, *, freeze_atoms: Optional[Sequence[int]] = None):
    """Instantiate mlmm calculator with freeze_atoms injected."""
    cfg_m = dict(cfg)                         # shallow copy
    cfg_m.update(out_hess_torch=True)
    if freeze_atoms is not None and "freeze_atoms" not in cfg_m:
        cfg_m["freeze_atoms"] = list(freeze_atoms)

    if cfg_m.get("vib_run", True) is False:
        WARNING = "[WARNING] vib_run is False, no vibrational analysis will be performed on MM region."
        print(WARNING, file=sys.stderr)
        print(WARNING, file=sys.stdout)

    return mlmm(**cfg_m)

# ---------------------------------------------------------------------
#   Hessian zero-padding for frozen DOF
# ---------------------------------------------------------------------
def _zero_frozen_dof(
    H: torch.Tensor, freeze_atoms: Sequence[int] | None
) -> torch.Tensor:
    """Set rows & cols of frozen atoms (3× per atom) to zero in-place."""
    if not freeze_atoms:
        return H
    n_atoms = H.shape[0] // 3
    bad = [i for i in freeze_atoms if i < 0 or i >= n_atoms]
    if bad:
        raise ValueError(f"freeze_atoms contains invalid index {bad} (n_atoms={n_atoms})")

    # 1D DOF indices to zero out
    dof_idx: List[int] = []
    for i in freeze_atoms:
        base = 3 * i
        dof_idx.extend((base, base + 1, base + 2))

    H[dof_idx, :] = 0.0
    H[:, dof_idx] = 0.0
    return H

# ---------------------------------------------------------------------
#   Single-structure analysis
# ---------------------------------------------------------------------
def analyze_single_structure(
    path: str,
    calc,
    *,
    temp: float = 298.15,
    vib_cutoff: float = 100.0,
    freeze_atoms: Optional[Sequence[int]] = None,
):
    atoms = read(path)

    res = calc.get_hessian(
        atoms.get_chemical_symbols(),
        atoms.get_positions() * ANG2BOHR,
    )
    E_h = float(res["energy"])
    H   = res["hessian"]
    del res; torch.cuda.empty_cache()

    # --- ZERO-OUT frozen DOF ---------------------------------------- 
    H = _zero_frozen_dof(H, freeze_atoms)

    # >>> TR-PROJ >>>  (prepare for TR projection)
    coords_bohr = atoms.get_positions() * ANG2BOHR
    masses_amu  = np.array([atomic_masses[z] for z in atoms.get_atomic_numbers()])
    # <<< TR-PROJ <<<

    freqs_cm, hnu_eV, _ = calc_freq_from_hessian(
        H,
        atoms.get_atomic_numbers(),
        # >>> TR-PROJ >>>
        project_tr=True,
        coords_bohr=coords_bohr,
        masses_amu=masses_amu,
        # <<< TR-PROJ <<<
        verbose=False,
    )
    del H, _; torch.cuda.empty_cache()

    # Insurance against 0 frequencies
    zero_tol = 1e-6            # cm⁻¹
    zero_mask = np.abs(freqs_cm) < zero_tol
    if np.any(zero_mask):
        print(f"[INFO] {zero_mask.sum()} mode(s) |ν|<{zero_tol} cm⁻¹ "
            "removed from thermodynamic corrections.", file=sys.stderr)
        freqs_cm = freqs_cm[~zero_mask]
        hnu_eV   = hnu_eV[~zero_mask]

    # discard 1 imaginary mode and raise ultra-low real modes to vib_cutoff (Truhlar type QRRHO)
    neg_big_mask = (freqs_cm < 0.0) & (np.abs(freqs_cm) > vib_cutoff)
    n_neg_big    = int(neg_big_mask.sum())
    if n_neg_big == 1:
        keep_mask = ~neg_big_mask
        freqs_cm = freqs_cm[keep_mask]
        hnu_eV   = hnu_eV[keep_mask]
    elif n_neg_big > 1:
        print(f"[WARNING] {n_neg_big} imaginary modes |ν|>{vib_cutoff} cm⁻¹ "
              "detected – all ignored in thermodynamic corrections.",
              file=sys.stderr)
        keep_mask = ~neg_big_mask
        freqs_cm = freqs_cm[keep_mask]
        hnu_eV   = hnu_eV[keep_mask]
    low_mask = np.abs(freqs_cm) < vib_cutoff
    if np.any(low_mask):
        scale = vib_cutoff / np.abs(freqs_cm[low_mask])
        hnu_eV[low_mask] = np.abs(hnu_eV[low_mask]) * scale
        freqs_cm[low_mask] = vib_cutoff

    zpe_eV, vib_th_eV, g_corr_eV = vib_corrections(hnu_eV, temp)

    return {
        "E_hartree":    E_h,
        "E_kcal":       E_h * H2KCAL,
        "ZPE_kcal":     zpe_eV    * EV2KCAL,
        "vib_th_kcal":  vib_th_eV * EV2KCAL,
        "G_kcal":       (E_h * H2KCAL) + (g_corr_eV * EV2KCAL),
        "freqs_cm":     freqs_cm.tolist(),
    }

# ---------------------------------------------------------------------
#   Three-structure analysis + I/O
# ---------------------------------------------------------------------
def analyze_three_structures(
    reac_path: str,
    ts_path: str,
    prod_path: str,
    calc,
    *,
    temp: float = 298.15,
    vib_cutoff: float = 100.0,
    freeze_atoms: Sequence[int] | None = None,
    out_stem: str | None = None,
    show_table: bool = True,
    save_json: bool = False,
    annotate: bool = True,
    ):
    start_time = time.time()
    print(f"Thermal Analysis on structures: {reac_path}, {ts_path}, {prod_path} ...")
    results = {
        "reac": analyze_single_structure(
            reac_path, calc, temp=temp, vib_cutoff=vib_cutoff,
            freeze_atoms=freeze_atoms),
        "ts":   analyze_single_structure(
            ts_path,   calc, temp=temp, vib_cutoff=vib_cutoff,
            freeze_atoms=freeze_atoms),
        "prod": analyze_single_structure(
            prod_path, calc, temp=temp, vib_cutoff=vib_cutoff,
            freeze_atoms=freeze_atoms),
    }
    elapsed = time.time() - start_time
    minutes, seconds = divmod(elapsed, 60)
    print(f"Thermal Analysis completed in {int(minutes)} min {seconds:.2f} sec")

    if show_table:
        hdr = f"{'Structure':<10} {'E (kcal)':>12} {'ZPE':>8} {'G_vib':>8} {'G (kcal)':>12}"
        print(hdr); print("-"*len(hdr))
        for k, r in results.items():
            g_vib = r['ZPE_kcal'] + r['vib_th_kcal']
            print(f"{k:<10} {r['E_kcal']:12.2f} {r['ZPE_kcal']:8.2f} "
                  f"{g_vib:8.2f} {r['G_kcal']:12.2f}")
        dg_act = results["ts"]["G_kcal"]   - results["reac"]["G_kcal"]
        dg_rxn = results["prod"]["G_kcal"] - results["reac"]["G_kcal"]
        print(f"\nΔG‡ = {dg_act:.2f} kcal/mol\nΔG  = {dg_rxn:.2f} kcal/mol")

    if out_stem:
        start_time = time.time()
        make_plots(results, out_stem, annotate=annotate)
        elapsed = time.time() - start_time
        minutes, seconds = divmod(elapsed, 60)
        print(f"Plots saved as {out_stem}_E.png and {out_stem}_G.png "
              f"in {int(minutes)} min {seconds:.2f} sec")
        if save_json:
            Path(f"{out_stem}_results.json").write_text(json.dumps(results, indent=2))

    return results

# ---------------------------------------------------------------------
#   Plot helpers  (unchanged)
# ---------------------------------------------------------------------
def _draw_energy_diagram(ax, y_vals, color, typ, annotate=True):
    """Draw one R-TS-P diagram (equal bar lengths, labelled x-axis)."""
    ER, ETS, EP = y_vals
    # bar positions
    pos = [0.0, 0.6, 1.2]          # left edges
    w   = 0.3                      # bar width (all equal)
    mid = [p + w/2 for p in pos]   # mid-points for dashed ramps

    height = max(y_vals) - min(y_vals)

    # Reactant, TS, Product bars
    ax.hlines(ER,  pos[0], pos[0]+w, colors=color, lw=4)
    ax.hlines(ETS, pos[1], pos[1]+w, colors=color, lw=4)
    ax.hlines(EP,  pos[2], pos[2]+w, colors=color, lw=4)

    # connecting ramps
    ax.plot([pos[0]+w, pos[1]], [ER,  ETS], ls='--', color=color, lw=2)
    ax.plot([pos[1]+w, pos[2]], [ETS, EP ], ls='--', color=color, lw=2)

    if annotate:
        # annotations
        ax.annotate(rf"$\Delta {typ}^{{‡}} = {ETS: .1f}$",
                    (mid[1], ETS + height/25), ha='center', fontsize=13, color=color)
        ax.annotate(rf"$\Delta {typ} = {EP : .1f}$",
                    (mid[2]+0.12, EP + height/25), ha='center', fontsize=13, color=color)

def make_plots(results: dict, stem: str, annotate: bool = True):
    refE = results["reac"]["E_kcal"]
    refG = results["reac"]["G_kcal"]
    relE = (0,
            results["ts"]["E_kcal"]  - refE,
            results["prod"]["E_kcal"]- refE)
    relG = (0,
            results["ts"]["G_kcal"]  - refG,
            results["prod"]["G_kcal"]- refG)

    for rel, ylabel, fname, typ in [
        (relE, "ΔE (kcal/mol)", f"{stem}_E.png", "E"),
        (relG, "ΔG (kcal/mol)", f"{stem}_G.png", "G"),
    ]:
        fig, ax = plt.subplots(figsize=(6, 4.5))
        _draw_energy_diagram(ax, rel, color='black', typ=typ, annotate=annotate)

        # x-axis
        ax.set_xticks([0.15, 0.75, 1.35])
        ax.set_xticklabels(["Reactant", "TS", "Product"], fontsize=14)

        # y-axis with ± energy-range * 0.2 kcal/mol margin
        height = max(rel) - min(rel)
        ymin = min(rel) - 0.2 * height
        ymax = max(rel) + 0.2 * height
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel(ylabel, fontsize=15)

        ax.tick_params(axis='both', which='major', labelsize=13)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

        plt.tight_layout()
        fig.savefig(fname, dpi=300)
        plt.close(fig)


# ---------------------------------------------------------------------
#   CLI
# ---------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Gibbs free-energy profile tool (freeze_atoms compatible)")
    p.add_argument("config", help="YAML file with 'geom' and 'calc'")
    p.add_argument("--temp", type=float, default=298.15, help="Temperature / K")
    p.add_argument("--no_annotation", action="store_true",
                   help="Do not write ΔG‡ / ΔG value on energy diagrams")

    args = p.parse_args()

    cfg_path = Path(args.config)
    cfg = yaml.safe_load(cfg_path.read_text())
    if "geom" not in cfg or "calc" not in cfg:
        sys.exit("[ERROR] YAML must contain 'geom' and 'calc' sections.")

    annotate = not args.no_annotation
    freeze_atoms = cfg["geom"].get("freeze_atoms", None)

    # --- build calculator (freeze passed automatically) -------------
    calc = build_calc(cfg["calc"], freeze_atoms=freeze_atoms)

    # --- run analysis -----------------------------------------------
    analyze_three_structures(
        cfg["geom"]["reac"], cfg["geom"]["ts"], cfg["geom"]["prod"],
        calc,
        temp=args.temp,
        vib_cutoff=100.0,
        freeze_atoms=freeze_atoms,
        out_stem=cfg_path.stem,
        annotate=annotate,
    )

if __name__ == "__main__":
    main()