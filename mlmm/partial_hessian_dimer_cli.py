"""
Run Partial-Hessian ML/MM Dimer search from a YAML configuration file
(updated for mass-scaled multi-mode flattening).

Usage
-----
$ partial_hessian_dimer config.yml
"""

from __future__ import annotations
import sys, argparse
from pathlib import Path
from typing import Any, Dict

import yaml
from mlmm.partial_hessian_dimer import PartialHessianDimer


# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------
def _as_dict(node: Any) -> Dict:
    """Cast *node* to dict (None → {}, mapping → dict, else TypeError)."""
    if node is None:
        return {}
    if isinstance(node, dict):
        return dict(node)
    raise TypeError(f"Expected a mapping but got {type(node).__name__}: {node!r}")


def build_runner(cfg: Dict) -> PartialHessianDimer:
    """Translate YAML → PartialHessianDimer(**kwargs)."""

    geom_cfg   = cfg.get("geom",   {})
    dimer_cfg  = cfg.get("dimer",  {})
    flatten_cfg= cfg.get("flatten",{})
    misc_cfg   = cfg.get("misc",   {})
    calc_cfg   = cfg.get("calc",   {})

    kwargs = dict(
        # geometry / IO
        # ---------------------------------------------------------------------
        fn           = geom_cfg.get("fn", "./coord/gs_hei.xyz"),
        freeze_atoms = geom_cfg.get("freeze_atoms", []),
        real_pdb     = geom_cfg.get(
                          "real_pdb",
                          calc_cfg.get("real_pdb", "./parm/complex.pdb")
                      ),
        out_dir      = misc_cfg.get("out_dir", "./dump/dimer/"),
        vib_dir      = misc_cfg.get("vib_dir", "./dump/vib/"),

        # thresholds & optimiser
        # ---------------------------------------------------------------------
        thresh_loose        = dimer_cfg.get("thresh_loose", "gau_loose"),
        thresh              = dimer_cfg.get("thresh", "baker"),
        update_interval_hessian = dimer_cfg.get("update_interval_hessian", 50),
        lobpcg              = dimer_cfg.get("lobpcg", True),
        max_cycles          = dimer_cfg.get("max_cycles", 100_000),
        max_cycles_preopt   = dimer_cfg.get("max_cycles_preopt"),
        mem                 = dimer_cfg.get("mem", 100_000),
        partial_mm_cutoff   = dimer_cfg.get("partial_mm_cutoff", 0.0),
        neg_freq_thresh     = dimer_cfg.get("neg_freq_thresh", 10.0),

        # flatten loop
        # ---------------------------------------------------------------------
        flatten_amp_ang     = flatten_cfg.get("amp_ang",   0.20),
        flatten_max_iter    = flatten_cfg.get("max_iter",  20),
        flatten_sep_cutoff  = flatten_cfg.get("sep_cutoff",2.0),
        flatten_k           = flatten_cfg.get("k", 10),

        # misc
        # ---------------------------------------------------------------------
        dump          = misc_cfg.get("dump", False),

        # low-level kwargs
        # ---------------------------------------------------------------------
        mlmm_kwargs   = calc_cfg,
        dimer_kwargs  = _as_dict(dimer_cfg.get("kwargs")),
    )

    return PartialHessianDimer(**kwargs)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(
        prog="partial_hessian_dimer",
        description="Run Partial-Hessian ML/MM Dimer search from a YAML config."
    )
    parser.add_argument("config", help="YAML configuration file.")
    args = parser.parse_args(argv)

    cfg_path = Path(args.config)
    if not cfg_path.is_file():
        sys.exit(f"[ERROR] YAML file not found: {cfg_path}")

    try:
        with cfg_path.open() as f:
            cfg = yaml.safe_load(f) or {}
    except yaml.YAMLError as e:
        sys.exit(f"[ERROR] Failed to parse YAML: {e}")

    try:
        runner = build_runner(cfg)
    except Exception as e:          # catch bad keys / types early
        sys.exit(f"[ERROR] Invalid configuration: {e}")

    runner.run()


if __name__ == "__main__":
    main()
