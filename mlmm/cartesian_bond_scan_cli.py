"""
CLI wrapper to run CartesianBondScan from a YAML input file.

Usage
-----
$ bond_scan my_scan.yml
"""

from __future__ import annotations

import sys
import argparse
from pathlib import Path
from typing import Any, Dict, Tuple, List

import yaml
from mlmm.cartesian_bond_scan import CartesianBondScan

# ---------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------
def _as_dict(x: Any) -> Dict:          # None → {}
    return {} if x is None else dict(x)


def build_scanner(cfg: Dict) -> CartesianBondScan:
    """Convert YAML dict → CartesianBondScan constructor kwargs."""
    geom_cfg = cfg.get("geom", {})
    scan_cfg = cfg.get("scan", {})
    calc_cfg = cfg.get("calc", {})
    misc_cfg = cfg.get("misc", {})

    scanner = CartesianBondScan(
        input_path = geom_cfg.get("fn", None),
        freeze_atoms = list(geom_cfg.get("freeze_atoms", [])),
        scan_bond  = tuple(scan_cfg.get("bond", (0, 1))),
        scan_step  = float(scan_cfg.get("step", 0.05)),
        scan_range= tuple(scan_cfg.get("scan_range", (1.4, 4.0))),
        init_thresh= scan_cfg.get("init_thresh", "gau"),
        thresh     = scan_cfg.get("thresh", "gau"),
        max_cycles = int(scan_cfg.get("max_cycles", 10_000)),
        out_dir    = scan_cfg.get("out_dir", "./dump/cart_scan/"),
        out_fn     = scan_cfg.get("out_fn", "final_geometries.trj"),
        mlmm_kwargs= calc_cfg,
    )
    return scanner


# ---------------------------------------------------------------------
# CLI Entry
# ---------------------------------------------------------------------
def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(
        prog="bond_scan",
        description="Run CartesianBondScan from YAML input",
    )
    parser.add_argument("input", help="YAML input file")
    args = parser.parse_args(argv)

    cfg_path = Path(args.input)
    if not cfg_path.exists():
        sys.exit(f"[ERROR] YAML file not found: {cfg_path}")

    with cfg_path.open() as f:
        cfg = yaml.safe_load(f)

    scanner = build_scanner(cfg)
    scanner.run()


if __name__ == "__main__":
    main()
