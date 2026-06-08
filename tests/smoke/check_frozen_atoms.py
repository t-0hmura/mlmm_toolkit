#!/usr/bin/env python
import math
import sys


def coords(path):
    result = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                result.append(
                    (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                )
    return result


def main(argv):
    ref_pdb, out_pdb, atoms_csv, label = argv[1:5]
    indices = [int(x) for x in atoms_csv.split(",") if x.strip()]
    ref = coords(ref_pdb)
    out = coords(out_pdb)
    if len(ref) != len(out):
        raise SystemExit(f"[freeze-check] {label}: atom count mismatch {len(ref)} != {len(out)}")

    max_disp = 0.0
    for idx_1based in indices:
        i = idx_1based - 1
        disp = math.sqrt(sum((out[i][j] - ref[i][j]) ** 2 for j in range(3)))
        max_disp = max(max_disp, disp)
        if disp > 1.0e-5:
            raise SystemExit(f"[freeze-check] {label}: atom {idx_1based} moved {disp:.6e} A")
    print(f"[freeze-check] {label}: max frozen displacement {max_disp:.3e} A")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit("usage: check_frozen_atoms.py REF_PDB OUT_PDB ATOMS_CSV LABEL")
    main(sys.argv)
