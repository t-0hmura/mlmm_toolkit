#!/usr/bin/env bash

# Methyltransferase — step-by-step ML/MM workflow.
# Methyl transfer: CS1 (SAM 359) – C8 (PHN 360) bond formation

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
work_dir=${1:-"$PWD/mlmm_methyltransferase_stepwise"}
mkdir -p -- "$(dirname -- "$work_dir")"
mkdir -- "$work_dir" || exit 1  # refuse to overwrite an earlier run
cd -- "$work_dir" || exit 1

# Step 1: Build topology and a topology-matched PDB in this work directory.
mlmm mm-parm -i "$script_dir/complex.pdb" -l "SAM:1,PHN:-1" -o complex || exit 1

# Step 2: Extract the ML region from LEaP's topology-matched PDB.
mlmm extract -i complex.pdb -c "SAM,PHN" -l "SAM:1,PHN:-1" -o pocket.pdb || exit 1

# Step 3: Layer that same topology-matched full-system PDB.
mlmm define-layer -i complex.pdb --model-pdb pocket.pdb -o r_layered.pdb || exit 1

# Step 4: Optimize initial (reactant) structure
mlmm opt -i r_layered.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 01_opt_init || exit 1

# Step 5: 1D bond scan (CS1 of SAM 359 – C8 of PHN 360, target 1.3 Å)
mlmm scan -i 01_opt_init/final_geometry.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 -s "[('SAM 359 CS1','PHN 360 C8',1.3)]" --out-dir 02_scan || exit 1

# Step 6: Path optimization (Growing String Method)
mlmm path-opt -i 01_opt_init/final_geometry.pdb 02_scan/stage_01/result.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 03_path_opt || exit 1

# Step 7: TS optimization
mlmm tsopt -i 03_path_opt/hei.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 04_tsopt || exit 1

# Step 8: IRC
mlmm irc -i 04_tsopt/final_geometry.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 05_irc || exit 1

# Step 9: Optimize the IRC forward endpoint (direction is not chemical identity).
mlmm opt -i 05_irc/forward_first.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 06_opt_forward --thresh baker || exit 1

# Step 10: Optimize the IRC backward endpoint.
mlmm opt -i 05_irc/backward_last.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 07_opt_backward --thresh baker || exit 1

# Step 11: Analyze forward, TS, and backward structures. Inspect the endpoint
# structures and directed bond changes before assigning reactant/product labels.
mlmm freq -i 06_opt_forward/final_geometry.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 08_freq_forward || exit 1
mlmm freq -i 04_tsopt/final_geometry.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 09_freq_ts || exit 1
mlmm freq -i 07_opt_backward/final_geometry.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 10_freq_backward || exit 1

# Step 12: Supply relative energies in forward/TS/backward order and the three
# validated chemical labels explicitly (for example ENERGIES='[0.0, 12.5, 4.3]'
# and ENDPOINT_LABELS='P TS R').
if [[ -n ${ENERGIES:-} && -n ${ENDPOINT_LABELS:-} ]]; then
  read -r -a endpoint_labels <<< "$ENDPOINT_LABELS"
  if [[ ${#endpoint_labels[@]} -ne 3 ]]; then
    printf '%s\n' 'ENDPOINT_LABELS must contain exactly three whitespace-separated labels.' >&2
    exit 2
  fi
  mlmm energy-diagram -i "$ENERGIES" -o energy_diagram.png \
    --label-x "${endpoint_labels[0]}" \
    --label-x "${endpoint_labels[1]}" \
    --label-x "${endpoint_labels[2]}" || exit 1
else
  printf '%s\n' 'After validating endpoint identities, set ENERGIES and ENDPOINT_LABELS to render the diagram.'
fi
