# Methyltransferase — step-by-step ML/MM workflow
# Methyl transfer: CS1 (SAM 359) – C8 (PHN 360) bond formation

# Step 1: Extract binding pocket
mlmm extract -i complex.pdb -c "SAM,PHN" -l "SAM:1,PHN:-1" -o pocket.pdb

# Step 2: Define ML/MM layers
mlmm define-layer -i complex.pdb --model-pdb pocket.pdb -o r_layered.pdb

# Step 3: Generate MM parameters
mlmm mm-parm -i complex.pdb -l "SAM:1,PHN:-1"

# Step 4: Optimize initial (reactant) structure
mlmm opt -i r_layered.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 01_opt_init

# Step 5: 1D bond scan (CS1 of SAM 359 – C8 of PHN 360, target 1.3 Å)
mlmm scan -i 01_opt_init/final_geometry.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 -s "[('SAM 359 CS1','PHN 360 C8',1.3)]" --out-dir 02_scan

# Step 6: Path optimization (Growing String Method)
mlmm path-opt -i 01_opt_init/final_geometry.pdb 02_scan/stage_01/result.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 03_path_opt

# Step 7: TS optimization
mlmm tsopt -i 03_path_opt/hei.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 04_tsopt

# Step 8: IRC
mlmm irc -i 04_tsopt/final_geometry.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 05_irc

# Step 9: Optimize IRC forward endpoint (use single-frame endpoint, not trajectory PDB)
mlmm opt -i 05_irc/forward_last.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 06_opt_fwd --thresh baker

# Step 10: Optimize IRC backward endpoint
mlmm opt -i 05_irc/backward_last.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 07_opt_bwd --thresh baker

# Step 11: Frequency analysis (R, TS, P)
mlmm freq -i 06_opt_fwd/final_geometry.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 08_freq_reac
mlmm freq -i 04_tsopt/final_geometry.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 09_freq_ts
mlmm freq -i 07_opt_bwd/final_geometry.pdb -l "SAM:1,PHN:-1" --parm complex.parm7 --out-dir 10_freq_prod

# Step 12: Energy diagram (replace placeholder values with actual relative energies)
mlmm energy-diagram -i "[0, 15.6, -43.7]" -o energy_diagram.png --label-x R TS P
