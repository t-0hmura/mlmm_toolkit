# Methyltransferase — automated end-to-end workflow via `mlmm all`
# Single PDB + scan mode: extract → scan → MEP → TS → IRC → freq

mlmm all -i complex.pdb -c "SAM,PHN" -l "SAM:1,PHN:-1" -s "[('SAM 359 CS1','PHN 360 C8',1.3)]" --tsopt True --thermo True -o result_all > all.log 2>&1
