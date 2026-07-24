# xTB routing in v0.3.3

Electronic embedding is unavailable in mlmm-toolkit v0.3.3.
`--embedcharge`, `--embedcharge-cutoff`, and `calc.embedcharge: true` are
compatibility inputs that fail before calculation. Do not install xTB to make
those options work; remove them and use the default mechanical embedding.

The retired experimental correction double-counted ML--MM electrostatics
already retained by the subtractive ONIOM expression and evaluated an uncapped
model inconsistent with the link-H high-level system. Results produced with
that path should be rerun.

xTB remains usable only as an independently defined ASE custom calculator.
Route that request through `--calc-file` and the custom-calculator guidance in
`SKILL.md`; it is not the retired embed-charge correction.
