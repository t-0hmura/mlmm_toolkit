# ONIOM Export (`oniom-export`)

## Overview

> **Summary:** Export an Amber-topology ML/MM system into external QM/MM input formats (Gaussian ONIOM or ORCA QM/MM).

### At a glance
- **Use when:** You want to run external QM/MM calculations on a system prepared with mlmm_toolkit.
- **Method:** Reads Amber parm7 topology; maps ML region from model PDB; writes Gaussian `.com`/`.gjf` or ORCA `.inp` based on `--mode` (or `-o` suffix when omitted).
- **Outputs:** Gaussian ONIOM input or ORCA QM/MM input with appropriate layer/connectivity annotations.
- **Next step:** Run the exported input in Gaussian or ORCA.

## Unified Command

```bash
mlmm oniom-export --parm7 real.parm7 -i pocket.pdb --model-pdb ml.pdb \
 -o out.<gjf|com|inp> --mode <g16|orca> -q 0 -m 1
```

- `--mode` is highest priority.
- If `--mode` is omitted, mode is inferred from `-o`:
  - `.gjf`/`.com` -> `g16`
  - `.inp` -> `orca`
- If `--mode` is omitted and `-o` suffix is unknown, command fails.

## Quick chooser

| If you need... | Use |
| --- | --- |
| Gaussian ONIOM input with link-atom annotations | `mlmm oniom-export --mode g16` |
| ORCA QM/MM input with ORCAFF handling | `mlmm oniom-export --mode orca` |

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [oniom_gaussian](oniom_gaussian.md) -- Gaussian-mode details (`--mode g16`)
- [oniom_orca](oniom_orca.md) -- ORCA-mode details (`--mode orca`)
- [mm_parm](mm_parm.md) -- Build Amber topology
- [define_layer](define_layer.md) -- Build/check layer annotations
