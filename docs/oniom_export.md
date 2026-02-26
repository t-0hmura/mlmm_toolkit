# ONIOM Export (`oniom-gaussian` / `oniom-orca`)

## Overview

> **Summary:** Export an Amber-topology ML/MM system into external QM/MM input formats (Gaussian ONIOM or ORCA QM/MM).

### At a glance
- **Use when:** You want to run external QM/MM calculations on a system prepared with mlmm_toolkit.
- **Method:** Reads Amber parm7 topology; maps ML region from model PDB; writes Gaussian `.com`/`.gjf` or ORCA `.inp`.
- **Outputs:** Gaussian ONIOM input or ORCA QM/MM input with appropriate layer/connectivity annotations.
- **Next step:** Run the exported input in Gaussian or ORCA.

## Command Guides

- [oniom-gaussian](oniom_gaussian.md) -- Gaussian-focused narrative, options, and examples
- [oniom-orca](oniom_orca.md) -- ORCA-focused narrative, options, and examples

## Quick chooser

| If you need... | Use |
| --- | --- |
| Gaussian ONIOM input with link-atom annotations | [`oniom-gaussian`](oniom_gaussian.md) |
| ORCA QM/MM input with ORCAFF handling | [`oniom-orca`](oniom_orca.md) |

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [oniom-gaussian](oniom_gaussian.md) -- Gaussian ONIOM exporter
- [oniom-orca](oniom_orca.md) -- ORCA QM/MM exporter
- [mm_parm](mm_parm.md) -- Build Amber topology
- [define_layer](define_layer.md) -- Build/check layer annotations
