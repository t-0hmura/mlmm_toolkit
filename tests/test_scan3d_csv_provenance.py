from __future__ import annotations

import json

import pandas as pd
from click.testing import CliRunner


def test_plot_only_csv_emits_complete_nullable_calculator_schema(
    tmp_path, monkeypatch
) -> None:
    from mlmm.cli import cli as root_cli
    from mlmm.workflows import scan3d

    rows = []
    for i, d1 in enumerate((1.0, 2.0)):
        for j, d2 in enumerate((1.5, 2.5)):
            for k, d3 in enumerate((2.0, 3.0)):
                rows.append(
                    {
                        "i": i,
                        "j": j,
                        "k": k,
                        "d1_A": d1,
                        "d2_A": d2,
                        "d3_A": d3,
                        "energy_hartree": float(i + 2 * j + 4 * k) / 1000.0,
                    }
                )
    csv_path = tmp_path / "surface.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    monkeypatch.setattr(scan3d, "_VOLUME_GRID_N", 3)
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(
        root_cli,
        ["scan3d", "--csv", str(csv_path), "--out-json", "--out-dir", str(out_dir)],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    for key in (
        "mlip_backend",
        "mlip_model",
        "mlip_precision",
        "mm_backend",
        "link_atom_method",
        "use_cmap",
        "charge",
        "spin",
    ):
        assert key in payload
        assert payload[key] is None


def test_fresh_scan3d_calculator_schema_uses_resolved_values() -> None:
    from mlmm.workflows.scan3d import _result_calculator_fields

    fields = _result_calculator_fields(
        {
            "backend": "orb",
            "orb_model": "orb_v3_conservative_omol",
            "orb_precision": "float64",
            "mm_backend": "openmm",
            "link_atom_method": "fixed",
            "use_cmap": True,
            "model_charge": -1,
            "model_mult": 2,
        }
    )

    assert fields == {
        "mlip_backend": "orb",
        "mlip_model": "orb_v3_conservative_omol",
        "mlip_precision": "fp64",
        "mm_backend": "openmm",
        "link_atom_method": "fixed",
        "use_cmap": True,
        "charge": -1,
        "spin": 2,
    }
