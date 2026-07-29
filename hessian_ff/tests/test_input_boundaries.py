from __future__ import annotations

from pathlib import Path
import re
from types import SimpleNamespace

import pytest
import torch

from hessian_ff import loaders, workflows


_DATA_DIR = Path(__file__).resolve().parent / "data" / "small"
_PRMTOP = _DATA_DIR / "complex.parm7"


def test_amber_restart_rejects_malformed_coordinate_field(tmp_path: Path) -> None:
    restart = tmp_path / "broken.rst7"
    restart.write_text(
        "title\n"
        "2\n"
        "  0.0000000  0.0000000BAD-FIELD!!"
        "  1.0000000  0.0000000  0.0000000\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Invalid coordinate field"):
        loaders._read_amber_inpcrd(restart)


@pytest.mark.parametrize("delta", [0.0, -1.0e-4, float("nan"), float("inf")])
def test_fd_hessian_rejects_invalid_delta_before_runtime_load(
    monkeypatch, delta,
) -> None:
    monkeypatch.setattr(
        workflows,
        "_load_runtime",
        lambda **_kwargs: (_ for _ in ()).throw(
            AssertionError("runtime load must not run")
        ),
    )

    with pytest.raises(ValueError, match="finite positive"):
        workflows.torch_hessian(
            "unused.parm7",
            torch.zeros((1, 3), dtype=torch.float64),
            partial_hessian=False,
            hessian_delta=delta,
        )


def test_public_batch_workflow_rejects_empty_tensor(monkeypatch) -> None:
    system = SimpleNamespace(natom=1)
    system.to = lambda **_kwargs: system
    monkeypatch.setattr(workflows, "load_system", lambda *_args, **_kwargs: system)

    with pytest.raises(ValueError, match="coords_batch is empty"):
        workflows._prepare_batch_ff(
            "unused.parm7",
            torch.empty((0, 1, 3), dtype=torch.float64),
            device="cpu",
            double=True,
            force_calc_mode="Analytical",
            with_force=False,
        )


@pytest.mark.skipif(
    not _PRMTOP.exists(),
    reason="benchmark/data/small topology is not available",
)
@pytest.mark.parametrize(
    ("missing_flag", "default_value", "preserved_attr"),
    [
        ("SCEE_SCALE_FACTOR", 1.0 / 1.2, "pair14_inv_scnb"),
        ("SCNB_SCALE_FACTOR", 1.0 / 2.0, "pair14_inv_scee"),
    ],
)
def test_legacy_scaling_defaults_only_the_missing_factor(
    tmp_path: Path,
    missing_flag: str,
    default_value: float,
    preserved_attr: str,
) -> None:
    original = loaders.load_system(_PRMTOP)
    text = _PRMTOP.read_text(encoding="utf-8")
    text, count = re.subn(
        rf"%FLAG {missing_flag}[^\n]*\n.*?(?=%FLAG|\Z)",
        "",
        text,
        flags=re.DOTALL,
    )
    assert count == 1
    modified_path = tmp_path / "missing-scale.parm7"
    modified_path.write_text(text, encoding="utf-8")

    modified = loaders.load_system(modified_path)
    default_attr = (
        "pair14_inv_scee"
        if missing_flag == "SCEE_SCALE_FACTOR"
        else "pair14_inv_scnb"
    )

    torch.testing.assert_close(
        getattr(modified, default_attr),
        torch.full_like(getattr(modified, default_attr), default_value),
    )
    torch.testing.assert_close(
        getattr(modified, preserved_attr),
        getattr(original, preserved_attr),
    )
