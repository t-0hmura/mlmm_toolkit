"""Regression tests for the opt-in IRC energy-stop bypass and the IRC workflow
output handling around it."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from mlmm.core.defaults import IRC_KW
from pysisyphus.irc.IRC import IRC


def _irc_stop_probe(*, never_stop: bool, increased: bool, converged: bool) -> IRC:
    irc = object.__new__(IRC)
    irc.never_stop = never_stop
    irc.energy_increased = increased
    irc.energy_converged = converged
    irc.never_stop_energy_increase_bypasses = 0
    irc.never_stop_energy_convergence_bypasses = 0
    return irc


def test_never_stop_is_opt_in() -> None:
    assert IRC_KW["never_stop"] is False


def test_default_irc_stops_on_energy_increase_or_plateau() -> None:
    assert _irc_stop_probe(
        never_stop=False, increased=True, converged=False
    )._energy_stop_message() == "Energy increased!"
    assert _irc_stop_probe(
        never_stop=False, increased=False, converged=True
    )._energy_stop_message() == "Energy converged!"


def test_never_stop_ignores_energy_only_stops() -> None:
    assert _irc_stop_probe(
        never_stop=True, increased=True, converged=False
    )._energy_stop_message() == ""
    assert _irc_stop_probe(
        never_stop=True, increased=False, converged=True
    )._energy_stop_message() == ""


def test_directional_endpoint_energy_fields_keep_legacy_aliases() -> None:
    from mlmm.workflows.irc import _directional_endpoint_energy_fields

    fields = _directional_endpoint_energy_fields([-10.0, -9.0, -11.0], -8.5)

    assert fields["energy_first_hartree"] == -10.0
    assert fields["energy_last_hartree"] == -11.0
    assert fields["energy_ts_hartree"] == -8.5
    assert fields["endpoint_energy_orientation"] == "finished_first_to_finished_last"
    assert fields["energy_reactant_hartree"] == fields["energy_first_hartree"]
    assert fields["energy_product_hartree"] == fields["energy_last_hartree"]


def test_downhill_endpoint_energy_orientation_is_explicit() -> None:
    from mlmm.workflows.irc import _directional_endpoint_energy_fields

    fields = _directional_endpoint_energy_fields(
        [-10.0, -11.0],
        -10.0,
        orientation="downhill_first_to_downhill_last",
    )

    assert (
        fields["endpoint_energy_orientation"]
        == "downhill_first_to_downhill_last"
    )


def test_irc_requires_a_direction_but_accepts_downhill() -> None:
    from mlmm.workflows.irc import _validate_irc_directions

    with pytest.raises(Exception, match="at least one IRC direction"):
        _validate_irc_directions(
            {"forward": False, "backward": False, "downhill": False}
        )
    _validate_irc_directions(
        {"forward": False, "backward": False, "downhill": True}
    )


def test_workflow_uses_engine_normalized_prefix(tmp_path) -> None:
    from mlmm.workflows.irc import _irc_output_path

    irc = object.__new__(IRC)
    irc.out_dir = tmp_path
    irc.prefix = "segment_"

    assert _irc_output_path(irc, "finished_irc_trj.xyz") == (
        tmp_path / "segment_finished_irc_trj.xyz"
    )


def test_real_irc_generation_invalidates_prefixed_directional_outputs(
    tmp_path,
) -> None:
    from mlmm.workflows.irc import _prepare_irc_output_dir

    stale = [
        tmp_path / "segment_forward_irc_trj.xyz",
        tmp_path / "segment_backward_irc.pdb",
        tmp_path / "segment_downhill_irc_trj.xyz",
        tmp_path / "segment_downhill_last.xyz",
        tmp_path / "segment_finished_first.xyz",
        tmp_path / "result.json",
        tmp_path / "summary.json",
    ]
    for path in stale:
        path.write_text("stale\n", encoding="utf-8")
    unrelated = tmp_path / "notes.txt"
    unrelated.write_text("keep\n", encoding="utf-8")

    assert _prepare_irc_output_dir(tmp_path, prefix="segment") == tmp_path.resolve()

    assert all(not path.exists() for path in stale)
    assert unrelated.read_text(encoding="utf-8") == "keep\n"


def test_irc_generation_refuses_to_delete_a_reserved_input(tmp_path) -> None:
    from mlmm.workflows.irc import _prepare_irc_output_dir

    source = tmp_path / "backward_last.cif"
    source.write_text("input\n", encoding="utf-8")

    with pytest.raises(Exception, match="collides"):
        _prepare_irc_output_dir(tmp_path, protected_inputs=(source,))

    assert source.read_text(encoding="utf-8") == "input\n"


@pytest.mark.parametrize(
    "device",
    ["cpu"] + (["cuda"] if torch.cuda.is_available() else []),
)
def test_terminal_hessian_conversion_consumes_tensor_in_place(device) -> None:
    from mlmm.workflows.irc import _consume_mw_hessian_to_cartesian_active

    hessian = torch.tensor(
        [[2.0, 0.5], [0.5, 3.0]],
        dtype=torch.float64,
        device=device,
    )
    original = hessian.detach().cpu().numpy().copy()
    data_ptr = hessian.data_ptr()
    masses = np.array([2.0, 3.0])

    result = _consume_mw_hessian_to_cartesian_active(hessian, masses)

    assert hessian.data_ptr() == data_ptr
    np.testing.assert_allclose(result, masses[:, None] * original * masses[None, :])
    np.testing.assert_allclose(hessian.detach().cpu().numpy(), result)


def test_terminal_hessian_conversion_consumes_numpy_in_place() -> None:
    from mlmm.workflows.irc import _consume_mw_hessian_to_cartesian_active

    hessian = np.array([[2.0, 0.5], [0.5, 3.0]])
    original = hessian.copy()
    masses = np.array([2.0, 3.0])

    result = _consume_mw_hessian_to_cartesian_active(hessian, masses)

    assert result is hessian
    np.testing.assert_allclose(result, masses[:, None] * original * masses[None, :])


def test_result_manifest_collects_prefixed_cif_companions(tmp_path) -> None:
    from mlmm.workflows.irc import _collect_irc_output_files

    class FakeIRC:
        def get_path_for_fn(self, filename: str) -> str:
            return str(tmp_path / f"segment_{filename}")

    for filename in (
        "finished_irc_trj.xyz",
        "finished_irc.pdb",
        "finished_irc.cif",
        "forward_last.xyz",
        "forward_last.cif",
        "forward_first.xyz",
        "forward_first.cif",
    ):
        (tmp_path / f"segment_{filename}").write_text("x", encoding="utf-8")

    assert _collect_irc_output_files(FakeIRC()) == {
        "finished_irc": "segment_finished_irc_trj.xyz",
        "finished_irc_pdb": "segment_finished_irc.pdb",
        "finished_irc_cif": "segment_finished_irc.cif",
        "forward_last": "segment_forward_last.xyz",
        "forward_last_cif": "segment_forward_last.cif",
        "forward_endpoint": "segment_forward_first.xyz",
        "forward_endpoint_cif": "segment_forward_first.cif",
    }
