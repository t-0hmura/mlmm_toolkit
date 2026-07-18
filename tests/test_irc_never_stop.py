"""Regression tests for opt-in IRC energy-stop bypass controls."""

from __future__ import annotations

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


def test_workflow_uses_engine_normalized_prefix(tmp_path) -> None:
    from mlmm.workflows.irc import _irc_output_path

    irc = object.__new__(IRC)
    irc.out_dir = tmp_path
    irc.prefix = "segment_"

    assert _irc_output_path(irc, "finished_irc_trj.xyz") == (
        tmp_path / "segment_finished_irc_trj.xyz"
    )


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
    ):
        (tmp_path / f"segment_{filename}").write_text("x", encoding="utf-8")

    assert _collect_irc_output_files(FakeIRC()) == {
        "finished_irc": "segment_finished_irc_trj.xyz",
        "finished_irc_pdb": "segment_finished_irc.pdb",
        "finished_irc_cif": "segment_finished_irc.cif",
        "forward_last": "segment_forward_last.xyz",
        "forward_last_cif": "segment_forward_last.cif",
    }
