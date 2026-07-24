from __future__ import annotations

import numpy as np
import pytest

from mlmm.backends import xtb_embedcharge_correction as xtb


def test_dormant_xtb_parsers_and_unit_sign_accept_d_exponents(tmp_path) -> None:
    engrad = tmp_path / "sample.engrad"
    engrad.write_text(
        "# The current total energy in Eh\n"
        "-1.250000000000D+00\n"
        "# The current gradient in Eh/bohr\n"
        "1.0D-02\n-2.0D-02\n3.0D-02\n",
        encoding="utf-8",
    )
    energy, gradient = xtb._parse_engrad(engrad, 1)
    assert energy == pytest.approx(-1.25)
    np.testing.assert_allclose(gradient, [[0.01, -0.02, 0.03]])

    pcgrad = tmp_path / "pcgrad"
    pcgrad.write_text("4.0D-02 -5.0D-02 6.0D-02\n", encoding="utf-8")
    np.testing.assert_allclose(
        xtb._parse_pcgrad(pcgrad, 1), [[0.04, -0.05, 0.06]],
    )

    energy_ev, forces = xtb._convert_units(
        energy_ha=energy,
        gradient_ha_bohr=gradient,
    )
    assert energy_ev == pytest.approx(energy * xtb.EV_PER_HARTREE)
    np.testing.assert_allclose(
        forces, -gradient * xtb.FORCE_EV_ANG_PER_HA_BOHR,
    )


def test_dormant_xtb_parsers_reject_truncated_records(tmp_path) -> None:
    engrad = tmp_path / "truncated.engrad"
    engrad.write_text(
        "# The current total energy in Eh\n-1.0\n"
        "# The current gradient in Eh/bohr\n0.1\n0.2\n",
        encoding="utf-8",
    )
    with pytest.raises(xtb.XTBEmbedError, match="expected 3, got 2"):
        xtb._parse_engrad(engrad, 1)

    pcgrad = tmp_path / "pcgrad"
    pcgrad.write_text("0.1 0.2 0.3\n", encoding="utf-8")
    with pytest.raises(xtb.XTBEmbedError, match="expected 2, got 1"):
        xtb._parse_pcgrad(pcgrad, 2)


def test_embedding_numerical_hessian_uses_negative_force_derivative() -> None:
    stiffness = np.diag([2.0, 3.0, 4.0])

    def forces(q_coords, _mm_coords):
        return (-stiffness @ np.asarray(q_coords).reshape(3)).reshape(1, 3)

    force0, hessian = xtb._numerical_hessian_from_forces(
        forces,
        np.array([[0.2, -0.3, 0.4]]),
        np.empty((0, 3)),
        step_ang=1.0e-4,
    )
    np.testing.assert_allclose(force0, [[-0.4, 0.9, -1.6]])
    np.testing.assert_allclose(hessian, stiffness, atol=1.0e-10)
