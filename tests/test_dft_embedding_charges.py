"""DFT electrostatic embedding must not double-count the ML-MM electrostatics.

With embedding, the DFT Hamiltonian holds the interaction between the ML-region
density and the MM point charges. The MM layers therefore run with the ML-region
charges set to zero, and each link-host charge is shifted onto its MM neighbours.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pmd = pytest.importorskip("parmed")

BOHR = 0.529177210903


def _structure(atoms, bonds=()):
    """Build a ParmEd structure from (name, element, Z, mass, charge, rmin, eps, xyz, resname)."""
    s = pmd.Structure()
    for resnum, (name, _el, z, mass, charge, rmin, eps, _xyz, resname) in enumerate(atoms, 1):
        atom_type = pmd.AtomType(name.upper(), None, mass, z)
        atom_type.set_lj_params(eps, rmin)
        atom = pmd.Atom(
            name=name, type=name.upper(), charge=charge, mass=mass, atomic_number=z
        )
        atom.atom_type = atom_type
        s.add_atom(atom, resname, resnum)
    if bonds:
        bond_type = pmd.BondType(300.0, 1.53)
        s.bond_types.append(bond_type)
        for i, j in bonds:
            s.bonds.append(pmd.Bond(s.atoms[i], s.atoms[j], type=bond_type))
    s.coordinates = np.array([a[7] for a in atoms], dtype=float)
    return s


def _chain():
    # C0-C1 | C2-C3-C4 : ML = {0, 1}, cut bond 1-2 (link host 2).
    atoms = [
        (f"C{i}", "C", 6, 12.011, q, 1.9, 0.086, (1.53 * i, 0.0, 0.0), "MOL")
        for i, q in enumerate([-0.1, 0.2, 0.3, -0.4, 0.15])
    ]
    return _structure(atoms, bonds=[(0, 1), (1, 2), (2, 3), (3, 4)])


def test_ml_charges_are_zeroed_and_link_host_charge_is_shifted() -> None:
    from mlmm.backends.mlmm_calc import apply_dft_embedding_charges

    s = _chain()
    mm_before = sum(a.charge for a in s.atoms[2:])
    summary = apply_dft_embedding_charges(s, [0, 1], [(2, 3)])

    charges = [a.charge for a in s.atoms]
    assert charges[:2] == [0.0, 0.0]
    assert charges[2] == 0.0  # link host emptied ...
    assert charges[3] == pytest.approx(-0.4 + 0.3)  # ... onto its only MM neighbour
    assert charges[4] == pytest.approx(0.15)
    assert sum(charges[2:]) == pytest.approx(mm_before)  # MM charge conserved
    assert summary["ml_charge_removed"] == pytest.approx(0.1)
    assert summary["link_hosts"] == 1


def test_link_host_without_mm_neighbour_is_rejected() -> None:
    from mlmm.backends.mlmm_calc import apply_dft_embedding_charges

    s = _chain()
    with pytest.raises(ValueError, match="link host"):
        # Atoms 2, 3 and 4 are all link hosts: host 2 has no MM receiver.
        apply_dft_embedding_charges(s, [0, 1], [(2, 3), (2, 4), (2, 5)])


def test_modified_charges_survive_a_parm7_round_trip(tmp_path) -> None:
    from parmed.amber import AmberParm

    from mlmm.backends.mlmm_calc import apply_dft_embedding_charges

    parm7 = tmp_path / "chain.parm7"
    AmberParm.from_structure(_chain()).save(str(parm7))
    loaded = pmd.load_file(str(parm7))
    apply_dft_embedding_charges(loaded, [0, 1], [(2, 3)])
    loaded.save(str(parm7), overwrite=True)

    reread = [a.charge for a in pmd.load_file(str(parm7)).atoms]
    assert reread == pytest.approx([0.0, 0.0, 0.0, -0.1, 0.15], abs=1.0e-6)


# --- end-to-end: a charged ML ion and one distant MM counter-charge ----------


def _ion_pair(tmp_path: Path, distance: float):
    from parmed.amber import AmberParm

    s = _structure([
        ("LI", "Li", 3, 6.94, 1.0, 1.137, 0.0183, (0.0, 0.0, 0.0), "LI"),
        ("CL", "Cl", 17, 35.45, -1.0, 2.513, 0.0356, (distance, 0.0, 0.0), "CL"),
    ])
    parm7 = tmp_path / f"pair_{distance:g}.parm7"
    pdb = tmp_path / f"pair_{distance:g}.pdb"
    model = tmp_path / f"model_{distance:g}.pdb"
    AmberParm.from_structure(s).save(str(parm7), overwrite=True)
    s.save(str(pdb), overwrite=True)
    s[":LI"].save(str(model), overwrite=True)
    return pdb, parm7, model, s.coordinates


def _total_energy_hartree(tmp_path: Path, distance: float, embed: bool) -> float:
    pytest.importorskip("pyscf")
    from pysisyphus.constants import AU2EV

    from mlmm.backends.mlmm_calc import MLMMCore
    from mlmm.core.dft_settings import resolve_dft_settings

    pdb, parm7, model, coords = _ion_pair(tmp_path, distance)
    settings = resolve_dft_settings({
        "model_charge": 1, "model_mult": 1,
        "embedcharge": embed, "embedcharge_cutoff": 100.0,
        "dft": {"func_basis": "hf/sto-3g", "engine": "cpu", "lowmem": False,
                "density_fit": False, "scf_stepwise_grid": False,
                "conv_tol": 1.0e-11},
    })
    core = MLMMCore(
        input_pdb=str(pdb), real_parm7=str(parm7), model_pdb=str(model),
        model_charge=1, model_mult=1, backend="dft",
        dft_settings=settings.to_dict(), embedcharge=embed,
        embedcharge_cutoff=100.0, ml_device="cpu", mm_device="cpu",
    )
    return float(core.compute(np.asarray(coords))["energy"]) / AU2EV


@pytest.mark.timeout(600)
@pytest.mark.parametrize("distance", [10.0, 20.0])
def test_embedding_does_not_double_count_the_ion_pair_coulomb(tmp_path, distance) -> None:
    """ME already holds the Li+...Cl- Coulomb through the MM charges; EE must
    replace it with the QM-charge interaction, not add a second copy.

    For a nearly unpolarizable Li+, E(EE) - E(ME) is a tiny polarization
    energy. With the double count it would be about -1/R (-0.053 Eh at 10 A).
    """
    e_me = _total_energy_hartree(tmp_path, distance, embed=False)
    e_ee = _total_energy_hartree(tmp_path, distance, embed=True)
    double_count = -1.0 / (distance / BOHR)

    assert abs(e_ee - e_me) < 1.0e-4
    assert abs(e_ee - e_me - double_count) > 0.9 * abs(double_count)
