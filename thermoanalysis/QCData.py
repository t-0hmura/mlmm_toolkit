import re
from collections.abc import Mapping

import numpy as np


def symmetry_number_from_point_group(point_group):
    """Return the external rotational symmetry number for a point group."""
    fixed = {
        "c1": 1,
        "ci": 1,
        "cs": 1,
        "cinf": 1,
        "cinfv": 1,
        "dinfh": 2,
        "t": 12,
        "th": 12,
        "td": 12,
        "o": 24,
        "oh": 24,
        "i": 60,
        "ih": 60,
        "kh": 1,
    }
    pg = str(point_group).lower()
    if pg in fixed:
        return fixed[pg]
    match = re.fullmatch(r"([cds])(\d+)[a-z]*", pg)
    if match is None:
        raise ValueError(f"Specified point group '{pg}' is invalid")
    family, order_text = match.groups()
    symmetry_number = int(order_text)
    if symmetry_number < 1:
        raise ValueError("Point-group order must be positive")
    if family == "d":
        symmetry_number *= 2
    elif family == "s":
        if symmetry_number % 2:
            raise ValueError(
                "Odd-order Sn labels must use the corresponding C point group"
            )
        symmetry_number //= 2
    return symmetry_number


def detect_point_group_and_symmetry_number(atomic_numbers, coords3d):
    """Detect molecular point group and external rotational symmetry number."""
    atomic_numbers = np.asarray(atomic_numbers, dtype=int).reshape(-1)
    coords3d = np.asarray(coords3d, dtype=float).reshape(-1, 3)
    if len(atomic_numbers) != len(coords3d):
        raise ValueError("atomic_numbers and coords3d must describe the same atoms")
    if len(atomic_numbers) == 0:
        raise ValueError("at least one atom is required for symmetry detection")
    if np.any(atomic_numbers < 1) or np.any(atomic_numbers > 92):
        return "C1", 1, "auto-fallback"
    if len(atomic_numbers) == 1:
        return "Kh", 1, "auto"

    try:
        import pymsym
        from ase.data import chemical_symbols
    except ImportError:
        return "C1", 1, "auto-fallback"

    try:
        elements = [
            pymsym.Element(
                name=chemical_symbols[int(atomic_number)],
                coordinates=coords,
            )
            for atomic_number, coords in zip(atomic_numbers, coords3d)
        ]
        with pymsym.Context(elements=elements) as context:
            point_group = str(context.find_symmetry())
        if point_group == "D0h":
            point_group = "Dinfh"
        elif point_group == "C0v":
            point_group = "Cinfv"
        symmetry_number = symmetry_number_from_point_group(point_group)
    except Exception:
        return "C1", 1, "auto-fallback"
    if symmetry_number < 1:
        return "C1", 1, "auto-fallback"
    return point_group, symmetry_number, "auto"

from thermoanalysis.constants import C, ANG2M, AMU2KG, PLANCK, KB


class QCData:
    def __init__(
        self,
        inp,
        point_group="c1",
        mult=None,
    ):
        self.point_group = point_group.lower()
        self.symmetry_number = self.get_symmetry_number()
        if not isinstance(inp, Mapping):
            raise TypeError(
                "QCData in mlmm expects a mapping with keys: "
                "coords3d, wavenumbers, scf_energy, masses, mult."
            )
        data = dict(inp)

        if mult is not None:
            data["mult"] = mult

        # Actually set data
        self.set_data(data)

        self.standard_orientation()
        I = self.inertia_tensor()
        w, v = np.linalg.eigh(I)
        self._linear = (abs(w[0]) < 1e-8) and (abs(w[1] - w[2]) < 1e-8)

    def set_data(self, data):
        expect = set("coords3d wavenumbers scf_energy masses mult".split())
        present = set(data.keys())
        missing = expect - present
        assert len(missing) == 0, f"Keys '{missing}' are missing!"

        self.masses = np.array(data["masses"], dtype=float)
        self.wavenumbers = np.array(data["wavenumbers"], dtype=float)
        self.coords3d = np.array(data["coords3d"], dtype=float).reshape(-1, 3)
        assert self.coords3d.size == 3 * len(self.masses)
        self.scf_energy = float(data["scf_energy"])
        self.mult = int(data["mult"])

    @property
    def atom_num(self):
        return len(self.masses)

    @property
    def M(self):
        """Molecular mass.

        Returns
        -------
        M : float
            Total molecular mass in amu.
        """
        return self.masses.sum()

    @property
    def mult(self):
        """Multiplicity.

        Returns
        -------
        2S+1 : int
            Multiplicity.
        """
        return self._mult

    @mult.setter
    def mult(self, mult):
        self._mult = mult

    @property
    def vib_frequencies(self):
        """Vibrational frequencies.

        Returns
        -------
        vibfreqs : np.array
            Vibrational frequencies in 1/s.
        """
        return C * self.wavenumbers * 100

    @property
    def is_linear(self):
        """Wether the molecule is linear.

        Returns
        -------
        is_linear : bool
            Wether the molecule is linear.
        """
        return self._linear

    @property
    def is_atom(self):
        """Wether the 'molecule' consists of only an atom.

        Returns
        -------
        is_atoms : bool
            Wether the molecule is only one atom.
        """
        return len(self.masses) == 1

    @property
    def rot_temperatures(self):
        """Rotational temperatures in K.

        Returns
        -------
        rot_temps : np.array
            Rotational temperatures in K.
        """
        if self.is_atom:
            return np.full(3, np.nan)

        self.standard_orientation()
        I = self.inertia_tensor() * ANG2M**2 * AMU2KG
        w, v = np.linalg.eigh(I)
        rot_temps = PLANCK**2 / (8 * np.pi**2 * w * KB)
        return rot_temps

    def get_symmetry_number(self, point_group=None):
        """Symmetry number for rotatioanl partiton function.

        Returns
        -------
        symmetry_number : int
            Symmetry number for calculation of rotational terms.
        """
        if point_group is None:
            point_group = self.point_group
        return symmetry_number_from_point_group(point_group)

    def inertia_tensor(self):
        """Inertita tensor.

                              | x² xy xz |
        (x y z)^T . (x y z) = | xy y² yz |
                              | xz yz z² |
        Returns
        -------
        I : np.array, shape (3, 3)
            Ineratia tensor  in units of Angstrom² * amu.
        """
        x, y, z = self.coords3d.T
        squares = np.sum(self.coords3d**2 * self.masses[:, None], axis=0)
        I_xx = squares[1] + squares[2]
        I_yy = squares[0] + squares[2]
        I_zz = squares[0] + squares[1]
        I_xy = -np.sum(self.masses * x * y)
        I_xz = -np.sum(self.masses * x * z)
        I_yz = -np.sum(self.masses * y * z)
        I = np.array(((I_xx, I_xy, I_xz), (I_xy, I_yy, I_yz), (I_xz, I_yz, I_zz)))
        return I

    @property
    def average_moment_of_inertia(self):
        """Average moment of inertia in Angstrom² * amu."""
        w, _ = np.linalg.eigh(self.inertia_tensor())
        I_avg = np.mean(w)
        return I_avg

    @property
    def center_of_mass(self):
        """Returns the center of mass.

        Returns
        -------
        R : np.array, shape (3, )
            Center of mass in Angstrom.
        """
        return 1 / self.M * np.sum(self.coords3d * self.masses[:, None], axis=0)

    def principal_axes_are_aligned(self):
        """Check if the principal axes are aligned with the cartesian axes.

        Returns
        -------
        aligned : bool
            Wether the principal axes are aligned or not.
        eigenvectors : np.array, shape (3, 3)
            Eigenvectors of the inertia tensor, as columns.
        """
        w, v = np.linalg.eigh(self.inertia_tensor())
        return np.allclose(v, np.eye(3)), v

    def align_principal_axes(self):
        """Align the principal axes to the cartesian axes.

        https://math.stackexchange.com/questions/145023
        """
        I = self.inertia_tensor()
        w, v = np.linalg.eigh(I)
        self.coords3d = v.T.dot(self.coords3d.T).T

    def standard_orientation(self):
        """Bring molecule in standard orientation."""
        # Translate center of mass to cartesian origin
        self.coords3d -= self.center_of_mass
        # Try to rotate the principal axes onto the cartesian axes
        for _ in range(5):
            self.align_principal_axes()
            aligned, vecs = self.principal_axes_are_aligned()
            if aligned:
                break
