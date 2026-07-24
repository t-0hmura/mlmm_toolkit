import numpy as np

from pysisyphus.intcoords.Torsion import Torsion


class DummyTorsion(Torsion):
    def __init__(self, indices, *args, fix_inner=True, **kwargs):
        self.fix_inner = fix_inner
        kwargs["calc_kwargs"] = ("fix_inner",)
        super().__init__(indices, *args, **kwargs)
        self.log("DummyTorsion is never checked for collinear atoms!")

    @staticmethod
    def get_fourth_coords(coords3d, indices, r=1.889, theta=90):
        r"""
        M       N <- add
         \     /
          u   v
           \ /
            P
        m, o, p, n = indices
        """
        _, O_, P_ = indices
        # Center
        P = coords3d[P_]
        # Bond, pointing away from O to M
        u = coords3d[O_] - P
        # Direction of u along x axis (left/right)
        sign = np.sign(u[0])
        # Polar coordinates
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        # Translate from center with correct orientation
        fourth_coords = P + (sign * x, sign * y, 0.0)
        return fourth_coords

    @staticmethod
    def get_coords3d_and_indices_ext(coords3d, indices):
        fourth_coords = DummyTorsion.get_fourth_coords(coords3d, indices)
        fourth_ind = len(coords3d)
        coords3d_ext = np.zeros((len(coords3d) + 1, 3))
        coords3d_ext[:fourth_ind] = coords3d
        coords3d_ext[fourth_ind] = fourth_coords
        indices_ext = indices + [fourth_ind]
        return coords3d_ext, indices_ext

    @staticmethod
    def _weight(*args, **kwargs):
        return 1

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, fix_inner=False):
        coords3d_ext, indices_ext = DummyTorsion.get_coords3d_and_indices_ext(
            coords3d, indices
        )

        results = Torsion._calculate(coords3d_ext, indices_ext, gradient=gradient)
        if gradient:
            val, grad = results

            # Remove entries that belong to the dummy atom.
            grad = grad[:-3]

            # Zero out contributions of the inner two atoms in the torsion.
            # So basically only the atom at the first index moves.
            #
            # This usually degrades optimization convergence.
            if fix_inner:
                grad.reshape(-1, 3)[indices_ext[1:3]] = 0.0
            return val, grad.flatten()
        else:
            return results

    @staticmethod
    def _jacobian(coords3d, indices, fix_inner=False):
        """Finite-difference the public three-atom gradient.

        The generated four-atom torsion Jacobian cannot be reused directly:
        the fourth atom is synthetic and coordinate dependent.
        """
        coords = np.asarray(coords3d, dtype=float)
        local_dofs = [
            3 * int(atom) + axis for atom in indices for axis in range(3)
        ]
        jac = np.zeros((len(local_dofs), len(local_dofs)), dtype=float)
        step = 1.0e-5
        for column, dof in enumerate(local_dofs):
            plus = coords.copy().reshape(-1)
            minus = coords.copy().reshape(-1)
            plus[dof] += step
            minus[dof] -= step
            grad_plus = DummyTorsion._calculate(
                plus.reshape(-1, 3),
                indices,
                gradient=True,
                fix_inner=fix_inner,
            )[1]
            grad_minus = DummyTorsion._calculate(
                minus.reshape(-1, 3),
                indices,
                gradient=True,
                fix_inner=fix_inner,
            )[1]
            jac[:, column] = (
                np.asarray(grad_plus)[local_dofs]
                - np.asarray(grad_minus)[local_dofs]
            ) / (2.0 * step)
        return jac.reshape(-1)
