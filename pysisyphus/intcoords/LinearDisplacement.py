import numpy as np

from pysisyphus.intcoords.Primitive import Primitive
from pysisyphus.intcoords.derivatives import dq_ld, d2q_ld
from pysisyphus.linalg import cross3, norm3


class LinearDisplacement(Primitive):
    _COMPLEMENT_FD_STEP = 1.0e-5

    def __init__(self, *args, complement=False, **kwargs):
        kwargs["calc_kwargs"] = ("complement", "cross_vec")
        super().__init__(*args, **kwargs)

        self.complement = complement
        self.cross_vec = None

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        raise Exception("Not yet implemented!")

    def calculate(self, coords3d, indices=None, gradient=False):
        if self.cross_vec is None:
            self.set_cross_vec(coords3d, indices)

        return super().calculate(coords3d, indices, gradient)

    @staticmethod
    def _complement_gradient(coords3d, indices, cross_vec):
        """Differentiate the complemented scalar coordinate, including its frame."""
        coords = np.asarray(coords3d, dtype=float)
        flat_coords = coords.reshape(-1)
        gradient = np.zeros_like(flat_coords)
        step = LinearDisplacement._COMPLEMENT_FD_STEP
        local_dofs = [
            3 * int(atom) + axis for atom in indices for axis in range(3)
        ]
        for dof in local_dofs:
            plus = flat_coords.copy()
            minus = flat_coords.copy()
            plus[dof] += step
            minus[dof] -= step
            value_plus = LinearDisplacement._calculate(
                plus.reshape(-1, 3),
                indices,
                complement=True,
                cross_vec=np.asarray(cross_vec, dtype=float).copy(),
            )
            value_minus = LinearDisplacement._calculate(
                minus.reshape(-1, 3),
                indices,
                complement=True,
                cross_vec=np.asarray(cross_vec, dtype=float).copy(),
            )
            gradient[dof] = (value_plus - value_minus) / (2.0 * step)
        return gradient

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, complement=False, cross_vec=None):
        m, o, n = indices
        w_dash = coords3d[n] - coords3d[m]
        w = w_dash / norm3(w_dash)

        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u = u_dash / norm3(u_dash)
        v = v_dash / norm3(v_dash)

        # Vector for cross product to determine first orthogonal direction
        if cross_vec is None:
            cross_vec = LinearDisplacement._get_cross_vec(coords3d, indices)

        base_cross_vec = np.asarray(cross_vec, dtype=float)
        if complement:
            effective_cross_vec = cross3(w, base_cross_vec)
        else:
            effective_cross_vec = base_cross_vec.copy()
        effective_cross_vec /= norm3(effective_cross_vec)

        # Orthogonal direction
        y = cross3(w, effective_cross_vec)
        y /= norm3(y)

        lin_disp = y.dot(u) + y.dot(v)

        if gradient:
            if complement:
                return lin_disp, LinearDisplacement._complement_gradient(
                    coords3d, indices, base_cross_vec,
                )
            row = np.zeros_like(coords3d)
            row[indices] = dq_ld(
                *coords3d[indices].flatten(), *effective_cross_vec,
            ).reshape(-1, 3)
            return lin_disp, row.flatten()

        return lin_disp

    def jacobian(self, coords3d, indices=None):
        if self.cross_vec is None:
            self.set_cross_vec(coords3d, indices)

        return super().jacobian(coords3d, indices)

    @staticmethod
    def _jacobian(coords3d, indices, complement=False, cross_vec=None):
        if cross_vec is None:
            cross_vec = LinearDisplacement._get_cross_vec(coords3d, indices)

        if complement:
            # The complemented direction depends on the terminal-atom axis.
            # Differentiate the total scalar gradient, including that moving
            # frame, rather than treating the generated vector parameter as
            # constant.
            coords = np.asarray(coords3d, dtype=float)
            local_dofs = [
                3 * int(atom) + axis for atom in indices for axis in range(3)
            ]
            jac = np.zeros((len(local_dofs), len(local_dofs)), dtype=float)
            step = LinearDisplacement._COMPLEMENT_FD_STEP
            for column, dof in enumerate(local_dofs):
                plus = coords.copy().reshape(-1)
                minus = coords.copy().reshape(-1)
                plus[dof] += step
                minus[dof] -= step
                grad_plus = LinearDisplacement._calculate(
                    plus.reshape(-1, 3),
                    indices,
                    gradient=True,
                    complement=True,
                    cross_vec=np.asarray(cross_vec, dtype=float).copy(),
                )[1]
                grad_minus = LinearDisplacement._calculate(
                    minus.reshape(-1, 3),
                    indices,
                    gradient=True,
                    complement=True,
                    cross_vec=np.asarray(cross_vec, dtype=float).copy(),
                )[1]
                jac[:, column] = (
                    np.asarray(grad_plus)[local_dofs]
                    - np.asarray(grad_minus)[local_dofs]
                ) / (2.0 * step)
            return jac.reshape(-1)
        cross_vec = cross_vec / norm3(cross_vec)

        return d2q_ld(*coords3d[indices].flatten(), *cross_vec)

    def __str__(self):
        return (
            f"LinearDisplacement({tuple(self.indices)}, complement={self.complement})"
        )
