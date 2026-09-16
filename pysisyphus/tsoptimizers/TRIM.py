# [1] https://doi.org/10.1016/0009-2614(91)90115-P
#     Helgaker, 1991


import numpy as np
import torch

from pysisyphus._array import as_numpy
from pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer


class TRIM(TSHessianOptimizer):

    def optimize(self):
        energy, gradient, H, eigvals, eigvecs, resetted = self.housekeeping()
        self.update_ts_mode(eigvals, eigvecs)

        # When the Hessian is a partial (active-block) Hessian — e.g. a frozen
        # active site — the eigvecs span only the active DOFs, but the gradient
        # may still be in full-coord space. Mirror RSIRFO/RSPRFO's reduction so
        # eigvecs.T @ gradient is shape-compatible. Without this, TRIM crashes
        # with a `coords(3N,) + step(3N_active,)` broadcast error under
        # freeze_atoms (partial Hessian is the default), which RSIRFO/RSPRFO
        # already guard against.
        if isinstance(H, torch.Tensor):
            if gradient.size(0) != eigvecs.size(0):
                gradient = self.active_from_full(gradient)
        else:
            if gradient.size != eigvecs.shape[0]:
                gradient = self.active_from_full(gradient)

        # The shared trust-region solver works with NumPy arrays.
        eigvals = as_numpy(eigvals)
        eigvecs = as_numpy(eigvecs)
        gradient = as_numpy(gradient)

        self.log(f"Signs of eigenvalue and -vector of root(s) {self.roots} "
                  "will be reversed!")
        # Transform gradient to basis of eigenvectors
        gradient_ = eigvecs.T.dot(gradient)

        # Construct image function by inverting the signs of the eigenvalue and
        # -vector of the mode to follow uphill.
        eigvals_ = eigvals.copy()
        eigvals_[self.roots] *= -1
        if self._physical_ts_mode is not None:
            residual_negative = (
                (eigvals_ < -self.small_eigval_thresh)
                & self._translation_mode_mask(eigvecs)
            )
            residual_negative[self.roots] = False
            residual_count = int(np.count_nonzero(residual_negative))
            if residual_count:
                eigvals_[residual_negative] *= -1
                self.log(
                    "Stabilized "
                    f"{residual_count} negative translational image-Hessian root(s)."
                )
        gradient_ = gradient_.copy()
        gradient_[self.roots] *= -1

        # Minimize the image quadratic even when complementary negative
        # curvature remains. The solver takes the gradient in coordinate space.
        step = self.get_newton_step_on_trust(
            eigvals_, eigvecs, eigvecs @ gradient_
        )

        step_norm = np.linalg.norm(step)
        self.log(f"norm(step)={step_norm:.6f}")

        gradient, step = self.validate_terminal_step_basis(gradient, step)
        step = self.apply_saddle_recovery_step(step)
        self.predicted_energy_changes.append(
            self.quadratic_model(gradient, as_numpy(self.cur_H), step)
        )

        # Expand the step back to full-coord space when the active subspace is in
        # use, so Optimizer.run() can do `geometry.coords + step` without a shape
        # mismatch (same convention as RSIRFOptimizer / RSPRFOptimizer).
        step = self.full_from_active(step)

        return step
