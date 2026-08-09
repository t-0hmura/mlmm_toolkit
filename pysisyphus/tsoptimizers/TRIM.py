# [1] https://doi.org/10.1016/0009-2614(91)90115-P
#     Helgaker, 1991


import numpy as np
import torch

from scipy.optimize import newton
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

        # TRIM uses scipy.optimize.newton + np.nan_to_num internally and is not
        # microiter-capable; coerce torch tensors from the MLIP Hessian path to
        # numpy so the legacy .dot / np.linalg.norm below stay valid.
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
            residual_negative = eigvals_ < -self.small_eigval_thresh
            residual_negative[self.roots] = False
            residual_count = int(np.count_nonzero(residual_negative))
            if residual_count:
                eigvals_[residual_negative] *= -1
                self.log(
                    "Stabilized "
                    f"{residual_count} residual negative image-Hessian root(s) "
                    "outside the PHVA-verified TS mode."
                )
        gradient_ = gradient_.copy()
        gradient_[self.roots] *= -1

        def get_step(mu):
            zetas = -gradient_ / (eigvals_ - mu)
            # Replace nan with 0.
            zetas = np.nan_to_num(zetas)
            # Transform to original basis
            step = eigvecs * zetas
            step = step.sum(axis=1)
            return step

        def get_step_norm(mu):
            return np.linalg.norm(get_step(mu))

        def func(mu):
            return get_step_norm(mu) - self.trust_radius

        mu = 0
        norm0 = get_step_norm(mu)
        if norm0 > self.trust_radius:
            try:
                mu, res = newton(func, x0=mu, full_output=True)
                if not res.converged:
                    raise RuntimeError("newton not converged")
                self.log(f"Using levelshift of μ={mu:.4f}")
                step = get_step(mu)
            except RuntimeError as exc:
                self.log(
                    f"Levelshift newton diverged ({exc}); falling back to "
                    "L_2 truncation at trust_radius."
                )
                step = get_step(0)
                step = step / np.linalg.norm(step) * self.trust_radius
        else:
            self.log("Took pure newton step without levelshift")
            step = get_step(mu)

        step_norm = np.linalg.norm(step)
        self.log(f"norm(step)={step_norm:.6f}")

        self.validate_terminal_saddle_for_step(step)
        step = self.apply_saddle_recovery_step(step)
        self.predicted_energy_changes.append(
            self.quadratic_model(gradient, as_numpy(self.cur_H), step)
        )

        # Expand the step back to full-coord space when the active subspace is in
        # use, so Optimizer.run() can do `geometry.coords + step` without a shape
        # mismatch (same convention as RSIRFOptimizer / RSPRFOptimizer).
        step = self.full_from_active(step)

        return step
