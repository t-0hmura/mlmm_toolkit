# [1] https://doi.org/10.1007/s002140050387
#     Bofill, 1998


import numpy as np

from pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer

import torch

class RSIRFOptimizer(TSHessianOptimizer):
    def optimize(self):
        energy, gradient, H, eigvals, eigvecs, resetted = self.housekeeping()
        self.update_ts_mode(eigvals, eigvecs)

        self.log(
            "Using projection to construct image potential gradient "
            f"and hessian for root(s) {self.roots}."
        )
        # Ensure gradient is in the same subspace as eigvecs (active DOFs)
        if isinstance(H, torch.Tensor):
            if gradient.size(0) != eigvecs.size(0):
                gradient = self.active_from_full(gradient)
        else:
            if gradient.size != eigvecs.shape[0]:
                gradient = self.active_from_full(gradient)
        # Projection matrix to construct g* and H*
        if isinstance(H, torch.Tensor):
            P = torch.eye(gradient.size(0), device=H.device, dtype=H.dtype)
            for root in self.roots:
                trans_vec = eigvecs[:, root]
                P -= 2 * torch.outer(trans_vec, trans_vec)
            H_star = P @ H
            eigvals_, eigvecs_ = torch.linalg.eigh(H_star)
        else:
            P = np.eye(gradient.size)
            for root in self.roots:
                trans_vec = eigvecs[:, root]
                P -= 2 * np.outer(trans_vec, trans_vec)
            H_star = P.dot(H)
            eigvals_, eigvecs_ = np.linalg.eigh(H_star)

        # Once PHVA has identified the physical first-order saddle mode, any
        # remaining negative roots of the unprojected Cartesian Hessian are
        # outside that target (commonly translation/rotation FD artifacts).
        # RS-I-RFO minimizes the image potential, so leaving those roots
        # negative makes that minimization indefinite and can drive a large
        # nonphysical step. Reflect only these residual image-Hessian roots to
        # positive curvature; the selected TS root has already been flipped by
        # P above.
        if self._physical_ts_mode is not None:
            residual_negative = eigvals_ < -self.small_eigval_thresh
            if isinstance(eigvals_, torch.Tensor):
                residual_count = int(residual_negative.sum().item())
                if residual_count:
                    eigvals_ = torch.where(
                        residual_negative, -eigvals_, eigvals_
                    )
            else:
                residual_count = int(np.count_nonzero(residual_negative))
                if residual_count:
                    eigvals_ = np.where(residual_negative, -eigvals_, eigvals_)
            if residual_count:
                self.log(
                    "Stabilized "
                    f"{residual_count} residual negative image-Hessian root(s) "
                    "outside the PHVA-verified TS mode."
                )
        # Neglect small eigenvalues
        eigvals_, eigvecs_ = self.filter_small_eigvals(eigvals_, eigvecs_)

        if isinstance(H, torch.Tensor):
            grad_star = P @ gradient
        else:
            grad_star = P.dot(gradient)
        step = self.get_rs_step(eigvals_, eigvecs_, grad_star, name="RS-I-RFO")

        step = self.apply_saddle_recovery_step(step)
        self.predicted_energy_changes.append(self.rfo_model(gradient, self.cur_H, step))

        step = self.full_from_active(step)
        if isinstance(step, torch.Tensor):
            step = step.cpu().numpy()
        return step
