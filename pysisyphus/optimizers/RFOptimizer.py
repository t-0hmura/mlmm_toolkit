# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import rms
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.poly_fit import poly_line_search
from pysisyphus.optimizers.gdiis import gdiis, gediis

import torch

class RFOptimizer(HessianOptimizer):
    def __init__(
        self,
        geometry: Geometry,
        line_search: bool = True,
        gediis: bool = False,
        gdiis: bool = True,
        gdiis_thresh: float = 2.5e-3,
        gediis_thresh: float = 1e-2,
        gdiis_test_direction: bool = True,
        max_micro_cycles: int = 25,
        adapt_step_func: bool = False,
        reject_uphill: bool = True,
        uphill_tolerance: float = 1e-8,
        rejection_trust_floor: float = 1e-7,
        max_rejections_at_floor: int = 3,
        **kwargs,
    ) -> None:
        """
        Rational function Optimizer.

        Parameters
        ----------
        geometry
            Geometry to be optimized.
        line_search
            Whether to carry out implicit line searches.
        gediis
            Whether to enable GEDIIS.
        gdiis
            Whether to enable GDIIS.
        gdiis_thresh
            Threshold for rms(forces) to enable GDIIS.
        gediis_thresh
            Threshold for rms(step) to enable GEDIIS.
        gdiis_test_direction
            Whether to the overlap of the RFO step and the GDIIS step.
        max_micro_cycles
            Number of restricted-step microcycles. Disabled by default.
        adapt_step_func
            Whether to switch between shifted Newton and RFO-steps.
        reject_uphill
            Reject energy-increasing minimization trials and retry from the
            previous accepted geometry.
        uphill_tolerance
            Positive energy change tolerated before rejecting a trial.
        rejection_trust_floor
            Emergency trust-radius floor for repeated rejected trials.
        max_rejections_at_floor
            Rejections allowed at the emergency floor before retaining the
            lower-energy point for one final convergence check.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments passed to the Optimizer/HessianOptimizer baseclass.
        """
        super().__init__(
            geometry,
            max_micro_cycles=max_micro_cycles,
            reject_uphill=reject_uphill,
            uphill_tolerance=uphill_tolerance,
            rejection_trust_floor=rejection_trust_floor,
            max_rejections_at_floor=max_rejections_at_floor,
            **kwargs,
        )

        self.line_search = line_search
        self.gediis = gediis
        self.gdiis = gdiis
        self.gdiis_thresh = gdiis_thresh  # Will be compared to rms(step)
        self.gediis_thresh = gediis_thresh  # Will be compared to rms(forces)
        self.gdiis_test_direction = gdiis_test_direction
        self.adapt_step_func = adapt_step_func

        self.successful_gediis = 0
        self.successful_gdiis = 0
        self.successful_line_search = 0

    def _accept_accelerated_step(self, step, ip_step, ref_step):
        """Return an accelerated step only when its full displacement is safe."""
        # The interpolation gradient can cross the NumPy/torch boundary even
        # when the reference RFO step lives on CUDA.  Keep the combined step in
        # the reference step's representation so neither NumPy's implicit
        # ``Tensor.__array__`` nor a device-changing copy is triggered.
        if isinstance(ref_step, torch.Tensor):
            candidate = (
                torch.as_tensor(step, dtype=ref_step.dtype, device=ref_step.device)
                + torch.as_tensor(
                    ip_step, dtype=ref_step.dtype, device=ref_step.device
                )
            )
            candidate_norm = float(torch.linalg.norm(candidate).detach().cpu())
            finite = bool(torch.isfinite(candidate).all().detach().cpu())
        else:
            def _as_numpy(value):
                if isinstance(value, torch.Tensor):
                    value = value.detach().cpu().numpy()
                return np.asarray(value)

            candidate = _as_numpy(step) + _as_numpy(ip_step)
            candidate_norm = float(np.linalg.norm(candidate))
            finite = bool(np.isfinite(candidate).all())

        if finite and candidate_norm <= self.trust_radius * (1.0 + 1e-12):
            return candidate

        reason = (
            "non-finite"
            if not finite
            else f"outside trust radius ({candidate_norm:.6f} > "
            f"{self.trust_radius:.6f})"
        )
        self.log(
            "Rejected the interpolated/GDIIS displacement because it is "
            f"{reason}; using the restricted reference step."
        )
        return ref_step

    def optimize(self):
        energy, gradient, H, big_eigvals, big_eigvecs, resetted = self.housekeeping()
        step_func, pred_func = self.get_step_func(big_eigvals, gradient)

        forces_act = self.active_list(self.forces)
        coords_act = self.active_list(self.coords)
        steps_act = self.active_list(self.steps)

        ref_gradient = gradient.copy() if isinstance(gradient, np.ndarray) else gradient.clone()
        # Reference step, used for judging the proposed GDIIS step
        ref_step = step_func(big_eigvals, big_eigvecs, gradient) # heavy-compute

        # Right everything is in place to check for convergence.  If all values are below
        # the thresholds, there is no need to do additional inter/extrapolations.
        # allow_stall=False: this is a provisional probe of a proposed step, not
        # the run-loop's terminal check, so a plateau here must not stall.
        if self.check_convergence(ref_step, allow_stall=False)[0]:  # Drop conv_info
            self.log("Convergence achieved! Skipping inter/extrapolation.")
            return ref_step

        # Try to interpolate an intermediate geometry, either from GDIIS or line search.
        #
        # Set some defaults
        ip_gradient = None
        ip_step = None
        diis_result = None

        # Check if we can do GDIIS or GEDIIS. If we (can) do a line search is decided
        # after trying GDIIS.
        rms_forces = rms(gradient)
        rms_step = rms(ref_step)
        can_diis = (rms_step <= self.gdiis_thresh) and (not resetted)
        can_gediis = (rms_forces <= self.gediis_thresh) and (not resetted)

        # GDIIS / GEDIIS, prefer GDIIS over GEDIIS
        if self.gdiis and can_diis:
            # Gradients as error vectors
            if isinstance(ref_step, torch.Tensor):
                err_vecs = -torch.from_numpy(np.array(forces_act)).to(ref_step.dtype).to(ref_step.device)
            else:
                err_vecs = -np.array(forces_act)
            diis_result = gdiis(
                err_vecs,
                coords_act,
                forces_act,
                ref_step,
                test_direction=self.gdiis_test_direction,
            )
            self.successful_gdiis += 1 if diis_result else 0
        # Don't try GEDIIS if GDIIS failed. If GEDIIS should be tried after GDIIS failed
        # comment the line below and uncomment the line following it.
        elif self.gediis and can_gediis:
            # if self.gediis and can_gediis and (diis_result == None):
            diis_result = gediis(coords_act, self.energies, forces_act, hessian=H)
            self.successful_gediis += 1 if diis_result else 0

        try:
            ip_coords = diis_result.coords
            if isinstance(ip_coords, torch.Tensor):
                ip_step = ip_coords - torch.from_numpy(coords_act[-1]).to(ip_coords.device, ip_coords.dtype)
            else:
                ip_step = ip_coords - coords_act[-1]
            ip_gradient = -diis_result.forces
        # When diis_result is None
        except AttributeError:
            self.log("GDIIS didn't succeed.")

        # Try line search if GDIIS failed or not requested
        if self.line_search and (diis_result is None) and (not resetted):
            ip_energy, ip_gradient, ip_step = poly_line_search(
                energy,
                self.energies[-2],
                gradient,
                -forces_act[-2],
                steps_act[-1],
                cubic_max_x=-1,
                quartic_max_x=2,
                logger=self.logger,
            )
            self.successful_line_search += 1 if ip_gradient is not None else 0

        # Use the interpolated gradient for the RFO step if interpolation succeeded
        if (ip_gradient is not None) and (ip_step is not None):
            gradient = ip_gradient
            step = step_func(big_eigvals, big_eigvecs, gradient) # heavy-compute
            # The trust region constrains the full displacement, including
            # the interpolation/GDIIS offset.
            step = self._accept_accelerated_step(step, ip_step, ref_step)
        # Keep the original gradient when the interpolation failed; reuse ref_step.
        else:
            step = ref_step

        # Use the original, actually calculated, gradient
        prediction = pred_func(ref_gradient, H, step)
        self.predicted_energy_changes.append(prediction)

        step = self.full_from_active(step)
        if isinstance(step, torch.Tensor):
            step = step.cpu().numpy()
        return step

    def postprocess_opt(self):
        if self.rejected_uphill_steps or self.skipped_bfgs_updates:
            print(
                "RFO safeguards: "
                f"rejected {self.rejected_uphill_steps} uphill trial(s), "
                f"skipped {self.skipped_bfgs_updates} unsafe BFGS update(s)."
            )
        if self.uphill_rejection_stalled:
            if self.is_converged:
                print(
                    "RFO reached the rejected-trial trust floor; "
                    "the retained lower-energy geometry satisfies the "
                    "convergence criteria."
                )
            else:
                print(
                    "RFO stopped at the rejected-trial trust floor; "
                    "the previous lower-energy geometry was retained without "
                    "claiming convergence."
                )
        msg = (
            f"Successful invocations:\n"
            f"\t     GEDIIS: {self.successful_gediis}\n"
            f"\t      GDIIS: {self.successful_gdiis}\n"
            f"\tLine Search: {self.successful_line_search}\n"
        )
        self.log(msg)
