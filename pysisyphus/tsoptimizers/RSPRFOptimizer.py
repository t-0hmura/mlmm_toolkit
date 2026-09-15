# See [1] https://pubs.acs.org/doi/pdf/10.1021/j100247a015
#         Banerjee, 1985
#     [2] https://aip.scitation.org/doi/abs/10.1063/1.2104507
#         Heyden, 2005
#     [3] https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540070402
#         Baker, 1985
#     [4] https://link.springer.com/article/10.1007/s002140050387
#         Besalu, 1998


import numpy as np
from scipy.optimize import root_scalar

from pysisyphus._array import as_numpy
from pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer


class RSPRFOptimizer(TSHessianOptimizer):
    supports_max_atom_trust = True
    def _defer_hosp_terminal_check(self, step):
        """Keep the physical refresh cadence while the model still rejects TS."""
        if (
            self.flatten_enabled
            or not self.verify_saddle
            or self._saddle_recovery_active
            or self.stop_requested
            or self._last_exact_validation != "higher_order"
            or self._last_exact_n_negative is None
            or self._last_exact_n_negative <= len(self.roots)
            or self._last_exact_frequencies_cm is None
            or not np.all(np.isfinite(self._last_exact_frequencies_cm))
            or self._last_exact_cart_coords is None
            or np.shape(self._last_exact_cart_coords) != np.shape(self.geometry.cart_coords)
            or self.hessian_recalc is None
            or not np.isfinite(self.hessian_recalc)
            or self.hessian_recalc <= 0
            or self.hessian_xtb
            or self._exact_phva_matches_current_geometry()
            or not self._all_configured_values_met(step)
        ):
            return False

        exact_projection = getattr(self, "_last_rigid_projection_info", None)
        from pysisyphus.normal_modes import _strict_negative_count
        if (
            not exact_projection
            or _strict_negative_count(self._last_exact_frequencies_cm, exact_projection)
            != self._last_exact_n_negative
        ):
            return False
        try:
            frequency_data = self._mw_frequencies_and_modes()
            model_projection = self._last_rigid_projection_info
            model_n_negative = (
                _strict_negative_count(frequency_data[0], model_projection)
                if frequency_data is not None else None
            )
        except Exception as err:
            self.log(f"Model PHVA screen unavailable; retaining exact check: {err}")
            return False
        finally:
            # Model screening must not replace the saved exact PHVA metadata.
            self._last_rigid_projection_info = exact_projection
        if frequency_data is None:
            return False
        scope_keys = (
            "active_atoms", "frozen_atoms", "treatment", "frequency_zero_cutoff_cm",
        )
        if any(
            key not in exact_projection or key not in model_projection
            or exact_projection[key] != model_projection[key]
            for key in scope_keys
        ):
            return False
        return bool(
            model_n_negative is not None and model_n_negative > len(self.roots)
        )

    def _exact_terminal_candidate_matches_current_geometry(self):
        if self.flatten_enabled:
            return super()._exact_terminal_candidate_matches_current_geometry()
        return self._exact_saddle_matches_current_geometry()

    def _image_trust_step(self):
        """Build a bounded image step from the current physical model."""
        gradient, _, eigvals, eigvecs = self._hessian_system(
            -np.asarray(self.forces[-1])
        )
        self.update_ts_mode(eigvals, eigvecs)
        eigvals, eigvecs, gradient = map(as_numpy, (eigvals, eigvecs, gradient))
        image_eigvals = eigvals.copy()
        image_eigvals[self.roots] *= -1
        image_gradient = eigvecs.T @ gradient
        image_gradient[self.roots] *= -1
        step = self.get_newton_step_on_trust(
            image_eigvals, eigvecs, eigvecs @ image_gradient
        )
        self.table.print(
            "RS-P-RFO: continuing with a restricted image-quadratic step."
        )
        return step, gradient

    @staticmethod
    def _partition_dstep2_dalpha(alpha, eigval, step, eigvals, gradient):
        """Derivative of a partitioned squared RFO step (Besalú Eq. 18)."""
        step2 = float(np.dot(step, step))
        if step2 == 0.0:
            return 0.0
        denom = (eigvals - eigval * alpha) ** 3
        return (
            2.0
            * eigval
            / (1.0 + step2 * alpha)
            * np.sum(gradient**2 / denom)
        )

    def _max_atom_prfo_step(
        self, eigvals, eigvecs, gradient_trans, ip_step_trans,
        max_indices, min_indices,
    ):
        """Select the existing PRFO scalar family by its Cartesian atom norm.

        No physical roots or Hessian signs are changed. The norm includes the
        line-search contribution. Fixed incoming augmented-root hints make
        bracket evaluations independent of the order in which they occur.
        This is a scalar restriction, not a product-of-balls quadratic solve.
        """
        radius = float(self.trust_radius)
        alpha0 = float(self.alpha0)
        if not np.isfinite(alpha0) or alpha0 <= 0:
            raise ValueError("max_atom PRFO requires a finite positive alpha0")
        cache = {}
        incoming = (self.prev_eigvec_max, self.prev_eigvec_min)

        class EvaluationBudget(Exception):
            pass

        def evaluate(alpha):
            alpha = float(alpha)
            if alpha in cache:
                return cache[alpha]
            if len(cache) >= self.max_micro_cycles:
                raise EvaluationBudget
            step = np.zeros_like(gradient_trans)
            returned = []
            for indices, kind, previous in zip(
                (max_indices, min_indices), ("max", "min"), incoming
            ):
                result = self.solve_rfo_secular(
                    eigvals[indices], gradient_trans[indices], alpha,
                    kind=kind, prev_eigvec=previous,
                )
                if result is None:
                    augmented = self.get_augmented_hessian(
                        eigvals[indices], gradient_trans[indices], alpha
                    )
                    result = self.solve_rfo(
                        augmented, kind, prev_eigvec=previous, alpha=alpha
                    )
                step[indices] = result[0]
                returned.append(result[3])
            cart = eigvecs @ (step + ip_step_trans)
            norm = self._trust_step_norm(cart)
            cache[alpha] = (step, norm, returned)
            return cache[alpha]

        if not np.any(gradient_trans):
            # Preserve the stationary-candidate terminal PHVA/recovery owner.
            self._last_atomic_trust = dict(
                alpha=alpha0, evaluations=0, termination="stationary",
                max_atom_bohr=self._trust_step_norm(eigvecs @ ip_step_trans),
                global_l2_bohr=float(np.linalg.norm(eigvecs @ ip_step_trans)),
            )
            return np.zeros_like(gradient_trans)
        step, norm, returned = evaluate(alpha0)
        alpha = alpha0
        termination = "interior"
        if norm > radius:
            lower, upper = alpha0, alpha0
            try:
                while norm > radius:
                    lower = upper
                    upper *= 2.0
                    if not np.isfinite(upper):
                        raise ValueError("max_atom PRFO could not bracket a finite restriction")
                    _, norm, _ = evaluate(upper)
                solution = root_scalar(
                    lambda value: evaluate(value)[1] - radius,
                    bracket=(lower, upper), method="brentq",
                    xtol=np.finfo(float).smallest_subnormal,
                    rtol=4 * np.finfo(float).eps,
                    maxiter=self.max_micro_cycles,
                )
                alpha = solution.root
                step, norm, returned = evaluate(alpha)
                termination = "boundary" if solution.converged else "feasible_budget"
            except EvaluationBudget:
                termination = "feasible_budget"
            # Keep an actually evaluated inside point if the scalar budget
            # is exhausted or the boundary rounds outside. No monotonicity
            # assumption or radial scaling is needed for this finite bound.
            if termination != "boundary" or norm > radius:
                feasible = [(value, row) for value, row in cache.items()
                            if row[1] <= radius]
                if not feasible:
                    raise ValueError("max_atom PRFO exhausted its budget without a feasible point")
                alpha, (step, norm, returned) = max(feasible, key=lambda item: item[1][1])
                if termination == "boundary":
                    termination = "feasible_roundoff"
        if norm > radius * (1.0 + 1e-12):
            raise ValueError("max_atom PRFO scalar solution exceeds its atomic bound")
        self.prev_eigvec_max, self.prev_eigvec_min = returned
        self._last_atomic_trust = dict(
            alpha=float(alpha), evaluations=len(cache), max_atom_bohr=float(norm),
            global_l2_bohr=float(np.linalg.norm(eigvecs @ (step + ip_step_trans))),
            termination=termination,
        )
        self.log(
            f"max_atom PRFO: alpha={alpha:.8g}, evaluations={len(cache)}, "
            f"max_atom={norm:.8g} Bohr, termination={termination}"
        )
        return step

    def optimize(self):
        energy, gradient, H, eigvals, eigvecs, resetted = self.housekeeping()
        self.update_ts_mode(eigvals, eigvecs)
        exact_negative_count = self._last_exact_n_negative
        if self._last_exact_frequencies_cm is None:
            # Preserve the legacy non-Cartesian step policy, not a PHVA proof.
            exact_negative_count = self._last_exact_n_imaginary

        # RS-PRFO uses np.linalg.norm and scalar Python loops; coerce torch
        # tensors from the MLIP Hessian path to numpy so the .dot and
        # fancy-indexing operations below stay valid.
        eigvals = as_numpy(eigvals)
        eigvecs = as_numpy(eigvecs)
        gradient = as_numpy(gradient)

        if (
            self._physical_ts_mode is not None
            and not (
                not self.flatten_enabled
                and exact_negative_count is not None
                and exact_negative_count > len(self.roots)
            )
        ):
            # A past PHVA result does not classify today's complementary
            # curvature. Only pure unconstrained translations are artifacts;
            # rotations at nonstationary points and mixed modes can be physical.
            residual_negative = (
                (eigvals < -self.small_eigval_thresh)
                & self._translation_mode_mask(eigvecs)
            )
            residual_negative[self.roots] = False
            residual_count = int(np.count_nonzero(residual_negative))
            if residual_count:
                eigvals = np.where(residual_negative, -eigvals, eigvals)
                self.log(
                    "Stabilized "
                    f"{residual_count} negative translational minimizing root(s)."
                )

        # Transform gradient to eigensystem of hessian
        gradient_trans = eigvecs.T.dot(gradient)
        # Minimize energy along all modes, except the TS-mode
        min_indices = [i for i in range(gradient_trans.size) if i not in self.roots]
        # Maximize energy along all requested modes.
        max_indices = [i for i in range(gradient_trans.size) if i in self.roots]
        # Get line search steps, if requested.
        ip_step_trans, gradient_trans = self.step_and_grad_from_line_search(
            energy,
            gradient_trans,
            eigvecs,
            min_indices,
            max_indices,
        )

        """In the RS-(P)RFO method we have to scale the matrices with alpha.
        Unscaled matrix (Eq. 8) in [1]:
            (H  g) (x)          (S 0) (x)
                       = lambda
            (g+ 0) (1)          (0 1) (1)
        with
            S = alpha * Identity matrix
        and multiplying from the left with the inverse of the scaling matrix
            (1/alpha 0)

            (0       1)
        we get
            (1/alpha 0) (H  g) (x)          (x)
                                   = lambda
            (0       1) (g+ 0) (1)          (1)
        eventually leading to the scaled matrix:
            (H/alpha  g/alpha) (x)          (x)
                                   = lambda     .
            (g+             0) (1)          (1)
        """

        alpha = self.alpha0
        image_step = False
        atomic_trust = getattr(self, "trust_norm", "l2") == "max_atom"
        if self.max_micro_cycles < 1:
            raise ValueError("RS-P-RFO requires at least one micro cycle.")
        if atomic_trust:
            try:
                step = self._max_atom_prfo_step(
                    eigvals, eigvecs, gradient_trans, ip_step_trans,
                    max_indices, min_indices,
                )
            except ZeroDivisionError:
                step, gradient = self._image_trust_step()
                image_step = True
        for mu in range(0 if atomic_trust else self.max_micro_cycles):
            self.log(f"RS-PRFO micro cycle {mu:02d}, alpha={alpha:.6f}")

            # A stationary candidate belongs to the terminal PHVA/recovery
            # gate below, including a minimum with an uncoupled uphill root.
            if not np.any(gradient_trans):
                step = np.zeros_like(gradient_trans)
                break

            # Maximize energy along the chosen TS modes.
            # Try secular equation solver first (O(N) vs O(N^3))
            secular_max = self.solve_rfo_secular(
                eigvals[max_indices], gradient_trans[max_indices], alpha,
                kind="max", prev_eigvec=self.prev_eigvec_max,
            )
            if secular_max is not None:
                step_max, eigval_max, nu_max, self.prev_eigvec_max = secular_max
            else:
                H_aug_max = self.get_augmented_hessian(
                    eigvals[max_indices], gradient_trans[max_indices], alpha
                )
                try:
                    step_max, eigval_max, nu_max, self.prev_eigvec_max = self.solve_rfo(
                        H_aug_max, "max", prev_eigvec=self.prev_eigvec_max, alpha=alpha
                    )
                except ZeroDivisionError:
                    step, gradient = self._image_trust_step()
                    image_step = True
                    break

            # Minimize energy along all modes, but the TS mode.
            secular_min = self.solve_rfo_secular(
                eigvals[min_indices], gradient_trans[min_indices], alpha,
                kind="min", prev_eigvec=self.prev_eigvec_min,
            )
            if secular_min is not None:
                step_min, eigval_min, nu_min, self.prev_eigvec_min = secular_min
            else:
                H_aug_min = self.get_augmented_hessian(
                    eigvals[min_indices], gradient_trans[min_indices], alpha
                )
                try:
                    step_min, eigval_min, nu_min, self.prev_eigvec_min = self.solve_rfo(
                        H_aug_min, "min", prev_eigvec=self.prev_eigvec_min, alpha=alpha
                    )
                except ZeroDivisionError:
                    step, gradient = self._image_trust_step()
                    image_step = True
                    break

            # Calculate overlap between directions over the course of the micro cycles
            # if mu == 0:
            # TODO: convert back to original space
            # ref_step_max = step_max.copy()
            # ref_step_min = step_min.copy()
            min_norm = np.linalg.norm(step_min)
            max_norm = np.linalg.norm(step_max)
            self.log(f"norm(step_max)={max_norm:.6f}")
            self.log(f"norm(step_min)={min_norm:.6f}")
            norm_ratio = max_norm / min_norm if min_norm > 0.0 else float("inf")
            self.log(f"norm(step_max)/norm(step_min)={norm_ratio:.2%}")
            # Calculate overlaps with originally proposed step in mu == 0
            # TODO: convert back to original space
            # max_ovlp = ref_step_max @ step_max
            # min_ovlp = ref_step_min @ step_min

            # As of Eq. (8a) of [4] max_eigval and min_eigval also
            # correspond to:
            # max_eigval = -forces_trans[max_indices].dot(max_step)
            # min_eigval = -forces_trans[min_indices].dot(min_step)

            # Create the full PRFO step
            step = np.zeros_like(gradient_trans)
            step[max_indices] = step_max
            step[min_indices] = step_min
            step_norm = np.linalg.norm(step)
            self.log(f"norm(step)={step_norm:.6f}")

            if not np.isfinite(step_norm):
                raise ValueError("RS-P-RFO produced a nonfinite step.")
            # Match the final displacement check's floating-point allowance.
            inside_trust = step_norm <= self.trust_radius * (1.0 + 1e-12)
            if inside_trust:
                self.log(
                    "Restricted step satisfies trust radius of "
                    f"{self.trust_radius:.6f}"
                )
                self.log(
                    f"Micro-cycles converged in cycle {mu:02d} with "
                    f"alpha={alpha:.6f}!"
                )
                break

            # One micro cycle deliberately requests unscaled P-RFO.
            if self.max_micro_cycles == 1:
                break
            if mu + 1 == self.max_micro_cycles:
                raise ValueError(
                    "RS-P-RFO exhausted its micro cycles outside the trust radius."
                )

            # Derivative of the squared step w.r.t. alpha for both
            # partitioned subspaces (Besalú and Bofill, Eq. 18).
            dstep2_dalpha_max = self._partition_dstep2_dalpha(
                alpha,
                eigval_max,
                step_max,
                eigvals[max_indices],
                gradient_trans[max_indices],
            )
            dstep2_dalpha_min = self._partition_dstep2_dalpha(
                alpha,
                eigval_min,
                step_min,
                eigvals[min_indices],
                gradient_trans[min_indices],
            )
            dstep2_dalpha = dstep2_dalpha_max + dstep2_dalpha_min
            if not np.isfinite(dstep2_dalpha) or dstep2_dalpha == 0.0:
                raise ValueError("RS-P-RFO alpha derivative is zero or nonfinite.")
            alpha_step = (
                2 * (self.trust_radius * step_norm - step_norm**2) / dstep2_dalpha
            )
            next_alpha = alpha + alpha_step
            if not np.isfinite(next_alpha) or next_alpha <= 0.0:
                raise ValueError("RS-P-RFO alpha update is not finite and positive.")
            alpha = next_alpha

        # Right now the step is still given in the Hessians eigensystem. We
        # transform it back now.
        if not image_step:
            step += ip_step_trans
            step = eigvecs.dot(step)
        step_norm = np.linalg.norm(step)
        if atomic_trust:
            step_norm = self._trust_step_norm(step)
        if not np.isfinite(step_norm):
            raise ValueError("RS-P-RFO combined step is nonfinite.")

        # With max_micro_cycles = 1 the RS part is disabled and the step
        # probably isn't scaled correctly in the one micro cycle.
        # In this case we use a naive scaling if the step is too big.
        if (self.max_micro_cycles == 1) and (step_norm > self.trust_radius):
            step = step / step_norm * self.trust_radius
        elif step_norm > self.trust_radius * (1.0 + 1e-12):
            raise ValueError("RS-P-RFO combined step exceeds the trust radius.")
        self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        # Eq. (6) from [4] seems erronous ... the prediction is usually only ~50%
        # of the actual change ...
        # predicted_energy_change = 1/2 * (eigval_max / nu_max**2 + eigval_min / nu_min**2)
        # self.predicted_energy_changes.append(predicted_energy_change)

        deferred_hosp_check = self._defer_hosp_terminal_check(step)
        if not deferred_hosp_check:
            self.validate_terminal_saddle_for_step(step)
        exact_negative_count = (
            self._last_exact_n_negative if self._last_exact_frequencies_cm is not None
            else self._last_exact_n_imaginary
        )
        if (
            not self.stop_requested
            and not self.flatten_enabled
            and (
                deferred_hosp_check
                or (
                    self._exact_phva_matches_current_geometry()
                    and exact_negative_count is not None
                    and exact_negative_count > len(self.roots)
                )
            )
        ):
            # Exact PHVA exposed complementary negative curvature that the
            # current physical model still carries.
            # Take one bounded image-quadratic trust step (as in TRIM), not a
            # finite-ratio P-RFO step at an uncoupled augmented root. Keep the
            # physical Hessian/gradient for updates and energy prediction.
            step, gradient = self._image_trust_step()
            image_step = True
        step = self.apply_saddle_recovery_step(step)
        if atomic_trust:
            # Recovery/terminal fallback can replace the initial proposal.
            # Bound the actual combined Cartesian proposal before prediction.
            step = self._bound_to_trust_radius(step, label="combined RS-PRFO step")
        prediction = self.quadratic_model if image_step else self.rfo_model
        self.predicted_energy_changes.append(prediction(gradient, as_numpy(self.cur_H), step))

        self.log("")
        step = self.full_from_active(step)
        return step
