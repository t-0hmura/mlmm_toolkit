# [1  ] https://aip.scitation.org/doi/pdf/10.1063/1.3514202?class=pdf
#       Original EulerPC
#       Hratchian, Schlegel, 2010
# [2  ] https://aip.scitation.org/doi/pdf/10.1063/1.1724823?class=pdf
#       Original HPC
#       Hratchian, Schlegel, 2004
# [3  ] https://pubs.rsc.org/en/content/articlepdf/2017/cp/c7cp03722h
#       EulerPC re-implementation
#       Meisner, Kästner, 2017
# [3.1] http://www.rsc.org/suppdata/c7/cp/c7cp03722h/c7cp03722h1.pdf
#       Corresponding SI
# [4  ] https://aip.scitation.org/doi/pdf/10.1063/1.3593456?class=pdf
#       Hratchian, Frisch, 2011
# [6  ] https://aip.scitation.org/doi/10.1063/1.3593456<Paste>
#       Hratchian, Frisch
# 	Further improvements for DWI; not implemented

import time
from collections import deque

import numpy as np

from pysisyphus._array import active_square
from pysisyphus.helpers_pure import rms
from pysisyphus.irc.DWI import DWI
from pysisyphus.irc.IRC import IRC
from pysisyphus.optimizers.hessian_updates import (
    bfgs_update,
    bofill_cpu_offload_enabled,
    bofill_rank2_factors,
    bofill_update,
)

import torch


class EulerPC(IRC):
    def __init__(
        self,
        *args,
        hessian_recalc=None,
        hessian_update="bofill",
        hessian_init="calc",
        max_pred_steps=500,
        loose_cycles=3,
        corr_func="mbs",
        **kwargs,
    ):
        super().__init__(*args, hessian_init=hessian_init, **kwargs)

        self.hessian_recalc = hessian_recalc
        self.hessian_update = {
            "bfgs": bfgs_update,
            "bofill": bofill_update,
        }
        self.hessian_update_func = self.hessian_update[hessian_update]
        self.max_pred_steps = int(max_pred_steps)
        self.loose_cycles = loose_cycles

        corr_funcs = {
            "mbs": self.corrector_step,
        }
        self.corr_func = corr_funcs[corr_func]

        self._to_active_vec = lambda v: v[self._act_dofs]

        def _to_full_vec(v_act, template):
            if isinstance(v_act, torch.Tensor):
                full = template.clone()
            else:
                full = template.copy()
            full[self._act_dofs] = v_act
            return full

        self._to_full_vec = _to_full_vec
        
    def prepare(self, *args, **kwargs):
        super().prepare(*args, **kwargs)

        # Initialize the distance weighted interpolator with the data
        # from the initial displacement.
        self.dwi = DWI()
        mw_grad_full = self.mw_gradient
        energy = self.energy

        mw_grad_act   = self._to_active_vec(mw_grad_full)
        mw_coords_act = self._to_active_vec(self.mw_coords)

        self.dwi.update(
            mw_coords_act,
            energy,
            mw_grad_act,
            self.mw_hessian,
            copy_hessian=True,
        )

        if self.downhill:
            return

        # Do a first Hessian update with the information between the TS
        # and the initially displaced geometry.
        dx = self._to_active_vec(self.mw_coords - self.ts_mw_coords)
        dg = mw_grad_act - self._to_active_vec(self.ts_mw_gradient)
        key = self._apply_hessian_update(dx, dg)
        self.log(f"Did {key} hessian update.")

    def _apply_hessian_update(self, dx, dg):
        """Apply a quasi-Newton update without retaining a dense temporary."""
        if (
            isinstance(self.mw_hessian, torch.Tensor)
            and self.hessian_update_func is bofill_update
            and not bofill_cpu_offload_enabled()
        ):
            U, C = bofill_rank2_factors(self.mw_hessian, dx, dg)
            self.mw_hessian.addmm_(U, C @ U.t())
            return "Bofill"

        dH, key = self.hessian_update_func(self.mw_hessian, dx, dg)
        if isinstance(self.mw_hessian, torch.Tensor):
            self.mw_hessian.add_(dH)
        else:
            self.mw_hessian += dH
        del dH
        return key

    def get_integration_length_func(self, init_mw_coords):
        """
        Return a closure that computes the un-mass-weighted path length Δs.
        """
        if init_mw_coords.shape[0] == self._m_sqrt.shape[0]:
            m_sqrt_vec = self._m_sqrt                    # full (3N)
        else:
            m_sqrt_vec = self._m_sqrt[self._act_dofs]    # active-only

        def get_integration_length(cur_mw_coords):
            return np.linalg.norm((cur_mw_coords - init_mw_coords) / m_sqrt_vec)

        return get_integration_length

    def step(self):
        ##################
        # PREDICTOR STEP #
        ##################

        mw_grad_full = self.mw_gradient            # 3N-vector, mass-weighted
        energy       = self.energy
        mw_grad_act  = self._to_active_vec(mw_grad_full)  # active slice

        if self.cur_cycle > 0:
            if self.hessian_recalc and (self.cur_cycle % self.hessian_recalc == 0):
                H_new = self.geometry.hessian
                act_n = len(self._act_dofs)
                if H_new.shape[0] != act_n:
                    H_new = active_square(H_new, self._act_dofs)
                _old_mw_hessian = self.mw_hessian
                self.mw_hessian = None
                del _old_mw_hessian
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                self.mw_hessian = self._mw_hessian_active(H_new)
                self.geometry.clear()
                self.log("Calculated exact hessian")
            else:
                dx = self._to_active_vec(self.mw_coords - self.irc_mw_coords[-2])
                dg = mw_grad_act - self._to_active_vec(self.irc_mw_gradients[-2])
                key = self._apply_hessian_update(dx, dg)
                self.log(f"Did {key} hessian update before predictor step.")
            self.dwi.update(
                self._to_active_vec(self.mw_coords).copy(),
                energy,
                mw_grad_act,
                self.mw_hessian,
                copy_hessian=True,
            )

        # Create a copy of the inital coordinates for the determination
        # of the actual step size in the predictor Euler integration.
        init_mw_coords = self.mw_coords.copy()

        get_integration_length = self.get_integration_length_func(init_mw_coords)

        # Calculate predictor Euler-integration step length. See get_conv_fact
        # method definition for a comment on this.
        conv_fact = self.get_conv_fact(mw_grad_act)
        euler_step_length = self.step_length / (self.max_pred_steps / conv_fact)

        def taylor_gradient(step_full):
            """Return full-length gradient via 2nd-order Taylor (active Hessian)."""
            step_act = self._to_active_vec(step_full)
            if isinstance(self.mw_hessian, torch.Tensor):
                step_t  = torch.tensor(step_act, dtype=self.mw_hessian.dtype,
                                       device=self.mw_hessian.device)
                hvp_act = (self.mw_hessian @ step_t).cpu().numpy()
            else:
                hvp_act = self.mw_hessian @ step_act
            grad_act = mw_grad_act + hvp_act
            return self._to_full_vec(grad_act, step_full)

        # These variables will hold the coordinates and gradients along
        # the Euler integration and will be updated frequently.
        euler_mw_coords = self.mw_coords.copy()
        euler_mw_grad   = mw_grad_full.copy()
        self.log(
            f"Predictor-Euler-integration with Δs={euler_step_length:.6f} "
            f"for up to {self.max_pred_steps} steps\n     #  |step|  d|step|"
        )
        prev_cur_length = 0.0
        for i in range(self.max_pred_steps):
            # Calculate step length in non-mass-weighted coordinates
            cur_length = get_integration_length(euler_mw_coords)
            if i % 50 == 0:
                diff = cur_length - prev_cur_length
                self.log(f"\t{i:03d}: {cur_length:.4f} Δ={diff:.4f}")
                prev_cur_length = cur_length

            # Check if we achieved the desired step length.
            if cur_length >= self.step_length:
                self.log(
                    "Predictor-Euler integration converged with "
                    f"Δs={cur_length:.4f} (desired Δs={self.step_length:.4f}) "
                    f"after {i+1} steps!"
                )
                break
            grad_norm = np.linalg.norm(euler_mw_grad)
            if not np.isfinite(grad_norm) or grad_norm == 0.0:
                self.log("Gradient norm is zero/NaN; using transition vector direction.")
                direction = getattr(self, "mw_transition_vector", euler_mw_grad)
                dir_norm = np.linalg.norm(direction)
                if not np.isfinite(dir_norm) or dir_norm == 0.0:
                    raise ValueError("Cannot determine IRC step direction (zero/NaN norms).")
                step_ = euler_step_length * -direction / dir_norm
            else:
                step_ = euler_step_length * -euler_mw_grad / grad_norm
            euler_mw_coords += step_
            # Determine actual step by comparing the current and the initial coordinates
            euler_step = euler_mw_coords - init_mw_coords
            euler_mw_grad = taylor_gradient(euler_step)
        else:
            # The final permitted microstep is applied at the end of the loop,
            # after the last pre-step convergence check.  Accept it when it
            # reached the requested integration length instead of reporting a
            # false exhaustion.
            cur_length = get_integration_length(euler_mw_coords)
            if cur_length >= self.step_length:
                self.log(
                    "Predictor-Euler integration converged with "
                    f"Δs={cur_length:.4f} (desired Δs={self.step_length:.4f}) "
                    f"after {self.max_pred_steps} steps!"
                )
            else:
                self.log(
                    f"Predictor-Euler integration did not converge in {i+1} "
                    f"steps. Δs={cur_length:.4f}."
                )

                self.mw_coords = euler_mw_coords

                if self.cur_cycle < self.loose_cycles:
                    self.log(
                        f"Current cycle {self.cur_cycle} is still in 'loose' mode.\n"
                        "Continuing IRC integration even though predictor integration "
                        f"did not succeed.\n{self.loose_cycles - self.cur_cycle - 1} loose "
                        "cycles remaining."
                    )
                else:
                    self.integration_stop_requested = True
                    self.integration_stop_reason = (
                        "Predictor integration exhausted before reaching its "
                        "requested step length."
                    )
                    return
        self.log("")

        # Calculate energy and gradient at new predicted geometry. Update the
        # hessian accordingly. These results will be added to the DWI for use
        # in the corrector step.
        self.mw_coords = euler_mw_coords
        self.log("Calculating energy and gradient at predictor step geometry.")
        mw_grad_full = self.mw_gradient
        energy       = self.energy
        mw_grad_act  = self._to_active_vec(mw_grad_full)

        # Hessian update
        dx = self._to_active_vec(self.mw_coords - self.irc_mw_coords[-1])
        dg = mw_grad_act - self._to_active_vec(self.irc_mw_gradients[-1])
        key = self._apply_hessian_update(dx, dg)
        self.log(f"Did {key} hessian update after predictor step.\n")
        self.dwi.update(
            self._to_active_vec(self.mw_coords).copy(),
            energy,
            mw_grad_act,
            self.mw_hessian,
            copy_hessian=True,
        )

        corrected_mw_coords_act = self.corr_func(self._to_active_vec(init_mw_coords), self.step_length, self.dwi)
        self.mw_coords = self._to_full_vec(corrected_mw_coords_act, self.mw_coords)
        corr_step_length = get_integration_length(self.mw_coords)
        self.log(f"Corrected unweighted step length: {corr_step_length:.6f}")

    def corrector_step(self, init_mw_coords, step_length, dwi):
        self.log("Corrector step using mBS integration (active dofs)")

        get_integration_length = self.get_integration_length_func(init_mw_coords)

        errors = list()
        richardson = dict()
        for k in range(15):
            points = 20 * (2**k)
            corr_step_length = step_length / (points - 1)
            cur_coords = init_mw_coords.copy()
            # Only keep the last 2 coords (needed for oscillation check)
            k_coords = deque(maxlen=2)
            cur_length = 0

            # A q-space microstep advances the unweighted path by a
            # mass-dependent amount. Bound by the heaviest active component,
            # then fail closed if the target still cannot be reached.
            m_sqrt_active = self._m_sqrt[self._act_dofs]
            max_euler_steps = max(
                100,
                int(np.ceil(3.0 * np.max(m_sqrt_active) * (points - 1))),
            )
            reached_target = False
            for _euler_step in range(max_euler_steps):
                k_coords.append(cur_coords.copy())
                if abs(step_length - cur_length) < 0.5 * corr_step_length:
                    self.log(
                        f"\tk={k:02d} points={points: >4d} "
                        f"step_length={corr_step_length:.4f} Δs={cur_length:.4f}"
                    )
                    reached_target = True
                    break

                energy, gradient = dwi.interpolate(cur_coords, gradient=True)
                grad_norm = np.linalg.norm(gradient)
                if not np.isfinite(grad_norm) or grad_norm < 1e-10:
                    raise RuntimeError(
                        "Corrector integration encountered a zero or "
                        f"non-finite gradient at level {k}."
                    )
                cur_coords += corr_step_length * -gradient / grad_norm
                # cur_length += corr_step_length
                cur_length = get_integration_length(cur_coords)

                # Check for oscillation
                if len(k_coords) >= 2:
                    prev_coords = k_coords[-2]
                    osc_norm = np.linalg.norm(cur_coords - prev_coords)
                    # TODO: Handle this by restarting everything with a smaller stepsize?
                    # Check 10.1039/c7cp03722h SI
                    if osc_norm <= corr_step_length:
                        # The corrector descends the DWI *interpolated* surface,
                        # not the real PES, so a reversal here is an artefact of
                        # the two-point interpolation, not a physical failure.
                        # Abort only the corrector and keep the last
                        # non-oscillating point, as upstream pysisyphus does:
                        # this branch is also the integration loop's escape
                        # hatch, so raising leaves the step budget below as the
                        # only exit and kills an otherwise healthy IRC.
                        msg = (
                            "Corrector-Euler integration oscillated at level "
                            f"{k} ({points} points); keeping the last "
                            "non-oscillating point for this IRC step."
                        )
                        print(f"WARNING: {msg}")
                        self.log(f"\t{msg}")
                        return prev_coords
            if not reached_target:
                raise RuntimeError(
                    "Corrector integration exhausted its mass-aware step "
                    f"budget at level {k} before reaching the target."
                )
            richardson[(k, 0)] = cur_coords

            # Refine using Richardson extrapolation
            # Set additional values using Richard extrapolation
            for j in range(1, k + 1):
                richardson[(k, j)] = (
                    (2**j) * richardson[(k, j - 1)] - richardson[(k - 1, j - 1)]
                ) / (2**j - 1)
            # Can only be done after the second successful integration
            if k > 0:
                # Error estimate according to Numerical Recipes Eq. (17.3.9).
                # We compare the last two entries/columns in the current row.
                # RMS error
                error = np.sqrt(
                    np.mean((richardson[(k, k)] - richardson[(k, k - 1)]) ** 2)
                )
                errors.append(error)
                if error <= 1e-5:
                    self.log(f"mBS integration converged (error={error:.4e})!")
                    break
        else:
            raise Exception("Richardson did not converge!")

        self.log(
            f"Returning corrected mass-weighted coordinates from richardson[({k},{k})]"
        )
        return richardson[(k, k)]
