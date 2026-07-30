from typing import Optional

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.intcoords.exceptions import NeedNewInternalsException
from pysisyphus.optimizers.closures import bfgs_multiply, get_update_mu_reg
from pysisyphus.optimizers.hessian_updates import double_damp
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.restrict_step import scale_by_max_step
from pysisyphus.optimizers.poly_fit import poly_line_search


class LBFGS(Optimizer):
    def __init__(
        self,
        geometry: Geometry,
        keep_last: int = 7,
        beta: float = 1,
        max_step: float = 0.2,
        double_damp: bool = True,
        gamma_mult: bool = False,
        line_search: bool = False,
        mu_reg: Optional[float] = None,
        max_mu_reg_adaptions: int = 10,
        control_step: bool = True,
        reject_uphill: bool = False,
        uphill_tolerance: float = 1e-3,
        rejection_step_floor: float = 1e-7,
        max_rejections_at_floor: int = 3,
        **kwargs,
    ) -> None:
        """Limited-memory BFGS optimizer.

        See [1] Nocedal, Wright - Numerical Optimization, 2006 for a general
        discussion of LBFGS. See pysisyphus.optimizers.hessian_updates for
        the references related to double damping and pysisyphus.optimizers.closures
        for references related to regularized LBFGS.

        Parameters
        ----------
        geometry
            Geometry to be optimized.
        keep_last
            History size. Keep last 'keep_last' steps and gradient differences.
        beta
            Force constant β in -(H + βI)⁻¹g.
        max_step
            Upper limit for the absolute component of the step vector in whatever
            unit the optimization is carried out.
        double_damp
            Use double damping procedure to modify steps s and gradient differences y
            to ensure sy > 0.
        gamma_mult
            Estimate β from previous cycle. Eq. (7.20) in [1]. See 'beta' argument.
        line_search
            Enable implicit linesearches.
        mu_reg
            Initial guess for regularization constant in regularized LBFGS.
        max_mu_reg_adaptions
            Maximum number of trial steps in regularized LBFGS.
        control_step
            Wheter to scale down the proposed step its biggest absolute component
            is equal to or below 'max_step'
        reject_uphill
            Reject energy-increasing minimization trials and retry from the
            previous accepted geometry with a smaller step limit.
        uphill_tolerance
            Positive energy change tolerated before a trial is rejected.
        rejection_step_floor
            Emergency maximum-component step floor for rejected trials.
        max_rejections_at_floor
            Rejections allowed at the emergency floor before stopping at the
            retained lower-energy geometry.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments passed to the Optimizer baseclass.
        """

        self.coord_diffs = list()
        self.grad_diffs = list()
        super().__init__(geometry, max_step=max_step, **kwargs)

        self.beta = beta
        self.keep_last = int(keep_last)
        self.double_damp = double_damp
        self.gamma_mult = gamma_mult
        self.mu_reg = mu_reg
        self.max_mu_reg_adaptions = max_mu_reg_adaptions
        self.line_search = (not self.is_cos) and line_search
        self.control_step = control_step
        # Chain-of-states and dimer-like objectives may legitimately increase
        # the physical energy.  Their callers leave this safeguard disabled.
        self.reject_uphill = bool(reject_uphill) and not self.is_cos
        self.uphill_tolerance = float(uphill_tolerance)
        if self.uphill_tolerance < 0.0:
            raise ValueError("uphill_tolerance must be non-negative")
        self.rejection_step_floor = float(rejection_step_floor)
        if self.rejection_step_floor <= 0.0:
            raise ValueError("rejection_step_floor must be positive")
        self.max_rejections_at_floor = int(max_rejections_at_floor)
        if self.max_rejections_at_floor < 1:
            raise ValueError("max_rejections_at_floor must be at least 1")
        self._trial_max_step = float(max_step)
        self.rejected_uphill_steps = 0
        self.rejections_at_floor = 0

        self.tot_adapt_mu_cycles = 0
        if self.mu_reg:
            self.mu_reg_0 = self.mu_reg  # Backup value
            self.update_mu_reg = get_update_mu_reg(logger=self.logger)
            self.control_step = False
            self.double_damp = False
            self.line_search = False
            self.log(
                f"Regularized-L-BFGS (μ_reg={self.mu_reg:.6f}) requested.\n"
                f"Disabling double damping, step control and line search."
            )

    def reset(self):
        self.coord_diffs = list()
        self.grad_diffs = list()

    # Subclass restart keys required for a bit-exact resume. The uphill-
    # rejection adaptive state (_trial_max_step, rejected_uphill_steps,
    # rejections_at_floor) is mutated every cycle by the reject/relax logic in
    # :meth:`optimize`; omitting it silently reset a resumed run to the
    # __init__ defaults and diverged the trajectory across a rejection
    # boundary.  Those three keys are new, so they are *restored tolerantly*
    # (see :meth:`_set_opt_restart_info`) and are kept out of the transactional
    # required set below so older checkpoints still load.
    required_opt_restart_keys = (
        "coord_diffs",
        "grad_diffs",
        "double_damp",
        "gamma_mult",
        "keep_last",
    )

    def _get_opt_restart_info(self):
        opt_restart_info = {
            "coord_diffs": np.array(self.coord_diffs).tolist(),
            "grad_diffs": np.array(self.grad_diffs).tolist(),
            "double_damp": self.double_damp,
            "gamma_mult": self.gamma_mult,
            "keep_last": self.keep_last,
            # Uphill-rejection adaptive state.
            "_trial_max_step": self._trial_max_step,
            "rejected_uphill_steps": self.rejected_uphill_steps,
            "rejections_at_floor": self.rejections_at_floor,
            # Regularized-L-BFGS adaptive state: mu_reg is adapted per cycle and
            # feeds get_lbfgs_step, so a resumed regularized run diverges without
            # it; tot_adapt_mu_cycles keeps the reported update count accurate.
            "mu_reg": self.mu_reg,
            "tot_adapt_mu_cycles": self.tot_adapt_mu_cycles,
        }
        return opt_restart_info

    def _set_opt_restart_info(self, opt_restart_info):
        self.coord_diffs = [np.array(cd) for cd in opt_restart_info["coord_diffs"]]
        self.grad_diffs = [np.array(gd) for gd in opt_restart_info["grad_diffs"]]
        for attr in ("double_damp", "gamma_mult", "keep_last"):
            setattr(self, attr, opt_restart_info[attr])
        # Backward-tolerant: checkpoints written before the adaptive state
        # was serialized keep the __init__ defaults already assigned above.
        if "_trial_max_step" in opt_restart_info:
            self._trial_max_step = float(opt_restart_info["_trial_max_step"])
        if "rejected_uphill_steps" in opt_restart_info:
            self.rejected_uphill_steps = int(opt_restart_info["rejected_uphill_steps"])
        if "rejections_at_floor" in opt_restart_info:
            self.rejections_at_floor = int(opt_restart_info["rejections_at_floor"])
        # mu_reg may be None (regularization off) or a float; restore verbatim.
        if "mu_reg" in opt_restart_info:
            self.mu_reg = opt_restart_info["mu_reg"]
        if "tot_adapt_mu_cycles" in opt_restart_info:
            self.tot_adapt_mu_cycles = int(opt_restart_info["tot_adapt_mu_cycles"])

    def get_lbfgs_step(self, forces):
        return bfgs_multiply(
            self.coord_diffs,
            self.grad_diffs,
            forces,
            beta=self.beta,
            gamma_mult=self.gamma_mult,
            mu_reg=self.mu_reg,
            logger=self.logger,
        )

    def optimize(self):
        if self.is_cos and self.align:
            rot_vecs, rot_vec_lists, _ = self.fit_rigid(
                vector_lists=(
                    self.steps,
                    self.forces,
                    self.coord_diffs,
                    self.grad_diffs,
                ),
            )
            rot_steps, rot_forces, rot_coord_diffs, rot_grad_diffs = rot_vec_lists
            self.steps = rot_steps
            self.forces = rot_forces
            self.coord_diffs = rot_coord_diffs
            self.grad_diffs = rot_grad_diffs

        forces = self.geometry.forces
        self.forces.append(forces)
        energy = self.geometry.energy
        self.energies.append(energy)
        norm = np.linalg.norm(forces)
        if not self.is_cos:
            self.log(f"      Energy={energy: >24.6f} au")
        self.log(f"norm(forces)={norm: >24.6f} au / bohr (rad)")

        trial_rejected = False
        if (
            self.reject_uphill
            and len(self.energies) > 1
            and energy - self.energies[-2] > self.uphill_tolerance
        ):
            actual_change = float(energy - self.energies[-2])
            self._trial_max_step = max(
                self._trial_max_step / 2.0,
                self.rejection_step_floor,
            )
            self.reject_current_trial()
            forces = self.forces[-1]
            energy = self.energies[-1]
            trial_rejected = True
            self.rejected_uphill_steps += 1
            at_floor = self._trial_max_step <= self.rejection_step_floor * (1.0 + 1e-12)
            self.rejections_at_floor = self.rejections_at_floor + 1 if at_floor else 0
            self.table.print(
                "Rejected uphill LBFGS trial and restored the previous geometry "
                f"(ΔE={actual_change:+.6f} au, max_step={self._trial_max_step:.3e})."
            )
            if self.rejections_at_floor >= self.max_rejections_at_floor:
                self.request_stop(
                    "repeated uphill LBFGS trials at the emergency step floor"
                )
                return np.zeros_like(forces)
        elif len(self.energies) > 1:
            self.rejections_at_floor = 0
            self._trial_max_step = min(
                self.max_step,
                self._trial_max_step * 2.0,
            )

        if (
            not trial_rejected
            and len(self.forces) > 1
            and self.forces[-2].size == forces.size
        ):
            y = self.forces[-2] - forces
            s = self.steps[-1]
            if self.double_damp:
                s, y = double_damp(
                    s, y, s_list=self.coord_diffs, y_list=self.grad_diffs
                )
            self.grad_diffs.append(y)
            self.coord_diffs.append(s)

            # Drop superfluous oldest vectors
            self.coord_diffs = self.coord_diffs[-self.keep_last :]
            self.grad_diffs = self.grad_diffs[-self.keep_last :]

        ###############
        # Line search #
        ###############

        ip_gradient, ip_step = None, None
        if self.line_search and len(self.energies) > 1:
            ip_energy, ip_gradient, ip_step = poly_line_search(
                energy, self.energies[-2], -forces, -self.forces[-2], self.steps[-1]
            )
        # Use the interpolated gradient for the step if interpolation succeeded
        if (ip_gradient is not None) and (ip_step is not None):
            forces = -ip_gradient
            self.log("Interpolation succeeded")
        # Keep the original gradient when the interpolation failed, but use
        # zero (interpolated) step.
        else:
            ip_step = np.zeros_like(forces)

        step = self.get_lbfgs_step(forces)

        # Skip adapation mu_reg in first cycle, because in this implementation
        # the first step will be steepest descent, which is independent of
        # self.mu_reg, so this loop would never break.
        adapt_mu_cycles = 0
        while self.mu_reg and (self.cur_cycle > 0):
            self.log(
                f"Adapt μ_reg={self.mu_reg:.6f}, norm(step)={np.linalg.norm(step):.6f}"
            )
            if adapt_mu_cycles == self.max_mu_reg_adaptions:
                raise Exception("Adapation of mu_reg failed! Breaking!")
            try:
                trial_energy = self.geometry.get_energy_at(self.geometry.coords + step)
                self.mu_reg, recompute_step = self.update_mu_reg(
                    self.mu_reg, energy, trial_energy, -forces, step
                )
            except NeedNewInternalsException:
                self.log("Internal coordinate breakdown in linesearch!")
                # Nothing further is done here, as the coordinates will probably also
                # breakdown when taking the step, which is then handled.
                recompute_step = False

            # Leave loop if step was accepted
            if not recompute_step:
                self.log(f"Next μ_reg={self.mu_reg:.6f}")
                break
            # Otherwise, recompute using updated μ_reg
            step = self.get_lbfgs_step(forces)
            adapt_mu_cycles += 1

        if self.mu_reg:
            self.tot_adapt_mu_cycles += adapt_mu_cycles + 1

        # Form full step. If we did not interpolate or interpolation failed,
        # ip_step will be zero.
        step = step + ip_step

        # Only try to scale down first step in regularized L-BFGS
        if (self.mu_reg and self.cur_cycle == 0) or self.control_step:
            step = scale_by_max_step(step, min(self.max_step, self._trial_max_step))

        return step

    def postprocess_opt(self):
        if self.rejected_uphill_steps:
            print(
                "LBFGS safeguards: "
                f"rejected {self.rejected_uphill_steps} uphill trial(s)."
            )
        if self.mu_reg:
            msg = f"\nNumber of μ updates: {self.tot_adapt_mu_cycles}"
            self.log(msg)
