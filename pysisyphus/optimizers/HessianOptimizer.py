from math import sqrt
from pathlib import Path
from typing import Literal, Optional

import numpy as np
from scipy.optimize import root_scalar

from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import rms
# from pysisyphus.io.hessian import save_hessian
from pysisyphus.optimizers.guess_hessians import (
    get_guess_hessian,
    xtb_hessian,
    HessInit,
)
from pysisyphus.optimizers.hessian_updates import (
    bfgs_update,
    flowchart_update,
    damped_bfgs_update,
    bofill_update,
    bofill_cpu_offload_enabled,
    ts_bfgs_update,
    ts_bfgs_update_org,
    ts_bfgs_update_revised,
)
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.exceptions import OptimizationError

from pysisyphus.helpers import array2string
import torch

# M4 (refinement_plan §M1+M4): one-place torch-vs-numpy dispatch. Only used at
# sites where the torch and numpy APIs share both the name and the semantics
# (`xp.linalg.eigh`, `xp.isfinite`, `xp.nan_to_num`, ...). Sites that use
# torch-specific ops (`index_select`, `matmul`, `.clone()`) keep their inline
# isinstance branches because the corresponding numpy code is structurally
# different, not just a different module.
from pysisyphus._array import get_xp, active_square

def dummy_hessian_update(H, dx, dg):
    return np.zeros_like(H), "no"


HESS_UPDATE_FUNCS = {
    "none": dummy_hessian_update,
    None: dummy_hessian_update,
    False: dummy_hessian_update,
    "bfgs": bfgs_update,
    "damped_bfgs": damped_bfgs_update,
    "flowchart": flowchart_update,
    "bofill": bofill_update,
    "ts_bfgs": ts_bfgs_update,
    "ts_bfgs_org": ts_bfgs_update_org,
    "ts_bfgs_rev": ts_bfgs_update_revised,
}
HessUpdate = Literal[
    "none",
    None,
    False,
    "bfgs",
    "damped_bfgs",
    "flowchart",
    "bofill",
    "ts_bfgs",
    "ts_bfgs_org",
    "ts_bfgs_rev",
]


class HessianOptimizer(Optimizer):
    rfo_dict = {
        "min": (0, "min"),
        "max": (-1, "max"),
    }

    def __init__(
        self,
        geometry: Geometry,
        trust_radius: float = 0.5,
        trust_update: bool = True,
        trust_min: float = 0.1,
        trust_max: float = 1,
        max_energy_incr: Optional[float] = None,
        hessian_update: HessUpdate = "bfgs",
        hessian_init: HessInit = "fischer",
        hessian_recalc: Optional[int] = None,
        hessian_recalc_adapt: Optional[float] = None,
        hessian_xtb: bool = False,
        hessian_recalc_reset: bool = False,
        small_eigval_thresh: float = 1e-8,
        line_search: bool = False,
        alpha0: float = 1.0,
        max_micro_cycles: int = 25,
        rfo_overlaps: bool = False,
        trust_band: bool = False,
        trust_band_sigma_inc: float = 1.15,
        trust_band_sigma_dec: float = 0.65,
        trust_band_rho_inc: float = 1.035,
        trust_band_rho_dec: float = 5.0,
        hessian_update_window: int = 1,
        # Opt-in: per-coord-type weighted L_inf trust check in lieu of
        # `np.linalg.norm(step)` (L2) when geometry has internals.
        weighted_trust: bool = False,
        weighted_trust_wb: float = 1.0,
        weighted_trust_wa: float = 0.5,
        weighted_trust_wd: float = 0.25,
        weighted_trust_wo: float = 0.5,
        weighted_trust_wx: float = 1.0,
        reject_uphill: bool = False,
        uphill_tolerance: float = 1e-8,
        rejection_trust_floor: float = 1e-7,
        max_rejections_at_floor: int = 3,
        **kwargs,
    ) -> None:
        """Baseclass for optimizers utilizing Hessian information.

        Parameters
        ----------
        geometry
            Geometry to be optimized.
        trust_radius
            Initial trust radius in whatever unit the optimization is carried out.
        trust_update
            Whether to update the trust radius throughout the optimization.
        trust_min
            Minimum trust radius.
        trust_max
            Maximum trust radius.
        max_energy_incr
            Maximum allowed energy increased after a faulty step. Optimization is
            aborted when the threshold is exceeded.
        hessian_update
            Type of Hessian update. Defaults to BFGS for minimizations and Bofill
            for saddle point searches.
        hessian_init
            Type of initial model Hessian.
        hessian_recalc
            Recalculate exact Hessian every n-th cycle instead of updating it.
        hessian_recalc_adapt
            Use a more flexible scheme to determine Hessian recalculation. Undocumented.
        hessian_xtb
            Recalculate the Hessian at the GFN2-XTB level of theory.
        hessian_recalc_reset
            Whether to skip Hessian recalculation after reset. Undocumented.
        small_eigval_thresh
            Threshold for small eigenvalues. Eigenvectors belonging to eigenvalues
            below this threshold are discardewd.
        line_search
            Whether to carry out a line search. Not implemented by a subclassing
            optimizers.
        alpha0
            Initial alpha for restricted-step (RS) procedure.
        max_micro_cycles
            Maximum number of RS iterations.
        rfo_overlaps
            Enable mode-following in RS procedure.
        reject_uphill
            Reject minimization trials whose actual energy rises.
        uphill_tolerance
            Positive energy change tolerated before a trial is rejected.
        rejection_trust_floor
            Emergency trust-radius floor used only while retrying a rejected
            uphill trial; it may be smaller than ``trust_min``.
        max_rejections_at_floor
            Repeated rejected trials allowed at the emergency floor before the
            lower-energy geometry is retained and checked once for convergence;
            the run stops as non-converged only if that check fails.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments passed to the Optimizer baseclass.
        """
        super().__init__(geometry, **kwargs)

        assert not issubclass(
            type(geometry), ChainOfStates
        ), "HessianOptimizer can't be used for and ChainOfStates objects!"

        self.trust_update = bool(trust_update)
        assert trust_min <= trust_max, "trust_min must be <= trust_max!"
        self.trust_min = float(trust_min)
        self.trust_max = float(trust_max)
        self.max_energy_incr = max_energy_incr
        self.reject_uphill = bool(reject_uphill)
        self.uphill_tolerance = float(uphill_tolerance)
        if self.uphill_tolerance < 0.0:
            raise ValueError("uphill_tolerance must be non-negative")
        self.rejection_trust_floor = float(rejection_trust_floor)
        if self.rejection_trust_floor <= 0.0:
            raise ValueError("rejection_trust_floor must be positive")
        self.max_rejections_at_floor = int(max_rejections_at_floor)
        if self.max_rejections_at_floor < 1:
            raise ValueError("max_rejections_at_floor must be at least 1")
        self.rejected_uphill_steps = 0
        self.rejections_at_floor = 0
        self.uphill_rejection_stalled = False
        self.skipped_bfgs_updates = 0
        # Constrain initial trust radius if trust_max > trust_radius
        self.trust_radius = min(trust_radius, trust_max)
        self.log(f"Initial trust radius: {self.trust_radius:.6f}")
        # Sella backport (opt-in): rho-band trust update parameters.
        # If trust_band=True, set_new_trust_radius uses gentler sigma_inc/sigma_dec
        # multipliers (1.15 / 0.65) inside an rho-band gate (1/rho_inc < rho < rho_inc
        # = "no update"; rho outside [1/rho_dec, rho_dec] = "shrink"). Default OFF
        # preserves Nocedal 4.1 binary x2/÷4 behaviour for reproducibility.
        self.trust_band = bool(trust_band)
        self.trust_band_sigma_inc = float(trust_band_sigma_inc)
        self.trust_band_sigma_dec = float(trust_band_sigma_dec)
        self.trust_band_rho_inc = float(trust_band_rho_inc)
        self.trust_band_rho_dec = float(trust_band_rho_dec)
        # window>=2 enables multistep_ts_bfgs_update from a sliding buffer of
        # recent (dx, dg) pairs; window=1 keeps single-step Bofill / BFGS / ts_bfgs.
        self.hessian_update_window = max(int(hessian_update_window), 1)
        # When weighted_trust is True AND geometry has internals, update_trust_radius
        # uses weighted_max_internal_step(step, typed_prims, **weights) instead of
        # plain L2 norm so that dihedrals don't burn the trust budget on small
        # cartesian moves (they get wd=0.25 by default vs. wb=1.0 for bonds).
        self.weighted_trust = bool(weighted_trust)
        self.weighted_trust_weights = dict(
            wb=float(weighted_trust_wb),
            wa=float(weighted_trust_wa),
            wd=float(weighted_trust_wd),
            wo=float(weighted_trust_wo),
            wx=float(weighted_trust_wx),
        )
        # Buffers populated lazily in update_hessian() when window>=2.
        self._sy_buffer_S: list = []
        self._sy_buffer_Y: list = []
        self.hessian_update = hessian_update
        self.hessian_update_func = HESS_UPDATE_FUNCS[hessian_update]
        self.hessian_init = hessian_init
        self.hessian_recalc = hessian_recalc
        self.hessian_recalc_adapt = hessian_recalc_adapt
        self.hessian_xtb = hessian_xtb
        self.hessian_recalc_reset = hessian_recalc_reset
        self.small_eigval_thresh = float(small_eigval_thresh)
        self.line_search = bool(line_search)
        # Restricted-step related
        self.alpha0 = alpha0
        self.max_micro_cycles = int(max_micro_cycles)
        assert max_micro_cycles >= 0
        self.rfo_overlaps = rfo_overlaps

        assert self.small_eigval_thresh > 0.0, "small_eigval_thresh must be > 0.!"
        if not self.restarted:
            self.hessian_recalc_in = None
            self.adapt_norm = None
            self.predicted_energy_changes = list()
        hessian_init_exists = Path(self.hessian_init).exists()
        if (
            # Allow actually calculated Hessians for all coordinate systems
            not hessian_init_exists
            and self.hessian_init not in ("calc", "xtb", "xtb1", "xtbff")
            # But disable model Hessian for Cartesian optimizations
            and self.geometry.coord_type in ("cart", "cartesian", "mwcartesian")
        ):
            self.hessian_init = "unit"
            self.log(
                f"Chosen initial (model) Hessian is incompatible with current "
                f"coord_type: {self.geometry.coord_type}!"
            )

        self._prev_eigvec_min = None
        self._prev_eigvec_max = None
        self._using_active_dofs = False
        self._active_dof_indices = None
        self.cur_H = None

    @property
    def using_active_dofs(self):
        return self._using_active_dofs

    @property
    def active_dof_indices(self):
        return self._active_dof_indices

    def _set_active_dofs(self, use_active):
        self._using_active_dofs = use_active
        if not use_active:
            self._active_dof_indices = None
            return
        if getattr(self.geometry, "within_partial_hessian", None) is not None:
            self._active_dof_indices = self.geometry.hess_active_dof_indices
            return
        # Fallback: infer active DOFs from calculator's Hessian-active atoms
        calc = getattr(self.geometry, "calculator", None)
        core = getattr(calc, "core", calc)
        hess_atoms = getattr(core, "hess_active_atoms", None)
        if hess_atoms is not None and len(hess_atoms) > 0:
            act = []
            for a in hess_atoms:
                base = 3 * int(a)
                act.extend([base, base + 1, base + 2])
            self._active_dof_indices = np.asarray(act, dtype=int)
            return
        self._active_dof_indices = self.geometry.active_dof_indices

    def active_from_full(self, vec):
        if not self.using_active_dofs:
            return vec
        inds = self._active_dof_indices
        if inds is None:
            return vec
        # Sanitize indices (drop negatives / out-of-bounds)
        try:
            inds = np.asarray(inds, dtype=int)
            if len(inds) > 0:
                if np.any(inds < 0):
                    inds = inds[inds >= 0]
                if len(inds) > 0:
                    max_valid = vec.shape[0] - 1
                    if np.any(inds > max_valid):
                        inds = inds[inds <= max_valid]
                if len(inds) > 0:
                    if np.min(inds) < 0 or np.max(inds) >= vec.shape[0]:
                        return vec
        except (ValueError, IndexError, TypeError):
            pass
        # Avoid double-slicing if vector is already in active space
        try:
            if vec.shape[0] == len(inds):
                return vec
        except (ValueError, IndexError, TypeError):
            pass
        try:
            if len(inds) > 0 and vec.shape[0] <= int(np.max(inds)):
                # Indices exceed vector length → already active (compact) vector.
                return vec
        except (ValueError, IndexError, TypeError):
            pass
        if isinstance(vec, torch.Tensor):
            # GPU-native: index_select avoids the CPU round-trip the legacy
            # CUDA branch used. Partial-DOF-correct: returns vec[inds] semantics.
            # Upstream sanitization at L226-251 already clamps inds, so the
            # deleted silent fallback (`except: return vec`) is unreachable in
            # practice and previously masked shape-mismatch bugs downstream.
            idx = torch.as_tensor(inds, dtype=torch.long, device=vec.device)
            return vec.index_select(0, idx)
        return vec[inds]

    def full_from_active(self, vec):
        if not self.using_active_dofs:
            return vec
        inds = self._active_dof_indices
        if inds is None:
            return vec
        if isinstance(vec, torch.Tensor):
            idx = torch.as_tensor(inds, dtype=torch.long, device=vec.device)
            full = torch.zeros(self.geometry.cart_coords.size, dtype=vec.dtype, device=vec.device)
            full.index_copy_(0, idx, vec)
            return full
        full = np.zeros(self.geometry.cart_coords.size, dtype=vec.dtype if hasattr(vec, "dtype") else float)
        full[inds] = vec
        return full

    def active_hessian(self, hessian):
        if not self.using_active_dofs:
            return hessian

        if getattr(self.geometry, "within_partial_hessian", None) is not None:
            act_n_dof = int(self.geometry.within_partial_hessian.get("active_n_dof", 0))
            if hessian.shape == (act_n_dof, act_n_dof):
                return hessian

        inds = self.active_dof_indices
        try:
            inds_arr = np.asarray(inds, dtype=int)
            if hessian.shape[0] == len(inds_arr) and len(inds_arr) > 0:
                if np.max(inds_arr) >= hessian.shape[0]:
                    # Likely already in active order; avoid double-slicing.
                    return hessian
        except (ValueError, IndexError, TypeError):
            pass
        try:
            inds = np.asarray(inds, dtype=int)
            if len(inds) > 0:
                max_valid = hessian.shape[0] - 1
                inds = inds[(inds >= 0) & (inds <= max_valid)]
        except (ValueError, IndexError, TypeError):
            pass
        if isinstance(hessian, torch.Tensor):
            # GPU-native: index_select avoids the CPU round-trip the legacy
            # CUDA branch used. Partial-Hessian-correct: returns (len(inds),)^2
            # identical to H[np.ix_(inds, inds)]. Upstream sanitization at
            # L289-302 already clamps inds to [0, n-1], so the deleted silent
            # fallback (`except: return hessian`) is unreachable in practice
            # and previously masked shape-mismatch bugs downstream.
            # C14 (H10): bounded row-chunk extraction avoids the full
            # (len(inds), n) row temporary of the chained index_select.
            idx = torch.as_tensor(inds, device=hessian.device, dtype=torch.int64)
            return active_square(hessian, idx)
        return hessian[np.ix_(inds, inds)]

    def active_list(self, seq):
        if not self.using_active_dofs:
            return seq
        return [self.active_from_full(item) for item in seq]

    @property
    def prev_eigvec_min(self):
        return self._prev_eigvec_min

    @prev_eigvec_min.setter
    def prev_eigvec_min(self, prev_eigvec_min):
        if self.rfo_overlaps:
            self._prev_eigvec_min = prev_eigvec_min

    @property
    def prev_eigvec_max(self):
        return self._prev_eigvec_max

    @prev_eigvec_max.setter
    def prev_eigvec_max(self, prev_eigvec_max):
        if self.rfo_overlaps:
            self._prev_eigvec_max = prev_eigvec_max

    def reset(self):
        # Don't recalculate the hessian if we have to reset the optimizer
        hessian_init = self.hessian_init
        if (
            (not self.hessian_recalc_reset)
            and hessian_init == "calc"
            and self.geometry.coord_type != "cart"
        ):
            hessian_init = "fischer"
        self.prepare_opt(hessian_init)

    # def save_hessian(self):
    #     # Don't try to save Hessians of analytical potentials
    #     if self.geometry.is_analytical_2d:
    #         return

    #     h5_fn = self.get_path_for_fn(f"hess_calc_cyc_{self.cur_cycle}.h5")
    #     # Save the cartesian hessian, as it is independent of the
    #     # actual coordinate system that is used.
    #     save_hessian(
    #         h5_fn,
    #         self.geometry,
    #         self.geometry.cart_hessian,
    #         self.geometry.energy,
    #         self.geometry.calculator.mult,
    #     )
    #     self.log(f"Wrote calculated cartesian Hessian to '{h5_fn}'")

    def prepare_opt(self, hessian_init=None):
        if hessian_init is None:
            hessian_init = self.hessian_init

        self.H, hess_str = get_guess_hessian(self.geometry, hessian_init)
        if self.hessian_init != "calc" and self.geometry.is_analytical_2d:
            assert self.H.shape == (3, 3)
            self.H[2, 2] = 0.0

        msg = f"Using {hess_str} Hessian"
        if hess_str == "saved":
            msg += f" from '{hessian_init}'"
        self.log(msg)

        # # Dump to disk if hessian was calculated
        # if self.hessian_init == "calc":
        #     self.save_hessian()

        if (
            hasattr(self.geometry, "coord_type")
            and self.geometry.coord_type == "dlc"
            # Calculated Hessian is already in DLC
            and hessian_init != "calc"
        ):
            U = self.geometry.internal.U
            self.H = U.T.dot(self.H).dot(U)

        if self.hessian_recalc_adapt:
            self.adapt_norm = np.linalg.norm(self.geometry.forces)

        if self.hessian_recalc:
            # Already substract one, as we don't do a hessian update in
            # the first cycle.
            self.hessian_recalc_in = self.hessian_recalc - 1

    # Subclass restart keys required for a bit-exact resume.  ``trust_radius``
    # is adapted every cycle in :meth:`update_trust_radius`; omitting it from
    # the restart payload silently reset a resumed run to
    # ``min(trust_radius, trust_max)`` and diverged the trajectory.  It is a
    # newly added key, so it is *restored tolerantly* (see
    # :meth:`_set_opt_restart_info`) and is deliberately absent from the
    # transactional-load required set below to keep pre-``trust_radius``
    # checkpoints loadable.  The same holds for the newer uphill-rejection
    # counters (``rejections_at_floor`` / ``rejected_uphill_steps`` /
    # ``uphill_rejection_stalled``), the multi-step Hessian buffer
    # (``_sy_buffer_S`` / ``_sy_buffer_Y``, empty unless
    # ``hessian_update_window >= 2``) and ``_prev_eigvec_min`` (``None`` unless
    # ``rfo_overlaps``): all are serialized unconditionally but restored
    # presence-guarded, so a checkpoint written before they existed still
    # loads.  Keeping them out of the required set below preserves that
    # backward tolerance; the presence-guard means an absent key never raises.
    required_opt_restart_keys = (
        "adapt_norm",
        "H",
        "hessian_recalc_in",
        "predicted_energy_changes",
    )

    def _get_opt_restart_info(self):
        opt_restart_info = {
            "adapt_norm": self.adapt_norm,
            "H": self.H.tolist(),
            "hessian_recalc_in": self.hessian_recalc_in,
            "predicted_energy_changes": self.predicted_energy_changes,
            # Per-cycle adaptive trust radius (see update_trust_radius).
            "trust_radius": self.trust_radius,
            # Uphill-rejection adaptive state (reject_current_uphill_step):
            # ``rejections_at_floor`` gates request_stop, so resetting it on a
            # resume changes the retry budget and the cycle a stalled run
            # terminates on.  The other two make the rejection state complete.
            "rejections_at_floor": self.rejections_at_floor,
            "rejected_uphill_steps": self.rejected_uphill_steps,
            "uphill_rejection_stalled": self.uphill_rejection_stalled,
            # Sliding (dx, dg) buffer for the multi-step TS-BFGS update used when
            # ``hessian_update_window >= 2``; the next Hessian update stacks
            # these columns, so an empty buffer on resume diverges the update.
            "_sy_buffer_S": [np.asarray(s).tolist() for s in self._sy_buffer_S],
            "_sy_buffer_Y": [np.asarray(y).tolist() for y in self._sy_buffer_Y],
            # Previous minimum-mode eigenvector for ``rfo_overlaps`` root
            # following; ``None`` unless overlap-based mode following is active.
            "_prev_eigvec_min": (
                None
                if self._prev_eigvec_min is None
                else np.asarray(self._prev_eigvec_min).tolist()
            ),
        }
        return opt_restart_info

    def _set_opt_restart_info(self, opt_restart_info):
        self.adapt_norm = opt_restart_info["adapt_norm"]
        self.H = np.array(opt_restart_info["H"])
        self.hessian_recalc_in = opt_restart_info["hessian_recalc_in"]
        self.predicted_energy_changes = opt_restart_info["predicted_energy_changes"]
        # Backward-tolerant: a checkpoint written before trust_radius was
        # serialized keeps the __init__-constrained value already in place.
        if "trust_radius" in opt_restart_info:
            self.trust_radius = float(opt_restart_info["trust_radius"])
        # Uphill-rejection adaptive state.  All presence-guarded so a checkpoint
        # written before these keys existed keeps the __init__ defaults.
        if "rejections_at_floor" in opt_restart_info:
            self.rejections_at_floor = int(opt_restart_info["rejections_at_floor"])
        if "rejected_uphill_steps" in opt_restart_info:
            self.rejected_uphill_steps = int(opt_restart_info["rejected_uphill_steps"])
        if "uphill_rejection_stalled" in opt_restart_info:
            self.uphill_rejection_stalled = bool(
                opt_restart_info["uphill_rejection_stalled"]
            )
        # Multi-step Hessian-update buffer; only non-empty for
        # ``hessian_update_window >= 2``.  Default [] when the key is absent.
        if "_sy_buffer_S" in opt_restart_info:
            self._sy_buffer_S = [np.array(s) for s in opt_restart_info["_sy_buffer_S"]]
        if "_sy_buffer_Y" in opt_restart_info:
            self._sy_buffer_Y = [np.array(y) for y in opt_restart_info["_sy_buffer_Y"]]
        # Previous minimum-mode eigenvector for ``rfo_overlaps``; default None.
        if "_prev_eigvec_min" in opt_restart_info:
            stored = opt_restart_info["_prev_eigvec_min"]
            self._prev_eigvec_min = None if stored is None else np.array(stored)

    def update_trust_radius(self):
        # The predicted change should be calculated at the end of optimize
        # of the previous cycle.
        assert (
            len(self.predicted_energy_changes) == len(self.forces) - 1
        ), "Did you forget to append to self.predicted_energy_changes?"
        self.log("Trust radius update")
        self.log(f"\tCurrent trust radius: {self.trust_radius:.6f}")
        predicted_change = self.predicted_energy_changes[-1]
        actual_change = self.energies[-1] - self.energies[-2]
        # Only report an unexpected increase if we actually predicted a
        # decrease.
        unexpected_increase = (actual_change > 0) and (predicted_change < 0)
        old_trust = self.trust_radius
        if unexpected_increase:
            self.log(f"Energy increased by {actual_change:.6f} au!")
            if self.max_energy_incr and (actual_change > self.max_energy_incr):
                raise OptimizationError("Actual energy change too high!")
        coeff = actual_change / predicted_change
        self.log(f"\tPredicted change: {predicted_change:.4e} au")
        self.log(f"\tActual change: {actual_change:.4e} au")
        self.log(f"\tCoefficient: {coeff:.2%}")
        step = self.steps[-1]
        # Weighted L_inf norm replaces L2 when weighted_trust is enabled AND
        # geometry has internals with matching typed_prims length.
        internal = getattr(self.geometry, "internal", None)
        typed_prims = getattr(internal, "typed_prims", None) if internal is not None else None
        if (
            self.weighted_trust
            and typed_prims is not None
            and len(typed_prims) == len(np.asarray(step).ravel())
        ):
            from pysisyphus.optimizers.restrict_step import weighted_max_internal_step
            last_step_norm = weighted_max_internal_step(
                step, typed_prims, **self.weighted_trust_weights
            )
        else:
            last_step_norm = np.linalg.norm(step)
        rejected_uphill = (
            unexpected_increase
            and self.reject_uphill
            and actual_change > self.uphill_tolerance
        )
        min_radius = self.rejection_trust_floor if rejected_uphill else None
        self.set_new_trust_radius(coeff, last_step_norm, min_radius=min_radius)
        if unexpected_increase:
            self.table.print(
                f"Unexpected energy increase ({actual_change:.6f} au)! "
                f"Trust radius: old={old_trust:.4}, new={self.trust_radius:.4}"
            )
        return unexpected_increase

    def _rollback_trial_state(self):
        if not self.predicted_energy_changes:
            raise RuntimeError("Missing model prediction for Hessian trial rejection.")
        self.predicted_energy_changes.pop()

    def reject_current_uphill_step(self):
        """Reject an uphill minimization trial through the shared transaction."""
        actual_change = float(self.energies[-1] - self.energies[-2])
        self.reject_current_trial()
        self.rejected_uphill_steps += 1
        self.table.print(
            "Rejected uphill RFO trial and restored the previous geometry "
            f"(ΔE={actual_change:+.6f} au)."
        )
        at_floor = self.trust_radius <= self.rejection_trust_floor * (1.0 + 1e-12)
        self.rejections_at_floor = self.rejections_at_floor + 1 if at_floor else 0
        if self.rejections_at_floor >= self.max_rejections_at_floor:
            self.uphill_rejection_stalled = True
            self.table.print(
                "Repeated uphill RFO trials at the emergency trust floor; "
                "retaining the lower-energy geometry and checking it for "
                "convergence."
            )
            self.request_stop(
                "repeated uphill RFO trials at the emergency trust floor"
            )

    def set_new_trust_radius(self, coeff, last_step_norm, min_radius=None):
        if min_radius is None:
            min_radius = self.trust_min
        # A rejected RFO trial may deliberately move the radius below the
        # ordinary user-facing trust_min.  Never make a small radius larger in
        # the branch whose purpose is to shrink it.
        min_radius = min(float(min_radius), self.trust_radius)
        if self.trust_band:
            # Sella backport (opt-in): rho-band trust update (sella/optimize/optimize.py:280-289).
            # `coeff` here is ΔE_actual / ΔE_predicted -- same definition as Sella's `rho`.
            # rho outside [1/rho_dec, rho_dec]      -> shrink with sigma_dec (default 0.65)
            # rho inside [1/rho_inc, rho_inc] band  -> grow with sigma_inc  (default 1.15)
            # otherwise                              -> keep current trust
            rho = float(coeff)
            if (rho < 1.0 / self.trust_band_rho_dec) or (rho > self.trust_band_rho_dec):
                new_tr = max(last_step_norm * self.trust_band_sigma_dec, min_radius)
                self.trust_radius = min(new_tr, self.trust_max)
                self.log(f"\tTrust band shrink (rho={rho:.3f}): {self.trust_radius:.6f}")
            elif (1.0 / self.trust_band_rho_inc) < rho < self.trust_band_rho_inc:
                new_tr = max(last_step_norm * self.trust_band_sigma_inc, self.trust_radius)
                self.trust_radius = min(new_tr, self.trust_max)
                self.log(f"\tTrust band grow (rho={rho:.3f}): {self.trust_radius:.6f}")
            else:
                self.log(f"\tTrust band keep (rho={rho:.3f}): {self.trust_radius:.6f}")
            return

        # Nocedal, Numerical optimization Chapter 4, Algorithm 4.1

        # If actual and predicted energy change have different signs
        # coeff will be negative and lead to a decreased trust radius,
        # which is fine.
        if coeff < 0.25:
            self.trust_radius = max(self.trust_radius / 4, min_radius)
            self.log("\tDecreasing trust radius.")
        # Only increase trust radius if last step norm was at least 80% of it
        # See [5], Appendix, step size and direction control
        # elif coeff > 0.75 and (last_step_norm >= .8*self.trust_radius):
        #
        # Only increase trust radius if last step norm corresponded approximately
        # to the trust radius.
        elif coeff > 0.75 and abs(self.trust_radius - last_step_norm) <= 1e-3:
            self.trust_radius = min(self.trust_radius * 2, self.trust_max)
            self.log("\tIncreasing trust radius.")
        else:
            self.log(f"\tKeeping current trust radius at {self.trust_radius:.6f}")
            return
        self.log(f"\tUpdated trust radius: {self.trust_radius:.6f}")

    def update_hessian(self):
        # Compare current forces to reference forces to see if we shall recalc the
        # hessian.
        try:
            cur_norm = np.linalg.norm(self.forces[-1])
            ref_norm = self.adapt_norm / self.hessian_recalc_adapt
            recalc_adapt = cur_norm <= ref_norm
            self.log(
                "Check for adaptive Hessian recalculation: "
                f"{cur_norm:.6f} <= {ref_norm:.6f}, {recalc_adapt}"
            )
        except TypeError:
            recalc_adapt = False

        try:
            self.hessian_recalc_in = max(self.hessian_recalc_in - 1, 0)
            self.log(f"Recalculation of Hessian in {self.hessian_recalc_in} cycle(s).")
        except TypeError:
            self.hessian_recalc_in = None

        # Update reference norm if needed
        # TODO: Decide on whether to update the norm when the recalculation is
        # initiated by 'recalc'.
        if recalc_adapt:
            self.adapt_norm = cur_norm

        recalc = self.hessian_recalc_in == 0

        if recalc or recalc_adapt:
            # Free old Hessian from GPU before recalculating
            H_old = self.H
            self.H = None
            self.cur_H = None
            del H_old
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            # Use xtb hessian
            self.log("Requested Hessian recalculation.")
            if self.hessian_xtb:
                self.H = xtb_hessian(self.geometry)
                key = "xtb"
            # Calculated hessian at actual level of theory
            else:
                self.H = self.geometry.hessian
                key = "exact"
                # self.save_hessian()
            if self.using_active_dofs:
                # Keep the optimizer Hessian in active-DOF space to avoid
                # shape mismatches during quasi-Newton updates.
                self.H = self.active_hessian(self.H)
            if not (self.cur_cycle == 0):
                self.log(f"Recalculated {key} Hessian in cycle {self.cur_cycle}.")
            # Reset counter. It is also reset when the recalculation was initiated
            # by the adaptive formulation.
            self.hessian_recalc_in = self.hessian_recalc
        # Simple hessian update
        else:
            dx = self.steps[-1]
            dg = -(self.forces[-1] - self.forces[-2])
            H_work = self.H
            if self.using_active_dofs:
                H_work = self.active_hessian(self.H)
                dx = self.active_from_full(dx)
                dg = self.active_from_full(dg)
            curv_cond = dx.dot(dg)
            curv_value = (
                float(curv_cond.detach().cpu())
                if isinstance(curv_cond, torch.Tensor)
                else float(curv_cond)
            )
            if curv_value <= 0.0:
                self.log(
                    f"Curvature condition (s·y = {curv_value:.4f} <= 0) not satisfied!"
                )
                # The BFGS Hessian update divides by s·y and only preserves a
                # minimization model when this curvature condition is positive.
                # Applying it anyway can inject an arbitrarily large negative
                # eigenvalue, which in turn sends RFO to the trust boundary.
                if self.hessian_update == "bfgs":
                    self.skipped_bfgs_updates += 1
                    self.log("Skipped unsafe BFGS Hessian update.")
                    return
            if self.hessian_update_window >= 2:
                # Multi-step TS-BFGS update from a sliding (dx, dg) buffer.
                # Multistep helper is numpy-only, so detach torch tensors to CPU
                # for the update and move the result back to H's device. Same
                # pattern as bofill_update's torch handling.
                from pysisyphus.optimizers.hessian_updates import multistep_ts_bfgs_update
                if isinstance(H_work, torch.Tensor):
                    H_np = H_work.detach().cpu().numpy()
                    dx_np = (
                        dx.detach().cpu().numpy()
                        if isinstance(dx, torch.Tensor) else np.asarray(dx)
                    )
                    dg_np = (
                        dg.detach().cpu().numpy()
                        if isinstance(dg, torch.Tensor) else np.asarray(dg)
                    )
                else:
                    H_np = np.asarray(H_work)
                    dx_np = np.asarray(dx)
                    dg_np = np.asarray(dg)
                self._sy_buffer_S.append(dx_np.copy())
                self._sy_buffer_Y.append(dg_np.copy())
                self._sy_buffer_S = self._sy_buffer_S[-self.hessian_update_window:]
                self._sy_buffer_Y = self._sy_buffer_Y[-self.hessian_update_window:]
                S = np.column_stack(self._sy_buffer_S)
                Y = np.column_stack(self._sy_buffer_Y)
                dH_np, key = multistep_ts_bfgs_update(H_np, S, Y)
                if isinstance(H_work, torch.Tensor):
                    dH = torch.as_tensor(dH_np, dtype=H_work.dtype, device=H_work.device)
                else:
                    dH = dH_np
                self.log(f"Did {key} Hessian update (window={len(self._sy_buffer_S)}).")
                self.H = H_work + dH
            elif (
                isinstance(H_work, torch.Tensor)
                and self.hessian_update == "bofill"
                and not bofill_cpu_offload_enabled()
            ):
                # C14: copy-on-write rank-two Bofill. Clone the accepted Hessian
                # once and apply the low-rank update into the clone with addmm_.
                # This preserves the replacement/rollback invariant (the trial
                # snapshot keeps the old H reference, which is never mutated)
                # while dropping the transient peak from H_work + dH + new_H to
                # H_work + new_H (one fewer dense N x N matrix). Numerically
                # identical to H_work + (U @ C @ U.T) to fp roundoff.
                from pysisyphus.optimizers.hessian_updates import bofill_rank2_factors

                U, C = bofill_rank2_factors(H_work, dx, dg)
                new_H = H_work.clone()
                new_H.addmm_(U, C @ U.t())
                self.H = new_H
                self.log("Did Bofill Hessian update (rank-2, copy-on-write).")
            else:
                dH, key = self.hessian_update_func(H_work, dx, dg)
                self.log(f"Did {key} Hessian update.")
                self.H = H_work + dH

    def solve_rfo(self, rfo_mat, kind="min", prev_eigvec=None):
        # When using the restricted step variant of RFO the RFO matrix
        # may not be symmetric. Thats why we can't use eigh here.
        # M4: torch/numpy dispatch via _array.get_xp — both backends share
        # isfinite / nan_to_num / allclose / linalg.{eigh,eig} / argsort by
        # name, so xp.foo(...) is equivalent to the prior if-is_torch branch.
        xp = get_xp(rfo_mat)
        is_torch = xp is torch  # still needed for the typed-exception in eigh fallback

        if not xp.isfinite(rfo_mat).all():
            self.log("RFO matrix contains NaN/inf; sanitizing entries.")
            rfo_mat = xp.nan_to_num(rfo_mat, nan=0.0, posinf=1e8, neginf=-1e8)

        is_sym = xp.allclose(rfo_mat, rfo_mat.T)

        if is_sym:
            try:
                eigenvalues, eigenvectors = xp.linalg.eigh(rfo_mat)
            except (torch._C._LinAlgError, np.linalg.LinAlgError):
                self.log("eigh failed; falling back to eig.")
                eigenvalues, eigenvectors = xp.linalg.eig(rfo_mat)
                eigenvalues = eigenvalues.real
                eigenvectors = eigenvectors.real
        else:
            eigenvalues, eigenvectors = xp.linalg.eig(rfo_mat)
            eigenvalues = eigenvalues.real
            eigenvectors = eigenvectors.real

        self.log("\tdiagonalized augmented Hessian")

        sorted_inds = get_xp(eigenvectors).argsort(eigenvalues)

        # Depending on wether we want to minimize (maximize) along
        # the mode(s) in the rfo mat we have to select the smallest
        # (biggest) eigenvalue and corresponding eigenvector.
        first_or_last, verbose = self.rfo_dict[kind]
        # Given sorted eigenvalue-indices (sorted_inds) use the first
        # (smallest eigenvalue) or the last (largest eigenvalue) index.
        if prev_eigvec is None:
            ind = sorted_inds[first_or_last]
        else:
            if isinstance(prev_eigvec, torch.Tensor):
                ovlps = prev_eigvec.matmul(eigenvectors)
            else:
                ovlps = np.array([prev_eigvec.dot(ev) for ev in eigenvectors.T])
            naive_ind = sorted_inds[first_or_last]
            ind = np.abs(ovlps).argmax() if isinstance(ovlps, np.ndarray) else torch.argmax(torch.abs(ovlps)).item()
            self.log(
                f"Overlap: {ind} ({eigenvalues[ind]:.6f}), "
                f"Naive: {naive_ind} ({eigenvalues[naive_ind]:.6f})"
            )
        follow_eigvec = eigenvectors.T[ind]
        if isinstance(follow_eigvec, torch.Tensor):
            step_nu = follow_eigvec.clone()
        else:
            step_nu = follow_eigvec.copy()
        nu = step_nu[-1]
        self.log(f"\tnu_{verbose}={nu:.8e}")
        # Scale eigenvector so that its last element equals 1. The
        # final is step is the scaled eigenvector without the last element.
        step = step_nu[:-1] / nu
        eigval = eigenvalues[ind]
        self.log(f"\teigenvalue_{verbose}={eigval:.8e}")
        return step, eigval, nu, follow_eigvec

    def solve_rfo_secular(self, eigvals, gradient, alpha=1.0, kind="min",
                          prev_eigvec=None, max_iter=50, tol=1e-12):
        """Solve the RFO eigenvalue problem via the secular equation.

        The augmented Hessian has arrowhead structure, so its eigenvalue
        problem reduces to: f(mu) = sum g_i^2/(alpha*mu - lam_i) - mu = 0,
        solvable in O(N) instead of O(N^3).

        Returns (step, eigval, nu, eigvec) on success, None on failure.
        """
        is_torch = isinstance(eigvals, torch.Tensor)

        # Convert all inputs to plain numpy/float for root-finding
        _np = lambda x: x.detach().cpu().numpy().copy() if isinstance(x, torch.Tensor) else np.asarray(x, dtype=np.float64)
        g = _np(gradient)
        lam = _np(eigvals)
        alpha = alpha.detach().cpu().item() if isinstance(alpha, torch.Tensor) else float(alpha)

        n = len(lam)
        g2 = g ** 2
        nz = g2 > 1e-30
        if not nz.any():
            step = np.zeros(n)
            eigvec = np.zeros(n + 1); eigvec[-1] = 1.0
            if is_torch:
                step = torch.zeros(n, device=eigvals.device, dtype=eigvals.dtype)
                eigvec = torch.zeros(n + 1, device=eigvals.device, dtype=eigvals.dtype)
                eigvec[-1] = 1.0
            return step, 0.0, 1.0, eigvec

        g2_nz, lam_nz = g2[nz], lam[nz]

        # Guard against alpha ≈ 0 (can arise from trust-radius adaptation)
        if abs(alpha) < 1e-14:
            return None

        def f_df(mu):
            d = alpha * mu - lam_nz
            return float(np.sum(g2_nz / d) - mu), float(-alpha * np.sum(g2_nz / d**2) - 1.0)

        _, verbose = self.rfo_dict[kind]

        # Bracket the root
        if kind == "min":
            pole = float(lam.min() / alpha)
            mu = pole - max(float(np.sqrt(g2_nz.sum())) / alpha, 1.0)
            for _ in range(20):
                if f_df(mu)[0] > 0: break
                mu = pole - 2.0 * (pole - mu)
            else:
                return None
            lo, hi = mu, pole - 1e-15 * max(abs(pole), 1.0)
        elif kind == "max":
            pole = float(lam.max() / alpha)
            mu = pole + max(float(np.sqrt(g2_nz.sum())) / alpha, 1.0)
            for _ in range(20):
                if f_df(mu)[0] < 0: break
                mu = pole + 2.0 * (mu - pole)
            else:
                return None
            lo, hi = pole + 1e-15 * max(abs(pole), 1.0), mu
        else:
            return None

        # Newton-Raphson with bisection safeguard
        mu_cur = (lo + hi) / 2.0
        for _ in range(max_iter):
            fval, dfval = f_df(mu_cur)
            if abs(fval) < tol:
                break
            mu_new = mu_cur - fval / dfval if abs(dfval) > 1e-30 else (lo + hi) / 2.0
            if mu_new <= lo or mu_new >= hi:
                mu_new = (lo + hi) / 2.0
            f_new = f_df(mu_new)[0]
            if f_new > 0: lo = mu_new
            else: hi = mu_new
            mu_cur = mu_new
        else:
            self.log(f"Secular equation did not converge in {max_iter} iters.")
            return None

        self.log(f"\teigenvalue_{verbose}={mu_cur:.8e} (secular)")
        self.log(f"\tnu_{verbose}={1.0:.8e}")

        # Compute step: s_i = g_i / (alpha * mu - lam_i)
        denom = alpha * mu_cur - lam
        step_np = np.where(nz, g / denom, 0.0)

        # Eigenvector for mode tracking
        eigvec_np = np.append(step_np, 1.0)
        eigvec_np /= np.linalg.norm(eigvec_np)

        # Mode tracking check
        if prev_eigvec is not None:
            prev_np = _np(prev_eigvec)
            if abs(float(np.dot(prev_np, eigvec_np))) < 0.3:
                self.log("Secular eigvec overlap too low; falling back.")
                return None

        if is_torch:
            step_np = torch.tensor(step_np, device=eigvals.device, dtype=eigvals.dtype)
            eigvec_np = torch.tensor(eigvec_np, device=eigvals.device, dtype=eigvals.dtype)

        return step_np, mu_cur, 1.0, eigvec_np

    def filter_small_eigvals(self, eigvals, eigvecs, mask=False):
        # M4: xp.abs handles both backends (torch.abs / np.abs share name + semantics).
        small_inds = get_xp(eigvals).abs(eigvals) < self.small_eigval_thresh
        eigvals = eigvals[~small_inds]
        eigvecs = eigvecs[:, ~small_inds]
        small_num = sum(small_inds)
        self.log(
            f"Found {small_num} small eigenvalues in Hessian. Removed "
            "corresponding eigenvalues and eigenvectors."
        )
        # assert small_num <= 6, (
        #     "Expected at most 6 small eigenvalues in cartesian hessian "
        #     f"but found {small_num}!"
        # )
        if mask:
            return eigvals, eigvecs, small_inds
        else:
            return eigvals, eigvecs

    def log_negative_eigenvalues(self, eigvals, pre_str=""):
        neg_inds = eigvals < -self.small_eigval_thresh
        neg_eigval_str = array2string(eigvals[neg_inds], precision=6)
        self.log(f"{pre_str}Hessian has {neg_inds.sum()} negative eigenvalue(s).")
        self.log(f"\t{neg_eigval_str}")

    def housekeeping(self):
        """Calculate gradient and energy. Update trust radius and hessian
        if needed. Return energy, gradient and hessian for the current cycle."""
        gradient_full = self.geometry.gradient
        energy = self.geometry.energy
        self.energies.append(energy)
        self.log(f"    Energy: {energy: >12.6f} au")
        self.log(
            f"norm(grad): {np.linalg.norm(gradient_full): >12.6f} au / bohr (rad)"
        )
        self.log(
            f" rms(grad): {np.sqrt(np.mean(gradient_full**2)): >12.6f} au / bohr (rad)"
        )
        self.forces.append(-gradient_full)

        can_update = (
            # Allows gradient differences
            len(self.forces) > 1
            and (self.forces[-2].shape == gradient_full.shape)
            and len(self.coords) > 1
            # Coordinates may have been rebuilt. Take care of that.
            and (self.coords[-2].shape == self.coords[1].shape)
            and len(self.energies) > 1
        )
        if can_update:
            unexpected_increase = False
            if self.trust_update:
                unexpected_increase = self.update_trust_radius()
            actual_change = self.energies[-1] - self.energies[-2]
            reject_trial = (
                self.reject_uphill
                and actual_change > self.uphill_tolerance
            )
            if reject_trial:
                # A minimizer must not accept an uphill point even if a faulty
                # model happened to predict an increase too.  Such a case does
                # not enter update_trust_radius's "unexpected" branch, so force
                # the radius down here to ensure that the retry differs.
                if not unexpected_increase:
                    self.trust_radius = max(
                        self.trust_radius / 4.0,
                        min(self.rejection_trust_floor, self.trust_radius),
                    )
                self.reject_current_uphill_step()
                energy = self.energies[-1]
                gradient_full = -self.forces[-1]
                can_update = False
            else:
                self.rejections_at_floor = 0
                self.update_hessian()

        gradient, H, eigvals, eigvecs = self._hessian_system(gradient_full)

        resetted = not can_update
        return energy, gradient, H, eigvals, eigvecs, resetted

    def _hessian_system(self, gradient_full):
        """Return the active Hessian eigensystem for ``gradient_full``.

        Keeping this transformation separate lets a safeguarded optimizer
        restore or exactly recalculate its Hessian and rebuild a consistent
        eigensystem without appending duplicate energy/force history entries.
        """
        # Convert gradient to match H device/dtype AFTER update_hessian(),
        # so that hessian_recalc (which may replace self.H with a new tensor
        # on a different device) is accounted for.
        if isinstance(self.H, torch.Tensor):
            gradient_full = torch.from_numpy(gradient_full).to(
                self.H.device, self.H.dtype
            )

        H = self.H
        if self.geometry.internal:
            # Shift eigenvalues of orthogonal part to high values, so they
            # don't contribute to the actual step.
            H_proj = self.geometry.internal.project_hessian(self.H)
            # Symmetrize hessian, as the projection may break it?!
            H = (H_proj + H_proj.T) / 2

        if self.geometry.internal:
            # Geometry.hessian expands a partial Cartesian block to 3N before
            # transforming it into the selected internal-coordinate basis.
            # The transformed Hessian and gradient are therefore already a
            # complete, self-consistent internal system. Cartesian PHVA DOF
            # indices cannot be applied to that basis (its dimension/order is
            # unrelated), and expanding its step with those indices causes a
            # shape mismatch in TS optimizers.
            use_active = False
        elif getattr(self.geometry, "within_partial_hessian", None) is not None:
            use_active = True
        elif (
            H.shape[0] != self.geometry.cart_coords.size
            and self.geometry.coord_type in ("cart", "cartesian", "mwcartesian")
        ):
            # Partial Hessian without explicit metadata: still use active slicing.
            # Only applies to Cartesian coordinate types; for internal coordinates
            # (e.g. DLC), the Hessian is naturally smaller than cart_coords.size.
            use_active = True
        else:
            use_active = (
                len(self.geometry.freeze_atoms) > 0
                and self.geometry.coord_type in ("cart", "cartesian", "mwcartesian")
                and H.shape[0] == self.geometry.cart_coords.size
            )
        self._set_active_dofs(use_active)

        H = self.active_hessian(H)
        if gradient_full.shape[0] == H.shape[0]:
            gradient = gradient_full
        else:
            gradient = self.active_from_full(gradient_full)

        # M4: xp.linalg.eigh handles both backends.
        eigvals, eigvecs = get_xp(H).linalg.eigh(H)
        # Neglect small eigenvalues
        eigvals, eigvecs = self.filter_small_eigvals(eigvals, eigvecs)

        self.cur_H = H
        return gradient, H, eigvals, eigvecs

    def get_augmented_hessian(self, eigvals, gradient, alpha=1.0):
        if isinstance(gradient, torch.Tensor):
            dim_ = eigvals.size(0) + 1
            H_aug = torch.zeros((dim_, dim_), device=gradient.device, dtype=gradient.dtype)
            H_aug[: dim_ - 1, : dim_ - 1] = torch.diag(eigvals / alpha)
        else:
            dim_ = eigvals.size + 1
            H_aug = np.zeros((dim_, dim_))
            H_aug[: dim_ - 1, : dim_ - 1] = np.diag(eigvals / alpha)
        H_aug[-1, :-1] = gradient
        H_aug[:-1, -1] = gradient

        H_aug[:-1, -1] /= alpha

        return H_aug

    def get_alpha_step(self, cur_alpha, rfo_eigval, step_norm, eigvals, gradient):
        # Derivative of the squared step w.r.t. alpha
        numer = gradient**2
        denom = (eigvals - rfo_eigval * cur_alpha) ** 3
        # M4: xp.sum / xp.isfinite / xp.abs all share name + semantics across torch+np.
        xp = get_xp(gradient)
        quot = xp.sum(numer / denom)
        self.log(f"quot={quot:.6f}")
        dstep2_dalpha = 2 * rfo_eigval / (1 + step_norm**2 * cur_alpha) * quot
        # bool(...) coerces both numpy scalar and torch 0-d tensor to Python bool.
        dstep2_valid = bool(
            xp.isfinite(dstep2_dalpha) & (xp.abs(dstep2_dalpha) > 1e-12)
        )
        if not dstep2_valid:
            self.log(
                "alpha update skipped due to invalid derivative "
                f"(dstep2_dalpha={dstep2_dalpha})"
            )
            return 0.0
        self.log(f"analytic deriv.={dstep2_dalpha:.6f}")
        # Update alpha
        alpha_step = (
            2 * (self.trust_radius * step_norm - step_norm**2) / dstep2_dalpha
        )
        self.log(f"alpha_step={alpha_step:.4f}")
        min_alpha = 1e-8
        if (cur_alpha + alpha_step) <= min_alpha:
            self.log(
                "alpha update would make alpha non-positive; "
                f"clamping to min_alpha={min_alpha:.1e}"
            )
            alpha_step = min_alpha - cur_alpha
        return alpha_step

    def get_rs_step(self, eigvals, eigvecs, gradient, name="RS"):
        # Transform gradient to basis of eigenvectors
        if isinstance(eigvecs, torch.Tensor):
            if not isinstance(gradient, torch.Tensor):
                gradient = torch.as_tensor(
                    gradient, device=eigvecs.device, dtype=eigvecs.dtype
                )
            elif gradient.device != eigvecs.device:
                gradient = gradient.to(device=eigvecs.device)
            gradient_ = eigvecs.T @ gradient
        else:
            gradient_ = eigvecs.T.dot(gradient)

        alpha = self.alpha0
        for mu in range(self.max_micro_cycles):
            self.log(f"{name} micro cycle {mu:02d}, alpha={alpha:.6f}")
            # Try secular equation solver first (O(N) vs O(N^3))
            secular_result = self.solve_rfo_secular(
                eigvals, gradient_, alpha, kind="min",
                prev_eigvec=self.prev_eigvec_min,
            )
            if secular_result is not None:
                rfo_step_, eigval_min, nu, self.prev_eigvec_min = secular_result
            else:
                # Fallback to full eigendecomposition
                self.log("Secular solver failed; using full eigendecomposition.")
                H_aug = self.get_augmented_hessian(eigvals, gradient_, alpha)
                rfo_step_, eigval_min, nu, self.prev_eigvec_min = self.solve_rfo(
                    H_aug, "min", prev_eigvec=self.prev_eigvec_min
                )
            if isinstance(rfo_step_, torch.Tensor):
                rfo_norm_ = torch.linalg.norm(rfo_step_)
            else:
                rfo_norm_ = np.linalg.norm(rfo_step_)
            self.log(f"norm(rfo step)={rfo_norm_:.6f}")

            if rfo_norm_ <= 0:
                self.log(
                    "RFO step length is zero; falling back to trust-region Newton step."
                )
                step_ = self.get_newton_step_on_trust(
                    eigvals, eigvecs, gradient, transform=False
                )
                break

            if (rfo_norm_ < self.trust_radius) or abs(
                rfo_norm_ - self.trust_radius
            ) <= 1e-3:
                step_ = rfo_step_
                break

            alpha_step = self.get_alpha_step(
                alpha, eigval_min, rfo_norm_, eigvals, gradient_
            )
            alpha += alpha_step
            self.log("")
        # Otherwise, use trust region newton step
        else:
            self.log(
                "RS algorithm did not produce a desired step length in "
                f"{self.max_micro_cycles} micro cycles. Trying RFO with α=1.0."
            )
            secular_result = self.solve_rfo_secular(
                eigvals, gradient_, alpha=1.0, kind="min"
            )
            if secular_result is not None:
                rfo_step_, eigval_min, nu, _ = secular_result
            else:
                H_aug = self.get_augmented_hessian(eigvals, gradient_, alpha=1.0)
                rfo_step_, eigval_min, nu, _ = self.solve_rfo(H_aug, "min")
            if isinstance(rfo_step_, torch.Tensor):
                rfo_norm_ = torch.linalg.norm(rfo_step_)
            else:
                rfo_norm_ = np.linalg.norm(rfo_step_)

            # This should always be True if the above algorithm failed but we
            # keep this line nonetheless,  to make it more obvious.
            if rfo_norm_ > self.trust_radius:
                self.log(
                    f"Proposed RFO step with norm {rfo_norm_:.4f} is outside trust "
                    f"radius Δ={self.trust_radius:.4f}. "
                )
                step_ = self.get_newton_step_on_trust(
                    eigvals, eigvecs, gradient, transform=False
                )
                # Simple, downscaled RFO step
                # step_ = rfo_step_ / rfo_norm_ * self.trust_radius
            else:
                step_ = rfo_step_

        # Transform step back to original basis
        if isinstance(eigvecs, torch.Tensor):
            if isinstance(step_, torch.Tensor):
                pass
            else:
                step_ = torch.tensor(step_, device=eigvecs.device, dtype=eigvecs.dtype)
            step = eigvecs @ step_
            step = step.cpu().numpy()
        else:
            step = eigvecs.dot(step_)
        min_step_norm = getattr(self, "min_step_norm", 0.0)
        step_norm = np.linalg.norm(step)
        if step_norm <= min_step_norm:
            self.log(
                "RFO step length below minimum threshold; "
                "falling back to trust-region Newton step."
            )
            step = self.get_newton_step_on_trust(eigvals, eigvecs, gradient)
        if not np.isfinite(step).all():
            self.log(
                "RFO step contains NaN/inf; falling back to trust-region Newton step."
            )
            step = self.get_newton_step_on_trust(eigvals, eigvecs, gradient)
            if not np.isfinite(step).all():
                raise ValueError(
                    "Fallback Newton step still contains NaN/inf; "
                    "aborting to avoid corrupting coordinates."
                )
        return step

    @staticmethod
    def get_shifted_step_trans(eigvals, gradient_trans, shift):
        return -gradient_trans / (eigvals + shift)

    @staticmethod
    def get_newton_step(eigvals, eigvecs, gradient):
        if isinstance(eigvecs, torch.Tensor):
            eigvals = eigvals.to(eigvecs.device, dtype=eigvecs.dtype)
            gradient = gradient.to(eigvecs.device, dtype=eigvecs.dtype)
            return (eigvecs @ (eigvecs.T @ gradient / eigvals)).cpu().numpy()
        else:
            return eigvecs.dot(eigvecs.T.dot(gradient) / eigvals)

    def get_newton_step_on_trust(self, eigvals, eigvecs, gradient, transform=True):
        """Step on trust-radius.

        See Nocedal 4.3 Iterative solutions of the subproblem
        """
        if isinstance(eigvals, torch.Tensor):
            eigvals = eigvals.cpu().numpy()

        min_ind = eigvals.argmin()
        min_eigval = eigvals[min_ind]
        pos_definite = bool((eigvals > 0.0).all())
        if isinstance(eigvecs, torch.Tensor):
            if not isinstance(gradient, torch.Tensor):
                gradient = torch.tensor(
                    gradient, device=eigvecs.device, dtype=eigvecs.dtype
                )
            else:
                gradient = gradient.to(device=eigvecs.device, dtype=eigvecs.dtype)
            gradient_trans = eigvecs.T @ gradient
            gradient_trans = gradient_trans.cpu().numpy()
        else:
            gradient_trans = eigvecs.T.dot(gradient)

        # This will be also be True when we come close to a minimizer,
        # but then the Hessian will also be positive definite and a
        # simple Newton step will be used.
        hard_case = abs(gradient_trans[min_ind]) <= 1e-6
        self.log(f"Smallest eigenvalue: {min_eigval:.6f}")
        self.log(f"Positive definite Hessian: {pos_definite}")
        self.log(f"Hard case: {hard_case}")

        def get_step(shift):
            return -gradient_trans / (eigvals + shift)

        # Unshifted Newton step
        newton_step_trans = get_step(0.0)
        newton_norm = np.linalg.norm(newton_step_trans)

        def on_trust_radius_lin(step):
            return 1 / self.trust_radius - 1 / np.linalg.norm(step)

        def finalize_step(shift):
            step = get_step(shift)
            if transform:
                if isinstance(eigvecs, torch.Tensor):
                    step = torch.tensor(step, device=eigvecs.device, dtype=eigvecs.dtype)
                    return (eigvecs @ step).cpu().numpy()
                else:
                    return eigvecs.dot(step)
            return step

        # Simplest case. Positive definite Hessian and predicted step is
        # already in trust radius.
        if pos_definite and newton_norm <= self.trust_radius:
            self.log("Using unshifted Newton step.")
            if isinstance(eigvecs, torch.Tensor):
                newton_step_trans = torch.tensor(
                    newton_step_trans, device=eigvecs.device, dtype=eigvecs.dtype
                )
                return (eigvecs @ newton_step_trans).cpu().numpy()
            else:
                return eigvecs.dot(newton_step_trans)

        # If the Hessian is not positive definite or if the step is too
        # long we have to determine the shift parameter lambda.
        rs_kwargs = {
            "f": lambda shift: on_trust_radius_lin(get_step(shift)),
            "xtol": 1e-3,
            # Would otherwise be chosen automatically, but we set it
            # here explicitly for verbosity.
            "method": "brentq",
        }

        def root_search(bracket):
            rs_kwargs.update(
                {
                    "bracket": bracket,
                    "x0": bracket[0] + 1e-3,
                }
            )
            res = root_scalar(**rs_kwargs)
            return res

        BRACKET_END = 1e10
        if not hard_case:
            bracket_start = 0.0 if pos_definite else -min_eigval + 1e-2
            bracket = (bracket_start, BRACKET_END)
            try:
                res = root_search(bracket)
                assert res.converged
                return finalize_step(res.root)
            # ValueError may be raised when the function values for the
            # initial bracket have the same sign. If so, we continue with
            # treating it as a hard case.
            except ValueError:
                pass

        # Now we would try the bracket (-b2, -b1). The resulting step should have
        # a suitable length, but the (shifted) Hessian would have an incorrect
        # eigenvalue spectrum (not positive definite). To solve this we use a
        # different formula to calculate the step.
        mask = np.ones_like(gradient_trans)
        mask[min_ind] = 0
        mask = mask.astype(bool)
        without_min = gradient_trans[mask] / (eigvals[mask] - min_eigval)
        tau_sq = self.trust_radius**2 - (without_min**2).sum()
        if tau_sq >= 0.0:
            tau = sqrt(tau_sq)
            step_trans = [tau] + (-without_min).tolist()
        else:
            # Hard case. Search in open interval (endpoints not included)
            # (-min_eigval, inf).
            bracket = (-min_eigval + 1e-6, BRACKET_END)
            try:
                res = root_search(bracket)
                if res.converged:
                    return finalize_step(res.root)
            except ValueError:
                pass
            # Fallback: clamp tau to 0 so the step excludes the
            # minimum-eigenvalue component but remains valid.
            self.log("Hard case fallback: tau clamped to 0.")
            tau = 0.0
            step_trans = [tau] + (-without_min).tolist()

        if not transform:
            return step_trans

        if isinstance(eigvecs, torch.Tensor):
            step_trans = torch.tensor(step_trans, device=eigvecs.device, dtype=eigvecs.dtype)
            return (eigvecs @ step_trans).cpu().numpy()
        else:
            return eigvecs.dot(step_trans)

    @staticmethod
    def quadratic_model(gradient, hessian, step):
        if isinstance(gradient, torch.Tensor):
            step = torch.tensor(step, device=gradient.device, dtype=gradient.dtype)
            return (step @ gradient + 0.5 * step @ hessian @ step).cpu().numpy()
        else:
            step = np.asarray(step).ravel()
            return step.dot(gradient) + 0.5 * step.dot(hessian).dot(step)

    @staticmethod
    def rfo_model(gradient, hessian, step):
        return HessianOptimizer.quadratic_model(gradient, hessian, step) / (
            1 + step.dot(step)
        )

    def get_step_func(self, eigvals, gradient, grad_rms_thresh=1e-2):
        positive_definite = (eigvals < 0).sum() == 0
        gradient_small = rms(gradient) < grad_rms_thresh

        if self.adapt_step_func and gradient_small and positive_definite:
            return self.get_newton_step_on_trust, self.quadratic_model
        # RFO fallback
        else:
            return self.get_rs_step, self.rfo_model
