from dataclasses import replace
from typing import List, Optional

import h5py
import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import log
from pysisyphus.intcoords.augment_bonds import augment_bonds
from pysisyphus.intcoords.PrimTypes import normalize_prim_input, normalize_prim_inputs
from pysisyphus.optimizers import poly_fit
from pysisyphus.optimizers.guess_hessians import ts_hessian, HessInit
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer, HessUpdate
from pysisyphus.optimizers.Optimizer import get_data_model, get_h5_group

from pysisyphus.helpers import array2string
import torch

class TSHessianOptimizer(HessianOptimizer):
    """Optimizer to find first-order saddle points."""

    valid_updates = ("bofill", "ts_bfgs", "ts_bfgs_org", "ts_bfgs_rev")

    def __init__(
        self,
        geometry: Geometry,
        roots: Optional[List[int]] = None,
        root: int = 0,
        hessian_ref: Optional[str] = None,
        reference_mode=None,
        rx_modes=None,
        prim_coord=None,
        rx_coords=None,
        hessian_init: HessInit = "calc",
        hessian_update: HessUpdate = "bofill",
        hessian_recalc_reset: bool = True,
        max_micro_cycles: int = 50,
        trust_radius: float = 0.3,
        trust_max: float = 0.5,
        augment_bonds: bool = False,
        min_line_search: bool = False,
        max_line_search: bool = False,
        assert_neg_eigval: bool = False,
        track_mode_by_overlap: bool = False,
        reject_mode_loss: bool = True,
        mode_loss_trust_floor: float = 1e-5,
        max_mode_loss_rejections: int = 5,
        verify_saddle: bool = True,
        saddle_imaginary_threshold_cm: float = 5.0,
        saddle_recovery_step: float = 0.01,
        saddle_recovery_check_interval: int = 50,
        saddle_recovery_max_cycles: int = 200,
        max_higher_order_checks: int = 3,
        **kwargs,
    ) -> None:
        """Baseclass for transition state optimizers utilizing Hessian information.

        Several arguments expect a typed primitive or an iterable of typed primitives.
        A typed primitive is specified as (PrimType, int, int, ...), e.g., for a bond
        between atoms 0 and 1: (BOND, 0, 1) or for a bend between the atom triple 0, 1, 2
        as (BEND, 0, 1, 2).

        Parameters
        ----------
        geometry
            Geometry to be optimized.
        roots
            Indices of modes to maximize along, e.g., to optimize saddle points of 2nd order.
            Overrides 'root'.
        root
            Index of imaginary mode to maximize along. Shortcut for 'roots' with only one root.
        hessian_ref
            Filename pointing to a pysisyphus HDF5 Hessian.
        reference_mode
            Cartesian reaction direction (full 3N or active-DOF vector).  The
            initial Hessian root is selected by overlap with this direction,
            and it is retained as the preferred recovery direction if exact
            saddle validation later finds n_imag=0.  Path workflows should
            normally provide the tangent through the HEI.
        rx_modes : iterable of (iterable of (typed_prim, phase_factor))
            Select initial root(s) by overlap with a modes constructed from the given
            typed primitives with respective phase factors.
        prim_coord : typed_prim
            Select initial root/imaginary mode by overlap with this internal coordinate.
            Shortcut for 'rx_modes' with only one internal coordinate.
        rx_coords : iterable of (typed_prim)
            Construct imaginary mode comprising the given typed prims by modifying
            a model Hessian.
        hessian_init
            Type of initial model Hessian.
        hessian_update
            Type of Hessian update. Defaults to BFGS for minimizations and Bofill
            for saddle point searches.
        hessian_recalc_reset
            Whether to skip Hessian recalculation after reset. Undocumented.
        max_micro_cycles
            Maximum number of RS iterations.
        trust_radius
            Initial trust radius in whatever unit the optimization is carried out.
        trust_max
            Maximum trust radius.
        augment_bonds
            Try to derive additional streching coordinates from the imaginary mode.
        min_line_search
            Carry out line search along the imaginary mode.
        max_line_search
            Carry out line search in the subspace that is minimized.
        assert_neg_eigval
            Check for the existences for at least one significant negative eigenvalue.
            If enabled and no negative eigenvalue is present the optimization will be
            aborted.
        reject_mode_loss
            Once a negative mode has been found, reject trials that leave the
            negative-curvature region and retry with a smaller trust radius.
        mode_loss_trust_floor
            Emergency trust-radius floor for mode-loss retries.
        max_mode_loss_rejections
            Rejections allowed at the emergency floor before stopping without
            reporting convergence.
        verify_saddle
            Recalculate the exact Hessian near apparent convergence and verify
            the saddle with a mass-weighted, translation/rotation-projected
            vibrational analysis.
        saddle_imaginary_threshold_cm
            Minimum absolute imaginary frequency required by exact saddle
            verification, in cm^-1.
        saddle_recovery_step
            Minimum uphill displacement used to escape a near-flat local minimum
            exposed by exact-Hessian validation.
        saddle_recovery_check_interval
            Recovery steps between exact physical-Hessian curvature checks.
        saddle_recovery_max_cycles
            Maximum bounded recovery steps allowed to regain negative curvature.
        max_higher_order_checks
            Consecutive terminal exact checks allowed at a saddle order above
            the requested order before stopping for an explicit flatten
            restart instead of repeating exact Hessians indefinitely.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments passed to the HessianOptimizer/Optimizer baseclass.
        """

        assert (
            hessian_update in self.valid_updates
        ), f"Invalid Hessian update. Please chose from: {self.valid_updates}!"

        super().__init__(
            geometry,
            hessian_init=hessian_init,
            hessian_update=hessian_update,
            trust_radius=trust_radius,
            trust_max=trust_max,
            hessian_recalc_reset=hessian_recalc_reset,
            **kwargs,
        )

        if (root is not None) and (roots is None):
            roots = [
                root,
            ]
        elif roots is None:
            roots = list()
        self.roots = roots
        self.track_mode_by_overlap = track_mode_by_overlap
        self.log(f"{self.roots=}")
        self.log(f"{self.track_mode_by_overlap=}")
        self.hessian_ref = hessian_ref
        try:
            with h5py.File(self.hessian_ref, "r") as handle:
                self.hessian_ref = handle["hessian"][:]
            _ = self.geometry.coords.size
            expected_shape = (_, _)
            shape = self.hessian_ref.shape
            # Hessian is not yet converted to the correct coordinate system if
            # coord_type != cart.
            assert (
                self.geometry.coord_type == "cart"
            ), "hessian_ref with internal coordinates are not yet supported."
            assert shape == expected_shape, (
                f"Shape of reference Hessian {shape} doesn't match the expected "
                f"shape {expected_shape} of the Hessian for the current coordinates!"
            )
        except OSError as err:
            self.log(
                f"Tried to load reference Hessian from '{self.hessian_ref}' "
                "but the file could not be found."
            )
            self.hessian_ref = None
        except (ValueError, TypeError) as err:
            self.log(f"No reference Hessian provided.")

        if reference_mode is None:
            self.reference_mode = None
        else:
            reference_mode = np.asarray(reference_mode, dtype=float).reshape(-1)
            mode_norm = float(np.linalg.norm(reference_mode))
            if (
                reference_mode.size == 0
                or not np.all(np.isfinite(reference_mode))
                or mode_norm <= 0.0
            ):
                raise ValueError("reference_mode must be a finite non-zero vector")
            self.reference_mode = reference_mode / mode_norm
            # A path tangent supplies the missing identity of the reaction
            # direction.  Preserve the selected mode by overlap as curvature
            # changes instead of resetting unconditionally to raw root 0.
            self.track_mode_by_overlap = True

        # Select initial root according to highest contribution of 'prim_coord'
        if prim_coord is not None:
            self.log("'prim_coord' given. Overwriting/setting 'rx_modes'.")
            rx_modes = [[[prim_coord, 1.0]]]
        self.prim_coord = prim_coord

        # Construct initial imaginary mode from the given inputs while respecting
        # the given weight/phase factors.
        self.rx_modes = rx_modes

        # Construct initial imaginary mode according the primitive internals in
        # 'rx_coords' by modifying a model Hessian.
        if rx_coords is not None:
            rx_coords = normalize_prim_inputs(rx_coords)
        self.rx_coords = rx_coords

        # Bond augmentation is only useful with calculated hessians
        self.augment_bonds = augment_bonds and (self.hessian_init == "calc")
        self.min_line_search = min_line_search
        self.max_line_search = max_line_search
        self.assert_neg_eigval = assert_neg_eigval
        self.reject_mode_loss = bool(reject_mode_loss)
        self.mode_loss_trust_floor = float(mode_loss_trust_floor)
        if self.mode_loss_trust_floor <= 0.0:
            raise ValueError("mode_loss_trust_floor must be positive")
        self.max_mode_loss_rejections = int(max_mode_loss_rejections)
        if self.max_mode_loss_rejections < 1:
            raise ValueError("max_mode_loss_rejections must be at least 1")
        self.verify_saddle = bool(verify_saddle)
        self.saddle_imaginary_threshold_cm = abs(
            float(saddle_imaginary_threshold_cm)
        )
        if self.saddle_imaginary_threshold_cm <= 0.0:
            raise ValueError("saddle_imaginary_threshold_cm must be positive")
        self.saddle_recovery_step = float(saddle_recovery_step)
        if self.saddle_recovery_step <= 0.0:
            raise ValueError("saddle_recovery_step must be positive")
        self.saddle_recovery_check_interval = int(
            saddle_recovery_check_interval
        )
        if self.saddle_recovery_check_interval < 1:
            raise ValueError(
                "saddle_recovery_check_interval must be at least 1"
            )
        self.saddle_recovery_max_cycles = int(saddle_recovery_max_cycles)
        if self.saddle_recovery_max_cycles < 1:
            raise ValueError("saddle_recovery_max_cycles must be at least 1")
        self.max_higher_order_checks = int(max_higher_order_checks)
        if self.max_higher_order_checks < 1:
            raise ValueError("max_higher_order_checks must be at least 1")
        self.negative_mode_seen = False
        self.rejected_mode_loss_steps = 0
        self.mode_loss_rejections_at_floor = 0
        self.exact_saddle_checks = 0
        self.saddle_recovery_cycles = 0
        self.saddle_recovery_steps = 0
        self.higher_order_saddle_checks = 0
        self._saddle_recovery_active = False
        self._saddle_recovery_mode = None
        self._saddle_recovery_sign = None
        self._physical_ts_mode = None
        self._last_recovery_mode_curvature = None
        self._last_recovery_mode_frequency_cm = None
        self._last_exact_frequencies_cm = None
        self._last_exact_n_imaginary = None
        self._last_exact_saddle_cycle = None
        self._last_exact_cart_coords = None
        self._last_exact_saddle_verified = False
        self._last_exact_target_mode_index = None
        self._last_exact_target_mode_overlap = None
        self._last_exact_target_mode_is_negative = None
        self._last_exact_target_mode_reanchored = False
        self._last_exact_physical_mode = None
        self._initial_reference_root_mode = None
        self._initial_reference_root_index = None
        self._initial_reference_root_overlap = None
        self._initial_reference_root_eigenvalue = None
        self._best_exact_saddle = None

        self.ts_modes = list()
        self.max_micro_cycles = max_micro_cycles
        self.prim_contrib_thresh = 0.05

        self.alpha0 = 1

    @property
    def root(self):
        return self.roots[0]

    @root.setter
    def root(self, root):
        raise Exception("Setting 'self.root' is deprecated. Set 'self.roots' instead.")

    @property
    def roots(self):
        return self._roots

    @roots.setter
    def roots(self, roots):
        self._roots = np.array(roots, dtype=int)

    def prepare_opt(self, *args, **kwargs):
        if self.augment_bonds:
            self.geometry = augment_bonds(self.geometry, root=self.root)
            # Update data model and HD5 shapes, as the number of coordinates
            # may have changed.
            if self.dump:
                self.data_model = get_data_model(
                    self.geometry, self.is_cos, self.max_cycles
                )
                # self.h5_group = get_h5_group(
                #     self.h5_fn, self.h5_group_name, self.data_model
                # )

        # Calculate/set initial hessian
        super().prepare_opt(*args, **kwargs)

        # Assume a guess hessian when not calculated. This hessian has to be
        # modified according to the assumed reaction coordinates.
        if self.hessian_init != "calc" and (self.rx_coords is not None):
            assert self.geometry.coord_type != "cart", (
                "Using a modified guess Hessian for TS-optimizations is "
                "only supported in redundand internal coordinates "
                "(coord_type=redund)"
            )
            prim_inds = [
                self.geometry.internal.get_index_of_typed_prim(rxc)
                for rxc in self.rx_coords
            ]
            missing_prim_inds = [
                self.rx_coords[i] for i, _ in enumerate(prim_inds) if _ is None
            ]
            assert len(missing_prim_inds) == 0, (
                "Some of the requested reaction coordinates are not defined: "
                f"{missing_prim_inds}"
            )
            _old_H = self.H
            self.H = ts_hessian(_old_H, coord_inds=prim_inds)
            del _old_H
            if torch.cuda.is_available() and isinstance(self.H, torch.Tensor):
                torch.cuda.empty_cache()

        # Determiniation of initial mode either by using a provided
        # reference hessian, or by using a supplied root.

        if isinstance(self.H, torch.Tensor):
            eigvals, eigvecs = torch.linalg.eigh(self.H)
            eigvals = eigvals.cpu().numpy()
        else:
            eigvals, eigvecs = np.linalg.eigh(self.H)
        neg_inds = eigvals < -self.small_eigval_thresh
        self.log_negative_eigenvalues(eigvals, "Initial ")

        self.log("Determining initial TS mode to follow uphill.")
        if self.reference_mode is not None:
            ref_mode = self._reference_mode_for_eigenspace(eigvecs)
            eigvecs_np = (
                eigvecs.detach().cpu().numpy()
                if isinstance(eigvecs, torch.Tensor)
                else np.asarray(eigvecs)
            )
            overlaps = eigvecs_np.T @ ref_mode
            selected_root = self._select_initial_reference_root(eigvals, overlaps)
            self.roots = [selected_root]
            self._initial_reference_root_index = int(selected_root)
            self._initial_reference_root_overlap = float(abs(overlaps[selected_root]))
            self._initial_reference_root_eigenvalue = float(eigvals[selected_root])
            used_str = (
                "overlap-band/curvature selection with Cartesian "
                "path/reference mode"
            )
        # Select an initial TS-mode by highest overlap with eigenvectors from
        # reference Hessian.
        elif self.hessian_ref is not None:
            eigvals_ref, eigvecs_ref = np.linalg.eigh(self.hessian_ref)
            self.log_negative_eigenvalues(eigvals_ref, "Reference ")
            assert eigvals_ref[0] < -self.small_eigval_thresh
            ref_mode = eigvecs_ref[:, 0]
            overlaps = np.einsum("ij,j->i", eigvecs.T, ref_mode)
            ovlp_str = array2string(overlaps, precision=4)
            self.log(
                "Overlaps between eigenvectors of current Hessian "
                f"TS mode from reference Hessian:"
            )
            self.log(f"\t{ovlp_str}")

            root = np.abs(overlaps).argmax()
            self.roots = [
                root,
            ]
            self.log(
                "Highest overlap between reference TS mode and "
                f"eigenvector {self.root:02d}."
            )
            used_str = "overlap with reference TS mode"
        # Construct an approximate initial mode according to user input
        # and calculate overlaps with the current eigenvectors.
        elif self.rx_modes is not None:
            self.log(f"Constructing reference mode, according to user input")
            assert self.geometry.coord_type != "cart"
            int_ = self.geometry.internal
            modes = list()
            for i, rx_mode in enumerate(self.rx_modes):
                mode = np.zeros_like(self.geometry.coords)
                for typed_prim, val in rx_mode:
                    typed_prim = normalize_prim_input(typed_prim)[0]
                    ind = int_.get_index_of_typed_prim(typed_prim)
                    mode[ind] = val
                    self.log(
                        f"\tIndex {ind} (coordinate {typed_prim}) set to {val:.4f}"
                    )
                mode /= np.linalg.norm(mode)
                modes.append(mode)
                mode_str = array2string(mode, precision=2)
                self.log(f"Normalized reference mode {i:02d}: {mode_str}")

            # Calculate overlaps in non-redundant subspace by zeroing overlaps
            # in the redundant subspace.
            small_inds = np.abs(eigvals) < self.small_eigval_thresh
            # Take absolute value, because sign of eigenvectors is ambiguous.
            ovlps = np.abs(np.einsum("ij,ki->kj", eigvecs, modes))
            ovlps[:, small_inds] = 0.0
            self.roots = ovlps.argmax(axis=1)
            used_str = "overlap with user-generated mode"
        else:
            used_str = f"root(s)={self.roots}"
        self.log(f"Used {used_str} to select inital TS mode.")

        # Below, some code is found, that checks if the chosen root(s) are a
        # sensible choice, i.e., if they are negative. Currently, it is commented out,
        # as we can also start from the convex region of the PES.
        #
        # Check if the selected mode (root) is a sensible choice.
        #
        # small_eigval_thresh is positive and we dont take the absolute value
        # of the eigenvalues. So we multiply small_eigval_thresh to get a
        # negative number.
        # assert (eigvals[self.roots] < -self.small_eigval_thresh).all(), (
        # "Expected negative eigenvalue(s)! Eigenvalues of selected TS-modes "
        # f"are above the the threshold of self.small_eigval_thresh:.6e}!"
        # )

        # Select an initial TS-mode by root index. self.root may have been
        # modified by using a reference hessian.
        self.ts_modes = eigvecs[:, self.roots].T
        self.ts_mode_eigvals = eigvals[self.roots]
        if self.reference_mode is not None and len(self.ts_modes) == 1:
            initial_mode = self.ts_modes[0]
            if isinstance(initial_mode, torch.Tensor):
                initial_mode = initial_mode.detach().cpu().numpy()
            initial_mode = np.asarray(initial_mode, dtype=float).reshape(-1)
            initial_norm = float(np.linalg.norm(initial_mode))
            if np.all(np.isfinite(initial_mode)) and initial_norm > 0.0:
                # Preserve the soft, path-correlated root selected at the
                # original HEI.  It is a bounded fallback when displacements
                # along a kinked discrete MEP tangent repeatedly find only an
                # unrelated imaginary mode.
                self._initial_reference_root_mode = (
                    initial_mode.copy() / initial_norm
                )
        # A raw negative root is not yet proof that the Cartesian path mode is
        # a physical vibration (TR/frozen-coordinate artifacts can be lower).
        # Path-guided runs arm mode-loss rollback only after exact PHVA or an
        # explicit recovery crossing confirms the requested mode.
        self.negative_mode_seen = (
            False
            if self.reference_mode is not None
            else self._has_required_negative_modes(eigvals)
        )
        self.log(
            f"Using root(s) {self.roots} with eigenvalues "
            f"{array2string(self.ts_mode_eigvals, precision=6)} as TS mode.\n"
        )

    def _select_initial_reference_root(self, eigvals, signed_overlaps):
        """Choose a soft root when an MEP tangent spans near-tied modes.

        A discrete Cartesian path tangent is generally a combination of normal
        modes. A strict ``argmax(overlap)`` can choose a stiff vibration whose
        overlap beats the relevant soft mode only marginally. Restrict the
        choice to roots carrying at least 90% of the maximum overlap, then take
        the lowest-curvature member. A clearly dominant overlap is unchanged.
        """
        eigvals = np.asarray(eigvals, dtype=float).reshape(-1)
        overlaps = np.abs(np.asarray(signed_overlaps, dtype=float).reshape(-1))
        if eigvals.size != overlaps.size or overlaps.size == 0:
            raise ValueError("reference overlaps must match the Hessian spectrum")
        maximum = float(np.max(overlaps))
        if not np.isfinite(maximum) or maximum <= 0.0:
            raise ValueError("reference mode has no finite Hessian-mode overlap")
        candidates = np.flatnonzero(overlaps >= 0.90 * maximum)
        selected = int(candidates[np.argmin(eigvals[candidates])])
        self.log(
            "Initial reference-root candidates within 90% of maximum "
            f"overlap: {candidates.tolist()}; selected root {selected}."
        )
        return selected

    def _reference_mode_for_eigenspace(self, eigvecs):
        """Return the configured Cartesian reference in ``eigvecs`` space."""
        if self.reference_mode is None:
            return None
        n_eig = int(eigvecs.shape[0])
        mode = np.asarray(self.reference_mode, dtype=float).reshape(-1)
        full_size = int(self.geometry.cart_coords.size)
        candidates = (
            getattr(self.geometry, "hess_active_dof_indices", None),
            getattr(self.geometry, "active_dof_indices", None),
        )

        # The public reference is Cartesian. Expand an explicitly active-only
        # vector first, then transform its differential into the optimizer's
        # coordinate space. For redundant/DLC/TRIC coordinates this is dq=B dx;
        # comparing Cartesian dx directly with internal-coordinate Hessian
        # eigenvectors is dimensionally and physically wrong.
        cart_mode = None
        if mode.size == full_size:
            cart_mode = mode
        else:
            for inds in candidates:
                if inds is None:
                    continue
                inds = np.asarray(inds, dtype=int).reshape(-1)
                if mode.size == inds.size and (
                    inds.size == 0 or int(inds.max()) < full_size
                ):
                    cart_mode = np.zeros(full_size, dtype=float)
                    cart_mode[inds] = mode
                    break

        internal = getattr(self.geometry, "internal", None)
        B = getattr(internal, "B", None) if internal is not None else None
        if cart_mode is not None and B is not None:
            B = np.asarray(B, dtype=float)
            if B.ndim == 2 and B.shape == (n_eig, full_size):
                mode = B @ cart_mode
            else:
                mode = np.empty(0, dtype=float)
        elif cart_mode is not None and cart_mode.size != n_eig:
            reduced = None
            for inds in candidates:
                if inds is None:
                    continue
                inds = np.asarray(inds, dtype=int).reshape(-1)
                if inds.size == n_eig and (
                    inds.size == 0 or int(inds.max()) < full_size
                ):
                    reduced = cart_mode[inds]
                    break
            mode = reduced if reduced is not None else np.empty(0, dtype=float)
        elif cart_mode is not None:
            mode = cart_mode

        if mode.size != n_eig:
            raise ValueError(
                "reference_mode cannot be mapped from Cartesian coordinates "
                f"to the Hessian space ({self.reference_mode.size} vs {n_eig})"
            )
        norm = float(np.linalg.norm(mode))
        if not np.isfinite(norm) or norm <= 0.0:
            raise ValueError("reference_mode is zero in the active Hessian space")
        return mode / norm

    def _transported_target_mode_for_eigenspace(self, eigvecs):
        """Return the continuously tracked reaction mode in this space.

        ``reference_mode`` identifies the initial root at the HEI. The mode can
        rotate substantially on the way to the stationary point, so exact
        validation must compare against the overlap-transported ``ts_modes``
        vector rather than repeatedly against the fixed initial tangent. If no
        transported mode is available (for example in a focused unit call
        before ``prepare_opt``), fall back to the Cartesian reference mapping.
        """
        n_eig = int(eigvecs.shape[0])
        if self.reference_mode is not None and len(self.ts_modes):
            mode = self.ts_modes[0]
            if isinstance(mode, torch.Tensor):
                mode = mode.detach().cpu().numpy()
            mode = np.asarray(mode, dtype=float).reshape(-1)
            norm = float(np.linalg.norm(mode))
            if (
                mode.size == n_eig
                and np.all(np.isfinite(mode))
                and norm > 0.0
            ):
                return mode / norm
        return self._reference_mode_for_eigenspace(eigvecs)

    def _path_recovery_mode_for_eigenspace(self, eigvecs):
        """Return the reaction direction appropriate for an ``n_imag=0`` escape.

        Before the requested mode has ever acquired negative curvature, an MEP
        tangent generally spans several positive Hessian eigenvectors. Replacing
        it by the single largest-overlap eigenvector discards most of that path
        information. Use the complete initial tangent for this first escape.
        Once negative curvature has actually been seen, the continuously
        transported mode is the safer identity through subsequent rotations.
        """
        if self.reference_mode is not None and not self.negative_mode_seen:
            return self._reference_mode_for_eigenspace(eigvecs)
        return self._transported_target_mode_for_eigenspace(eigvecs)

    def _exact_identity_reference_for_eigenspace(self, eigvecs, *, reanchor=False):
        """Return the authoritative identity for an exact saddle check.

        Before exact PHVA has ever confirmed the requested negative mode, the
        raw Hessian root followed by the optimizer is only a numerical uphill
        direction. A coarse MEP tangent can span several positive roots, so
        transporting that single initial root would discard the path identity
        even though recovery correctly keeps using the complete tangent. Keep
        exact validation on the immutable path direction until the first
        physical crossing; only then does continuous mode transport become
        authoritative. Higher-order saddles also re-anchor to the path because
        individual roots can exchange inside the negative subspace.
        """
        if reanchor or not self.negative_mode_seen:
            return self._reference_mode_for_eigenspace(eigvecs)
        return self._transported_target_mode_for_eigenspace(eigvecs)

    def _selected_ts_modes_are_negative(self, eigvals):
        """Whether the currently selected reaction roots are all negative."""
        if len(self.roots) == 0:
            return False
        selected = eigvals[self.roots]
        negative = selected < -self.small_eigval_thresh
        if isinstance(negative, torch.Tensor):
            return bool(torch.all(negative).item())
        return bool(np.all(negative))

    def _has_required_negative_modes(self, eigvals):
        neg = eigvals < -self.small_eigval_thresh
        if isinstance(neg, torch.Tensor):
            neg_num = int(neg.sum().item())
        else:
            neg_num = int(np.count_nonzero(neg))
        return neg_num >= len(self.roots)

    def _mw_frequencies_and_modes(self):
        """Return PHVA frequencies/modes using the final ``freq`` criterion.

        The frequency helper assumes a Cartesian full or atom-wise active
        Hessian.  Other coordinate systems retain the legacy raw-Hessian
        criterion because converting their Hessian back to Cartesian space is
        outside this optimizer's responsibility.
        """
        if self.geometry.coord_type not in ("cart", "cartesian"):
            return None

        H = self.cur_H
        if H is None:
            return None

        n_atoms = len(self.geometry.atoms)
        n_full = 3 * n_atoms
        frozen = {
            int(index)
            for index in np.asarray(self.geometry.freeze_atoms, dtype=int)
            if 0 <= int(index) < n_atoms
        }
        n_active = n_full - 3 * len(frozen)
        if tuple(H.shape) not in ((n_full, n_full), (n_active, n_active)):
            self.log(
                "Skipping PHVA saddle verification for a non atom-wise "
                f"partial Hessian with shape {tuple(H.shape)}."
            )
            return None

        from pysisyphus.normal_modes import _frequencies_cm_and_modes
        from pysisyphus.tr_projection import DEFAULT_TR_PROJECTION

        if isinstance(H, torch.Tensor):
            H_in = H.detach().clone().to(dtype=torch.float64)
            device = H_in.device
        else:
            H_in = torch.as_tensor(
                np.array(H, dtype=np.float64, copy=True), dtype=torch.float64
            )
            device = H_in.device

        projection_info = {}
        freqs_cm, modes = _frequencies_cm_and_modes(
            H_in,
            list(self.geometry.atomic_numbers),
            self.geometry.cart_coords.reshape(-1, 3).copy(),
            device,
            freeze_idx=sorted(frozen) or None,
            tr_projection=getattr(
                self.geometry, "tr_projection", DEFAULT_TR_PROJECTION
            ),
            projection_info=projection_info,
        )
        self._last_rigid_projection_info = projection_info
        return np.asarray(freqs_cm, dtype=float), modes

    def _recovery_mode_from_mw(self, modes, mode_index):
        """Convert a full mass-weighted vibration to active Cartesian DOFs."""
        if modes is None or len(modes) == 0:
            return None

        from pysisyphus.constants import AMU2AU
        from pysisyphus.normal_modes import _mw_mode_to_cart, _safe_masses_amu

        mode_mw = modes[int(mode_index)]
        masses_au = torch.as_tensor(
            _safe_masses_amu(self.geometry.atomic_numbers) * AMU2AU,
            dtype=mode_mw.dtype,
            device=mode_mw.device,
        )
        mode_full = _mw_mode_to_cart(mode_mw, masses_au)
        mode = np.asarray(self.active_from_full(mode_full), dtype=float)
        expected = int(self.cur_H.shape[0])
        norm = float(np.linalg.norm(mode))
        if mode.size != expected or not np.isfinite(norm) or norm <= 0.0:
            return None
        return mode / norm

    def _fallback_recovery_mode(self, eigvecs):
        root = int(self.roots[0]) if len(self.roots) else 0
        root = min(max(root, 0), eigvecs.shape[1] - 1)
        mode = eigvecs[:, root]
        if isinstance(mode, torch.Tensor):
            return mode.detach().clone()
        return np.asarray(mode, dtype=float).copy()

    def _verify_exact_vibrational_structure(self, eigvals, eigvecs):
        """Verify exact curvature in the same space used by final TS freq."""
        try:
            frequency_data = self._mw_frequencies_and_modes()
        except Exception as err:
            self.table.print(
                "Exact PHVA saddle verification failed; stopping without "
                f"convergence ({type(err).__name__}: {err})."
            )
            self.request_stop("exact PHVA saddle verification failed")
            return False, None, True

        if frequency_data is None:
            # Internal/custom coordinate spaces cannot use the Cartesian PHVA
            # helper, but a supplied path tangent can still be transformed to
            # this Hessian space. Match it against the complete spectrum so an
            # unrelated negative root cannot stand in for a positive path mode.
            neg = eigvals < -self.small_eigval_thresh
            if isinstance(neg, torch.Tensor):
                negative_count = int(neg.sum().item())
                eigvals_np = eigvals.detach().cpu().numpy()
                eigvecs_np = eigvecs.detach().cpu().numpy()
            else:
                negative_count = int(np.count_nonzero(neg))
                eigvals_np = np.asarray(eigvals)
                eigvecs_np = np.asarray(eigvecs)

            physical_mode = None
            target_is_negative = None
            self._last_exact_target_mode_index = None
            self._last_exact_target_mode_overlap = None
            self._last_exact_target_mode_is_negative = None
            self._last_exact_target_mode_reanchored = False
            if self.reference_mode is not None:
                # Individual eigenvectors are not uniquely transportable inside
                # a multi-negative subspace: near-degenerate roots can exchange
                # identity even when cycle-to-cycle overlap is continuous.  At
                # a higher-order saddle, re-anchor the mode that must survive
                # flattening to the immutable MEP tangent. At order zero/one,
                # retain the path itself until the first physical crossing,
                # then use overlap transport so genuine rotation is followed.
                reanchor = negative_count > len(self.roots)
                reference = self._exact_identity_reference_for_eigenspace(
                    eigvecs,
                    reanchor=reanchor,
                )
                overlaps = np.abs(eigvecs_np.T @ reference)
                target_index = int(np.argmax(overlaps))
                target_overlap = float(overlaps[target_index])
                physical_mode = eigvecs_np[:, target_index].copy()
                target_is_negative = bool(
                    eigvals_np[target_index] < -self.small_eigval_thresh
                )
                if not target_is_negative:
                    physical_mode = self._path_recovery_mode_for_eigenspace(
                        eigvecs
                    )
                self._last_exact_target_mode_index = target_index
                self._last_exact_target_mode_overlap = target_overlap
                self._last_exact_target_mode_is_negative = target_is_negative
                self._last_exact_target_mode_reanchored = reanchor

            has_saddle_modes = (
                self._has_required_negative_modes(eigvals)
                if target_is_negative is None
                else target_is_negative
            )
            exact_order = negative_count == len(self.roots) and has_saddle_modes
            self._last_exact_n_imaginary = negative_count
            self._last_exact_cart_coords = self.geometry.cart_coords.copy()
            self._last_exact_saddle_verified = exact_order
            self._last_exact_saddle_cycle = (
                self.cur_cycle if exact_order else None
            )
            if exact_order:
                self._record_exact_saddle_candidate()
            if negative_count > len(self.roots) and has_saddle_modes:
                self.higher_order_saddle_checks += 1
                if self.higher_order_saddle_checks >= self.max_higher_order_checks:
                    self.request_stop(
                        "persistent higher-order saddle requires a flatten restart"
                    )
            else:
                self.higher_order_saddle_checks = 0
            self.table.print(
                "Exact optimizer-space saddle validation: "
                f"n_negative={negative_count}."
            )
            if self._last_exact_target_mode_reanchored:
                self.table.print(
                    "Re-anchored the reaction-mode identity to the MEP "
                    "tangent inside the higher-order negative subspace."
                )
            if target_is_negative is True and physical_mode is not None:
                self._remember_exact_physical_mode(physical_mode)
            return has_saddle_modes, physical_mode, physical_mode is not None

        freqs_cm, modes = frequency_data
        self._last_exact_frequencies_cm = freqs_cm.copy()
        neg_mask = freqs_cm < -self.saddle_imaginary_threshold_cm
        n_imaginary = int(np.count_nonzero(neg_mask))
        self._last_exact_n_imaginary = n_imaginary
        self._last_exact_cart_coords = self.geometry.cart_coords.copy()
        self._last_exact_target_mode_index = None
        self._last_exact_target_mode_overlap = None
        self._last_exact_target_mode_is_negative = None
        self._last_exact_target_mode_reanchored = False

        physical_mode = None
        target_is_negative = None
        if self.reference_mode is not None:
            # Do not select only among negative modes. If the intended path
            # mode has become positive while an unrelated mode remains
            # negative, an n_imag=1 count alone is the wrong reaction saddle.
            # Match the reference against the complete physical spectrum and
            # require that same matched mode to carry negative curvature.
            reanchor = n_imaginary > len(self.roots)
            reference = self._exact_identity_reference_for_eigenspace(
                eigvecs,
                reanchor=reanchor,
            )
            candidate_modes = []
            candidate_indices = []
            candidate_overlaps = []
            for mode_index, frequency in enumerate(freqs_cm):
                # Projected translations/rotations occupy a numerically null
                # subspace with arbitrary eigenvectors and cannot be a
                # chemically meaningful transported path mode.
                if abs(float(frequency)) <= 1.0e-6:
                    continue
                candidate = self._recovery_mode_from_mw(modes, mode_index)
                if candidate is None:
                    continue
                candidate = np.asarray(candidate, dtype=float)
                candidate_modes.append(candidate)
                candidate_indices.append(int(mode_index))
                candidate_overlaps.append(
                    abs(float(np.dot(reference, candidate)))
                )
            if candidate_modes:
                selected = int(np.argmax(candidate_overlaps))
                target_overlap = float(candidate_overlaps[selected])
                self._last_exact_target_mode_overlap = target_overlap
                self._last_exact_target_mode_reanchored = reanchor
                if target_overlap <= 1.0e-12:
                    # A custom/partial mode conversion may not span the
                    # supplied path vector. Preserve the explicit direction
                    # for recovery, but never let an unrelated negative root
                    # authorize convergence.
                    physical_mode = self._path_recovery_mode_for_eigenspace(
                        eigvecs
                    )
                    target_is_negative = False
                    self._last_exact_target_mode_is_negative = False
                else:
                    target_index = candidate_indices[selected]
                    target_is_negative = bool(neg_mask[target_index])
                    physical_mode = (
                        candidate_modes[selected]
                        if target_is_negative
                        else self._path_recovery_mode_for_eigenspace(eigvecs)
                    )
                    self._last_exact_target_mode_index = target_index
                    self._last_exact_target_mode_is_negative = target_is_negative

        if target_is_negative is None:
            has_saddle_modes = n_imaginary >= len(self.roots)
        else:
            has_saddle_modes = target_is_negative
        exact_order = n_imaginary == len(self.roots) and has_saddle_modes
        # Several unrelated negative modes do not make a higher-order saddle
        # for the requested reaction if the path-correlated mode itself has
        # become positive. Recover that target first; flattening at this point
        # would preserve an arbitrary wrong negative mode.
        if n_imaginary > len(self.roots) and has_saddle_modes:
            self.higher_order_saddle_checks += 1
            if self.higher_order_saddle_checks >= self.max_higher_order_checks:
                self.request_stop(
                    "persistent higher-order saddle requires a flatten restart"
                )
        else:
            self.higher_order_saddle_checks = 0
        self._last_exact_saddle_verified = exact_order
        self._last_exact_saddle_cycle = (
            self.cur_cycle if exact_order else None
        )
        if exact_order:
            self._record_exact_saddle_candidate()
        lowest = f"{float(freqs_cm[0]):+.2f} cm^-1" if freqs_cm.size else "n/a"
        self.table.print(
            "Exact PHVA saddle validation: "
            f"n_imag={n_imaginary}, lowest={lowest}."
        )
        if self._last_exact_target_mode_reanchored:
            self.table.print(
                "Re-anchored the reaction-mode identity to the MEP tangent "
                "inside the higher-order imaginary-mode subspace."
            )

        negative_indices = np.flatnonzero(neg_mask)
        if physical_mode is None and negative_indices.size:
            candidate_modes = []
            candidate_indices = []
            for mode_index in negative_indices:
                candidate = self._recovery_mode_from_mw(modes, int(mode_index))
                if candidate is not None:
                    candidate_modes.append(np.asarray(candidate, dtype=float))
                    candidate_indices.append(int(mode_index))
            if candidate_modes:
                selected = 0
                selected_overlap = None
                if self.reference_mode is not None:
                    reference = self._transported_target_mode_for_eigenspace(eigvecs)
                    overlaps = [
                        abs(float(np.dot(reference, candidate)))
                        for candidate in candidate_modes
                    ]
                    selected = int(np.argmax(overlaps))
                    selected_overlap = float(overlaps[selected])
                physical_mode = candidate_modes[selected]
                self._last_exact_target_mode_index = candidate_indices[selected]
                self._last_exact_target_mode_overlap = selected_overlap
        elif physical_mode is None and freqs_cm.size:
            physical_mode = self._recovery_mode_from_mw(
                modes, int(np.argmin(freqs_cm))
            )
        if (
            not has_saddle_modes
            and physical_mode is None
            and self.reference_mode is not None
        ):
            # At a local minimum the lowest vibration is not necessarily the
            # intended reaction coordinate.  A path tangent carries that
            # otherwise missing information and therefore takes precedence
            # for saddle recovery.
            physical_mode = self._transported_target_mode_for_eigenspace(eigvecs)
        if not has_saddle_modes and physical_mode is None:
            physical_mode = self._fallback_recovery_mode(eigvecs)
        if target_is_negative is True and physical_mode is not None:
            self._remember_exact_physical_mode(physical_mode)
        return has_saddle_modes, physical_mode, True

    def _remember_exact_physical_mode(self, mode):
        """Snapshot a PHVA/physical path mode for cross-optimizer handoff."""
        if isinstance(mode, torch.Tensor):
            mode = mode.detach().cpu().numpy()
        mode = np.asarray(mode, dtype=float).reshape(-1)
        norm = float(np.linalg.norm(mode))
        if np.all(np.isfinite(mode)) and norm > 0.0:
            self._last_exact_physical_mode = mode.copy() / norm

    def _record_exact_saddle_candidate(self):
        """Keep the best exact first-order saddle seen during this run."""
        if self.forces:
            force = np.asarray(self.forces[-1], dtype=float).reshape(-1)
        else:
            force = np.asarray(self.geometry.cart_forces, dtype=float).reshape(-1)
        score = (
            float(np.max(np.abs(force))) if force.size else float("inf"),
            float(np.sqrt(np.mean(force**2))) if force.size else float("inf"),
        )
        if self._best_exact_saddle is None or score < self._best_exact_saddle["score"]:
            self._best_exact_saddle = {
                "score": score,
                "cart_coords": self.geometry.cart_coords.copy(),
                "cycle": int(self.cur_cycle),
            }

    def _restore_best_exact_saddle(self):
        """Restore a verified n_imag=order point before a guarded failure."""
        if self._best_exact_saddle is None:
            return False
        current = np.asarray(self.geometry.cart_coords)
        best = np.asarray(self._best_exact_saddle["cart_coords"])
        if current.shape == best.shape and np.allclose(
            current, best, rtol=0.0, atol=1e-12
        ):
            return False
        self.geometry.cart_coords = best.copy()
        self.table.print(
            "Restored the best exact first-order saddle candidate from cycle "
            f"{self._best_exact_saddle['cycle']} after guarded TS failure "
            "(run remains non-converged)."
        )
        return True

    def _exact_saddle_matches_current_geometry(self):
        """Whether exact PHVA verified the geometry that is current now.

        An optimizer cycle can verify a trial and later roll it back.  A cycle
        number therefore does not identify the verified coordinates.  Keep the
        coordinate identity explicit so a stale ``n_imag=1`` result can never
        authorize convergence at a different (possibly minimum) geometry.
        """
        if (
            self._last_exact_cart_coords is None
            or not self._last_exact_saddle_verified
        ):
            return False
        current = np.asarray(self.geometry.cart_coords)
        verified = np.asarray(self._last_exact_cart_coords)
        return current.shape == verified.shape and np.allclose(
            current, verified, rtol=0.0, atol=1e-12
        )

    def _recovery_mode_has_negative_curvature(self, H, mode=None):
        """Test the target recovery direction, ignoring unrelated TR artifacts."""
        if mode is None:
            mode = self._saddle_recovery_mode
        if mode is None:
            return False
        if isinstance(H, torch.Tensor):
            mode_t = torch.as_tensor(mode, dtype=H.dtype, device=H.device)
            mode_t = mode_t / torch.linalg.norm(mode_t)
            curvature = float((mode_t @ H @ mode_t).detach().cpu())
            mode_np = mode_t.detach().cpu().numpy()
        else:
            mode_arr = np.asarray(mode, dtype=float)
            mode_arr = mode_arr / np.linalg.norm(mode_arr)
            curvature = float(mode_arr @ np.asarray(H) @ mode_arr)
            mode_np = mode_arr
        self._last_recovery_mode_curvature = curvature

        # Use the same physical frequency threshold as exact PHVA.  Testing
        # only curvature < 0 would immediately accept a sub-threshold soft
        # mode that final frequency analysis intentionally discards.
        if self.geometry.coord_type in ("cart", "cartesian"):
            from ase import units
            from pysisyphus.constants import AU2EV, BOHR2ANG
            from pysisyphus.normal_modes import _safe_masses_amu

            masses_3n = np.repeat(
                _safe_masses_amu(self.geometry.atomic_numbers), 3
            )
            masses = np.asarray(self.active_from_full(masses_3n), dtype=float)
            if masses.size == mode_np.size:
                mass_norm = float(np.dot(masses, mode_np**2))
                if mass_norm > 0.0:
                    omega2 = curvature / mass_norm
                    conversion = (
                        units._hbar
                        * 1e10
                        / np.sqrt(units._e * units._amu)
                        * np.sqrt(AU2EV)
                        / BOHR2ANG
                        / units.invcm
                    )
                    frequency_cm = float(
                        np.copysign(conversion * np.sqrt(abs(omega2)), omega2)
                    )
                    self._last_recovery_mode_frequency_cm = frequency_cm
                    return frequency_cm < -self.saddle_imaginary_threshold_cm

        return curvature < -self.small_eigval_thresh

    def _near_terminal_without_eigval_check(self):
        """Whether force/energy or plateau criteria are nearly terminal."""
        if not self.forces:
            return False
        forces = np.asarray(self.forces[-1])
        force_ok = True
        if "max_force_thresh" in self.convergence:
            force_ok &= np.abs(forces).max() <= self.max_force_thresh
        if "rms_force_thresh" in self.convergence:
            force_ok &= np.sqrt(np.mean(forces**2)) <= self.rms_force_thresh

        plateau = False
        if self.energy_plateau and len(self.energies) >= self.energy_plateau_window:
            window = self.energies[-self.energy_plateau_window :]
            plateau = float(np.max(window) - np.min(window)) < self.energy_plateau_thresh
        # Exact PHVA is the authority for saddle order and the gradient is the
        # stationarity criterion.  A quasi-Newton model may still propose a
        # large next step along a pseudo-zero mode even at a valid stationary
        # saddle, so revalidate whenever the current forces are terminal.
        return bool(force_ok or plateau)

    def _restore_hessian_trial_state(self, snapshot):
        self.H = snapshot["H"]
        self.hessian_recalc_in = snapshot["hessian_recalc_in"]
        self.adapt_norm = snapshot["adapt_norm"]
        self._sy_buffer_S = snapshot["sy_S"]
        self._sy_buffer_Y = snapshot["sy_Y"]

    def _reject_lost_mode_trial(self, snapshot):
        self.trust_radius = max(
            self.trust_radius / 4.0,
            min(self.mode_loss_trust_floor, self.trust_radius),
        )
        self._restore_hessian_trial_state(snapshot)
        self.reject_current_trial()
        self.rejected_mode_loss_steps += 1
        at_floor = self.trust_radius <= self.mode_loss_trust_floor * (1.0 + 1e-12)
        self.mode_loss_rejections_at_floor = (
            self.mode_loss_rejections_at_floor + 1 if at_floor else 0
        )
        self.table.print(
            "Rejected TS trial that lost the required saddle mode; restored the "
            f"previous saddle-side geometry (trust={self.trust_radius:.3e})."
        )
        if self.mode_loss_rejections_at_floor >= self.max_mode_loss_rejections:
            self._restore_best_exact_saddle()
            self.request_stop(
                "repeated TS mode-loss trials at the emergency trust floor"
            )

        energy = self.energies[-1]
        gradient_full = -np.asarray(self.forces[-1])
        gradient, H, eigvals, eigvecs = self._hessian_system(gradient_full)
        return energy, gradient, H, eigvals, eigvecs, True

    def _refresh_exact_saddle_model(self, gradient_full):
        self.H = self.geometry.hessian
        if self.using_active_dofs:
            self.H = self.active_hessian(self.H)
        if self.hessian_recalc is not None:
            self.hessian_recalc_in = self.hessian_recalc
        self.exact_saddle_checks += 1
        return self._hessian_system(np.asarray(gradient_full))

    def housekeeping(self):
        recovery_active_at_entry = self._saddle_recovery_active
        can_reject_trial = (
            self.reject_mode_loss
            and self.negative_mode_seen
            and not self._saddle_recovery_active
            and len(self.coords) > 1
            and len(self.steps) > 0
            and len(self.predicted_energy_changes) > 0
        )
        snapshot = None
        if can_reject_trial:
            snapshot = {
                # update_hessian replaces self.H instead of mutating it, so a
                # reference is sufficient and avoids an extra NxN copy.
                "H": self.H,
                "hessian_recalc_in": self.hessian_recalc_in,
                "adapt_norm": self.adapt_norm,
                "sy_S": list(self._sy_buffer_S),
                "sy_Y": list(self._sy_buffer_Y),
            }

        energy, gradient, H, eigvals, eigvecs, resetted = super().housekeeping()
        has_negative = (
            self._selected_ts_modes_are_negative(eigvals)
            if self.reference_mode is not None
            else self._has_required_negative_modes(eigvals)
        )
        if recovery_active_at_entry:
            has_negative = self._recovery_mode_has_negative_curvature(H)
        elif self._physical_ts_mode is not None:
            # Once a physical saddle mode has been recovered/verified, follow
            # that mode rather than allowing an unrelated raw-Hessian
            # translation/rotation artifact to satisfy the mode guard.
            has_negative = self._recovery_mode_has_negative_curvature(
                H, self._physical_ts_mode
            )

        # Validate an apparent terminal point with the exact Hessian before it
        # can be reported as a transition state.
        exact_checked = False
        if (
            self.verify_saddle
            and self._near_terminal_without_eigval_check()
            and not self._saddle_recovery_active
        ):
            gradient, H, eigvals, eigvecs = self._refresh_exact_saddle_model(
                -np.asarray(self.forces[-1])
            )
            has_negative, physical_mode, physical_checked = (
                self._verify_exact_vibrational_structure(eigvals, eigvecs)
            )
            exact_checked = True
            resetted = True
            if not has_negative:
                self._saddle_recovery_active = True
                self._saddle_recovery_mode = physical_mode
                self._saddle_recovery_sign = None
                self._physical_ts_mode = None
                self.saddle_recovery_cycles = 0
                self.table.print(
                    "Exact saddle validation found no physical imaginary mode; "
                    "continuing with saddle-recovery steps instead of accepting "
                    "a local minimum."
                )
            elif physical_checked:
                self._saddle_recovery_mode = None
                self._physical_ts_mode = physical_mode

        if has_negative:
            # For a path-guided run that started in a convex region, a
            # quasi-Newton eigenvalue can flicker negative before the physical
            # reaction mode is established. Arming mode-loss rollback on that
            # transient traps the optimizer at the emergency trust floor.
            # Require exact PHVA or an explicit recovery crossing first.
            if (
                self.reference_mode is None
                or self.negative_mode_seen
                or exact_checked
                or recovery_active_at_entry
                or self._physical_ts_mode is not None
            ):
                self.negative_mode_seen = True
            self.mode_loss_rejections_at_floor = 0
            if self._saddle_recovery_active:
                self.table.print("Recovered negative curvature; resuming TS optimization.")
                mode = self._saddle_recovery_mode
                if isinstance(mode, torch.Tensor):
                    mode = mode.detach().clone()
                elif mode is not None:
                    mode = np.asarray(mode, dtype=float).copy()
                self._physical_ts_mode = mode
            self._saddle_recovery_active = False
            self._saddle_recovery_mode = None
            self._saddle_recovery_sign = None
            self.saddle_recovery_cycles = 0
        elif self._saddle_recovery_active:
            if recovery_active_at_entry:
                self.saddle_recovery_cycles += 1
            recovery_at_limit = (
                self.saddle_recovery_cycles
                >= self.saddle_recovery_max_cycles
            )
            recovery_checkpoint = (
                recovery_active_at_entry
                and (
                    recovery_at_limit
                    or self.saddle_recovery_cycles
                    % self.saddle_recovery_check_interval
                    == 0
                )
            )
            if recovery_checkpoint:
                # Bofill updates can lag behind the actual curvature while a
                # recovery displacement crosses out of a local-minimum basin.
                # Recheck at bounded intervals so the optimizer neither misses
                # a real crossing nor runs indefinitely along a stale model.
                recovered_on_exact = False
                physical_mode = None
                if self.verify_saddle:
                    gradient, H, eigvals, eigvecs = (
                        self._refresh_exact_saddle_model(
                            -np.asarray(self.forces[-1])
                        )
                    )
                    has_negative, physical_mode, physical_checked = (
                        self._verify_exact_vibrational_structure(
                            eigvals, eigvecs
                        )
                    )
                    exact_checked = True
                    resetted = True
                    if has_negative and not self.stop_requested:
                        self.negative_mode_seen = True
                        self.mode_loss_rejections_at_floor = 0
                        self._saddle_recovery_active = False
                        self._saddle_recovery_mode = None
                        self._saddle_recovery_sign = None
                        self.saddle_recovery_cycles = 0
                        if physical_checked:
                            self._physical_ts_mode = physical_mode
                        self.table.print(
                            "Exact recovery checkpoint found negative "
                            "curvature; resuming TS optimization."
                        )
                        recovered_on_exact = True
                if (
                    not recovered_on_exact
                    and not recovery_at_limit
                    and not self.stop_requested
                ):
                    if physical_mode is not None:
                        if isinstance(physical_mode, torch.Tensor):
                            physical_mode = (
                                physical_mode.detach().cpu().numpy()
                            )
                        updated_mode = np.asarray(
                            physical_mode, dtype=float
                        ).reshape(-1)
                        previous_mode = self._saddle_recovery_mode
                        if previous_mode is not None:
                            previous_mode = np.asarray(
                                previous_mode, dtype=float
                            ).reshape(-1)
                            if (
                                previous_mode.size == updated_mode.size
                                and float(
                                    np.dot(previous_mode, updated_mode)
                                )
                                < 0.0
                            ):
                                updated_mode = -updated_mode
                        self._saddle_recovery_mode = updated_mode
                    self.table.print(
                        "Exact recovery checkpoint still has no requested "
                        "negative mode; continuing bounded recovery "
                        f"({self.saddle_recovery_cycles}/"
                        f"{self.saddle_recovery_max_cycles})."
                    )
                elif not recovered_on_exact:
                    self._restore_best_exact_saddle()
                    if not self.stop_requested:
                        self.request_stop(
                            "exact-Hessian saddle recovery did not regain "
                            "negative curvature"
                        )
        elif can_reject_trial and not exact_checked:
            return self._reject_lost_mode_trial(snapshot)

        return energy, gradient, H, eigvals, eigvecs, resetted

    def check_convergence(self, *args, **kwargs):
        # The TS optimizer owns saddle-recovery vs. energy-plateau precedence,
        # so it always suppresses the base plateau stall (allow_stall=False),
        # decides exact-saddle convergence/recovery first, and only then
        # finalizes a stall if the requested curvature is present while the
        # configured current force/step criteria fail.
        kwargs = dict(kwargs)
        outer_allow_stall = kwargs.pop("allow_stall", True)
        converged, conv_info = super().check_convergence(
            *args, allow_stall=False, **kwargs
        )
        if self._saddle_recovery_active or self.stop_requested:
            # A wrong/missing physical mode arms bounded saddle recovery, and a
            # pre-existing bounded stop already owns the outcome; neither is an
            # energy-plateau stall.
            return False, conv_info
        if self.verify_saddle:
            # A plateau or a raw/quasi-Newton root is insufficient.  Exact
            # PHVA at these exact coordinates plus converged current forces is
            # the physical definition of a stationary saddle.  The proposed
            # outgoing quasi-Newton step is not allowed to drive an already
            # verified saddle into a neighboring basin merely because a
            # pseudo-zero model mode makes that proposal large.
            exact_current_saddle = self._exact_saddle_matches_current_geometry()
            physical_criteria = all(
                (
                    conv_info.max_force_converged,
                    conv_info.rms_force_converged,
                )
            )
            converged = bool(exact_current_saddle and physical_criteria)
            if converged:
                # Exact curvature plus current force/step criteria makes a
                # repeated zero-displacement cycle unnecessary.
                conv_info = replace(conv_info, energy_converged=True)
            elif outer_allow_stall and exact_current_saddle:
                # Required curvature is present (a verified first-order saddle
                # at these exact coordinates) but the energy is flat while the
                # configured current force/step criteria fail: this is a
                # stall, never a converged saddle.
                self._maybe_request_energy_plateau_stall(
                    conv_info, curvature_ok=True
                )
        return converged, conv_info

    def apply_saddle_recovery_step(self, step):
        """Ensure a finite uphill component while escaping a local minimum."""
        if not self._saddle_recovery_active or self.stop_requested:
            return step
        mode = self._saddle_recovery_mode
        if mode is None:
            mode = self.ts_modes[0]
        amplitude = min(self.trust_radius, self.saddle_recovery_step)
        if isinstance(step, torch.Tensor):
            if not isinstance(mode, torch.Tensor):
                mode = torch.as_tensor(mode, dtype=step.dtype, device=step.device)
            mode = mode / torch.linalg.norm(mode)
            component = torch.dot(step, mode)
            comp_value = float(component.detach().cpu())
            if self._saddle_recovery_sign is None and self.steps:
                incoming = self.active_from_full(np.asarray(self.steps[-1]))
                incoming = torch.as_tensor(
                    incoming, dtype=step.dtype, device=step.device
                )
                incoming_component = float(torch.dot(incoming, mode).detach().cpu())
                if abs(incoming_component) > 1e-14:
                    sign = -1.0 if incoming_component > 0.0 else 1.0
                elif abs(comp_value) > 1e-14:
                    sign = 1.0 if comp_value > 0.0 else -1.0
                else:
                    pivot = int(torch.argmax(torch.abs(mode)).item())
                    sign = 1.0 if float(mode[pivot].detach().cpu()) >= 0.0 else -1.0
                self._saddle_recovery_sign = sign
            elif self._saddle_recovery_sign is None:
                if abs(comp_value) > 1e-14:
                    sign = 1.0 if comp_value > 0.0 else -1.0
                else:
                    pivot = int(torch.argmax(torch.abs(mode)).item())
                    sign = 1.0 if float(mode[pivot].detach().cpu()) >= 0.0 else -1.0
                self._saddle_recovery_sign = sign
            sign = self._saddle_recovery_sign
            if sign * comp_value < amplitude:
                step = step + (sign * amplitude - component) * mode
                component = torch.dot(step, mode)
                orthogonal = step - component * mode
                orthogonal_norm = torch.linalg.norm(orthogonal)
                orthogonal_limit = np.sqrt(
                    max(self.trust_radius**2 - amplitude**2, 0.0)
                )
                if float(orthogonal_norm.detach().cpu()) > orthogonal_limit:
                    if orthogonal_limit == 0.0:
                        orthogonal = torch.zeros_like(orthogonal)
                    else:
                        orthogonal = orthogonal * (
                            orthogonal_limit / orthogonal_norm
                        )
                    step = component * mode + orthogonal
        else:
            mode = np.asarray(mode, dtype=float)
            mode = mode / np.linalg.norm(mode)
            component = float(np.dot(step, mode))
            if self._saddle_recovery_sign is None and self.steps:
                incoming = np.asarray(self.active_from_full(self.steps[-1]), dtype=float)
                incoming_component = float(np.dot(incoming, mode))
                if abs(incoming_component) > 1e-14:
                    sign = -1.0 if incoming_component > 0.0 else 1.0
                elif abs(component) > 1e-14:
                    sign = 1.0 if component > 0.0 else -1.0
                else:
                    pivot = int(np.argmax(np.abs(mode)))
                    sign = 1.0 if mode[pivot] >= 0.0 else -1.0
                self._saddle_recovery_sign = sign
            elif self._saddle_recovery_sign is None:
                if abs(component) > 1e-14:
                    sign = 1.0 if component > 0.0 else -1.0
                else:
                    pivot = int(np.argmax(np.abs(mode)))
                    sign = 1.0 if mode[pivot] >= 0.0 else -1.0
                self._saddle_recovery_sign = sign
            sign = self._saddle_recovery_sign
            if sign * component < amplitude:
                step = step + (sign * amplitude - component) * mode
                component = float(np.dot(step, mode))
                orthogonal = step - component * mode
                orthogonal_norm = float(np.linalg.norm(orthogonal))
                orthogonal_limit = np.sqrt(
                    max(self.trust_radius**2 - amplitude**2, 0.0)
                )
                if orthogonal_norm > orthogonal_limit:
                    if orthogonal_limit == 0.0:
                        orthogonal = np.zeros_like(orthogonal)
                    else:
                        orthogonal *= orthogonal_limit / orthogonal_norm
                    step = component * mode + orthogonal
        self.saddle_recovery_steps += 1
        return step

    def update_ts_mode(self, eigvals, eigvecs):
        neg_eigval_inds = eigvals < -self.small_eigval_thresh
        neg_num = neg_eigval_inds.sum()
        self.log_negative_eigenvalues(eigvals)

        # A recovered/PHVA-verified vibration is the chemically meaningful TS
        # direction. Select the raw-Hessian eigenvector with maximum overlap
        # among significant negative roots, so a spurious TR root cannot take
        # over merely because it has the lowest unprojected eigenvalue.
        if self._physical_ts_mode is not None and len(self.roots) == 1:
            reference = self._physical_ts_mode
            if isinstance(eigvecs, torch.Tensor):
                reference = torch.as_tensor(
                    reference, dtype=eigvecs.dtype, device=eigvecs.device
                )
                if reference.shape[0] != eigvecs.shape[0]:
                    reference = self.active_from_full(reference)
                overlaps = torch.abs(eigvecs.T @ reference)
                negative = eigvals < -self.small_eigval_thresh
                if self.reference_mode is not None:
                    selected = int(torch.argmax(overlaps).item())
                elif bool(negative.any()):
                    overlaps = torch.where(
                        negative, overlaps, torch.full_like(overlaps, -1.0)
                    )
                    selected = int(torch.argmax(overlaps).item())
                else:
                    selected = None
            else:
                reference = np.asarray(reference, dtype=float)
                if reference.shape[0] != eigvecs.shape[0]:
                    reference = self.active_from_full(reference)
                overlaps = np.abs(eigvecs.T @ reference)
                negative = eigvals < -self.small_eigval_thresh
                if self.reference_mode is not None:
                    selected = int(np.argmax(overlaps))
                else:
                    selected = (
                        int(np.argmax(np.where(negative, overlaps, -1.0)))
                        if np.any(negative)
                        else None
                    )
            if selected is not None:
                self.roots = np.array([selected], dtype=int)
                self.ts_modes = eigvecs[:, self.roots].T
                self.ts_mode_eigvals = eigvals[self.roots]
                if self.reference_mode is not None:
                    mode = self.ts_modes[0]
                    self._physical_ts_mode = (
                        mode.detach().clone()
                        if isinstance(mode, torch.Tensor)
                        else np.asarray(mode, dtype=float).copy()
                    )
                self.log(
                    "Following PHVA-verified TS mode at raw-Hessian root "
                    f"{selected}."
                )
                return

        # --- DEBUG: dump all negative eigvals + overlaps with prev TS mode ---
        # all_freqs = self._all_mw_freqs_cm()
        # if isinstance(eigvals, torch.Tensor):
        #     evs_np = eigvals.cpu().numpy()
        # else:
        #     evs_np = np.asarray(eigvals)
        # n_show = min(int(neg_num) + 2, len(evs_np))
        # neg_evs_str = ", ".join(f"{float(evs_np[i]):+.3e}" for i in range(n_show))
        # neg_cm_str = ", ".join(f"{all_freqs[i]:+.1f}" for i in range(min(n_show, len(all_freqs))))
        # print(f"[ts-mode] cycle neg_count={neg_num}")
        # print(f"[ts-mode]   eigvals (au)  [0:{n_show}] = [{neg_evs_str}]")
        # print(f"[ts-mode]   freqs   (cm⁻¹)[0:{n_show}] = [{neg_cm_str}]")
        #
        # # Overlaps with previous TS mode (if exists) for all candidates
        # if self.ts_modes is not None:
        #     try:
        #         prev_mode = self.ts_modes[0]
        #         if isinstance(eigvecs, torch.Tensor):
        #             if isinstance(prev_mode, torch.Tensor):
        #                 pm = prev_mode
        #                 if pm.shape[0] != eigvecs.shape[0]:
        #                     pm = torch.from_numpy(self.active_from_full(pm.cpu().numpy())).to(eigvecs.device, eigvecs.dtype)
        #             else:
        #                 pm_arr = prev_mode if prev_mode.shape[0] == eigvecs.shape[0] else self.active_from_full(prev_mode)
        #                 pm = torch.from_numpy(np.asarray(pm_arr, dtype=float)).to(eigvecs.device, eigvecs.dtype)
        #             ovlps_all = torch.abs(eigvecs.T @ pm).cpu().numpy()
        #         else:
        #             pm = prev_mode if prev_mode.shape[0] == eigvecs.shape[0] else self.active_from_full(prev_mode)
        #             ovlps_all = np.abs(eigvecs.T @ pm)
        #         n_ovlp = min(n_show, len(ovlps_all))
        #         ovlp_str = ", ".join(f"{float(ovlps_all[i]):.3f}" for i in range(n_ovlp))
        #         print(f"[ts-mode]   overlap w/ prev[0:{n_ovlp}] = [{ovlp_str}]")
        #         best = int(np.argmax(ovlps_all))
        #         print(f"[ts-mode]   best overlap at root={best} ({ovlps_all[best]:.4f}), "
        #               f"freq={all_freqs[best] if best < len(all_freqs) else 0:+.1f} cm⁻¹")
        #     except Exception as _e:
        #         print(f"[ts-mode]   overlap calc failed: {_e}")
        # --- END DEBUG ---

        if not self.track_mode_by_overlap:
            # Fixed-root mode: always follow root=0 (most negative eigenvalue).
            # Robust for Cartesian-coordinate TS optimization with MLIPs.
            self.roots = np.array([0] * len(self.roots), dtype=int)
            self.ts_modes = eigvecs[:, self.roots].T
            self.ts_mode_eigvals = eigvals[self.roots]
            # _cm1 = self._lowest_mw_freq_cm()
            # print(f"[ts-mode] SELECTED root=0  eigval={float(eigvals[0]):+.6e}  {_cm1:+.1f} cm⁻¹")
            return

        # --- Overlap-based mode tracking (track_mode_by_overlap=True) ---
        if (
            self.reference_mode is None
            and (self.ts_mode_eigvals < 0).all()
            and neg_num > 0
        ):
            infix = "imaginary "
            ovlp_eigvecs = eigvecs[:, :neg_num]
            eigvals = eigvals[:neg_num]
        elif self.assert_neg_eigval and neg_num == 0:
            raise AssertionError(
                "Need at least 1 negative eigenvalue for TS optimization."
            )
        else:
            infix = ""
            ovlp_eigvecs = eigvecs

        if (
            self.using_active_dofs
            and self.ts_modes is not None
            and self.ts_modes.shape[1] != ovlp_eigvecs.shape[0]
        ):
            self.log("Projecting previous TS mode(s) onto active DOFs.")
            if isinstance(self.ts_modes, torch.Tensor):
                self.ts_modes = torch.stack(
                    [self.active_from_full(mode) for mode in self.ts_modes]
                )
            else:
                self.ts_modes = np.stack(
                    [self.active_from_full(mode) for mode in self.ts_modes]
                )

        self.log(f"Overlaps of previous TS mode with current {infix}mode(s):")
        if isinstance(eigvecs, torch.Tensor):
            ovlps = torch.abs(torch.einsum("ij,ki->kj", ovlp_eigvecs, self.ts_modes))
            ovlps = ovlps.cpu().numpy()
        else:
            ovlps = np.abs(np.einsum("ij,ki->kj", ovlp_eigvecs, self.ts_modes))
        for i, ovlp in enumerate(ovlps):
            self.log(f"\tTS mode {i:02d}: {array2string(ovlp, precision=3)}")
        max_ovlp_inds = np.argmax(ovlps, axis=1)
        for i in range(len(self.ts_modes)):
            max_ovlp_ind = int(max_ovlp_inds[i])
            self.log(
                f"Mode {i}: highest overlap: {ovlps[i, max_ovlp_ind]:.6f} with mode "
                f"{max_ovlp_ind}"
            )
        self.roots = max_ovlp_inds
        self.ts_modes = ovlp_eigvecs.T[self.roots]
        self.ts_mode_eigvals = eigvals[self.roots]
        # for i, ev in enumerate(self.ts_mode_eigvals):
        #     selected_root = int(self.roots[i])
        #     selected_freq = all_freqs[selected_root] if selected_root < len(all_freqs) else 0.0
        #     print(
        #         f"[ts-mode] SELECTED root={selected_root:3d}  eigval={float(ev):+.6e}  "
        #         f"freq={selected_freq:+.1f} cm⁻¹  "
        #         f"overlap={float(ovlps[i, int(max_ovlp_inds[i])]):.4f}"
        #     )

    def _lowest_mw_freq_cm(self):
        """Compute the lowest frequency (cm⁻¹) from current Hessian using
        the same _frequencies_cm_and_modes as the freq command."""
        freqs = self._all_mw_freqs_cm()
        return freqs[0] if len(freqs) > 0 else 0.0

    def _all_mw_freqs_cm(self):
        """Return all mass-weighted frequencies (cm⁻¹) from current Hessian."""
        try:
            result = self._mw_frequencies_and_modes()
            if result is None:
                return []
            freqs_cm, _ = result
            return [float(f) for f in freqs_cm]
        except Exception as err:
            self.log(f"Mass-weighted frequency analysis failed: {err}")
            return []

    @staticmethod
    def do_line_search(e0, e1, g0, g1, prev_step, maximize, logger=None):
        poly_fit_kwargs = {
            "e0": e0,
            "e1": e1,
            "g0": g0,
            "g1": g1,
            "maximize": maximize,
        }
        if not maximize:
            poly_fit_kwargs.update(
                {
                    "g0": prev_step.dot(g0),
                    "g1": prev_step.dot(g1),
                }
            )
        prefix = "Max" if maximize else "Min"

        fit_result = poly_fit.quartic_fit(**poly_fit_kwargs)
        fit_energy = None
        fit_grad = None
        fit_step = None
        if fit_result and (0.0 < fit_result.x <= 2.0):
            x = fit_result.x
            log(logger, f"{prefix}-subpsace interpolation succeeded: x={x:.6f}")
            fit_energy = fit_result.y
            fit_step = (1 - x) * -prev_step
            fit_grad = (1 - x) * g0 + x * g1
        return fit_energy, fit_grad, fit_step

    def step_and_grad_from_line_search(
        self,
        energy,
        gradient_trans,
        eigvecs,
        min_indices,
        max_indices,
    ):
        ip_step = np.zeros_like(gradient_trans)
        ip_gradient_trans = gradient_trans.copy()

        # ``prev_*`` are consumed by BOTH the max and the min line search below,
        # so compute them once above both blocks; otherwise min_line_search=True
        # with max_line_search=False references undefined names (NameError).
        if (self.max_line_search or self.min_line_search) and self.cur_cycle > 0:
            prev_energy = self.energies[-2]
            prev_gradient = -self.forces[-2]
            try:
                prev_gradient_trans = eigvecs.T.dot(prev_gradient)
                prev_step = self.steps[-1]
                prev_step_trans = eigvecs.T.dot(prev_step)
            # Will be raised when coordinates were rebuilt and the array shapes differe.
            except ValueError:
                return ip_step, ip_gradient_trans

        if self.max_line_search and self.cur_cycle > 0:
            # Max subspace
            # max_energy, max_gradient, max_step = self.do_max_line_search(
            max_energy, max_gradient, max_step = self.do_line_search(
                prev_energy,
                energy,
                prev_gradient_trans[max_indices],
                gradient_trans[max_indices],
                prev_step=prev_step_trans[max_indices],
                maximize=True,
                logger=self.logger,
            )
            if max_gradient is not None:
                ip_gradient_trans[max_indices] = max_gradient
                ip_step[max_indices] = max_step

        if self.min_line_search and self.cur_cycle > 0:
            # Min subspace
            # min_energy, min_gradient, min_step = self.do_min_line_search(
            min_energy, min_gradient, min_step = self.do_line_search(
                prev_energy,
                energy,
                prev_gradient_trans[min_indices],
                gradient_trans[min_indices],
                prev_step=prev_step_trans[min_indices],
                maximize=False,
                logger=self.logger,
            )
            if min_gradient is not None:
                ip_gradient_trans[min_indices] = min_gradient
                ip_step[min_indices] = min_step
        return ip_step, ip_gradient_trans
