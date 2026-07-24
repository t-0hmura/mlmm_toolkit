"""Lightweight canonicalization of calculator method names."""

from __future__ import annotations

from typing import Optional


def normalize_hessian_calc_mode(value: Optional[str]) -> str:
    """Return the canonical Hessian method or reject an unknown token."""

    text = "FiniteDifference" if value is None else str(value).strip()
    canonical = {
        "finitedifference": "FiniteDifference",
        "analytical": "Analytical",
    }.get(text.casefold())
    if canonical is None:
        raise ValueError(
            "Unsupported hessian_calc_mode "
            f"{value!r}. Choose from: FiniteDifference, Analytical."
        )
    return canonical


def normalize_mm_hessian_mode(
    value: Optional[str], *, mm_fd: bool = True
) -> str:
    """Resolve the low-level MM Hessian policy."""

    if value is None:
        return "finite_difference" if bool(mm_fd) else "analytical"
    text = str(value).strip().casefold().replace("-", "_")
    canonical = {
        "finitedifference": "finite_difference",
        "finite_difference": "finite_difference",
        "fd": "finite_difference",
        "analytical": "analytical",
        "analytic": "analytical",
        "none": "none",
    }.get(text)
    if canonical is None:
        raise ValueError(
            "Unsupported mm_hessian_mode "
            f"{value!r}. Choose from: finite_difference, analytical, none."
        )
    return canonical


def normalize_link_atom_method(value: Optional[str]) -> str:
    """Return the canonical link-atom placement method or reject a typo."""

    text = "scaled" if value is None else str(value).strip()
    canonical = {"scaled": "scaled", "fixed": "fixed"}.get(text.casefold())
    if canonical is None:
        raise ValueError(
            "Unsupported link_atom_method "
            f"{value!r}. Choose from: scaled, fixed."
        )
    return canonical


__all__ = [
    "normalize_hessian_calc_mode",
    "normalize_link_atom_method",
    "normalize_mm_hessian_mode",
]
