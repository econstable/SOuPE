"""Coefficient conversion for one explicit real-Wybourne c/s convention.

Even ranks 2, 4 and 6 follow the Amelin et al. 2020, Appendix B, and
Alexanian et al. 2023 supplementary material, Eq. (16), with supporting
thesis Eqs. (III.18-19).
Published factors and the general extension are distinguished in
docs/CONVERSIONS.md.
The relation is B[k,q] = lambda[k,q] * theta[k] * L[k,q].
"""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite, sqrt
from numbers import Integral, Real

__all__ = ["real_wybourne_to_stevens"]

# Each row is ordered by q = 0, ..., k. These multiply COEFFICIENTS;
# they must not be inverted as if they were operator-basis conversions.
_LAMBDA = {
    2: (1 / 2, sqrt(6), sqrt(6) / 2),
    4: (1 / 8, sqrt(5) / 2, sqrt(10) / 4, sqrt(35) / 2, sqrt(70) / 8),
    6: (
        1 / 16,
        sqrt(42) / 8,
        sqrt(105) / 16,
        sqrt(105) / 8,
        3 * sqrt(14) / 16,
        3 * sqrt(77) / 8,
        sqrt(231) / 16,
    ),
}


def real_wybourne_to_stevens(
    coefficients: Mapping[tuple[int, int, str | None], Real],
    *,
    theta: Mapping[int, Real],
) -> dict[tuple[int, int, str | None], float]:
    """Convert real-Wybourne L coefficients to tesseral Stevens B.

    Parameters
    ----------
    coefficients : mapping of (int, int, str or None) to real
        Real-Wybourne L coefficients in the convention defined below.
        Keys are
        (k, q, component): k in {2, 4, 6}, 0 <= q <= k, with None for
        q = 0 and 'c' or 's' otherwise. Values must be finite.
    theta : mapping of int to real
        Explicit Stevens factors theta_2, theta_4, theta_6 for every rank
        present in coefficients. Values must be finite; zero is allowed.
        For input L values already multiplied by theta, explicitly pass
        1 for those ranks. No ion factors are inferred from J or ion names.

    Returns
    -------
    dict of (int, int, str or None) to float
        B[k,q,component] = lambda[k,q] * theta[k] * L[k,q,component].
        Keys, zero entries, and energy units are preserved. Neither input
        is modified. Omitted components remain omitted.

    Raises
    ------
    TypeError
        Inputs are not mappings, keys have invalid types, or values are
        not real numbers. Booleans and complex values are not accepted.
    ValueError
        A rank, order, or component is unsupported, a rank factor is
        missing, or a value or converted result is nonfinite.

    Notes
    -----
    This convention uses real cosine/sine combinations without 1/sqrt(2).
    For q > 0, W_c = C_-q + (-1)**q C_q and
    W_s = 1j * (C_-q - (-1)**q C_q), with C_q in the Condon-Shortley
    phase convention. For q = 0, W_0 = C_0. These follow the combinations
    of Alexanian et al. 2023 supplementary material Eq. (16), expressed using
    explicit c/s components. External signed labels require phase checks.
    It does not accept complex spherical Wybourne D coefficients or the
    alternative orthonormal real combinations. Both expansions refer to
    the same coordinate frame and fixed-J manifold; no rotation or radial
    integral calculation is performed. A zero theta removes that rank's
    contribution, so this projection need not be invertible.

    Examples
    --------
    Illustrative factors, not parameters for a particular ion:

    >>> real_wybourne_to_stevens({(2, 0, None): 10.0}, theta={2: -0.01})
    {(2, 0, None): -0.05}
    """
    if not isinstance(coefficients, Mapping) or not isinstance(theta, Mapping):
        raise TypeError("coefficients and theta must be mappings.")
    factors = {}
    for rank, value in theta.items():
        _validate_rank(rank)
        factors[rank] = _finite_real(value, f"theta[{rank}]")

    converted = {}
    for key, value in coefficients.items():
        if not isinstance(key, tuple) or len(key) != 3:
            raise TypeError("Coefficient keys must be (k, q, component).")
        rank, order, component = key
        _validate_rank(rank)
        if not isinstance(order, Integral) or isinstance(order, bool):
            raise TypeError("q must be an integer.")
        if not 0 <= order <= rank:
            raise ValueError("q must be between 0 and k.")
        if order == 0 and component is not None:
            raise ValueError("component must be None for q = 0.")
        if order > 0 and component not in ("c", "s"):
            raise ValueError("component must be 'c' or 's' for q > 0.")
        if rank not in factors:
            raise ValueError(f"Missing theta factor for rank {rank}.")
        coefficient = _finite_real(value, f"coefficient {key}")
        result = _LAMBDA[rank][order] * factors[rank] * coefficient
        if not isfinite(result):
            raise ValueError(f"Converted coefficient {key} is nonfinite.")
        converted[key] = result
    return converted


def _validate_rank(rank: int) -> None:
    """Limit conversion to the documented even crystal-field ranks."""
    if not isinstance(rank, Integral) or isinstance(rank, bool):
        raise TypeError("k must be an integer.")
    if rank not in _LAMBDA:
        raise ValueError("Conversion supports only ranks 2, 4, and 6.")


def _finite_real(value: Real, label: str) -> float:
    """Reject ambiguous complex input and invalid numerical values."""
    if not isinstance(value, Real) or isinstance(value, bool):
        raise TypeError(f"{label} must be a real number.")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{label} must be finite.")
    return result
