"""
Spin-J generators and cached polynomial calculations.
JSet retains its Qobjs; OperatorSet returns independent copies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isfinite
from numbers import Real

from qutip import Qobj, jmat, qeye, qeye_like

__all__ = ["JSet", "jset", "jset_from_ops"]


@dataclass(frozen=True, eq=False)
class JSet:
    """Container for total-J operators in one coordinate frame.

    Construct with `jset` or `jset_from_ops`. Operators use units where
    hbar = 1. Treat the operators and scalar attributes as read-only;
    fields cannot be reassigned. QuTiP objects themselves remain mutable,
    so copy them before editing to avoid invalidating cached results.

    Attributes
    ----------
    J : float
        Nonnegative integer or half-integer angular momentum.
    Jx, Jy, Jz : qutip.Qobj
        Cartesian angular momentum operators.
    Jp, Jm : qutip.Qobj
        Raising and lowering operators, Jx +/- 1j * Jy.
    JJ : float
        J * (J + 1).
    I : qutip.Qobj
        Identity with the same dimensions as the supplied operators.
    """

    J: float
    Jx: Qobj
    Jy: Qobj
    Jz: Qobj
    Jp: Qobj
    Jm: Qobj
    JJ: float  # J(J+1) as a scalar
    I: Qobj  # noqa: E741 -- established mathematical identity attribute
    _cache: dict[tuple, Qobj] = field(default_factory=dict, repr=False)


def jset(J: float, *, dtype=None) -> JSet:
    """Build angular momentum operators using QuTiP's `jmat`.

    Parameters
    ----------
    J : float
        Nonnegative integer or half-integer angular momentum.
    dtype : str or type, optional
        QuTiP storage backend passed to jmat or qeye.

    Returns
    -------
    JSet
        Operators in the ordered basis m = J, J-1, ..., -J, with hbar = 1.
    """
    _validate_j(J)
    # The J=0 representation has zero generators and unit identity.
    if J == 0:
        zero = 0 * qeye(1, dtype=dtype)
        return jset_from_ops(J, zero, zero.copy(), zero.copy())
    return jset_from_ops(J, *jmat(J, dtype=dtype))


def jset_from_ops(J: float, Jx: Qobj, Jy: Qobj, Jz: Qobj) -> JSet:
    """Build a container from existing angular momentum operators.

    Parameters
    ----------
    J : float
        Nonnegative integer or half-integer angular momentum.
    Jx, Jy, Jz : qutip.Qobj
        Hermitian operators for the same spin-J representation, each with
        shape (2*J+1, 2*J+1) and matching QuTiP dimensions.

    Returns
    -------
    JSet
        Container preserving the supplied Qobj instances and dimensions,
        with a fresh operator cache.

    Notes
    -----
    Shape, dimensions, and Hermiticity are validated. The caller is
    responsible for the angular momentum commutators and Casimir identity.
    """
    _validate_j(J)
    d = int(2 * J + 1)
    for op in (Jx, Jy, Jz):
        if not isinstance(op, Qobj):
            raise TypeError("Jx, Jy, and Jz must be QuTiP Qobj operators.")
        if not op.isoper or op.shape != (d, d):
            raise ValueError("Operators must have shape (2*J+1, 2*J+1).")
        if op.dims != Jx.dims:
            raise ValueError("Operators must have matching QuTiP dimensions.")
        if not op.isherm:
            raise ValueError("Angular momentum operators must be Hermitian.")
    identity = qeye_like(Jx)
    Jp = Jx + 1j * Jy
    Jm = Jx - 1j * Jy
    return JSet(
        J=J,
        Jx=Jx,
        Jy=Jy,
        Jz=Jz,
        Jp=Jp,
        Jm=Jm,
        JJ=J * (J + 1),
        I=identity,
        _cache={},
    )


def _validate_j(J: float) -> None:
    """Require a finite, nonnegative integer or half-integer spin."""
    if not isinstance(J, Real) or isinstance(J, bool):
        raise TypeError("J must be a real number.")
    if not isfinite(J) or J < 0 or not float(2 * J).is_integer():
        raise ValueError("J must be a nonnegative integer or half-integer.")
