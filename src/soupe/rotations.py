"""Euler rotations using U=exp(-i alpha Jz) exp(-i beta Jy) exp(-i gamma Jz).

Angles describe intrinsic ZYZ rotations. Apply U O U.dag() to an operator
or U ket to a state. Every exponential uses the original spin generators.
"""

from __future__ import annotations

from math import isfinite, radians
from numbers import Real

from qutip import Qobj

from ._jset import JSet, jset

__all__ = ["rotation_operator"]


def _angles(alpha, beta, gamma, degrees):
    """Reject ambiguous or nonfinite angles before constructing a unitary."""
    if not isinstance(degrees, bool):
        raise TypeError("degrees must be True or False.")
    values = []
    for name, angle in zip(("alpha", "beta", "gamma"), (alpha, beta, gamma)):
        if not isinstance(angle, Real) or isinstance(angle, bool):
            raise TypeError(f"{name} must be a real angle.")
        if not isfinite(angle):
            raise ValueError(f"{name} must be finite.")
        values.append(radians(angle) if degrees else float(angle))
    return values


def _unitary(ops: JSet, alpha, beta, gamma, *, degrees=False):
    """Use this set's generators, preserving its basis and dimensions."""
    alpha, beta, gamma = _angles(alpha, beta, gamma, degrees)
    # QuTiP scalar Qobjs need no exponential in the trivial J=0 manifold.
    if ops.J == 0:
        return ops.I.copy()
    return (
        (-1j * alpha * ops.Jz).expm()
        * (-1j * beta * ops.Jy).expm()
        * (-1j * gamma * ops.Jz).expm()
    )


def rotation_operator(J, alpha, beta, gamma, *, degrees=False) -> Qobj:
    """Construct the intrinsic ZYZ rotation unitary for spin J.

    Parameters
    ----------
    J : real
        Finite nonnegative integer or half-integer angular momentum.
    alpha, beta, gamma : real
        Angles in U=exp(-i alpha Jz) exp(-i beta Jy) exp(-i gamma Jz).
        Radians unless degrees=True. All generators are unrotated.
    degrees : bool, optional
        Interpret all three angles as degrees when True.

    Returns
    -------
    qutip.Qobj
        Unitary in descending m=J,...,-J order, hbar=1. Apply U ket
        to a state, U O U.dag() to an operator, or U.dag() for the inverse.

    Raises
    ------
    TypeError
        Angles are not real numbers or degrees is not a boolean.
    ValueError
        J is invalid or an angle is nonfinite.

    Notes
    -----
    At beta=0 or pi the first and third angles are not independent
    (gimbal lock). Supplied angles still give a valid rotation; no angle
    extraction or automatic warning is performed.
    """
    return _unitary(jset(J), alpha, beta, gamma, degrees=degrees)
