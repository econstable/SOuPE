"""Explicit extended Stevens polynomials, ranks 0 through 6.

The expressions follow Rafik Ballou's
[Stevens operator table](Stevens_Operators-Ballou.pdf); see docs/CONVENTIONS.md
q is nonnegative; c/s label real tesseral cosine/sine components,
and hbar = 1. X below denotes J(J+1) I. Ion-dependent Stevens factors belong
in coefficients, not these operators.

Keep expressions explicit so the LaTeX exporter can inspect their syntax.
JSet caches calculations; OperatorSet returns independent Qobj copies.
"""

from __future__ import annotations

from numbers import Integral

from qutip import Qobj

from ._jset import JSet

__all__ = ["stevens_operator"]


def _sym(A: Qobj, B: Qobj) -> Qobj:
    """Return half the anticommutator: (AB + BA)/2.

    A tabulated (AB + BA)/4 therefore requires 0.5 * _sym(A, B).
    """
    return 0.5 * (A * B + B * A)


# k = 0
def O_0_0(ops: JSet) -> Qobj:
    return ops.I


# k = 1
def O_1_0(ops: JSet) -> Qobj:
    return ops.Jz


def O_1_1_c(ops: JSet) -> Qobj:
    return 0.5 * (ops.Jp + ops.Jm)


def O_1_1_s(ops: JSet) -> Qobj:
    return (-0.5j) * (ops.Jp - ops.Jm)


# k = 2
def O_2_0(ops: JSet) -> Qobj:
    return 3 * (ops.Jz**2) - ops.JJ * ops.I


def O_2_1_c(ops: JSet) -> Qobj:
    return _sym(ops.Jz, ops.Jx)


def O_2_1_s(ops: JSet) -> Qobj:
    return _sym(ops.Jz, ops.Jy)


def O_2_2_c(ops: JSet) -> Qobj:
    return ops.Jx**2 - ops.Jy**2


def O_2_2_s(ops: JSet) -> Qobj:
    return 2 * _sym(ops.Jx, ops.Jy)


# k = 3
def O_3_0(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    return 5 * (ops.Jz**3) - (3 * X - ops.I) * ops.Jz


def O_3_1_c(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 5 * (ops.Jz**2) - X - 0.5 * ops.I
    Jsum1 = ops.Jp + ops.Jm
    return 0.5 * _sym(T, Jsum1)


def O_3_1_s(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 5 * (ops.Jz**2) - X - 0.5 * ops.I
    Jdif1 = ops.Jp - ops.Jm
    return (-0.5j) * _sym(T, Jdif1)


def O_3_2_c(ops: JSet) -> Qobj:
    Jsum2 = (ops.Jp**2) + (ops.Jm**2)
    return 0.5 * _sym(ops.Jz, Jsum2)


def O_3_2_s(ops: JSet) -> Qobj:
    Jdif2 = (ops.Jp**2) - (ops.Jm**2)
    return (-0.5j) * _sym(ops.Jz, Jdif2)


def O_3_3_c(ops: JSet) -> Qobj:
    return 0.5 * ((ops.Jp**3) + (ops.Jm**3))


def O_3_3_s(ops: JSet) -> Qobj:
    return (-0.5j) * ((ops.Jp**3) - (ops.Jm**3))


# k = 4
def O_4_0(ops: JSet) -> Qobj:
    Jz2 = ops.Jz**2
    JJ = ops.JJ
    return (
        35 * (Jz2**2)
        - (30 * JJ - 25) * Jz2
        + 3 * (JJ**2) * ops.I
        - 6 * JJ * ops.I
    )


def O_4_1_c(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 7 * (ops.Jz**3) - (3 * X + ops.I) * ops.Jz
    Jsum1 = ops.Jp + ops.Jm
    return 0.5 * _sym(T, Jsum1)


def O_4_1_s(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 7 * (ops.Jz**3) - (3 * X + ops.I) * ops.Jz
    Jdif1 = ops.Jp - ops.Jm
    return (-0.5j) * _sym(T, Jdif1)


def O_4_2_c(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 7 * (ops.Jz**2) - X - 5 * ops.I
    Jsum2 = (ops.Jp**2) + (ops.Jm**2)
    return 0.5 * _sym(T, Jsum2)


def O_4_2_s(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 7 * (ops.Jz**2) - X - 5 * ops.I
    Jdif2 = (ops.Jp**2) - (ops.Jm**2)
    return (-0.5j) * _sym(T, Jdif2)


def O_4_3_c(ops: JSet) -> Qobj:
    Jsum3 = (ops.Jp**3) + (ops.Jm**3)
    return 0.5 * _sym(ops.Jz, Jsum3)


def O_4_3_s(ops: JSet) -> Qobj:
    Jdif3 = (ops.Jp**3) - (ops.Jm**3)
    return (-0.5j) * _sym(ops.Jz, Jdif3)


def O_4_4_c(ops: JSet) -> Qobj:
    return 0.5 * ((ops.Jp**4) + (ops.Jm**4))


def O_4_4_s(ops: JSet) -> Qobj:
    return (-0.5j) * ((ops.Jp**4) - (ops.Jm**4))


# k = 5
def O_5_0(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    return (
        63 * (ops.Jz**5)
        - (70 * X - 105 * ops.I) * (ops.Jz**3)
        + (15 * (X**2) - 50 * X + 12 * ops.I) * ops.Jz
    )


def O_5_1_c(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 21 * (ops.Jz**4) - 14 * X * (ops.Jz**2) + (X**2 - X + 1.5 * ops.I)
    Jsum1 = ops.Jp + ops.Jm
    return 0.5 * _sym(T, Jsum1)


def O_5_1_s(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 21 * (ops.Jz**4) - 14 * X * (ops.Jz**2) + (X**2 - X + 1.5 * ops.I)
    Jdif1 = ops.Jp - ops.Jm
    return (-0.5j) * _sym(T, Jdif1)


def O_5_2_c(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 3 * (ops.Jz**3) - (X + 6 * ops.I) * ops.Jz
    Jsum2 = (ops.Jp**2) + (ops.Jm**2)
    return 0.5 * _sym(T, Jsum2)


def O_5_2_s(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 3 * (ops.Jz**3) - (X + 6 * ops.I) * ops.Jz
    Jdif2 = (ops.Jp**2) - (ops.Jm**2)
    return (-0.5j) * _sym(T, Jdif2)


def O_5_3_c(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 9 * (ops.Jz**2) - X - 16.5 * ops.I
    Jsum3 = (ops.Jp**3) + (ops.Jm**3)
    return 0.5 * _sym(T, Jsum3)


def O_5_3_s(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 9 * (ops.Jz**2) - X - 16.5 * ops.I
    Jdif3 = (ops.Jp**3) - (ops.Jm**3)
    return (-0.5j) * _sym(T, Jdif3)


def O_5_4_c(ops: JSet) -> Qobj:
    Jsum4 = (ops.Jp**4) + (ops.Jm**4)
    return 0.5 * _sym(ops.Jz, Jsum4)


def O_5_4_s(ops: JSet) -> Qobj:
    Jdif4 = (ops.Jp**4) - (ops.Jm**4)
    return (-0.5j) * _sym(ops.Jz, Jdif4)


def O_5_5_c(ops: JSet) -> Qobj:
    return 0.5 * ((ops.Jp**5) + (ops.Jm**5))


def O_5_5_s(ops: JSet) -> Qobj:
    return (-0.5j) * ((ops.Jp**5) - (ops.Jm**5))


# k = 6
def O_6_0(ops: JSet) -> Qobj:
    Jz2 = ops.Jz**2
    Jz4 = Jz2**2
    Jz6 = Jz4 * Jz2
    JJ = ops.JJ
    return (
        231 * Jz6
        - (315 * JJ - 735) * Jz4
        + (105 * JJ**2 - 525 * JJ + 294) * Jz2
        - 5 * JJ**3 * ops.I
        + 40 * JJ**2 * ops.I
        - 60 * JJ * ops.I
    )


def O_6_3_c(ops: JSet) -> Qobj:
    const = 3 * ops.JJ + 59.0
    T = 11 * (ops.Jz**3) - const * ops.Jz
    Jsum3 = (ops.Jp**3) + (ops.Jm**3)
    return 0.5 * _sym(T, Jsum3)


def O_6_3_s(ops: JSet) -> Qobj:
    const = 3 * ops.JJ + 59.0
    T = 11 * (ops.Jz**3) - const * ops.Jz
    Jdif3 = (ops.Jp**3) - (ops.Jm**3)
    return (-0.5j) * _sym(T, Jdif3)


def O_6_1_c(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = (
        33 * (ops.Jz**5)
        - (30 * X - 15 * ops.I) * (ops.Jz**3)
        + (5 * (X**2) - 10 * X + 12 * ops.I) * ops.Jz
    )
    Jsum1 = ops.Jp + ops.Jm
    return 0.5 * _sym(T, Jsum1)


def O_6_1_s(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = (
        33 * (ops.Jz**5)
        - (30 * X - 15 * ops.I) * (ops.Jz**3)
        + (5 * (X**2) - 10 * X + 12 * ops.I) * ops.Jz
    )
    Jdif1 = ops.Jp - ops.Jm
    return (-0.5j) * _sym(T, Jdif1)


def O_6_2_c(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = (
        33 * (ops.Jz**4)
        - (18 * X + 123 * ops.I) * (ops.Jz**2)
        + ((X**2) + 10 * X + 102 * ops.I)
    )
    Jsum2 = (ops.Jp**2) + (ops.Jm**2)
    return 0.5 * _sym(T, Jsum2)


def O_6_2_s(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = (
        33 * (ops.Jz**4)
        - (18 * X + 123 * ops.I) * (ops.Jz**2)
        + ((X**2) + 10 * X + 102 * ops.I)
    )
    Jdif2 = (ops.Jp**2) - (ops.Jm**2)
    return (-0.5j) * _sym(T, Jdif2)


def O_6_4_c(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 11 * (ops.Jz**2) - X - 38 * ops.I
    Jsum4 = (ops.Jp**4) + (ops.Jm**4)
    return 0.5 * _sym(T, Jsum4)


def O_6_4_s(ops: JSet) -> Qobj:
    X = ops.JJ * ops.I
    T = 11 * (ops.Jz**2) - X - 38 * ops.I
    Jdif4 = (ops.Jp**4) - (ops.Jm**4)
    return (-0.5j) * _sym(T, Jdif4)


def O_6_5_c(ops: JSet) -> Qobj:
    Jsum5 = (ops.Jp**5) + (ops.Jm**5)
    return 0.5 * _sym(ops.Jz, Jsum5)


def O_6_5_s(ops: JSet) -> Qobj:
    Jdif5 = (ops.Jp**5) - (ops.Jm**5)
    return (-0.5j) * _sym(ops.Jz, Jdif5)


def O_6_6_c(ops: JSet) -> Qobj:
    return 0.5 * ((ops.Jp**6) + (ops.Jm**6))


def O_6_6_s(ops: JSet) -> Qobj:
    return (-0.5j) * ((ops.Jp**6) - (ops.Jm**6))


def stevens_operator(
    ops: JSet, k: int, q: int, component: str | None = None
) -> Qobj:
    """Return a cached real tesseral Stevens operator."""
    if not isinstance(ops, JSet):
        raise TypeError("ops must be a JSet.")
    for name, value in (("k", k), ("q", q)):
        if not isinstance(value, Integral) or isinstance(value, bool):
            raise TypeError(f"{name} must be an integer.")
    if not 0 <= k <= 6:
        raise ValueError("k must be between 0 and 6.")
    if not 0 <= q <= k:
        raise ValueError("q must be between 0 and k.")
    if q == 0 and component is not None:
        raise ValueError("component must be None for q = 0.")
    if q > 0 and component not in ("c", "s"):
        raise ValueError("component must be 'c' or 's' for q > 0.")
    key = ("O", k, q, component)
    if key in ops._cache:
        return ops._cache[key]

    if q == 0:
        fn = {
            0: O_0_0,
            1: O_1_0,
            2: O_2_0,
            3: O_3_0,
            4: O_4_0,
            5: O_5_0,
            6: O_6_0,
        }.get(k)
        if fn is None:
            raise NotImplementedError(f"O_{k}_0 not implemented")
    else:
        if component not in ("c", "s"):
            raise ValueError(f"O_{k}_{q}: component must be 'c' or 's'")
        table = {
            (1, 1, "c"): O_1_1_c,
            (1, 1, "s"): O_1_1_s,
            (2, 1, "c"): O_2_1_c,
            (2, 1, "s"): O_2_1_s,
            (2, 2, "c"): O_2_2_c,
            (2, 2, "s"): O_2_2_s,
            (3, 1, "c"): O_3_1_c,
            (3, 1, "s"): O_3_1_s,
            (3, 2, "c"): O_3_2_c,
            (3, 2, "s"): O_3_2_s,
            (3, 3, "c"): O_3_3_c,
            (3, 3, "s"): O_3_3_s,
            (4, 1, "c"): O_4_1_c,
            (4, 1, "s"): O_4_1_s,
            (4, 2, "c"): O_4_2_c,
            (4, 2, "s"): O_4_2_s,
            (4, 3, "c"): O_4_3_c,
            (4, 3, "s"): O_4_3_s,
            (4, 4, "c"): O_4_4_c,
            (4, 4, "s"): O_4_4_s,
            (5, 1, "c"): O_5_1_c,
            (5, 1, "s"): O_5_1_s,
            (5, 2, "c"): O_5_2_c,
            (5, 2, "s"): O_5_2_s,
            (5, 3, "c"): O_5_3_c,
            (5, 3, "s"): O_5_3_s,
            (5, 4, "c"): O_5_4_c,
            (5, 4, "s"): O_5_4_s,
            (5, 5, "c"): O_5_5_c,
            (5, 5, "s"): O_5_5_s,
            (6, 1, "c"): O_6_1_c,
            (6, 1, "s"): O_6_1_s,
            (6, 2, "c"): O_6_2_c,
            (6, 2, "s"): O_6_2_s,
            (6, 3, "c"): O_6_3_c,
            (6, 3, "s"): O_6_3_s,
            (6, 4, "c"): O_6_4_c,
            (6, 4, "s"): O_6_4_s,
            (6, 5, "c"): O_6_5_c,
            (6, 5, "s"): O_6_5_s,
            (6, 6, "c"): O_6_6_c,
            (6, 6, "s"): O_6_6_s,
        }
        fn = table.get((k, q, component))
        if fn is None:
            raise NotImplementedError(f"O_{k}_{q}_{component} not implemented")

    # QuTiP scalar Qobjs do not support powers; the spin-zero
    # representation has only the rank-zero identity and zero multipoles.
    op = 0 * ops.I if ops.J == 0 and k > 0 else fn(ops)
    ops._cache[key] = op
    return op
