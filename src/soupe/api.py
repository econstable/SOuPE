"""QuTiP Stevens operators and Euler rotations."""

from __future__ import annotations

from dataclasses import dataclass, field

from qutip import Qobj

from ._jset import JSet, jset, jset_from_ops
from ._stevens import stevens_operator
from .rotations import _unitary

__all__ = ["OperatorSet", "operator", "operators"]


@dataclass(frozen=True, eq=False)
class OperatorSet:
    """Reusable Stevens operator factory for one angular momentum manifold.

    Construct using `operators` or `OperatorSet.from_qutip`. Results and
    angular momentum properties return independent QuTiP objects; changing
    one cannot alter the factory's internal cache. Rotation returns a new
    set. The public convention is real tesseral extended Stevens throughout.

    Attributes
    ----------
    J : float
        Nonnegative integer or half-integer angular momentum.
    Jx, Jy, Jz, Jp, Jm : qutip.Qobj
        Copies of Cartesian, raising, and lowering operators, with hbar = 1.
    identity : qutip.Qobj
        Copy of the identity with matching QuTiP dimensions and backend.
    dtype : type
        Storage backend of the angular momentum operators.
    """

    _ops: JSet = field(repr=False)

    def __repr__(self) -> str:
        return f"OperatorSet(J={self.J}, dtype={self.dtype.__name__!r})"

    @property
    def J(self) -> float:
        """Angular momentum quantum number."""
        return self._ops.J

    @property
    def dtype(self) -> type:
        """QuTiP storage backend of the angular momentum operators."""
        return self._ops.Jx.dtype

    @property
    def Jx(self) -> Qobj:
        """Independent copy of the x angular momentum operator."""
        return self._ops.Jx.copy()

    @property
    def Jy(self) -> Qobj:
        """Independent copy of the y angular momentum operator."""
        return self._ops.Jy.copy()

    @property
    def Jz(self) -> Qobj:
        """Independent copy of the z angular momentum operator."""
        return self._ops.Jz.copy()

    @property
    def Jp(self) -> Qobj:
        """Independent copy of the raising operator."""
        return self._ops.Jp.copy()

    @property
    def Jm(self) -> Qobj:
        """Independent copy of the lowering operator."""
        return self._ops.Jm.copy()

    @property
    def identity(self) -> Qobj:
        """Independent identity with the set's dimensions and backend."""
        return self._ops.I.copy()

    @classmethod
    def from_qutip(
        cls, J: float, Jx: Qobj, Jy: Qobj, Jz: Qobj, *, dtype=None
    ) -> OperatorSet:
        """Construct from existing operators in a spin-J representation.

        Parameters
        ----------
        J : float
            Nonnegative integer or half-integer angular momentum.
        Jx, Jy, Jz : qutip.Qobj
            Hermitian matrices with shape (2J+1, 2J+1) and matching dims.
            Inputs are copied. Their spin algebra is the caller's contract.
        dtype : str or type, optional
            QuTiP backend. If omitted, use Jx's backend for all operators.

        Returns
        -------
        OperatorSet
            Independent set preserving the supplied QuTiP dimensions.
        """
        checked = jset_from_ops(J, Jx, Jy, Jz)
        backend = Jx.dtype if dtype is None else dtype
        copied = [op.to(backend).copy() for op in (Jx, Jy, Jz)]
        return cls(jset_from_ops(checked.J, *copied))

    def operator(
        self, k: int, q: int, component: str | None = None, *, dtype=None
    ) -> Qobj:
        """Construct a real tesseral Stevens operator.

        Parameters
        ----------
        k : int
            Rank from 0 to 6 inclusive.
        q : int
            Order from 0 to k inclusive.
        component : {None, 'c', 's'}, optional
            Required for q > 0; must be None for q = 0.
        dtype : str or type, optional
            Output backend; defaults to this set's backend.

        Returns
        -------
        qutip.Qobj
            Independent Hermitian operator. Internal caching never exposes
            a shared mutable result. No ion-dependent factors are included.
        """
        result = stevens_operator(self._ops, k, q, component)
        # High-rank polynomials after rotation can acquire roundoff-sized
        # anti-Hermitian entries. Tesseral Stevens operators are Hermitian;
        # restore that property numerically before handing them to QuTiP.
        result = 0.5 * (result + result.dag())
        return result.to(self.dtype if dtype is None else dtype).copy()

    def rotated(self, alpha, beta, gamma, *, degrees=False) -> OperatorSet:
        """Return a set transformed by the intrinsic ZYZ Euler unitary.

        Parameters
        ----------
        alpha, beta, gamma : real
            Angles in U=exp(-i alpha Jz) exp(-i beta Jy) exp(-i gamma Jz).
            All generators belong to the current set and are unrotated.
        degrees : bool, optional
            Interpret angles as degrees when True; otherwise radians.

        Returns
        -------
        OperatorSet
            New set using U J_i U.dag(), preserving dimensions and storage.
            States and coefficients are untouched. Each call transforms
            the current set, rather than assigning an absolute orientation.
        """
        unitary = _unitary(self._ops, alpha, beta, gamma, degrees=degrees)
        transformed = [
            unitary * op * unitary.dag()
            for op in (self._ops.Jx, self._ops.Jy, self._ops.Jz)
        ]
        return type(self).from_qutip(self.J, *transformed, dtype=self.dtype)


def operators(J: float, *, dtype=None) -> OperatorSet:
    """Create a reusable Stevens operator set from QuTiP angular momentum.

    Parameters
    ----------
    J : float
        Nonnegative integer or half-integer angular momentum.
    dtype : str or type, optional
        QuTiP storage backend, e.g. 'CSR', 'Dense', or 'Dia'.

    Returns
    -------
    OperatorSet
        Operator factory in QuTiP's descending m = J, ..., -J basis.
        Euler rotations do not relabel that Hilbert-space basis.
    """
    return OperatorSet(jset(J, dtype=dtype))


def operator(
    J: float,
    k: int,
    q: int,
    component: str | None = None,
    *,
    dtype=None,
) -> Qobj:
    """Construct one Stevens operator as an ordinary QuTiP Qobj.

    Parameters
    ----------
    J : float
        Nonnegative integer or half-integer angular momentum.
    k : int
        Rank from 0 to 6 inclusive.
    q : int
        Order from 0 to k inclusive.
    component : {None, 'c', 's'}, optional
        Required for q > 0; must be None for q = 0.
    dtype : str or type, optional
        QuTiP output backend, e.g. 'CSR', 'Dense', or 'Dia'.

    Returns
    -------
    qutip.Qobj
        Independent Hermitian operator, with hbar = 1 and extended Stevens
        normalisation. No radial or ion-dependent factors are included.
    """
    return operators(J, dtype=dtype).operator(k, q, component)
