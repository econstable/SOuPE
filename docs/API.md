# Public API

The supported top-level names are `operator`, `operators`, `OperatorSet`,
`rotation_operator`, `real_wybourne_to_stevens` and `convert_energy`.

## Operators

| Call | Purpose |
| --- | --- |
| `operator(J,k,q,component=None, *, dtype=None)` | Construct one independent QuTiP Qobj. |
| `operators(J, *, dtype=None)` | Construct a reusable OperatorSet. |
| `ops.operator(k,q,component=None, *, dtype=None)` | Obtain an independent operator, reusing calculations. |
| `OperatorSet.from_qutip(J,Jx,Jy,Jz, *, dtype=None)` | Use copies of an existing spin-J representation. |

J is a finite nonnegative integer or half-integer. Integer k and q satisfy
0<=k<=6 and 0<=q<=k. For q=0 use component=None; otherwise specify `c` or
`s`. Boolean quantum numbers are rejected. Invalid types raise TypeError;
invalid values raise ValueError. Ranks above 2J vanish mathematically, subject
to floating-point cancellation. No ion-dependent factors are included.

The basis is descending m=J,...,-J and hbar=1. `dtype` may specify a QuTiP
storage type such as CSR, Dense or Dia. The public factory returns copies;
changing returned operators cannot corrupt its internal cache. `ops.Jx`,
`Jy`, `Jz`, `Jp`, `Jm` and `identity` are also independent Qobjs.

`from_qutip` checks dimensions, shape and Hermiticity. Supplying the correct
spin algebra is the caller's responsibility. Construct single-manifold
operators before embedding with `qutip.tensor`.

## Euler rotations

`rotation_operator(J,alpha,beta,gamma, *, degrees=False)` returns the ordinary
QuTiP unitary

```text
U = exp(-i alpha Jz) exp(-i beta Jy) exp(-i gamma Jz)
```

This is intrinsic ZYZ; all three generators are the original, unrotated
operators. The corresponding spatial matrix is Rz(alpha) Ry(beta) Rz(gamma).
Angles default to radians; degrees=True applies to all three angles.

`ops.rotated(alpha,beta,gamma, *, degrees=False)` returns a new operator set
using U J_i U.dag(). Its Stevens operators are therefore U O U.dag().
States and coefficients are untouched. Each call rotates the current set;
it is not an absolute orientation assignment. Dimensions and storage are
preserved. Use U ket to rotate a state, and U.dag() for the inverse:

```python
U = rotation_operator(6, 30, 40, 50, degrees=True)
rotated_state = U * state
rotated_operator = U * O * U.dag()
inverse_operator = U.dag() * O * U
```

Basic validation rejects nonfinite or non-real angles, boolean angles and a
non-boolean degrees flag. The operator factory already validates J.

**Gimbal lock:** at beta=0 or pi (modulo equivalent rotations), the first and
third Euler angles cannot be determined independently. The rotation and its
unitary remain valid. SOuPE accepts supplied angles without extracting them
from matrices, so it does not emit a gimbal-lock warning.


## Conversions

`real_wybourne_to_stevens(coefficients, *, theta)` returns B=lambda*theta*L
with unchanged `(k,q,component)` keys, units and frame. Supported ranks are
2,4,6. Values and factors must be finite real numbers. Every supplied rank
requires an explicit factor; zero factors are permitted. The exact real
combinations and source relationships are in [CONVERSIONS.md](CONVERSIONS.md).

`convert_energy(value, *, from_unit,to_unit)` converts scalars or arrays
between `meV`, `cm^-1` and `K`. K denotes E/k_B. It returns a float for scalar
input or a NumPy array with unchanged shape; finite real input is required.
It does not alter a coefficient convention.

## LaTeX

`soupe.latex.render_latex(selection="ALL")` returns a complete document;
`write_latex(path,selection="ALL")` writes it. The CLI `soupe latex --ops`
uses the same selection syntax, such as `[O_2_0,O_4_3_c]`. Compact forms
such as O20 are accepted; omitting c/s for nonzero q selects both. Canonical
identifiers are rank/order/component; the mathematical label is O_k^q.
Export follows the implemented syntax and is not an independent reference.

## Example inputs

Example configs and shared helpers are documented in examples/README.md
