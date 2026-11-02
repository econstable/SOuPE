# Stevens operator conventions and scientific origins

The polynomial definitions follow the
[Stevens operator table](Stevens_Operators-Ballou.pdf) derived by R. Ballou.
The package constructs its 49 real tesseral components, ranks 0-6, as QuTiP
matrices. Rank k is the subscript and order q is the superscript: O_k^q.
Identifiers and API arguments use `(k, q, component)`; c/s distinguish cosine
and sine for q>0. Signed negative orders in the publications correspond to
sine components here. The identity O_0^0=I is included.

The basis is |J,m> in descending m order, hbar=1, J_+=J_x+iJ_y and
J_-=J_x-iJ_y. X=J(J+1)I. Ion-dependent Stevens factors and radial expectation
values are not part of the polynomials. For example,

```text
O_2^0 = 3 J_z^2 - X
O_2^1(c) = (J_z J_x + J_x J_z)/2
O_2^2(s) = J_x J_y + J_y J_x
```

These definitions also appear in Amelin et al.,
[Phys. Rev. B 102, 134428 (2020)](https://doi.org/10.1103/PhysRevB.102.134428),
Appendix A. Supporting derivations appear in
[Alexanian’s thesis](https://hal.science/tel-03668858v1), Eqs. (III.20-21). The [McPhase operator
 table](https://mcphase.github.io/webpage/manual/node159.html) was used as a secondary
correctness check and is not the origin of the package.

## Notes on Stevens_Operators-Ballou.pdf

Stevens_Operators-Ballou.pdf sometimes labels axial components c; 
the API uses None because there is only one q=0 component. This is notation only.

## Euler rotations

Amelin et al. (2020), Eqs. (5-7), expresses local crystal-field operators in the
common cubic coordinate system. Alexanian et al.,
[Phys. Rev. B 107, 224404 (2023)](https://doi.org/10.1103/PhysRevB.107.224404),
supplementary Table VII, provides local axes. SOuPE uses intrinsic ZYZ:
U=exp(-i alpha Jz) exp(-i beta Jy) exp(-i gamma Jz), with unrotated generators.
Operators transform as U O U.dag(); the inverse uses U.dag() O U.
Example Euler angles are derived from the published axes; the example guide
shows their relationship explicitly. Supporting rotation mathematics appears
in Alexanian's thesis, section III.2.1.2.

At beta=0 or pi, Euler angles are non-unique (gimbal lock), but U is valid.
This limitation should be noted when reporting or fitting orientations;
see API.md. A future update could accept local axes
or rotation matrices directly, avoiding the need to extract Euler angles
from a specified orientation.

Changing the signs of both order-three crystal-field coefficients together
with reversing local x and y describes the same orientation. The supplementary 
material of Alexanian et al., p. 2, explicitly discusses its frame/sign 
change relative to Amelin et al. The examples preserve each source's pairing. 
Equal zero-field spectra alone do not establish that a field-dependent orientation 
is correct.

Returned public operators are independent copies. High-rank floating-point
roundoff is symmetrised to preserve Hermiticity. This numerical treatment
and storage selection are software choices, not additional physical models.
