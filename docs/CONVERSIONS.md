# Coefficient conventions and energy units

These are two separate operations: coefficient conversion changes the
operator basis normalisation; energy conversion changes the unit of a value.

## Real-Wybourne to Stevens

Amelin et al., [Phys. Rev. B 102, 134428 (2020)](https://doi.org/10.1103/PhysRevB.102.134428),
Appendix B, Eqs. (B1-B5), defines the operators and coefficient relation:

```text
B[k,q,component] = lambda[k,q] * theta[k] * L[k,q,component]
H = sum(B[k,q,component] * O[k,q,component])
```

Here L names Wybourne input within the API; publications use different labels.
Alexanian et al., [Phys. Rev. B 107, 224404 (2023)](https://doi.org/10.1103/PhysRevB.107.224404),
supplement Eq. (1), Tables I-III, applies theta and lambda explicitly
to the supplied crystal-field parameters. Radial expectation values are already
included in those parameters.

Use the Condon-Shortley phase convention. For positive order p,

```text
W[k,0] = C[k,0]
W[k,p,c] = C[k,-p] + (-1)**p C[k,p]
W[k,p,s] = i * (C[k,-p] - (-1)**p C[k,p])
P_J W P_J = theta[k] * lambda[k,p] * O[k,p,component]
```

There is no 1/sqrt(2) in these combinations. This is the signed-order
construction in the supplementary material of Alexanian et al. Eq. (16), 
expressed with explicit c/s keys instead of signed orders. Sine phases 
must be checked before importing coefficients from another convention.

| k | q=0 | q=1 | q=2 | q=3 | q=4 | q=5 | q=6 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 1/2 | sqrt(6) | sqrt(6)/2 | | | | |
| 4 | 1/8 | sqrt(5)/2 | sqrt(10)/4 | sqrt(35)/2 | sqrt(70)/8 | | |
| 6 | 1/16 | sqrt(42)/8 | sqrt(105)/16 | sqrt(105)/8 | 3 sqrt(14)/16 | 3 sqrt(77)/8 | sqrt(231)/16 |

The Tables II-III of Amelin et al (2020). and Table III of the 
Alexanian et al. (2023) supplementary material tabulate the pyrochlore factors. 
Supporting [thesis](https://hal.science/tel-03668858v1)
Table III.3 also gives all rank-2 entries and the pyrochlore entries
(4,0), (4,3), (6,0), (6,3), (6,6). Table III of the Alexanian et al.
supplementary material supplies the six crystal-field entries. 
The table has been extended in SOuPE including additional orders to 
support lower-symmetry crystal fields.The additional normalisation factors 
follow those tabulated in the [McPhase manual](https://mcphase.github.io/webpage/manual/node157.html),
using the real cosine/sine combinations defined above.

The converter accepts 27 components across even ranks 2,4,6. Supply theta
explicitly for each input rank; no ion is inferred. Use theta=1 only when the
input already contains that factor. Zero theta is permitted, so a general
inverse is not provided. Units and coordinate frame are preserved. Complex
spherical coefficients and alternative orthonormal real combinations are
not interchangeable inputs.

## Energy units

```python
from soupe import convert_energy

wavenumber = convert_energy(1.67, from_unit="meV", to_unit="cm^-1")
kelvin = convert_energy([0, 1.67], from_unit="meV", to_unit="K")
```

`meV` means millielectronvolt, `cm^-1` spectroscopic wavenumber, and `K`
energy expressed as E/k_B.Conversion follows E=hc*wavenumber and E=k_B*T. 
Exact SI constants are h=6.62607015e-34 J s, c=299792458 m/s, e=1.602176634e-19 C and
k_B=1.380649e-23 J/K. Thus 1 meV is approximately 8.065543937 cm^-1 or
11.604518122 K. Arrays retain their shape, negative coefficients retain
sign, and inputs are not modified. Unknown units, complex data and nonfinite
values are rejected. No automatic unit inference is performed.
