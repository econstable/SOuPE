# Spectrum examples

These examples demonstrate the SOuPE operator factory, coefficient convention
conversion, energy units and Euler rotations using Tb2Ti2O7 parameters.
They calculate independent-ion crystal-field and Zeeman Hamiltonians with
QuTiP.

Install the example project environment with `uv sync`.

## Inputs

- `tb2ti2o7_Amelin2020.json`: Amelin et al.,
  [Phys. Rev. B 102, 134428 (2020)](https://doi.org/10.1103/PhysRevB.102.134428),
  Table I, Eqs. (5-6), section IV. The coefficients are Stevens energies;
  g=1.4 is the effective g-factor stated in this work.
- `tb2ti2o7_Alexanian2023.json`: Alexanian et al.,
  [Phys. Rev. B 107, 224404 (2023)](https://doi.org/10.1103/PhysRevB.107.224404),
  supplementary Eq. (1), Tables I-III and VII. The coefficients are Wybourne
  energies with explicit Tb3+ factors. The default static field is zero.

The Alexanian et al. configuration includes an effective g-factor of g=1.35 for 
illustrative nonzero-field use. That value is sourced separately from 
[Alexanian's thesis](https://hal.science/tel-03668858v1), section V.3.2.1. 
Changing the field therefore makes the extension and its chosen g explicit. 
The same crystal-field values and axes occur in the thesis; a user interested in 
its field calculation can use this setting.

Each coefficient row specifies rank, order, component and value. These map
directly to `ops.operator(rank, order, component)`. For zero order the
component is null; otherwise use c or s. `crystal_field.unit` is meV, cm^-1
or K; K means E/k_B, not the calculation temperature. The examples first
convert energies to meV, then apply lambda*theta to Wybourne inputs once.
Stevens inputs multiply operators directly. Factors are not inferred from
the ion label or J.

`static_field.direction` is a nonzero global cubic direction, normalised by
the example. `magnitude_T` is B=mu_0 H in tesla; signed values are allowed.
Edit these fields to calculate a different applied field. A sweep uses the
same direction and its own start_T, stop_T and inclusive number of points.
`excited_levels` selects E1-E0 through En-E0; output_unit sets their unit.
Omitting field_sweep requests one static-field sample.

## Static spectrum

```bash
uv run python examples/static_spectrum.py examples/tb2ti2o7_Amelin2020.json
uv run python examples/static_spectrum.py examples/tb2ti2o7_Alexanian2023.json --output static.txt
```

The report lists input and converted coefficients, all relative levels in
meV, cm^-1 and K, and dominant local |J,m> components of the lowest four
states for every configured site. The basis is descending m=J,...,-J.

At zero field, diagonalising projected local Jz within each numerically
degenerate group selects a reporting basis. These states have definite
projected Jz, not necessarily definite full Jz. A residual degeneracy of the
projected operator can still leave its basis non-unique. The largest amplitude
is displayed real positive; overall phases have no physical significance.
At nonzero field the energy eigenstates are reported in local coordinates
without mixing split levels. State amplitudes are recovered using U.dag()
from the laboratory eigenvectors.

## Field sweep

```bash
uv run python examples/field_sweep_spectrum.py examples/tb2ti2o7_Amelin2020.json --output field_spectrum.tsv
```

The TSV contains B_T, then each site's selected relative excitation energies.
Headers identify the site, level and energy unit. Levels are sorted at every
field; the script does not track individual eigenstate identities through
crossings. All four ions are calculated independently; the [111] field leaves
sites 2-4 equivalent. The file can be analysed with tools of the user's choice.

## LaTeX expressions

```bash
uv run python examples/export_operators.py examples/tb2ti2o7_Amelin2020.json --output configured_operators.tex
```

The script exports the unrotated Stevens polynomials selected by the config.

## Euler angles and published axes

The executable inputs use intrinsic ZYZ, with
U=exp(-i alpha Jz) exp(-i beta Jy) exp(-i gamma Jz).
The spatial matrix is Rz(alpha) Ry(beta) Rz(gamma). Each site's crystal-field
polynomial uses U J_i U.dag(); the common Zeeman term uses laboratory J_i.
Coefficients and their local xy orientation must remain paired.

The following angles in degrees are derived from the local axes in the
supplementary material of Alexanian et al., Table VII.
Let t=acos(1/sqrt(3))=54.735610317245 degrees.

|Site | sqrt(6) x | sqrt(2) y | sqrt(3) z | alpha | beta | gamma|
| --- | --- | --- | --- | --- | --- | --- |
| 1 | (-2,1,1) | (0,-1,1) | (1,1,1) | 45 | t | 120 |
| 2 | (-2,-1,-1) | (0,1,-1) | (1,-1,-1) | -45 | 180-t | -60 |
| 3 | (2,1,-1) | (0,-1,-1) | (-1,1,-1) | 135 | 180-t | -60 |
| 4 | (2,-1,1) | (0,1,1) | (-1,-1,1) | -135 | t | 120 |

Normalise the integer vectors by the indicated square roots. To audit the
angle derivation, form z=(cos(alpha)sin(beta),sin(alpha)sin(beta),cos(beta)),
e_theta=(cos(alpha)cos(beta),sin(alpha)cos(beta),-sin(beta)) and
e_phi=(-sin(alpha),cos(alpha),0). Then x=cos(gamma)e_theta+sin(gamma)e_phi
and y=-sin(gamma)e_theta+cos(gamma)e_phi reproduce the table.

For the Amelin et al. config, add 180 degrees to gamma: x/y reverse while z is
unchanged, paired with its order-three coefficient signs. The
supplementary material of Alexanian et al., p. 2, describes the sign/frame 
convention difference between these works. Angles outside the range [-180,180] 
are valid inputs.


## How the scripts use SOuPE

`example_utilities.py` contains the small shared config reader and calculations.
`coefficients_mev` calls `convert_energy` and `real_wybourne_to_stevens`;
`build_site_hamiltonians` calls `operators`, `rotated` and `operator`, then
adds -g*mu_B*B_direction.J per tesla. The two spectrum scripts call QuTiP's
eigensolvers. These helpers are shared by the example scripts.
