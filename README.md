# SOuPE

A **Stevens Operator Utility Package for Everyone** that constructs Stevens 
operators as ordinary QuTiP objects, making rare-earth crystal-field 
Hamiltonians easy to write without coding the operator polynomials yourself. 
Its scientific basis includes Amelin et al., 
[Phys. Rev. B 102, 134428 (2020)](https://doi.org/10.1103/PhysRevB.102.134428)
and Alexanian et al., [Phys. Rev. B 107, 224404 (2023)](https://doi.org/10.1103/PhysRevB.107.224404).
[Alexanian's thesis](https://hal.science/tel-03668858v1) provides supporting
mathematical detail. The tools are general-purpose building blocks for other
ions and systems, including QuTiP projects with time-dependent Hamiltonians.


## Construct an operator

```python
import soupe

O20 = soupe.operator(J=6, k=2, q=0)
O43 = soupe.operator(J=6, k=4, q=3, component="c")
ops = soupe.operators(J=6)
H = -0.26 * ops.operator(2, 0) + 0.0045 * ops.operator(4, 0)
energies = H.eigenenergies()
```

J may be any finite nonnegative integer or half-integer. All 49 real tesseral
components through rank 6 are available in descending m=J,...,-J order,
hbar=1. For nonzero q select c or s; q=0 has no component label.
Ion factors and energy units belong in coefficients, rather than matrices.

## Rotations and conversions

```python
oriented = ops.rotated(30, 40, 50, degrees=True)
U = soupe.rotation_operator(6, 30, 40, 50, degrees=True)
B = soupe.real_wybourne_to_stevens(
    {(2, 0, None): 51.4}, theta={2: -1 / 99}
)
gap_cm = soupe.convert_energy(1.67, from_unit="meV", to_unit="cm^-1")
gap_K = soupe.convert_energy(1.67, from_unit="meV", to_unit="K")
```

Intrinsic ZYZ rotations use U=exp(-i alpha Jz) exp(-i beta Jy) exp(-i gamma Jz)
and transform operators as U O U.dag(). Angles default to radians. At beta=0
or pi the first and third angles are not independent (gimbal lock), although
the rotation remains valid.

Coefficient conversion applies B=lambda*theta*L for an explicitly defined
real-Wybourne convention across ranks 2,4,6. Energy-unit conversion is
separate and accepts scalars/arrays in meV, cm^-1 and K (E/k_B).
See [API](docs/API.md), [operator conventions](docs/CONVENTIONS.md) and
[conversions](docs/CONVERSIONS.md) for precise definitions.

## Calculations from publication-based configs

```bash
uv sync
uv run python examples/static_spectrum.py examples/tb2ti2o7_Amelin2020.json
uv run python examples/static_spectrum.py examples/tb2ti2o7_Alexanian2023.json
uv run python examples/field_sweep_spectrum.py examples/tb2ti2o7_Amelin2020.json --output field_spectrum.tsv
uv run python examples/export_operators.py examples/tb2ti2o7_Alexanian2023.json --output operators.tex
```

The configs specify Stevens or Wybourne input, units, Euler angles and field
settings. The static report shows coefficients, energies and local state
components. The sweep writes independent-site transition energies as TSV.
These examples need only the package's NumPy/QuTiP dependencies.
See the [example guide](examples/README.md) for sources and input meanings.

Example configuration and Hamiltonian helpers are described in examples/README.md.

## LaTeX export

```bash
uv run soupe latex --ops "[O_2_0,O_4_3_c]" -o selected.tex
uv run soupe latex --ops ALL -o stevens_operators.tex
```

Identifiers use rank, order, component; output uses O_k^q. Expressions are
exported from code and do not provide independent scientific validation.

## Development

```bash
uv run ruff check src examples
uv run ruff format --check src examples
uv build
```

See [contributing](CONTRIBUTING.md) and [source traceability](docs/TRACEABILITY.md).
The [development and review disclosure](CONTRIBUTING.md#development-and-review)
describes the use of generative coding agents.

## Licence

SOuPE's original code, documentation and examples are licensed under the
[BSD 3-Clause License](LICENSE). Third-party material retains its own
licensing terms. Scientific references are listed in
[source traceability](docs/TRACEABILITY.md).