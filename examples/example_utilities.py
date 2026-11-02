"""Small shared calculations for the examples.

The Hamiltonian is built in laboratory coordinates. Reported state components
use each site's local descending-m basis. All internal energies are meV.
"""

# Read transparent example inputs without a configuration library.
import json

# Accept command-line config paths.
from pathlib import Path

# Validate fields and resolve degenerate reporting bases.
import numpy as np

# Report local angular-momentum expectation values.
from qutip import expect

# Obtain U for reporting local states.
from soupe import (
    convert_energy,
    operators,
    real_wybourne_to_stevens,
    rotation_operator,
)

# CODATA 2022 Bohr magneton, converted from J/T to meV/T.
MU_B_MEV_T = 9.2740100657e-24 / 1.602176634e-22
# Numerical degeneracy threshold in meV, far below published precision.
DEGENERACY_TOL = 1e-9


def finite_number(value, label):
    """Reject ambiguous or nonfinite values with an input-specific message."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a real number.")
    # JSON may accept NaN, but physics inputs cannot.
    if not np.isfinite(value):
        raise ValueError(f"{label} must be finite.")
    # Use the same numerical type throughout the examples.
    return float(value)


def load_config(path):
    """Read the example schema, retaining ordinary dictionaries."""
    # Keep the JSON visible.
    config = json.loads(Path(path).read_text())
    try:
        # SOuPE validates the angular momentum.
        base = operators(config["J"])
        # Required Zeeman factor.
        finite_number(config["g_factor"], "g_factor")
        # Source descriptions are explicit inputs, not inferred metadata.
        for key in ("citation", "url", "g_factor"):
            config["source"][key]
        # A label only; numerical factors are always explicit.
        config["ion"]
        # One explicit example convention.
        if config["euler_unit"] != "degrees":
            raise ValueError("euler_unit must be 'degrees'.")
        # There must be at least one oriented ion.
        if not config["sites"]:
            raise ValueError("sites must not be empty.")
        # Avoid ambiguous site columns in exported tables.
        labels = set()
        for site in config["sites"]:
            # Keep labels suitable for TSV headers.
            label = str(site["label"])
            if label in labels or not label.replace("_", "").isalnum():
                raise ValueError(
                    "Site labels must be unique alphanumeric names."
                )
            # Preserve the configured site ordering.
            labels.add(label)
            # Intrinsic ZYZ, in degrees.
            angles = site["euler_angles"]
            if len(angles) != 3:
                raise ValueError("Each site requires three Euler angles.")
            # Validate angles.
            rotation_operator(base.J, *angles, degrees=True)
        # Validate terms and apply explicit units.
        coefficients_mev(config)
        # Reject an invalid laboratory direction.
        field_direction(config)
        finite_number(config["static_field"]["magnitude_T"], "magnitude_T")
        # Check optional sweep settings when supplied.
        sweep_fields(config)
    except (KeyError, TypeError) as error:
        raise ValueError(f"Invalid example configuration: {error}") from error
    # No package-level config classes or inferred ion factors.
    return config


def coefficients_mev(config):
    """Return supplied and final Stevens mappings, with unambiguous keys."""
    # Units/convention are explicit inputs.
    crystal = config["crystal_field"]
    # Translate JSON rows into the SOuPE key convention.
    supplied = {}
    # Validate operators through the public API.
    base = operators(config["J"])
    for term in crystal["coefficients"]:
        key = (term["rank"], term["order"], term["component"])
        # Reject invalid rank, order or c/s combinations.
        base.operator(*key)
        if key in supplied:
            raise ValueError(f"Duplicate crystal-field term {key}.")
        supplied[key] = finite_number(term["value"], f"coefficient {key}")
    # An empty model would not demonstrate crystal fields.
    if not supplied:
        raise ValueError("At least one crystal-field coefficient is required.")
    # Unit scaling precedes convention conversion; the two steps are separate.
    mev = {
        key: convert_energy(value, from_unit=crystal["unit"], to_unit="meV")
        for key, value in supplied.items()
    }
    # Direct coefficients multiply O.
    if crystal["convention"] == "stevens":
        return supplied, mev
    # Do not guess the input convention.
    if crystal["convention"] != "wybourne":
        raise ValueError("convention must be 'stevens' or 'wybourne'.")
    theta = {int(k): value for k, value in crystal["theta"].items()}
    # Apply lambda and the explicitly supplied ion factors exactly once.
    return supplied, real_wybourne_to_stevens(mev, theta=theta)


def field_direction(config):
    """Normalise the global field direction, leaving its magnitude separate."""
    direction = np.asarray(config["static_field"]["direction"], dtype=float)
    if direction.shape != (3,) or not np.all(np.isfinite(direction)):
        raise ValueError("Field direction must be a finite three-vector.")
    # [111] becomes [1,1,1]/sqrt(3).
    norm = np.linalg.norm(direction)
    if norm == 0:
        raise ValueError(
            "Field direction must be nonzero, even at zero field."
        )
    # Changing direction does not change field units.
    return direction / norm


def sweep_fields(config):
    """Use the configured sweep or a single static-field value if omitted."""
    # A static-only config remains useful.
    sweep = config.get("field_sweep")
    if sweep is None:
        return np.array(
            [
                finite_number(
                    config["static_field"]["magnitude_T"], "magnitude_T"
                )
            ]
        )
    # Number of inclusive, equally spaced samples.
    points = sweep["points"]
    if isinstance(points, bool) or not isinstance(points, int) or points < 2:
        raise ValueError("Sweep points must be an integer of at least two.")
    # Signed tesla allowed.
    start = finite_number(sweep["start_T"], "start_T")
    # Descending sweeps allowed.
    stop = finite_number(sweep["stop_T"], "stop_T")
    if start == stop:
        raise ValueError("Sweep start_T and stop_T must differ.")
    # Validate energy units without guessing.
    unit = sweep["output_unit"]
    convert_energy(0, from_unit="meV", to_unit=unit)
    # Ground state is separately excluded.
    levels = sweep["excited_levels"]
    if (
        isinstance(levels, bool)
        or not isinstance(levels, int)
        or not 1 <= levels <= int(2 * config["J"])
    ):
        raise ValueError("excited_levels must lie between 1 and 2J.")
    # Includes both field endpoints.
    return np.linspace(start, stop, points)


def build_site_hamiltonians(config):
    """Construct each site's CEF and a common laboratory Zeeman term per T."""
    # Final Stevens energies, meV.
    _, coefficients = coefficients_mev(config)
    # Unrotated laboratory spin generators.
    base = operators(config["J"])
    # Same field direction at every site.
    direction = field_direction(config)
    j_field = sum(
        (v * op for v, op in zip(direction, (base.Jx, base.Jy, base.Jz))),
        0 * base.identity,
    )
    # Published sign.
    zeeman_per_t = -config["g_factor"] * MU_B_MEV_T * j_field
    # Independent ions, not a tensor-product interacting model.
    sites = []
    for site in config["sites"]:
        # Keep publication axes/coefficients paired.
        angles = site["euler_angles"]
        # U J_i U.dag().
        oriented = base.rotated(*angles, degrees=True)
        h_cf = sum(
            (
                value * oriented.operator(*key)
                for key, value in coefficients.items()
            ),
            0 * base.identity,
        )
        unitary = rotation_operator(config["J"], *angles, degrees=True)
        # U maps local to lab.
        sites.append((str(site["label"]), h_cf, unitary))
    # Show the separate physical terms.
    return base, sites, zeeman_per_t


def energy_groups(energies):
    """Group numerically degenerate levels without altering their energies."""
    # Energy order is maintained for the spectrum table.
    groups = []
    for index, energy in enumerate(energies):
        if (
            not groups
            or abs(energy - energies[groups[-1][0]]) > DEGENERACY_TOL
        ):
            # Start a new physical energy level.
            groups.append([index])
        else:
            # Same level within numerical precision.
            groups[-1].append(index)
    # Compare to the first member to avoid chained grouping.
    return groups


def reporting_states(energies, states, unitary, local_jz, *, zero_field):
    """Report local components; resolve zero-field degeneracies with P Jz P.

    This changes only the basis within a degenerate subspace. No artificial
    field or perturbation is added, and no nondegenerate eigenstates are mixed.
    """
    # Local |J,m> basis.
    local = [unitary.dag() * state for state in states]
    if zero_field:
        # Choose a basis only if needed.
        for group in energy_groups(energies):
            if len(group) == 1:
                # Leave isolated energy eigenstates unchanged.
                continue
            columns = np.column_stack([local[i].full().ravel() for i in group])
            # P Jz P.
            projected = columns.conj().T @ local_jz.full() @ columns
            # Resolve projected moments.
            _, mixing = np.linalg.eigh(projected)
            # Preserve the whole subspace.
            original = [local[i] for i in group]
            for column, index in enumerate(group):
                local[index] = sum(
                    (
                        mixing[row, column] * state
                        for row, state in enumerate(original)
                    ),
                    0 * original[0],
                ).unit()
    # Fix only an arbitrary global phase.
    for index, state in enumerate(local):
        # Descending m=J,...,-J order.
        amplitudes = state.full().ravel()
        largest = amplitudes[np.argmax(np.abs(amplitudes))]
        local[index] = state * np.exp(-1j * np.angle(largest))
    return local


def state_summary(state, J):
    """Display the five largest local amplitudes and local <Jz>."""
    # Basis indices correspond to m=J-index.
    amplitudes = state.full().ravel()
    # Keep output readable.
    dominant = np.argsort(-np.abs(amplitudes))[:5]
    # Complex amplitudes, rather than probabilities.
    terms = []
    for index in dominant:
        value = amplitudes[index]
        # Suppress entries below display precision.
        if abs(value) >= 1e-3:
            terms.append(
                f"({value.real:+.4f}{value.imag:+.4f}i)|{J - index:g}>"
            )
    return " + ".join(terms), float(np.real(expect(operators(J).Jz, state)))
