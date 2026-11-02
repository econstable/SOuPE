"""Calculate a static CEF + Zeeman spectrum from a publication-based config."""

# Accept a JSON file and an optional text report destination.
import argparse

# Write the same readable report shown in the terminal.
from pathlib import (
    Path,
)

# Shared example code, not a package config API.
from example_utilities import (
    build_site_hamiltonians,
    coefficients_mev,
    energy_groups,
    load_config,
    reporting_states,
    state_summary,
)

# Show calculated energies in community units.
from soupe import (
    convert_energy,
)


def static_report(config):
    """Return input details and spectra for every configured site."""
    # Separate unit/convention steps.
    supplied, stevens = coefficients_mev(config)
    # Describe what the input numbers mean.
    crystal = config["crystal_field"]
    # Field magnitude is in tesla, B=mu_0 H.
    field = config["static_field"]
    lines = [
        config["source"]["citation"],
        config["source"]["url"],
        f"Ion: {config['ion']}; J={config['J']}; g={config['g_factor']}",
        f"Field: {field['magnitude_T']:g} T along {field['direction']}",
        f"g source: {config['source']['g_factor']}",
        "",
        "======== Crystal-field coefficients ========",
        f"Input: {crystal['convention']}; units: {crystal['unit']}",
        " rank order component       input       Stevens (meV)",
    ]
    # Include every term, even a zero value.
    for (
        key,
        value,
    ) in supplied.items():
        # Same ordering as soupe.operator arguments.
        k, q, component = key
        lines.append(
            f"{k:5d} {q:5d} {component or '-':>9} "
            f"{value:12.7g} {stevens[key]:19.10g}"
        )
    # Explain the conversion used here.
    if crystal["convention"] == "wybourne":
        lines.append(
            "B = lambda * theta * L after converting input units to meV."
        )
        lines.append(f"Explicit theta: {crystal['theta']}")
    else:
        lines.append("Stevens input multiplies O directly after unit scaling.")
    # QuTiP Hamiltonians.
    base, sites, zeeman = build_site_hamiltonians(config)
    # Independently oriented local ions.
    for label, h_cf, unitary in sites:
        # Both terms in meV.
        hamiltonian = h_cf + field["magnitude_T"] * zeeman
        # Sorted energy eigenstates.
        energies, states = hamiltonian.eigenstates()
        # Remove each site's arbitrary origin.
        gaps = energies - energies[0]
        cm = convert_energy(gaps, from_unit="meV", to_unit="cm^-1")
        kelvin = convert_energy(gaps, from_unit="meV", to_unit="K")
        # Label numerical degeneracies clearly.
        groups = energy_groups(energies)
        group_id = {
            i: group for group, indices in enumerate(groups) for i in indices
        }
        lines.extend(
            [
                "",
                f"======== Site {label}: all relative levels ========",
                "state group       meV        cm^-1            K",
            ]
        )
        # All 2J+1 states are reported.
        for index, gap in enumerate(gaps):
            lines.append(
                f"{index:5d} {group_id[index]:5d} {gap:11.6f} "
                f"{cm[index]:12.6f} {kelvin[index]:12.6f}"
            )
        # Do not mix split, nonzero-field states.
        zero = field["magnitude_T"] == 0
        local = reporting_states(
            energies, states, unitary, base.Jz, zero_field=zero
        )
        lines.extend(
            [
                "",
                "Lowest four states: local |J,m> components",
                (
                    "Zero-field degenerate groups: projected local Jz basis."
                    if zero
                    else "Nonzero field: local energy eigenstates."
                ),
                "Largest amplitude real positive; overall phase arbitrary.",
            ]
        )
        # Keep the state report compact.
        for index, state in enumerate(local[:4]):
            terms, moment = state_summary(state, base.J)
            lines.append(f"psi{index}: <Jz_local>={moment:+.6f}; {terms}")
    # Terminal and saved report are identical.
    return "\n".join(lines) + "\n"


def main():
    """Run the calculation without requiring an installed example command."""
    # Standalone script.
    parser = argparse.ArgumentParser(description=__doc__)
    # Either supplied publication config.
    parser.add_argument("config", type=Path)
    # Optional plain-text report.
    parser.add_argument("--output", type=Path)
    # Use the specified inputs without hidden profiles.
    args = parser.parse_args()
    try:
        # Validate then solve.
        report = static_report(load_config(args.config))
        if args.output:
            # No output-directory framework.
            args.output.write_text(report)
    except (ValueError, OSError) as error:
        # Give a concise, useful input error.
        parser.error(str(error))
    # Present only the calculated results and their inputs.
    print(report, end="")


# Importing the example does not run the calculation.
if __name__ == "__main__":
    # Direct invocation provides the CLI.
    main()
