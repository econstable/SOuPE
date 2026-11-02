"""Export independent-site CEF + Zeeman transition energies as a TSV table."""

# Accept config/output paths without adding a package command.
import argparse

# Write a tab-separated numerical result.
from pathlib import Path

# Tabulate all site/level gaps at each field.
import numpy as np

# Reuse the same inputs and Hamiltonian logic.
from example_utilities import (
    build_site_hamiltonians,
    load_config,
    sweep_fields,
)

# Keep output units separate from convention.
from soupe import convert_energy


def sweep_spectrum(config):
    """Return fields, site labels and gaps with shape (B,site,level)."""
    # Inclusive configured field grid, in tesla.
    fields = sweep_fields(config)
    # Static fallback if no sweep given.
    sweep = config.get("field_sweep") or {}
    count = sweep.get("excited_levels", min(4, int(2 * config["J"])))
    # Default community energy unit.
    unit = sweep.get("output_unit", "cm^-1")
    # Build CEFs once.
    _, sites, zeeman = build_site_hamiltonians(config)
    # Every site independently.
    values = np.empty((len(fields), len(sites), count))
    # Sweep magnitude along fixed direction.
    stop = count + 1  # Slice excludes the ground state and includes E_count.
    for row, field in enumerate(fields):
        for column, (_, h_cf, _) in enumerate(sites):
            # Sorted each field.
            energy = (h_cf + field * zeeman).eigenenergies()
            # E_n-E_0.
            values[row, column] = energy[1:stop] - energy[0]
    gaps = convert_energy(values, from_unit="meV", to_unit=unit)
    return fields, gaps, [label for label, _, _ in sites], unit


def write_sweep(config, path):
    """Write a TSV with units and site/level identities."""
    fields, gaps, labels, unit = sweep_spectrum(config)
    columns = ["B_T"] + [
        f"site{site}_dE{level}_{unit}"
        for site in labels
        for level in range(1, gaps.shape[2] + 1)
    ]
    table = np.column_stack((fields, gaps.reshape(len(fields), -1)))
    # The header describes sorted levels, without eigenstate tracking.
    header = (
        f"{config['source']['citation']}\n"
        f"direction={config['static_field']['direction']}; "
        f"g={config['g_factor']}; energy unit={unit}\n"
        f"g source: {config['source']['g_factor']}\n"
        "E_n-E_0; levels sorted at every field.\n" + "\t".join(columns)
    )
    np.savetxt(path, table, delimiter="\t", header=header, fmt="%.12g")
    # A short confirmation of output dimensions.
    return len(fields), len(labels)


def main():
    """Run the field sweep and save the calculated energies."""
    # Script-only interface.
    parser = argparse.ArgumentParser(description=__doc__)
    # Same JSON used by static_spectrum.
    parser.add_argument("config", type=Path)
    # Explicit TSV path.
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        fields, sites = write_sweep(load_config(args.config), args.output)
    except (ValueError, OSError) as error:
        # Reject invalid inputs without a long traceback.
        parser.error(str(error))
    print(
        f"Wrote {fields} fields for {sites} independent sites to {args.output}"
    )


if __name__ == "__main__":
    main()
