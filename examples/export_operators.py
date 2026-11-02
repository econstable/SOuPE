"""Export the Stevens expressions selected by a JSON example config."""

# Read config and output paths without changing the package CLI.
import argparse

# Keep file destinations explicit.
from pathlib import Path

# Reuse the small example validation.
from example_utilities import (
    load_config,
)

# Existing exporter accepts operator names.
from soupe.latex import (
    write_latex,
)


def main():
    """Translate config terms into canonical exporter identifiers."""
    # Standalone wrapper.
    parser = argparse.ArgumentParser(description=__doc__)
    # Select either publication config.
    parser.add_argument("config", type=Path)
    # Complete TeX file.
    parser.add_argument("--output", type=Path, required=True)
    # No new JSON-aware behaviour in soupe.latex.
    args = parser.parse_args()
    try:
        # Validated rank/order/component rows.
        config = load_config(args.config)
        # Include configured zero-valued terms as well.
        names = []
        for term in config["crystal_field"]["coefficients"]:
            # Rank before order.
            name = f"O_{term['rank']}_{term['order']}"
            if term["component"] is not None:
                # Explicit cosine/sine selection.
                name += f"_{term['component']}"
            # Export unrotated definitions, not coefficients.
            names.append(name)
        # Delegate all rendering.
        write_latex(args.output, ",".join(names))
    except (ValueError, OSError) as error:
        # Show actionable configuration/file errors.
        parser.error(str(error))
    print(f"Wrote {len(names)} Stevens operator expressions to {args.output}")


# Imports leave the example inert.
if __name__ == "__main__":
    # Run the config-to-export demonstration.
    main()
