"""Command-line entry point for exporting operator definitions as LaTeX."""

from __future__ import annotations

import argparse
from pathlib import Path

from .latex import write_latex


def main(argv: list[str] | None = None) -> int:
    """Run the standalone LaTeX export command."""
    parser = argparse.ArgumentParser(prog="soupe")
    subparsers = parser.add_subparsers(dest="command")

    latex_parser = subparsers.add_parser(
        "latex", help="write Stevens operator definitions as LaTeX"
    )
    latex_parser.add_argument(
        "--ops",
        default="ALL",
        help="operator list, e.g. ALL or '[O_2_0,O_4_3_c]'",
    )
    latex_parser.add_argument(
        "-o",
        "--output",
        default="stevens_operators.tex",
        type=Path,
        help="output .tex file",
    )

    args = parser.parse_args(argv)
    if args.command == "latex":
        out = write_latex(args.output, args.ops)
        print(out)
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
