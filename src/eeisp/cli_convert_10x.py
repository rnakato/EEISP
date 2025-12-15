from __future__ import annotations

from pathlib import Path
from typing import Literal

import typer

from eeisp.convert_10x import convert_10x_to_matrix
from eeisp.version import get_version

app = typer.Typer(add_completion=False, help="Convert 10X CellRanger output to EEISP matrix format.")


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"eeisp-convert-10x {get_version()}")
        raise typer.Exit()


@app.command()
def cli(
    input_path: Path = typer.Argument(
        ...,
        help="Path to 10X data (directory for mtx format, .h5 file for h5 format)",
    ),
    output: Path = typer.Argument(..., help="Output matrix file path"),
    format: Literal["mtx", "h5"] = typer.Option(
        "mtx",
        "--format",
        "-f",
        help="Input format: 'mtx' for directory with matrix.mtx, or 'h5' for .h5 file",
    ),
    var_names: str = typer.Option(
        "gene_symbols",
        "--var-names",
        help="Variable names to use: 'gene_symbols' or 'gene_ids'",
    ),
    sep: str = typer.Option(
        ",",
        "--sep",
        help="Separator for output file (default: ',')",
    ),
    version: bool = typer.Option(
        False,
        "-v",
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    """Convert 10X CellRanger output to EEISP input matrix format using scanpy."""
    _ = version
    
    convert_10x_to_matrix(
        input_path=input_path,
        output_path=output,
        input_format=format,
        var_names=var_names,
        sep=sep,
    )


def main() -> None:
    app()


if __name__ == "__main__":
    main()

