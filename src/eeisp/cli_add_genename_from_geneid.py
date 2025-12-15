from __future__ import annotations

from pathlib import Path

import typer

from eeisp.add_genename import add_genenames_from_geneids
from eeisp.version import get_version

app = typer.Typer(add_completion=False, help="Add gene names to EEISP score files.")


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"eeisp_add_genename_from_geneid {get_version()}")
        raise typer.Exit()


@app.command()
def cli(
    input: Path = typer.Argument(..., help="Input score file."),
    output: Path = typer.Argument(..., help="Output path."),
    genelist: Path = typer.Argument(..., help="Gene list (id<tab>name)."),
    i_id: int = typer.Option(0, "--i_id", help="column number of gene id (default: 0)"),
    i_name: int = typer.Option(1, "--i_name", help="column number of gene name (default: 1)"),
    version: bool = typer.Option(
        False, "-v", "--version", callback=_version_callback, is_eager=True, help="Show version and exit."
    ),
) -> None:
    _ = version
    add_genenames_from_geneids(
        input_path=input,
        output_path=output,
        genelist_path=genelist,
        gene_id_col=i_id,
        gene_name_col=i_name,
    )


def main() -> None:
    app()


if __name__ == "__main__":
    main()


