from __future__ import annotations

from pathlib import Path

import typer

from eeisp.heatmap import plot_heatmap
from eeisp.version import get_version

app = typer.Typer(add_completion=False, help="Plot a heatmap from an EEISP score file.")


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"eeisp.heatmap {get_version()}")
        raise typer.Exit()


@app.command()
def cli(
    input: Path = typer.Argument(..., help="Score_data (*_[CDI|EEI]_score_data_*.txt)"),
    output: Path = typer.Argument(..., help="Output file name"),
    version: bool = typer.Option(
        False, "-v", "--version", callback=_version_callback, is_eager=True, help="Show version and exit."
    ),
) -> None:
    _ = version
    plot_heatmap(input_path=input, output_path=output)


def main() -> None:
    app()


if __name__ == "__main__":
    main()


