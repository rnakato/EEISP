from __future__ import annotations

from pathlib import Path

import typer

from eeisp.louvain_signed import LouvainSigned
from eeisp.network import load_graph_from_file
from eeisp.version import get_version

app = typer.Typer(add_completion=False, help="Community detection from a CDI score file (unsigned).")


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"eeisp.Louvain {get_version()}")
        raise typer.Exit()


@app.command()
def cli(
    cdi_score_data: Path = typer.Argument(..., help="CDI_score_data", metavar="CDI_score_data"),
    thre_CDI: float = typer.Option(10.0, "--thre_CDI", help="Threshold of CDI (default: 10)"),
    resolution: float = typer.Option(1.0, "--resolution", help="resolution for louvain (default: 1.0)"),
    seed: int | None = typer.Option(None, "--seed", help="seed for LouvainSigned"),
    version: bool = typer.Option(
        False, "-v", "--version", callback=_version_callback, is_eager=True, help="Show version and exit."
    ),
) -> None:
    _ = version
    g = load_graph_from_file(cdi_score_data, threshold=thre_CDI)
    print(f"G: Nodes={g.number_of_nodes()}, Edges={g.number_of_edges()}")
    l = LouvainSigned(g, g)
    partition = l.best_partition(alpha=1.0, resolution=resolution, seed=seed)
    print(f"partition: {partition}")


def main() -> None:
    app()


if __name__ == "__main__":
    main()


