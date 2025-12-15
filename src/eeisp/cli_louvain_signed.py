from __future__ import annotations

from pathlib import Path

import networkx as nx
import typer

from eeisp.louvain_signed import LouvainSigned
from eeisp.network import load_graph_from_file
from eeisp.version import get_version

app = typer.Typer(add_completion=False, help="Signed community detection from CDI (positive) and EEI (negative).")


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"eeisp.LouvainSigned {get_version()}")
        raise typer.Exit()


@app.command()
def cli(
    cdi_score_data: Path = typer.Argument(..., help="CDI_score_data", metavar="CDI_score_data"),
    eei_score_data: Path = typer.Argument(..., help="EEI_score_data", metavar="EEI_score_data"),
    thre_CDI: float = typer.Option(10.0, "--thre_CDI", help="Threshold of CDI (default: 10)"),
    thre_EEI: float = typer.Option(5.0, "--thre_EEI", help="Threshold of EEI (default: 5)"),
    alpha: float = typer.Option(0.5, "--alpha", help="alpha parameter (from 0 to 1, default: 0.5)"),
    resolution: float = typer.Option(1.0, "--resolution", help="resolution for louvain (default: 1.0)"),
    seed: int | None = typer.Option(None, "--seed", help="seed for LouvainSigned"),
    version: bool = typer.Option(
        False, "-v", "--version", callback=_version_callback, is_eager=True, help="Show version and exit."
    ),
) -> None:
    _ = version
    g_positive = load_graph_from_file(cdi_score_data, threshold=thre_CDI)
    g_negative = load_graph_from_file(eei_score_data, threshold=thre_EEI)

    print(f"G_positive: Nodes={g_positive.number_of_nodes()}, Edges={g_positive.number_of_edges()}")
    print(f"G_negative: Nodes={g_negative.number_of_nodes()}, Edges={g_negative.number_of_edges()}")
    print(
        f"G: Nodes={nx.compose(g_positive, g_negative).number_of_nodes()}, "
        f"Edges={nx.compose(g_positive, g_negative).number_of_edges()}"
    )

    l = LouvainSigned(g_positive, g_negative)
    partition = l.best_partition(alpha=alpha, resolution=resolution, seed=seed)
    print(f"partition: {partition}")


def main() -> None:
    app()


if __name__ == "__main__":
    main()


