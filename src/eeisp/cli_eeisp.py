from __future__ import annotations

from pathlib import Path

import typer

from eeisp.runner import run_eeisp
from eeisp.version import get_version

app = typer.Typer(add_completion=False, help="EEISP: compute CDI and EEI from an expression matrix.")


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"eeisp {get_version()}")
        raise typer.Exit()


@app.command()
def cli(
    matrix: Path = typer.Argument(..., help="Input matrix (genes x cells)."),
    output: str = typer.Argument(..., help="Output prefix. If no path is given, outputs go to data/output/."),
    threCDI: float = typer.Option(10.0, "--threCDI", help="Threshold for CDI (default: 10.0)"),
    threEEI: float = typer.Option(5.0, "--threEEI", help="Threshold for EEI (default: 5.0)"),
    tsv: bool = typer.Option(False, "--tsv", help="Specify when the input file is tab-delimited (.tsv)"),
    gpu: bool = typer.Option(False, "--gpu", help="GPU mode (requires cupy)"),
    CDIonly: bool = typer.Option(False, "--CDIonly", help="Calculate CDI only"),
    EEIonly: bool = typer.Option(False, "--EEIonly", help="Calculate EEI only"),
    threads: int = typer.Option(2, "-p", "--threads", help="number of threads (default: 2)"),
    version: bool = typer.Option(
        False, "-v", "--version", callback=_version_callback, is_eager=True, help="Show version and exit."
    ),
) -> None:
    _ = version
    output_path = Path(output)
    # Default behavior: if user passes just a prefix (no directory), write into data/output/<prefix>_*.*
    if output_path.parent == Path("."):
        output_path = Path("data") / "output" / output_path.name

    try:
        run_eeisp(
            matrix,
            output_prefix=output_path,
            thre_cdi=threCDI,
            thre_eei=threEEI,
            tsv=tsv,
            gpu=gpu,
            cdi_only=CDIonly,
            eei_only=EEIonly,
            threads=threads,
        )
    except ImportError as e:
        typer.echo(f"Import error: {e}", err=True)
        typer.echo("If you used --gpu, install a compatible cupy build for your CUDA setup.", err=True)
        raise typer.Exit(code=1) from e


def main() -> None:
    app()


if __name__ == "__main__":
    main()


