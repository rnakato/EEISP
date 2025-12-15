from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from numpy.typing import NDArray


def read_expression_matrix(path: Path, *, tsv: bool) -> tuple[NDArray[np.number[Any]], NDArray[np.str_]]:
    """
    Read an expression matrix where rows are genes and columns are cells.

    The first column is treated as the gene identifier column (like pandas `index_col=0`).
    """
    sep = "\t" if tsv else ","
    df = pl.read_csv(path, separator=sep, has_header=True)
    if not df.columns:
        raise ValueError("Input matrix has no columns")
    gene_col = df.columns[0]

    gene_ids = df.select(pl.col(gene_col).cast(pl.Utf8)).to_series().to_numpy()

    values_df = df.drop(gene_col)
    # Ensure numeric (count matrix). If there are empty strings, treat as null then fill with 0.
    values_df = values_df.with_columns(
        [pl.all().cast(pl.Float64, strict=False).fill_null(0.0).fill_nan(0.0)]
    )
    values = values_df.to_numpy()

    return values, gene_ids


def filter_nonzero_genes(
    values: NDArray[np.number[Any]],
    gene_ids: NDArray[np.str_],
) -> tuple[NDArray[np.number[Any]], NDArray[np.str_]]:
    mask = np.any(values > 0, axis=1)
    return values[mask], gene_ids[mask]


def read_gene_id_name_map(
    path: Path,
    *,
    sep: str = "\t",
    gene_id_col: int = 0,
    gene_name_col: int = 1,
) -> dict[str, str]:
    """
    Reads a two-column (or multi-column) table mapping gene ID to gene name.
    """
    df = pl.read_csv(path, separator=sep, has_header=False)
    if df.width <= max(gene_id_col, gene_name_col):
        raise ValueError("genelist file has fewer columns than requested")

    id_series = df.to_series(gene_id_col).cast(pl.Utf8)
    name_series = df.to_series(gene_name_col).cast(pl.Utf8)
    return dict(zip(id_series.to_list(), name_series.to_list(), strict=False))


