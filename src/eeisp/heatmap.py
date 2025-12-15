from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
from numpy.typing import NDArray


def plot_heatmap(*, input_path: Path, output_path: Path) -> None:
    df = pl.read_csv(input_path, separator="\t", has_header=False)
    if df.width < 5:
        raise ValueError("Score file must have at least 5 columns (i, j, label1, label2, value)")

    df = df.select([pl.col(c) for c in df.columns[:5]]).rename(
        {
            df.columns[0]: "row_index",
            df.columns[1]: "col_index",
            df.columns[2]: "row_label",
            df.columns[3]: "col_label",
            df.columns[4]: "value",
        }
    )

    row_idx = df.select(pl.col("row_index").cast(pl.Int64)).to_series().to_numpy()
    col_idx = df.select(pl.col("col_index").cast(pl.Int64)).to_series().to_numpy()
    values = df.select(pl.col("value").cast(pl.Float64)).to_series().to_numpy()

    n_genes = int(max(int(row_idx.max(initial=0)), int(col_idx.max(initial=0))) + 1)
    matrix: NDArray[np.float64] = np.zeros((n_genes, n_genes), dtype=np.float64)

    matrix[row_idx, col_idx] = values
    matrix[col_idx, row_idx] = values

    plt.figure(figsize=(12, 10))
    sns.heatmap(matrix, cmap="Blues")
    plt.savefig(output_path)
    plt.close()


