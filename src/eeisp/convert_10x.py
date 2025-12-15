from __future__ import annotations

import csv
from pathlib import Path
from typing import Literal

import numpy as np


def convert_10x_to_matrix(
    input_path: Path,
    output_path: Path,
    input_format: Literal["mtx", "h5"] = "mtx",
    var_names: str = "gene_symbols",
    sep: str = ",",
) -> None:
    """
    Convert 10X CellRanger output to EEISP input matrix format.
    
    Args:
        input_path: Path to 10X data (directory for mtx, file for h5)
        output_path: Output file path
        input_format: Format of input data ("mtx" or "h5")
        var_names: Variable names to use ('gene_symbols' or 'gene_ids')
        sep: Separator for output file (default: ',')
    """
    try:
        import scanpy as sc
    except ImportError as e:
        raise ImportError(
            "scanpy is required for 10X conversion. Install it with:\n"
            "  uv sync --extra scanpy\n"
            "or:\n"
            "  pip install 'eeisp[scanpy]'"
        ) from e
    
    # Read 10X data based on format
    if input_format == "h5":
        adata = sc.read_10x_h5(input_path)
    else:
        adata = sc.read_10x_mtx(input_path, var_names=var_names, cache=True)

    # adata.X is (cells x genes). EEISP expects (genes x cells).
    x = adata.X
    x_t = x.T

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=sep)
        writer.writerow([""] + [str(c) for c in adata.obs_names])

        # Iterate genes (rows) to avoid materializing the full dense matrix at once.
        gene_names = [str(g) for g in adata.var_names]
        for i, gene in enumerate(gene_names):
            row = x_t[i]
            if hasattr(row, "toarray"):
                arr = np.asarray(row.toarray()).ravel()
            else:
                arr = np.asarray(row).ravel()
            writer.writerow([gene] + arr.tolist())

    print(f"Converted {int(adata.n_vars)} genes × {int(adata.n_obs)} cells")
    print(f"Output saved to: {output_path}")

