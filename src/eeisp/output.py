from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray


def write_score_edges(
    matrix: NDArray[np.floating[Any]],
    *,
    threshold: float,
    score_name: str,
    output_prefix: Path,
    gene_ids: NDArray[np.str_],
) -> list[int]:
    """
    Writes `<output_prefix>_<score_name>_score_data_thre<threshold>.txt` (tab-delimited, no header).

    Each line: i, j, gene_id_i, gene_id_j, value
    """
    filename = f"{score_name}_score_data_thre{threshold}"
    prefix = str(output_prefix)
    pdf_path = Path(f"{prefix}_{filename}.pdf")
    data_path = Path(f"{prefix}_{filename}.txt")

    degree = (matrix > threshold).sum(axis=0).astype(int).tolist()

    upper = np.triu(matrix, k=1)
    ii, jj = np.where(upper > threshold)
    vals = matrix[ii, jj]

    if vals.size > 0:
        order = np.lexsort((jj, ii, -vals))
        ii = ii[order]
        jj = jj[order]
        vals = vals[order]

    plt.figure(figsize=(10, 6))
    if vals.size > 0:
        plt.hist(vals, bins=50, color="blue", alpha=0.7, edgecolor="black")
    else:
        plt.hist([], bins=50, color="blue", alpha=0.7, edgecolor="black")
    plt.title(f"Distribution of {output_prefix}")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.grid(True)
    plt.savefig(pdf_path)
    plt.close()

    print(f"output degree data in {data_path}")
    print(f"number of gene pairs over threshold (>{threshold}): {int(vals.size)}")
    with data_path.open("w", encoding="utf-8") as f:
        for i, j, v in zip(ii.tolist(), jj.tolist(), vals.tolist(), strict=False):
            f.write(f"{i}\t{j}\t{gene_ids[i]}\t{gene_ids[j]}\t{v}\n")

    return degree


def write_degree_distribution(
    degree: list[int],
    *,
    score_name: str,
    output_prefix: Path,
) -> None:
    max_value = max(degree)
    min_value = min(degree)
    value_range = max_value - min_value
    print(f"max degree:{max_value:.3f} min degree:{min_value:.3f} value_width={value_range:.3f}")

    # Preserve the legacy behavior (skip min_value, keep only degrees with non-zero counts).
    freq: list[tuple[int, int]] = []
    for a in range(min_value + 1, max_value + 1):
        cnt = degree.count(a)
        if cnt > 0:
            freq.append((a, cnt))

    out_path = Path(f"{str(output_prefix)}_{score_name}_degree_distribution.tsv")
    with out_path.open("w", encoding="utf-8") as f:
        f.write("\tLog_Degree\tLog_The number of genes\tDegree\tThe number of genes\n")
        for idx, (deg, cnt) in enumerate(freq):
            log_deg = math.log(deg)
            log_cnt = math.log(cnt)
            f.write(f"{idx}\t{log_deg}\t{log_cnt}\t{deg}\t{cnt}\n")


