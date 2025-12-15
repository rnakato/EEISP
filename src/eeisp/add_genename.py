from __future__ import annotations

from pathlib import Path

import polars as pl

from eeisp.io import read_gene_id_name_map


def add_genenames_from_geneids(
    *,
    input_path: Path,
    output_path: Path,
    genelist_path: Path,
    gene_id_col: int = 0,
    gene_name_col: int = 1,
) -> None:
    mapping = read_gene_id_name_map(
        genelist_path, sep="\t", gene_id_col=gene_id_col, gene_name_col=gene_name_col
    )

    df = pl.read_csv(input_path, separator="\t", has_header=False)
    if df.width < 5:
        raise ValueError("Input score file must have at least 5 tab-delimited columns")

    df = df.select([pl.col(c) for c in df.columns[:5]]).rename(
        {
            df.columns[0]: "i",
            df.columns[1]: "j",
            df.columns[2]: "geneid1",
            df.columns[3]: "geneid2",
            df.columns[4]: "val",
        }
    )

    df = df.with_columns(
        [
            pl.col("geneid1").cast(pl.Utf8),
            pl.col("geneid2").cast(pl.Utf8),
        ]
    ).with_columns(
        [
            pl.col("geneid1").replace(mapping, default=None).fill_null("").alias("genename1"),
            pl.col("geneid2").replace(mapping, default=None).fill_null("").alias("genename2"),
        ]
    )

    df = df.select(["i", "j", "geneid1", "geneid2", "genename1", "genename2", "val"])
    df.write_csv(output_path, separator="\t", include_header=False)


