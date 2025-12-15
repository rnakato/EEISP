from __future__ import annotations

import gc
import time
from pathlib import Path

from eeisp.io import filter_nonzero_genes, read_expression_matrix
from eeisp.matrix import generate_cdi_matrix, generate_eei_matrix
from eeisp.output import write_degree_distribution, write_score_edges


def run_eeisp(
    matrix_path: Path,
    *,
    output_prefix: Path,
    thre_cdi: float,
    thre_eei: float,
    tsv: bool,
    gpu: bool,
    cdi_only: bool,
    eei_only: bool,
    threads: int,
) -> None:
    startt = time.time()

    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    values, gene_ids = read_expression_matrix(matrix_path, tsv=tsv)

    n_genes_total = int(values.shape[0])
    n_cells = int(values.shape[1])
    print("number of cells: ", n_cells)
    print("number of genes: ", n_genes_total)

    values, gene_ids = filter_nonzero_genes(values, gene_ids)
    gc.collect()

    ngene = int(values.shape[0])
    print("number of nonzero genes: ", ngene)
    print("-----------------------------------------------")

    if not eei_only:
        print("Calculating CDI...")
        print("using GPU for CDI calculation." if gpu else "using CPU for CDI calculation.")
        cdi = generate_cdi_matrix(values, threads=threads, use_gpu=gpu)
        degree_cdi = write_score_edges(
            cdi,
            threshold=thre_cdi,
            score_name="CDI",
            output_prefix=output_prefix,
            gene_ids=gene_ids,
        )
        write_degree_distribution(degree_cdi, score_name="CDI", output_prefix=output_prefix)
        del cdi, degree_cdi
        gc.collect()
        print("Finish to calculate CDI!")

    if not cdi_only:
        print("Calculating EEI...")
        print("using GPU for EEI calculation." if gpu else "using CPU for EEI calculation.")
        eei = generate_eei_matrix(values, threads=threads, use_gpu=gpu)
        del values
        degree_eei = write_score_edges(
            eei,
            threshold=thre_eei,
            score_name="EEI",
            output_prefix=output_prefix,
            gene_ids=gene_ids,
        )
        write_degree_distribution(degree_eei, score_name="EEI", output_prefix=output_prefix)
        del eei, degree_eei
        gc.collect()
        print("Finish to calculate EEI!")

    elapsed_time = time.time() - startt
    print(f"Elapsed_time:{elapsed_time}[sec]")
    print("*************************************************************")


