# EEISP

[![PyPI](https://img.shields.io/pypi/v/eeisp?logo=pypi&logoColor=white)](https://pypi.org/project/eeisp/)

EEISP identifies gene pairs that are codependent and mutually exclusive from single-cell RNA-seq data. 

## 0. Changelog

See [Changelog](CHANGELOG.md)

## 1. Installation

```bash
pip3 install -U eeisp

# Optional: Install with additional features
pip install 'eeisp[scanpy]'  # For 10X CellRanger conversion
pip install 'eeisp[gpu]'      # For GPU acceleration (requires CUDA)
```

### Install as a `uv` tool (global)

If you use `uv`, you can install EEISP as an isolated tool (recommended for CLI usage):

```
uv tool install eeisp
```

This will expose the commands on your `PATH` (e.g. `eeisp`, `add-names`, `heatmap`, `louvain`, `louvain-signed`).

### Run once without installing (uvx)

```
uvx eeisp --help
uvx --from eeisp add-names --help
```

## 2. Usage

EEISP takes a read count matrix as an input, in which rows and columns represent genes and cells, respectively. A gzipped file (.gz) is also acceptable.

   0. (Optional) Convert CellRanger output to an input matrix (requires scanpy: `pip install 'eeisp[scanpy]'`)
       ```bash
       # From directory format (matrix.mtx)
       eeisp-convert-10x outs/filtered_feature_bc_matrix/ matrix.txt
       
       # Or from .h5 format
       eeisp-convert-10x filtered_feature_bc_matrix.h5 matrix.txt --format h5
       ```

   1.  `eeisp` calculates the CDI and EEI scores for all gene pairs. The output contains lists of gene pairs that have CDI or EEI values above the specified threshold and the tables of degree distribution.
       ```
         usage: eeisp [-h] [--threCDI THRECDI] [--threEEI THREEEI] [--tsv] [--gpu] [-p THREADS] [-v] matrix output

         positional arguments:
           matrix                Input matrix
           output                Output prefix

         optional arguments:
           -h, --help            show this help message and exit
           --threCDI THRECDI     Threshold for CDI (default: 20.0)
           --threEEI THREEEI     Threshold for EEI (default: 10.0)
           --tsv                 Specify when the input file is tab-delimited (.tsv)
           --gpu                 GPU mode
           -p THREADS, --threads THREADS  number of threads (default: 2)
           -v, --version         show program's version number and exit
       ```  
   2. `eeisp_add_genename_from_geneid` add Gene Names (Symbols) to the output files of `eeisp`.
        ```
         usage: eeisp_add_genename_from_geneid [-h] [--i_id I_ID] [--i_name I_NAME] input output genelist

         positional arguments:
           input            Input matrix
           output           Output prefix
           genelist         Gene list

         optional arguments:
           -h, --help       show this help message and exit
           --i_id I_ID      column number of gene id (default: 0)
           --i_name I_NAME  column number of gene name (default: 1)
       ```

## 3. Tutorial

The sample data is included in `data/sample` directory.
   * `data.txt`: the input matrix of scRNA-seq data.
   * `genelidlist.txt`: the gene list for `eeisp_add_genename_from_geneid`.

### Converting 10X CellRanger data (optional)

If you have 10X CellRanger output, convert it first:

```bash
# Install with scanpy support
uv sync --extra scanpy

# Convert CellRanger output
uv run eeisp-convert-10x outs/filtered_feature_bc_matrix/ data/input/matrix.txt
```

### Running EEISP analysis


Run with `uv` (recommended for this repo checkout). Outputs go to `data/output/` by default if you pass a plain prefix like `Sample`:

```
uv sync
uv run eeisp data/sample/data.txt Sample --threCDI 0.5 --threEEI 0.5 -p 8
uv run add-names data/output/Sample_CDI_score_data_thre0.5.txt data/output/Sample_CDI_score_data_thre0.5.addgenename.txt data/sample/geneidlist.txt
uv run add-names data/output/Sample_EEI_score_data_thre0.5.txt data/output/Sample_EEI_score_data_thre0.5.addgenename.txt data/sample/geneidlist.txt
```

Supply `--gpu` option for GPU computation:

```bash
# Install with GPU support
uv sync --extra gpu
# or: pip install 'eeisp[gpu]'

# Run with GPU
uv run eeisp data/sample/data.txt Sample --threCDI 0.5 --threEEI 0.5 -p 8 --gpu
```
    
(Note: Since GPU computation covers a part of eeisp, it is better to use multiple CPUs even in `--gpu` mode for the fast computation.)

Output files are:
```
   data/output/Sample_CDI_score_data_thre0.5.txt            # A list of gene pairs with CDI score.
   data/output/Sample_CDI_degree_distribution.tsv           # A table of the number of CDI degree and genes.
   data/output/Sample_EEI_score_data_thre0.5.txt            # A list of gene pairs with EEI scores.
   data/output/Sample_EEI_degree_distribution.tsv           # A table of the number of EEI degree and genes.
```
The output files might include gene ids only. 

```
   $ head Sample_CDI_score_data_thre0.5.txt
   2       7       ESG000003       ESG000008       0.96384320244841
   0       1       ESG000001       ESG000002       0.6852891560232545
   0       6       ESG000001       ESG000007       0.6852891560232545
   7       8       ESG000008       ESG000009       0.6852891560232545
   3       9       ESG000004       ESG000010       0.6469554204484568
   4       6       ESG100005       ESG000007       0.5258703930217091
```

If you want to add gene names (Symbols), use `eeisp_add_genename_from_geneid` with `geneidlist.txt`, which contains the pairs of gene ids and names.

```
uv run add-names \
    data/output/Sample_CDI_score_data_thre0.5.txt \
    data/output/Sample_CDI_score_data_thre0.5.addgenename.txt \
    data/sample/geneidlist.txt
uv run add-names \
    data/output/Sample_EEI_score_data_thre0.5.txt \
    data/output/Sample_EEI_score_data_thre0.5.addgenename.txt \
    data/sample/geneidlist.txt
```

The output files include gene names.

```
   $ head Sample_CDI_score_data_thre0.5.addgenename.txt
   2       7       ESG000003       ESG000008       OR4F5   FO538757.3      0.96384320244841
   0       1       ESG000001       ESG000002       RP11-34P13.3    FAM138A 0.6852891560232545
   0       6       ESG000001       ESG000007       RP11-34P13.3    RP11-34P13.9    0.6852891560232545
   7       8       ESG000008       ESG000009       FO538757.3      FO538757.2      0.6852891560232545
   3       9       ESG000004       ESG000010       RP11-34P13.7    AP006222.2      0.6469554204484568
   4       6       ESG100005       ESG000007       RP11-34P13.8    RP11-34P13.9    0.5258703930217091
```

## 4. Reference

Nakajima N., Hayashi T., Fujiki K., Shirahige K., Akiyama T., Akutsu T. and Nakato R., [Codependency and mutual exclusivity for gene community detection from sparse single-cell transcriptome data](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab601/6324613), *Nucleic Acids Research*, 2021.
