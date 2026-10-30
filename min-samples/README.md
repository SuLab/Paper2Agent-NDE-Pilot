# How few samples or cells can each tool handle?

These tests set the size filters in the NDE searches.

## decoupler server (PyDESeq2 + ulm)

[`deseq_small_n.py`](deseq_small_n.py) simulates negative-binomial counts (2,000 genes, 200 of them 4-fold up) at small group sizes. It runs them through PyDESeq2, called as the server calls it, and through decoupler's ulm. Output: [`output/deseq_small_n.txt`](output/deseq_small_n.txt).

- **1 vs 1:** fails. With as many samples as design terms, PyDESeq2 raises "no replicates to estimate the dispersion" (`pydeseq2/utils.py`).
- **1–3 residual degrees of freedom:** runs, with a warning (`dds.py:859` warns when samples − design variables ≤ 3).
- **3 vs 3 unpaired:** the smallest design that runs cleanly.
- **Outliers:** Cook's-distance outlier flagging only applies to groups with at least 3 replicates (`dds.py:1077–1096`), and refitting needs 7.
- **ulm:** the significance test uses n_genes − 2 degrees of freedom (`mt/_ulm.py:95`), so the number of samples does not limit it.

The NDE bulk search therefore requires at least 3 samples (1 vs 2); 3 vs 3 or more is better.

## CellTypist server

[`celltypist_small_n.py`](celltypist_small_n.py) annotates random subsets of PBMC 3k (5 to 200 cells), with and without majority voting, and trains on small subsets. Output: [`output/celltypist_small_n.txt`](output/celltypist_small_n.txt).

- **Annotation:** there is no minimum.
- **Majority voting:** skipped at 50 cells or fewer (`annotate.py:88`). The MCP server reports when this happens.
- **Training:** needs at least 2 cell types. It trained on 2 types × 1 cell and failed on 1 type.
- **From the source code, not tested here:** feature selection needs more genes than `top_genes`, and mini-batch training needs more cells than `batch_size`.

Versions: PyDESeq2 0.5.4, decoupler 2.3.1, CellTypist 1.7.1. File and line references are to those releases.
