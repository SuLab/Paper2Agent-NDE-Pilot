# decoupler bulk-enrichment MCP server

This MCP server turns the decoupler **Bulk Enrichment Analysis** tutorial into five tools:
- filter genes;
- explore sample variability;
- test differential expression with PyDESeq2;
- infer transcription-factor (TF), pathway and Hallmark activities with `ulm`;
- inspect one TF, pathway or gene set through its target genes.

Each tool calls [scverse/decoupler](https://github.com/scverse/decoupler) 2.3.1, pinned to commit `27a1be7c75d88306bf3c2e4d0e7fb4b42d48f78a`, and PyDESeq2 0.5.4. None of their methods are reimplemented.

- **Tutorial:** `docs/notebooks/bulk/rna.ipynb` in [scverse/decoupler-tutorials](https://github.com/scverse/decoupler-tutorials), commit `306be6673a5df5baf780e44ea929391616dbe5d8`.
- **Paper:** Badia-i-Mompel *et al.*, "decoupleR: ensemble of computational methods to infer biological activities from omics data", *Bioinformatics Advances* 2, vbac016 (2022). <https://doi.org/10.1093/bioadv/vbac016>

## What it covers and what it doesn't

**Covered.** The bulk RNA-seq tutorial's workflow, run on *your* samples and metadata. Human networks:
- **CollecTRI:** TF regulons.
- **PROGENy:** 14 pathway footprints.
- **MSigDB Hallmark:** 50 gene sets.

**Not covered:**
- The single-cell, pseudobulk, spatial, benchmarking and OmniPath tutorials.
- Non-human organisms, which would need `dc.op.translate`.
- Enrichment methods other than `ulm`.
- Arbitrary custom gene-set collections. `network_path` only accepts tables in the CollecTRI, PROGENy or Hallmark format; see below.
- Gene-identifier mapping. The networks use **HGNC gene symbols**, so convert Ensembl or Entrez IDs to symbols before using these tools.

**Input contract.** An AnnData `.h5ad` with:
- samples as observations, genes as variables;
- **raw integer counts** in `.X`;
- sample metadata in `.obs`.

Conditions are always read from a metadata column you name. They are never inferred from sample order.

## Install

Requirements:
- Linux x86_64, the tested platform;
- [uv](https://docs.astral.sh/uv/);
- git;
- network access to GitHub and PyPI.

Run these from the extracted `decoupler-mcp/` folder:

```bash
uv venv .venv --python 3.13
uv pip install --python .venv/bin/python -r src/requirements.txt
uv pip check --python .venv/bin/python
```

`src/requirements.txt` pins the full tested dependency closure:
- CPython 3.13.12
- fastmcp 4.0.3
- PyDESeq2 0.5.4
- scanpy 1.12.4
- anndata 0.13.4
- pandas 3.1.0
- numpy 2.5.4
- numba 0.68.1
- igraph 1.0.0
- decoupler, from the pinned git commit

## Start the server (stdio)

```bash
.venv/bin/python src/decoupler_mcp.py
```

Client configuration (use your absolute paths):

```json
{
  "mcpServers": {
    "decoupler": {
      "command": "/ABSOLUTE/PATH/decoupler-mcp/.venv/bin/python",
      "args": ["/ABSOLUTE/PATH/decoupler-mcp/src/decoupler_mcp.py"]
    }
  }
}
```

Claude Code: the syntax below comes from `fastmcp install claude-code --help` (fastmcp 4.0.3). It was not executed during validation.

```bash
.venv/bin/fastmcp install claude-code src/decoupler_mcp.py --name decoupler \
  --python 3.13 --with-requirements src/requirements.txt
```

## Networks: download vs. saved copy

By default, the networks are downloaded on each call:
- CollecTRI from Zenodo;
- PROGENy from omnipathdb.org;
- Hallmark from static.omnipathdb.org.

These resources can change over time. Every call to the two network tools saves the table it actually used as `<network>_network_used.csv`. Pass that file back as `network_path` (CSV, TSV or Parquet with `source`, `target` and an optional `weight` column) to reproduce a result exactly or to work offline.

## Tools and how they chain

```
filter_genes_by_expression ─► pydeseq2_test_differential_expression ─► score_activities_ulm ─► inspect_source_targets
            └─► rank_pcs_by_metadata (optional: which metadata drives variation? e.g. a batch/donor term for the design)
```

Each call writes into a fresh subdirectory:
- under `output_dir` if you pass one;
- otherwise under `<package>/tmp/outputs/bulk_enrichment/`.

It returns `message`, `reference`, absolute `artifacts` paths, and compact summaries such as shapes, counts and top sources.

### `decoupler_filter_genes_by_expression`
Removes lowly expressed genes using edgeR `filterByExpr` logic: `dc.pp.filter_by_expr`, plus the diagnostic plot `dc.pl.filter_by_expr`.

- **Required:**
  - `data_path`
  - `group_key`: the obs column of experimental groups.
- **Options** (defaults are the tutorial's):
  - `min_count=10`
  - `min_total_count=15`
  - `large_n=10`
  - `min_prop=0.7`
- **Outputs:** filtered `.h5ad`, diagnostic PNG, and gene counts before and after.
- **Note:** upstream also allows `group=None`, which treats all samples as one group. This tool requires a group column.

### `decoupler_rank_pcs_by_metadata`
Runs the tutorial's PCA: normalise to 1e4, log1p, scale (max 10), PCA. It then calls `dc.tl.rankby_obsm` to test each principal component (PC) against each metadata column. Use this to see whether condition, donor or batch drives variation, and so whether to add a covariate to the DE design.

- **Required:**
  - `data_path`
- **Options:**
  - `obs_keys`: columns to test. Default: all.
  - `color_keys`: columns for the PCA scatter.
  - `nvar=5`
- **Outputs:**
  - association table CSV;
  - PC scores;
  - variance-ratio plot, association heatmap and PCA scatter PNGs.

### `pydeseq2_test_differential_expression`
Runs the tutorial's PyDESeq2 Wald test: `DeseqDataSet(..., refit_cooks=True)`, then `.deseq2()`, then `DeseqStats(contrast=[condition_key, group, reference]).summary()`.

- **Required:**
  - `data_path`
  - `condition_key`
  - `group`: the numerator of the fold change.
  - `reference`
- **Options:**
  - `design`: default `~<condition_key>`. Add covariates, e.g. `~donor + condition` for paired designs.
  - `n_cpus=8`: results do not depend on it.
- **Outputs:**
  - gene-indexed results CSV, with columns `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`;
  - volcano PNG.
- **Input checks:** the tool rejects a missing column, a missing group level, and counts that are negative or non-integer.

### `decoupler_score_activities_ulm`
Infers activities for one contrast by running `dc.mt.ulm` on a DE statistic.

- **Required:**
  - `de_results_path`: the CSV from the DE tool.
  - `network`: `collectri`, `progeny` or `hallmark`.
- **Options:**
  - `stat_column="stat"`: the Wald statistic, as the tutorial recommends.
  - `padj_threshold=0.05`
  - `contrast_name`
  - `barplot_top`
  - `dotplot_top=10`
  - `network_path`
- **Outputs:**
  - scores and adjusted p-values for **all** sources;
  - the significant subset;
  - the network used;
  - bar plot and dot plot PNGs.

### `decoupler_inspect_source_targets`
Explains one source's activity through its target genes, using:
- a target volcano (`dc.pl.volcano`);
- a weight-vs-statistic scatter (`dc.pl.source_targets`);
- leading-edge plots (`dc.pl.leading_edge`), split by weight sign for weighted networks, as the tutorial does for PROGENy;
- a network plot (`dc.pl.network`).

- **Required:**
  - `de_results_path`
  - `network`
  - `source`: e.g. a top hit from the scoring tool.
- **Options:**
  - `scores_path`: the scores CSV from the scoring tool. Required for the network plot.
  - `network_sources`
  - `stat_column`
  - `network_path`
- **Outputs:** leading-edge gene lists, the network used, and PNGs.

## Validation performed

The tools' correctness was checked against independent upstream evidence: the real tutorial notebook, run end to end with live downloads, and direct decoupler and PyDESeq2 calls on the same inputs. A separate verifier agent wrote the tests; it did not write the wrappers. The tutorial dataset is GSE151251: hepatic stellate cells, 3 controls vs 3 TGF-β treated.

**Filtering**
- 58,674 genes reduce to 19,510, identical to the tutorial and in the same order.

**PCA association**
- The association table and PC scores match. The comparison ignores sign, because PCA axis signs are arbitrary.

**Differential expression**
- The full PyDESeq2 results table is bit-for-bit equal to the tutorial's.
- Swapping the contrast matched a direct PyDESeq2 call bit for bit.
- Relabelling and reordering the metadata matched a direct call bit for bit. PyDESeq2 itself is not exactly invariant to sample order, so the comparison with the tutorial's table uses a correlation threshold.
- A two-factor design (`~batch + condition`) also matched a direct call bit for bit.

**Activity scoring (ulm)**
- Scores and adjusted p-values are bit-identical for all three networks.
- 143 TFs, 12 pathways and 30 Hallmark sets are significant.
- Top TFs are SRF, MYOCD and SMAD3 (up) and IRF1, SKIL and NR1H4 (down). TGF-β has the top pathway score, 59.36.

**Leading edges**
- TGF-β PROGENy: 581 (+) and 292 (−) genes.
- Hallmark TGF_BETA_SIGNALING: 19 genes.

**Errors, repeat calls and composition**
- Missing files, columns and groups, invalid counts, an unknown source and malformed networks all return MCP tool errors.
- Repeated calls write fresh directories.
- The filter → DE → scoring → inspect chain was tested end to end.

**Packaged server over real stdio**
- 27 acceptance cases (16 successful-call cases and 11 error cases) pass in the development environment.
- They also pass in a fresh environment installed only from `src/requirements.txt`.

### Limits

- **Live network downloads are not used for exact checks.** Exact comparisons use network snapshots saved on 2026-10-27. The live-download tests check structure and rough agreement with the snapshots.
- **Two figures vary between runs.** In the volcano and source–target plots, label positions move slightly from run to run because of adjustText. Their underlying data is identical.
- **`sc.set_figure_params(figsize=(3, 3), frameon=False)` is set once at import**, as in the tutorial. It applies to the whole server process.
- **Only Linux x86_64 was tested.**
- **These checks only show that the selected tools faithfully wrap decoupler and PyDESeq2.** They do not re-validate the methods, the network resources, or the paper's benchmark results.
