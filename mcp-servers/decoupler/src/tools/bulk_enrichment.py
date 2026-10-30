"""Tools extracted from scverse/decoupler docs/notebooks/bulk/rna.ipynb (Bulk Enrichment Analysis).

Workflow of the tutorial, split into composable tools:
filter lowly expressed genes (edgeR filterByExpr logic) -> explore PCA / metadata associations ->
PyDESeq2 Wald test between two conditions -> ulm enrichment of the gene-level statistic with
CollecTRI TF regulons, PROGENy pathway footprints or MSigDB Hallmark gene sets -> inspect the
target genes behind one source.

Composition contract:
- ``decoupler_filter_genes_by_expression`` writes a filtered raw-count ``.h5ad`` that is the input
  of ``decoupler_rank_pcs_by_metadata`` and ``pydeseq2_test_differential_expression``.
- ``pydeseq2_test_differential_expression`` writes the PyDESeq2 ``results_df`` as CSV with the gene
  index in the first column; that CSV is the ``de_results_path`` input of
  ``decoupler_score_activities_ulm`` and ``decoupler_inspect_source_targets``.
- ``decoupler_score_activities_ulm`` writes the wide ulm score table (one contrast row, one column
  per source); it is the optional ``scores_path`` input of ``decoupler_inspect_source_targets``.
"""

import threading
import uuid
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Annotated, Literal

import matplotlib

matplotlib.use("Agg")

import anndata as ad  # noqa: E402
import decoupler as dc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sps  # noqa: E402
from fastmcp import FastMCP  # noqa: E402
from pydeseq2.dds import DefaultInference, DeseqDataSet  # noqa: E402
from pydeseq2.ds import DeseqStats  # noqa: E402

REFERENCE = "https://github.com/scverse/decoupler-tutorials/blob/306be6673a5df5baf780e44ea929391616dbe5d8/bulk/rna.ipynb"
bulk_enrichment_mcp = FastMCP(name="bulk_enrichment")

# Tutorial setup (cell 2): scanpy figure defaults used for every figure of the notebook.
sc.set_figure_params(figsize=(3, 3), frameon=False)

_DEFAULT_OUTPUT_BASE = Path(__file__).resolve().parents[2] / "tmp" / "outputs" / "bulk_enrichment"

# pyplot keeps process-global figure state and FastMCP runs sync tools in worker threads,
# so figure creation, saving and closing are serialised.
_PLOT_LOCK = threading.Lock()

NetworkName = Literal["collectri", "progeny", "hallmark"]

# Prior-knowledge resources used in the tutorial (cells 30, 40, 52), always for human.
_NETWORK_LOADERS = {
    "collectri": dc.op.collectri,
    "progeny": dc.op.progeny,
    "hallmark": dc.op.hallmark,
}

# Plot settings taken from the tutorial cells: barplot cells 34 (CollecTRI), 44 (PROGENy),
# 56 (Hallmark); dotplot cells 46 (PROGENy) and 58 (Hallmark). The tutorial draws no CollecTRI
# dotplot, so CollecTRI uses the cell-46 dotplot settings.
_PLOT_PRESETS = {
    "collectri": {"bar_top": 25, "bar_figsize": (3, 5), "dot_scale": 0.15, "dot_figsize": (4, 4)},
    "progeny": {"bar_top": 25, "bar_figsize": (3, 3), "dot_scale": 0.15, "dot_figsize": (4, 4)},
    "hallmark": {"bar_top": 10, "bar_figsize": (6, 3), "dot_scale": 0.05, "dot_figsize": (6.5, 4)},
}


def _new_output_dir(output_dir: str | None, prefix: str) -> Path:
    base = Path(output_dir).expanduser() if output_dir else _DEFAULT_OUTPUT_BASE
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out = base / f"{prefix}_{stamp}_{uuid.uuid4().hex[:8]}"
    out.mkdir(parents=True, exist_ok=False)
    return out.resolve()


def _read_h5ad(path: str) -> ad.AnnData:
    p = Path(path).expanduser()
    if not p.is_file():
        raise FileNotFoundError(f"data_path does not exist: {path}")
    return ad.read_h5ad(p)


def _require_obs_columns(adata: ad.AnnData, keys: list[str], param: str) -> None:
    missing = [k for k in keys if k not in adata.obs.columns]
    if missing:
        raise ValueError(f"{param} {missing} not found in adata.obs columns {list(adata.obs.columns)}")


def _read_de_results(path: str, columns: list[str]) -> pd.DataFrame:
    p = Path(path).expanduser()
    if not p.is_file():
        raise FileNotFoundError(f"de_results_path does not exist: {path}")
    de = pd.read_csv(p, index_col=0)
    if de.index.hasnans or not de.index.is_unique:
        raise ValueError("DE results must have unique, non-missing gene names in the first (index) column")
    missing = [c for c in columns if c not in de.columns]
    if missing:
        raise ValueError(f"DE results lack column(s) {missing}; available columns: {list(de.columns)}")
    for c in columns:
        if not pd.api.types.is_numeric_dtype(de[c]):
            raise ValueError(f"DE results column {c!r} is not numeric")
    return de


def _load_network(network: str, network_path: str | None) -> tuple[pd.DataFrame, str]:
    if network_path is None:
        net = _NETWORK_LOADERS[network](organism="human")
        origin = f"downloaded with decoupler.op.{network}(organism='human')"
    else:
        p = Path(network_path).expanduser()
        if not p.is_file():
            raise FileNotFoundError(f"network_path does not exist: {network_path}")
        name = p.name.lower()
        if name.endswith(".parquet"):
            net = pd.read_parquet(p)
        elif name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
            net = pd.read_csv(p, sep="\t")
        else:
            net = pd.read_csv(p)
        origin = str(p.resolve())
    if not {"source", "target"}.issubset(net.columns):
        raise ValueError(f"network table needs 'source' and 'target' columns; found {list(net.columns)}")
    return net, origin


def _num(x) -> float | None:
    if x is None:
        return None
    x = float(x)
    return x if np.isfinite(x) else None


def _top_sources(scores: pd.Series, padj: pd.Series, n: int = 5) -> dict:
    order = scores.sort_values()

    def rows(idx):
        return [{"source": str(s), "score": _num(scores[s]), "padj": _num(padj[s])} for s in idx]

    return {
        "positive": rows([s for s in order.index[::-1][:n] if scores[s] > 0]),
        "negative": rows([s for s in order.index[:n] if scores[s] < 0]),
    }


@bulk_enrichment_mcp.tool()
def decoupler_filter_genes_by_expression(
    data_path: Annotated[str, "AnnData .h5ad of raw counts with samples as observations and genes as variables"],
    group_key: Annotated[str, "adata.obs column with the experimental groups; its smallest group sets the minimum sample size"],
    min_count: Annotated[int, "Minimum count a gene needs in the minimum number of samples (applied as a CPM cutoff)"] = 10,
    min_total_count: Annotated[int, "Minimum total count of a gene across all samples"] = 15,
    large_n: Annotated[int, "Number of samples per group considered large"] = 10,
    min_prop: Annotated[float, "Minimum proportion of samples in the smallest group expressing the gene when groups are large"] = 0.7,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Remove lowly expressed genes from a bulk RNA-seq count matrix with edgeR filterByExpr logic.
    Input is a raw-count .h5ad with a group column; output is the filtered .h5ad and a threshold diagnostic plot.
    """
    adata = _read_h5ad(data_path)
    _require_obs_columns(adata, [group_key], "group_key")
    n_before = int(adata.n_vars)
    out = _new_output_dir(output_dir, "filter_by_expr")
    plot_path = out / "filter_by_expr_thresholds.png"
    h5ad_path = out / "filtered_counts.h5ad"

    # Diagnostic plot (cell 11) on the unfiltered matrix, then the filter itself (cell 13).
    with _PLOT_LOCK:
        fig = dc.pl.filter_by_expr(
            adata=adata,
            group=group_key,
            min_count=min_count,
            min_total_count=min_total_count,
            large_n=large_n,
            min_prop=min_prop,
            return_fig=True,
            save=str(plot_path),
        )
        plt.close(fig)
    dc.pp.filter_by_expr(
        adata=adata,
        group=group_key,
        min_count=min_count,
        min_total_count=min_total_count,
        large_n=large_n,
        min_prop=min_prop,
    )
    adata.write_h5ad(h5ad_path)

    n_after = int(adata.n_vars)
    group_sizes = {str(k): int(v) for k, v in adata.obs[group_key].value_counts().sort_index().items()}
    message = f"Kept {n_after} of {n_before} genes ({n_before - n_after} removed) across {adata.n_obs} samples."
    if n_after == 0:
        message += " No gene passed the thresholds; relax them after inspecting the diagnostic plot."
    return {
        "message": message,
        "reference": REFERENCE,
        "artifacts": [
            {"description": "filtered raw-count AnnData (obs kept, genes subset)", "path": str(h5ad_path)},
            {"description": "log10 total counts vs number of samples with threshold lines", "path": str(plot_path)},
        ],
        "n_samples": int(adata.n_obs),
        "n_genes_before": n_before,
        "n_genes_after": n_after,
        "n_genes_removed": n_before - n_after,
        "group_sizes": group_sizes,
        "parameters": {
            "group_key": group_key,
            "min_count": min_count,
            "min_total_count": min_total_count,
            "large_n": large_n,
            "min_prop": min_prop,
        },
    }


@bulk_enrichment_mcp.tool()
def decoupler_rank_pcs_by_metadata(
    data_path: Annotated[str, "AnnData .h5ad of bulk counts (samples x genes) with sample metadata, e.g. the filtered .h5ad from decoupler_filter_genes_by_expression"],
    color_keys: Annotated[list[str] | None, "adata.obs columns used to colour the PCA scatter; defaults to the tested metadata columns"] = None,
    obs_keys: Annotated[list[str] | None, "adata.obs columns tested for association with each PC; defaults to all obs columns"] = None,
    nvar: Annotated[int, "Number of leading PCs shown in the association heatmap"] = 5,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Explore sample variability: PCA of normalised counts and association of each PC with sample metadata.
    Input is a count .h5ad with metadata columns; output is the PC-metadata association table, PC scores and PCA plots.
    """
    adata = _read_h5ad(data_path)
    if adata.obs.shape[1] == 0:
        raise ValueError("adata.obs has no metadata columns to test against the PCs")
    if obs_keys is not None:
        _require_obs_columns(adata, obs_keys, "obs_keys")
    tested = list(obs_keys) if obs_keys is not None else [str(c) for c in adata.obs.columns]
    color_keys = list(color_keys) if color_keys is not None else tested
    _require_obs_columns(adata, color_keys, "color_keys")

    # Work on a writeable floating copy of the counts (the tutorial's counts are float64).
    if np.issubdtype(adata.X.dtype, np.floating):
        adata.X = adata.X.copy()
    else:
        adata.X = adata.X.astype(np.float64)

    # Normalize, scale and compute PCA (cell 15); test PC-metadata associations (cell 17).
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)
    dc.tl.rankby_obsm(adata, key="X_pca", obs_keys=obs_keys)

    out = _new_output_dir(output_dir, "rank_pcs")
    assoc = adata.uns["rank_obsm"]
    assoc_path = out / "pc_metadata_associations.csv"
    assoc.to_csv(assoc_path, index=False)
    n_pcs = adata.obsm["X_pca"].shape[1]
    pc_names = [f"PC{i}" for i in range(1, n_pcs + 1)]
    var_path = out / "pca_variance_ratio.csv"
    pd.DataFrame(
        {
            "PC": pc_names,
            "variance": adata.uns["pca"]["variance"],
            "variance_ratio": adata.uns["pca"]["variance_ratio"],
        }
    ).to_csv(var_path, index=False)
    scores_path = out / "pca_scores.csv"
    pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names, columns=pc_names).to_csv(scores_path)

    var_plot = out / "pca_variance_ratio.png"
    obsm_plot = out / "pc_metadata_association_heatmap.png"
    pca_plot = out / "pca_scatter.png"
    # Plots of cells 19 and 21.
    with _PLOT_LOCK:
        sc.pl.pca_variance_ratio(adata, show=False)
        fig = plt.gcf()
        fig.savefig(var_plot, bbox_inches="tight")
        plt.close(fig)
        fig = dc.pl.obsm(
            adata=adata,
            return_fig=True,
            nvar=nvar,
            titles=["PC scores", "Adjusted p-values"],
            figsize=(6, 3),
            save=str(obsm_plot),
        )
        plt.close(fig)
        fig = sc.pl.pca(
            adata,
            color=color_keys,
            ncols=1,
            size=300,
            frameon=True,
            show=False,
            return_fig=True,
        )
        fig.savefig(pca_plot, bbox_inches="tight")
        plt.close(fig)

    sig = assoc[assoc["padj"] < 0.05]
    significant = [
        {"obsm": str(r.obsm), "obs": str(r.obs), "stat": _num(r.stat), "padj": _num(r.padj)}
        for r in sig.head(20).itertuples()
    ]
    return {
        "message": f"PCA with {n_pcs} PCs on {adata.n_obs} samples x {adata.n_vars} genes; "
        f"{len(sig)} of {len(assoc)} PC-metadata tests have padj < 0.05.",
        "reference": REFERENCE,
        "artifacts": [
            {"description": "PC x metadata association table (obsm, obs, stat, pval, padj; ANOVA or Spearman, BH)", "path": str(assoc_path)},
            {"description": "PCA explained variance and variance ratio per PC", "path": str(var_path)},
            {"description": "PC scores per sample", "path": str(scores_path)},
            {"description": "PCA variance ratio plot", "path": str(var_plot)},
            {"description": "PC scores and -log10 adjusted p-value heatmap", "path": str(obsm_plot)},
            {"description": "PCA scatter coloured by metadata", "path": str(pca_plot)},
        ],
        "n_samples": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_pcs": int(n_pcs),
        "variance_ratio": [_num(v) for v in adata.uns["pca"]["variance_ratio"][:10]],
        "tested_obs_keys": tested,
        "color_keys": color_keys,
        "n_tests": int(len(assoc)),
        "significant_associations": significant,
    }


@bulk_enrichment_mcp.tool()
def pydeseq2_test_differential_expression(
    data_path: Annotated[str, "AnnData .h5ad of raw integer counts (samples x genes) with sample metadata, e.g. the filtered .h5ad from decoupler_filter_genes_by_expression"],
    condition_key: Annotated[str, "adata.obs column holding the condition labels to contrast"],
    group: Annotated[str, "Condition level tested against the reference (numerator of the fold change)"],
    reference: Annotated[str, "Condition level used as the reference (denominator of the fold change)"],
    design: Annotated[str | None, "Formulaic design formula such as '~batch + condition'; defaults to '~<condition_key>'"] = None,
    n_cpus: Annotated[int, "CPU workers for PyDESeq2 inference; results do not depend on it"] = 8,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Test differential expression between two conditions with the PyDESeq2 Wald test.
    Input is a raw-count .h5ad with a condition column; output is the gene-indexed results CSV and a volcano plot.
    """
    adata = _read_h5ad(data_path)
    _require_obs_columns(adata, [condition_key], "condition_key")
    labels = adata.obs[condition_key]
    if pd.api.types.is_numeric_dtype(labels) or pd.api.types.is_bool_dtype(labels):
        raise ValueError(
            f"condition_key {condition_key!r} is numeric/boolean; PyDESeq2 treats it as continuous. "
            "Store the condition labels as strings."
        )
    levels = sorted(str(v) for v in labels.dropna().unique())
    if group == reference:
        raise ValueError("group and reference must be different condition levels")
    missing = [g for g in (group, reference) if g not in levels]
    if missing:
        raise ValueError(f"level(s) {missing} not found in adata.obs[{condition_key!r}]; available levels: {levels}")
    if n_cpus < 1:
        raise ValueError("n_cpus must be >= 1")

    # PyDESeq2 precondition: a dense matrix of non-negative integer counts.
    if sps.issparse(adata.X):
        adata.X = adata.X.toarray()
    counts = np.asarray(adata.X)
    if not np.issubdtype(counts.dtype, np.number) or not np.isfinite(counts).all():
        raise ValueError("counts must be finite numbers (raw counts expected)")
    if (counts < 0).any() or (counts % 1 != 0).any():
        raise ValueError("counts must be non-negative integers (raw counts); normalised or log values are not valid input")

    design = design if design is not None else f"~{condition_key}"
    # DE cells 23-24. PyDESeq2 shares and mutates the AnnData it receives, so it gets a copy;
    # quiet=True only silences its stdout prints (stdout carries the MCP protocol).
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        adata=adata.copy(),
        design=design,
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    if condition_key not in dds.variables:
        raise ValueError(f"design {design!r} does not include condition_key {condition_key!r}")
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=[condition_key, group, reference], inference=inference, quiet=True)
    stat_res.summary()
    results_df = stat_res.results_df

    out = _new_output_dir(output_dir, "deseq2")
    csv_path = out / "deseq2_results.csv"
    results_df.to_csv(csv_path)
    volcano_path = out / "deseq2_volcano.png"
    # Volcano plot (cell 26).
    with _PLOT_LOCK:
        fig = dc.pl.volcano(
            data=results_df,
            x="log2FoldChange",
            y="padj",
            top=10,
            return_fig=True,
            save=str(volcano_path),
        )
        plt.close(fig)

    sig = results_df[results_df["padj"] < 0.05]

    def top(df):
        return [
            {"gene": str(r.Index), "log2FoldChange": _num(r.log2FoldChange), "stat": _num(r.stat), "padj": _num(r.padj)}
            for r in df.itertuples()
        ]

    up = sig[sig["log2FoldChange"] > 0].sort_values("padj").head(5)
    down = sig[sig["log2FoldChange"] < 0].sort_values("padj").head(5)
    return {
        "message": f"PyDESeq2 Wald test {condition_key}: {group} vs {reference} on {adata.n_obs} samples x "
        f"{results_df.shape[0]} genes; {len(sig)} genes with padj < 0.05.",
        "reference": REFERENCE,
        "artifacts": [
            {"description": "PyDESeq2 results_df (gene index; baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)", "path": str(csv_path)},
            {"description": "volcano plot of log2FoldChange vs -log10 padj", "path": str(volcano_path)},
        ],
        "design": design,
        "contrast": [condition_key, group, reference],
        "n_samples": int(adata.n_obs),
        "samples_per_level": {str(k): int(v) for k, v in labels.astype(str).value_counts().sort_index().items()},
        "n_genes": int(results_df.shape[0]),
        "n_missing": {c: int(results_df[c].isna().sum()) for c in results_df.columns},
        "n_padj_lt_0.05": int(len(sig)),
        "n_padj_lt_0.05_up": int((sig["log2FoldChange"] > 0).sum()),
        "n_padj_lt_0.05_down": int((sig["log2FoldChange"] < 0).sum()),
        "top_up": top(up),
        "top_down": top(down),
        "pydeseq2_version": version("pydeseq2"),
        "n_cpus": n_cpus,
    }


@bulk_enrichment_mcp.tool()
def decoupler_score_activities_ulm(
    de_results_path: Annotated[str, "CSV of gene-level DE results with genes in the first column, as written by pydeseq2_test_differential_expression"],
    network: Annotated[NetworkName, "Human prior-knowledge network: CollecTRI TF regulons, PROGENy pathway footprints or MSigDB Hallmark gene sets"],
    stat_column: Annotated[str, "DE column used as the gene-level statistic; the tutorial recommends the Wald 'stat' over log2FoldChange"] = "stat",
    contrast_name: Annotated[str | None, "Row label of the contrast in output tables and plots; defaults to the stat column name"] = None,
    padj_threshold: Annotated[float, "Adjusted p-value cutoff defining significant sources"] = 0.05,
    barplot_top: Annotated[int | None, "Number of top sources in the bar plot; defaults to the tutorial value (25 CollecTRI/PROGENy, 10 Hallmark)"] = None,
    dotplot_top: Annotated[int, "Number of top significant sources in the dot plot"] = 10,
    network_path: Annotated[str | None, "Saved network table (CSV, TSV or Parquet with source, target and optional weight) used instead of downloading"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Infer TF, pathway or Hallmark activities for a contrast by running ulm on the gene-level DE statistic.
    Input is the DE results CSV; output is ulm scores and adjusted p-values for all sources, the significant subset, the network used and bar/dot plots.
    """
    if not 0 < padj_threshold <= 1:
        raise ValueError("padj_threshold must be in (0, 1]")
    de = _read_de_results(de_results_path, [stat_column])
    n_nonfinite = int((~np.isfinite(de[stat_column].to_numpy(dtype=float))).sum())
    if n_nonfinite:
        raise ValueError(
            f"{n_nonfinite} genes have a non-finite {stat_column!r}; ulm requires finite values "
            "(filter lowly expressed genes before DE, e.g. with decoupler_filter_genes_by_expression)"
        )
    # Wide 1 x genes matrix of the statistic (cell 28).
    data = de[[stat_column]].T
    if contrast_name is not None:
        data = data.rename(index={stat_column: contrast_name})
    row = str(data.index[0])
    net, origin = _load_network(network, network_path)

    # Scoring and significance filter (cells 32 / 42 / 54).
    acts, padj = dc.mt.ulm(data=data, net=net)
    msk = (padj.T < padj_threshold).iloc[:, 0]
    sig_acts = acts.loc[:, msk]

    out = _new_output_dir(output_dir, f"ulm_{network}")
    net_path = out / f"{network}_network_used.csv"
    net.to_csv(net_path, index=False)
    scores_path = out / f"{network}_ulm_scores.csv"
    acts.to_csv(scores_path)
    padj_path = out / f"{network}_ulm_padj.csv"
    padj.to_csv(padj_path)
    sig_path = out / f"{network}_ulm_scores_significant.csv"
    sig_acts.to_csv(sig_path)
    artifacts = [
        {"description": f"ulm scores for all sources (row {row!r}, one column per source)", "path": str(scores_path)},
        {"description": "BH-adjusted ulm p-values for all sources (same layout)", "path": str(padj_path)},
        {"description": f"ulm scores of sources with padj < {padj_threshold}", "path": str(sig_path)},
        {"description": "network table used for scoring (provenance)", "path": str(net_path)},
    ]

    preset = _PLOT_PRESETS[network]
    n_sig = int(sig_acts.shape[1])
    if n_sig > 0:
        bar_path = out / f"{network}_ulm_barplot.png"
        dot_path = out / f"{network}_ulm_dotplot.png"
        # Dot plot table (cells 46 / 58): significant scores with -log10 padj clipped at 2.22e-16.
        df = sig_acts.melt(value_name="score").merge(
            padj.melt(value_name="padj")
            .assign(padj=lambda x: x["padj"].clip(2.22e-16, 1))
            .assign(padj=lambda x: -np.log10(x["padj"]))
        )
        with _PLOT_LOCK:
            # Bar plot (cells 34 / 44 / 56).
            fig = dc.pl.barplot(
                data=sig_acts,
                name=row,
                top=barplot_top if barplot_top is not None else preset["bar_top"],
                figsize=preset["bar_figsize"],
                return_fig=True,
                save=str(bar_path),
            )
            plt.close(fig)
            fig = dc.pl.dotplot(
                df=df,
                x="score",
                y="variable",
                s="padj",
                c="score",
                top=dotplot_top,
                scale=preset["dot_scale"],
                figsize=preset["dot_figsize"],
                return_fig=True,
                save=str(dot_path),
            )
            plt.close(fig)
        artifacts += [
            {"description": "bar plot of the top significant source scores", "path": str(bar_path)},
            {"description": "dot plot of significant scores sized by -log10 padj", "path": str(dot_path)},
        ]

    scores = acts.iloc[0]
    padjs = padj.iloc[0]
    sig_scores = sig_acts.iloc[0]
    message = f"ulm with {network} ({origin}): {acts.shape[1]} sources scored, {n_sig} with padj < {padj_threshold}."
    if n_sig == 0:
        message += " No significant source, so no bar/dot plot was drawn."
    return {
        "message": message,
        "reference": REFERENCE,
        "artifacts": artifacts,
        "network": network,
        "network_origin": origin,
        "network_rows": int(net.shape[0]),
        "network_n_sources": int(net["source"].nunique()),
        "stat_column": stat_column,
        "contrast_row": row,
        "n_genes": int(data.shape[1]),
        "n_sources_scored": int(acts.shape[1]),
        "n_significant": n_sig,
        "n_significant_positive": int((sig_scores > 0).sum()),
        "n_significant_negative": int((sig_scores < 0).sum()),
        "top_significant_sources": _top_sources(sig_scores, padjs[sig_scores.index]),
        "top_sources_all": _top_sources(scores, padjs),
    }


@bulk_enrichment_mcp.tool()
def decoupler_inspect_source_targets(
    de_results_path: Annotated[str, "CSV of gene-level DE results (gene index; log2FoldChange, pvalue and the stat column), as written by pydeseq2_test_differential_expression"],
    network: Annotated[NetworkName, "Human prior-knowledge network containing the source: CollecTRI, PROGENy or MSigDB Hallmark"],
    source: Annotated[str, "Source to inspect (TF, pathway or gene set), e.g. a top hit of decoupler_score_activities_ulm"],
    stat_column: Annotated[str, "DE column used as the ranking statistic for the leading edge and target plots"] = "stat",
    scores_path: Annotated[str | None, "Wide ulm scores CSV from decoupler_score_activities_ulm; needed for the network plot, which is skipped without it"] = None,
    network_sources: Annotated[list[str] | None, "Sources drawn in the network plot; defaults to [source]"] = None,
    network_path: Annotated[str | None, "Saved network table (CSV, TSV or Parquet with source, target and optional weight) used instead of downloading"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Explain one source's activity through its target genes: target volcano, weight-vs-statistic scatter, leading edge and network plot.
    Input is the DE results CSV (plus optional ulm scores CSV); output is the leading-edge gene lists, the network used and the plots.
    """
    de = _read_de_results(de_results_path, ["log2FoldChange", "pvalue", stat_column])
    net, origin = _load_network(network, network_path)
    src_net = net[net["source"].astype(str) == source]
    if src_net.empty:
        raise ValueError(f"source {source!r} not found in the {network} network")
    n_measured = int(src_net["target"].astype(str).isin(de.index.astype(str)).sum())
    if n_measured == 0:
        raise ValueError(f"none of the {len(src_net)} targets of {source!r} is present in the DE results")
    weighted = "weight" in net.columns

    score = None
    if scores_path is not None:
        p = Path(scores_path).expanduser()
        if not p.is_file():
            raise FileNotFoundError(f"scores_path does not exist: {scores_path}")
        score = pd.read_csv(p, index_col=0)
        if score.shape[0] != 1:
            raise ValueError(f"scores table must have exactly one contrast row; found {score.shape[0]}")
        sources = list(network_sources) if network_sources is not None else [source]
        absent = [s for s in sources if s not in score.columns]
        if absent:
            raise ValueError(f"source(s) {absent} not in the scores table; pass the full scores CSV or other sources")

    out = _new_output_dir(output_dir, f"inspect_{network}")
    safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in source)
    net_path = out / f"{network}_network_used.csv"
    net.to_csv(net_path, index=False)
    artifacts = [{"description": "network table used (provenance)", "path": str(net_path)}]
    skipped = []
    leading_edges = {}

    with _PLOT_LOCK:
        # Volcano of the source's targets (cell 38).
        path = out / f"{safe}_targets_volcano.png"
        fig = dc.pl.volcano(
            data=de, x="log2FoldChange", y="pvalue", net=net, name=source, top=10, return_fig=True, save=str(path)
        )
        plt.close(fig)
        artifacts.append({"description": f"volcano of {source} targets (log2FoldChange vs -log10 pvalue)", "path": str(path)})

        # Weight vs statistic scatter (cell 48), only for networks with interaction weights.
        if weighted:
            path = out / f"{safe}_source_targets.png"
            fig = dc.pl.source_targets(
                data=de, x="weight", y=stat_column, net=net, name=source, top=15, figsize=(4, 4),
                return_fig=True, save=str(path),
            )
            plt.close(fig)
            artifacts.append({"description": f"{source} target weights vs {stat_column}", "path": str(path)})
        else:
            skipped.append("source_targets plot: the network has no weight column")

        # Leading edge (cells 50 and 60): split by weight sign for weighted networks, as for PROGENy.
        if weighted:
            parts = [("positive_weights", net[net["weight"] > 0]), ("negative_weights", net[net["weight"] < 0])]
        else:
            parts = [("all_targets", net)]
        for label, sub in parts:
            sub_targets = sub.loc[sub["source"].astype(str) == source, "target"].astype(str)
            if not sub_targets.isin(de.index.astype(str)).any():
                skipped.append(f"leading edge {label}: {source} has no measured targets in this component")
                continue
            path = out / f"{safe}_leading_edge_{label}.png"
            fig, le = dc.pl.leading_edge(de, stat=stat_column, net=sub, name=source, return_fig=True, save=str(path))
            plt.close(fig)
            genes_path = out / f"{safe}_leading_edge_{label}.txt"
            genes_path.write_text("".join(f"{g}\n" for g in le))
            leading_edges[label] = {"n_genes": int(len(le)), "first_genes": [str(g) for g in le[:10]], "path": str(genes_path)}
            artifacts.append({"description": f"leading-edge plot ({label})", "path": str(path)})
            artifacts.append({"description": f"leading-edge genes ({label}), one per line in rank order", "path": str(genes_path)})

        # Regulator-target network (cell 36) coloured by ulm score and gene-level statistic.
        if score is not None:
            data = de[[stat_column]].T
            data.index = score.index
            path = out / f"{safe}_network.png"
            fig = dc.pl.network(
                net=net,
                data=data,
                score=score,
                sources=sources,
                targets=10,
                figsize=(5, 5),
                vcenter=True,
                by_abs=True,
                size_node=10,
                t_label="Gene-level\nstatistic",
                return_fig=True,
                save=str(path),
            )
            plt.close(fig)
            artifacts.append({"description": f"network of {sources} and their top targets", "path": str(path)})
        else:
            skipped.append("network plot: no scores_path from decoupler_score_activities_ulm was given")

    result = {
        "message": f"Inspected {source} in {network} ({origin}): {len(src_net)} targets, {n_measured} measured; "
        f"leading edges: " + ", ".join(f"{k} {v['n_genes']} genes" for k, v in leading_edges.items()) + ".",
        "reference": REFERENCE,
        "artifacts": artifacts,
        "source": source,
        "network": network,
        "network_origin": origin,
        "network_rows": int(net.shape[0]),
        "weighted_network": weighted,
        "n_targets_in_network": int(len(src_net)),
        "n_targets_measured": n_measured,
        "leading_edges": leading_edges,
        "skipped": skipped,
    }
    if score is not None:
        result["network_plot_sources"] = sources
        if source in score.columns:
            result["source_score"] = _num(score[source].iloc[0])
    return result
