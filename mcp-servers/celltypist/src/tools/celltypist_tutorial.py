"""Tools extracted from Teichlab/celltypist docs/notebook/celltypist_tutorial.ipynb.

CellTypist assigns cell type labels to log1p-normalised (10,000 counts/cell) scRNA-seq data with
logistic-regression models, optionally refining per-cell predictions by majority voting within
over-clustered subclusters. Every tool below calls the installed celltypist package directly.

Built-in models are cached under ``$CELLTYPIST_FOLDER/data/models`` (default ``~/.celltypist``);
celltypist fixes that location when it is first imported, so set the variable before the server starts.
"""
import os
import tempfile
import time
from pathlib import Path
from typing import Annotated, Literal

import matplotlib

matplotlib.use("Agg")  # headless server: never open a display
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from matplotlib import pyplot as plt  # noqa: E402

import celltypist  # noqa: E402
from celltypist import models  # noqa: E402
from fastmcp import FastMCP  # noqa: E402

REFERENCE = (
    "https://github.com/Teichlab/celltypist/blob/fe357564a6625d3b1732a022fd39f18e55696e80/"
    "docs/notebook/celltypist_tutorial.ipynb"
)
DEFAULT_OUTPUT_BASE = Path(tempfile.gettempdir()) / "celltypist_mcp_outputs"

celltypist_tutorial_mcp = FastMCP(name="celltypist_tutorial")


def _fresh_output_dir(output_dir: str | None, stem: str) -> Path:
    """Create a new, uniquely named directory for one tool invocation."""
    base = Path(output_dir).expanduser().resolve() if output_dir else DEFAULT_OUTPUT_BASE
    base.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=f"{stem}_{time.strftime('%Y%m%d-%H%M%S')}_", dir=base)).resolve()


def _h5ad_path(path: str) -> Path:
    """Resolve a user AnnData path and check it is an existing .h5ad file."""
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"No such file: {p}")
    if p.suffix != ".h5ad":
        raise ValueError(f"Expected an AnnData .h5ad file, got: {p.name}")
    return p


def _resolve_model(model: str) -> tuple[str, str]:
    """Return (argument for celltypist, model file path) for a built-in model name or a .pkl path.

    celltypist's Model.load treats a string without '/' as a built-in model name and, when the model
    cache holds no .pkl at all, get_all_models() downloads every model in the index. Fetching only the
    requested model first (models.download_models(model=name), as in tutorial cell 16) prevents that.
    """
    if "/" in model:
        p = Path(model).expanduser().resolve()
        if not p.is_file():
            raise FileNotFoundError(f"No such model file: {p}")
        return str(p), str(p)
    model_file = models.get_model_path(model)
    if not os.path.isfile(model_file):
        models.download_models(model=model)  # raises ValueError if the name is not in the model index
        if not os.path.isfile(model_file):
            raise RuntimeError(f"Built-in model {model!r} could not be downloaded to {models.models_path}")
    return model, model_file


def _frequency_records(freq: pd.DataFrame) -> list[dict]:
    return [{"celltype": str(r.celltype), "counts": int(r.counts)} for r in freq.itertuples(index=False)]


@celltypist_tutorial_mcp.tool()
def celltypist_list_models(
    refresh_index: Annotated[bool, "Re-download the CellTypist model index (models.json) before listing"] = False,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """List the built-in CellTypist models and what each one represents.
    Cached (or refreshed) CellTypist model index -> CSV of model filenames and descriptions."""
    if refresh_index:
        models.get_models_index(force_update=True)
    description = models.models_description()
    out = _fresh_output_dir(output_dir, "celltypist_list_models")
    csv_path = out / "models_description.csv"
    description.to_csv(csv_path, index=False)
    return {
        "message": f"{len(description)} built-in CellTypist models listed",
        "reference": REFERENCE,
        "artifacts": [{"description": "model filenames and descriptions", "path": str(csv_path)}],
        "output_dir": str(out),
        "n_models": int(len(description)),
        "default_model": models.get_default_model(),
        "index_path": str(Path(models.get_model_path("models.json")).resolve()),
        "models": description.to_dict("records"),
    }


@celltypist_tutorial_mcp.tool()
def celltypist_annotate_cells(
    adata_path: Annotated[str, "AnnData .h5ad with log1p-normalised (10,000 counts/cell) expression in .X or .raw.X and gene symbols as var_names"],
    model: Annotated[str, "Built-in model filename (see celltypist_list_models) or path to a custom .pkl model"] = "Immune_All_Low.pkl",
    mode: Annotated[Literal["best match", "prob match"], "Single most likely label per cell, or multi-label assignment by probability threshold"] = "best match",
    p_thres: Annotated[float, "Probability threshold for 'prob match' mode; ignored in 'best match' mode"] = 0.5,
    majority_voting: Annotated[bool, "Refine per-cell labels by majority voting within over-clustered subclusters"] = False,
    over_clustering_key: Annotated[str | None, "Existing obs column with user over-clustering for majority voting; default builds a Leiden over-clustering"] = None,
    min_prop: Annotated[float, "Minimum proportion of the dominant label for naming a subcluster, else 'Heterogeneous' (majority voting only)"] = 0.0,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Annotate cells of a normalised scRNA-seq AnnData with a built-in or custom CellTypist model.
    AnnData .h5ad + model -> annotated.h5ad plus predicted-label, decision-score and probability CSV tables."""
    data_file = _h5ad_path(adata_path)
    model_arg, model_file = _resolve_model(model)
    adata = sc.read(str(data_file))
    if majority_voting and over_clustering_key is not None and over_clustering_key not in adata.obs:
        raise KeyError(f"over_clustering_key {over_clustering_key!r} is not an obs column of {data_file.name}; "
                       f"available columns: {list(adata.obs.columns)}")

    predictions = celltypist.annotate(
        adata,
        model=model_arg,
        mode=mode,
        p_thres=p_thres,
        majority_voting=majority_voting,
        over_clustering=over_clustering_key,
        min_prop=min_prop,
    )

    out = _fresh_output_dir(output_dir, "celltypist_annotate_cells")
    predictions.to_table(str(out))  # predicted_labels.csv, decision_matrix.csv, probability_matrix.csv
    # Probabilities stay in the CSV tables: several model cell-type names contain '/', which h5ad keys reject.
    annotated = predictions.to_adata(insert_labels=True, insert_conf=True)
    h5ad_path = out / "annotated.h5ad"
    annotated.write_h5ad(h5ad_path)

    labels = predictions.predicted_labels
    mv_done = "majority_voting" in labels.columns
    frequencies = {by: _frequency_records(predictions.summary_frequency(by=by))
                   for by in ("predicted_labels", "majority_voting") if by in labels.columns}
    message = f"Annotated {predictions.cell_count} cells with {Path(model_file).name} ({mode})"
    if majority_voting and not mv_done:
        message += "; majority voting was skipped by CellTypist because the dataset has 50 or fewer cells"
    result = {
        "message": message,
        "reference": REFERENCE,
        "artifacts": [
            {"description": "per-cell predicted labels (and over_clustering, majority_voting when run)", "path": str(out / "predicted_labels.csv")},
            {"description": "cell-by-cell-type decision scores", "path": str(out / "decision_matrix.csv")},
            {"description": "cell-by-cell-type probabilities", "path": str(out / "probability_matrix.csv")},
            {"description": "input AnnData with predicted labels and conf_score in obs", "path": str(h5ad_path)},
        ],
        "output_dir": str(out),
        "n_cells": int(predictions.cell_count),
        "n_model_cell_types": int(predictions.probability_matrix.shape[1]),
        "model": model_arg,
        "model_file": model_file,
        "mode": mode,
        "majority_voting_requested": majority_voting,
        "majority_voting_performed": mv_done,
        "label_frequencies": frequencies,
    }
    if mode == "prob match":
        result["p_thres"] = p_thres
    if mv_done:
        result["over_clustering_source"] = over_clustering_key or "celltypist Leiden over-clustering"
        result["n_over_clusters"] = int(labels["over_clustering"].nunique())
        result["min_prop"] = min_prop
    return result


@celltypist_tutorial_mcp.tool()
def celltypist_train_model(
    adata_path: Annotated[str, "Reference AnnData .h5ad with log1p-normalised (10,000 counts/cell) expression in .X (or .raw.X if .X is scaled)"],
    label_key: Annotated[str, "obs column holding the cell type labels to learn"],
    feature_selection: Annotated[bool, "Two-pass training: select top genes per cell type with SGD, then retrain on them"] = False,
    top_genes: Annotated[int, "Genes kept per cell type by absolute coefficient when feature_selection is on"] = 300,
    n_jobs: Annotated[int | None, "CPUs for training (-1 for all; default one); seeded results are bit-identical only at the same n_jobs"] = None,
    seed: Annotated[int | None, "random_state passed to the SGD and logistic-regression classifiers for reproducible training"] = None,
    model_name: Annotated[str, "File name (without directory) for the saved .pkl model"] = "celltypist_model",
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Train a custom CellTypist logistic-regression model from a labelled reference AnnData.
    AnnData .h5ad + label column -> saved .pkl model usable by celltypist_annotate_cells and celltypist_inspect_model."""
    data_file = _h5ad_path(adata_path)
    if Path(model_name).name != model_name or not model_name:
        raise ValueError(f"model_name must be a plain file name without directories, got {model_name!r}")
    if top_genes < 1:
        raise ValueError(f"top_genes must be a positive integer, got {top_genes}")
    adata = sc.read(str(data_file))
    if label_key not in adata.obs:
        raise KeyError(f"label_key {label_key!r} is not an obs column of {data_file.name}; "
                       f"available columns: {list(adata.obs.columns)}")

    extra = {} if seed is None else {"random_state": seed}  # forwarded by train(**kwargs) to the sklearn classifiers
    new_model = celltypist.train(adata, labels=label_key, n_jobs=n_jobs,
                                 feature_selection=feature_selection, top_genes=top_genes, **extra)

    out = _fresh_output_dir(output_dir, "celltypist_train_model")
    stem = model_name[:-4] if model_name.endswith(".pkl") else model_name
    model_path = out / f"{stem}.pkl"
    new_model.write(str(model_path))
    return {
        "message": f"Trained a CellTypist model with {len(new_model.cell_types)} cell types and {len(new_model.features)} features",
        "reference": REFERENCE,
        "artifacts": [{"description": "trained CellTypist model (.pkl)", "path": str(model_path)}],
        "output_dir": str(out),
        "model_path": str(model_path),
        "n_cells": int(adata.n_obs),
        "n_cell_types": int(len(new_model.cell_types)),
        "cell_types": [str(x) for x in new_model.cell_types],
        "n_features": int(len(new_model.features)),
        "label_key": label_key,
        "feature_selection": feature_selection,
        "top_genes": top_genes,
        "n_jobs": n_jobs,
        "seed": seed,
        "description": {str(k): str(v) for k, v in new_model.description.items()},
    }


@celltypist_tutorial_mcp.tool()
def celltypist_inspect_model(
    model: Annotated[str, "Built-in model filename (see celltypist_list_models) or path to a custom .pkl model"],
    cell_types: Annotated[list[str] | None, "Model cell types to extract driving genes for; default is every cell type in the model"] = None,
    top_n: Annotated[int, "Number of top driving genes per cell type"] = 10,
    only_positive: Annotated[bool, "Rank by positive coefficients only; False ranks by absolute coefficient"] = True,
    adata_path: Annotated[str | None, "Optional AnnData .h5ad in which to plot expression of the driving genes"] = None,
    groupby_key: Annotated[str | None, "obs column used to group cells in the violin plots; required with adata_path"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Report a CellTypist model's metadata and the top driving genes of its cell types, optionally with violin plots.
    Model name or .pkl (+ optional AnnData .h5ad and obs grouping column) -> markers CSV and per-cell-type violin PNGs."""
    if top_n < 1:
        raise ValueError(f"top_n must be a positive integer, got {top_n}")
    adata = None
    if adata_path is not None:
        if groupby_key is None:
            raise ValueError("groupby_key is required when adata_path is given")
        if not cell_types:
            raise ValueError("cell_types must list the model cell types to plot when adata_path is given")
        data_file = _h5ad_path(adata_path)
    model_arg, model_file = _resolve_model(model)
    loaded = models.Model.load(model_arg)
    selected = list(cell_types) if cell_types else [str(x) for x in loaded.cell_types]
    markers = {ct: [str(g) for g in loaded.extract_top_markers(ct, top_n, only_positive)] for ct in selected}

    if adata_path is not None:
        adata = sc.read(str(data_file))
        if groupby_key not in adata.obs:
            raise KeyError(f"groupby_key {groupby_key!r} is not an obs column of {data_file.name}; "
                           f"available columns: {list(adata.obs.columns)}")

    out = _fresh_output_dir(output_dir, "celltypist_inspect_model")
    markers_csv = out / "top_markers.csv"
    pd.DataFrame(
        [{"cell_type": ct, "rank": i + 1, "gene": g} for ct, genes in markers.items() for i, g in enumerate(genes)],
        columns=["cell_type", "rank", "gene"],
    ).to_csv(markers_csv, index=False)
    artifacts = [{"description": "top driving genes per cell type", "path": str(markers_csv)}]

    missing_genes = {}
    if adata is not None:
        # Tutorial cells 67-68: sc.pl.violin(adata, top_genes, groupby=..., rotation=90); scanpy plots .raw when present.
        available = adata.raw.var_names if adata.raw is not None else adata.var_names
        for i, (ct, genes) in enumerate(markers.items(), start=1):
            present = [g for g in genes if g in available]
            absent = [g for g in genes if g not in available]
            if absent:
                missing_genes[ct] = absent
            if not present:
                continue
            sc.pl.violin(adata, present, groupby=groupby_key, rotation=90, show=False)
            png = out / f"violin_{i:02d}_{''.join(c if c.isalnum() else '_' for c in ct)}.png"
            plt.savefig(png, bbox_inches="tight")
            plt.close("all")
            artifacts.append({"description": f"expression of {ct} driving genes by {groupby_key}", "path": str(png)})

    result = {
        "message": f"Model with {len(loaded.cell_types)} cell types and {len(loaded.features)} features; "
                   f"top {top_n} driving genes extracted for {len(markers)} cell types",
        "reference": REFERENCE,
        "artifacts": artifacts,
        "output_dir": str(out),
        "model": model_arg,
        "model_file": model_file,
        "description": {str(k): str(v) for k, v in loaded.description.items()},
        "n_cell_types": int(len(loaded.cell_types)),
        "n_features": int(len(loaded.features)),
        "cell_types": [str(x) for x in loaded.cell_types],
        "top_markers": markers,
    }
    if adata is not None:
        result["groupby_key"] = groupby_key
        result["genes_missing_from_adata"] = missing_genes
    return result
