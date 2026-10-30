"""Tools extracted from celltypist/docs/notebook/celltypist_tutorial_ml.ipynb.

Compares CellTypist predictions with an existing reference label/cluster column,
as in the tutorial: celltypist.dotplot (ml cells 49/51; main tutorial cells 41/60),
pd.crosstab of reference vs prediction (ml cell 47) and per-reference-group
probability boxplots (ml cell 82).
"""
import contextlib
import re
import sys
import tempfile
import time
from pathlib import Path
from typing import Annotated, Literal

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import anndata as ad  # noqa: E402
import pandas as pd  # noqa: E402
import celltypist  # noqa: E402
from celltypist.classifier import AnnotationResult  # noqa: E402
from fastmcp import FastMCP  # noqa: E402

REFERENCE = "https://github.com/Teichlab/celltypist/blob/fe357564a6625d3b1732a022fd39f18e55696e80/docs/notebook/celltypist_tutorial_ml.ipynb"

celltypist_tutorial_ml_mcp = FastMCP(name="celltypist_tutorial_ml")

# Files written by celltypist_annotate_cells: AnnotationResult.to_table(folder) (no prefix, no xlsx)
# plus AnnotationResult.to_adata().write_h5ad(...).
_LABELS_CSV = "predicted_labels.csv"
_DECISION_CSV = "decision_matrix.csv"
_PROBABILITY_CSV = "probability_matrix.csv"
_ANNOTATED_H5AD = "annotated.h5ad"


def _check_same_cells(expected: pd.Index, other: pd.Index, other_name: str) -> None:
    if len(other) != len(expected) or not other.equals(expected):
        detail = f"{len(expected)} vs {len(other)} cells"
        for pos, (a, b) in enumerate(zip(expected, other)):
            if a != b:
                detail += f"; first mismatch at row {pos}: {a!r} vs {b!r}"
                break
        raise ValueError(
            f"Cell IDs in {_LABELS_CSV} and {other_name} do not align ({detail}). "
            "All four files must come from the same celltypist_annotate_cells call."
        )


def _load_annotation_result(annotation_dir: Path) -> AnnotationResult:
    """Rebuild celltypist.classifier.AnnotationResult from a celltypist_annotate_cells output directory."""
    if not annotation_dir.is_dir():
        raise FileNotFoundError(f"annotation_dir does not exist or is not a directory: {annotation_dir}")
    missing = [f for f in (_LABELS_CSV, _DECISION_CSV, _PROBABILITY_CSV, _ANNOTATED_H5AD) if not (annotation_dir / f).is_file()]
    if missing:
        raise FileNotFoundError(f"annotation_dir {annotation_dir} lacks {missing}; expected the output of celltypist_annotate_cells.")

    # to_table writes the cell IDs as the unnamed first column; keep IDs and labels as literal strings.
    labels = pd.read_csv(annotation_dir / _LABELS_CSV, index_col=0, dtype=str, keep_default_na=False)
    if "predicted_labels" not in labels.columns:
        raise ValueError(f"{_LABELS_CSV} has no 'predicted_labels' column (columns: {list(labels.columns)}).")
    # Upstream builds these as pandas 'category' (classifier.py:389, :498); celltypist.dotplot orders rows by
    # pred.cat.categories (plot.py:48). astype('category') reproduces the same sorted categories.
    for col in ("predicted_labels", "majority_voting"):
        if col in labels.columns:
            labels[col] = labels[col].astype("category")
    decision = pd.read_csv(annotation_dir / _DECISION_CSV, index_col=0, converters={0: str})
    probability = pd.read_csv(annotation_dir / _PROBABILITY_CSV, index_col=0, converters={0: str})
    # Only obs is used for the comparison, so the expression matrix stays on disk.
    adata = ad.read_h5ad(annotation_dir / _ANNOTATED_H5AD, backed="r")

    try:
        if not labels.index.is_unique:
            raise ValueError(f"Cell IDs in {_LABELS_CSV} are not unique.")
        _check_same_cells(labels.index, decision.index, _DECISION_CSV)
        _check_same_cells(labels.index, probability.index, _PROBABILITY_CSV)
        _check_same_cells(labels.index, adata.obs_names, f"{_ANNOTATED_H5AD} obs_names")
        if not decision.columns.equals(probability.columns):
            raise ValueError(f"{_DECISION_CSV} and {_PROBABILITY_CSV} list different model cell types.")
    except Exception:
        adata.file.close()
        raise
    return AnnotationResult(labels, decision, probability, adata)


def _safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_")[:60] or "group"


@celltypist_tutorial_ml_mcp.tool()
def celltypist_compare_with_reference(
    annotation_dir: Annotated[str, "Output directory of celltypist_annotate_cells (predicted_labels.csv, decision_matrix.csv, probability_matrix.csv, annotated.h5ad)"],
    reference_key: Annotated[str, "obs column of annotated.h5ad holding the existing reference cell types or clusters to compare against"],
    use_as_prediction: Annotated[Literal["majority_voting", "predicted_labels"], "Prediction to assess: majority-voting labels (needs annotation with majority voting) or raw per-cell labels"] = "majority_voting",
    reference_groups: Annotated[list[str] | None, "Reference groups (values of reference_key) whose crosstab counts are returned and which get probability boxplots; None for none"] = None,
    plot_group_probabilities: Annotated[bool, "Draw a boxplot of every model cell type's probability for each reference group"] = True,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Compare CellTypist predictions with an existing reference label or cluster column via a label-transfer dot plot, crosstab and per-group probability boxplots.
    Input is a celltypist_annotate_cells output directory plus a reference obs column → dot-plot PNG with its fraction/mean-probability CSVs, crosstab CSV and optional boxplot PNGs."""
    with contextlib.redirect_stdout(sys.stderr):
        predictions = _load_annotation_result(Path(annotation_dir).expanduser().resolve())
        try:
            return _compare(predictions, reference_key, use_as_prediction, reference_groups,
                            plot_group_probabilities, output_dir)
        finally:
            predictions.adata.file.close()


def _compare(predictions: AnnotationResult, reference_key: str, use_as_prediction: str,
             reference_groups: list[str] | None, plot_group_probabilities: bool,
             output_dir: str | None) -> dict:
    obs = predictions.adata.obs
    if reference_key not in obs.columns:
        raise ValueError(f"reference_key {reference_key!r} is not an obs column of {_ANNOTATED_H5AD}; available columns: {list(obs.columns)}")
    if use_as_prediction not in predictions.predicted_labels.columns:
        raise ValueError(
            f"{_LABELS_CSV} has no {use_as_prediction!r} column (columns: {list(predictions.predicted_labels.columns)}). "
            "Majority voting was not run (or was skipped for <= 50 cells); use use_as_prediction='predicted_labels' or re-annotate with majority voting."
        )
    reference = obs[reference_key]
    prediction = predictions.predicted_labels[use_as_prediction]

    # Map user-supplied group names onto the reference values (compared as strings, so numeric cluster IDs work).
    present = {str(v): v for v in reference.dropna().unique()}
    groups = list(dict.fromkeys(reference_groups or []))
    absent = [g for g in groups if g not in present]
    if absent:
        shown = sorted(present)[:50]
        raise ValueError(f"Reference group(s) {absent} not found in obs[{reference_key!r}]; available groups: {shown}{' ...' if len(present) > 50 else ''}")

    base = Path(output_dir).expanduser() if output_dir else Path(tempfile.gettempdir()) / "celltypist_mcp_outputs"
    base.mkdir(parents=True, exist_ok=True)
    out = Path(tempfile.mkdtemp(prefix=f"compare_with_reference_{time.strftime('%Y%m%d_%H%M%S')}_", dir=base)).resolve()
    artifacts = []

    # Label-transfer dot plot: tutorial `celltypist.dotplot(predictions, use_as_reference=..., use_as_prediction=...)`.
    dp = celltypist.dotplot(predictions, use_as_reference=reference_key, use_as_prediction=use_as_prediction, return_fig=True)
    dp.make_figure()
    dotplot_png = out / "dotplot.png"
    dp.fig.savefig(dotplot_png, bbox_inches="tight", dpi=150)
    plt.close(dp.fig)
    fraction_csv = out / "dotplot_fraction.csv"
    mean_prob_csv = out / "dotplot_mean_probability.csv"
    dp.dot_size_df.to_csv(fraction_csv, index_label=use_as_prediction)
    dp.dot_color_df.to_csv(mean_prob_csv, index_label=use_as_prediction)
    artifacts += [
        {"description": f"CellTypist label-transfer dot plot ({use_as_prediction} rows vs {reference_key} columns)", "path": str(dotplot_png)},
        {"description": "Dot sizes: fraction of each reference group assigned to each predicted label", "path": str(fraction_csv)},
        {"description": "Dot colours: mean probability of the predicted label within each reference group", "path": str(mean_prob_csv)},
    ]

    # Contingency table: tutorial `pd.crosstab(adata.obs.cell_type, adata.obs.majority_voting)` (ml cell 47).
    crosstab = pd.crosstab(reference, prediction)
    crosstab_csv = out / "crosstab.csv"
    crosstab.to_csv(crosstab_csv)
    artifacts.append({"description": f"Cell counts: {reference_key} (rows) x {use_as_prediction} (columns)", "path": str(crosstab_csv)})

    selected = {}
    for i, name in enumerate(groups, start=1):
        value = present[name]
        row = crosstab.loc[value]
        row = row[row > 0].sort_values(ascending=False, kind="stable")
        selected[name] = {"n_cells": int(row.sum()), "prediction_counts": {str(k): int(v) for k, v in row.items()}}
        if plot_group_probabilities:
            # Tutorial ml cell 82: obs.loc[obs.cell_type == group, model.cell_types].plot(kind='box', rot=90, figsize=[8, 4], title=group);
            # probability_matrix columns are the model cell types in model order (equivalent to to_adata(insert_prob=True) columns).
            # The tutorial's 8-inch width (10-type model) is kept up to 32 cell types and widened beyond that so labels stay legible.
            probs = predictions.probability_matrix
            width = max(8.0, 0.25 * probs.shape[1])
            ax = probs.loc[(reference == value).to_numpy()].plot(kind="box", rot=90, figsize=[width, 4], title=name)
            box_png = out / f"probability_boxplot_{i:02d}_{_safe_name(name)}.png"
            ax.figure.savefig(box_png, bbox_inches="tight", dpi=150)
            plt.close(ax.figure)
            artifacts.append({"description": f"Probability of each model cell type across {reference_key} == {name!r} cells", "path": str(box_png)})

    result = {
        "message": (f"Compared {use_as_prediction} with obs[{reference_key!r}] for {predictions.cell_count} cells: "
                    f"{crosstab.shape[0]} reference groups x {crosstab.shape[1]} predicted labels."),
        "reference": REFERENCE,
        "artifacts": artifacts,
        "output_dir": str(out),
        "n_cells": int(predictions.cell_count),
        "reference_key": reference_key,
        "use_as_prediction": use_as_prediction,
        "n_reference_groups": int(crosstab.shape[0]),
        "n_predicted_labels": int(crosstab.shape[1]),
        "n_cells_missing_reference": int(reference.isna().sum()),
        "reference_group_sizes": {str(k): int(v) for k, v in crosstab.sum(axis=1).items()},
    }
    if groups:
        result["selected_groups"] = selected
    return result
