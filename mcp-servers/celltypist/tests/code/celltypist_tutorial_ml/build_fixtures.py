"""Build fixtures and direct-upstream expected values for the celltypist_compare_with_reference tests.

Every annotation directory is produced by direct upstream calls, which is the tool's frozen input
contract (the output of celltypist_annotate_cells):

    predictions = celltypist.annotate(adata, model=..., mode=..., p_thres=..., majority_voting=...)
    predictions.to_table(dir)                                   # predicted_labels / decision / probability CSVs
    predictions.to_adata(insert_labels=True, insert_conf=True).write_h5ad(dir / "annotated.h5ad")

Expected values come from the same in-memory AnnotationResult, with no CSV round-trip:

    celltypist.dotplot(predictions, use_as_reference=ref, use_as_prediction=pred, return_fig=True)
        -> .dot_size_df / .dot_color_df            (tutorial ml cells 49/51)
    pd.crosstab(adata.obs[ref], adata.obs[pred])   (tutorial ml cell 47)
    to_adata(insert_prob=True).obs[[ref] + list(model.cell_types)]   (tutorial ml cell 82 data)

The wrapper under test (src/tools/celltypist_tutorial_ml.py) is never imported here.

Run from the project root:
    CELLTYPIST_FOLDER=$PWD/tests/data/celltypist_tutorial_ml/celltypist_home \
        celltypist-env/bin/python tests/code/celltypist_tutorial_ml/build_fixtures.py
"""
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "tests/data/celltypist_tutorial_ml"
CT_HOME = DATA / "celltypist_home"
os.environ["CELLTYPIST_FOLDER"] = str(CT_HOME)  # must precede the celltypist import

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

import celltypist  # noqa: E402
from celltypist import models  # noqa: E402
from celltypist.classifier import AnnotationResult  # noqa: E402

FIXTURES = DATA / "fixtures"
EXPECTED = DATA / "expected"
DEMO_500 = DATA / "inputs/demo_500_cells.h5ad"
CUSTOM_MODEL = DATA / "inputs/model_from_immune2000_seed0_run1_njobs10.pkl"
# Main tutorial (celltypist_tutorial.ipynb cells 41/60) dataset, read in place (read-only executor evidence).
DEMO_2000 = ROOT / "notebooks/celltypist_tutorial_ml/data/demo_2000_cells.h5ad"

SPECS = {
    # Tutorial ml cell 34 (and main tutorial): best match + majority voting with the built-in model.
    "builtin_bestmatch": dict(model="Immune_All_Low.pkl", mode="best match", p_thres=0.5, majority_voting=True),
    # No majority voting: majority_voting column absent (error case) and predicted_labels-only use.
    "builtin_bestmatch_no_mv": dict(model="Immune_All_Low.pkl", mode="best match", p_thres=0.5, majority_voting=False),
    # Changed input: prob match with a stricter threshold than the tutorial's 0.5 (multi-label and Unassigned labels).
    "builtin_probmatch_p07": dict(model="Immune_All_Low.pkl", mode="prob match", p_thres=0.7, majority_voting=True),
    # Tutorial ml cell 77 with the executor's seeded custom model (deterministic stand-in for the unseeded as-run model).
    "custom_seed0_probmatch": dict(model=str(CUSTOM_MODEL), mode="prob match", p_thres=0.5, majority_voting=True),
    # Main tutorial cell 41: demo_2000 + Immune_All_Low, best match + majority voting (second dataset).
    "main_demo2000_bestmatch": dict(model="Immune_All_Low.pkl", mode="best match", p_thres=0.5, majority_voting=True),
}
INPUT_OF = {"main_demo2000_bestmatch": DEMO_2000}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_expected(name: str, pred: AnnotationResult, refs: list[str], model_cell_types: list[str]) -> dict:
    out = EXPECTED / name
    out.mkdir(parents=True)
    pred_cols = [c for c in ("predicted_labels", "majority_voting") if c in pred.predicted_labels.columns]
    obs = pred.adata.obs
    record = {"dotplot": {}, "crosstab": {}}
    for ref in refs:
        for p in pred_cols:
            dp = celltypist.dotplot(pred, use_as_reference=ref, use_as_prediction=p, return_fig=True)
            dp.make_figure()
            dp.fig.savefig(out / f"dotplot_{ref}_{p}.png", bbox_inches="tight", dpi=100)
            plt.close(dp.fig)
            pd.to_pickle(dp.dot_size_df, out / f"dotplot_{ref}_{p}_size.pkl")
            pd.to_pickle(dp.dot_color_df, out / f"dotplot_{ref}_{p}_color.pkl")
            record["dotplot"][f"{ref}|{p}"] = {"shape": list(dp.dot_size_df.shape),
                                               "rows": [str(x) for x in dp.dot_size_df.index],
                                               "columns": [str(x) for x in dp.dot_size_df.columns]}
            ct = pd.crosstab(obs[ref], obs[p])
            pd.to_pickle(ct, out / f"crosstab_{ref}_{p}.pkl")
            record["crosstab"][f"{ref}|{p}"] = {"shape": list(ct.shape), "total": int(ct.to_numpy().sum())}
    # Tutorial ml cell 82 data: probabilities inserted into obs by to_adata(insert_prob=True), model.cell_types order.
    prob_view = AnnotationResult(pred.predicted_labels, pred.decision_matrix, pred.probability_matrix, pred.adata.copy())
    prob_obs = prob_view.to_adata(insert_prob=True).obs[refs + model_cell_types]
    pd.to_pickle(prob_obs, out / "obs_probabilities.pkl")
    record["model_cell_types"] = model_cell_types
    return record


def build_annotation_dir(name: str, pred: AnnotationResult) -> dict:
    d = FIXTURES / name
    d.mkdir(parents=True)
    pred.to_table(str(d))
    pred.to_adata(insert_labels=True, insert_conf=True).write_h5ad(d / "annotated.h5ad")
    return {f.name: sha256(f) for f in sorted(d.iterdir())}


def main() -> None:
    assert Path(models.models_path).resolve() == (CT_HOME / "data/models").resolve(), models.models_path
    for p in (FIXTURES, EXPECTED):
        if p.exists():
            shutil.rmtree(p)
    manifest = {"celltypist": celltypist.__version__, "scanpy": sc.__version__, "pandas": pd.__version__,
                "models_path": models.models_path, "demo_500_sha256": sha256(DEMO_500),
                "custom_model_sha256": sha256(CUSTOM_MODEL), "demo_2000_sha256": sha256(DEMO_2000), "fixtures": {}}
    in_memory = {}
    for name, spec in SPECS.items():
        source = INPUT_OF.get(name, DEMO_500)
        adata = sc.read(str(source))
        pred = celltypist.annotate(adata, **spec)
        model_cell_types = [str(x) for x in models.Model.load(spec["model"]).cell_types]
        files = build_annotation_dir(name, pred)
        record = write_expected(name, pred, ["cell_type"], model_cell_types)
        manifest["fixtures"][name] = {"call": f"celltypist.annotate({source.name}, **{spec})",
                                      "files_sha256": files, **record}
        in_memory[name] = pred

    # Changed reference column: the best-match result with the author column renamed and recoded
    # (reverse-alphabetical codes, so any reliance on the tutorial names or their order is visible),
    # plus an integer-coded cluster column taken from the over-clustering. 'cell_type' is removed.
    base = in_memory["builtin_bestmatch"]
    adata = base.adata.copy()
    types = sorted(adata.obs["cell_type"].astype(str).unique())
    recode = {t: f"pop_{i:02d}" for i, t in enumerate(reversed(types))}
    adata.obs["expert_annotation"] = pd.Categorical(adata.obs["cell_type"].astype(str).map(recode))
    adata.obs["leiden_int"] = base.predicted_labels["over_clustering"].astype(int).to_numpy()
    del adata.obs["cell_type"]
    relabelled = AnnotationResult(base.predicted_labels, base.decision_matrix, base.probability_matrix, adata)
    files = build_annotation_dir("relabelled_reference", relabelled)
    model_cell_types = [str(x) for x in models.Model.load("Immune_All_Low.pkl").cell_types]
    record = write_expected("relabelled_reference", relabelled, ["expert_annotation", "leiden_int"], model_cell_types)
    manifest["fixtures"]["relabelled_reference"] = {
        "call": "builtin_bestmatch AnnotationResult; obs 'cell_type' recoded to 'expert_annotation' and removed; "
                "'leiden_int' = int(predicted_labels.over_clustering)",
        "recode": recode, "files_sha256": files, **record}

    (DATA / "fixtures_manifest.json").write_text(json.dumps(manifest, indent=1))
    print(json.dumps({k: v["files_sha256"] for k, v in manifest["fixtures"].items()}, indent=1), file=sys.stderr)


if __name__ == "__main__":
    main()
