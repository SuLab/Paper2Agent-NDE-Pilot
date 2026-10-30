"""Verification of celltypist_train_model (tutorial cells 49-54: celltypist.train + Model.write/Model.load).

Expectations: the executor's seeded reference notebooks/celltypist_tutorial/reference/train_seeded_rs0/run1/model_arrays.npz
(celltypist.train(demo_2000, labels='cell_type', n_jobs=10, feature_selection=True, random_state=0), bit-identical across
three runs in two processes) and new direct celltypist.train calls on the same fixture files for changed inputs.
Model arrays are compared exactly (same machine, same n_jobs; see executor seed_provenance).
"""
import inspect
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import scanpy as sc

import celltypist
from celltypist import models as ct_models
from fastmcp import Client

from _celltypist_tutorial_support import (INPUTS, REF, assert_annotation_matches, assert_artifacts, call_tool,
                                          call_tool_error, celltypist_tutorial_mcp, read_labels, run_tools_in_subprocess,
                                          seeded_cache_copy, fresh_tmp, sha256)

TOOL = "celltypist_train_model"
D2000 = INPUTS / "demo_2000_cells.h5ad"
D400 = INPUTS / "demo_400_cells.h5ad"
USER_META = INPUTS / "demo_400_user_meta.h5ad"
COUNTS_LIKE = INPUTS / "demo_400_counts_like.h5ad"
D2000_SHA = "167d038d57bc7aee2854f12c3337cc1f58846f4ae7382745639a0957831bae9b"

_SEEDED: dict = {}


def assert_same_model(got: ct_models.Model, exp: ct_models.Model):
    np.testing.assert_array_equal(np.asarray(got.cell_types).astype(str), np.asarray(exp.cell_types).astype(str))
    np.testing.assert_array_equal(np.asarray(got.features).astype(str), np.asarray(exp.features).astype(str))
    np.testing.assert_array_equal(got.classifier.coef_, exp.classifier.coef_)
    np.testing.assert_array_equal(got.classifier.intercept_, exp.classifier.intercept_)
    np.testing.assert_array_equal(got.scaler.mean_, exp.scaler.mean_)
    np.testing.assert_array_equal(got.scaler.scale_, exp.scaler.scale_)


async def seeded_tutorial_training() -> dict:
    """Tutorial cell 49 with random_state=0 through the tool, once per session (about 12 s)."""
    if not _SEEDED:
        base = Path(__file__).resolve().parents[2] / "results" / "celltypist_tutorial" / "shared_train_seeded"
        base.mkdir(parents=True, exist_ok=True)
        before = sha256(D2000)
        _SEEDED["result"] = await call_tool(TOOL, {"adata_path": str(D2000), "label_key": "cell_type", "feature_selection": True,
                                                   "n_jobs": 10, "seed": 0, "model_name": "model_from_immune2000",
                                                   "output_dir": str(base)})
        _SEEDED["base"], _SEEDED["sha"] = base, (before, sha256(D2000))
    return _SEEDED


async def test_schema_defaults_follow_upstream():
    async with Client(celltypist_tutorial_mcp) as client:
        tool = {t.name: t for t in await client.list_tools()}[TOOL]
    schema = tool.input_schema
    props = schema["properties"]
    assert set(schema["required"]) == {"adata_path", "label_key"}
    assert set(props) == {"adata_path", "label_key", "feature_selection", "top_genes", "n_jobs", "seed", "model_name", "output_dir"}
    assert all(p.get("description") for p in props.values())
    up = inspect.signature(celltypist.train).parameters
    assert props["feature_selection"]["default"] == up["feature_selection"].default is False
    assert props["top_genes"]["default"] == up["top_genes"].default == 300
    assert props["n_jobs"]["default"] is None and up["n_jobs"].default is None
    assert props["seed"]["default"] is None
    assert len(tool.description.strip().splitlines()) == 2


async def test_seeded_tutorial_training_matches_executor_reference():
    shared = await seeded_tutorial_training()
    result = shared["result"]
    (pkl,) = assert_artifacts(result, shared["base"], 1)
    assert pkl.name == "model_from_immune2000.pkl" and Path(result["model_path"]) == pkl
    model = ct_models.Model.load(str(pkl))
    ref = np.load(REF / "train_seeded_rs0" / "run1" / "model_arrays.npz")
    np.testing.assert_array_equal(np.asarray(model.features).astype(str), ref["features"])   # exact gene list and order
    np.testing.assert_array_equal(np.asarray(model.cell_types).astype(str), ref["cell_types"])
    np.testing.assert_array_equal(model.classifier.coef_, ref["coef"])
    np.testing.assert_array_equal(model.classifier.intercept_, ref["intercept"])
    np.testing.assert_array_equal(model.scaler.mean_, ref["scaler_mean"])
    np.testing.assert_array_equal(model.scaler.scale_, ref["scaler_scale"])
    assert len(model.features) == 2652 and model.classifier.coef_.shape == (10, 2652)

    # The executor's saved seeded .pkl (same call) loads to the same arrays.
    assert_same_model(model, ct_models.Model.load(str(Path(__file__).resolve().parents[2] / "data" / "celltypist_tutorial" /
                                                      "models" / "model_seeded_rs0_run1.pkl")))
    assert result["n_features"] == 2652 and result["n_cell_types"] == 10 and result["n_cells"] == 2000
    assert result["cell_types"] == [str(x) for x in ref["cell_types"]]
    assert result["seed"] == 0 and result["n_jobs"] == 10 and result["feature_selection"] is True and result["top_genes"] == 300
    assert result["description"]["number_celltypes"] == "10" and result["description"]["details"] == ""
    assert shared["sha"][0] == shared["sha"][1] == D2000_SHA


async def test_trained_model_composes_with_annotate(out_base):
    """train -> annotate (tutorial cell 56 with the seeded model): equals the executor's direct-upstream reference."""
    shared = await seeded_tutorial_training()
    res = await call_tool("celltypist_annotate_cells", {"adata_path": str(D400), "model": shared["result"]["model_path"],
                                                         "majority_voting": True, "output_dir": str(out_base)})
    folder = REF / "annotate_demo400_seeded_rs0_mv"
    import pandas as pd
    z = np.load(folder / "matrices.npz")
    cells, classes = z["cells"].astype(str), z["classes"].astype(str)
    assert_annotation_matches(Path(res["output_dir"]), read_labels(folder / "predicted_labels.csv"),
                              pd.DataFrame(z["decision"], index=cells, columns=classes),
                              pd.DataFrame(z["probability"], index=cells, columns=classes))


async def test_renamed_label_column_drives_training(out_base):
    """User labels in a non-tutorial column (fixture has no 'cell_type'); expectation from a direct upstream train."""
    before = sha256(USER_META)
    result = await call_tool(TOOL, {"adata_path": str(USER_META), "label_key": "coarse_group", "feature_selection": True,
                                    "top_genes": 100, "n_jobs": 10, "seed": 1, "output_dir": str(out_base)})
    model = ct_models.Model.load(result["model_path"])
    direct = celltypist.train(sc.read(str(USER_META)), labels="coarse_group", feature_selection=True, top_genes=100,
                              n_jobs=10, random_state=1)
    assert_same_model(model, direct)
    assert result["cell_types"] == ["B lineage", "Myeloid", "Stromal", "T cells"]
    assert result["label_key"] == "coarse_group" and result["top_genes"] == 100
    assert sha256(USER_META) == before


async def test_plain_logistic_regression_without_feature_selection(out_base):
    result = await call_tool(TOOL, {"adata_path": str(USER_META), "label_key": "author_label", "output_dir": str(out_base)})
    model = ct_models.Model.load(result["model_path"])
    direct = celltypist.train(sc.read(str(USER_META)), labels="author_label")
    assert_same_model(model, direct)
    # all expressed genes are kept when feature_selection is off
    assert result["n_features"] == len(direct.features) and result["feature_selection"] is False
    assert Path(result["model_path"]).name == "celltypist_model.pkl"


async def test_repeated_seeded_calls_fresh_dirs_identical_models(out_base):
    args = {"adata_path": str(USER_META), "label_key": "coarse_group", "feature_selection": True, "top_genes": 50,
            "seed": 3, "model_name": "coarse.pkl", "output_dir": str(out_base)}
    r1, r2 = await call_tool(TOOL, args), await call_tool(TOOL, args)
    assert r1["output_dir"] != r2["output_dir"] and r1["model_path"] != r2["model_path"]
    assert Path(r1["model_path"]).name == "coarse.pkl"
    assert_same_model(ct_models.Model.load(r1["model_path"]), ct_models.Model.load(r2["model_path"]))


# ------------------------------------------------------------------------------------------------------ errors
async def test_missing_label_column(out_base):
    """A hard-coded tutorial 'cell_type' would fail on this fixture: the error names the user's columns."""
    msg = await call_tool_error(TOOL, {"adata_path": str(USER_META), "label_key": "cell_type", "output_dir": str(out_base)})
    assert "label_key 'cell_type' is not an obs column" in msg and "author_label" in msg
    assert list(out_base.iterdir()) == []


async def test_missing_input_file(out_base):
    msg = await call_tool_error(TOOL, {"adata_path": str(INPUTS / "absent.h5ad"), "label_key": "cell_type",
                                       "output_dir": str(out_base)})
    assert "No such file" in msg


async def test_non_normalised_input_surfaces_upstream_check(out_base):
    with pytest.raises(ValueError) as direct:
        celltypist.train(sc.read(str(COUNTS_LIKE)), labels="cell_type")
    before = sha256(COUNTS_LIKE)
    msg = await call_tool_error(TOOL, {"adata_path": str(COUNTS_LIKE), "label_key": "cell_type", "output_dir": str(out_base)})
    assert "expect log1p normalized expression to 10000 counts per cell" in str(direct.value)
    assert str(direct.value) in msg
    assert list(out_base.iterdir()) == [] and sha256(COUNTS_LIKE) == before


async def test_top_genes_larger_than_gene_count_surfaces_upstream_error(out_base):
    msg = await call_tool_error(TOOL, {"adata_path": str(USER_META), "label_key": "coarse_group", "feature_selection": True,
                                       "top_genes": 50000, "output_dir": str(out_base)})
    assert "is fewer than the `top_genes` (50000)" in msg


async def test_invalid_top_genes_and_model_name(out_base):
    assert "top_genes must be a positive integer" in await call_tool_error(
        TOOL, {"adata_path": str(USER_META), "label_key": "coarse_group", "top_genes": 0, "output_dir": str(out_base)})
    assert "model_name must be a plain file name" in await call_tool_error(
        TOOL, {"adata_path": str(USER_META), "label_key": "coarse_group", "model_name": "../escape", "output_dir": str(out_base)})
    assert list(out_base.iterdir()) == []


# ---------------------------------------------------------------------------------- protocol stdout (fd level)
def test_no_stdout_from_parallel_training_annotation_and_plotting():
    """fd-1 capture of a fresh process running train (n_jobs=4, joblib workers), annotate with Leiden majority voting,
    and inspect with violin plots: an MCP stdio server must not write anything but protocol frames to stdout."""
    home = seeded_cache_copy("stdout_cache")
    out = fresh_tmp("stdout_out")
    calls = [
        {"tool": TOOL, "arguments": {"adata_path": str(USER_META), "label_key": "coarse_group", "feature_selection": True,
                                     "top_genes": 50, "n_jobs": 4, "seed": 0, "output_dir": str(out)}},
        {"tool": "celltypist_annotate_cells", "arguments": {"adata_path": str(D400), "model": "Immune_All_Low.pkl",
                                                            "majority_voting": True, "output_dir": str(out)}},
        {"tool": "celltypist_inspect_model", "arguments": {"model": "Immune_All_Low.pkl", "cell_types": ["Mast cells"], "top_n": 3,
                                                           "adata_path": str(D400), "groupby_key": "cell_type", "output_dir": str(out)}},
        {"tool": "celltypist_list_models", "arguments": {"output_dir": str(out)}},
        {"tool": "celltypist_inspect_model", "arguments": {"model": "No_Such_Model.pkl", "output_dir": str(out)}},
    ]
    records, stdout, stderr, rc = run_tools_in_subprocess(home, calls, "stdout_run")
    assert rc == 0, stderr[-3000:]
    assert stdout == b"", stdout[:1000]
    assert all("error" not in r for r in records[1:5]), records
    assert "No models match" in records[5]["error"]
    assert all(r["cache_after"] == ["Immune_All_Low.pkl", "models.json"] for r in records[1:])
