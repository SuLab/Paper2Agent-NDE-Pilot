"""Verification of celltypist_annotate_cells (tutorial cells 28-35 and 56-57; prob match from celltypist_tutorial_ml 55/77).

Expectations: the executor's saved direct-upstream outputs under notebooks/celltypist_tutorial/reference/
(annotate_demo2000_ImmuneAllLow_mv / _nomv, annotate_demo400_seeded_rs0_mv, notebook_run/) and new direct
celltypist.annotate calls on the same fixture files for changed inputs.
Labels / cluster IDs: exact. Decision scores: atol 1e-9, probabilities / conf_score: atol 2.5e-10 (executor replay policy).
"""
import asyncio
import inspect
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import celltypist
from celltypist import models as ct_models
from fastmcp import Client
from fastmcp.client.transports import StdioTransport

from _celltypist_tutorial_support import (CACHE_MODELS, CONTRACT_FILES, INPUTS, MODELS, PROBABILITY_ATOL, PROJECT_PYTHON,
                                          PROJECT_ROOT, REF, assert_annotation_matches, assert_artifacts, call_tool,
                                          call_tool_error, celltypist_tutorial_mcp, fresh_tmp, network_available,
                                          read_labels, read_matrix, run_tools_in_subprocess, seeded_cache_copy, sha256)

TOOL = "celltypist_annotate_cells"
D2000 = INPUTS / "demo_2000_cells.h5ad"
D400 = INPUTS / "demo_400_cells.h5ad"
USER_META = INPUTS / "demo_400_user_meta.h5ad"
FIRST40 = INPUTS / "demo_400_first40.h5ad"
COUNTS_LIKE = INPUTS / "demo_400_counts_like.h5ad"
SEEDED = MODELS / "model_seeded_rs0_run1.pkl"
IAL_CACHED = CACHE_MODELS / "Immune_All_Low.pkl"
D2000_SHA = "167d038d57bc7aee2854f12c3337cc1f58846f4ae7382745639a0957831bae9b"
D400_SHA = "7f1efb7242cf043b6e91b6fdeb99ad025cba152887c5a4567d54f6a748a6a119"

_MV2000: dict = {}


async def annotate_demo2000_mv() -> dict:
    """Tutorial cell 28 through the tool, computed once per session (about 35 s) and shared by several tests."""
    if not _MV2000:
        base = PROJECT_ROOT / "tests" / "results" / "celltypist_tutorial" / "shared_demo2000_mv"
        base.mkdir(parents=True, exist_ok=True)
        before = sha256(D2000)
        _MV2000["result"] = await call_tool(TOOL, {"adata_path": str(D2000), "model": "Immune_All_Low.pkl",
                                                   "majority_voting": True, "output_dir": str(base)})
        _MV2000["base"] = base
        _MV2000["sha_before"], _MV2000["sha_after"] = before, sha256(D2000)
    return _MV2000


def reference_tables(folder: Path):
    labels = read_labels(folder / "predicted_labels.csv")
    z = np.load(folder / "matrices.npz")
    cells, classes = z["cells"].astype(str), z["classes"].astype(str)
    dec = pd.DataFrame(z["decision"], index=cells, columns=classes)
    prob = pd.DataFrame(z["probability"], index=cells, columns=classes)
    return labels, dec, prob


def direct_tables(pred):
    return pred.predicted_labels.astype(str), pred.decision_matrix, pred.probability_matrix


def frequencies_from_csv(path) -> dict:
    df = pd.read_csv(path, dtype={"celltype": str})
    return {r.celltype: int(r.counts) for r in df.itertuples()}


# --------------------------------------------------------------------------------------------- schema / format
async def test_schema_defaults_follow_upstream():
    async with Client(celltypist_tutorial_mcp) as client:
        tool = {t.name: t for t in await client.list_tools()}[TOOL]
    schema = tool.input_schema if hasattr(tool, "input_schema") else tool.inputSchema
    props = schema["properties"]
    assert set(schema["required"]) == {"adata_path"}
    assert set(props) == {"adata_path", "model", "mode", "p_thres", "majority_voting", "over_clustering_key", "min_prop", "output_dir"}
    assert all(p.get("description") for p in props.values())
    assert props["mode"]["enum"] == ["best match", "prob match"]
    up = inspect.signature(celltypist.annotate).parameters
    assert props["mode"]["default"] == up["mode"].default
    assert props["p_thres"]["default"] == up["p_thres"].default
    assert props["majority_voting"]["default"] == up["majority_voting"].default
    assert props["min_prop"]["default"] == up["min_prop"].default
    assert props["over_clustering_key"]["default"] is None and up["over_clustering"].default is None
    # upstream model=None resolves to the index's default model
    assert props["model"]["default"] == ct_models.get_default_model() == "Immune_All_Low.pkl"
    assert len(tool.description.strip().splitlines()) == 2


# ------------------------------------------------------------------------------- tutorial reference (cell 28)
async def test_demo2000_majority_voting_matches_executor_reference():
    shared = await annotate_demo2000_mv()
    result = shared["result"]
    paths = assert_artifacts(result, shared["base"], 4)
    out = Path(result["output_dir"])
    assert {p.name for p in paths} == CONTRACT_FILES
    assert {p.name for p in out.iterdir()} == CONTRACT_FILES

    labels, dec, prob = reference_tables(REF / "annotate_demo2000_ImmuneAllLow_mv")
    assert list(labels.columns) == ["predicted_labels", "over_clustering", "majority_voting"]
    assert_annotation_matches(out, labels, dec, prob)

    # The notebook's own cell-28 labels (captured from the kernel's Out cache) are reproduced as well.
    nb_labels = read_labels(REF / "notebook_run" / "demo2000_ImmuneAllLow_mv" / "predicted_labels.csv")
    pd.testing.assert_frame_equal(read_labels(out / "predicted_labels.csv"), nb_labels)

    assert result["n_cells"] == 2000 and result["n_model_cell_types"] == 98
    assert result["majority_voting_requested"] is True and result["majority_voting_performed"] is True
    assert result["n_over_clusters"] == 46 == labels["over_clustering"].nunique()
    assert result["model"] == "Immune_All_Low.pkl" and Path(result["model_file"]) == IAL_CACHED
    assert result["over_clustering_source"] == "celltypist Leiden over-clustering"
    for by in ("predicted_labels", "majority_voting"):
        got = {r["celltype"]: r["counts"] for r in result["label_frequencies"][by]}
        assert got == frequencies_from_csv(REF / "annotate_demo2000_ImmuneAllLow_mv" / f"summary_frequency_{by}.csv")


async def test_demo2000_annotated_h5ad_matches_to_adata_reference():
    shared = await annotate_demo2000_mv()
    out = Path(shared["result"]["output_dir"])
    got = ad.read_h5ad(out / "annotated.h5ad")
    ref_obs = pd.read_csv(REF / "annotate_demo2000_ImmuneAllLow_mv" / "obs_after_to_adata.csv", index_col=0,
                          dtype={"over_clustering": str}, keep_default_na=False)
    assert list(got.obs.columns) == list(ref_obs.columns) == ["cell_type", "predicted_labels", "over_clustering",
                                                              "majority_voting", "conf_score"]
    assert got.obs_names.equals(ref_obs.index)
    for col in ["cell_type", "predicted_labels", "over_clustering", "majority_voting"]:
        assert (got.obs[col].astype(str).to_numpy() == ref_obs[col].astype(str).to_numpy()).all(), col
    np.testing.assert_allclose(got.obs["conf_score"].to_numpy(), ref_obs["conf_score"].to_numpy(), rtol=0, atol=PROBABILITY_ATOL)
    # conf_score is max probability (to_adata insert_conf_by='predicted_labels'); probabilities are not inserted into obs.
    prob = read_matrix(out / "probability_matrix.csv")
    np.testing.assert_array_equal(got.obs["conf_score"].to_numpy(), prob.max(axis=1).to_numpy())
    # Expression matrix is the user's, unchanged; upstream's over-clustering side effects are kept.
    src = ad.read_h5ad(D2000)
    assert got.var_names.equals(src.var_names)
    assert (got.X != src.X).nnz == 0
    assert "connectivities" in got.obsp and "X_pca" in got.obsm


async def test_demo2000_input_file_unchanged():
    shared = await annotate_demo2000_mv()
    assert shared["sha_before"] == shared["sha_after"] == D2000_SHA


async def test_composition_with_compare_with_reference_tool():
    """Frozen cross-module contract: celltypist_tutorial_ml reads this tool's output_dir (read-only use of that module)."""
    from tools.celltypist_tutorial_ml import celltypist_tutorial_ml_mcp
    shared = await annotate_demo2000_mv()
    base = fresh_tmp("composition")
    async with Client(celltypist_tutorial_ml_mcp) as client:
        res = (await client.call_tool("celltypist_compare_with_reference", {
            "annotation_dir": shared["result"]["output_dir"], "reference_key": "cell_type",
            "use_as_prediction": "majority_voting", "output_dir": str(base)})).data
    out = Path(res["output_dir"])
    frac = pd.read_csv(out / "dotplot_fraction.csv", index_col=0)
    ref_frac = pd.read_csv(REF / "annotate_demo2000_ImmuneAllLow_mv" / "dotplot_majority_voting_fraction.csv", index_col=0)
    assert list(frac.index) == list(ref_frac.index) and list(frac.columns) == list(ref_frac.columns)
    np.testing.assert_allclose(frac.to_numpy(), ref_frac.to_numpy(), rtol=0, atol=1e-12)
    xt = pd.read_csv(out / "crosstab.csv", index_col=0)
    ref_xt = pd.read_csv(REF / "annotate_demo2000_ImmuneAllLow_mv" / "crosstab_cell_type_vs_majority_voting.csv", index_col=0)
    pd.testing.assert_frame_equal(xt, ref_xt, check_names=False)


# ----------------------------------------------------------------------------- reference without majority voting
async def test_demo2000_default_model_without_majority_voting(out_base):
    # model omitted -> tool default Immune_All_Low.pkl (= upstream default model); majority_voting default False.
    result = await call_tool(TOOL, {"adata_path": str(D2000), "output_dir": str(out_base)})
    assert_artifacts(result, out_base, 4)
    out = Path(result["output_dir"])
    labels, dec, prob = reference_tables(REF / "annotate_demo2000_ImmuneAllLow_nomv")
    assert list(labels.columns) == ["predicted_labels"]
    assert_annotation_matches(out, labels, dec, prob)
    assert result["majority_voting_performed"] is False and "n_over_clusters" not in result
    assert set(result["label_frequencies"]) == {"predicted_labels"}
    assert sha256(D2000) == D2000_SHA


# --------------------------------------------------------------------------- changed inputs: custom .pkl model
async def test_custom_seeded_model_on_demo400_matches_executor_reference(out_base):
    result = await call_tool(TOOL, {"adata_path": str(D400), "model": str(SEEDED), "majority_voting": True,
                                    "output_dir": str(out_base)})
    assert_artifacts(result, out_base, 4)
    labels, dec, prob = reference_tables(REF / "annotate_demo400_seeded_rs0_mv")
    assert_annotation_matches(Path(result["output_dir"]), labels, dec, prob)
    assert result["n_model_cell_types"] == 10 and Path(result["model_file"]) == SEEDED
    assert sha256(D400) == D400_SHA


async def test_custom_model_csv_format_is_upstream_to_table(out_base):
    """The three CSVs are AnnotationResult.to_table(folder) output: compare with a direct upstream to_table."""
    result = await call_tool(TOOL, {"adata_path": str(D400), "model": str(SEEDED), "output_dir": str(out_base)})
    out = Path(result["output_dir"])
    direct_dir = fresh_tmp("direct_to_table")
    pred = celltypist.annotate(sc.read(str(D400)), model=str(SEEDED))
    pred.to_table(str(direct_dir))
    assert (out / "predicted_labels.csv").read_bytes() == (direct_dir / "predicted_labels.csv").read_bytes()
    for name in ("decision_matrix.csv", "probability_matrix.csv"):
        got_lines = (out / name).read_text().splitlines()
        exp_lines = (direct_dir / name).read_text().splitlines()
        assert got_lines[0] == exp_lines[0] and len(got_lines) == len(exp_lines)
        assert [l.split(",", 1)[0] for l in got_lines] == [l.split(",", 1)[0] for l in exp_lines]
    assert_annotation_matches(out, *direct_tables(pred))


# ------------------------------------------------------------------- changed inputs: prob match, non-default p_thres
async def test_prob_match_builtin_model_nondefault_threshold(out_base):
    result = await call_tool(TOOL, {"adata_path": str(D2000), "model": "Immune_All_Low.pkl", "mode": "prob match",
                                    "p_thres": 0.7, "output_dir": str(out_base)})
    pred = celltypist.annotate(sc.read(str(D2000)), model=str(IAL_CACHED), mode="prob match", p_thres=0.7)
    assert_annotation_matches(Path(result["output_dir"]), *direct_tables(pred))
    got = read_labels(Path(result["output_dir"]) / "predicted_labels.csv")["predicted_labels"]
    assert (got == "Unassigned").sum() == 149 and got.str.contains("|", regex=False).sum() == 32
    # p_thres is honoured: the default threshold gives different multi-label calls.
    at_default = celltypist.annotate(sc.read(str(D2000)), model=str(IAL_CACHED), mode="prob match").predicted_labels
    assert (at_default["predicted_labels"].astype(str).to_numpy() != got.to_numpy()).sum() > 0
    assert result["mode"] == "prob match" and result["p_thres"] == 0.7


async def test_prob_match_custom_model_with_majority_voting(out_base):
    result = await call_tool(TOOL, {"adata_path": str(D400), "model": str(SEEDED), "mode": "prob match", "p_thres": 0.3,
                                    "majority_voting": True, "output_dir": str(out_base)})
    pred = celltypist.annotate(sc.read(str(D400)), model=str(SEEDED), mode="prob match", p_thres=0.3, majority_voting=True)
    assert list(pred.predicted_labels.columns) == ["predicted_labels", "over_clustering", "majority_voting"]
    assert_annotation_matches(Path(result["output_dir"]), *direct_tables(pred))
    assert result["majority_voting_performed"] is True


# --------------------------------------------------------- changed inputs: user over-clustering key + min_prop
async def test_user_over_clustering_key_and_min_prop(out_base):
    before = sha256(USER_META)
    result = await call_tool(TOOL, {"adata_path": str(USER_META), "model": str(SEEDED), "majority_voting": True,
                                    "over_clustering_key": "my_clusters", "min_prop": 0.6, "output_dir": str(out_base)})
    pred = celltypist.annotate(sc.read(str(USER_META)), model=str(SEEDED), majority_voting=True,
                               over_clustering="my_clusters", min_prop=0.6)
    out = Path(result["output_dir"])
    assert_annotation_matches(out, *direct_tables(pred))
    got = read_labels(out / "predicted_labels.csv")
    user = ad.read_h5ad(USER_META).obs["my_clusters"].astype(str)
    assert (got["over_clustering"].to_numpy() == user.to_numpy()).all()   # the user's clusters drove the vote
    assert (got["majority_voting"] == "Heterogeneous").sum() == 210       # min_prop honoured (3 of 7 clusters)
    assert result["over_clustering_source"] == "my_clusters" and result["n_over_clusters"] == 7 and result["min_prop"] == 0.6
    annotated = ad.read_h5ad(out / "annotated.h5ad")
    assert "connectivities" not in annotated.obsp          # no Leiden graph built when the user supplies clusters
    assert list(annotated.obs.columns[:3]) == ["author_label", "coarse_group", "my_clusters"]
    assert sha256(USER_META) == before


async def test_majority_voting_skipped_for_50_or_fewer_cells(out_base):
    result = await call_tool(TOOL, {"adata_path": str(FIRST40), "model": str(SEEDED), "majority_voting": True,
                                    "output_dir": str(out_base)})
    pred = celltypist.annotate(sc.read(str(FIRST40)), model=str(SEEDED), majority_voting=True)
    assert list(pred.predicted_labels.columns) == ["predicted_labels"]
    assert_annotation_matches(Path(result["output_dir"]), *direct_tables(pred))
    assert result["majority_voting_requested"] is True and result["majority_voting_performed"] is False
    assert "skipped" in result["message"]
    assert {p.name for p in Path(result["output_dir"]).iterdir()} == CONTRACT_FILES


async def test_repeated_calls_write_fresh_directories(out_base):
    args = {"adata_path": str(FIRST40), "model": str(SEEDED), "output_dir": str(out_base)}
    r1, r2 = await call_tool(TOOL, args), await call_tool(TOOL, args)
    assert r1["output_dir"] != r2["output_dir"]
    for a1, a2 in zip(r1["artifacts"], r2["artifacts"]):
        assert a1["path"] != a2["path"]
    for name in ("predicted_labels.csv", "decision_matrix.csv", "probability_matrix.csv"):
        assert (Path(r1["output_dir"]) / name).read_bytes() == (Path(r2["output_dir"]) / name).read_bytes()


# ------------------------------------------------------------------------------------------------------ errors
async def test_missing_input_file(out_base):
    msg = await call_tool_error(TOOL, {"adata_path": str(INPUTS / "absent.h5ad"), "output_dir": str(out_base)})
    assert "No such file" in msg
    assert list(out_base.iterdir()) == []


async def test_non_h5ad_input_rejected(out_base):
    msg = await call_tool_error(TOOL, {"adata_path": str(INPUTS / "not_anndata.csv"), "output_dir": str(out_base)})
    assert "Expected an AnnData .h5ad file" in msg


async def test_missing_custom_model_file(out_base):
    msg = await call_tool_error(TOOL, {"adata_path": str(FIRST40), "model": str(MODELS / "absent.pkl"), "output_dir": str(out_base)})
    assert "No such model file" in msg


async def test_missing_over_clustering_key(out_base):
    before = sha256(USER_META)
    msg = await call_tool_error(TOOL, {"adata_path": str(USER_META), "model": str(SEEDED), "majority_voting": True,
                                       "over_clustering_key": "cell_type", "output_dir": str(out_base)})
    assert "over_clustering_key 'cell_type' is not an obs column" in msg
    assert list(out_base.iterdir()) == [] and sha256(USER_META) == before


async def test_non_normalised_input_surfaces_upstream_check(out_base):
    with pytest.raises(ValueError) as direct:
        celltypist.annotate(sc.read(str(COUNTS_LIKE)), model=str(SEEDED))
    before = sha256(COUNTS_LIKE)
    msg = await call_tool_error(TOOL, {"adata_path": str(COUNTS_LIKE), "model": str(SEEDED), "output_dir": str(out_base)})
    assert "Invalid expression matrix in `.X`" in str(direct.value)
    assert str(direct.value) in msg
    assert list(out_base.iterdir()) == [] and sha256(COUNTS_LIKE) == before


async def test_invalid_mode_rejected_by_schema(out_base):
    async with Client(celltypist_tutorial_mcp) as client:
        with pytest.raises(Exception) as info:
            await client.call_tool(TOOL, {"adata_path": str(FIRST40), "mode": "best", "output_dir": str(out_base)})
    assert "best match" in str(info.value) and "prob match" in str(info.value)


class _RequestRecorder:
    """Stands in for celltypist.models._requests_get: records URLs and fails like an unreachable server."""

    def __init__(self):
        self.urls = []

    def __call__(self, url, timeout=30):
        self.urls.append(url)
        raise Exception(f"blocked in test: {url}")


async def test_unknown_model_name_errors_without_downloads(out_base, monkeypatch):
    rec = _RequestRecorder()
    monkeypatch.setattr(ct_models, "_requests_get", rec)
    msg = await call_tool_error(TOOL, {"adata_path": str(FIRST40), "model": "No_Such_Model.pkl", "output_dir": str(out_base)})
    assert "No models match the celltypist model repertoire" in msg
    assert rec.urls == []


async def test_unknown_model_on_cache_without_pkl_does_not_download_all_models(out_base, monkeypatch):
    """Offline proof of the _resolve_model guard: an index-only cache (no .pkl) makes bare upstream Model.load attempt
    to download every model in the index (get_all_models -> download_if_required), while the tool attempts none."""
    home = fresh_tmp("index_only_cache") / "celltypist_home"
    mdir = home / "data" / "models"
    mdir.mkdir(parents=True)
    (mdir / "models.json").write_bytes((CACHE_MODELS / "models.json").read_bytes())
    monkeypatch.setattr(ct_models, "models_path", str(mdir))
    n_index = len(json.loads((mdir / "models.json").read_text())["models"])

    rec = _RequestRecorder()
    monkeypatch.setattr(ct_models, "_requests_get", rec)
    msg = await call_tool_error(TOOL, {"adata_path": str(FIRST40), "model": "No_Such_Model.pkl", "output_dir": str(out_base)})
    assert "No models match the celltypist model repertoire" in msg
    assert rec.urls == [] and sorted(os.listdir(mdir)) == ["models.json"]

    # Contrast: the bare upstream name lookup tries to fetch all models in the index.
    bare = _RequestRecorder()
    monkeypatch.setattr(ct_models, "_requests_get", bare)
    with pytest.raises(FileNotFoundError):
        ct_models.Model.load("No_Such_Model.pkl")
    assert len(bare.urls) == n_index == 61


# ---------------------------------------------------------------------------------- network (separate caches)
@pytest.mark.network
def test_by_name_load_on_empty_cache_downloads_only_requested_model():
    if not network_available():
        pytest.skip("celltypist.cog.sanger.ac.uk unreachable")
    home = fresh_tmp("net_empty_cache") / "celltypist_home"
    out = fresh_tmp("net_empty_cache_out")
    calls = [
        {"tool": TOOL, "arguments": {"adata_path": str(FIRST40), "model": "Immune_All_Low.pkl", "output_dir": str(out)}},
        {"tool": TOOL, "arguments": {"adata_path": str(FIRST40), "model": "No_Such_Model.pkl", "output_dir": str(out)}},
    ]
    records, stdout, stderr, rc = run_tools_in_subprocess(home, calls, "net_empty_cache_run")
    assert rc == 0, stderr[-3000:]
    assert stdout == b""
    assert Path(records[0]["models_path"]) == home / "data" / "models"
    ok, bad = records[1], records[2]
    assert "error" not in ok, ok
    assert ok["cache_after"] == ["Immune_All_Low.pkl", "models.json"]
    assert "No models match the celltypist model repertoire" in bad["error"]
    assert bad["cache_after"] == ["Immune_All_Low.pkl", "models.json"]
    # Results on the freshly downloaded model: compare with a direct upstream call on that same downloaded file.
    downloaded = home / "data" / "models" / "Immune_All_Low.pkl"
    pred = celltypist.annotate(sc.read(str(FIRST40)), model=str(downloaded))
    assert_annotation_matches(Path(ok["data"]["output_dir"]), *direct_tables(pred))
    print("downloaded Immune_All_Low.pkl sha256:", sha256(downloaded))


# ------------------------------------------------------------------------------- protocol: real stdio round trip
async def test_stdio_server_round_trip():
    """Launch the module as a stdio MCP server; any stdout pollution would corrupt the JSON-RPC stream."""
    home = seeded_cache_copy("stdio_cache")
    out = fresh_tmp("stdio_out")
    launcher = Path(__file__).resolve().parent / "_stdio_server.py"
    env = {k: v for k, v in os.environ.items()}
    env["CELLTYPIST_FOLDER"] = str(home)
    transport = StdioTransport(command=str(PROJECT_PYTHON), args=["-I", str(launcher)], env=env, cwd=str(out))
    async with Client(transport) as client:
        names = {t.name for t in await client.list_tools()}
        assert names == {"celltypist_list_models", "celltypist_annotate_cells", "celltypist_train_model", "celltypist_inspect_model"}
        res = (await client.call_tool(TOOL, {"adata_path": str(D400), "model": str(SEEDED), "majority_voting": True,
                                             "output_dir": str(out)})).data
    labels, dec, prob = reference_tables(REF / "annotate_demo400_seeded_rs0_mv")
    assert_annotation_matches(Path(res["output_dir"]), labels, dec, prob)
