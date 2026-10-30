"""Verification of celltypist_inspect_model (tutorial cells 22-26 model meta information, 63-68 top driving genes + violins).

Expectations: executor references notebooks/celltypist_tutorial/reference/{immune_all_low_model, top_markers.json,
notebook_run/top_markers_mast_cells_3.json}, the tutorial figures images/figure_5.png (cell 67) and figure_6.png
(cell 68), and new direct Model.load / extract_top_markers / sc.pl.violin calls.

Violin figures: scanpy's strip-plot jitter is random, so two direct renders of the same sc.pl.violin call differ in about
0.27 % of pixels (measured on this machine; tool and tutorial figure_5 differ by the same amount). Figures are compared by
exact image size plus a differing-pixel fraction below 1 % (channel difference > 0.1).
"""
import inspect
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from matplotlib import image as mpimg
from matplotlib import pyplot as plt

from celltypist import models as ct_models
from fastmcp import Client

from _celltypist_tutorial_support import (CACHE_MODELS, INPUTS, MODELS, NB, REF, assert_artifacts, call_tool,
                                          call_tool_error, celltypist_tutorial_mcp, fresh_tmp, sha256)

TOOL = "celltypist_inspect_model"
D2000 = INPUTS / "demo_2000_cells.h5ad"
D400 = INPUTS / "demo_400_cells.h5ad"
SEEDED = MODELS / "model_seeded_rs0_run1.pkl"
ASRUN = MODELS / "model_from_immune2000.pkl"
IAL_CACHED = CACHE_MODELS / "Immune_All_Low.pkl"
PIXEL_FRACTION_MAX = 0.01


def differing_pixel_fraction(a: Path, b: Path) -> float:
    x, y = mpimg.imread(a), mpimg.imread(b)
    assert x.shape == y.shape, (x.shape, y.shape)
    return float((np.abs(x - y).max(axis=2) > 0.1).mean())


def markers_csv(path) -> dict:
    df = pd.read_csv(path, dtype={"cell_type": str, "gene": str}, keep_default_na=False)
    assert list(df.columns) == ["cell_type", "rank", "gene"]
    out = {}
    for ct, sub in df.groupby("cell_type", sort=False):
        assert sub["rank"].tolist() == list(range(1, len(sub) + 1))
        out[ct] = sub["gene"].tolist()
    return out


async def test_schema_defaults_follow_upstream():
    async with Client(celltypist_tutorial_mcp) as client:
        tool = {t.name: t for t in await client.list_tools()}[TOOL]
    schema = tool.input_schema
    props = schema["properties"]
    assert set(schema["required"]) == {"model"}
    assert set(props) == {"model", "cell_types", "top_n", "only_positive", "adata_path", "groupby_key", "output_dir"}
    assert all(p.get("description") for p in props.values())
    up = inspect.signature(ct_models.Model.extract_top_markers).parameters
    assert props["top_n"]["default"] == up["top_n"].default == 10
    assert props["only_positive"]["default"] is up["only_positive"].default is True
    assert len(tool.description.strip().splitlines()) == 2


async def test_builtin_model_metadata_and_all_markers(out_base):
    """Tutorial cells 22-26 (Model.load('Immune_All_Low.pkl'), repr, cell_types) plus default top-10 markers of every type."""
    result = await call_tool(TOOL, {"model": "Immune_All_Low.pkl", "output_dir": str(out_base)})
    (csv,) = assert_artifacts(result, out_base, 1)
    ref_types = pd.read_csv(REF / "immune_all_low_model" / "cell_types.csv", dtype=str)["cell_type"].tolist()
    ref_summary = json.loads((REF / "immune_all_low_model" / "summary.json").read_text())
    assert result["cell_types"] == ref_types and result["n_cell_types"] == 98 == ref_summary["n_cell_types"]
    assert result["n_features"] == 6639 == ref_summary["n_features"]
    assert result["description"] == {k: str(v) for k, v in ref_summary["description"].items()}
    assert result["model"] == "Immune_All_Low.pkl" and Path(result["model_file"]) == IAL_CACHED

    direct = ct_models.Model.load(str(IAL_CACHED))
    expected = {str(ct): [str(g) for g in direct.extract_top_markers(ct, 10, True)] for ct in direct.cell_types}
    assert result["top_markers"] == expected
    assert markers_csv(csv) == expected and len(pd.read_csv(csv)) == 980


async def test_seeded_model_top3_all_types_match_executor_reference(out_base):
    result = await call_tool(TOOL, {"model": str(SEEDED), "top_n": 3, "output_dir": str(out_base)})
    ref = json.loads((REF / "top_markers.json").read_text())["seeded_rs0_run1"]
    assert result["top_markers"] == ref["all_cell_types_top3"]
    assert result["top_markers"]["Mast cells"] == ref["Mast cells top3 (extract_top_markers('Mast cells', 3))"] == ["TPSB2", "KRT1", "SLC45A3"]
    assert markers_csv(result["artifacts"][0]["path"]) == ref["all_cell_types_top3"]
    assert result["n_features"] == 2652 and result["n_cell_types"] == 10


async def test_nondefault_top_n_and_absolute_ranking(out_base):
    result = await call_tool(TOOL, {"model": str(SEEDED), "cell_types": ["Plasma cells", "pDC"], "top_n": 7,
                                    "only_positive": False, "output_dir": str(out_base)})
    direct = ct_models.Model.load(str(SEEDED))
    expected = {ct: [str(g) for g in direct.extract_top_markers(ct, 7, False)] for ct in ["Plasma cells", "pDC"]}
    assert result["top_markers"] == expected
    # only_positive is honoured: for Plasma cells the absolute-coefficient top 7 differs from the positive top 7
    # (for Mast cells / pDC both rankings coincide at top_n=7, so they cannot reveal an ignored flag).
    positive = [str(g) for g in direct.extract_top_markers("Plasma cells", 7, True)]
    assert expected["Plasma cells"] != positive
    assert markers_csv(result["artifacts"][0]["path"]) == expected


async def test_tutorial_cell67_violin_as_run_model_on_training_data(out_base):
    """Cells 63-67: as-run model, extract_top_markers('Mast cells', 3), sc.pl.violin(adata_2000, genes, groupby='cell_type', rotation=90)."""
    before = sha256(D2000)
    result = await call_tool(TOOL, {"model": str(ASRUN), "cell_types": ["Mast cells"], "top_n": 3, "adata_path": str(D2000),
                                    "groupby_key": "cell_type", "output_dir": str(out_base)})
    paths = assert_artifacts(result, out_base, 2)
    nb = json.loads((REF / "notebook_run" / "top_markers_mast_cells_3.json").read_text())["result"]
    assert result["top_markers"] == {"Mast cells": nb} and nb == ["TPSB2", "SLC45A3", "CPA3"]
    assert result["genes_missing_from_adata"] == {} and result["groupby_key"] == "cell_type"
    png = paths[1]
    assert png.name == "violin_01_Mast_cells.png"
    assert differing_pixel_fraction(png, NB / "images" / "figure_5.png") < PIXEL_FRACTION_MAX

    # Same upstream call made directly here, saved the same way.
    direct_png = fresh_tmp("direct_violin") / "direct.png"
    sc.pl.violin(sc.read(str(D2000)), nb, groupby="cell_type", rotation=90, show=False)
    axes = plt.gcf().axes
    assert [ax.get_ylabel() for ax in axes] == nb
    plt.savefig(direct_png, bbox_inches="tight")
    plt.close("all")
    assert differing_pixel_fraction(png, direct_png) < PIXEL_FRACTION_MAX
    assert sha256(D2000) == before


async def test_tutorial_cell68_violin_on_annotated_query_by_majority_voting(out_base):
    """annotate -> inspect composition: cell 56 (as-run model on demo_400, majority voting) then cell 68 violin by majority_voting."""
    ann = await call_tool("celltypist_annotate_cells", {"adata_path": str(D400), "model": str(ASRUN), "majority_voting": True,
                                                         "output_dir": str(out_base)})
    # majority-voting labels equal the notebook's own cell-56 run (replayed by the executor)
    got_mv = pd.read_csv(Path(ann["output_dir"]) / "predicted_labels.csv", index_col=0, dtype=str)["majority_voting"]
    ref_mv = pd.read_csv(REF / "notebook_run" / "demo400_custom_asrun_mv" / "predicted_labels.csv", index_col=0, dtype=str)["majority_voting"]
    assert got_mv.equals(ref_mv)
    h5 = Path(ann["output_dir"]) / "annotated.h5ad"
    before = sha256(h5)
    result = await call_tool(TOOL, {"model": str(ASRUN), "cell_types": ["Mast cells"], "top_n": 3, "adata_path": str(h5),
                                    "groupby_key": "majority_voting", "output_dir": str(out_base)})
    paths = assert_artifacts(result, out_base, 2)
    assert differing_pixel_fraction(paths[1], NB / "images" / "figure_6.png") < PIXEL_FRACTION_MAX
    assert sha256(h5) == before


async def test_genes_absent_from_adata_are_reported_and_plot_uses_the_rest(out_base):
    direct = ct_models.Model.load(str(SEEDED))
    genes = [str(g) for g in direct.extract_top_markers("pDC", 3)]
    sub = ad.read_h5ad(D400)
    sub = sub[:, [g for g in sub.var_names if g != genes[1]]].copy()
    path = fresh_tmp("gene_subset") / "demo_400_minus_one_gene.h5ad"
    sub.write_h5ad(path)
    result = await call_tool(TOOL, {"model": str(SEEDED), "cell_types": ["pDC", "Mast cells"], "top_n": 3,
                                    "adata_path": str(path), "groupby_key": "cell_type", "output_dir": str(out_base)})
    assert result["genes_missing_from_adata"] == {"pDC": [genes[1]]}
    assert [Path(a["path"]).name for a in result["artifacts"]] == ["top_markers.csv", "violin_01_pDC.png", "violin_02_Mast_cells.png"]
    assert result["top_markers"]["pDC"] == genes        # the reported markers are unaffected by the user's gene set


async def test_repeated_calls_write_fresh_directories(out_base):
    args = {"model": str(SEEDED), "top_n": 4, "output_dir": str(out_base)}
    r1, r2 = await call_tool(TOOL, args), await call_tool(TOOL, args)
    assert r1["output_dir"] != r2["output_dir"]
    assert r1["artifacts"][0]["path"] != r2["artifacts"][0]["path"]
    assert Path(r1["artifacts"][0]["path"]).read_bytes() == Path(r2["artifacts"][0]["path"]).read_bytes()


# ------------------------------------------------------------------------------------------------------ errors
async def test_unknown_cell_type_surfaces_upstream_error(out_base):
    with pytest.raises(ValueError) as direct:
        ct_models.Model.load(str(SEEDED)).extract_top_markers("Unicorn cells", 3)
    msg = await call_tool_error(TOOL, {"model": str(SEEDED), "cell_types": ["Unicorn cells"], "output_dir": str(out_base)})
    assert str(direct.value) in msg and "'Unicorn cells' is not found" in msg
    assert list(out_base.iterdir()) == []


async def test_unknown_builtin_model_name_without_downloads(out_base, monkeypatch):
    urls = []

    def blocked(url, timeout=30):
        urls.append(url)
        raise Exception("blocked in test")

    monkeypatch.setattr(ct_models, "_requests_get", blocked)
    msg = await call_tool_error(TOOL, {"model": "No_Such_Model.pkl", "output_dir": str(out_base)})
    assert "No models match the celltypist model repertoire" in msg and urls == []


async def test_input_errors(out_base):
    cases = [
        ({"model": str(MODELS / "absent.pkl")}, "No such model file"),
        ({"model": str(SEEDED), "top_n": 0}, "top_n must be a positive integer"),
        ({"model": str(SEEDED), "cell_types": ["Mast cells"], "adata_path": str(D400)}, "groupby_key is required"),
        ({"model": str(SEEDED), "adata_path": str(D400), "groupby_key": "cell_type"}, "cell_types must list the model cell types"),
        ({"model": str(SEEDED), "cell_types": ["Mast cells"], "adata_path": str(INPUTS / "absent.h5ad"), "groupby_key": "cell_type"},
         "No such file"),
        ({"model": str(SEEDED), "cell_types": ["Mast cells"], "adata_path": str(INPUTS / "demo_400_user_meta.h5ad"),
          "groupby_key": "cell_type"}, "groupby_key 'cell_type' is not an obs column"),
    ]
    for args, fragment in cases:
        msg = await call_tool_error(TOOL, {**args, "output_dir": str(out_base)})
        assert fragment in msg, (args, msg)
    assert list(out_base.iterdir()) == []
