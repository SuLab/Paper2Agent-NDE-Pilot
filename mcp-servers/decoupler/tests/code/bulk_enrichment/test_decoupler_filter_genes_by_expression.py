"""decoupler_filter_genes_by_expression vs tutorial cells 11/13 and direct dc.pp.filter_by_expr calls."""

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import pytest

from conftest import (
    DER, RAW_H5AD, REF, Spy, artifact, assert_png, call_tool, call_tool_error, check_common_return,
    check_tool_format, read_lines, sha256,
)

TOOL = "decoupler_filter_genes_by_expression"


def direct_filter(path, group, **params) -> ad.AnnData:
    a = ad.read_h5ad(path)
    dc.pp.filter_by_expr(adata=a, group=group, **params)
    return a


def test_tool_format():
    check_tool_format(TOOL, {"data_path", "group_key"})


def test_reference_filter(out_dir, monkeypatch):
    before = sha256(RAW_H5AD)
    spy = Spy()
    spy.wrap(monkeypatch, dc.pl, "filter_by_expr")
    data = call_tool(TOOL, {"data_path": str(RAW_H5AD), "group_key": "condition", "output_dir": str(out_dir)})
    check_common_return(data, out_dir, 2)
    res = ad.read_h5ad(artifact(data, "filtered_counts.h5ad"))
    # Exact gene set and order from the notebook (cell 13).
    kept = read_lines(REF / "filter_by_expr" / "kept_genes.txt")
    assert list(res.var_names) == kept and len(kept) == 19510
    ref = ad.read_h5ad(REF / "filter_by_expr" / "adata_filtered_by_expr.h5ad")
    np.testing.assert_array_equal(res.X, ref.X)
    pd.testing.assert_frame_equal(res.obs, ref.obs)
    # Independent direct call on the same saved input.
    direct = direct_filter(RAW_H5AD, "condition", min_count=10, min_total_count=15, large_n=10, min_prop=0.7)
    assert list(direct.var_names) == list(res.var_names)
    assert data["n_genes_before"] == 58674 and data["n_genes_after"] == 19510 and data["n_genes_removed"] == 39164
    assert data["n_samples"] == 6 and data["group_sizes"] == {"control": 3, "treatment": 3}
    # Diagnostic plot (cell 11) drawn on the unfiltered matrix with the same thresholds.
    (call,) = spy.calls["filter_by_expr"]
    assert call["adata"].n_vars == 58674
    assert {k: call[k] for k in ("group", "min_count", "min_total_count", "large_n", "min_prop")} == {
        "group": "condition", "min_count": 10, "min_total_count": 15, "large_n": 10, "min_prop": 0.7}
    assert_png(artifact(data, "filter_by_expr_thresholds.png"))
    assert sha256(RAW_H5AD) == before


def test_relabelled_permuted_metadata(out_dir):
    """User labels in a renamed column drive grouping; no 'condition' fallback; input untouched."""
    path = DER / "relabelled_permuted_raw.h5ad"
    before = sha256(path)
    data = call_tool(TOOL, {"data_path": str(path), "group_key": "treatment_status", "output_dir": str(out_dir)})
    res = ad.read_h5ad(artifact(data, "filtered_counts.h5ad"))
    direct = direct_filter(path, "treatment_status")
    assert list(res.var_names) == list(direct.var_names)
    np.testing.assert_array_equal(res.X, direct.X)
    src = ad.read_h5ad(path)
    pd.testing.assert_frame_equal(res.obs, src.obs)  # labels and interleaved order preserved
    assert list(res.obs["treatment_status"].astype(str)) == [
        "vehicle", "TGFb_treated", "vehicle", "TGFb_treated", "vehicle", "TGFb_treated"]
    assert data["group_sizes"] == {"TGFb_treated": 3, "vehicle": 3}
    # Same group sizes as the tutorial, so the same genes survive (gene order unchanged by row permutation).
    assert list(res.var_names) == read_lines(REF / "filter_by_expr" / "kept_genes.txt")
    assert sha256(path) == before
    err = call_tool_error(TOOL, {"data_path": str(path), "group_key": "condition", "output_dir": str(out_dir)})
    assert "group_key" in err and "condition" in err


def test_non_default_thresholds(out_dir, monkeypatch):
    params = {"min_count": 50, "min_total_count": 100, "large_n": 2, "min_prop": 0.5}
    spy = Spy()
    spy.wrap(monkeypatch, dc.pl, "filter_by_expr")
    data = call_tool(TOOL, {"data_path": str(RAW_H5AD), "group_key": "condition", "output_dir": str(out_dir), **params})
    res = ad.read_h5ad(artifact(data, "filtered_counts.h5ad"))
    direct = direct_filter(RAW_H5AD, "condition", **params)
    assert list(res.var_names) == list(direct.var_names)
    assert data["n_genes_after"] == direct.n_vars != 19510
    assert data["parameters"] == {"group_key": "condition", **params}
    (call,) = spy.calls["filter_by_expr"]
    assert {k: call[k] for k in params} == params


@pytest.mark.parametrize(
    "args, fragment",
    [
        ({"group_key": "treatment_status"}, "group_key"),
        ({"data_path": "/nonexistent/counts.h5ad", "group_key": "condition"}, "does not exist"),
    ],
)
def test_errors(out_dir, args, fragment):
    full = {"data_path": str(RAW_H5AD), "output_dir": str(out_dir), **args}
    assert fragment in call_tool_error(TOOL, full)


def test_repeated_calls_fresh_dirs(out_dir):
    args = {"data_path": str(RAW_H5AD), "group_key": "condition", "output_dir": str(out_dir)}
    a, b = call_tool(TOOL, args), call_tool(TOOL, args)
    pa, pb = check_common_return(a, out_dir, 2), check_common_return(b, out_dir, 2)
    assert pa[0].parent != pb[0].parent and not set(pa) & set(pb)
    ha, hb = ad.read_h5ad(pa[0]), ad.read_h5ad(pb[0])
    assert list(ha.var_names) == list(hb.var_names)
    np.testing.assert_array_equal(ha.X, hb.X)
