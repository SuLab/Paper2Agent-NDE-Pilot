"""decoupler_rank_pcs_by_metadata vs tutorial cells 15-21 and direct scanpy/decoupler calls.

Tolerances: against a direct upstream run on the same saved input the wrapper must be bitwise
equal. Against the notebook reference, float64 variance ratios may differ by ~2e-16 (probe
tmp/verify_bulk_enrichment/probe_varratio.py: a direct run on the h5ad-reloaded, C-contiguous
filtered matrix gives the same 2.2e-16 offset as the wrapper; the in-memory non-contiguous
matrix of the notebook gives the reference bits), so atol=1e-12 is used there. X_pca and the
rank_obsm statistics derive from float32 PCA scores, so they are compared with rtol=1e-5 and
sign-invariantly (|stat|, per-column sign alignment of PC scores).
"""

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from conftest import (
    DER, FILTERED_H5AD, REF, Spy, artifact, assert_png, call_tool, call_tool_error, check_common_return,
    check_tool_format, sha256,
)

TOOL = "decoupler_rank_pcs_by_metadata"
KEYS = ["obsm", "obs"]


def direct_rank(path, obs_keys=None, cast_float64=False) -> ad.AnnData:
    """Tutorial cells 15 and 17 verbatim (including the counts layer and swap_layer)."""
    a = ad.read_h5ad(path)
    if cast_float64:
        a.X = a.X.astype(np.float64)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a)
    dc.pp.swap_layer(adata=a, key="counts", inplace=True)
    dc.tl.rankby_obsm(a, key="X_pca", obs_keys=obs_keys)
    return a


def sorted_assoc(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["obsm"], df["obs"] = df["obsm"].astype(str), df["obs"].astype(str)
    return df.sort_values(KEYS).reset_index(drop=True)


def assert_assoc_exact(actual: pd.DataFrame, expected: pd.DataFrame) -> None:
    a, e = sorted_assoc(actual), sorted_assoc(expected)
    assert a[KEYS].values.tolist() == e[KEYS].values.tolist()
    for c in ("stat", "pval", "padj"):
        # float32 columns are written by pandas in shortest float32 form; compare in that dtype.
        dtype = np.float32 if e[c].dtype == np.float32 else np.float64
        np.testing.assert_array_equal(a[c].to_numpy(dtype), e[c].to_numpy(dtype))


def assert_assoc_sign_invariant(actual: pd.DataFrame, expected: pd.DataFrame, rtol=1e-5, cols=("pval", "padj")) -> None:
    a, e = sorted_assoc(actual), sorted_assoc(expected)
    assert a[KEYS].values.tolist() == e[KEYS].values.tolist()
    np.testing.assert_allclose(np.abs(a["stat"].to_numpy(float)), np.abs(e["stat"].to_numpy(float)), rtol=rtol)
    for c in cols:
        np.testing.assert_allclose(a[c].to_numpy(float), e[c].to_numpy(float), rtol=rtol)


def assert_pcs_sign_invariant(actual: np.ndarray, expected: np.ndarray, rtol=1e-5, atol=1e-5) -> None:
    assert actual.shape == expected.shape
    for j in range(actual.shape[1]):
        sign = 1.0 if np.dot(actual[:, j], expected[:, j]) >= 0 else -1.0
        np.testing.assert_allclose(sign * actual[:, j], expected[:, j], rtol=rtol, atol=atol)


def read_wrapper_outputs(data):
    assoc = pd.read_csv(artifact(data, "pc_metadata_associations.csv"))
    var = pd.read_csv(artifact(data, "pca_variance_ratio.csv"))
    pcs = pd.read_csv(artifact(data, "pca_scores.csv"), index_col=0)
    return assoc, var, pcs


def test_tool_format():
    check_tool_format(TOOL, {"data_path"})


def test_reference_pca_and_associations(out_dir, monkeypatch):
    spy = Spy()
    spy.wrap(monkeypatch, dc.pl, "obsm")
    spy.wrap(monkeypatch, sc.pl, "pca")
    before = sha256(FILTERED_H5AD)
    data = call_tool(TOOL, {"data_path": str(FILTERED_H5AD), "output_dir": str(out_dir)})
    paths = check_common_return(data, out_dir, 6)
    for p in paths:
        if p.suffix == ".png":
            assert_png(p)
    assoc, var, pcs = read_wrapper_outputs(data)
    direct = direct_rank(FILTERED_H5AD)
    # Direct upstream on the same input: exact.
    assert_assoc_exact(assoc, direct.uns["rank_obsm"])
    np.testing.assert_array_equal(var["variance_ratio"].to_numpy(), direct.uns["pca"]["variance_ratio"])
    np.testing.assert_array_equal(var["variance"].to_numpy(), direct.uns["pca"]["variance"])
    np.testing.assert_array_equal(pcs.to_numpy(np.float32), direct.obsm["X_pca"])
    assert list(pcs.index) == list(direct.obs_names)
    # Notebook reference (cells 15/17): sign-invariant, float32 tolerance; variance ratio atol 1e-12.
    ref_assoc = pd.read_parquet(REF / "pca_rankby_obsm" / "rank_obsm.parquet")
    assert_assoc_sign_invariant(assoc, ref_assoc)
    ref_var = pd.read_csv(REF / "pca_rankby_obsm" / "pca_variance_ratio.csv")
    np.testing.assert_allclose(var["variance_ratio"], ref_var["variance_ratio"], rtol=0, atol=1e-12)
    ref_pcs = pd.read_csv(REF / "pca_rankby_obsm" / "X_pca.csv", index_col=0)
    assert list(ref_pcs.index) == list(pcs.index)
    assert_pcs_sign_invariant(pcs.to_numpy(np.float32), ref_pcs.to_numpy(np.float32))
    # Summary: only PC1 ~ condition is significant (padj 8.18e-4).
    assert data["n_pcs"] == 5 and data["n_tests"] == 10 and data["n_genes"] == 19510
    assert data["tested_obs_keys"] == ["condition", "sample_id"] == data["color_keys"]
    (sig,) = data["significant_associations"]
    assert (sig["obsm"], sig["obs"]) == ("PC1", "condition")
    assert sig["padj"] == pytest.approx(0.0008177597192116082, rel=1e-5)
    # Plot calls carry the tutorial settings (cells 19 and 21).
    (obsm_call,) = spy.calls["obsm"]
    assert obsm_call["nvar"] == 5 and obsm_call["titles"] == ["PC scores", "Adjusted p-values"]
    assert obsm_call["figsize"] == (6, 3)
    (pca_call,) = spy.calls["pca"]
    assert pca_call["color"] == ["condition", "sample_id"]
    assert (pca_call["ncols"], pca_call["size"], pca_call["frameon"]) == (1, 300, True)
    assert sha256(FILTERED_H5AD) == before


def test_relabelled_permuted_metadata_and_obs_keys(out_dir, monkeypatch):
    path = DER / "relabelled_permuted_filtered.h5ad"
    before = sha256(path)
    spy = Spy()
    spy.wrap(monkeypatch, dc.pl, "obsm")
    spy.wrap(monkeypatch, sc.pl, "pca")
    data = call_tool(TOOL, {"data_path": str(path), "obs_keys": ["treatment_status"],
                            "color_keys": ["treatment_status"], "nvar": 3, "output_dir": str(out_dir)})
    assoc, var, pcs = read_wrapper_outputs(data)
    direct = direct_rank(path, obs_keys=["treatment_status"])
    assert_assoc_exact(assoc, direct.uns["rank_obsm"])
    np.testing.assert_array_equal(var["variance_ratio"].to_numpy(), direct.uns["pca"]["variance_ratio"])
    assert data["n_tests"] == 5 and set(assoc["obs"]) == {"treatment_status"}
    # Relabelling and row permutation do not change the science: PC1 ~ labels as in the tutorial.
    ref_assoc = pd.read_parquet(REF / "pca_rankby_obsm" / "rank_obsm.parquet")
    ref_cond = ref_assoc[ref_assoc["obs"] == "condition"].assign(obs="treatment_status")
    # padj is not comparable: BH runs over 5 tests here vs 10 in the tutorial (obs_keys subset).
    assert_assoc_sign_invariant(assoc, ref_cond, cols=("pval",))
    ref_pcs = pd.read_csv(REF / "pca_rankby_obsm" / "X_pca.csv", index_col=0).loc[pcs.index]
    assert_pcs_sign_invariant(pcs.to_numpy(np.float32), ref_pcs.to_numpy(np.float32))
    assert spy.calls["obsm"][0]["nvar"] == 3
    assert spy.calls["pca"][0]["color"] == ["treatment_status"]
    assert sha256(path) == before


def test_integer_counts_are_processed_as_float64(out_dir):
    """Flagged choice: integer X is cast to float64 (the tutorial loader's dtype) before scanpy."""
    path = DER / "filtered_int64_counts.h5ad"
    data = call_tool(TOOL, {"data_path": str(path), "output_dir": str(out_dir)})
    assoc, var, _ = read_wrapper_outputs(data)
    direct64 = direct_rank(path, cast_float64=True)
    assert_assoc_exact(assoc, direct64.uns["rank_obsm"])
    np.testing.assert_array_equal(var["variance_ratio"].to_numpy(), direct64.uns["pca"]["variance_ratio"])
    # scanpy alone would normalise integer counts in float32; the difference is float32-sized only.
    direct_native = direct_rank(path)
    diff = np.max(np.abs(var["variance_ratio"].to_numpy() - direct_native.uns["pca"]["variance_ratio"]))
    assert diff < 1e-5


@pytest.mark.parametrize(
    "args, fragment",
    [
        ({"obs_keys": ["treatment_status"]}, "obs_keys"),
        ({"color_keys": ["batch"]}, "color_keys"),
        ({"data_path": "/nonexistent/x.h5ad"}, "does not exist"),
    ],
)
def test_errors(out_dir, args, fragment):
    full = {"data_path": str(FILTERED_H5AD), "output_dir": str(out_dir), **args}
    assert fragment in call_tool_error(TOOL, full)


def test_repeated_calls_fresh_dirs(out_dir):
    args = {"data_path": str(FILTERED_H5AD), "output_dir": str(out_dir)}
    a, b = call_tool(TOOL, args), call_tool(TOOL, args)
    pa, pb = check_common_return(a, out_dir, 6), check_common_return(b, out_dir, 6)
    assert pa[0].parent != pb[0].parent and not set(pa) & set(pb)
    pd.testing.assert_frame_equal(pd.read_csv(pa[0]), pd.read_csv(pb[0]))
