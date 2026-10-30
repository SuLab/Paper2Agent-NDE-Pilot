"""pydeseq2_test_differential_expression vs tutorial cells 23-26 and direct PyDESeq2 calls.

Tolerance: 0 (bitwise) everywhere a wrapper result is compared with a direct PyDESeq2 run or
with the notebook reference. Justification: the stage-2 executor showed bitwise-identical
results across 4 replays (n_cpus 1 vs 8, notebook-chain vs filtered-only input), and the CSV
written by the wrapper reloads float64 values exactly with pandas 3 read_csv. The swapped
contrast is checked as an exact sign flip: the contrast vector is negated, and IEEE rounding is
sign-symmetric, so LFC and stat flip exactly while lfcSE, pvalue and padj are unchanged.
"""

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import pytest
from pydeseq2.dds import DefaultInference, DeseqDataSet
from pydeseq2.ds import DeseqStats

from conftest import (
    DE_TIMEOUT, DER, FILTERED_H5AD, REF, Spy, artifact, assert_frame_exact, assert_png, call_tool,
    call_tool_error, check_common_return, check_tool_format, module_out_dir, read_de_csv, sha256,
)

TOOL = "pydeseq2_test_differential_expression"
COLS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]


def direct_fit(path, condition_key="condition", design=None, n_cpus=8):
    """Tutorial cell 23 verbatim (no quiet flag) on the saved input."""
    adata = ad.read_h5ad(path)
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(adata=adata, design=design or f"~{condition_key}", refit_cooks=True, inference=inference)
    dds.deseq2()
    return dds, inference


def direct_stats(dds, inference, contrast):
    stat_res = DeseqStats(dds, contrast=contrast, inference=inference)
    stat_res.summary()
    return stat_res.results_df


@pytest.fixture(scope="module")
def reference_call():
    out = module_out_dir("de_reference_call")
    with pytest.MonkeyPatch.context() as mp:
        spy = Spy()
        spy.wrap(mp, dc.pl, "volcano")
        data = call_tool(TOOL, {"data_path": str(FILTERED_H5AD), "condition_key": "condition", "group": "treatment",
                                "reference": "control", "n_cpus": 8, "output_dir": str(out)}, timeout=DE_TIMEOUT)
    return data, out, spy


@pytest.fixture(scope="module")
def reference_direct():
    return direct_fit(FILTERED_H5AD)


@pytest.fixture(scope="module")
def reference_df():
    return pd.read_parquet(REF / "deseq2" / "results_df.parquet")


def test_tool_format():
    check_tool_format(TOOL, {"data_path", "condition_key", "group", "reference"})


def test_reference_results_bitwise(reference_call, reference_df):
    data, out, _ = reference_call
    check_common_return(data, out, 2)
    res = read_de_csv(artifact(data, "deseq2_results.csv"))
    assert_frame_exact(res, reference_df)
    assert data["design"] == "~condition" and data["contrast"] == ["condition", "treatment", "control"]
    assert data["n_genes"] == 19510 and data["n_samples"] == 6
    assert data["samples_per_level"] == {"control": 3, "treatment": 3}
    assert (data["n_padj_lt_0.05"], data["n_padj_lt_0.05_up"], data["n_padj_lt_0.05_down"]) == (8517, 4328, 4189)
    assert data["n_missing"] == {c: 0 for c in COLS}
    expected_up = reference_df[(reference_df.padj < 0.05) & (reference_df.log2FoldChange > 0)].sort_values("padj")
    assert [r["gene"] for r in data["top_up"]] == list(expected_up.index[:5])


def test_direct_call_reproduces_reference(reference_direct, reference_df):
    dds, inference = reference_direct
    res = direct_stats(dds, inference, ["condition", "treatment", "control"])
    assert_frame_exact(res, reference_df)


def test_volcano_plot_data(reference_call, reference_df):
    data, _, spy = reference_call
    (call,) = spy.calls["volcano"]
    assert (call["x"], call["y"], call["top"]) == ("log2FoldChange", "padj", 10)
    assert "net" not in call
    assert_frame_exact(call["data"], reference_df)
    assert_png(artifact(data, "deseq2_volcano.png"))


def test_swapped_contrast(reference_direct, reference_df):
    out = module_out_dir("de_swapped")
    data = call_tool(TOOL, {"data_path": str(FILTERED_H5AD), "condition_key": "condition", "group": "control",
                            "reference": "treatment", "output_dir": str(out)}, timeout=DE_TIMEOUT)
    res = read_de_csv(artifact(data, "deseq2_results.csv"))
    dds, inference = reference_direct
    assert_frame_exact(res, direct_stats(dds, inference, ["condition", "control", "treatment"]))
    np.testing.assert_array_equal(res["log2FoldChange"], -reference_df["log2FoldChange"])
    np.testing.assert_array_equal(res["stat"], -reference_df["stat"])
    for c in ("baseMean", "lfcSE", "pvalue", "padj"):
        np.testing.assert_array_equal(res[c], reference_df[c])
    assert data["contrast"] == ["condition", "control", "treatment"]
    assert (data["n_padj_lt_0.05_up"], data["n_padj_lt_0.05_down"]) == (4189, 4328)


def test_relabelled_permuted_metadata():
    """Labels come from the user's column, not row order or hard-coded names; input untouched."""
    path = DER / "relabelled_permuted_filtered.h5ad"
    before = sha256(path)
    out = module_out_dir("de_relabelled")
    data = call_tool(TOOL, {"data_path": str(path), "condition_key": "treatment_status", "group": "TGFb_treated",
                            "reference": "vehicle", "output_dir": str(out)}, timeout=DE_TIMEOUT)
    res = read_de_csv(artifact(data, "deseq2_results.csv"))
    dds, inference = direct_fit(path, condition_key="treatment_status")
    expected = direct_stats(dds, inference, ["treatment_status", "TGFb_treated", "vehicle"])
    assert_frame_exact(res, expected)
    assert data["samples_per_level"] == {"TGFb_treated": 3, "vehicle": 3}
    assert data["design"] == "~treatment_status"
    assert sha256(path) == before
    # Scientific property: same samples and grouping as the tutorial (only labels, baseline level and
    # row order differ). PyDESeq2 itself is not bitwise invariant to these (attempt-2 inspection:
    # baseMean max |d| 7e-10, LFC max |d| 5.8e-4, lfcSE up to 0.13 for a few high-count genes, so a
    # Pearson r of stat 0.99990), hence rank/LFC/significant-set criteria. A row-order split
    # ({Ctrl1,TGFb1,Ctrl2} vs rest) would give unrelated statistics and fail all of them.
    ref = pd.read_parquet(REF / "deseq2" / "results_df.parquet")
    assert list(res.index) == list(ref.index)
    np.testing.assert_allclose(res["baseMean"], ref["baseMean"], rtol=1e-9)
    np.testing.assert_allclose(res["log2FoldChange"], ref["log2FoldChange"], rtol=0, atol=1e-3)
    assert res["stat"].corr(ref["stat"], method="spearman") > 0.9999
    sig_new, sig_ref = set(res.index[res.padj < 0.05]), set(ref.index[ref.padj < 0.05])
    assert len(sig_new ^ sig_ref) <= 0.005 * len(sig_ref)


def test_two_factor_design_reaches_pydeseq2():
    path = DER / "filtered_with_batch.h5ad"
    out = module_out_dir("de_two_factor")
    data = call_tool(TOOL, {"data_path": str(path), "condition_key": "condition", "group": "treatment",
                            "reference": "control", "design": "~batch + condition", "n_cpus": 4,
                            "output_dir": str(out)}, timeout=DE_TIMEOUT)
    res = read_de_csv(artifact(data, "deseq2_results.csv"))
    dds, inference = direct_fit(path, design="~batch + condition", n_cpus=8)
    assert list(dds.obsm["design_matrix"].columns) == ["Intercept", "batch[T.b2]", "batch[T.b3]", "condition[T.treatment]"]
    assert_frame_exact(res, direct_stats(dds, inference, ["condition", "treatment", "control"]))
    assert data["design"] == "~batch + condition" and data["n_cpus"] == 4
    ref = pd.read_parquet(REF / "deseq2" / "results_df.parquet")
    assert np.max(np.abs(res["stat"].to_numpy() - ref["stat"].to_numpy())) > 1e-2  # the covariate changed the fit


def test_quiet_and_copy_choices(capsys):
    """Flagged choices: quiet=True only silences stdout; PyDESeq2 mutates the AnnData it is given."""
    path = DER / "filtered_small_2000genes.h5ad"
    results = {}
    for quiet in (False, True):
        adata = ad.read_h5ad(path)
        inference = DefaultInference(n_cpus=8)
        dds = DeseqDataSet(adata=adata, design="~condition", refit_cooks=True, inference=inference, quiet=quiet)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["condition", "treatment", "control"], inference=inference, quiet=quiet)
        st.summary()
        printed = capsys.readouterr().out
        results[quiet] = st.results_df
        if quiet:
            assert printed == ""
        else:
            assert "Log2 fold change & Wald test p-value" in printed
        assert set(adata.uns) >= {"trend_coeffs", "disp_function_type"}  # input object mutated upstream
    assert_frame_exact(results[True], results[False])
    out = module_out_dir("de_small_wrapper")
    before = sha256(path)
    data = call_tool(TOOL, {"data_path": str(path), "condition_key": "condition", "group": "treatment",
                            "reference": "control", "output_dir": str(out)}, timeout=DE_TIMEOUT)
    assert capsys.readouterr().out == ""
    assert_frame_exact(read_de_csv(artifact(data, "deseq2_results.csv")), results[False])
    assert sha256(path) == before


@pytest.mark.parametrize(
    "args, fragment",
    [
        ({"condition_key": "treatment_status"}, "condition_key"),
        ({"group": "TGFb"}, "not found"),
        ({"group": "control"}, "different"),
        ({"data_path": str(DER / "non_integer_counts.h5ad")}, "non-negative integers"),
        ({"data_path": str(DER / "negative_counts.h5ad")}, "non-negative integers"),
        ({"data_path": str(DER / "numeric_condition.h5ad"), "condition_key": "cond_num", "group": "1",
          "reference": "0"}, "numeric"),
        ({"data_path": str(DER / "filtered_with_batch.h5ad"), "design": "~batch"}, "does not include condition_key"),
        ({"design": "~not_a_column + condition"}, "not_a_column"),
        ({"data_path": "/nonexistent/counts.h5ad"}, "does not exist"),
    ],
)
def test_errors(args, fragment):
    out = module_out_dir("de_errors")
    full = {"data_path": str(FILTERED_H5AD), "condition_key": "condition", "group": "treatment",
            "reference": "control", "output_dir": str(out), **args}
    assert fragment in call_tool_error(TOOL, full, timeout=DE_TIMEOUT)


def test_repeated_calls_fresh_dirs():
    out = module_out_dir("de_repeated")
    args = {"data_path": str(DER / "filtered_small_2000genes.h5ad"), "condition_key": "condition",
            "group": "treatment", "reference": "control", "output_dir": str(out)}
    a, b = call_tool(TOOL, args, timeout=DE_TIMEOUT), call_tool(TOOL, args, timeout=DE_TIMEOUT)
    pa, pb = check_common_return(a, out, 2), check_common_return(b, out, 2)
    assert pa[0].parent != pb[0].parent and not set(pa) & set(pb)
    assert sha256(pa[0]) == sha256(pb[0])
