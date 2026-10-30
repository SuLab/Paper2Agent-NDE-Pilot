"""decoupler_score_activities_ulm vs tutorial cells 28-34, 40-46, 52-58 and direct dc.mt.ulm calls.

Exact comparisons use the executor's network snapshots through network_path. ulm is deterministic
(executor replays were bitwise identical) and the wrapper's CSVs reload float64 values exactly, so
tolerance is 0. Live-download checks are marked `network` and only test structure and rough
agreement, because CollecTRI/PROGENy/Hallmark are remote resources that can drift.
"""

import json

import decoupler as dc
import numpy as np
import pandas as pd
import pytest

from conftest import (
    CONTRAST, DE_CSV, DER, NETS, REF, Spy, artifact, assert_frame_exact, assert_png, call_tool, call_tool_error,
    check_common_return, check_tool_format, read_de_csv,
)

TOOL = "decoupler_score_activities_ulm"
EXPECT = {
    "collectri": {"n_src": 655, "n_sig": 143, "pos": 64, "neg": 79, "rows": 42990, "top_pos": "SRF", "top_neg": "IRF1",
                  "bar_top": 25, "bar_figsize": (3, 5), "dot_scale": 0.15, "dot_figsize": (4, 4)},
    "progeny": {"n_src": 14, "n_sig": 12, "pos": 8, "neg": 4, "rows": 62421, "top_pos": "TGFb", "top_neg": "JAK-STAT",
                "bar_top": 25, "bar_figsize": (3, 3), "dot_scale": 0.15, "dot_figsize": (4, 4)},
    "hallmark": {"n_src": 50, "n_sig": 30, "pos": 13, "neg": 17, "rows": 7318,
                 "top_pos": "EPITHELIAL_MESENCHYMAL_TRANSITION", "top_neg": "INTERFERON_ALPHA_RESPONSE",
                 "bar_top": 10, "bar_figsize": (6, 3), "dot_scale": 0.05, "dot_figsize": (6.5, 4)},
}
LIVE_URLS = {
    "collectri": "https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1",
    "progeny": "https://omnipathdb.org/",
    "hallmark": "https://static.omnipathdb.org/tables/msigdb-hallmark.tsv.gz",
}


def ref_data() -> pd.DataFrame:
    return pd.read_parquet(REF / "ulm" / "data_stat_matrix.parquet")


def snapshot(net: str) -> pd.DataFrame:
    return pd.read_parquet(NETS / f"{net}.parquet")


def dotplot_frame(acts_sig: pd.DataFrame, padj: pd.DataFrame) -> pd.DataFrame:
    """Tutorial cell 46 transformation, applied directly to reference objects."""
    return acts_sig.melt(value_name="score").merge(
        padj.melt(value_name="padj")
        .assign(padj=lambda x: x["padj"].clip(2.22e-16, 1))
        .assign(padj=lambda x: -np.log10(x["padj"]))
    )


def read_tables(data, net):
    scores = pd.read_csv(artifact(data, f"{net}_ulm_scores.csv"), index_col=0)
    padj = pd.read_csv(artifact(data, f"{net}_ulm_padj.csv"), index_col=0)
    sig = pd.read_csv(artifact(data, f"{net}_ulm_scores_significant.csv"), index_col=0)
    return scores, padj, sig


def test_tool_format():
    check_tool_format(TOOL, {"de_results_path", "network"})


@pytest.mark.parametrize("net", ["collectri", "progeny", "hallmark"])
def test_reference_scores(net, out_dir, monkeypatch):
    spy = Spy()
    spy.wrap(monkeypatch, dc.pl, "barplot")
    spy.wrap(monkeypatch, dc.pl, "dotplot")
    data = call_tool(TOOL, {"de_results_path": str(DE_CSV), "network": net, "network_path": str(NETS / f"{net}.parquet"),
                            "contrast_name": CONTRAST, "output_dir": str(out_dir)})
    paths = check_common_return(data, out_dir, 6)
    scores, padj, sig = read_tables(data, net)
    # Notebook reference (full scores from the driver's identical call; padj and filtered acts from the kernel).
    assert_frame_exact(scores, pd.read_parquet(REF / "ulm" / f"{net}_scores_full.parquet"))
    assert_frame_exact(padj, pd.read_parquet(REF / "ulm" / f"{net}_padj_full.parquet"))
    ref_sig = pd.read_parquet(REF / "ulm" / f"{net}_scores_padj_lt_0.05.parquet")
    assert_frame_exact(sig, ref_sig)
    # Fresh direct upstream call on the saved stat matrix and snapshot.
    acts, pv = dc.mt.ulm(data=ref_data(), net=snapshot(net))
    assert_frame_exact(scores, acts)
    assert_frame_exact(padj, pv)
    e = EXPECT[net]
    assert sig.shape[1] == e["n_sig"] == data["n_significant"]
    assert (data["n_significant_positive"], data["n_significant_negative"]) == (e["pos"], e["neg"])
    assert data["n_sources_scored"] == e["n_src"] and data["network_rows"] == e["rows"]
    assert data["contrast_row"] == CONTRAST and data["n_genes"] == 19510
    assert data["network_origin"] == str((NETS / f"{net}.parquet").resolve())
    assert data["top_significant_sources"]["positive"][0]["source"] == e["top_pos"]
    assert data["top_significant_sources"]["negative"][0]["source"] == e["top_neg"]
    used = pd.read_csv(artifact(data, f"{net}_network_used.csv"))
    snap = snapshot(net)
    assert list(used.columns) == list(snap.columns) and len(used) == len(snap)
    assert used["source"].astype(str).tolist() == snap["source"].astype(str).tolist()
    # Bar plot (cells 34/44/56): significant scores, tutorial top and figsize.
    (bar,) = spy.calls["barplot"]
    assert_frame_exact(bar["data"], ref_sig)
    assert (bar["name"], bar["top"], bar["figsize"]) == (CONTRAST, e["bar_top"], e["bar_figsize"])
    # Dot plot (cells 46/58; CollecTRI has no tutorial dot plot and uses the cell-46 settings).
    (dot,) = spy.calls["dotplot"]
    assert (dot["x"], dot["y"], dot["s"], dot["c"], dot["top"]) == ("score", "variable", "padj", "score", 10)
    assert (dot["scale"], dot["figsize"]) == (e["dot_scale"], e["dot_figsize"])
    got = dot["df"][["variable", "score", "padj"]].sort_values("variable").reset_index(drop=True)
    if net == "collectri":
        exp = dotplot_frame(ref_sig, pd.read_parquet(REF / "ulm" / "collectri_padj_full.parquet"))
    else:
        exp = pd.read_csv(REF / "ulm" / f"{net}_dotplot_df.csv")
    exp = exp[["variable", "score", "padj"]].sort_values("variable").reset_index(drop=True)
    assert got["variable"].astype(str).tolist() == exp["variable"].astype(str).tolist()
    np.testing.assert_array_equal(got[["score", "padj"]].to_numpy(float), exp[["score", "padj"]].to_numpy(float))
    for p in paths:
        if p.suffix == ".png":
            assert_png(p)


def test_padj_threshold_non_default(out_dir, monkeypatch):
    spy = Spy()
    spy.wrap(monkeypatch, dc.pl, "barplot")
    data = call_tool(TOOL, {"de_results_path": str(DE_CSV), "network": "collectri", "padj_threshold": 0.001,
                            "network_path": str(NETS / "collectri.parquet"), "contrast_name": CONTRAST,
                            "output_dir": str(out_dir)})
    _, _, sig = read_tables(data, "collectri")
    acts, pv = dc.mt.ulm(data=ref_data(), net=snapshot("collectri"))
    msk = (pv.T < 0.001).iloc[:, 0]
    assert_frame_exact(sig, acts.loc[:, msk])
    assert data["n_significant"] == int(msk.sum()) < 143
    assert_frame_exact(spy.calls["barplot"][0]["data"], acts.loc[:, msk])


def test_stat_column_log2foldchange(out_dir):
    data = call_tool(TOOL, {"de_results_path": str(DE_CSV), "network": "progeny", "stat_column": "log2FoldChange",
                            "network_path": str(NETS / "progeny.parquet"), "output_dir": str(out_dir)})
    scores, padj, _ = read_tables(data, "progeny")
    de = read_de_csv(DE_CSV)
    acts, pv = dc.mt.ulm(data=de[["log2FoldChange"]].T, net=snapshot("progeny"))
    assert_frame_exact(scores, acts)
    assert_frame_exact(padj, pv)
    assert data["contrast_row"] == "log2FoldChange" == scores.index[0]
    ref = pd.read_parquet(REF / "ulm" / "progeny_scores_full.parquet")
    assert not np.allclose(scores.to_numpy(), ref.to_numpy())


def test_default_contrast_name_keeps_values(out_dir):
    data = call_tool(TOOL, {"de_results_path": str(DE_CSV), "network": "hallmark",
                            "network_path": str(NETS / "hallmark.parquet"), "output_dir": str(out_dir)})
    scores, _, _ = read_tables(data, "hallmark")
    assert list(scores.index) == ["stat"]
    ref = pd.read_parquet(REF / "ulm" / "hallmark_scores_full.parquet")
    np.testing.assert_array_equal(scores.to_numpy(), ref.to_numpy())


def test_plot_top_parameters(out_dir, monkeypatch):
    spy = Spy()
    spy.wrap(monkeypatch, dc.pl, "barplot")
    spy.wrap(monkeypatch, dc.pl, "dotplot")
    call_tool(TOOL, {"de_results_path": str(DE_CSV), "network": "hallmark", "barplot_top": 5, "dotplot_top": 3,
                     "network_path": str(NETS / "hallmark.parquet"), "output_dir": str(out_dir)})
    assert spy.calls["barplot"][0]["top"] == 5 and spy.calls["dotplot"][0]["top"] == 3


@pytest.mark.parametrize("path", [NETS / "hallmark.csv", DER / "hallmark.tsv"])
def test_network_file_formats(path, out_dir):
    data = call_tool(TOOL, {"de_results_path": str(DE_CSV), "network": "hallmark", "network_path": str(path),
                            "contrast_name": CONTRAST, "output_dir": str(out_dir)})
    scores, padj, _ = read_tables(data, "hallmark")
    assert_frame_exact(scores, pd.read_parquet(REF / "ulm" / "hallmark_scores_full.parquet"))
    assert_frame_exact(padj, pd.read_parquet(REF / "ulm" / "hallmark_padj_full.parquet"))


@pytest.mark.parametrize(
    "args, fragment",
    [
        ({"de_results_path": str(DER / "deseq2_results_no_stat.csv")}, "lack column"),
        ({"de_results_path": str(DER / "deseq2_results_nonfinite_stat.csv")}, "non-finite"),
        ({"de_results_path": "/nonexistent/de.csv"}, "does not exist"),
        ({"padj_threshold": 0.0}, "padj_threshold"),
        ({"padj_threshold": 1.5}, "padj_threshold"),
        ({"network_path": "/nonexistent/net.parquet"}, "does not exist"),
        ({"network_path": str(DER / "network_bad_columns.csv")}, "'source' and 'target'"),
        ({"network": "dorothea"}, "network"),
    ],
)
def test_errors(args, fragment, out_dir):
    full = {"de_results_path": str(DE_CSV), "network": "hallmark", "network_path": str(NETS / "hallmark.parquet"),
            "output_dir": str(out_dir), **args}
    assert fragment in call_tool_error(TOOL, full)


def test_repeated_calls_fresh_dirs(out_dir):
    args = {"de_results_path": str(DE_CSV), "network": "progeny", "network_path": str(NETS / "progeny.parquet"),
            "contrast_name": CONTRAST, "output_dir": str(out_dir)}
    a, b = call_tool(TOOL, args), call_tool(TOOL, args)
    pa, pb = check_common_return(a, out_dir, 6), check_common_return(b, out_dir, 6)
    assert pa[0].parent != pb[0].parent and not set(pa) & set(pb)
    assert pa[0].read_bytes() == pb[0].read_bytes()


@pytest.mark.network
@pytest.mark.parametrize("net", ["collectri", "progeny", "hallmark"])
def test_live_download_rough_agreement(net, out_dir):
    import requests

    try:
        requests.head(LIVE_URLS[net], timeout=15, allow_redirects=True).raise_for_status()
    except Exception as exc:  # pragma: no cover - depends on connectivity
        pytest.skip(f"network resource unreachable: {exc}")
    data = call_tool(TOOL, {"de_results_path": str(DE_CSV), "network": net, "contrast_name": CONTRAST,
                            "output_dir": str(out_dir)}, timeout=600)
    assert data["network_origin"] == f"downloaded with decoupler.op.{net}(organism='human')"
    used = pd.read_csv(artifact(data, f"{net}_network_used.csv"))
    snap = snapshot(net)
    assert {"source", "target"} <= set(used.columns)
    assert ("weight" in used.columns) == ("weight" in snap.columns)
    assert 0.8 <= len(used) / len(snap) <= 1.25
    scores, _, _ = read_tables(data, net)
    ref = pd.read_parquet(REF / "ulm" / f"{net}_scores_full.parquet")
    common = scores.columns.intersection(ref.columns)
    frac = len(common) / ref.shape[1]
    r = float(np.corrcoef(scores[common].iloc[0], ref[common].iloc[0])[0, 1])
    summary = {"network": net, "rows_live": int(len(used)), "rows_snapshot": int(len(snap)),
               "sources_live": int(scores.shape[1]), "sources_snapshot": int(ref.shape[1]),
               "common_source_fraction": frac, "pearson_r_common": r, "n_significant_live": data["n_significant"],
               "n_significant_snapshot": EXPECT[net]["n_sig"]}
    (out_dir / "live_vs_snapshot.json").write_text(json.dumps(summary, indent=2))
    assert frac >= 0.9 and r >= 0.95, summary
