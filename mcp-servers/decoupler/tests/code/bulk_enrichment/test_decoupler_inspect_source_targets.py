"""decoupler_inspect_source_targets vs tutorial cells 36, 38, 48, 50, 60 and direct dc.pl calls.

Leading edges are compared as exact ordered gene lists with the kernel-captured references
(cells 50 and 60) and with fresh direct dc.pl.leading_edge calls. The volcano (cell 38) and
source_targets (cell 48) figures have adjustText label jitter, so the data and arguments handed
to the upstream plotters are compared instead of pixels.
"""

import decoupler as dc
import pandas as pd
import pytest

from conftest import (
    CONTRAST, DE_CSV, DER, NETS, REF, SIG_SCORES_CSV, Spy, artifact, assert_frame_exact, assert_png, call_tool,
    call_tool_error, check_common_return, check_tool_format, direct_leading_edge, read_de_csv, read_lines,
)

TOOL = "decoupler_inspect_source_targets"
CELL36_SOURCES = ["SRF", "MYOCD", "IRF1", "SKIL"]


def snapshot(net: str) -> pd.DataFrame:
    return pd.read_parquet(NETS / f"{net}.parquet")


def spies(monkeypatch) -> Spy:
    spy = Spy()
    for name in ("volcano", "source_targets", "leading_edge", "network"):
        spy.wrap(monkeypatch, dc.pl, name)
    return spy


def call(net, source, out_dir, **extra):
    return call_tool(TOOL, {"de_results_path": str(DE_CSV), "network": net, "source": source,
                            "network_path": str(NETS / f"{net}.parquet"), "output_dir": str(out_dir), **extra})


def le_files(data) -> dict:
    return {label: read_lines(info["path"]) for label, info in data["leading_edges"].items()}


def test_tool_format():
    check_tool_format(TOOL, {"de_results_path", "network", "source"})


def test_progeny_tgfb_cells_48_50(out_dir, monkeypatch):
    spy = spies(monkeypatch)
    data = call("progeny", "TGFb", out_dir)
    paths = check_common_return(data, out_dir, 6)
    de, net = read_de_csv(DE_CSV), snapshot("progeny")
    got = le_files(data)
    assert got["positive_weights"] == read_lines(REF / "leading_edge" / "progeny_TGFb_positive_weights.txt")
    assert got["negative_weights"] == read_lines(REF / "leading_edge" / "progeny_TGFb_negative_weights.txt")
    assert (len(got["positive_weights"]), len(got["negative_weights"])) == (581, 292)
    assert got == direct_leading_edge(de, net, "TGFb")
    # Cell 48 arguments and data.
    (st,) = spy.calls["source_targets"]
    assert (st["x"], st["y"], st["name"], st["top"], st["figsize"]) == ("weight", "stat", "TGFb", 15, (4, 4))
    assert_frame_exact(st["data"], de)
    pd.testing.assert_frame_equal(st["net"], net)
    # Leading-edge calls use the weight-sign split of cell 50.
    le_nets = [c["net"] for c in spy.calls["leading_edge"]]
    pd.testing.assert_frame_equal(le_nets[0], net[net["weight"] > 0])
    pd.testing.assert_frame_equal(le_nets[1], net[net["weight"] < 0])
    assert all(c["stat"] == "stat" and c["name"] == "TGFb" for c in spy.calls["leading_edge"])
    (vol,) = spy.calls["volcano"]
    assert (vol["x"], vol["y"], vol["name"], vol["top"]) == ("log2FoldChange", "pvalue", "TGFb", 10)
    assert "network" not in spy.calls or not spy.calls["network"]
    assert any("network plot" in s for s in data["skipped"])
    tgfb = net[net["source"] == "TGFb"]
    assert data["n_targets_in_network"] == len(tgfb)
    assert data["n_targets_measured"] == int(tgfb["target"].isin(de.index).sum())
    assert data["weighted_network"] is True
    for p in paths:
        if p.suffix == ".png":
            assert_png(p)


def test_hallmark_tgf_beta_cell_60(out_dir, monkeypatch):
    spy = spies(monkeypatch)
    data = call("hallmark", "TGF_BETA_SIGNALING", out_dir)
    got = le_files(data)
    assert list(got) == ["all_targets"]
    assert got["all_targets"] == read_lines(REF / "leading_edge" / "hallmark_TGF_BETA_SIGNALING.txt")
    assert len(got["all_targets"]) == 19
    pd.testing.assert_frame_equal(spy.calls["leading_edge"][0]["net"], snapshot("hallmark"))
    assert "source_targets" not in spy.calls or not spy.calls["source_targets"]
    assert any("source_targets" in s for s in data["skipped"]) and data["weighted_network"] is False


def test_collectri_srf_network_and_volcano_cells_36_38(out_dir, monkeypatch):
    spy = spies(monkeypatch)
    data = call("collectri", "SRF", out_dir, scores_path=str(SIG_SCORES_CSV), network_sources=CELL36_SOURCES)
    de, net = read_de_csv(DE_CSV), snapshot("collectri")
    # Cell 38 volcano.
    (vol,) = spy.calls["volcano"]
    assert (vol["x"], vol["y"], vol["name"], vol["top"]) == ("log2FoldChange", "pvalue", "SRF", 10)
    assert_frame_exact(vol["data"], de)
    pd.testing.assert_frame_equal(vol["net"], net)
    # Cell 36 network: data = the cell-28 stat matrix, score = the cell-32 significant tf_acts.
    (nw,) = spy.calls["network"]
    assert_frame_exact(nw["data"], pd.read_parquet(REF / "ulm" / "data_stat_matrix.parquet"))
    assert_frame_exact(nw["score"], pd.read_parquet(REF / "ulm" / "collectri_scores_padj_lt_0.05.parquet"))
    assert nw["data"].index[0] == CONTRAST
    pd.testing.assert_frame_equal(nw["net"], net)
    assert nw["sources"] == CELL36_SOURCES
    assert (nw["targets"], nw["figsize"], nw["vcenter"], nw["by_abs"], nw["size_node"], nw["t_label"]) == (
        10, (5, 5), True, True, 10, "Gene-level\nstatistic")
    assert_png(artifact(data, "SRF_network.png"))
    # Leading edge with the same weight-sign split (not demonstrated for CollecTRI in the tutorial).
    assert le_files(data) == direct_leading_edge(de, net, "SRF")
    assert data["source_score"] == pytest.approx(7.383742530399846, rel=1e-12)
    assert data["network_plot_sources"] == CELL36_SOURCES


@pytest.mark.parametrize("net, source", [("progeny", "JAK-STAT"), ("hallmark", "EPITHELIAL_MESENCHYMAL_TRANSITION"),
                                         ("collectri", "IRF1")])
def test_other_sources(net, source, out_dir):
    data = call(net, source, out_dir)
    expected = direct_leading_edge(read_de_csv(DE_CSV), snapshot(net), source)
    assert le_files(data) == expected and expected
    assert {k: v["n_genes"] for k, v in data["leading_edges"].items()} == {k: len(v) for k, v in expected.items()}
    assert data["source"] == source


def test_stat_column_changed(out_dir, monkeypatch):
    spy = spies(monkeypatch)
    data = call("progeny", "TGFb", out_dir, stat_column="log2FoldChange")
    expected = direct_leading_edge(read_de_csv(DE_CSV), snapshot("progeny"), "TGFb", stat="log2FoldChange")
    assert le_files(data) == expected
    assert le_files(data)["positive_weights"] != read_lines(REF / "leading_edge" / "progeny_TGFb_positive_weights.txt")
    assert spy.calls["source_targets"][0]["y"] == "log2FoldChange"


@pytest.mark.parametrize(
    "args, fragment",
    [
        ({"source": "NOT_A_PATHWAY"}, "not found"),
        ({"de_results_path": str(DER / "deseq2_results_no_stat.csv")}, "lack column"),
        ({"de_results_path": "/nonexistent/de.csv"}, "does not exist"),
        ({"network_path": "/nonexistent/net.parquet"}, "does not exist"),
        ({"scores_path": "/nonexistent/scores.csv"}, "does not exist"),
        ({"scores_path": str(SIG_SCORES_CSV), "network_sources": ["TGFb", "NOT_A_TF"]}, "not in the scores table"),
        ({"scores_path": str(DER / "scores_two_rows.csv")}, "exactly one"),
    ],
)
def test_errors(args, fragment, out_dir):
    full = {"de_results_path": str(DE_CSV), "network": "progeny", "source": "TGFb",
            "network_path": str(NETS / "progeny.parquet"), "output_dir": str(out_dir), **args}
    assert fragment in call_tool_error(TOOL, full)


def test_repeated_calls_fresh_dirs(out_dir):
    a, b = call("hallmark", "TGF_BETA_SIGNALING", out_dir), call("hallmark", "TGF_BETA_SIGNALING", out_dir)
    pa, pb = check_common_return(a, out_dir, 4), check_common_return(b, out_dir, 4)
    assert pa[0].parent != pb[0].parent and not set(pa) & set(pb)
    assert le_files(a) == le_files(b)
