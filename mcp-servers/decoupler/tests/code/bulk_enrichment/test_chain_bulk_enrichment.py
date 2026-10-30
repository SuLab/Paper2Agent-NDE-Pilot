"""Composition, concurrency and stdio-transport checks for the bulk_enrichment tools.

The chain passes each tool's own output file to the next tool explicitly (no reliance on test
order) and compares every stage with the notebook reference.
"""

import asyncio
import json
import queue
import subprocess
import threading
import time

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd

from conftest import (
    CONTRAST, DE_CSV, DE_TIMEOUT, DER, FILTERED_H5AD, NETS, NETWORKS, PROJECT_PYTHON, RAW_H5AD, REF, SERVER, SRC,
    Spy, artifact, assert_frame_exact, assert_png, call_tool, check_common_return, read_de_csv, read_lines,
)
from fastmcp import Client

CELL36_SOURCES = ["SRF", "MYOCD", "IRF1", "SKIL"]


def test_chain_filter_rank_de_ulm_inspect(out_dir, monkeypatch):
    spy = Spy()
    spy.wrap(monkeypatch, dc.pl, "network")
    f = call_tool("decoupler_filter_genes_by_expression",
                  {"data_path": str(RAW_H5AD), "group_key": "condition", "output_dir": str(out_dir)})
    filtered = artifact(f, "filtered_counts.h5ad")
    assert list(ad.read_h5ad(filtered).var_names) == read_lines(REF / "filter_by_expr" / "kept_genes.txt")

    r = call_tool("decoupler_rank_pcs_by_metadata", {"data_path": str(filtered), "output_dir": str(out_dir)})
    assoc = pd.read_csv(artifact(r, "pc_metadata_associations.csv")).sort_values(["obsm", "obs"])
    ref_assoc = pd.read_parquet(REF / "pca_rankby_obsm" / "rank_obsm.parquet").sort_values(["obsm", "obs"])
    np.testing.assert_allclose(assoc["stat"].abs(), ref_assoc["stat"].abs(), rtol=1e-5)
    np.testing.assert_allclose(assoc["padj"], ref_assoc["padj"], rtol=1e-5)

    d = call_tool("pydeseq2_test_differential_expression",
                  {"data_path": str(filtered), "condition_key": "condition", "group": "treatment",
                   "reference": "control", "output_dir": str(out_dir)}, timeout=DE_TIMEOUT)
    de_csv = artifact(d, "deseq2_results.csv")
    assert_frame_exact(read_de_csv(de_csv), pd.read_parquet(REF / "deseq2" / "results_df.parquet"))

    sig_paths = {}
    for net in NETWORKS:
        u = call_tool("decoupler_score_activities_ulm",
                      {"de_results_path": str(de_csv), "network": net, "network_path": str(NETS / f"{net}.parquet"),
                       "contrast_name": CONTRAST, "output_dir": str(out_dir)})
        scores = pd.read_csv(artifact(u, f"{net}_ulm_scores.csv"), index_col=0)
        padj = pd.read_csv(artifact(u, f"{net}_ulm_padj.csv"), index_col=0)
        assert_frame_exact(scores, pd.read_parquet(REF / "ulm" / f"{net}_scores_full.parquet"))
        assert_frame_exact(padj, pd.read_parquet(REF / "ulm" / f"{net}_padj_full.parquet"))
        sig_paths[net] = artifact(u, f"{net}_ulm_scores_significant.csv")
        assert u["n_significant"] == {"collectri": 143, "progeny": 12, "hallmark": 30}[net]

    p = call_tool("decoupler_inspect_source_targets",
                  {"de_results_path": str(de_csv), "network": "progeny", "source": "TGFb",
                   "network_path": str(NETS / "progeny.parquet"), "output_dir": str(out_dir)})
    assert read_lines(p["leading_edges"]["positive_weights"]["path"]) == read_lines(
        REF / "leading_edge" / "progeny_TGFb_positive_weights.txt")
    assert read_lines(p["leading_edges"]["negative_weights"]["path"]) == read_lines(
        REF / "leading_edge" / "progeny_TGFb_negative_weights.txt")
    h = call_tool("decoupler_inspect_source_targets",
                  {"de_results_path": str(de_csv), "network": "hallmark", "source": "TGF_BETA_SIGNALING",
                   "network_path": str(NETS / "hallmark.parquet"), "output_dir": str(out_dir)})
    assert read_lines(h["leading_edges"]["all_targets"]["path"]) == read_lines(
        REF / "leading_edge" / "hallmark_TGF_BETA_SIGNALING.txt")
    c = call_tool("decoupler_inspect_source_targets",
                  {"de_results_path": str(de_csv), "network": "collectri", "source": "SRF",
                   "scores_path": str(sig_paths["collectri"]), "network_sources": CELL36_SOURCES,
                   "network_path": str(NETS / "collectri.parquet"), "output_dir": str(out_dir)})
    (nw,) = spy.calls["network"]
    assert_frame_exact(nw["score"], pd.read_parquet(REF / "ulm" / "collectri_scores_padj_lt_0.05.parquet"))
    assert_frame_exact(nw["data"], pd.read_parquet(REF / "ulm" / "data_stat_matrix.parquet"))
    assert_png(artifact(c, "SRF_network.png"))


def test_concurrent_calls_share_one_server(out_dir, monkeypatch):
    """FastMCP runs sync tools in worker threads; plotting must stay serialised and results correct."""
    spy_log = []
    lock = threading.Lock()
    for name in ("obsm", "barplot", "dotplot", "volcano", "leading_edge", "network", "source_targets"):
        orig = getattr(dc.pl, name)

        def make(orig=orig, name=name):
            def wrapped(*a, **k):
                t0 = time.perf_counter()
                try:
                    return orig(*a, **k)
                finally:
                    with lock:
                        spy_log.append((name, threading.get_ident(), t0, time.perf_counter()))
            return wrapped

        monkeypatch.setattr(dc.pl, name, make())

    calls = [
        ("decoupler_rank_pcs_by_metadata", {"data_path": str(FILTERED_H5AD)}),
        ("decoupler_score_activities_ulm", {"de_results_path": str(DE_CSV), "network": "progeny",
                                            "network_path": str(NETS / "progeny.parquet"), "contrast_name": CONTRAST}),
        ("decoupler_score_activities_ulm", {"de_results_path": str(DE_CSV), "network": "hallmark",
                                            "network_path": str(NETS / "hallmark.parquet"), "contrast_name": CONTRAST}),
        ("decoupler_inspect_source_targets", {"de_results_path": str(DE_CSV), "network": "progeny", "source": "TGFb",
                                              "network_path": str(NETS / "progeny.parquet")}),
        ("decoupler_inspect_source_targets", {"de_results_path": str(DE_CSV), "network": "collectri", "source": "SRF",
                                              "scores_path": str(DE_CSV.parent / "collectri_scores_significant_reference.csv"),
                                              "network_sources": CELL36_SOURCES,
                                              "network_path": str(NETS / "collectri.parquet")}),
    ]

    async def go():
        async with Client(SERVER) as client:
            return await asyncio.gather(*[client.call_tool(n, {**a, "output_dir": str(out_dir)}, timeout=600,
                                                           raise_on_error=False) for n, a in calls])

    results = asyncio.run(go())
    assert not any(r.is_error for r in results), [r.content for r in results if r.is_error]
    data = [r.data for r in results]
    for d in data:
        for p in check_common_return(d, out_dir, 1):
            if p.suffix == ".png":
                assert_png(p)
    assert_frame_exact(pd.read_csv(artifact(data[1], "progeny_ulm_scores.csv"), index_col=0),
                       pd.read_parquet(REF / "ulm" / "progeny_scores_full.parquet"))
    assert_frame_exact(pd.read_csv(artifact(data[2], "hallmark_ulm_scores.csv"), index_col=0),
                       pd.read_parquet(REF / "ulm" / "hallmark_scores_full.parquet"))
    assert read_lines(data[3]["leading_edges"]["positive_weights"]["path"]) == read_lines(
        REF / "leading_edge" / "progeny_TGFb_positive_weights.txt")
    # The lock serialises upstream plotting: no two plot calls overlap in time.
    spans = sorted((s, e) for _, _, s, e in spy_log)
    assert all(spans[i][1] <= spans[i + 1][0] for i in range(len(spans) - 1))
    threads = {t for _, t, _, _ in spy_log}
    (out_dir / "concurrency.json").write_text(json.dumps({"plot_calls": len(spy_log), "worker_threads": len(threads)}))


class StdioSession:
    """Minimal newline-delimited JSON-RPC client that keeps every raw stdout line."""

    def __init__(self, argv, cwd):
        self.proc = subprocess.Popen(argv, cwd=cwd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, text=True, bufsize=1)
        self.lines, self.stderr, self.q = [], [], queue.Queue()
        threading.Thread(target=self._read_out, daemon=True).start()
        threading.Thread(target=lambda: self.stderr.extend(self.proc.stderr), daemon=True).start()

    def _read_out(self):
        for line in self.proc.stdout:
            self.lines.append(line)
            self.q.put(line)

    def send(self, msg):
        self.proc.stdin.write(json.dumps(msg) + "\n")
        self.proc.stdin.flush()

    def request(self, rid, method, params, timeout):
        self.send({"jsonrpc": "2.0", "id": rid, "method": method, "params": params})
        deadline = time.time() + timeout
        while time.time() < deadline:
            try:
                line = self.q.get(timeout=max(0.1, deadline - time.time()))
            except queue.Empty:
                break
            msg = json.loads(line)
            if msg.get("id") == rid:
                return msg
        raise TimeoutError(f"no response to {method} (id {rid}); stderr tail: {''.join(self.stderr[-20:])}")

    def close(self):
        self.proc.stdin.close()
        try:
            self.proc.wait(timeout=30)
        except subprocess.TimeoutExpired:
            self.proc.kill()


def test_stdio_transport_stdout_is_pure_jsonrpc(out_dir):
    launcher = out_dir / "serve_bulk_enrichment.py"
    launcher.write_text(
        "import sys\n"
        f"sys.path.insert(0, {str(SRC)!r})\n"
        "from tools.bulk_enrichment import bulk_enrichment_mcp\n"
        "if __name__ == '__main__':\n"
        "    bulk_enrichment_mcp.run()\n"
    )
    s = StdioSession([str(PROJECT_PYTHON), "-I", str(launcher)], cwd=out_dir)
    try:
        init = s.request(1, "initialize", {"protocolVersion": "2025-06-18", "capabilities": {},
                                           "clientInfo": {"name": "bulk-enrichment-verifier", "version": "2"}}, 120)
        assert "result" in init
        s.send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        listed = s.request(2, "tools/list", {}, 60)
        assert {t["name"] for t in listed["result"]["tools"]} == {
            "decoupler_filter_genes_by_expression", "decoupler_rank_pcs_by_metadata",
            "pydeseq2_test_differential_expression", "decoupler_score_activities_ulm",
            "decoupler_inspect_source_targets"}
        calls = [
            ("decoupler_filter_genes_by_expression", {"data_path": str(RAW_H5AD), "group_key": "condition"}),
            ("decoupler_rank_pcs_by_metadata", {"data_path": str(FILTERED_H5AD)}),
            ("pydeseq2_test_differential_expression", {"data_path": str(DER / "filtered_small_2000genes.h5ad"),
                                                       "condition_key": "condition", "group": "treatment",
                                                       "reference": "control"}),
            ("decoupler_score_activities_ulm", {"de_results_path": str(DE_CSV), "network": "hallmark",
                                                "network_path": str(NETS / "hallmark.parquet")}),
            ("decoupler_inspect_source_targets", {"de_results_path": str(DE_CSV), "network": "hallmark",
                                                  "source": "TGF_BETA_SIGNALING",
                                                  "network_path": str(NETS / "hallmark.parquet")}),
        ]
        for i, (name, args) in enumerate(calls, start=3):
            resp = s.request(i, "tools/call", {"name": name, "arguments": {**args, "output_dir": str(out_dir)}},
                             DE_TIMEOUT)
            assert not resp["result"].get("isError"), resp
            assert resp["result"]["structuredContent"]["artifacts"]
    finally:
        s.close()
    # Every byte on stdout must be a JSON-RPC 2.0 message (no prints from PyDESeq2, scanpy or decoupler).
    assert s.lines
    for line in s.lines:
        msg = json.loads(line)
        assert isinstance(msg, dict) and msg.get("jsonrpc") == "2.0", line
    (out_dir / "stdio_transcript_summary.json").write_text(json.dumps(
        {"stdout_lines": len(s.lines), "all_jsonrpc": True, "stderr_lines": len(s.stderr)}, indent=2))
