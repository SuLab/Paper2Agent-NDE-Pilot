"""Score the decoupler exemplar runs against their pre-declared expectations.

Usage: python -I evaluate_decoupler.py EXEMPLARS_DECOUPLER_DIR
Reads each <GSE>/results*.json written by mcp_call.py and writes <GSE>/expectation_results.json.
"""
import json
import sys
from pathlib import Path

import pandas as pd

root = Path(sys.argv[1])


def load(results_json):
    calls = {c["name"]: c["result"] for c in json.load(open(results_json))["calls"]}
    run = {"de": pd.read_csv(calls["de"]["artifacts"][0]["path"], index_col=0),
           "n_genes_kept": calls["filter"]["n_genes_after"]}
    for net in ["collectri", "progeny", "hallmark"]:
        arts = calls[f"ulm_{net}"]["artifacts"]
        run[net] = pd.DataFrame({"score": pd.read_csv(arts[0]["path"], index_col=0).iloc[0],
                                 "padj": pd.read_csv(arts[1]["path"], index_col=0).iloc[0]})
    pcs = calls["pcs"]["artifacts"]
    run["pc_assoc"] = pd.read_csv(pcs[0]["path"])
    run["pc_var"] = pd.read_csv(pcs[1]["path"])
    return run


def up(df, name):
    """score/log2FC > 0 and padj < 0.05; a name absent from df counts as not up."""
    if name not in df.index:
        return False, None
    row = df.loc[name]
    val = row["score"] if "score" in row else row["log2FoldChange"]
    return bool(val > 0 and row["padj"] < 0.05), {"value": round(float(val), 3), "padj": float(row["padj"])}


def check_up(df, names):
    detail = {n: up(df, n) for n in names}
    return all(v[0] for v in detail.values()), {n: v[1] for n, v in detail.items()}


def top(df, k=10):
    d = df.sort_values("score", ascending=False)
    return {"up": [f"{s} ({r.score:+.2f})" for s, r in d.head(k).iterrows()],
            "down": [f"{s} ({r.score:+.2f})" for s, r in d.tail(k)[::-1].iterrows()]}


def rank(df, name):
    return int(df["score"].rank(ascending=False)[name]) if name in df.index else None


def report(run):
    sig = run["pc_assoc"][run["pc_assoc"]["padj"] < 0.05]
    return {"genes_kept": run["n_genes_kept"],
            "de_genes_padj_0.05": int((run["de"]["padj"] < 0.05).sum()),
            "pc_variance_ratio_top4": run["pc_var"].iloc[:4].round(3).to_dict("records"),
            "pc_metadata_associations_padj_0.05": sig.round(4).to_dict("records"),
            "top_collectri": top(run["collectri"]), "n_tfs_scored": len(run["collectri"]),
            "top_progeny": top(run["progeny"], 5), "top_hallmark": top(run["hallmark"], 5)}


def lfc(run, g):
    return float(run["de"].loc[g, "log2FoldChange"]) if g in run["de"].index else 0.0


out = {}

# GSE164287: Mtb vs uninfected macrophages
r = load(root / "GSE164287/results.json")
res = {}
res["E1"] = check_up(r["collectri"], ["NFKB1", "RELA"])
res["E2"] = check_up(r["progeny"], ["NFkB", "TNFa"])
res["E3"] = check_up(r["hallmark"], ["TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE"])
res["E4"] = check_up(r["de"], ["IL1B", "IL6", "CXCL8", "PTGS2"])
rep = report(r)
rep["ranks_of_NFKB1_RELA"] = {g: rank(r["collectri"], g) for g in ["NFKB1", "RELA"]}
rep["type_I_IFN"] = {"collectri": check_up(r["collectri"], ["STAT1", "STAT2", "IRF9"])[1],
                     "hallmark": check_up(r["hallmark"], ["INTERFERON_ALPHA_RESPONSE"])[1]}
out["GSE164287"] = (res, rep)

# GSE153970: SARS-CoV-2 vs mock airway epithelium
r = load(root / "GSE153970/results.json")
res = {}
res["E1"] = check_up(r["collectri"], ["NFKB1", "RELA"])
res["E2"] = check_up(r["collectri"], ["ATF4"])
res["E3"] = check_up(r["de"], ["IL6", "TNF", "CXCL8"])
ifn = {"hallmark INTERFERON_ALPHA_RESPONSE": up(r["hallmark"], "INTERFERON_ALPHA_RESPONSE"),
       "collectri STAT2": up(r["collectri"], "STAT2"), "collectri IRF9": up(r["collectri"], "IRF9"),
       "IFNB1": up(r["de"], "IFNB1"), "IFNL1": up(r["de"], "IFNL1")}
res["E4"] = (not any(v[0] for v in ifn.values()),
             {k: (v[1] if v[1] is not None else "removed by filtering") for k, v in ifn.items()})
rep = report(r)
rep["ranks_of_NFKB1_RELA_ATF4"] = {g: rank(r["collectri"], g) for g in ["NFKB1", "RELA", "ATF4"]}
rep["superseded_NDE_hypothesis_interferon_STAT_IRF"] = {
    "collectri": check_up(r["collectri"], ["STAT1", "STAT2", "IRF1", "IRF3", "IRF7", "IRF9"])[1],
    "progeny JAK-STAT": up(r["progeny"], "JAK-STAT")[1],
    "hallmark": check_up(r["hallmark"], ["INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE"])[1]}
out["GSE153970"] = (res, rep)

# GSE198542: ZIKV vs mock, WT and IFNAR2-deficient (MUT) iPSC-macrophages
wt, mut = load(root / "GSE198542/results_WT.json"), load(root / "GSE198542/results_MUT.json")
res = {}
res["E1"] = check_up(wt["hallmark"], ["INTERFERON_ALPHA_RESPONSE"])
res["E2"] = check_up(wt["collectri"], ["STAT1", "STAT2", "IRF9"])
res["E3"] = check_up(wt["de"], ["MX1", "OAS1"])
isgf3 = ["STAT1", "STAT2", "IRF9"]
h_wt, h_mut = wt["hallmark"].loc["INTERFERON_ALPHA_RESPONSE", "score"], mut["hallmark"].loc["INTERFERON_ALPHA_RESPONSE", "score"]
c_wt, c_mut = wt["collectri"].loc[isgf3, "score"].mean(), mut["collectri"].loc[isgf3, "score"].mean()
g = {x: (lfc(wt, x), lfc(mut, x)) for x in ["MX1", "OAS1"]}
e4 = {"a_hallmark_ifna": (bool(h_mut < 0.5 * h_wt), {"WT": round(h_wt, 3), "MUT": round(h_mut, 3)}),
      "b_isgf3_mean": (bool(c_mut < 0.5 * c_wt), {"WT": round(c_wt, 3), "MUT": round(c_mut, 3)}),
      "c_MX1_OAS1_lfc": (all(m < 0.5 * w for w, m in g.values()),
                          {x: {"WT": round(w, 3), "MUT": round(m, 3)} for x, (w, m) in g.items()})}
res["E4"] = (all(v[0] for v in e4.values()), e4)
rep = {"WT": report(wt), "MUT": report(mut)}
for k, run in [("WT", wt), ("MUT", mut)]:
    rep[k]["IRF3_axis"] = {"IFNB1": up(run["de"], "IFNB1")[1], "ISG15": up(run["de"], "ISG15")[1],
                           "IFIT1": up(run["de"], "IFIT1")[1], "CXCL10": up(run["de"], "CXCL10")[1],
                           "collectri IRF3": up(run["collectri"], "IRF3")[1]}
    rep[k]["NFkB"] = {"collectri RELA": up(run["collectri"], "RELA")[1], "progeny NFkB": up(run["progeny"], "NFkB")[1]}
out["GSE198542"] = (res, rep)

for gse, (res, rep) in out.items():
    doc = {"expectations": {k: {"pass": bool(v[0]), "evidence": v[1]} for k, v in res.items()},
           "verified_sensible": all(v[0] for v in res.values()), "report_only": rep}
    json.dump(doc, open(root / gse / "expectation_results.json", "w"), indent=1, default=str)
    print(gse, {k: v[0] for k, v in res.items()}, "-> sensible" if doc["verified_sensible"] else "-> NOT all met")
