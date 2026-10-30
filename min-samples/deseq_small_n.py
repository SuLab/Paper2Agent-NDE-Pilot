"""Run PyDESeq2 (as the decoupler MCP calls it) and ulm on simulated counts at small sample sizes."""
import warnings
import numpy as np, pandas as pd, anndata as ad
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
import decoupler as dc

rng = np.random.default_rng(0)
G, TRUE = 2000, 200                      # genes; first 200 are 4-fold up in "trt"
base = rng.lognormal(mean=5, sigma=1.5, size=G)
disp = 0.05 + 1 / base                   # realistic-ish mean-dispersion trend
lfc = np.zeros(G); lfc[:TRUE] = 2.0
donor_eff = lambda: rng.normal(0, 0.3, size=G)

def sim(n_ctrl, n_trt, paired=False):
    rows, obs = [], []
    for cond, n in (("ctrl", n_ctrl), ("trt", n_trt)):
        for i in range(n):
            obs.append({"condition": cond, "donor": f"d{i}"})
    effs = {f"d{i}": donor_eff() for i in range(max(n_ctrl, n_trt))}
    for o in obs:
        mu = base * np.exp(np.log(2) * lfc * (o["condition"] == "trt"))
        if paired:
            mu = mu * np.exp(effs[o["donor"]])
        p = 1 / (1 + mu * disp)
        rows.append(rng.negative_binomial(1 / disp, p))
    a = ad.AnnData(np.array(rows, dtype=np.int64), obs=pd.DataFrame(obs, index=[f"s{i}" for i in range(len(obs))]))
    a.var_names = [f"g{i}" for i in range(G)]
    return a

# ulm network: one source whose targets are the true DE genes, one decoy
net = pd.DataFrame({"source": ["UP"] * 50 + ["DECOY"] * 50,
                    "target": [f"g{i}" for i in range(50)] + [f"g{i}" for i in range(1000, 1050)],
                    "weight": 1.0})

cases = [("1 vs 1", 1, 1, False), ("1 vs 2", 1, 2, False), ("2 vs 2", 2, 2, False),
         ("3 vs 3", 3, 3, False), ("4 vs 4", 4, 4, False),
         ("2 donors paired (~donor+condition)", 2, 2, True), ("3 donors paired", 3, 3, True),
         ("4 donors paired", 4, 4, True)]
for name, nc, nt, paired in cases:
    a = sim(nc, nt, paired)
    design = "~donor + condition" if paired else "~condition"
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            inf = DefaultInference(n_cpus=1)
            dds = DeseqDataSet(adata=a.copy(), design=design, refit_cooks=True, inference=inf, quiet=True)
            dds.deseq2()
            st = DeseqStats(dds, contrast=["condition", "trt", "ctrl"], inference=inf, quiet=True)
            st.summary()
            r = st.results_df
            sig = r["padj"] < 0.05
            tp = int(sig.iloc[:TRUE].sum()); fp = int(sig.iloc[TRUE:].sum())
            stat = r[["stat"]].T.fillna(0)
            acts, padj = dc.mt.ulm(data=stat, net=net)
            msg = (f"OK  true DE found {tp}/{TRUE}, false positives {fp}; "
                   f"ulm UP score {acts['UP'].iloc[0]:.1f} (padj {padj['UP'].iloc[0]:.1e}), DECOY {acts['DECOY'].iloc[0]:.1f}")
        except Exception as e:
            msg = f"ERROR {type(e).__name__}: {str(e)[:160]}"
    ws = sorted({str(x.message)[:70] for x in w if "residual degrees" in str(x.message) or "full rank" in str(x.message)})
    print(f"{name:<38} n={a.n_obs:<2} resid_df={a.n_obs - (1 + 1 + ((nc if paired else 1) - 1))}  {msg}  {('WARN: ' + '; '.join(ws)) if ws else ''}")
