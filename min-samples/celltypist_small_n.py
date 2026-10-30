"""Run celltypist.annotate (with and without majority voting) and celltypist.train on small subsets of PBMC 3k."""
import sys, logging
import numpy as np, scanpy as sc, celltypist
logging.disable(logging.WARNING)
h5ad, model = sys.argv[1], sys.argv[2]
full = sc.read_h5ad(h5ad)
lab = [c for c in full.obs.columns if "louvain" in c or "label" in c or "cell" in c][0]
print("label column:", lab, "| cells:", full.n_obs)
rng = np.random.default_rng(0)
for n in (5, 12, 30, 51, 60, 200):
    sub = full[rng.choice(full.n_obs, n, replace=False)].copy()
    for mv in (False, True):
        try:
            p = celltypist.annotate(sub.copy(), model=model, majority_voting=mv)
            col = "majority_voting" if mv else "predicted_labels"
            out = f"OK ({p.predicted_labels[col].nunique()} labels)"
        except Exception as e:
            out = f"ERROR {type(e).__name__}: {str(e).splitlines()[0][:110]}"
        print(f"annotate n={n:<4} majority_voting={mv!s:<5} {out}")
# training: per-class cell counts
for per_class, k in ((1, 2), (3, 2), (10, 3)):
    cats = full.obs[lab].value_counts().index[:k]
    idx = np.concatenate([np.where(full.obs[lab] == c)[0][:per_class] for c in cats])
    sub = full[idx].copy()
    for fs in (False, True):
        try:
            m = celltypist.train(sub, labels=lab, feature_selection=fs, top_genes=300, random_state=0)
            out = f"OK ({len(m.cell_types)} types, {len(m.features)} features)"
        except Exception as e:
            out = f"ERROR {type(e).__name__}: {str(e).splitlines()[0][:110]}"
        print(f"train {k} types x {per_class} cells feature_selection={fs!s:<5} {out}")
one = full[np.where(full.obs[lab] == full.obs[lab].value_counts().index[0])[0][:50]].copy()
try:
    celltypist.train(one, labels=lab, random_state=0); print("train 1 type x 50 cells OK")
except Exception as e:
    print(f"train 1 type x 50 cells ERROR {type(e).__name__}: {str(e).splitlines()[0][:110]}")
