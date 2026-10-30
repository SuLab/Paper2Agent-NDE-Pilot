"""Build a CellTypist-ready PBMC 3k AnnData: log1p expression normalised to 10,000 counts per cell,
all 13,714 retained genes, with the scanpy-tutorial cluster labels kept in obs['louvain'].

The processed file's .X is scaled HVG data and its .raw.X is log1p(raw counts) (not depth-normalised),
so counts are recovered as expm1(raw.X), checked to be integers matching obs['n_counts'], then
normalised as CellTypist's README requires.
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

src, out = Path(sys.argv[1]), Path(sys.argv[2])
a = ad.read_h5ad(src)
raw = a.raw.to_adata()
X = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
X = X.astype(np.float64)
X.data = np.expm1(X.data)
dev = np.abs(X.data - np.round(X.data)).max()
print(f"max deviation of expm1(raw.X) from integers: {dev:.2e}")
assert dev < 1e-3
X.data = np.round(X.data)
totals = np.asarray(X.sum(1)).ravel()
print("recovered totals == obs n_counts:", np.allclose(totals, a.obs["n_counts"].to_numpy()))

counts = ad.AnnData(X=X.astype(np.float32), obs=a.obs[["n_genes", "percent_mito", "n_counts", "louvain"]].copy(),
                    var=pd.DataFrame(index=raw.var_names))
counts.layers["counts"] = counts.X.copy()
# log1p normalised to 10,000 counts per cell (CellTypist input contract)
lib = np.asarray(counts.X.sum(1)).ravel()
norm = sp.diags(1e4 / lib) @ counts.X
norm = sp.csr_matrix(norm)
norm.data = np.log1p(norm.data)
counts.X = norm.astype(np.float32)
counts.obs = counts.obs.rename(columns={"louvain": "louvain"})
counts.obsm["X_umap"] = a.obsm["X_umap"]
counts.uns["source"] = {
    "dataset": "10x Genomics PBMC 3k (Satija/Seurat & scanpy tutorial), processed h5ad from chanzuckerberg/cellxgene example-dataset/pbmc3k.h5ad",
    "labels": "obs['louvain']: cluster annotations from the scanpy PBMC3k tutorial (marker-gene based)",
    "X": "log1p(counts normalised to 1e4 per cell); layers['counts'] = recovered integer counts",
}
out.parent.mkdir(parents=True, exist_ok=True)
counts.write_h5ad(out)
check = np.expm1(counts.X[:5].toarray()).sum(1)
print(counts)
print("expm1 row sums (first 5):", np.round(check, 1).tolist())
print(counts.obs["louvain"].value_counts().to_dict())
