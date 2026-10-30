"""Build an airway (GSE52778) AnnData in the shape of decoupler's dc.ds.hsctgfb():
samples x genes raw counts (float), gene symbols as var_names, categorical obs columns.

Keeps the 8 untreated / dexamethasone samples (4 airway smooth-muscle cell lines);
albuterol arms are dropped. Sample annotations come from the GEO SOFT file, not row order.
"""
import gzip
import re
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

download, out = Path(sys.argv[1]), Path(sys.argv[2])

# GEO sample annotations
samples, cur = {}, None
with gzip.open(download / "GSE52778_family.soft.gz", "rt") as fh:
    for line in fh:
        if line.startswith("^SAMPLE"):
            cur = line.split("=", 1)[1].strip()
            samples[cur] = {}
        elif cur and line.startswith("!Sample_title"):
            samples[cur]["title"] = line.split("=", 1)[1].strip()
        elif cur and line.startswith("!Sample_characteristics_ch1"):
            key, _, val = line.split("=", 1)[1].strip().partition(": ")
            samples[cur][key] = val
meta = pd.DataFrame.from_dict(samples, orient="index")
meta = meta[meta["treatment"].isin(["Untreated", "Dexamethasone"])]
assert len(meta) == 8 and meta.groupby("treatment").size().eq(4).all()

# NCBI-generated raw counts (GeneID rows x GSM columns)
counts = pd.read_csv(download / "GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="\t", index_col=0)
counts = counts[meta.index]

# GeneID -> current official symbol
info = pd.read_csv(download / "Homo_sapiens.gene_info.gz", sep="\t", usecols=["GeneID", "Symbol"])
sym = info.set_index("GeneID")["Symbol"]
mapped = counts.index.isin(sym.index)
print(f"GeneIDs: {len(counts)}; mapped to symbols: {mapped.sum()}; dropped (no current symbol): {(~mapped).sum()}")
counts = counts[mapped]
counts.index = sym.loc[counts.index].values
dups = counts.index.duplicated(keep="first")
print(f"duplicate symbols dropped (first kept, as hsctgfb does): {dups.sum()}")
counts = counts[~dups]

adata = ad.AnnData(counts.T.astype(float))
adata.obs["condition"] = pd.Categorical(
    meta.loc[adata.obs_names, "treatment"].map({"Untreated": "untreated", "Dexamethasone": "dexamethasone"}),
    categories=["untreated", "dexamethasone"],
)
adata.obs["cell_line"] = pd.Categorical(meta.loc[adata.obs_names, "cell line"])
adata.obs["sample_id"] = pd.Categorical(meta.loc[adata.obs_names, "title"])
adata.uns["source"] = {
    "geo": "GSE52778",
    "counts": "NCBI GEO RNA-seq raw counts, GRCh38.p13 (GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz)",
    "symbols": "NCBI Homo_sapiens.gene_info (GeneID -> Symbol)",
    "subset": "Untreated and Dexamethasone samples only",
}
out.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(out)
print(adata)
print(adata.obs)
print("library sizes:", adata.X.sum(1).astype(int).tolist())
for g in ["FKBP5", "TSC22D3", "PER1", "KLF15", "DUSP1", "ZBTB16", "NR3C1"]:
    print(g, adata[:, g].X.ravel().astype(int).tolist() if g in adata.var_names else "MISSING")
