"""Build a decoupler-ready AnnData from a GEO series with NCBI-generated raw counts, in the shape of
decoupler's dc.ds.hsctgfb(): samples x genes raw counts (float), gene symbols as var_names, categorical obs.

Usage: python -I prep_geo_bulk.py CONFIG_JSON DOWNLOAD_DIR GENE_INFO_GZ OUT_H5AD

CONFIG_JSON:
  {"gse": "GSE...",
   "condition": {"field": "<characteristic or 'title'>", "levels": {"<GEO value>": "<label>", ...}},  # first label = reference
   "covariates": {"<obs name>": "<characteristic or 'title' or 'extraction date'>", ...},
   "subset": {"<characteristic>": "<value>"}}                                                          # optional
Sample annotations come from the GEO SOFT file, not column order.
"""
import gzip
import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

cfg_path, download, gene_info, out = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
cfg = json.load(open(cfg_path))
gse = cfg["gse"]

samples, cur = {}, None
with gzip.open(download / f"{gse}_family.soft.gz", "rt") as fh:
    for line in fh:
        if line.startswith("^SAMPLE"):
            cur = line.split("=", 1)[1].strip()
            samples[cur] = {}
        elif cur and line.startswith("!Sample_title"):
            samples[cur]["title"] = line.split("=", 1)[1].strip()
        elif cur and line.startswith("!Sample_organism_ch1"):
            samples[cur]["organism"] = line.split("=", 1)[1].strip()
        elif cur and line.startswith("!Sample_characteristics_ch1"):
            key, _, val = line.split("=", 1)[1].strip().partition(": ")
            samples[cur][key] = val
meta = pd.DataFrame.from_dict(samples, orient="index")
assert meta["organism"].eq("Homo sapiens").all(), meta["organism"].value_counts()
for key, val in cfg.get("subset", {}).items():
    meta = meta[meta[key] == val]
cond = cfg["condition"]
meta = meta[meta[cond["field"]].isin(cond["levels"])]
print(f"{gse}: {len(meta)} samples kept;", meta[cond["field"]].value_counts().to_dict())

counts = pd.read_csv(download / f"{gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="\t", index_col=0)
counts = counts[meta.index]

info = pd.read_csv(gene_info, sep="\t", usecols=["GeneID", "Symbol"])
sym = info.set_index("GeneID")["Symbol"]
mapped = counts.index.isin(sym.index)
print(f"GeneIDs: {len(counts)}; mapped to symbols: {mapped.sum()}; dropped (no current symbol): {(~mapped).sum()}")
counts = counts[mapped]
counts.index = sym.loc[counts.index].values
dups = counts.index.duplicated(keep="first")
print(f"duplicate symbols dropped (first kept, as hsctgfb does): {dups.sum()}")
counts = counts[~dups]

adata = ad.AnnData(counts.T.astype(float))
labels = list(dict.fromkeys(cond["levels"].values()))
adata.obs["condition"] = pd.Categorical(meta.loc[adata.obs_names, cond["field"]].map(cond["levels"]), categories=labels)
for name, field in cfg.get("covariates", {}).items():
    adata.obs[name] = pd.Categorical(meta.loc[adata.obs_names, field])
adata.obs["sample_title"] = pd.Categorical(meta.loc[adata.obs_names, "title"])
adata.uns["source"] = {
    "geo": gse,
    "counts": f"NCBI GEO RNA-seq raw counts, GRCh38.p13 ({gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz)",
    "symbols": f"NCBI Homo_sapiens.gene_info (GeneID -> Symbol), {gene_info.name}",
    "config": json.dumps(cfg),
}
out.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(out)
print(adata)
print(adata.obs)
print("library sizes:", adata.X.sum(1).astype(int).tolist())
