# NDE searches and spot checks

All queries go to the NDE Discovery API, `https://api.data.niaid.nih.gov/v1/query`. Counts and hit lists are from 2026-10-28.

## Queries

**decoupler: human bulk RNA-seq series with at least 3 samples (48,684 hits)**

```
@type:"Dataset" AND measurementTechnique.name:"high-throughput expression assay"
AND species.name:"homo sapiens" AND sample.numberOfItems.value:>=3
NOT ("single-cell" OR "single cell" OR "single-nucleus" OR "single nucleus" OR scRNA-seq OR snRNA-seq OR scRNAseq OR snRNAseq)
```

**CellTypist: single-cell RNA-seq, any species (117,912 hits; 43,409 tagged human)**

```
@type:"Dataset" AND ("single-cell" OR "single cell" OR "single-nucleus" OR "single nucleus" OR scRNA-seq OR snRNA-seq OR scRNAseq OR snRNAseq)
AND (RNA OR transcriptome OR transcriptomic OR transcriptomics OR "gene expression" OR scRNA-seq OR snRNA-seq OR scRNAseq OR snRNAseq)
```

`docs/nde-pilot.md` explains each filter, and why the queries leave out species, tissue, disease and file-format filters.

## Scripts

| Script | What it does |
|---|---|
| `scripts/nde_query.py` | Queries the API, with facets and with scrolling past the 1,000-hit page limit. Uses the standard library only. |
| `scripts/geo_check.py OUT_DIR IDS_JSON` | For each GEO series: reads the SOFT sample sheet (organism, characteristics), lists the supplementary files and the contents of `RAW.tar`, and checks whether NCBI has generated raw counts. |
| `scripts/ncbi_to_gse.py RECORDS_JSON OUT_JSON` | Maps SRA and BioProject records to their GEO series with NCBI E-utilities. |

## Results

| File | Contents |
|---|---|
| `results/bulk_human_ge3_all.json.gz` | All 48,684 bulk hits: identifier, name, sample count. |
| `results/sc_human_all.json.gz` | The 43,409 human-tagged single-cell hits: identifier, name, catalog, URL. |
| `results/bulk_sample.json`, `results/sc_sample.json` | The 20 records drawn at random from each list for the spot check. The sampling code was run inline and not saved. |
| `results/sc_gse_map.json` | The GEO series found for each sampled SRA or BioProject record. |
| `results/bulk_spotcheck.csv`, `results/sc_spotcheck.csv` | Verdict and reason for each spot-checked record. |
| `results/tool-records/` | The NDE records `biotools_decoupler` and `biotools_celltypist`, and the bio.tools source records behind them. |

## How far the spot checks went

- **Formats:** inferred from file names.
- **Not opened:** the authors' count files and the per-sample archives.
- **Species:** the organism was taken from the GEO sample sheet, not from NDE's species tag. That tag was wrong for 3 of the 40 records; in each case it had been inferred from the linked paper.
