# Applying Paper2Agent tools to NIAID Data Ecosystem datasets

*Follow-on to [Exploring Paper2Agent: from papers to working analysis tools](https://docs.google.com/document/d/13MywPte5tHFky83whzLX7COnULGPHVTFhtb36zRfR7I/edit?tab=t.0#heading=h.ssj1d4wdjbj2). October 28, 2026.*

> This is a Markdown snapshot of the working Google Doc as of 2026-10-28. The exemplar runs in [`runs/exemplars`](../runs/exemplars/STATUS.md) postdate it, so its statements that no exemplar has been run are out of date.

## Summary

The companion doc shows that Paper2Agent's Paper2MCP workflow can turn a published method's code into tested tools an AI agent can call. It built two:

- **decoupleR**, which finds the transcription factors and pathways that change between two groups of bulk RNA-seq samples;
- **CellTypist**, which assigns a cell type to every cell in a single-cell RNA-seq dataset.

Each tool comes with an explicit input contract: the kind of data it needs. Matching those contracts against the NIAID Data Ecosystem (NDE) lets NDE add a layer of annotation to every dataset that fits. That annotation would say three things:

- the dataset can be analysed with this tool;
- the tool has been run on it;
- what the run found, such as regulator activity for the dataset's comparison or a cell-type label for every cell.

The two vignettes below follow each tool from its record in NDE's Computational Tools to the NDE datasets it could analyse. The approach scales with the number of compatible datasets, so the searches use only the filters each tool truly requires.

- **decoupleR:** 48,684 human bulk RNA-seq series in NDE have at least 3 samples. Another 44,159 mouse series could be added: decoupleR can translate its networks to mouse, and the server already accepts a saved network file.
  - In a random sample of 20, 14 could be analysed now, 1 was unclear and 5 could not.
  - If the sample is representative, that suggests roughly 23,000–42,000 runnable human series.
- **CellTypist:** 117,912 NDE records describe single-cell RNA-seq.
  - In a random sample of 20 human-tagged records, 16 had a downloadable expression matrix in a standard format.
  - None was a ready-to-use AnnData file, so a routine conversion step is needed.
- **Exemplars:** we picked a handful of datasets by hand for each tool. All of them fall inside these searches, and they confirm that the searches return datasets the tools can use. None of them has been run yet.

NDE's tool records and dataset records are barely linked today. MCP servers with input contracts could link them, and file-level descriptions in the style of Croissant would make the matching more accurate.

## The idea: tools with input contracts, matched against NDE metadata

A dataset's annotation can be built up in three tiers:

1. **Predicted compatible.** The dataset's metadata satisfies the tool's input contract. This can be computed across the whole catalog.
2. **Verified executable.** An agent downloads the dataset, adapts it to the contract, and runs the tool to completion. The adaptation steps are saved so they can be reused.
3. **Verified sensible.** The output meets expectations written down before the run, such as a known regulator being activated by a treatment, or cell labels agreeing with the authors' annotation.

Each vignette goes through the same phases:

1. Find the tool in NDE's Computational Tools.
2. Wrap it as an MCP server with Paper2MCP, which fixes its input contract.
3. Search NDE's datasets for that contract.
4. Spot-check a random sample of the hits at their source repositories.
5. Confirm with hand-picked exemplars.

Each vignette starts where an NDE user would, at the tool's NDE record. Our own path started from the paper and its code.

All searches use the NDE Discovery API (api.data.niaid.nih.gov). Dataset searches are restricted to dataset records; NDE also holds per-sample records, which make up 8.9 million of its 14.6 million records. Counts are as of October 28, 2026. Disease, pathogen, tissue and treatment are not requirements of either tool, so the searches don't filter on them.

## Vignette 1: decoupleR on bulk RNA-seq

### Find a tool in NDE's Computational Tools

NDE holds 34,275 Computational Tool records, all harvested from bio.tools. Searching them for "transcription factor" and "activity" returns 75 tools. decoupleR is one of them (63rd by relevance), with the record [biotools\_decoupler](https://data.niaid.nih.gov/resources?id=biotools_decoupler).

- **What the record gives:** the paper ([Badia-i-Mompel et al. 2022](https://doi.org/10.1093/bioadv/vbac016)), the authors, the license and a link to the code.
- **Which implementation:** the record describes the R/Bioconductor package (version 2.17.0). The server in the companion doc wraps the Python implementation from the same lab and paper (decoupler 2.3.1), which has no NDE record. Its differential-expression step uses PyDESeq2, which has no record either.
- **What it doesn't say:** what data the tool takes. Its listed functions are "Visualisation" and "Data retrieval", and its input and output fields are empty, so the record alone can't say which datasets the tool fits.

### Wrap it as an MCP server

The companion doc built five tools from decoupler's bulk RNA-seq tutorial: filter genes, rank principal components against sample metadata, test differential expression with PyDESeq2, score regulator and pathway activity, and inspect a regulator's targets. Run on airway smooth-muscle cells treated with dexamethasone (GSE52778), it ranked the glucocorticoid receptor 2nd of 657 transcription factors.

Its input contract, and where each requirement can be checked:

| Requirement | What the server needs | Where it is checked |
|---|---|---|
| Assay | Bulk RNA-seq | NDE search |
| Organism | Human (its regulator networks are human) | NDE search, then the GEO sample sheet, because NDE species tags can be wrong |
| Sample size | At least 3 samples. 1 vs 1 fails; 3 vs 3 or more is better. | NDE search (sample count) |
| Design | Two groups to compare, such as treated vs control | GEO sample sheet |
| Data files | Raw integer counts per gene | File listing, then file contents |

### Search NDE for compatible datasets

```
@type:"Dataset" AND measurementTechnique.name:"high-throughput expression assay"
AND species.name:"homo sapiens" AND sample.numberOfItems.value:>=3
NOT ("single-cell" OR "single cell" OR "single-nucleus" OR "single nucleus" OR scRNA-seq OR snRNA-seq OR scRNAseq OR snRNAseq)
```

Why each filter is there:

- **"High-throughput expression assay"** is NDE's term for GEO's "expression profiling by high-throughput sequencing", so count data exist. GEO microarray series are indexed separately, as "expression assay"; they would need a different differential-expression step.
- **Human,** because the server's regulator networks are human.
- **At least 3 samples,** the minimum for a two-group comparison (1 vs 2).
- **Not single-cell,** because single-cell series would first need to be pooled into pseudo-bulk samples. decoupleR has a separate tutorial for that, which this server doesn't include.

Whether a series actually contains two groups to compare is not recorded in NDE metadata, so that is read from the GEO sample sheet instead.

| Step | Datasets |
|---|---|
| GEO high-throughput expression assay, all species | 146,115 |
| Human | 68,577 |
| Excluding single-cell series | 53,765 |
| At least 3 samples | **48,684** |
| (of which at least 6 samples, enough for 3 vs 3) | 44,427 |
| Mouse equivalent of the final step (needs translated networks) | 44,159 |
| For comparison: our earlier infection-only search | 326 |

The earlier search required a pathogen tag and words like "mock" and "infected". It found useful datasets, but it was about 150 times narrower than the tool requires. It also excluded GSE52778, the dexamethasone dataset the server was actually run on.

### Spot-check a random sample

For 20 random series, we read the GEO sample sheet and the list of supplementary files, and checked whether NCBI has generated raw counts.

- **14 usable now:** two groups with at least 2 samples each, plus raw counts.
  - 9 have NCBI-generated counts.
  - 5 have authors' count files, identified by file name. These were not opened.
- **1 unclear:** the design fits, but the authors' expression file may not contain raw counts.
- **5 not usable:**
  - 2 are mouse-only series that NDE tagged as human based on the linked paper. NDE marks that tag as not curated.
  - 1 has no comparison.
  - 1 has no group information in its sample sheet.
  - 1 has no RNA-seq count matrix among its GEO files.

The usable designs were varied: drug vs vehicle, gene knockdown vs control, stem vs non-stem cancer cells, disease stages, environmental exposures. None was an infection study. The 95% interval for 14 of 20 is about 48–85%.

### Confirm with exemplars

These five were picked by hand from the earlier infection search. All fall inside the broader search above, and all have NCBI-generated raw counts. None has been run yet.

| Dataset | System and comparison | Samples and design | Expected signals (hypotheses) |
|---|---|---|---|
| **GSE164287** (top pick) | *M. tuberculosis* H37Rv vs uninfected, primary human monocyte-derived macrophages, 48 h | 8: 4 donors × 2 conditions, paired | NF-κB transcription factors (NFKB1, RELA); PROGENy NFkB and TNFa; Hallmark TNFα-via-NF-κB and inflammatory response up |
| **GSE193977** | Chikungunya virus vs uninfected, monocyte-derived macrophages, 6 h and 24 h | 12: 4 donors × 3 conditions; 24 h vs uninfected uses 8, paired | STAT1, STAT2 and IRF factors; PROGENy JAK-STAT; Hallmark interferon-α and -γ responses up |
| **GSE128661** | *Aspergillus fumigatus* vs uninfected, monocyte-derived macrophages, 2 h and 6 h | 9: 3 donors × 3 conditions; 6 h vs uninfected uses 6, paired | Early NF-κB/TNF response, stronger at 6 h than 2 h |
| **GSE153970** | SARS-CoV-2 vs mock, primary human airway epithelial cultures, 48 h | 6: 3 vs 3, unpaired | Interferon (STAT/IRF) and inflammatory signalling. The GEO summary reports cytokine and chemokine induction, including IL-6. |
| **GSE198542** | Zika virus vs mock, stem-cell-derived macrophages from an IFNAR2-deficient patient and from healthy controls | 12: 2 genotypes × 2 conditions × 3; run as two separate 3 vs 3 comparisons | Interferon/STAT activity in control cells, absent or much weaker in IFNAR2-deficient cells. This is a built-in negative control. |

**Set aside:**

- **GSE98665:** one sample per condition and time point, and no raw counts.
- **GSE279463:** no NCBI raw counts, so its raw files would need custom processing.

## Vignette 2: CellTypist on single-cell RNA-seq

### Find a tool in NDE's Computational Tools

Searching NDE's tool records for "cell type" together with annotation, classification or identification returns 153 tools. CellTypist is one of them (51st by relevance), with the record [biotools\_celltypist](https://data.niaid.nih.gov/resources?id=biotools_celltypist).

- **Which implementation:** the record points to github.com/Teichlab/celltypist, the same repository the server wraps (version 1.7.1).
- **What it doesn't say:** almost everything else. The record has a one-line description ("A tool for semi-automatic cell type classification") and no paper, version, license, topics, or input and output data types. The repository itself cites [Domínguez Conde et al. 2022](https://doi.org/10.1126/science.abl5197).

### Wrap it as an MCP server

The companion doc built five tools from CellTypist's main and multi-label tutorials: list the built-in models, annotate cells, train a custom model, inspect a model, and compare predicted labels with reference labels. Run on the 10x PBMC 3k dataset, its labels agreed with the scanpy tutorial's clusters for 97.7% of 2,638 cells.

Its input contract, and where each requirement can be checked:

| Requirement | What the server needs | Where it is checked |
|---|---|---|
| Assay | Single-cell or single-nucleus RNA-seq | NDE search |
| Organism | Any species with a built-in model (48 human, 10 mouse, 2 macaque, 1 pig), or any species via a custom-trained model | GEO sample sheet, because NDE species tags can be wrong or missing |
| Cell count | More than 50 cells for majority voting; no minimum otherwise | File contents |
| Labels | None needed to annotate; author cell labels needed to score agreement | Sample sheet or file contents |
| Data files | Expression matrix with gene symbols, log-normalised to 10,000 counts per cell (raw counts can be converted) | File listing, then file contents |

### Search NDE for compatible datasets

```
@type:"Dataset" AND ("single-cell" OR "single cell" OR "single-nucleus" OR "single nucleus" OR scRNA-seq OR snRNA-seq OR scRNAseq OR snRNAseq)
AND (RNA OR transcriptome OR transcriptomic OR transcriptomics OR "gene expression" OR scRNA-seq OR snRNA-seq OR scRNAseq OR snRNAseq)
```

Why each filter is there, and what is deliberately left out:

- **Single-cell terms:** NDE has no assay term for single-cell sequencing, so the search uses words in the record.
- **RNA terms:** these keep RNA studies and drop most single-cell chromatin, methylation or cytometry studies.
- **No species filter:**
  - CellTypist's built-in models cover human, mouse, macaque and pig, and the server can train a model for any other species from labelled data.
  - 18,208 of the hits have no species tag at all. Both Zenodo exemplars below are among them.
- **No tissue or disease filter:** the built-in models cover immune cells, lung, airway, gut, brain, heart, liver, pancreas, skin, endometrium and more.
- **No file-format filter:** most single-cell matrices are published as 10x Genomics or plain-text files that one standard function call converts to the AnnData format CellTypist reads. Records that mention AnnData or .h5ad files are the ready-to-run subset.

| Slice | Records |
|---|---|
| All single-cell RNA-seq records | **117,912** |
| Tagged human | 43,409 |
| Tagged mouse | 49,532 |
| No species tag | 18,208 |
| Mentioning AnnData or .h5ad files | 361 |

Of the human-tagged records, 17,895 are from GEO, 11,654 from BioProject, 9,230 from SRA, 2,590 from Figshare and 548 from Zenodo; the rest come from a dozen other repositories. GEO, BioProject and SRA records often describe the same study, so these counts overstate the number of distinct studies.

### Spot-check a random sample

Of 20 random human-tagged records, 12 were SRA or BioProject records. 10 of those 12 link to a GEO series, which we checked in their place.

- **16 have a downloadable expression matrix:**
  - 10x Genomics files (.h5 or .mtx), CSV or plain-text tables;
  - 4 inside per-sample archives that we did not open;
  - none as AnnData.
  - One of the 16 is a mouse-only series tagged as human.
- **2 had sequencing reads only:** no linked GEO series was found.
- **2 were not single-cell RNA-seq:** one is a DNA-methylation and chromatin study; the other is a bulk RNA-modification study whose description mentions a related single-cell series.

Most of the 16 contain immune cells or come from a tissue with a built-in model. A few, such as stem-cell-derived cultures, ear cartilage and mouse retina, would need a custom model. Only one series had author cell-type labels visible in its GEO sample sheet; we did not open the matrices to look for label files.

That gap is the point. CellTypist would add consistent cell-type labels to datasets that lack them. Scoring agreement, however, needs the labelled minority.

### Confirm with exemplars

These four were picked by hand from an earlier search for labelled immune single-cell data. All fall inside the broader search above. None has been run yet.

| Dataset | Contents | Labels | Status |
|---|---|---|---|
| **Zenodo 20988712** (top pick; CC BY 4.0) | COVID-19 and healthy PBMC from two source studies. 29,989 cells × 14,881 genes, 8 donors (4 COVID-19, 4 healthy). Raw integer counts with gene symbols. 94 MB. | 5 lineages and 26 cluster labels; donor and severity columns | **Opened and checked** |
| **Zenodo 19581815** (CC BY 4.0) | Healthy-donor PBMC from 10x Genomics' Multiome PBMC 10k. 9,631 cells × 29,095 genes. Raw integer counts. 51 MB. | 19 fine-grained labels from label transfer against the Hao et al. 2021 PBMC atlas | **Opened and checked** |
| **GEO GSE206289** | PBMC from SARS-CoV-2-infected infants (26) and matched controls (14), per the GEO description. A 1.0 GB .h5ad. | Not yet checked | Not opened |
| **GEO GSE233321** | PBMC from 167 healthy people aged 2 months to 105 years, per the GEO description. One file per cell population. | Not yet checked | Not opened |

**Set aside:**

- one dataset with no labels or sample metadata;
- one with only normalised values;
- three germinal-centre B-cell sets of 0.6–9.6 GB per file;
- one 12.7 GB archive;
- two lung and airway atlases that likely overlap with CellTypist's training data.

## What it means for NDE

### Linking Computational Tools to MCP servers to datasets

NDE holds 34,275 tool records and 5.7 million dataset records, and almost nothing connects them. On October 28, 21 dataset records pointed to a tool and no tool record pointed to a dataset. Most tool records don't say what data the tool takes (4,107 list input data types), and dataset records don't say what could be run on them. Both vignettes started from a tool record that could not answer "which datasets does this fit?"

A Paper2MCP server fills that gap. It is a runnable, tested form of the tool, and its input contract is the missing description of the data the tool needs. That gives a chain NDE could record and display: **tool record → MCP server → compatible datasets → results**.

- **On a tool record:** the servers that wrap it, and how many datasets are predicted compatible, verified executable and verified sensible.
- **On a dataset record:** the tools predicted to fit it, the ones that have been run, and what they found.
- **Results as searchable metadata:** cell-type labels and regulator activities would let people search for datasets containing a given cell type, or comparisons in which a given regulator is activated.
- **More precise tool records:** the decoupleR record describes the R package, while the server wraps the Python one; neither the Python package nor PyDESeq2 has a record, and the CellTypist record lacks a paper, version and license. Links from servers to tools need records for the exact implementation and version that was run.
- **More tools:** the two tool searches above returned 75 tools for transcription-factor activity and 153 for cell-type annotation. Each one with public code and tutorials is a candidate for the same treatment.

### Croissant-style annotations for more accurate matching

Today "predicted compatible" rests on dataset-level metadata, and the spot checks show where that falls short:

- species tags inferred from linked papers were wrong for 3 of the 40 records checked;
- which samples form the groups to compare is not recorded, so it had to be read from GEO sample sheets;
- what the files contain (raw or normalised values, which gene identifiers, whether cell labels are included) is unknown until the files are listed or opened.

[Croissant](https://docs.mlcommons.org/croissant/docs/croissant-spec-1.1.html), the MLCommons metadata format, describes the files inside a dataset and the fields in each, with checksums. With descriptions like that, most of the checks we ran by hand would become metadata queries, and NDE could compute compatibility across the whole catalog with fewer false matches.

- **Gaps in Croissant 1.1:** it can describe columns in CSV or JSON files, which covers count tables like NCBI's, but it has no standard way to describe what is inside an AnnData or HDF5 file, where single-cell matrices, cell labels and gene names live. Its data types also don't say whether values are raw or normalised, or which gene-identifier system is used; shared vocabulary terms would be needed for that.
- **Declared is not verified:** a description can be wrong, just as species tags are now. Descriptions generated by reading the files are more trustworthy than ones drafted from a landing page.
- **The opportunity:** an agent that runs a tool reads the files anyway. Each run could write back a file-level description of what it found, plus the steps used to adapt the data. Those descriptions would come from the files themselves, and each one would improve matching for every other tool.
- **The tool side:** NDE's tool records already have input and output fields, filled in for 4,107 of 34,275 tools with terms such as the EDAM data type "Transcription factor name". Both records here leave them empty. Filling them from each server's input contract, in the same vocabulary as the dataset descriptions, would make matching a tool to datasets a metadata query.
