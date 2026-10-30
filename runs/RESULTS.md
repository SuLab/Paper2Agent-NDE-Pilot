# Paper2MCP: decoupler and CellTypist, re-run on new data

Date: 2026-10-27. Each project was built with the Paper2Agent `paper2mcp` workflow. Each one:
- passed every workflow gate, from setup through `complete`;
- had independent verifier agents check its tools against real tutorial runs;
- was delivered as a ZIP that a separate agent re-validated after extracting it to a new location.

## Deliverables

| Server | ZIP | Tools |
|---|---|---|
| decoupler (bulk RNA-seq: filter, DE, TF/pathway activity) | `decoupler/dist/decoupler-mcp.zip` (sha256 94ff548c…) | `decoupler_filter_genes_by_expression`, `decoupler_rank_pcs_by_metadata`, `pydeseq2_test_differential_expression`, `decoupler_score_activities_ulm`, `decoupler_inspect_source_targets` |
| CellTypist (scRNA-seq cell-type annotation) | `celltypist/dist/celltypist-mcp.zip` (sha256 a52110da…) | `celltypist_list_models`, `celltypist_annotate_cells`, `celltypist_train_model`, `celltypist_inspect_model`, `celltypist_compare_with_reference` |

Each ZIP contains a `USAGE.md` covering installation, startup, every tool's inputs and outputs, and the validation limits.

### Validation

| | decoupler | CellTypist |
|---|---|---|
| Tutorial reproduced end to end, with replay identical | `bulk/rna.ipynb` | main and multi-label tutorials |
| Independent verifier tests | 74/74 | 90/90 (40 + 50) |
| Code repairs needed | 0 | 0 |
| Real-stdio acceptance cases (project env / fresh env / extracted ZIP) | 27/27 in each | 33/33 in each |

## New dataset 1: airway smooth muscle + dexamethasone → decoupler MCP

**Data.** GEO GSE52778, using NCBI's raw counts mapped to gene symbols. I kept the 4 untreated and 4 dexamethasone samples from 4 donor cell lines; the albuterol arms were dropped.

**Run.** All seven MCP calls went over stdio to the extracted, delivered server, with live network downloads:
- filter genes;
- rank PCs against metadata;
- PyDESeq2 with `~cell_line + condition`;
- ulm scoring against CollecTRI, PROGENy and Hallmark;
- inspect NR3C1.

**Expectations.** I declared these in `new-data/airway/predeclared_expectations.json` before running. Two sources:
- The study that produced the data. According to PubMed, Himes *et al.* 2014, *PLoS One*, [DOI](https://doi.org/10.1371/journal.pone.0099625).
- Glucocorticoid-receptor biology.

| Expectation | Result |
|---|---|
| E1. DUSP1, KLF15, PER1, TSC22D3, CRISPLD2, C7 and CCDC69 have padj < 0.05 (the paper's named genes) | **Pass**: all 7; DUSP1 ranks #5 genome-wide, KLF15 #9, PER1 #11, CCDC69 #14 |
| E2. DUSP1, KLF15, PER1, TSC22D3, CRISPLD2, FKBP5 and ZBTB16 are induced | **Pass**: log2FC from +2.6 (CRISPLD2) to +7.4 (ZBTB16) |
| E3. The glucocorticoid receptor NR3C1 is significantly *activated* | **Pass**: ulm score +5.85, padj 1.6e-6, **rank 2 of 657 TFs**. NR3C1's own mRNA *falls* (log2FC −0.99), so the activity call comes from its targets. The positive leading edge includes MAOA, DUSP1, KLF15, PER1, MT2A, IRS2, FOXO1, KLF9 and CEBPD. |

**Reported only (no pass/fail):**
- **PCA:** PC1 (33% of variance) tracks treatment, and PCs 2–4 track donor. That supports putting `cell_line` in the design.
- **DE count:** 4,441 genes at padj < 0.05, against 316 in the paper. The pipelines differ (aligner, annotation, test and filtering), so the counts aren't comparable.
- **TF ranking:** PGR is #1 and NR3C1 #2. *Interpretation, not tested:* PGR and GR share the same response elements. In the CollecTRI network used here, 44 of PGR's 133 targets are also NR3C1 targets. Activity inferred from expression alone may not separate these steroid receptors.
- **PROGENy:** the top pathway is "Androgen" (same caveat).
- **Hallmark:** the top terms are ADIPOGENESIS and TNFA_SIGNALING_VIA_NFKB, both positive. The TNFA set contains classic GR-induced genes (DUSP1, KLF9, PER1, SGK1, CEBPD, KLF6), so a positive score here does not mean NF-κB activation.

## New dataset 2: PBMC 3k → CellTypist MCP

**Data.** 10x PBMC 3k: 2,638 cells with the 8 cluster labels from the scanpy tutorial. Raw counts were recovered from the public file and log1p-normalised to 10k counts per cell, which is CellTypist's input contract.

**Run.** Four MCP calls to the extracted, delivered server, starting from an empty model cache. Only the two models requested were downloaded.

**Mapping.** I declared which CellTypist label counts as a match for each author label before annotating (`new-data/pbmc3k/predeclared_mapping.json`).

| Author label | Cells | Immune_All_High + majority voting agrees |
|---|---|---|
| CD4 T cells | 1,144 | 99.4% (T cells) |
| CD14+ Monocytes | 480 | 99.0% (Monocytes) |
| B cells | 342 | 98.5% |
| CD8 T cells | 316 | 91.8% (25 called ILC) |
| NK cells | 154 | 89.6% (ILC; 16 called T cells) |
| FCGR3A+ Monocytes | 150 | 100% |
| Dendritic cells | 37 | 100% |
| Megakaryocytes | 15 | 14/15 |
| **All** | **2,638** | **97.7%** (per-cell labels without voting: 96.4%) |

**Fine-grained model (Immune_All_Low, 98 types).** Most clusters map cleanly:
- CD4 T cells → naive/central-memory helper T cells
- NK cells → CD16+ NK cells
- CD14+ monocytes → classical monocytes (404) or non-classical monocytes (71)
- FCGR3A+ monocytes → non-classical monocytes
- DCs → DC
- Megakaryocytes → megakaryocytes/platelets

The author "CD8 T cells" cluster is the exception. It splits into "NK cells" (92), effector-memory/resident cytotoxic T cells (69), MAIT cells (61) and naive/central-memory helper T cells (58), all with low probabilities.

## Limits

- **Tool scope.** The tools wrap the selected tutorials only. They do not cover all of decoupler or CellTypist, and they do not reproduce the papers' benchmarks.
- **Paper skills.** None were built, so these are MCP servers only.
- **Tested environment.**
  - Only Linux x86_64 (WSL2) was tested.
  - Installs used a warm uv cache.
  - Built-in CellTypist models and the decoupler networks are downloaded live, and can change.
- **Majority voting.** Exact majority-voting labels depend on the pinned scanpy, leidenalg and igraph versions.

## Where things are

- Projects, with all evidence, reports and workflow state: `/home/asu/Science/paper2agent-runs/{decoupler,celltypist}/`
- New-data inputs, preparation scripts, MCP call logs and outputs: `/home/asu/Science/paper2agent-runs/new-data/`
