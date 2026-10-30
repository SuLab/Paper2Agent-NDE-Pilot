# Step 5 exemplar runs (cap: 3 per tool) — status

Paused 2026-10-29 at the user's request (computer going offline). Nothing is running.

## Setup
- Servers re-extracted from the delivered ZIPs (sha256 94ff548c… decoupler, a52110da… CellTypist) into `servers/`, installed with uv from `src/requirements.txt`; `uv pip check` clean (128 / 122 packages).
- All tool calls go over stdio to those servers via `scripts/mcp_call.py` (same harness as `../new-data`).
- Expectations are written before each run (`*/predeclared_expectations.json`; decoupler hashes in `decoupler/EXPECTATIONS.sha256`).

## decoupler — DONE (3/3 verified executable, 2/3 verified sensible)
| Exemplar | Design | Result |
|---|---|---|
| GSE164287 Mtb vs uninfected MDM | ~donor + condition, 4 v 4 | E1–E4 pass. NFKB1 #1, RELA #3 of 624 TFs; IL1B/IL6/CXCL8/PTGS2 up; 4,295 DE genes |
| GSE153970 SARS-CoV-2 vs mock airway epithelium | ~condition, 3 v 3 | E1–E4 pass. NFKB1 #1, RELA #9, ATF4 #20 of 664; IL6/TNF/CXCL8 up; no IFN response, as the paper reports. The NDE doc's interferon hypothesis was replaced before running, because the paper's abstract reports no IFN response. |
| GSE198542 ZIKV vs mock, WT and IFNAR2-deficient iPSC macrophages | ~condition per genotype, 3 v 3 each | E1–E3 pass (WT ISGF3 strongly active). E4 negative control **fails as declared**: the ISGF3 mean score is 55% of WT (the threshold was <50%), and MX1's log2FC is as large in MUT as in WT. Post hoc: infected-cell MX1 is 83 vs 3,097 CPM, and MUT baselines are near zero, so the defect is real but a fold-change readout misses it. See `post_hoc_not_scored` in `expectation_results.json`. |

Dropped from the NDE doc's 5 candidates: GSE193977 (CHIKV, which duplicates the viral/macrophage/IFN case) and GSE128661 (Aspergillus, which duplicates the NF-κB/macrophage case with only 3 donors).

## CellTypist — IN PROGRESS (data downloaded, nothing run)
Chosen: Zenodo 20988712 (`test_RNA.h5ad`, md5 OK), Zenodo 19581815 (`10x-Multiome-Pbmc10k-RNA.h5ad`, md5 OK), and GEO GSE206289 (`GSE206289_pCoV40_cleaned_12062021.h5ad.gz`, 1,109,227,079 bytes, matches GEO's Content-Length).
Dropped: GSE233321. Its all-PBMC file is 12.5 GB, and the per-population files are 1–2 GB each.

Already inspected (labels only, no predictions):
- **Zenodo 20988712:** 29,989 cells × 14,881 genes of raw integer counts with symbols. obs `large` has 5 lineages (T, Mono, B, Megakaryocyte, Neutrophil), with no NK or DC lineage. obs `celltype` has 26 lineage-prefixed clusters. Also `sample` (8 donors: 4 CoV, 4 HD), `sev.level` and `batch` (Wilk/Aruna).
- **Zenodo 19581815:** 9,631 cells × 29,095 genes of raw integer counts with symbols. obs `cell_type` has 19 Azimuth level-2 labels (CD14 Mono 2,554, CD4 Naive 1,382, … Plasma 12).

Next steps:
1. `gunzip -k` the GSE206289 file (it was interrupted at pause; the partial output was removed). Inspect its obs and labels with `python -I`, and read the GEO description.
2. List the exact Immune_All_High and Immune_All_Low label names: run the `celltypist_inspect_model` tool, or read `../new-data/pbmc3k/results_step1.json`.
3. Write `predeclared_mapping.json` for each dataset before annotating. Planned: score at High level for all three; for 19581815 also score at Low level; for 20988712 the T lineage accepts T cells or ILC, because there is no NK lineage. Possible E2 for 20988712: more plasmablasts in CoV than HD donors (Wilk 2020).
4. Normalise each to log1p CP10k (`scripts/prep_*`), then run annotate (majority voting) and compare over MCP, and score.
5. Write up the tier counts for both tools.
