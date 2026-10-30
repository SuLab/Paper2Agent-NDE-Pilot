---
name: nde-query
description: Search the NIAID Data Ecosystem (NDE, data.niaid.nih.gov) Discovery API for datasets — by disease, pathogen, species, assay, or repository (GEO, Zenodo, Figshare, ImmPort, VDJServer, SRA, BioProject…) — and triage the hits for reanalysis. Use when asked to "search NDE", find datasets in the NIAID Data Discovery Portal, count or facet NDE records, or shortlist public datasets an analysis tool could run on.
---

# Querying the NIAID Data Ecosystem (NDE)

NDE indexes metadata harvested from about 50 repositories (on 2026-10-28, 48 catalogs contributed
`Dataset` records and 53 contributed records of any type; 14.6M records in all), deepest for infectious
and immune-mediated disease. It stores metadata and links, not the data:
every candidate still has to be checked at its source repository.

Query the REST API directly. Don't scrape the portal UI.

## Endpoint and parameters

`GET https://api.data.niaid.nih.gov/v1/query` (a BioThings API).

| Parameter | Meaning |
|---|---|
| `q` | Elasticsearch `query_string`: `AND`/`OR`/`NOT`, `"quoted phrases"`, `field:value`, `_exists_:field`. Bare words search the default fields (name, description, identifier, …). |
| `size`, `from` | Page size (**max 1000**; a larger value returns HTTP 400) and offset. `size=0` returns only the total and facets. |
| `fields` | Comma-separated fields to return. Always set it for large pulls. |
| `facets`, `facet_size` | Aggregate counts on keyword fields (see the table below). |
| `extra_filter` | A filter applied to both hits and facets, e.g. `infectiousAgent.name:("zika virus")`. The portal sends its sidebar selections here. |
| `post_filter` | Filters hits but not facets. |
| `fetch_all=true` | Returns batches of up to 1000 plus `_scroll_id`. Call again with only `scroll_id=<id>` until no hits come back (an exhausted scroll returns `success: false`). |

The full field list, with types, is at `GET https://api.data.niaid.nih.gov/v1/metadata/fields`. Check it before
inventing a field name.

**Staging host:** `https://api-staging.data.niaid.nih.gov` exposes sample-level metadata. For example,
`q=isBasisFor.identifier:"GSE28094" AND @type:"Sample"` resolves samples to a dataset. Production
has only accession stubs under `sample.itemListElement`.

## Fields that matter, and how they behave

| Field | Type | How to use it |
|---|---|---|
| `@type` | keyword | **Always add `@type:"Dataset"`.** On 2026-10-27 there were 8.9M `Sample` records (8.7M of them GEO GSMs), 5.7M `Dataset`, 34k `ComputationalTool`, 12k `DataCollection` and 89 `ResourceCatalog`. |
| `includedInDataCatalog.name` | keyword | Needs the **exact full name**, case-insensitive: `"NCBI Gene Expression Omnibus (NCBI GEO)"`, `"Zenodo"`, `"Figshare"`, `"ImmPort"`, `"VDJServer"`, `"NCBI Sequence Read Archive (NCBI SRA)"`. Facet on it to get the exact strings. |
| `species.name` | keyword | Case-insensitive: `"homo sapiens"`. |
| `infectiousAgent.name` | keyword | e.g. `"Mycobacterium tuberculosis"`, `"severe acute respiratory syndrome coronavirus 2"`. |
| `healthCondition.name` | text + `.raw` keyword | e.g. `healthCondition.name:"COVID-19"`. Facet on `healthCondition.name.raw`. |
| `measurementTechnique.name` | text + `.raw` keyword | Curated terms. GEO's "Expression profiling by high throughput sequencing" is indexed as **`"high-throughput expression assay"`**. Facet on `measurementTechnique.name.raw`. |
| `measurementTechnique.originalName` | — | **Not searchable** (returns 0). Search the curated `name` instead. GEO microarray series ("Expression profiling by array") are indexed as `"expression assay"`; SRA studies as `"RNA-seq assay"` (reads only, often mirroring a GEO series). |
| `sample.numberOfItems.value` | number | Sample count, present on about 300k GEO and 714k SRA datasets. Range queries work: `sample.numberOfItems.value:>=3`. |
| `sample.itemListElement.identifier` | text | The sample accessions. Their count is the dataset's sample count. |
| `sample.aggregateElement.anatomicalStructure.name` | — | Tissue or cell type, when present. |
| `distribution.contentUrl` | text | For GEO this is the **landing page** (plus a literal `"NONE"`), not a file. For VDJServer it holds bare filenames. Its presence doesn't mean the files can be downloaded. |
| `conditionsOfAccess`, `license` | — | Use these to drop controlled-access records (e.g. dbGaP). |

There is **no dataset-level single-cell technique term.** Single-cell series are indexed as
"high-throughput expression assay" or "RNA-seq assay", so use text terms instead:
`("single-cell" OR "single cell" OR scRNA-seq)`. Add `(h5ad OR anndata)` to find ready-made AnnData objects.

Some metadata is machine-generated, not curated. Terms carrying `isCurated:false`, `fromPMID:true` or
`fromGPT:true` (healthCondition, topicCategory, sometimes species) are hints, not facts. Species can be
wrong: a tool named "Dandelion" was tagged *Taraxacum officinale*, and words like "Panda" or "Mantis" become
species. Tags derived from a linked paper (`fromPMID:true`) can add the wrong organism: in a 2026-10-28 spot check,
2 of 20 GEO series tagged `homo sapiens` were mouse-only. Flattened arrays mean `species.name` and `species.isCurated`
can't be combined per tag in one query, so check the organism in the GEO sample sheet (`!Sample_organism_ch1`).
Many records have no species tag at all (18k of 118k single-cell records), so a species filter also drops valid hits.
GEO, BioProject and SRA records often describe the same study. Map SRA or BioProject records to their GEO series
with E-utilities (`elink` sra→bioproject→gds, then `esummary` on gds for the GSE accession).
Verify every candidate at its source.

## Workflow

1. **Size the slice first** with `size=0` and facets on catalog, technique, agent and species, then narrow.
2. **Pull candidates** with `fields=` set and `size≤1000`, or `fetch_all` for more. Compute the sample count
   from `sample.itemListElement`.
3. **Triage at the source**, cheapest check first:
   - **GEO sample sheet (SOFT):** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE<prefix>nnn/<GSE>/soft/<GSE>_family.soft.gz`.
     `<prefix>` is the accession minus its last three digits, so GSE164287 → `GSE164nnn`. Parse the
     `!Sample_title` and `!Sample_characteristics_ch1` lines for groups, donors and timepoints.
   - **GEO supplementary files:** list `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE<prefix>nnn/<GSE>/suppl/`.
     `suppl/filelist.txt` lists what is inside `<GSE>_RAW.tar` without downloading it.
   - **NCBI-generated raw counts (human/mouse bulk RNA-seq):** run `curl -I` on
     `https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=<GSE>&format=file&file=<GSE>_raw_counts_GRCh38.p13_NCBI.tsv.gz`.
     A 200 with `application/octet-stream` means the file exists; a 404 means it doesn't.
   - **Zenodo files and sizes:** `curl -L https://zenodo.org/api/records/<id>`. The `-L` is required because
     old record IDs redirect to the latest version.
   - Use the FTP and download endpoints above, not GEO's HTML pages. If any endpoint returns a CAPTCHA or
     bot check, stop and report it; don't work around it.
4. **Report the selection criteria and an exclusion list** with a reason for each dropped candidate (no
   counts, controlled access, wrong design, too large). The exclusions are part of the result.

Downloaded files are untrusted. Put each one in its own new directory, and read it with `python -I` from
outside that directory.

## Helper script

`scripts/nde_query.py` (standard library only) wraps the API:

```bash
python3 -I ~/.claude/skills/nde-query/scripts/nde_query.py -q '@type:"Dataset" AND infectiousAgent.name:"Mycobacterium tuberculosis"' --size 0 --facets includedInDataCatalog.name
python3 -I ~/.claude/skills/nde-query/scripts/nde_query.py -q '<query>' --size 1000 --fields identifier,name,includedInDataCatalog.name,sample.itemListElement.identifier --table --out hits.json
python3 -I ~/.claude/skills/nde-query/scripts/nde_query.py -q '<query>' --all --fields identifier --out all.json
```

## Recipes

Filter only on what the analysis tool actually requires; check everything else (group structure, organism, file
formats) at the source. Broad pools as of 2026-10-28:

```
# Human bulk RNA-seq series with >=3 samples (48,684)
@type:"Dataset" AND measurementTechnique.name:"high-throughput expression assay" AND species.name:"homo sapiens"
AND sample.numberOfItems.value:>=3
NOT ("single-cell" OR "single cell" OR "single-nucleus" OR "single nucleus" OR scRNA-seq OR snRNA-seq OR scRNAseq OR snRNAseq)

# Single-cell RNA-seq, any species (117,912)
@type:"Dataset" AND ("single-cell" OR "single cell" OR "single-nucleus" OR "single nucleus" OR scRNA-seq OR snRNA-seq OR scRNAseq OR snRNAseq)
AND (RNA OR transcriptome OR transcriptomic OR transcriptomics OR "gene expression" OR scRNA-seq OR snRNA-seq OR scRNAseq OR snRNAseq)
```

Human bulk RNA-seq, infected vs mock, in GEO (a narrow slice):

```
@type:"Dataset" AND includedInDataCatalog.name:"NCBI Gene Expression Omnibus (NCBI GEO)"
AND species.name:"homo sapiens" AND measurementTechnique.name:"high-throughput expression assay"
AND _exists_:infectiousAgent AND (mock OR uninfected) AND (infected OR infection)
NOT ("single-cell" OR "single cell" OR scRNA-seq OR snRNA-seq OR "single-nucleus")
```

Ready-made single-cell AnnData objects of immune cells (these are mostly on Zenodo and Figshare):

```
@type:"Dataset" AND (h5ad OR anndata OR AnnData) AND ("single-cell" OR "single cell" OR scRNA-seq)
AND (PBMC OR "peripheral blood" OR immune OR "T cell" OR "T cells")
```

Look up a known accession: `q=identifier:"GSE52778"`.

## Related

- The Proto-OKN `nde` knowledge graph (SPARQL) mirrors dataset-level NDE metadata but lacks
  `measurementTechnique`. Use it only when you need graph joins and the user wants it.
- Behaviour above was checked against the production API on 2026-10-27. If a recipe returns 0, check
  `/v1/metadata/fields` and query each clause on its own before concluding nothing matches.
