# Paper2Agent × NIAID Data Ecosystem pilot

Can tools built by [Paper2Agent](https://github.com/jmiao24/Paper2Agent) add a layer of analysis-based annotation to datasets in the [NIAID Data Ecosystem](https://data.niaid.nih.gov) (NDE)?

This pilot used Paper2Agent's Paper2MCP workflow to turn two published methods into MCP servers that an AI agent can call:

- **decoupler** (Python, [Badia-i-Mompel et al. 2022](https://doi.org/10.1093/bioadv/vbac016)): transcription-factor and pathway activity between two groups of bulk RNA-seq samples;
- **CellTypist** ([Domínguez Conde et al. 2022](https://doi.org/10.1126/science.abl5197)): cell-type labels for single-cell RNA-seq.

It then:

1. ran each server on a new public dataset;
2. searched NDE for every dataset that matches each server's input contract;
3. spot-checked random hits at their source repositories;
4. ran hand-picked exemplars from NDE through the servers.

The write-up is in [`docs/nde-pilot.md`](docs/nde-pilot.md), a snapshot of the working Google Doc as of 2026-10-28.

## Status (2026-10-29)

| Step | decoupler | CellTypist |
|---|---|---|
| MCP server built and validated | Done: 74/74 verifier tests; 27/27 stdio acceptance cases in each of three environments | Done: 90/90 verifier tests; 33/33 stdio acceptance cases in each of three environments |
| New dataset | GSE52778 (airway, dexamethasone): 3/3 pre-declared expectations pass | PBMC 3k: 97.7% agreement with author labels |
| NDE search | 48,684 human bulk RNA-seq series with at least 3 samples | 117,912 single-cell RNA-seq records |
| Random spot check (20 each) | 14 usable, 1 unclear, 5 not usable | 16 with a downloadable matrix |
| NDE exemplars | 3 run: all 3 verified executable; 2 of 3 met every pre-declared expectation (see [`runs/exemplars/STATUS.md`](runs/exemplars/STATUS.md)) | 3 downloaded, none annotated yet (paused) |

The exemplar runs postdate the doc snapshot, which still says none has been run.

## Using one of the MCP servers

You don't need to clone decoupler or CellTypist. Installing a server's requirements fetches the tool itself at the tested commit. You need this repo, [uv](https://docs.astral.sh/uv/), git, and an MCP client such as [Claude Code](https://code.claude.com/docs/en/overview). Linux x86_64 is the only tested platform.

### 1. Clone this repo and pick a server

```bash
git clone https://github.com/SuLab/Paper2Agent-NDE-Pilot.git
cd Paper2Agent-NDE-Pilot/mcp-servers/decoupler     # or mcp-servers/celltypist
```

Run every later command from that folder.

### 2. Create the server's environment

The two servers need different Python versions. If yours is missing, uv downloads it. The first install can take a few minutes; each environment is under 1 GB.

**decoupler** (Python 3.13):

```bash
uv venv .venv --python 3.13
uv pip install --python .venv/bin/python -r src/requirements.txt
uv pip check --python .venv/bin/python
```

**CellTypist** (Python 3.12):

```bash
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python -r src/requirements.txt
uv pip check --python .venv/bin/python
```

### 3. Register it with Claude Code

**decoupler:**

```bash
claude mcp add decoupler --scope user -- "$PWD/.venv/bin/python" "$PWD/src/decoupler_mcp.py"
```

**CellTypist:**

```bash
claude mcp add celltypist --scope user -e CELLTYPIST_FOLDER="$PWD/celltypist_home" -- "$PWD/.venv/bin/python" "$PWD/src/celltypist_mcp.py"
```

- **`--scope user`** makes the server available in every project. Leave it out to register the server for the current directory only.
- **`CELLTYPIST_FOLDER`** is where CellTypist caches the models it downloads.
- **Checking it worked:** `claude mcp list` should show the server as connected. Inside Claude Code, `/mcp` shows the same.
- **Other MCP clients:** each server's `USAGE.md` has a JSON configuration block.

### 4. Put your data in the format the server expects

Both servers read an AnnData `.h5ad` file. Run the conversion with the server's own Python (`.venv/bin/python`), which already has the packages it needs.

**decoupler** expects samples as rows and genes (HGNC symbols) as columns, raw integer counts, and the sample metadata in `.obs`: a condition column, plus donor or batch columns for paired designs. From a counts table and a sample sheet:

```python
import anndata as ad
import pandas as pd

counts = pd.read_csv("counts.csv", index_col=0)   # genes x samples, raw integer counts, gene symbols as row names
meta = pd.read_csv("samples.csv", index_col=0)    # one row per sample: condition, donor, ...
adata = ad.AnnData(counts.T.loc[meta.index], obs=meta)
adata.write_h5ad("my_counts.h5ad")
```

**CellTypist** expects cells as rows, gene symbols as columns, and values log-normalised to 10,000 counts per cell. From 10x Genomics output:

```python
import scanpy as sc

adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")   # or sc.read_10x_mtx("filtered_feature_bc_matrix/")
adata.var_names_make_unique()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad("my_cells.h5ad")
```

For worked examples, see [`runs/new-data/scripts/prep_airway.py`](runs/new-data/scripts/prep_airway.py) (NCBI GEO counts) and [`prep_pbmc3k.py`](runs/new-data/scripts/prep_pbmc3k.py) (10x PBMC 3k). [`min-samples/`](min-samples) gives the smallest designs each tool handles.

### 5. Ask Claude

Name the file, the comparison and an output folder. For example:

> Use the decoupler tools on /data/my_counts.h5ad. Filter lowly expressed genes, test treated vs control in the `condition` column with the design `~donor + condition`, then score transcription-factor activity with CollecTRI and pathway activity with PROGENy. Write all outputs under /data/decoupler_out.

> Use the celltypist tools to annotate /data/my_cells.h5ad with the Immune_All_Low model and majority voting, then compare the predictions with the `cluster` column. Write all outputs under /data/celltypist_out.

**Where outputs go:**
- Each tool call writes a fresh folder of tables and plots, and returns their paths.
- Without an output folder, decoupler writes under `mcp-servers/decoupler/tmp/outputs/`, and CellTypist under the system temp folder.

**Downloads during the first calls:** decoupler's regulator networks and CellTypist's models are fetched from the internet on first use.

**Reference:** each server's `USAGE.md` documents every tool's arguments and outputs, and the validation behind it.

**Without Claude:** [`runs/new-data/scripts/mcp_call.py`](runs/new-data/scripts/mcp_call.py) calls a server over stdio from a JSON plan, the way an MCP client would; the plans under `runs/` are examples.

**What was tested:** on 2026-10-29, steps 1, 2 and 4 were run from a fresh clone. Each server then ran a short analysis over stdio:
- **decoupler on GSE52778:** 4,441 differentially expressed genes, as in `runs/RESULTS.md`.
- **CellTypist on PBMC 3k:** all 2,638 cells annotated.

The `claude mcp add` commands follow `claude mcp add --help` but were not run.

## Layout

| Path | Contents |
|---|---|
| [`docs/`](docs) | The NDE pilot write-up (Markdown snapshot). |
| [`mcp-servers/decoupler/`](mcp-servers/decoupler), [`mcp-servers/celltypist/`](mcp-servers/celltypist) | Paper2MCP output for each tool: server source (`src/`), the delivered ZIP (`dist/`), `USAGE.md`, workflow and validation reports (`reports/`), verifier test code and logs (`tests/`), and the executed tutorial notebooks. |
| [`mcp-servers/build/`](mcp-servers/build) | Helpers used while coordinating the Paper2MCP build, and its agent lifecycle log. |
| [`runs/RESULTS.md`](runs/RESULTS.md) | Results of the first new-data runs (GSE52778 and PBMC 3k). |
| [`runs/new-data/`](runs/new-data) | Prep scripts, the MCP stdio client (`scripts/mcp_call.py`), call plans, pre-declared expectations, results and tool outputs for those runs. |
| [`runs/exemplars/`](runs/exemplars) | The same for the NDE exemplar runs, plus `STATUS.md`. |
| [`nde-search/`](nde-search) | NDE queries, search and source-check scripts, the full hit lists, the spot-check samples and verdicts, and the NDE and bio.tools records for both tools. |
| [`min-samples/`](min-samples) | Tests of how few samples or cells each tool can handle. |
| [`skills/nde-query/`](skills/nde-query) | A Claude Code skill for querying the NDE Discovery API. |

## Not included

- **Python environments:** rebuild each from `mcp-servers/*/src/requirements.txt`; the pinned closures are in `mcp-servers/*/reports/environment-requirements.txt`.
- **Upstream code:** Paper2MCP cloned [scverse/decoupler](https://github.com/scverse/decoupler) at v2.3.1 (commit 27a1be7) and [Teichlab/celltypist](https://github.com/Teichlab/celltypist) at 1.7.1 (commit fe35756). The servers install these as dependencies rather than copying them.
- **Data:**
  - downloaded datasets, intermediate `.h5ad` files and tutorial data;
  - the copies of the CollecTRI, PROGENy and Hallmark networks saved with each run;
  - the verifiers' full result folders (about 580 MB).

  Data sources are public. Their accessions or URLs are in the plans, `STATUS.md` and the NDE results.
- **Paths:** reports, plans and logs keep absolute paths from the machine where the work was done. Scripts take their paths as arguments.

## License

[MIT](LICENSE), for the code and documentation written for this pilot. Three exceptions keep their own terms:

- **Tutorial-derived code:** each server's tools adapt that tool's own tutorials, and the executed tutorial notebooks are upstream material. Their notices are in `mcp-servers/<tool>/THIRD_PARTY/` and inside each ZIP:
  - decoupler: BSD 3-Clause, © scverse;
  - CellTypist: MIT, © Teichmann Lab.
- **Records:** those retrieved from NDE and bio.tools remain under those services' terms.
- **Run outputs:** these derive from public datasets (GEO, 10x Genomics), whose terms apply.
