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

## Running a server

Unzip `mcp-servers/<tool>/dist/<tool>-mcp.zip`. The `USAGE.md` inside covers installation, startup, every tool's inputs and outputs, and the validation limits. [`runs/new-data/scripts/mcp_call.py`](runs/new-data/scripts/mcp_call.py) calls a server over stdio from a JSON plan, the way an MCP client would.
