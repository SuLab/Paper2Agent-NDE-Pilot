# Environment manager results, attempt 2: install python-igraph

**Status: succeeded.** The coordinator asked for this dependency so that `dc.pl.network` (bulk/rna.ipynb, cell 36) can draw its plot. decoupler lists igraph in its own `test` and `full` dependency groups. This file adds to `reports/environment-manager_results.md` and does not replace it; that file was left unchanged.

## Environment

- PROJECT_PYTHON: `/home/asu/Science/paper2agent-runs/decoupler/decoupler-env/bin/python`
- Python: 3.13.12 (unchanged)
- Kernel `decoupler-p2a`: unchanged. argv[0] is still PROJECT_PYTHON.

## What changed

The resolver accepted igraph with the current stack and installed two new packages. The environment went from 162 to 164 packages, and nothing that was already installed was upgraded or downgraded.

| Package | Version |
|---|---|
| igraph (python-igraph) | 1.0.0 |
| texttable | 1.7.2 (dependency of igraph) |

These versions were unchanged: decoupler 2.3.1 (still `file:///home/asu/Science/paper2agent-runs/decoupler/repo/decoupler`), fastmcp 4.0.3, pandas 3.1.0, scanpy 1.12.4 and pydeseq2 0.5.4.

## Commands (run from PROJECT_ROOT with the explicit interpreter) and exit statuses

| # | Command | Exit |
|---|---|---|
| 1 | `uv pip install --python $PROJECT_PYTHON igraph`: installed igraph 1.0.0 and texttable 1.7.2 | 0 |
| 2 | `uv pip check --python $PROJECT_PYTHON`: "All installed packages are compatible" (164 packages) | 0 |
| 3 | `MPLBACKEND=Agg $PROJECT_PYTHON -I smoke.py` (details below) | 0 |
| 4 | `uv pip freeze --python $PROJECT_PYTHON > reports/environment-requirements.txt`: regenerated, 164 lines, includes `igraph==1.0.0` and `texttable==1.7.2` | 0 |
| 5 | `git -C repo/decoupler status --short --ignored`: empty, so `repo/` is untouched | 0 |

## What the smoke test (command 3) checked

- **Import:** `import igraph` worked, giving version 1.0.0 from `decoupler-env/lib/python3.13/site-packages/igraph/__init__.py`. `importlib.metadata.version('igraph')` also returned 1.0.0.
- **Input data:** the test used `dc.ds.toy(nobs=40, nvar=20, seed=42)`. Because pandas 3 makes `X` read-only, it ran `adata.X = adata.X.copy()` first; `X.flags.writeable` was then `True`. It built a one-row `data` DataFrame of per-gene means and ran `dc.mt.ulm(data, net, tmin=3)`, which gave acts of shape (1, 5).
- **The plot call:** `dc.pl.network(net=net, data=data, score=acts, sources=['T1','T2'], targets=5, figsize=(5,5), vcenter=True, by_abs=True, size_node=10, t_label='Gene-level\nstatistic', return_fig=True)`. These arguments match the tutorial's cell 36 settings, with toy sources in place of the tutorial's.
- **Result:** the call returned a matplotlib Figure with 4 axes. Saved as PNG it was 41,358 bytes. I looked at the image: it shows two source nodes (T1, T2), their target genes, positive and negative edge colours, and the enrichment-score and gene-level-statistic colourbars.

## Outcome and unresolved requirements

- Nothing is unresolved. igraph's matplotlib backend works with the Agg backend, so no cairo library was needed for `dc.pl.network`.
- This was a plotting smoke test only. The real cell 36, using CollecTRI and the tutorial's results, belongs to stage 2.
