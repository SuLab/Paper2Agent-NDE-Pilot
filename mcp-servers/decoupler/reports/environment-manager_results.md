# Environment manager results: decoupler (stage 1, Python route)

**Status: succeeded.** I found no blockers. There is one behavioural caveat for downstream workers: pandas 3 makes `AnnData(DataFrame).X` read-only. It is described below.

## Interpreter

| Item | Value |
|---|---|
| PROJECT_ENV | `/home/asu/Science/paper2agent-runs/decoupler/decoupler-env` |
| PROJECT_PYTHON | `/home/asu/Science/paper2agent-runs/decoupler/decoupler-env/bin/python` |
| Python | CPython 3.13.12 (uv-managed: `~/.local/share/uv/python/cpython-3.13.12-linux-x86_64-gnu`), not the system Python |
| Why 3.13 | decoupler declares `requires-python >=3.12` and classifiers for 3.12 to 3.15. CPython 3.13 was already installed and supported by numba 0.68.1, pydeseq2, scanpy and fastmcp 4.0.3. Everything resolved together in one pass. |
| venv tool | uv 0.10.11 (`uv venv --seed`, so pip 26.2.2 is also present) |

## Key package versions (resolved together in one `uv pip install`)

decoupler 2.3.1 (local checkout) · scanpy 1.12.4 · pydeseq2 0.5.4 · anndata 0.13.4 · numba 0.68.1 / llvmlite 0.50.0 · numpy 2.5.4 · scipy 1.18.1 · pandas 3.1.0 · scikit-learn 1.9.1 · marsilea 0.8.1 · adjusttext 1.4.0 · matplotlib 3.11.2 · requests 2.34.2 · **fastmcp 4.0.3** (+ fastmcp-slim 4.0.3, mcp 2.3.0, pydantic 2.14.1) · pytest 9.2.0 · pytest-asyncio 1.4.1 · papermill 2.7.0 · nbclient 0.11.0 · ipykernel 7.4.0 · jupyter-client 8.10.0 · jupytext 1.19.6.

FastMCP 4.0.3 did **not** conflict with the scientific stack. It resolved in the same pass, and `uv pip check` reports no conflicts. The full snapshot of 162 packages is in `reports/environment-requirements.txt`.

## decoupler source identity evidence

- Checkout: `repo/decoupler` at HEAD `27a1be7c75d88306bf3c2e4d0e7fb4b42d48f78a`. The `docs/notebooks` submodule is at `306be6673a5df5baf780e44ea929391616dbe5d8`. Both were re-checked after the install.
- `importlib.metadata.version('decoupler')` returns `2.3.1`, and `decoupler.__version__` also returns `2.3.1`.
- `decoupler.__file__` is `/home/asu/Science/paper2agent-runs/decoupler/decoupler-env/lib/python3.13/site-packages/decoupler/__init__.py`. This is a non-editable install.
- `direct_url.json` contains `{"url":"file:///home/asu/Science/paper2agent-runs/decoupler/repo/decoupler","dir_info":{}}`, and INSTALLER is `uv`. The freeze line reads `decoupler @ file:///home/asu/Science/paper2agent-runs/decoupler/repo/decoupler`.
- `diff -rq --exclude=__pycache__ repo/decoupler/src/decoupler <site-packages>/decoupler` exited 0. The installed files are byte-identical to the checkout.
- `git -C repo/decoupler status --short --ignored` is empty after the install. Nothing under `repo/` was modified or left behind as a build artifact.

## Commands run (from PROJECT_ROOT, explicit interpreter) and exit statuses

| # | Command | Exit |
|---|---|---|
| 1 | `uv venv decoupler-env --python ~/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 --seed` | 0 |
| 2 | `uv pip install --python $PROJECT_PYTHON ./repo/decoupler scanpy pydeseq2 'fastmcp==4.0.3' pytest pytest-asyncio papermill nbclient ipykernel jupytext` | 0 |
| 3 | `uv pip check --python $PROJECT_PYTHON`: "All installed packages are compatible" (162 packages) | 0 |
| 4 | `$PROJECT_PYTHON -I envcheck.py`: identity checks plus imports of decoupler, scanpy, pydeseq2 (`DeseqDataSet`, `DefaultInference`, `DeseqStats`), fastmcp (`FastMCP`, `Client`), anndata, numba, numpy, scipy, pandas, marsilea, adjustText, matplotlib, pytest, pytest_asyncio, papermill, nbclient, ipykernel, jupytext, pydantic, mcp and requests. All imports succeeded. All decoupler submodules are present: pp, mt, op, ds, pl, tl, bm. | 0 |
| 5 | `diff -rq` of the installed package against `repo/decoupler/src/decoupler` | 0 |
| 6 | Network HEAD checks (see below) with `curl -I` and Python `requests.head` | 0 |
| 7 | `$PROJECT_PYTHON -I -m ipykernel install --prefix $PROJECT_ENV --name decoupler-p2a --display-name "Python (decoupler-p2a)"` | 0 |
| 8 | `decoupler-env/bin/jupyter kernelspec list` | 0 |
| 9 | A live kernel test via `jupyter_client.start_new_kernel(kernel_name="decoupler-p2a")`. In the kernel, `sys.executable` was PROJECT_PYTHON, decoupler 2.3.1 imported from site-packages, and `dc.mt.ulm` on `dc.ds.toy()` worked. The scanpy step failed on the toy data with a read-only array error (see the caveat). | 0 (script); cell raised ValueError |
| 10 | An offline simulation of the tutorial chain, cells 4 to 32. It used a synthetic counts DataFrame built the same way as `dc.ds.hsctgfb()`: `AnnData(df)` then `X.astype(float)`. It ran `dc.pp.filter_by_expr`, `sc.pp.normalize_total/log1p/scale`, `sc.tl.pca`, `dc.pp.swap_layer`, `dc.tl.rankby_obsm`, pydeseq2 `DeseqDataSet(...).deseq2()` + `DeseqStats(...).summary()` (with `DefaultInference(n_cpus=8)`), and `dc.mt.ulm` on a DataFrame. All steps succeeded. | 0 |
| 11 | A plotting smoke test with `MPLBACKEND=Agg`: `dc.pl.obsm` (marsilea), `dc.pl.obsbar`, `dc.pl.volcano` (adjustText) and `dc.pl.barplot` each returned a Figure. | 0 |
| 12 | The `fastmcp version` CLI reported FastMCP 4.0.3 and MCP 2.3.0. | 0 |
| 13 | A FastMCP API smoke test: a `FastMCP` server with one `@mcp.tool()`, then `Client(server)`, `list_tools`, and `call_tool('add', {a:2,b:3})`. `result.data` was `{'sum': 5}`. | 0 |
| 14 | `papermill --version` (2.7.0), `jupytext --version` (1.19.6), `pytest --version` (9.2.0) | 0 / 0 / 0 |
| 15 | `$PROJECT_PYTHON -m pytest --collect-only`. pytest read `configfile: pytest.ini` and `testpaths: tests/code`. The run exited 5 because `tests/code/` is still empty, which is expected at this stage. | 5 (no tests yet) |
| 16 | `uv pip freeze --python $PROJECT_PYTHON > reports/environment-requirements.txt` | 0 |
| 17 | `$PROJECT_PYTHON --version > reports/python-version.txt` (`Python 3.13.12`) | 0 |
| 18 | `uv pip check` again after the snapshot | 0 |

The helper scripts for steps 4, 9, 10, 11 and 13 were kept in the session scratchpad, outside the project.

## Kernel

- Name: `decoupler-p2a`, display name "Python (decoupler-p2a)"
- Path: `/home/asu/Science/paper2agent-runs/decoupler/decoupler-env/share/jupyter/kernels/decoupler-p2a/kernel.json`
- argv: `["/home/asu/Science/paper2agent-runs/decoupler/decoupler-env/bin/python", "-Xfrozen_modules=off", "-m", "ipykernel_launcher", "-f", "{connection_file}"]`. argv[0] is PROJECT_PYTHON, which I verified.
- `jupyter kernelspec list`, run from the env, shows no other `decoupler-p2a` (for example in the user-level `~/.local/share/jupyter/kernels`), so nothing shadows it. A live kernel launch confirmed that `sys.executable` is PROJECT_PYTHON.
- Use `"$PROJECT_PYTHON" -m papermill IN OUT --kernel decoupler-p2a`.

## Device check

- CPU only. `os.cpu_count()` is 20 and `numba.config.NUMBA_NUM_THREADS` is 20.
- There is no GPU: `nvidia-smi` is not on PATH, `/dev/nvidia*` does not exist, and there is no WSL `libcuda`. No torch or cupy is installed.
- The tutorial needs no GPU. pydeseq2 uses joblib CPU workers (`DefaultInference(n_cpus=8)` fits the 20-core host).

## Network check (HEAD requests only, no resource bodies downloaded)

| URL (used by the tutorial) | Result |
|---|---|
| `https://omnipathdb.org/` (OmniPath base; `dc.op.progeny` / resource queries) | 200 text/html. Python `requests.head` from PROJECT_PYTHON also returned 200. |
| `https://static.omnipathdb.org/tables/msigdb-hallmark.tsv.gz` (`dc.op.hallmark`) | 200 application/gzip, 90,854 B |
| `https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1` (`dc.op.collectri`) | 200 text/plain, 4,345,649 B |
| `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE151251&...counts.tsv.gz` (`dc.ds.hsctgfb`) | 200 application/x-gzip, 1,578,642 B |

These checks show that the hosts were reachable at the time of the check. They do not show that the full downloads will succeed later during execution.

## Caveat for downstream workers (not a blocker for rna.ipynb)

- **What I observed directly:** with pandas 3.1.0, `pd.DataFrame(...).values.flags.writeable` is `False`, and `anndata.AnnData(df).X.flags.writeable` is `False`.
  - `dc.ds.toy()` builds its AnnData this way and does not copy, so its `X` is read-only.
  - On that data, in-place scanpy steps fail: `sc.pp.log1p` stops at `np.log1p(x, out=x)` and `sc.pp.scale` stops at `x -= mean`, both with `ValueError: output array is read-only`.
- **Not affected:** `dc.ds.hsctgfb()`, which the selected tutorial uses, calls `adata.X = adata.X.astype(float)`. That makes a writeable copy, and my offline simulation of the tutorial chain passed.
- **Hypothesis, not tested:** this is pandas 3 copy-on-write behaviour, and pandas below 3 would return a writeable view. I did not build a pandas<3 environment to confirm it.
- **What I did not do:** I did not pin pandas<3. decoupler 2.3.1 declares no pandas cap, and the selected tutorial path works.
- **Advice for later work:** a tool or test fixture that builds AnnData from a DataFrame and then calls in-place scanpy preprocessing should copy `X` first, for example `adata.X = adata.X.copy()` or `np.array(..., copy=True)`.

## Notes and unresolved requirements

- **Unresolved requirements: none.** No credentials or native libraries are needed, and R is not used on this route.
- I did not execute the real tutorial or its full network downloads. That is stage 2. The checks above are environment and import smoke tests, not scientific validation.
- `pytest.ini` contains `[pytest]`, `testpaths = tests/code`, `asyncio_mode = auto` and `asyncio_default_fixture_loop_scope = function`.
- The `fastmcp version` CLI reports that 4.1.4 is available. I left the tested baseline at 4.0.3 as instructed.
- I did not install the repo's other `test` / `full` group extras (gseapy, igraph, xgboost, dcor, ipywidgets, memory-profiler, coverage, pytest-cov). The selected tutorial does not use them. statsmodels 0.15.0 is installed anyway as a transitive dependency. Later package additions requested by workers should go through this environment owner one at a time.
