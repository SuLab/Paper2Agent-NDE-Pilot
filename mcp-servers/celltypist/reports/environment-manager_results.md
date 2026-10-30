# Environment manager results: celltypist (stage 1, Python route)

Status: **READY**. One dependency constraint was required: `scikit-learn==1.7.2`. The reason is below.

## Interpreter
- PROJECT_ENV: `/home/asu/Science/paper2agent-runs/celltypist/celltypist-env`
- PROJECT_PYTHON: `/home/asu/Science/paper2agent-runs/celltypist/celltypist-env/bin/python` (`Python 3.12.13`; uv-managed CPython at `/home/asu/.local/share/uv/python/cpython-3.12.13-linux-x86_64-gnu/bin/python3.12`)
- Why 3.12: celltypist declares `python_requires>=3.6`. The full stack, including fastmcp 4.0.3, resolved cleanly for 3.12 with `uv pip compile --python-version 3.12` (exit 0).
- Snapshot: `reports/environment-requirements.txt` (164 packages, `uv pip freeze`). Python version is in `reports/python-version.txt`.

## Key package versions (resolved and installed)
| package | version |
|---|---|
| celltypist | 1.7.1 (local checkout, see below) |
| scanpy | 1.12.4 |
| anndata | 0.13.4 |
| scikit-learn | **1.7.2** (constrained; see unresolved/notes) |
| numpy | 2.5.4 |
| pandas | 3.1.0 |
| scipy | 1.18.1 |
| matplotlib | 3.11.2 |
| leidenalg | 0.12.0 |
| igraph | 1.0.0 |
| umap-learn | 0.5.12 / numba 0.68.1 |
| openpyxl / click / requests | 3.1.5 / 8.5.0 / 2.34.2 |
| fastmcp | 4.0.3 (fastmcp-slim 4.0.3, mcp 2.3.0, pydantic 2.14.1) |
| pytest / pytest-asyncio | 9.2.0 / 1.4.1 |
| papermill / nbclient / ipykernel / jupytext | 2.7.0 / 0.11.0 / 7.4.0 / 1.19.6 |

## Source identity evidence
- Checkout: `repo/celltypist` at `git rev-parse HEAD` = `fe357564a6625d3b1732a022fd39f18e55696e80`
- `celltypist.__version__` = `1.7.1`
- `celltypist.__file__` = `/home/asu/Science/paper2agent-runs/celltypist/celltypist-env/lib/python3.12/site-packages/celltypist/__init__.py` (non-editable install)
- `direct_url.json` = `{"url":"file:///home/asu/Science/paper2agent-runs/celltypist/repo/celltypist","dir_info":{}}`, INSTALLER = `uv`
- `uv pip freeze` line: `celltypist @ file:///home/asu/Science/paper2agent-runs/celltypist/repo/celltypist`
- `diff -rq --exclude=__pycache__ repo/celltypist/celltypist <site-packages>/celltypist` → no differences (exit 0). The installed package is byte-identical to the pinned checkout.
- Repo cleanliness: the setuptools build wrote the git-ignored `build/` and `celltypist.egg-info/` into `repo/celltypist` at 15:45:40, three minutes after the 15:42 clone. I removed both, and afterwards `git status --short --ignored` is empty. No tracked file was modified.

## Commands and exit status
(cwd = PROJECT_ROOT; `uv` = `~/.local/bin/uv`)
| # | command | exit |
|---|---|---|
| 1 | `uv pip compile <reqs: ./repo/celltypist, fastmcp==4.0.3, pytest, pytest-asyncio, papermill, nbclient, ipykernel, jupytext> --python-version 3.12` (dry run) | 0 (picked scikit-learn 1.9.1) |
| 2 | `uv venv celltypist-env --python cpython-3.12.13 --managed-python` | 0 |
| 3 | `uv pip install --python celltypist-env/bin/python ./repo/celltypist` | 0 |
| 4 | `rm -rf repo/celltypist/build repo/celltypist/celltypist.egg-info` (cleanup of my own build artifacts) | 0 |
| 5 | `uv pip install --python celltypist-env/bin/python 'fastmcp==4.0.3' pytest pytest-asyncio papermill nbclient ipykernel jupytext` | 0 |
| 6 | `uv pip check` (before the sklearn change) | 0 |
| 7 | `LogisticRegression(multi_class='ovr').fit(...)` under scikit-learn 1.9.1 | **1**: `TypeError: unexpected keyword argument 'multi_class'` |
| 8 | Same probe in scratch venvs: sklearn 1.8.0 → TypeError; sklearn 1.7.2 → fit OK, with a FutureWarning ("deprecated in 1.5, removed in 1.8") | 1 / 0 |
| 9 | `uv pip install --python celltypist-env/bin/python 'scikit-learn==1.7.2'` | 0 |
| 10 | `uv pip check --python celltypist-env/bin/python` | 0 ("All installed packages are compatible", 164 pkgs) |
| 11 | Import check: celltypist, scanpy, anndata, leidenalg, igraph, sklearn, numpy, pandas, scipy, matplotlib, fastmcp, mcp, pydantic, papermill, nbclient, ipykernel, jupytext, pytest, pytest_asyncio, openpyxl, click, requests, umap, numba | 0 |
| 12 | Synthetic API smoke test (scratch script, random data, no downloads; listed below) | 0 |
| 13 | `celltypist-env/bin/celltypist --help` (shebang = PROJECT_PYTHON) | 0 |
| 14 | `GET https://celltypist.cog.sanger.ac.uk/models/models.json` (requests, timeout 30) | 200, 25,568 bytes, 61 models; Immune_All_Low.pkl and Immune_All_High.pkl listed |
| 15 | `curl -sI` HEAD for the tutorial data URLs (no body download) | all 200 (see network section) |
| 16 | `python -m ipykernel install --prefix celltypist-env --name celltypist-p2a --display-name "Python (celltypist-p2a)"` | 0 |
| 17 | `python -m jupyter kernelspec list` | 0 |
| 18 | `python -m papermill kernel_probe.ipynb kernel_probe_out.ipynb --kernel celltypist-p2a` (scratch probe notebook, CELLTYPIST_FOLDER set) | 0 |
| 19 | FastMCP probe: `FastMCP` + `@mcp.tool()` + `Client(server)` + `call_tool` → `result.data == 5` | 0 |
| 20 | Wrote `pytest.ini`; `python -m pytest --co -q` | 5 (expected: no tests exist yet in tests/code; config parsed with no warnings) |
| 21 | `uv pip freeze > reports/environment-requirements.txt`; `python --version > reports/python-version.txt` | 0 / 0 |

### Synthetic API smoke test (#12)
This checks API compatibility only. It is not scientific validation. Setup: 1600 cells × 1000 genes of random Poisson counts, 4 planted types, `normalize_total(1e4)` + `log1p`, with `CELLTYPIST_FOLDER` set. Everything ran without error:
- `celltypist.train(adata, 'cell_type', feature_selection=True, top_genes=100, n_jobs=10)` ran the SGD pass, then the second LR pass with `multi_class='ovr'`, in 12.5 s and kept 330 features.
- `model.extract_top_markers('T0', 5)`, `model.write(path)`, and `models.Model.load(path)` all worked.
- `celltypist.annotate(adata, model=..., majority_voting=True)` built the neighbor graph, ran Leiden over-clustering (`flavor='igraph'` path, since scanpy ≥ 1.10), and ran majority voting.
- `predictions.to_adata(insert_prob=True)` worked.
- `celltypist.dotplot`, `sc.tl.umap`, `sc.pl.umap`, `sc.pl.violin`, and `sc.pl.dotplot` all saved PNGs (Agg backend).

## Kernel
- Name: `celltypist-p2a`; display name `Python (celltypist-p2a)`
- Path: `/home/asu/Science/paper2agent-runs/celltypist/celltypist-env/share/jupyter/kernels/celltypist-p2a/kernel.json`
- argv: `["/home/asu/Science/paper2agent-runs/celltypist/celltypist-env/bin/python", "-Xfrozen_modules=off", "-m", "ipykernel_launcher", "-f", "{connection_file}"]`
- `KernelSpecManager().get_kernel_spec('celltypist-p2a')` resolves to the path above, and `argv[0] == PROJECT_PYTHON` is True. No other `celltypist-p2a` kernelspec exists (`~/.local/share/jupyter/kernels` does not exist).
- Papermill probe: inside the kernel, `sys.executable` = PROJECT_PYTHON, celltypist 1.7.1 imports from the env, and `models_path` honours CELLTYPIST_FOLDER.
- Papermill must be run as `"$PROJECT_PYTHON" -m papermill ... --kernel celltypist-p2a` so that it discovers the env-prefix kernelspec.

## Device check
- There is no usable GPU: `nvidia-smi` is not found. `/dev/dxg` (the WSL paravirtual device) exists, but there are no CUDA tools or libraries in the env.
- `cuml`, `rapids_singlecell`, and `torch` are not installed (`find_spec` → None).
- celltypist's `use_GPU` defaults to False. If it were set True, celltypist would log a warning and fall back to CPU for over-clustering. Training with `use_GPU=True` checks for `cuml` in `sys.modules`. **Executors must keep `use_GPU=False`.**
- CPU: `os.cpu_count()` = 20, about 31 GB RAM. The tutorial's `n_jobs=10` fits.

## Network check
- `https://celltypist.cog.sanger.ac.uk/models/models.json`: HTTP 200, 61 models. Its URL for Immune_All_Low is `https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl`, and a HEAD request returned 200 (2,824,990 bytes).
- Tutorial `backup_url` data (HEAD only, nothing downloaded):
  - `Notebook_demo_data/demo_2000_cells.h5ad`: 200, 35.7 MB
  - `demo_400_cells.h5ad`: 200, 8.0 MB
  - `demo_500_cells.h5ad`: 200, 9.9 MB
- I downloaded no models or datasets.

## CELLTYPIST_FOLDER handling (required runtime variable)
- **Required:** `export CELLTYPIST_FOLDER=/home/asu/Science/paper2agent-runs/celltypist/tmp/celltypist_home` for every executor, papermill, pytest, and MCP-server process. celltypist reads it once at `import celltypist.models` (`os.getenv('CELLTYPIST_FOLDER', ~/.celltypist)`). Set it before Python starts, not after import.
- Verified with the variable set: `models_path` = `/home/asu/Science/paper2agent-runs/celltypist/tmp/celltypist_home/data/models`, and `celltypist_path` = `.../tmp/celltypist_home`. The same holds inside the papermill kernel.
- On import, celltypist creates `tmp/celltypist_home/data/models/`. That directory currently holds no models.
- I wrote nothing to `~/.celltypist`. It pre-existed (dated Mar 1, 2026), and `find ~/.celltypist -newermt "2026-10-27 15:40"` returns nothing.
- Optional: matplotlib writes its font cache to `~/.cache/matplotlib` unless `MPLCONFIGDIR` is set. I set it to scratch for my smoke run only.

## pytest
`/home/asu/Science/paper2agent-runs/celltypist/pytest.ini`:
```
[pytest]
asyncio_mode = auto
asyncio_default_fixture_loop_scope = function
testpaths = tests/code
```
Run tests with `"$PROJECT_PYTHON" -m pytest tests/code/`.

## Unresolved requirements and notes
1. **scikit-learn is constrained to 1.7.2.** This is required, not cosmetic. `celltypist/train.py` (`_LRClassifier`, line 126, and `_cuLRClassifier`) calls `LogisticRegression(..., multi_class='ovr')`. sklearn 1.8.0 and 1.9.1 (the unconstrained resolver pick) reject that keyword with `TypeError`, which would break `celltypist.train(..., feature_selection=True)` in the ML tutorial. 1.7.2 is the newest version I verified that accepts it; it emits a FutureWarning that celltypist's own filters suppress. I did not modify the scientific source. **Stage 5 `src/requirements.txt` must pin `scikit-learn==1.7.2`, or at least `<1.8`.**
2. **Real-model unpickling is not yet verified.** Pretrained models were pickled with older sklearn and numpy, so loading them under sklearn 1.7.2 / numpy 2.5.4 may emit `InconsistentVersionWarning`. Prediction uses `decision_function` + `expit`, not `predict_proba`. Executors will be the first to load a real downloaded model and should confirm there are no errors.
3. **pandas 3.1.0 / numpy 2.5.4 / scanpy 1.12.4 are much newer than celltypist 1.7.1's release (June 2025).** The synthetic smoke covered every API path the tutorials use without error, but numerical parity with the published notebook outputs is for stage 2 and 3 to establish.
4. **Harmless warnings.** Every import of `celltypist` emits a FutureWarning, ``scanpy `__version__` is deprecated`` (from `classifier.py:11`). It still works under scanpy 1.12.4. With `pytest -W error` it would become an error, so do not escalate FutureWarnings globally.
5. **No blockers:** no credentials, native libraries, or R are needed. The network is reachable. There is no GPU, and none is needed (CPU route).
