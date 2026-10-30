# CellTypist MCP server

This MCP server exposes CellTypist's automated cell-type annotation as five tools. Each tool calls the [Teichlab/celltypist](https://github.com/Teichlab/celltypist) Python package, pinned to commit `fe357564a6625d3b1732a022fd39f18e55696e80` (celltypist 1.7.1). None of the classification, training, majority voting or plotting is reimplemented.

Paper: Domínguez Conde *et al.*, "Cross-tissue immune cell analysis reveals tissue-specific features in humans", *Science* 376, eabl5197 (2022). <https://doi.org/10.1126/science.abl5197>

## What it covers

The tools come from two CellTypist tutorials:
- `docs/notebook/celltypist_tutorial.ipynb`: annotation with built-in or custom models, majority voting, model training and marker inspection.
- `docs/notebook/celltypist_tutorial_ml.ipynb`: multi-label "prob match" annotation and comparison with existing labels.

Not included:
- **The third tutorial** (`celltypist_tutorial_cv.ipynb`, cross-dataset label transfer). It needs a multi-GB gut-atlas download.
- **CellTypist functions the two tutorials never use**: `convert`, `subset`, `to_plots`, `downsample_adata` and the command-line interface.
- **UMAP plotting.** It is generic scanpy, not CellTypist.
- **`dotplot`'s `filter_prediction` option.** With the pinned scanpy 1.12.4 and pandas 3.1.0, upstream `celltypist.dotplot` raises a `KeyError` whenever a non-zero threshold actually removes a label. This was reproduced with direct upstream calls.

### Input expectations

CellTypist requires an AnnData `.h5ad` with these properties:
- **Expression:** log1p-normalised to 10,000 counts per cell, in `.X`, or in `.raw.X` when `.X` is scaled.
- **Genes:** human gene symbols as `var_names`.

CellTypist checks this itself and the tools report its error. Data stored as raw counts, or as log1p of un-normalised counts, must be normalised first: `scanpy.pp.normalize_total(adata, target_sum=1e4)` followed by `scanpy.pp.log1p(adata)`.

## Install

You need Linux x86_64 (the tested platform), [uv](https://docs.astral.sh/uv/), git, and network access to GitHub and PyPI.

From the extracted `celltypist-mcp/` folder:

```bash
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python -r src/requirements.txt
uv pip check --python .venv/bin/python
```

`src/requirements.txt` pins the full tested dependency closure:
- CPython 3.12.13, fastmcp 4.0.3, scanpy 1.12.4, anndata 0.13.4, pandas 3.1.0, numpy 2.5.4.
- leidenalg 0.12.0 and igraph 1.0.0.
- scikit-learn 1.7.2. Keep this pin: `celltypist.train` passes `LogisticRegression(multi_class='ovr')`, which scikit-learn ≥ 1.8 rejects.
- CellTypist itself, installed from the pinned git commit.

## Models and the model cache

Built-in models are downloaded on first use from `https://celltypist.cog.sanger.ac.uk`. They are cached under `$CELLTYPIST_FOLDER/data/models`, which defaults to `~/.celltypist`. Set `CELLTYPIST_FOLDER` **before** the server starts, because celltypist reads it at import time.

When a tool is given a built-in model name, it downloads only that one model. A bare `celltypist.models.Model.load(name)` on an empty cache would instead download all 61 models. Custom models are passed as a path to a `.pkl` file.

## Start the server (stdio)

```bash
CELLTYPIST_FOLDER="$PWD/celltypist_home" .venv/bin/python src/celltypist_mcp.py
```

**Client configuration.** Use the absolute paths of your extracted folder.

```json
{
  "mcpServers": {
    "celltypist": {
      "command": "/ABSOLUTE/PATH/celltypist-mcp/.venv/bin/python",
      "args": ["/ABSOLUTE/PATH/celltypist-mcp/src/celltypist_mcp.py"],
      "env": {"CELLTYPIST_FOLDER": "/ABSOLUTE/PATH/celltypist-mcp/celltypist_home"}
    }
  }
}
```

**Claude Code.** The syntax below is taken from `fastmcp install claude-code --help` (fastmcp 4.0.3). This command was not run during validation.

```bash
.venv/bin/fastmcp install claude-code src/celltypist_mcp.py --name celltypist \
  --python 3.12 --with-requirements src/requirements.txt \
  --env CELLTYPIST_FOLDER="$PWD/celltypist_home"
```

## Tools

Every tool creates a fresh subdirectory for each call, under `output_dir` if you pass one, otherwise under `<system temp>/celltypist_mcp_outputs`. It returns a `message`, the source `reference`, and `artifacts` as absolute paths, plus compact summaries.

### `celltypist_list_models`
Lists the built-in models and what each one represents.
- `refresh_index` (default `false`): re-download `models.json` first.
- Output: a CSV of model file names and descriptions. The result also gives the default model.

### `celltypist_annotate_cells`
Labels each cell with a built-in or custom model, using `celltypist.annotate`.
- **Required:**
  - `adata_path`: a normalised `.h5ad`.
- **Options:**
  - `model`: default `Immune_All_Low.pkl`, or a `.pkl` path.
  - `mode`: `best match` (default) or `prob match`.
  - `p_thres`: probability threshold for prob match, default 0.5.
  - `majority_voting`: default `false`.
  - `over_clustering_key`: an existing obs column to use instead of CellTypist's own Leiden over-clustering.
  - `min_prop`: default 0.
- **Outputs:** written into the returned `output_dir`.
  - `predicted_labels.csv`, `decision_matrix.csv` and `probability_matrix.csv`, exactly as `AnnotationResult.to_table` writes them.
  - `annotated.h5ad`: the input plus label and confidence columns. Probabilities stay in CSV because some cell-type names contain `/`.
- The result also includes a label-frequency summary.

### `celltypist_train_model`
Trains a custom logistic-regression model from a labelled reference, using `celltypist.train` and `Model.write`.
- **Required:**
  - `adata_path`.
  - `label_key`: the obs column of labels.
- **Options:**
  - `feature_selection`: default `false`.
  - `top_genes`: default 300.
  - `n_jobs`.
  - `seed`: passed as `random_state` to the SGD and logistic-regression classifiers.
  - `model_name`.
- **Output:** a `.pkl` model usable by the annotate and inspect tools.
- **Reproducibility:** without `seed`, feature selection is not reproducible. With `seed`, results are bit-identical at the same `n_jobs`.

### `celltypist_inspect_model`
Shows a model's metadata and the top genes driving each cell type, using `Model.extract_top_markers`.
- **Required:**
  - `model`.
- **Options:**
  - `cell_types`: default all.
  - `top_n`: default 10.
  - `only_positive`: default `true`.
  - `adata_path` and `groupby_key`: when both are given, the tool also draws the tutorial's `scanpy.pl.violin` of those genes for each group.
- **Output:** a markers CSV and, optionally, violin PNGs.

### `celltypist_compare_with_reference`
Compares CellTypist predictions with labels or clusters already present in the data.
- **Required:**
  - `annotation_dir`: the `output_dir` returned by `celltypist_annotate_cells`.
  - `reference_key`: an obs column of `annotated.h5ad`.
- **Options:**
  - `use_as_prediction`: `majority_voting` (default) or `predicted_labels`.
  - `reference_groups`: groups to report and plot.
  - `plot_group_probabilities`: default `true`.
- **Outputs:**
  - the tutorial's `celltypist.dotplot` PNG, with its size and colour tables as CSV;
  - the full reference × prediction crosstab as CSV;
  - per-group probability boxplots.
- **Composition:** run `celltypist_annotate_cells` first, then pass its `output_dir` here as `annotation_dir`.

## Validation performed

Tool correctness was checked against independent upstream evidence: the real tutorial notebooks, executed end to end, and direct `celltypist` calls on the same inputs. Separate verifier agents wrote the tests; they did not write the wrappers.

**Annotation**
- On the tutorial's 2,000-cell demo with `Immune_All_Low` (98 cell types), labels, over-clustering and majority-voting labels match exactly.
- Decision scores match within 1e-9 and probabilities within 2.5e-10.
- 5,645 features are used, as in the published tutorial.

**Training**
- Seeded training reproduces the reference's 2,652 genes and identical coefficients.

**Comparison**
- The multi-label tutorial's crosstab rows, and 13 full crosstab/dot-plot tables, match direct `celltypist.dotplot` calls exactly.

**Changed inputs and errors**
- Tests covered a renamed label column, a custom model, prob match at a different threshold, a user-supplied over-clustering, different groups, and both prediction columns.
- Each was compared with a new direct upstream call.
- Missing files, columns and groups, unknown models, non-normalised input and misaligned annotation files all return MCP tool errors.

**Packaged server over real stdio**
- 33 acceptance cases (18 successful-call cases and 15 error cases) pass in the development environment.
- They also pass in a fresh environment installed only from `src/requirements.txt`.

### Limits

- **Majority voting:** the over-clustering, and therefore the majority-voting labels, are exact only for the pinned scanpy, leidenalg and igraph versions. Their cluster IDs differ from the tutorial's stored scanpy 1.7.1 output.
- **Multi-label tutorial and model builds:** the multi-label tutorial's stored outputs came from an older 90-type `Immune_All_Low` build. Under the current 98-type build:
  - "Microglia → Unassigned" still holds.
  - The built-in model's "Macrophages|pDC" call does not reproduce. The new model splits macrophage subtypes.
  - Expect similar shifts whenever a model build changes.
- **Not validated:**
  - **Model downloads:** they depend on the Sanger server. One validation download timed out at celltypist's fixed 30 s limit and succeeded on retry.
  - **GPU paths:** not exposed and not tested.
  - **Platforms other than Linux x86_64:** not tested.
- **Scope:** these checks establish faithful wrappers for the selected tools only. They do not re-validate CellTypist's models or the paper's results.
