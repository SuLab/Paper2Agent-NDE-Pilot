"""Independent verification of celltypist_compare_with_reference (src/tools/celltypist_tutorial_ml.py).

Expected values never come from the wrapper. Sources:
  * the executed ml notebook (cell 47 printed output, parsed from the notebook JSON);
  * executor reference outputs in notebooks/celltypist_tutorial_ml/reference/ and
    notebooks/celltypist_tutorial/reference/ (direct upstream calls, same environment);
  * direct in-memory upstream calls saved by build_fixtures.py on the same AnnotationResult that wrote
    each fixture: celltypist.dotplot(..., return_fig=True).dot_size_df/.dot_color_df, pd.crosstab of the
    to_adata() obs (ml cell 47 form) and to_adata(insert_prob=True).obs[model.cell_types] (ml cell 82 data).
All tool calls go through the FastMCP client (Client(server).call_tool -> result.data).

Tolerances: crosstab counts are integers and compared exactly. Dot-plot fractions / mean probabilities and
boxplot probabilities are compared exactly (tolerance 0): the tool reads probabilities written by pandas
to_csv (shortest round-trip float repr, lossless) and upstream performs the identical computation on them.
No nonfinite values are expected; the comparisons assert finiteness of the reference values.
"""
import json
import re
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from fastmcp import Client
from fastmcp.exceptions import ToolError

from src.tools.celltypist_tutorial_ml import celltypist_tutorial_ml_mcp

TOOL = "celltypist_compare_with_reference"
ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "tests/data/celltypist_tutorial_ml"
FIXTURES = DATA / "fixtures"
EXPECTED = DATA / "expected"
INPUTS = DATA / "inputs"
ML_REF = ROOT / "notebooks/celltypist_tutorial_ml/reference"
MAIN_REF = ROOT / "notebooks/celltypist_tutorial/reference/annotate_demo2000_ImmuneAllLow_mv"
ML_NOTEBOOK = ROOT / "notebooks/celltypist_tutorial_ml/celltypist_tutorial_ml_execution_final.ipynb"
RESULTS = ROOT / "tests/results/celltypist_tutorial_ml"
MANIFEST = json.loads((DATA / "fixtures_manifest.json").read_text())
CONTRACT_FILES = ("predicted_labels.csv", "decision_matrix.csv", "probability_matrix.csv", "annotated.h5ad")


# ----------------------------------------------------------------------------------------------- helpers
async def call(args: dict) -> dict:
    async with Client(celltypist_tutorial_ml_mcp) as client:
        result = await client.call_tool(TOOL, args)
    return result.data


async def call_error(args: dict) -> str:
    async with Client(celltypist_tutorial_ml_mcp) as client:
        with pytest.raises(ToolError) as excinfo:
            await client.call_tool(TOOL, args)
    return str(excinfo.value)


def art(result: dict, filename: str) -> Path:
    paths = [Path(a["path"]) for a in result["artifacts"] if Path(a["path"]).name == filename]
    assert len(paths) == 1, (filename, result["artifacts"])
    assert paths[0].is_absolute() and paths[0].is_file()
    return paths[0]


def read_out_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0, keep_default_na=False)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def assert_same_table(got: pd.DataFrame, exp: pd.DataFrame, integer: bool = False) -> None:
    assert list(got.index) == [str(x) for x in exp.index]
    assert list(got.columns) == [str(x) for x in exp.columns]
    e = exp.to_numpy(dtype=float)
    assert np.isfinite(e).all()
    np.testing.assert_array_equal(got.to_numpy(dtype=float), e)
    if integer:
        assert all(np.issubdtype(t, np.integer) for t in got.dtypes)


def expected(fixture: str, name: str) -> pd.DataFrame:
    return pd.read_pickle(EXPECTED / fixture / f"{name}.pkl")


def keep(result: dict, tag: str) -> None:
    """Copy a tool output directory into tests/results for visual review."""
    dest = RESULTS / tag
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(result["output_dir"], dest)


def tree_hashes(folder: Path) -> dict:
    import hashlib
    return {str(p.relative_to(folder)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(folder).rglob("*")) if p.is_file()}


def is_png(path: Path) -> bool:
    return path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n" and path.stat().st_size > 5000


def notebook_cell47() -> pd.DataFrame:
    """Parse the printed DataFrame of executed ml notebook cell 47 (wrapped pandas repr)."""
    nb = json.loads(ML_NOTEBOOK.read_text())
    cell = nb["cells"][47]
    assert "pd.crosstab(adata.obs.cell_type, adata.obs.majority_voting)" in "".join(cell["source"])
    text = "".join(next(o for o in cell["outputs"] if "data" in o)["data"]["text/plain"])
    columns, rows = [], {}
    for block in [b for b in text.split("\n\n") if b.strip()]:
        lines = [ln.rstrip().rstrip("\\").rstrip() for ln in block.splitlines() if ln.strip()]
        header = re.split(r"\s{2,}", lines[0].strip())
        assert header[0] == "majority_voting" and lines[1].strip() == "cell_type"
        columns += header[1:]
        for ln in lines[2:]:
            name, *vals = re.split(r"\s{2,}", ln.strip())
            rows.setdefault(name, []).extend(int(v) for v in vals)
    df = pd.DataFrame.from_dict(rows, orient="index", columns=columns)
    assert list(df.index) == ["Microglia", "Macro_pDC"]
    return df


# ------------------------------------------------------------------------------------- fixture anchoring
@pytest.mark.parametrize("fixture,ref_dir", [
    ("builtin_bestmatch", ML_REF / "builtin_bestmatch"),
    ("custom_seed0_probmatch", ML_REF / "custom_seed0_probmatch"),
    ("main_demo2000_bestmatch", MAIN_REF),
])
def test_fixtures_are_byte_identical_to_executor_reference(fixture, ref_dir):
    files = CONTRACT_FILES if (ref_dir / "annotated.h5ad").exists() else CONTRACT_FILES[:3]
    for f in files:
        assert (FIXTURES / fixture / f).read_bytes() == (ref_dir / f).read_bytes(), f


# ---------------------------------------------------------------------------- tutorial reference values
async def test_ml_cell47_rows_exact(out_base):
    """ml cell 47: pd.crosstab(adata.obs.cell_type, adata.obs.majority_voting).loc[['Microglia','Macro_pDC']]."""
    res = await call({"annotation_dir": str(FIXTURES / "builtin_bestmatch"), "reference_key": "cell_type",
                      "use_as_prediction": "majority_voting", "reference_groups": ["Microglia", "Macro_pDC"],
                      "output_dir": str(out_base)})
    keep(res, "builtin_bestmatch_majority_voting_cell47")
    crosstab = read_out_table(art(res, "crosstab.csv"))
    nb47 = notebook_cell47()
    got = crosstab.loc[["Microglia", "Macro_pDC"]]
    assert_same_table(got, nb47, integer=True)
    executor = read_out_table(ML_REF / "builtin_bestmatch/crosstab_cell_type_x_majority_voting_Microglia_Macro_pDC.csv")
    assert_same_table(got, executor, integer=True)
    for g in ("Microglia", "Macro_pDC"):
        row = nb47.loc[g]
        assert res["selected_groups"][g]["prediction_counts"] == {k: int(v) for k, v in row.items() if v > 0}
        assert res["selected_groups"][g]["n_cells"] == int(row.sum())
        counts = list(res["selected_groups"][g]["prediction_counts"].values())
        assert counts == sorted(counts, reverse=True)
    assert res["n_cells"] == 500 and res["n_reference_groups"] == 12 and res["n_predicted_labels"] == 11
    assert res["reference_group_sizes"]["Microglia"] == 60 and res["reference_group_sizes"]["Macro_pDC"] == 40
    assert res["n_cells_missing_reference"] == 0
    assert res["reference_key"] == "cell_type" and res["use_as_prediction"] == "majority_voting"
    assert res["reference"].startswith("https://github.com/Teichlab/celltypist/blob/fe357564a6625d3b1732a022fd39f18e55696e80/")


CASES = [  # (fixture, reference_key, use_as_prediction, executor reference dir or None)
    ("builtin_bestmatch", "cell_type", "majority_voting", ML_REF / "builtin_bestmatch"),
    ("builtin_bestmatch", "cell_type", "predicted_labels", ML_REF / "builtin_bestmatch"),
    ("builtin_probmatch_p07", "cell_type", "majority_voting", None),
    ("builtin_probmatch_p07", "cell_type", "predicted_labels", None),
    ("custom_seed0_probmatch", "cell_type", "majority_voting", ML_REF / "custom_seed0_probmatch"),
    ("custom_seed0_probmatch", "cell_type", "predicted_labels", ML_REF / "custom_seed0_probmatch"),
    ("builtin_bestmatch_no_mv", "cell_type", "predicted_labels", None),
    ("main_demo2000_bestmatch", "cell_type", "majority_voting", MAIN_REF),
    ("main_demo2000_bestmatch", "cell_type", "predicted_labels", None),
    ("relabelled_reference", "expert_annotation", "majority_voting", None),
    ("relabelled_reference", "expert_annotation", "predicted_labels", None),
    ("relabelled_reference", "leiden_int", "majority_voting", None),
    ("relabelled_reference", "leiden_int", "predicted_labels", None),
]


@pytest.mark.parametrize("fixture,ref,pred,executor_dir", CASES,
                         ids=[f"{c[0]}-{c[1]}-{c[2]}" for c in CASES])
async def test_full_crosstab_and_dotplot_tables(out_base, fixture, ref, pred, executor_dir):
    res = await call({"annotation_dir": str(FIXTURES / fixture), "reference_key": ref,
                      "use_as_prediction": pred, "output_dir": str(out_base)})
    keep(res, f"{fixture}_{ref}_{pred}")
    # Full crosstab vs pd.crosstab on the in-memory to_adata() obs (ml cell 47 without .loc).
    crosstab = read_out_table(art(res, "crosstab.csv"))
    assert_same_table(crosstab, expected(fixture, f"crosstab_{ref}_{pred}"), integer=True)
    # Dot-plot size/colour data vs in-memory celltypist.dotplot(predictions, ..., return_fig=True).
    size = read_out_table(art(res, "dotplot_fraction.csv"))
    color = read_out_table(art(res, "dotplot_mean_probability.csv"))
    exp_size = expected(fixture, f"dotplot_{ref}_{pred}_size")
    exp_color = expected(fixture, f"dotplot_{ref}_{pred}_color")
    assert_same_table(size, exp_size)
    assert_same_table(color, exp_color)
    np.testing.assert_allclose(size.sum(axis=0).to_numpy(), 1.0, rtol=0, atol=1e-12)  # columns are fractions
    assert is_png(art(res, "dotplot.png"))
    if executor_dir is not None:
        if (executor_dir / f"crosstab_cell_type_x_{pred}.csv").exists():
            assert_same_table(crosstab, read_out_table(executor_dir / f"crosstab_cell_type_x_{pred}.csv"), integer=True)
        elif (executor_dir / f"crosstab_cell_type_vs_{pred}.csv").exists():
            assert_same_table(crosstab, read_out_table(executor_dir / f"crosstab_cell_type_vs_{pred}.csv"), integer=True)
        assert_same_table(size, read_out_table(executor_dir / f"dotplot_{pred}_fraction.csv"))
        assert_same_table(color, read_out_table(executor_dir / f"dotplot_{pred}_mean_probability.csv"))
    # Summary fields agree with the reference crosstab.
    exp_ct = expected(fixture, f"crosstab_{ref}_{pred}")
    assert res["n_reference_groups"] == exp_ct.shape[0] and res["n_predicted_labels"] == exp_ct.shape[1]
    assert res["reference_group_sizes"] == {str(k): int(v) for k, v in exp_ct.sum(axis=1).items()}
    assert "selected_groups" not in res
    assert [Path(a["path"]).name for a in res["artifacts"]] == [
        "dotplot.png", "dotplot_fraction.csv", "dotplot_mean_probability.csv", "crosstab.csv"]


# ------------------------------------------------------------------------------ boxplot group probabilities
@pytest.fixture
def box_spy(monkeypatch):
    """Record the DataFrame and arguments of every DataFrame.plot(...) call made by the tool."""
    from pandas.plotting import PlotAccessor
    captured = []
    original = PlotAccessor.__call__

    def spy(self, *args, **kwargs):
        captured.append({"data": self._parent.copy(), "args": args, "kwargs": dict(kwargs)})
        return original(self, *args, **kwargs)

    monkeypatch.setattr(PlotAccessor, "__call__", spy)
    return captured


async def test_boxplot_values_custom_model_cell82(out_base, box_spy):
    """ml cell 82: obs.loc[obs.cell_type == g, new_model.cell_types].plot(kind='box', rot=90, figsize=[8, 4], title=g)."""
    groups = ["Microglia", "Macro_pDC"]
    res = await call({"annotation_dir": str(FIXTURES / "custom_seed0_probmatch"), "reference_key": "cell_type",
                      "use_as_prediction": "majority_voting", "reference_groups": groups, "output_dir": str(out_base)})
    keep(res, "custom_seed0_probmatch_boxplots_cell82")
    model_types = MANIFEST["fixtures"]["custom_seed0_probmatch"]["model_cell_types"]
    assert len(model_types) == 10
    obs_prob = expected("custom_seed0_probmatch", "obs_probabilities")
    executor_values = pd.read_csv(ML_REF / "custom_seed0_probmatch/boxplot_values_Microglia_Macro_pDC.csv", index_col=0)
    assert len(box_spy) == 2
    for g, cap in zip(groups, box_spy):
        exp = obs_prob.loc[obs_prob["cell_type"] == g, model_types]
        assert list(cap["data"].columns) == model_types
        assert list(cap["data"].index) == list(exp.index)
        np.testing.assert_array_equal(cap["data"].to_numpy(), exp.to_numpy())
        ex = executor_values.loc[executor_values["cell_type"] == g, model_types]
        np.testing.assert_array_equal(cap["data"].to_numpy(), ex.to_numpy())
        assert cap["kwargs"] == {"kind": "box", "rot": 90, "figsize": [8.0, 4], "title": g}
    boxes = [Path(a["path"]) for a in res["artifacts"] if Path(a["path"]).name.startswith("probability_boxplot_")]
    assert [p.name for p in boxes] == ["probability_boxplot_01_Microglia.png", "probability_boxplot_02_Macro_pDC.png"]
    assert all(is_png(p) for p in boxes)
    # ml cell 47 form for the custom model: executor's direct crosstab rows (seed-0 model; the as-run notebook
    # model's rows in reference/notebook_asrun are identical for these two groups).
    rows = read_out_table(ML_REF / "custom_seed0_probmatch/crosstab_cell_type_x_majority_voting_Microglia_Macro_pDC.csv")
    assert_same_table(read_out_table(art(res, "crosstab.csv")).loc[groups], rows, integer=True)
    for g in groups:
        assert res["selected_groups"][g]["prediction_counts"] == {k: int(v) for k, v in rows.loc[g].items() if v > 0}


async def test_boxplot_values_builtin_model_other_groups(out_base, box_spy):
    """98-type model, groups not used by the tutorial; only the figure width differs from the tutorial call."""
    groups = ["Kupffer cells", "pDC", "Mast cells"]
    res = await call({"annotation_dir": str(FIXTURES / "builtin_bestmatch"), "reference_key": "cell_type",
                      "use_as_prediction": "predicted_labels", "reference_groups": groups, "output_dir": str(out_base)})
    keep(res, "builtin_bestmatch_boxplots_other_groups")
    model_types = MANIFEST["fixtures"]["builtin_bestmatch"]["model_cell_types"]
    assert len(model_types) == 98
    obs_prob = expected("builtin_bestmatch", "obs_probabilities")
    exp_ct = expected("builtin_bestmatch", "crosstab_cell_type_predicted_labels")
    assert len(box_spy) == 3
    for g, cap in zip(groups, box_spy):
        exp = obs_prob.loc[obs_prob["cell_type"] == g, model_types]
        assert list(cap["data"].columns) == model_types and list(cap["data"].index) == list(exp.index)
        np.testing.assert_array_equal(cap["data"].to_numpy(), exp.to_numpy())
        assert cap["kwargs"] == {"kind": "box", "rot": 90, "figsize": [24.5, 4], "title": g}
        row = exp_ct.loc[g]
        assert res["selected_groups"][g] == {"n_cells": int(row.sum()),
                                             "prediction_counts": {str(k): int(v) for k, v in row.items() if v > 0}}


async def test_plot_group_probabilities_false(out_base, box_spy):
    res = await call({"annotation_dir": str(FIXTURES / "builtin_bestmatch"), "reference_key": "cell_type",
                      "reference_groups": ["Microglia"], "plot_group_probabilities": False, "output_dir": str(out_base)})
    assert box_spy == []
    assert not any("boxplot" in a["path"] for a in res["artifacts"])
    assert res["selected_groups"]["Microglia"]["n_cells"] == 60


# --------------------------------------------------------------------------------------- changed inputs
async def test_relabelled_reference_column_and_groups(out_base):
    """Renamed + recoded reference column: no dependence on 'cell_type' or the tutorial group names."""
    recode = MANIFEST["fixtures"]["relabelled_reference"]["recode"]
    groups = [recode["Microglia"], recode["Macro_pDC"]]
    res = await call({"annotation_dir": str(FIXTURES / "relabelled_reference"), "reference_key": "expert_annotation",
                      "use_as_prediction": "majority_voting", "reference_groups": groups, "output_dir": str(out_base)})
    keep(res, "relabelled_reference_expert_annotation")
    nb47 = notebook_cell47()  # relabelling the same cells cannot change their counts
    for g, original in zip(groups, ["Microglia", "Macro_pDC"]):
        assert res["selected_groups"][g]["prediction_counts"] == {k: int(v) for k, v in nb47.loc[original].items() if v > 0}
    assert any(Path(a["path"]).name == f"probability_boxplot_01_{groups[0]}.png" for a in res["artifacts"])
    # the original column is gone from this fixture: asking for it must fail
    msg = await call_error({"annotation_dir": str(FIXTURES / "relabelled_reference"), "reference_key": "cell_type",
                            "output_dir": str(out_base)})
    assert "reference_key 'cell_type' is not an obs column" in msg
    msg = await call_error({"annotation_dir": str(FIXTURES / "relabelled_reference"), "reference_key": "expert_annotation",
                            "reference_groups": ["Microglia"], "output_dir": str(out_base)})
    assert "['Microglia'] not found" in msg


async def test_integer_cluster_reference(out_base):
    """Integer-coded reference column; groups are given as strings."""
    res = await call({"annotation_dir": str(FIXTURES / "relabelled_reference"), "reference_key": "leiden_int",
                      "use_as_prediction": "predicted_labels", "reference_groups": ["3", "17"],
                      "plot_group_probabilities": False, "output_dir": str(out_base)})
    exp_ct = expected("relabelled_reference", "crosstab_leiden_int_predicted_labels")
    for g in ("3", "17"):
        row = exp_ct.loc[int(g)]
        assert res["selected_groups"][g] == {"n_cells": int(row.sum()),
                                             "prediction_counts": {str(k): int(v) for k, v in row.items() if v > 0}}


async def test_prob_match_p07_selected_groups(out_base):
    """Different annotation (prob match, p_thres=0.7): multi-label/Unassigned labels pass through the CSV round-trip."""
    exp_ct = expected("builtin_probmatch_p07", "crosstab_cell_type_predicted_labels")
    assert any("|" in str(c) for c in exp_ct.columns) and "Unassigned" in map(str, exp_ct.columns)
    groups = ["Macro_pDC", "Neutrophil-myeloid progenitor"]
    res = await call({"annotation_dir": str(FIXTURES / "builtin_probmatch_p07"), "reference_key": "cell_type",
                      "use_as_prediction": "predicted_labels", "reference_groups": groups,
                      "plot_group_probabilities": False, "output_dir": str(out_base)})
    for g in groups:
        row = exp_ct.loc[g]
        assert res["selected_groups"][g]["prediction_counts"] == {str(k): int(v) for k, v in row.items() if v > 0}


async def test_without_majority_voting(out_base):
    fx = str(FIXTURES / "builtin_bestmatch_no_mv")
    msg = await call_error({"annotation_dir": fx, "reference_key": "cell_type", "output_dir": str(out_base)})
    assert "has no 'majority_voting' column" in msg
    assert list(out_base.iterdir()) == []
    res = await call({"annotation_dir": fx, "reference_key": "cell_type", "use_as_prediction": "predicted_labels",
                      "output_dir": str(out_base)})
    # best-match per-cell labels do not depend on majority voting
    assert_same_table(read_out_table(art(res, "crosstab.csv")),
                      expected("builtin_bestmatch", "crosstab_cell_type_predicted_labels"), integer=True)


# ------------------------------------------------------------------------------------------------ errors
def _copy_fixture(tmp_path: Path, name: str = "builtin_bestmatch") -> Path:
    dest = tmp_path / f"annot_{name}"
    shutil.copytree(FIXTURES / name, dest)
    return dest


def _reverse_probability_rows(d):
    df = pd.read_csv(d / "probability_matrix.csv", index_col=0)
    df.iloc[::-1].to_csv(d / "probability_matrix.csv")


def _permute_h5ad_obs(d):
    a = ad.read_h5ad(d / "annotated.h5ad")
    a[np.random.default_rng(0).permutation(a.n_obs)].copy().write_h5ad(d / "annotated.h5ad")


def _drop_decision_row(d):
    df = pd.read_csv(d / "decision_matrix.csv", index_col=0)
    df.iloc[:-1].to_csv(d / "decision_matrix.csv")


def _duplicate_label_id(d):
    df = pd.read_csv(d / "predicted_labels.csv", index_col=0)
    df.index = [df.index[0]] + list(df.index[:-1])
    df.to_csv(d / "predicted_labels.csv")


def _drop_decision_column(d):
    df = pd.read_csv(d / "decision_matrix.csv", index_col=0)
    df.drop(columns=df.columns[0]).to_csv(d / "decision_matrix.csv")


def _remove_probability_file(d):
    (d / "probability_matrix.csv").unlink()


def _remove_predicted_labels_column(d):
    df = pd.read_csv(d / "predicted_labels.csv", index_col=0)
    df.rename(columns={"predicted_labels": "labels"}).to_csv(d / "predicted_labels.csv")


CORRUPTIONS = [
    (_reverse_probability_rows, "probability_matrix.csv do not align"),
    (_permute_h5ad_obs, "annotated.h5ad obs_names do not align"),
    (_drop_decision_row, "decision_matrix.csv do not align"),
    (_duplicate_label_id, "are not unique"),
    (_drop_decision_column, "list different model cell types"),
    (_remove_probability_file, "lacks ['probability_matrix.csv']"),
    (_remove_predicted_labels_column, "has no 'predicted_labels' column"),
]


@pytest.mark.parametrize("corrupt,fragment", CORRUPTIONS, ids=[c[0].__name__.strip("_") for c in CORRUPTIONS])
async def test_misaligned_or_incomplete_annotation_dir(tmp_path, out_base, corrupt, fragment):
    d = _copy_fixture(tmp_path)
    corrupt(d)
    msg = await call_error({"annotation_dir": str(d), "reference_key": "cell_type",
                            "use_as_prediction": "predicted_labels", "output_dir": str(out_base)})
    assert fragment in msg
    assert list(out_base.iterdir()) == []


@pytest.mark.parametrize("args,fragment", [
    ({"reference_key": "author_cell_type"}, "reference_key 'author_cell_type' is not an obs column"),
    ({"reference_key": "cell_type", "reference_groups": ["Microglia", "Not_a_type"]}, "['Not_a_type'] not found"),
    ({"reference_key": "cell_type", "use_as_prediction": "over_clustering"}, "use_as_prediction"),
], ids=["missing_reference_column", "missing_group", "unsupported_prediction_column"])
async def test_input_errors(out_base, args, fragment):
    msg = await call_error({"annotation_dir": str(FIXTURES / "builtin_bestmatch"), "output_dir": str(out_base), **args})
    assert fragment in msg
    assert list(out_base.iterdir()) == []


async def test_missing_directory(tmp_path, out_base):
    msg = await call_error({"annotation_dir": str(tmp_path / "does_not_exist"), "reference_key": "cell_type",
                            "output_dir": str(out_base)})
    assert "annotation_dir does not exist" in msg
    assert list(out_base.iterdir()) == []


# ------------------------------------------------------------------ input immutability and repeated calls
async def test_inputs_unchanged_and_repeated_calls_isolated(out_base):
    fx = FIXTURES / "custom_seed0_probmatch"
    before = tree_hashes(fx)
    args = {"annotation_dir": str(fx), "reference_key": "cell_type", "reference_groups": ["Microglia"],
            "output_dir": str(out_base)}
    r1 = await call(args)
    r2 = await call(args)
    assert tree_hashes(fx) == before
    assert r1["output_dir"] != r2["output_dir"]
    for r in (r1, r2):
        assert Path(r["output_dir"]).parent == out_base.resolve()
    p1 = {a["path"] for a in r1["artifacts"]}
    p2 = {a["path"] for a in r2["artifacts"]}
    assert len(p1) == len(r1["artifacts"]) == 5 and not p1 & p2
    assert all(Path(p).is_file() and Path(p).parent == Path(r1["output_dir"]) for p in p1)
    for name in ("crosstab.csv", "dotplot_fraction.csv", "dotplot_mean_probability.csv"):
        assert art(r1, name).read_bytes() == art(r2, name).read_bytes()


async def test_default_output_base(monkeypatch, tmp_path):
    import tempfile
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path / "systmp"))
    (tmp_path / "systmp").mkdir()
    res = await call({"annotation_dir": str(FIXTURES / "builtin_bestmatch"), "reference_key": "cell_type"})
    assert Path(res["output_dir"]).parent == (tmp_path / "systmp/celltypist_mcp_outputs").resolve()


# ------------------------------------------------------------------------------ tool format and transport
async def test_tool_schema_and_format():
    async with Client(celltypist_tutorial_ml_mcp) as client:
        tools = {t.name: t for t in await client.list_tools()}
    assert set(tools) == {TOOL}
    schema = tools[TOOL].input_schema
    assert set(schema["required"]) == {"annotation_dir", "reference_key"}
    props = schema["properties"]
    assert set(props) == {"annotation_dir", "reference_key", "use_as_prediction", "reference_groups",
                          "plot_group_probabilities", "output_dir"}
    assert all(p.get("description") for p in props.values())
    assert props["use_as_prediction"]["enum"] == ["majority_voting", "predicted_labels"]
    assert props["use_as_prediction"]["default"] == "majority_voting"
    assert props["plot_group_probabilities"]["type"] == "boolean" and props["plot_group_probabilities"]["default"] is True
    assert len(tools[TOOL].description.strip().splitlines()) == 2


async def test_no_stdout_output(capfd, out_base):
    await call({"annotation_dir": str(FIXTURES / "builtin_bestmatch"), "reference_key": "cell_type",
                "reference_groups": ["Microglia"], "output_dir": str(out_base)})
    assert capfd.readouterr().out == ""


# ------------------------------------------------------------------------------------------- composition
async def test_composition_with_celltypist_annotate_cells(tmp_path, out_base):
    """celltypist_annotate_cells (other module, read-only) -> its output_dir -> celltypist_compare_with_reference."""
    from src.tools.celltypist_tutorial import celltypist_tutorial_mcp
    model = DATA / "celltypist_home/data/models/Immune_All_Low.pkl"
    async with Client(celltypist_tutorial_mcp) as client:
        annot = (await client.call_tool("celltypist_annotate_cells", {
            "adata_path": str(INPUTS / "demo_500_cells.h5ad"), "model": str(model), "mode": "best match",
            "majority_voting": True, "output_dir": str(tmp_path / "annotate")})).data
    annotation_dir = Path(annot["output_dir"])
    assert all((annotation_dir / f).is_file() for f in CONTRACT_FILES)
    res = await call({"annotation_dir": str(annotation_dir), "reference_key": "cell_type",
                      "use_as_prediction": "majority_voting", "reference_groups": ["Microglia", "Macro_pDC"],
                      "output_dir": str(out_base)})
    keep(res, "composition_annotate_then_compare")
    crosstab = read_out_table(art(res, "crosstab.csv"))
    assert_same_table(crosstab.loc[["Microglia", "Macro_pDC"]], notebook_cell47(), integer=True)
    assert_same_table(crosstab, expected("builtin_bestmatch", "crosstab_cell_type_majority_voting"), integer=True)
    assert_same_table(read_out_table(art(res, "dotplot_fraction.csv")),
                      expected("builtin_bestmatch", "dotplot_cell_type_majority_voting_size"))
    assert_same_table(read_out_table(art(res, "dotplot_mean_probability.csv")),
                      expected("builtin_bestmatch", "dotplot_cell_type_majority_voting_color"))
