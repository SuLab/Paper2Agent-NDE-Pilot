"""Shared helpers/fixtures for the independent verification tests of module celltypist_tutorial (imported by conftest.py).

Expected values come from (a) the stage-2 executor's saved reference under notebooks/celltypist_tutorial/reference/
(produced by direct upstream CellTypist calls) or (b) new direct upstream calls made inside the tests. The wrapper's
own outputs are never used as expectations.

CellTypist fixes its model cache path when it is first imported (celltypist/models.py:16-20), so CELLTYPIST_FOLDER
is exported here before anything imports celltypist, and an autouse fixture also pins celltypist.models.models_path
in case another conftest imported celltypist first.
"""
import hashlib
import os
import shutil
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA = PROJECT_ROOT / "tests" / "data" / "celltypist_tutorial"
INPUTS = DATA / "inputs"
MODELS = DATA / "models"
CACHE_HOME = DATA / "celltypist_home"
CACHE_MODELS = CACHE_HOME / "data" / "models"
RESULTS = PROJECT_ROOT / "tests" / "results" / "celltypist_tutorial"
NB = PROJECT_ROOT / "notebooks" / "celltypist_tutorial"
REF = NB / "reference"
VERIFY_TMP = PROJECT_ROOT / "tmp" / "verify_celltypist_tutorial"
PROJECT_PYTHON = PROJECT_ROOT / "celltypist-env" / "bin" / "python"

os.environ["CELLTYPIST_FOLDER"] = str(CACHE_HOME)
os.environ.setdefault("MPLCONFIGDIR", str(VERIFY_TMP / ".mpl"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(VERIFY_TMP / ".numba"))
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from celltypist import models as ct_models  # noqa: E402
from fastmcp import Client  # noqa: E402
from fastmcp.exceptions import ToolError  # noqa: E402

from tools.celltypist_tutorial import celltypist_tutorial_mcp  # noqa: E402

# Seeded cache content (byte copies of notebooks/celltypist_tutorial/data/, see tests/data/.../fixture_provenance.json).
SEEDED_CACHE = {
    "Immune_All_Low.pkl": "290874d35dac039d4c9218c343fde4aac1077709b72a331ce7266f6828c36502",
    "models.json": "f617fcd7863b08e6ba2dca1db959d60d6f82f06a0833897e8968c9164a18b132",
}

# Tolerances justified by the executor's replay (reports/executed_notebook_celltypist_tutorial.json -> replay):
# decision atol = max(1e-9, float64 forward-error bound of the model dot products; 1.79e-10 for demo_2000 x Immune_All_Low,
# 3.8e-11 for demo_400 x seeded model); probability = expit(decision), |expit'| <= 1/4 -> atol/4.
DECISION_ATOL = 1e-9
PROBABILITY_ATOL = 2.5e-10

CONTRACT_FILES = {"predicted_labels.csv", "decision_matrix.csv", "probability_matrix.csv", "annotated.h5ad"}


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cache_listing() -> dict:
    return {p.name: sha256(p) for p in sorted(CACHE_MODELS.iterdir()) if p.is_file()}


@pytest.fixture(autouse=True)
def _pinned_seeded_cache(monkeypatch):
    """Point celltypist at the seeded test cache and fail any test that adds/changes cache files (e.g. downloads)."""
    monkeypatch.setattr(ct_models, "celltypist_path", str(CACHE_HOME))
    monkeypatch.setattr(ct_models, "data_path", str(CACHE_HOME / "data"))
    monkeypatch.setattr(ct_models, "models_path", str(CACHE_MODELS))
    before = cache_listing()
    assert before == SEEDED_CACHE, f"seeded test cache is not pristine: {before}"
    yield
    after = cache_listing()
    assert after == before, f"the model cache changed during the test: {before} -> {after}"


@pytest.fixture
def out_base(request) -> Path:
    """Fresh base output directory for this test under tests/results/celltypist_tutorial/."""
    p = RESULTS / request.node.name.replace("[", "_").replace("]", "").replace("/", "_")
    if p.exists():
        shutil.rmtree(p)
    p.mkdir(parents=True)
    return p


async def call_tool(name: str, arguments: dict) -> dict:
    async with Client(celltypist_tutorial_mcp) as client:
        result = await client.call_tool(name, arguments)
    assert isinstance(result.data, dict), result
    return result.data


async def call_tool_error(name: str, arguments: dict) -> str:
    """Call a tool that must fail; return the MCP tool-error text."""
    async with Client(celltypist_tutorial_mcp) as client:
        with pytest.raises(ToolError) as info:
            await client.call_tool(name, arguments)
    return str(info.value)


def assert_artifacts(result: dict, root: Path, min_count: int = 1) -> list[Path]:
    """Return-format contract: message/reference/artifacts with absolute existing paths inside a fresh output_dir."""
    assert isinstance(result["message"], str) and result["message"]
    assert result["reference"].startswith(
        "https://github.com/Teichlab/celltypist/blob/fe357564a6625d3b1732a022fd39f18e55696e80/docs/notebook/celltypist_tutorial.ipynb")
    arts = result["artifacts"]
    assert len(arts) >= min_count
    paths = []
    out = Path(result["output_dir"])
    assert out.is_absolute() and out.is_dir() and out.parent == root.resolve()
    for a in arts:
        assert set(a) == {"description", "path"} and a["description"]
        p = Path(a["path"])
        assert p.is_absolute() and p.is_file() and p.parent == out, a
        paths.append(p)
    return paths


def read_labels(path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0, dtype=str, keep_default_na=False)


def read_matrix(path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0, converters={0: str})


def assert_annotation_matches(tool_dir: Path, labels: pd.DataFrame, decision: pd.DataFrame, probability: pd.DataFrame):
    """Compare the tool's to_table CSVs with expected upstream tables: labels exact, matrices within the replay tolerances."""
    got_labels = read_labels(tool_dir / "predicted_labels.csv")
    exp_labels = labels.astype(str)
    assert list(got_labels.columns) == list(exp_labels.columns)
    assert got_labels.index.equals(exp_labels.index.astype(str))
    for col in exp_labels.columns:
        mism = (got_labels[col].to_numpy() != exp_labels[col].to_numpy()).sum()
        assert mism == 0, f"{col}: {mism} label mismatches"
    for name, exp, atol in (("decision_matrix.csv", decision, DECISION_ATOL), ("probability_matrix.csv", probability, PROBABILITY_ATOL)):
        got = read_matrix(tool_dir / name)
        assert list(got.columns) == [str(c) for c in exp.columns], name
        assert got.index.equals(exp.index.astype(str)), name
        g, e = got.to_numpy(dtype=float), exp.to_numpy(dtype=float)
        assert np.isfinite(g).all() and np.isfinite(e).all(), name
        np.testing.assert_allclose(g, e, rtol=0, atol=atol, err_msg=name)


def fresh_tmp(tag: str) -> Path:
    """New empty directory under tmp/verify_celltypist_tutorial/."""
    import tempfile
    VERIFY_TMP.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=f"{tag}_", dir=VERIFY_TMP)).resolve()


def seeded_cache_copy(tag: str) -> Path:
    """A private copy of the seeded cache (for subprocess runs that must not touch the shared test cache)."""
    home = fresh_tmp(tag) / "celltypist_home"
    shutil.copytree(CACHE_HOME, home)
    return home


def run_tools_in_subprocess(cache_home: Path, calls: list, tag: str, timeout: int = 900):
    """Run tool calls in a fresh interpreter with CELLTYPIST_FOLDER=cache_home; returns (records, stdout_bytes, stderr, rc)."""
    import json
    import subprocess
    work = fresh_tmp(tag)
    calls_path, out_path = work / "calls.json", work / "result.json"
    calls_path.write_text(json.dumps(calls))
    env = dict(os.environ)
    env["CELLTYPIST_FOLDER"] = str(cache_home)
    proc = subprocess.run(
        [str(PROJECT_PYTHON), "-I", str(Path(__file__).resolve().parent / "_tool_subprocess.py"),
         str(cache_home), str(calls_path), str(out_path)],
        capture_output=True, timeout=timeout, env=env, cwd=str(work))
    (work / "stderr.txt").write_bytes(proc.stderr)
    records = json.loads(out_path.read_text()) if out_path.is_file() else None
    return records, proc.stdout, proc.stderr.decode(errors="replace"), proc.returncode


def network_available() -> bool:
    import requests
    try:
        r = requests.head("https://celltypist.cog.sanger.ac.uk/models/models.json", timeout=15)
        return r.status_code == 200
    except Exception:
        return False
