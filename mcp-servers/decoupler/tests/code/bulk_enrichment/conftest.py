"""Shared helpers for the bulk_enrichment verification tests (verifier attempt 2).

Expectations come from the stage-2 executor reference files and from direct upstream calls
made inside the tests; wrapper outputs are never used as their own expectation.
"""

import asyncio
import hashlib
import json
import re
import sys
import time
import uuid
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SRC = PROJECT_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from fastmcp import Client  # noqa: E402

from tools import bulk_enrichment as be  # noqa: E402

SERVER = be.bulk_enrichment_mcp
PROJECT_PYTHON = PROJECT_ROOT / "decoupler-env" / "bin" / "python"
DATA = PROJECT_ROOT / "tests" / "data" / "bulk_enrichment"
INP = DATA / "inputs"
DER = DATA / "derived"
NETS = INP / "networks"
NB = PROJECT_ROOT / "notebooks" / "bulk_rna"
REF = NB / "reference"
RESULTS = PROJECT_ROOT / "tests" / "results" / "bulk_enrichment"
RUN_ID = time.strftime("%Y%m%d-%H%M%S") + "_" + uuid.uuid4().hex[:6]

RAW_H5AD = INP / "hsctgfb_raw.h5ad"
FILTERED_H5AD = INP / "hsctgfb_filtered_reference.h5ad"
DE_CSV = INP / "deseq2_results_reference.csv"
SIG_SCORES_CSV = INP / "collectri_scores_significant_reference.csv"
NETWORKS = ("collectri", "progeny", "hallmark")
CONTRAST = "treatment.vs.control"
DE_TIMEOUT = 900  # seconds; PyDESeq2 took 28-125 s per fit in stage 2
CALL_TIMEOUT = 300


def pytest_configure(config):
    config.addinivalue_line("markers", "network: downloads live prior-knowledge networks (Zenodo/OmniPath)")


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _check_evidence() -> None:
    manifest = json.loads((REF / "manifest.json").read_text())
    for entry in manifest["files"]:
        p = PROJECT_ROOT / entry["path"]
        if p.suffix in {".h5ad", ".parquet", ".csv", ".txt"} and "sha256" in entry:
            assert sha256(p) == entry["sha256"], f"reference evidence changed: {entry['path']}"
    fixtures = json.loads((DATA / "fixtures_manifest.json").read_text())
    for section in fixtures.values():
        for rel, info in section.items():
            assert sha256(DATA / rel) == info["sha256"], f"fixture changed: {rel}"


@pytest.fixture(scope="session", autouse=True)
def evidence_integrity():
    """Reference files and fixtures match their recorded hashes before and after all tool calls
    (no tool may modify a user input, network snapshot or reference file)."""
    _check_evidence()
    yield
    _check_evidence()


def _sanitize(nodeid: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", nodeid.split("::", 1)[-1])[:120]


@pytest.fixture
def out_dir(request) -> Path:
    d = RESULTS / RUN_ID / _sanitize(request.node.nodeid)
    d.mkdir(parents=True, exist_ok=True)
    return d


def module_out_dir(name: str) -> Path:
    d = RESULTS / RUN_ID / name
    d.mkdir(parents=True, exist_ok=True)
    return d


async def _call(name: str, args: dict, timeout: float):
    async with Client(SERVER) as client:
        return await client.call_tool(name, args, timeout=timeout, raise_on_error=False)


def error_text(result) -> str:
    return "\n".join(getattr(block, "text", "") for block in result.content)


def call_tool(name: str, args: dict, timeout: float = CALL_TIMEOUT) -> dict:
    """Successful call through the FastMCP client; returns result.data."""
    result = asyncio.run(_call(name, args, timeout))
    assert not result.is_error, f"{name} returned a tool error: {error_text(result)}"
    data = result.data
    assert isinstance(data, dict)
    return data


def call_tool_error(name: str, args: dict, timeout: float = CALL_TIMEOUT) -> str:
    """Expected failure: must come back as an MCP tool error (is_error); returns its text."""
    result = asyncio.run(_call(name, args, timeout))
    assert result.is_error, f"{name} succeeded although an error was expected: {result.data}"
    return error_text(result)


def list_tools():
    async def go():
        async with Client(SERVER) as client:
            return await client.list_tools()

    return asyncio.run(go())


def check_tool_format(name: str, required: set[str]) -> None:
    """Python tool format: described params, Literal enums, two-line docstring, exact required set."""
    tool = {t.name: t for t in list_tools()}[name]
    schema = tool.input_schema
    assert set(schema.get("required", [])) == required
    for pname, prop in schema["properties"].items():
        assert isinstance(prop.get("description"), str) and prop["description"].strip(), pname
    lines = [ln for ln in (tool.description or "").strip().splitlines() if ln.strip()]
    assert len(lines) == 2, tool.description


def check_common_return(data: dict, root: Path, min_artifacts: int) -> list[Path]:
    assert {"message", "reference", "artifacts"} <= set(data)
    assert data["reference"] == be.REFERENCE
    assert isinstance(data["message"], str) and data["message"]
    paths = []
    for art in data["artifacts"]:
        assert set(art) >= {"description", "path"} and art["description"]
        p = Path(art["path"])
        assert p.is_absolute() and p.exists(), p
        assert p.resolve().is_relative_to(root.resolve()), p
        paths.append(p)
    assert len(paths) >= min_artifacts
    assert len(set(paths)) == len(paths)
    json.dumps(data)  # JSON-serialisable
    return paths


def artifact(data: dict, suffix: str) -> Path:
    hits = [Path(a["path"]) for a in data["artifacts"] if a["path"].endswith(suffix)]
    assert len(hits) == 1, (suffix, [a["path"] for a in data["artifacts"]])
    return hits[0]


def assert_png(path: Path, min_bytes: int = 2000) -> None:
    raw = Path(path).read_bytes()
    assert raw[:8] == b"\x89PNG\r\n\x1a\n" and len(raw) > min_bytes, path


def assert_frame_exact(actual: pd.DataFrame, expected: pd.DataFrame, check_names: bool = True) -> None:
    """Identifiers, order and float values bitwise equal (NaN positions must match)."""
    assert list(map(str, actual.index)) == list(map(str, expected.index))
    assert list(map(str, actual.columns)) == list(map(str, expected.columns))
    a = actual.to_numpy(dtype=float)
    e = expected.to_numpy(dtype=float)
    assert np.array_equal(np.isnan(a), np.isnan(e))
    np.testing.assert_array_equal(a, e)


class Spy:
    """Record calls to an upstream function while still running it."""

    def __init__(self):
        self.calls: dict[str, list[dict]] = {}

    def wrap(self, mp: pytest.MonkeyPatch, owner, attr: str) -> None:
        orig = getattr(owner, attr)
        record = self.calls.setdefault(attr, [])

        def spy(*args, **kwargs):
            # Snapshot mutable inputs at call time (e.g. AnnData later subset in place).
            def snap_value(v):
                return v.copy() if isinstance(v, (pd.DataFrame, ad.AnnData)) else v

            snap = {k: snap_value(v) for k, v in kwargs.items()}
            snap["_args"] = [snap_value(a) for a in args]
            record.append(snap)
            return orig(*args, **kwargs)

        mp.setattr(owner, attr, spy)


def read_de_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def read_lines(path: Path) -> list[str]:
    return [ln for ln in Path(path).read_text().splitlines() if ln]


def direct_leading_edge(de: pd.DataFrame, net: pd.DataFrame, source: str, stat: str = "stat"):
    """Direct dc.pl.leading_edge calls (cells 50 / 60) used as independent expectations."""
    import decoupler as dc
    import matplotlib.pyplot as plt

    out = {}
    if "weight" in net.columns:
        parts = {"positive_weights": net[net["weight"] > 0], "negative_weights": net[net["weight"] < 0]}
    else:
        parts = {"all_targets": net}
    for label, sub in parts.items():
        if not sub.loc[sub["source"] == source, "target"].isin(de.index).any():
            continue
        fig, le = dc.pl.leading_edge(de, stat=stat, net=sub, name=source, return_fig=True)
        plt.close(fig)
        out[label] = [str(g) for g in le]
    return out
