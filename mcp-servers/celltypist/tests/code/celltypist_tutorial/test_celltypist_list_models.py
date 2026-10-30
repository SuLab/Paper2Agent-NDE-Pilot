"""Verification of celltypist_list_models (tutorial cells 16-20: models.models_description()).

Expected values: notebooks/celltypist_tutorial/reference/models_description.csv (executor's direct upstream call on the
same cached models.json) and a new direct celltypist.models call in this process.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

from _celltypist_tutorial_support import (CACHE_MODELS, REF, assert_artifacts, call_tool, celltypist_tutorial_mcp,
                                          fresh_tmp, network_available, run_tools_in_subprocess)
from celltypist import models as ct_models
from fastmcp import Client

TOOL = "celltypist_list_models"


async def test_schema_and_required_inputs():
    async with Client(celltypist_tutorial_mcp) as client:
        tools = {t.name: t for t in await client.list_tools()}
    schema = tools[TOOL].input_schema
    assert set(schema.get("required", [])) == set()
    assert set(schema["properties"]) == {"refresh_index", "output_dir"}
    assert schema["properties"]["refresh_index"]["default"] is False
    assert all(p.get("description") for p in schema["properties"].values())
    doc = tools[TOOL].description.strip().splitlines()
    assert len(doc) == 2


async def test_listing_matches_executor_reference_and_direct_upstream(out_base):
    result = await call_tool(TOOL, {"output_dir": str(out_base)})
    (csv,) = assert_artifacts(result, out_base, 1)
    assert csv.name == "models_description.csv"

    got = pd.read_csv(csv, dtype=str, keep_default_na=False)
    ref = pd.read_csv(REF / "models_description.csv", dtype=str, keep_default_na=False)
    direct = ct_models.models_description().astype(str)
    assert list(got.columns) == ["model", "description"]
    assert got.shape == (61, 2) == ref.shape
    pd.testing.assert_frame_equal(got, ref)
    pd.testing.assert_frame_equal(got, direct.reset_index(drop=True))

    # Returned summary fields are the same upstream values.
    assert result["n_models"] == 61
    assert result["default_model"] == ct_models.get_default_model() == "Immune_All_Low.pkl"
    assert result["models"] == direct.to_dict("records")
    assert Path(result["index_path"]) == (CACHE_MODELS / "models.json").resolve()
    row = got[got.model == "Immune_All_Low.pkl"]
    assert row.description.tolist() == ["immune sub-populations combined from 20 tissues of 18 studies"]


async def test_repeated_calls_write_fresh_directories(out_base):
    r1 = await call_tool(TOOL, {"output_dir": str(out_base)})
    r2 = await call_tool(TOOL, {"output_dir": str(out_base)})
    assert r1["output_dir"] != r2["output_dir"]
    p1, p2 = r1["artifacts"][0]["path"], r2["artifacts"][0]["path"]
    assert p1 != p2 and Path(p1).read_bytes() == Path(p2).read_bytes()
    assert len([p for p in out_base.iterdir() if p.is_dir()]) == 2


@pytest.mark.network
def test_refresh_index_downloads_only_the_index():
    """refresh_index=True re-downloads models.json (get_models_index(force_update=True)) and never any model file."""
    if not network_available():
        pytest.skip("celltypist.cog.sanger.ac.uk unreachable")
    home = fresh_tmp("net_list_models") / "celltypist_home"
    out = fresh_tmp("net_list_models_out")
    records, stdout, stderr, rc = run_tools_in_subprocess(
        home, [{"tool": TOOL, "arguments": {"refresh_index": True, "output_dir": str(out)}}], "net_list_models_run")
    assert rc == 0, stderr[-3000:]
    assert stdout == b"", stdout[:500]
    rec = records[1]
    assert "error" not in rec, rec
    assert rec["cache_after"] == ["models.json"]
    index = json.loads((home / "data" / "models" / "models.json").read_text())
    assert rec["data"]["n_models"] == len(index["models"])
    listed = pd.read_csv(rec["data"]["artifacts"][0]["path"], dtype=str, keep_default_na=False)
    assert listed.model.tolist() == [m["filename"] for m in index["models"]]
    assert listed.description.tolist() == [m["details"] for m in index["models"]]
