"""Run celltypist_tutorial tools in a fresh process with its own CELLTYPIST_FOLDER (verifier helper, not a test).

Usage: python -I _tool_subprocess.py <celltypist_home> <calls.json> <result.json>
calls.json: [{"tool": name, "arguments": {...}}, ...]. Each call is made through fastmcp.Client(celltypist_tutorial_mcp);
the result (or MCP tool-error text) and the model-cache listing after the call are written to result.json.
Nothing is printed to stdout on purpose: the parent asserts that the process wrote zero bytes to fd 1.
"""
import asyncio
import json
import os
import sys
from pathlib import Path

home, calls_path, out_path = sys.argv[1:4]
os.environ["CELLTYPIST_FOLDER"] = home  # must precede the first celltypist import (models.py:16)
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from fastmcp import Client  # noqa: E402
from fastmcp.exceptions import ToolError  # noqa: E402

from tools.celltypist_tutorial import celltypist_tutorial_mcp  # noqa: E402
from celltypist import models  # noqa: E402

models_dir = Path(home) / "data" / "models"


async def main():
    calls = json.loads(Path(calls_path).read_text())
    records = [{"models_path": models.models_path}]
    async with Client(celltypist_tutorial_mcp) as client:
        for call in calls:
            rec = {"tool": call["tool"], "arguments": call["arguments"]}
            try:
                res = await client.call_tool(call["tool"], call["arguments"])
                rec["data"] = res.data
            except ToolError as e:
                rec["error"] = str(e)
            rec["cache_after"] = sorted(p.name for p in models_dir.iterdir()) if models_dir.is_dir() else None
            records.append(rec)
    Path(out_path).write_text(json.dumps(records, indent=1, default=str))


asyncio.run(main())
