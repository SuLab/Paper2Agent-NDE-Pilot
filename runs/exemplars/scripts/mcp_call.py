"""Call tools on an MCP server over stdio, the way an MCP client would, and save every result.

Usage: python mcp_call.py SERVER_PY SERVER_PYTHON PLAN_JSON RESULTS_JSON
PLAN_JSON is a list of {"name", "tool", "arguments"}; an argument string "$<name>.<key>" is replaced by
that key of an earlier step's result (e.g. "$de.output_dir"). Extra env vars come from PLAN env "MCP_ENV_*".
"""
import asyncio
import json
import os
import sys
import time

from fastmcp import Client
from fastmcp.client.transports import StdioTransport

server, python, plan_path, out_path = sys.argv[1:5]
plan = json.load(open(plan_path))
env = {k[len("MCP_ENV_"):]: v for k, v in os.environ.items() if k.startswith("MCP_ENV_")}
env["PATH"] = os.environ.get("PATH", "")


def resolve(value, results):
    if isinstance(value, str) and value.startswith("$"):
        step, _, key = value[1:].partition(".")
        node = results[step]
        for part in key.split("."):
            node = node[int(part)] if isinstance(node, list) else node[part]
        return node
    if isinstance(value, list):
        return [resolve(v, results) for v in value]
    return value


async def main():
    results, log = {}, []
    transport = StdioTransport(command=python, args=[server], env=env, cwd=os.path.dirname(os.path.dirname(server)))
    async with Client(transport, timeout=1800) as client:
        tools = sorted(t.name for t in await client.list_tools())
        print("tools:", tools, file=sys.stderr)
        for step in plan:
            args = {k: resolve(v, results) for k, v in step["arguments"].items()}
            t0 = time.time()
            res = await client.call_tool(step["tool"], args)
            results[step["name"]] = res.data
            log.append({"name": step["name"], "tool": step["tool"], "arguments": args,
                        "seconds": round(time.time() - t0, 1), "result": res.data})
            print(f"{step['name']}: {step['tool']} ok in {time.time() - t0:.1f}s", file=sys.stderr)
    json.dump({"server": server, "python": python, "tools": tools, "calls": log}, open(out_path, "w"), indent=2, default=str)


asyncio.run(main())
