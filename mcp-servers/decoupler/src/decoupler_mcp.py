"""decoupler MCP server: bulk RNA-seq filtering, PyDESeq2 differential expression and enrichment
(TF / pathway / Hallmark activity) tools bound to scverse/decoupler (commit 27a1be7)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fastmcp import FastMCP  # noqa: E402

from tools.bulk_enrichment import bulk_enrichment_mcp  # noqa: E402

mcp = FastMCP(name="decoupler")
mcp.mount(bulk_enrichment_mcp)

if __name__ == "__main__":
    mcp.run(show_banner=False)
