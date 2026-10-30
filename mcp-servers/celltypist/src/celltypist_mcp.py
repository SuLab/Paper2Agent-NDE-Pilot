"""CellTypist MCP server: cell type annotation tools bound to Teichlab/celltypist (commit fe357564).

Built-in models are cached under $CELLTYPIST_FOLDER/data/models (default ~/.celltypist); set the
variable before starting the server to choose the cache location.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fastmcp import FastMCP  # noqa: E402

from tools.celltypist_tutorial import celltypist_tutorial_mcp  # noqa: E402
from tools.celltypist_tutorial_ml import celltypist_tutorial_ml_mcp  # noqa: E402

mcp = FastMCP(name="celltypist")
mcp.mount(celltypist_tutorial_mcp)
mcp.mount(celltypist_tutorial_ml_mcp)

if __name__ == "__main__":
    mcp.run(show_banner=False)
