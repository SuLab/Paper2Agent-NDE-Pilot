"""Stdio launcher for the celltypist_tutorial tool module (verifier helper for the protocol round-trip test).

CELLTYPIST_FOLDER comes from the launching environment, as it would for the packaged server.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from tools.celltypist_tutorial import celltypist_tutorial_mcp  # noqa: E402

if __name__ == "__main__":
    celltypist_tutorial_mcp.run(transport="stdio", show_banner=False)
