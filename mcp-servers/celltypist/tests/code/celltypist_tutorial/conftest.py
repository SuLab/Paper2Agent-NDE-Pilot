"""pytest wiring for tests/code/celltypist_tutorial (helpers live in _celltypist_tutorial_support.py)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _celltypist_tutorial_support import _pinned_seeded_cache, out_base  # noqa: E402,F401  (fixtures)


def pytest_configure(config):
    config.addinivalue_line("markers", "network: needs https://celltypist.cog.sanger.ac.uk (skipped when unreachable)")
