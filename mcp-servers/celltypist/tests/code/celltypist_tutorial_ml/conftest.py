"""Environment setup for the celltypist_tutorial_ml verification tests.

Fixtures and expected values are produced by build_fixtures.py (direct upstream calls only);
the tests never derive expectations from the wrapper's own output.
"""
import os
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[3]
_DATA = _ROOT / "tests/data/celltypist_tutorial_ml"
# celltypist fixes its model directory on first import. The composition test passes an absolute model
# path, so it does not depend on which test module imported celltypist first.
os.environ.setdefault("CELLTYPIST_FOLDER", str(_DATA / "celltypist_home"))
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# Test-strength check only: CTML_MUTANT=<file> makes the tests import a mutated copy of the tool module
# (kept under tmp/verify_celltypist_tutorial_ml/mutants/) instead of the production file. Unset in normal runs.
if os.environ.get("CTML_MUTANT"):
    import importlib.util

    import src.tools  # noqa: F401
    _spec = importlib.util.spec_from_file_location("src.tools.celltypist_tutorial_ml", os.environ["CTML_MUTANT"])
    _mod = importlib.util.module_from_spec(_spec)
    sys.modules[_spec.name] = _mod
    _spec.loader.exec_module(_mod)


@pytest.fixture(scope="session", autouse=True)
def _ctml_fixtures_built():
    names = ("builtin_bestmatch", "builtin_bestmatch_no_mv", "builtin_probmatch_p07",
             "custom_seed0_probmatch", "relabelled_reference")
    files = ("predicted_labels.csv", "decision_matrix.csv", "probability_matrix.csv", "annotated.h5ad")
    missing = [f"{n}/{f}" for n in names for f in files if not (_DATA / "fixtures" / n / f).is_file()]
    if missing or not (_DATA / "expected").is_dir():
        pytest.fail(f"Fixtures missing ({missing[:3]}...). Run tests/code/celltypist_tutorial_ml/build_fixtures.py first.")


@pytest.fixture
def out_base(tmp_path):
    base = tmp_path / "outputs"
    base.mkdir()
    return base
