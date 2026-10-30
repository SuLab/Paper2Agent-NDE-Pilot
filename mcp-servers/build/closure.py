"""Print the pinned runtime dependency closure of the given root distributions, as installed in this interpreter.

Usage: <project-python> closure.py ROOT [ROOT ...] [--git NAME=URL@COMMIT ...]
Extras markers are evaluated against the installed environment; dev-only tools are not roots and so are excluded.
"""
import re
import sys
from importlib.metadata import PackageNotFoundError, distribution

from packaging.markers import default_environment
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

args = sys.argv[1:]
git = dict(a.split("=", 1) for a in args if "=" in a and not a.startswith("--"))
roots = [a for a in args if "=" not in a and not a.startswith("--")]
env = default_environment()
seen, order, queue = {}, [], [(r, ()) for r in roots]
while queue:
    name, extras = queue.pop()
    key = canonicalize_name(name)
    try:
        dist = distribution(name)
    except PackageNotFoundError:
        continue
    prev = seen.get(key)
    if prev is not None and set(extras) <= prev:
        continue
    seen[key] = (prev or set()) | set(extras)
    if prev is None:
        order.append((key, dist.version))
    for spec in dist.requires or []:
        req = Requirement(spec)
        if req.marker:
            ok = any(req.marker.evaluate({**env, "extra": e}) for e in (extras or ("",)))
            if not ok:
                continue
        queue.append((req.name, tuple(req.extras)))
for key, version in sorted(order):
    print(git[key] if key in git else f"{key}=={version}")
