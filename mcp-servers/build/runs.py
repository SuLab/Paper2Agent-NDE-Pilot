"""Coordinator helper: maintain <project>/reports/agent-runs.json (workflow schema 1) from observed lifecycle events.

  runs.py init  PROJECT URL COMMIT TUTORIAL_JSON
  runs.py add   PROJECT --id ID --role ROLE --task TASK --agent AGENT --start T --finish T --status S
                --owned p1,p2 --report R [--produced f1,f2] [--impl RUN_ID] [--tested f1,f2]
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(project):
    return json.loads((project / "reports/agent-runs.json").read_text())


def save(project, state):
    (project / "reports/agent-runs.json").write_text(json.dumps(state, indent=2) + "\n")


def split(value):
    return [v for v in (value or "").split(",") if v]


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("init")
    i.add_argument("project", type=Path)
    i.add_argument("url")
    i.add_argument("commit")
    i.add_argument("tutorials", help="JSON list of {id, module, source}")
    a = sub.add_parser("add")
    a.add_argument("project", type=Path)
    for f in ("id", "role", "task", "agent", "start", "finish", "status", "owned", "report"):
        a.add_argument(f"--{f}", required=f != "finish")
    for f in ("produced", "impl", "tested"):
        a.add_argument(f"--{f}")
    args = p.parse_args()
    root = args.project
    if args.cmd == "init":
        tutorials = json.loads(args.tutorials)
        for t in tutorials:
            t["sha256"] = sha(root / t["source"])
        state = {"schema_version": 1, "route": "python",
                 "repository": {"url": args.url, "commit": args.commit},
                 "tutorials": tutorials, "runs": [], "exclusions": []}
        save(root, state)
        return
    state = load(root)
    assert args.id not in {r["id"] for r in state["runs"]}, f"duplicate run id {args.id}"
    run = {"id": args.id, "role": args.role, "task": args.task, "agent_id": args.agent,
           "status": args.status, "started_at": args.start, "finished_at": args.finish,
           "owned_paths": split(args.owned), "report": args.report,
           "report_sha256": sha(root / args.report)}
    if args.produced:
        run["produced_files"] = {f: sha(root / f) for f in split(args.produced)}
    if args.impl:
        run["implementation_run_id"] = args.impl
    if args.tested:
        run["tested_files"] = {f: sha(root / f) for f in split(args.tested)}
    state["runs"].append(run)
    save(root, state)
    print(json.dumps(run, indent=1)[:600])


if __name__ == "__main__":
    sys.exit(main())
