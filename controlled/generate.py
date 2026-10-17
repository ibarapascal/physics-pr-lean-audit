#!/usr/bin/env python3
"""generate.py: sample generator for the controlled comparison tasks.

Description:
    For every task package (tasks/<id>/request.md + snapshot/), calls the Codex CLI K times independently. Each call
    copies the snapshot into a fresh temporary working directory; the agent may read the repository with read-only
    commands (read-only sandbox, no web search, isolated CODEX_HOME) and replies with a unified diff against the
    repository root. All events and executed commands are recorded; finished samples are skipped, so the run can be
    resumed. Only text generation happens here; no simulation is run.
Usage: GEN_CODEX_HOME=<isolated home> GEN_WORK_ROOT=<temporary directory> [CODEX_BIN=<path to codex>] \
       python3 generate.py [--k 10] [--par 2] [--only <id>]
       Outputs are written to samples/<task>/codex/sNN.{md,json}.
Dependencies: Python standard library; Codex CLI (`codex` on PATH, or CODEX_BIN)
Change Log: 2026-10-16 first version
"""
import argparse, hashlib, json, os, pathlib, re, shutil, subprocess, tempfile, time
from concurrent.futures import ThreadPoolExecutor

ROOT = pathlib.Path(__file__).resolve().parent
CODEX = os.environ.get("CODEX_BIN", "codex")
MODEL = "gpt-6-astra"
TASK = (
    "\n\nThe repository checkout in the current directory is the library source at the time of this request. "
    "You may read files with read-only shell commands (ls, cat, grep, sed -n). You cannot modify files, run simulations, "
    "or access the network.\n\nReply with:\n"
    "1. A single unified diff (```diff fenced) against the repository root, with paths relative to the repository root "
    "in the --- / +++ headers.\n2. At most 10 lines explaining the change.\n"
)


def extract_diff(text):
    m = re.findall(r"```diff\n(.*?)```", text, flags=re.S)
    return m[0] if m else None


def run_codex(prompt, workdir, home, timeout_s=3600):
    env = dict(os.environ, CODEX_HOME=str(home))
    cmd = [CODEX, "exec", "--skip-git-repo-check", "--ephemeral", "--json", "-m", MODEL,
           "-c", "model_reasoning_effort=medium", "-s", "read-only", "-C", str(workdir), prompt]
    r = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=workdir, env=env, timeout=timeout_s)
    events, text, cmds, usage, err = [], "", [], None, None
    for line in r.stdout.splitlines():
        try:
            e = json.loads(line)
        except json.JSONDecodeError:
            continue
        events.append(e)
        it = e.get("item") or {}
        if e.get("type") == "item.completed" and it.get("type") == "agent_message":
            text = it.get("text") or text
        if e.get("type") == "item.completed" and it.get("type") not in (None, "agent_message", "reasoning"):
            cmds.append({"type": it.get("type"), "command": it.get("command"), "query": it.get("query")})
        if e.get("type") == "turn.completed":
            usage = e.get("usage")
        if e.get("type") in ("error", "turn.failed"):
            err = e
    events = [{k: v for k, v in e.items() if k != "thread_id"} for e in events]  # drop per-session identifiers
    return r.returncode, text, cmds, usage, err, events, r.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--par", type=int, default=2)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    home = pathlib.Path(os.environ["GEN_CODEX_HOME"])
    assert (home / "config.toml").exists() and not (home / "AGENTS.md").exists(), "an isolated CODEX_HOME (with config.toml and without AGENTS.md) is required"
    work_root = pathlib.Path(os.environ["GEN_WORK_ROOT"]); work_root.mkdir(parents=True, exist_ok=True)
    jobs = []
    for t in sorted((ROOT / "tasks").iterdir()):
        if not (t / "request.md").exists() or (a.only and t.name != a.only):
            continue
        for k in range(1, a.k + 1):
            stem = ROOT / "samples" / t.name / "codex" / f"s{k:02d}"
            if not stem.with_suffix(".json").exists():
                jobs.append((t, k, stem))
    print(f"{len(jobs)} samples to generate", flush=True)

    def do(job):
        t, k, stem = job
        req = (t / "request.md").read_text()
        wd = pathlib.Path(tempfile.mkdtemp(prefix=f"{t.name}_s{k:02d}_", dir=work_root))
        shutil.copytree(t / "snapshot", wd, dirs_exist_ok=True)
        t0 = time.time()
        code, text, cmds, usage, err, events, stderr = run_codex(req + TASK, wd, home)
        shutil.rmtree(wd, ignore_errors=True)
        if code != 0 or err or not text:
            print(f"[failed] {t.name} s{k:02d} code={code} {str(err)[:200]} {stderr[-200:]}", flush=True)
            return
        net = [c for c in cmds if c["type"] not in ("command_execution",)]
        rec = {"task": t.name, "model": MODEL, "effort": "medium", "sample": k,
               "request_sha256": hashlib.sha256(req.encode()).hexdigest(), "wall_s": round(time.time() - t0, 1),
               "commands": cmds, "non_command_tools": net, "usage": usage, "events": events,
               "diff": extract_diff(text), "valid_format": extract_diff(text) is not None}
        stem.parent.mkdir(parents=True, exist_ok=True)
        stem.with_suffix(".md").write_text(text)
        stem.with_suffix(".json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
        print(f"[done] {t.name} s{k:02d} {rec['wall_s']}s diff={'yes' if rec['valid_format'] else 'no'} commands {len(cmds)}"
              f"{' non-command tools ' + str(len(net)) if net else ''}", flush=True)

    with ThreadPoolExecutor(max_workers=a.par) as ex:
        list(ex.map(do, jobs))


if __name__ == "__main__":
    main()
