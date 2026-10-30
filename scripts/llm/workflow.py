#!/usr/bin/env python3
"""The end-to-end LLM workflow (docs/llm_test.md, docs/readability.md).

A fresh model session gets one request and writes one input file. The driver
checks the file, gives the checker output back to the model until the check
passes or the round limit, then runs the file and records the result.

    scripts/llm/workflow.py skein --out DIR   # the model gets `skein help --llm`
    scripts/llm/workflow.py inp --out DIR     # a CalculiX deck for femahoy, no card

Path `skein`: the check is `skein check`, the run is `skein run`.
Path `inp`: the check and the run are one `fh run` (femahoy gives the deck to
ccx, which is the parser of the deck).

The model is the claude CLI in print mode, with no tools, no settings and no
memory. Each run is one new session. The repair rounds resume that session.
DIR gets the files of each round and results.toon.
"""

import argparse
import json
import re
import subprocess
import sys
import time
import uuid
from datetime import date
from pathlib import Path

REQUEST = ("a steel cantilever 1 m x 50 mm x 100 mm, fixed at one end, "
           "1 kN tip load, report the tip deflection")
SYSTEM = ("You write simulation input files. Answer with the complete file "
          "in one fenced code block.")
ASK = {
    "skein": "Write the complete Skein model.",
    "inp": ("Write the complete CalculiX input deck. femahoy solves it with "
            "`fh run deck.inp`, which gives the deck to ccx. Put the mesh in "
            "the deck."),
}
FIX = {
    "skein": "skein check printed this. Fix the model.",
    "inp": "fh run printed this. Fix the deck.",
}
SUFFIX = {"skein": ".skn", "inp": ".inp"}
ROOT = Path(__file__).resolve().parents[2]


def ask_model(model, cwd, session, first, prompt):
    """Sends one prompt. Returns the reply text and the usage record."""
    cmd = ["claude", "-p", "--output-format", "json", "--model", model,
           "--tools", "", "--setting-sources", "", "--strict-mcp-config",
           "--disable-slash-commands", "--system-prompt", SYSTEM]
    cmd += ["--session-id", session] if first else ["--resume", session]
    out = subprocess.run(cmd + [prompt], cwd=cwd, capture_output=True,
                         text=True)
    try:
        reply = json.loads(out.stdout)
    except ValueError:
        reply = {}
    if out.returncode or reply.get("is_error") or "result" not in reply:
        sys.exit(f"workflow: the model call failed (exit {out.returncode}): "
                 f"{reply.get('result') or out.stderr or out.stdout}")
    return reply["result"], reply


def last_code_block(text):
    blocks = re.findall(r"```[^\n]*\n(.*?)```", text, re.S)
    return blocks[-1] if blocks else text


def tokens(reply):
    u = reply["usage"]
    inp = (u["input_tokens"] + u.get("cache_creation_input_tokens", 0)
           + u.get("cache_read_input_tokens", 0))
    return inp, u["output_tokens"]


def run(cmd, cwd):
    out = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return out.returncode, out.stdout + out.stderr


def skein_value(text):
    """The first scalar length and the first reaction in the report records."""
    value = reaction = ""
    for line in text.splitlines():
        m = re.match(r'\s+[^,"]+,("[^"]*"|[^,]*),(.*)$', line)
        if not m:
            continue
        expr, val = (x.strip('"') for x in m.groups())
        if not reaction and expr.startswith("reaction("):
            reaction = val
        elif not value and re.fullmatch(r"-?[\d.eE+-]+ (mm|m|cm|um)", val):
            value = val
    return value, reaction


def inp_value(text):
    """The largest displacement component in job.dat, and the sum of the
    first block of forces (the reaction of the first set that prints RF)."""
    m = re.search(r"text results (\S+\.dat)", text)
    if not m or not Path(m.group(1)).exists():
        return "", ""
    value, force, block, done = 0.0, [0.0, 0.0, 0.0], "", False
    for line in Path(m.group(1)).read_text().splitlines():
        if " for set " in line and "time" in line:
            done = done or block == "f"
            block = ("u" if "displacements" in line
                     else "f" if "force" in line else "")
            continue
        nums = line.split()
        try:
            v = [float(x) for x in nums[-3:]]
        except ValueError:
            continue
        if len(v) < 3:
            continue
        if block == "u":
            value = max(value, *v, key=abs)
        elif block == "f" and not done:
            force = [a + b for a, b in zip(force, v)]
    reaction = "(" + ", ".join(f"{f:.6g}" for f in force) + ") deck units"
    return (f"{value:.6g} deck units" if value else ""), reaction


def one_run(path, n, args, card):
    work = args.out / f"{path}-{n}"
    work.mkdir(parents=True, exist_ok=True)
    session = str(uuid.uuid4())
    prompt = "\n\n".join(x for x in (card, REQUEST, ASK[path]) if x)
    record = {"path": path, "run": n, "first_try_pass": "no", "passed": "no",
              "repair_rounds": 0, "input_tokens": 0, "output_tokens": 0}
    start = time.monotonic()
    output = ""
    for k in range(1, args.rounds + 1):
        text, reply = ask_model(args.model, work, session, k == 1, prompt)
        record["model"] = ",".join(reply.get("modelUsage", {})) or args.model
        i, o = tokens(reply)
        record["input_tokens"] += i
        record["output_tokens"] += o
        file = work / f"round-{k}{SUFFIX[path]}"
        file.write_text(last_code_block(text))
        if path == "skein":
            status, output = run([str(args.skein), "check", file.name], work)
        else:
            status, output = run(["fh", "run", file.name, "--timeout",
                                  str(args.timeout)], work)
        (work / f"check-{k}.txt").write_text(output)
        if status == 0:
            record["passed"] = "yes"
            record["first_try_pass"] = "yes" if k == 1 else "no"
            break
        record["repair_rounds"] = k
        prompt = f"{FIX[path]}\n\n{output}"
    if record["passed"] == "yes":
        if path == "skein":
            status, output = run([str(args.skein), "run", file.name], work)
            (work / "run.txt").write_text(output)
            record["value"], record["reaction"] = skein_value(output)
            record["run_status"] = status
        else:
            record["value"], record["reaction"] = inp_value(output)
            record["run_status"] = 0
    record["wall_s"] = round(time.monotonic() - start, 1)
    return record


FIELDS = ["path", "run", "model", "first_try_pass", "passed", "repair_rounds",
          "input_tokens", "output_tokens", "wall_s", "run_status", "value",
          "reaction"]


def toon(records):
    rows = [",".join(json.dumps(r.get(f, "")) if isinstance(r.get(f, ""), str)
                     else str(r.get(f, "")) for f in FIELDS) for r in records]
    head = f"runs[{len(records)}]{{{','.join(FIELDS)}}}:"
    return "\n".join([f"date: {date.today()}", f"request: {json.dumps(REQUEST)}",
                      head] + ["  " + r for r in rows]) + "\n"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("path", choices=["skein", "inp"])
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--runs", type=int, default=3)
    p.add_argument("--rounds", type=int, default=5)
    p.add_argument("--model", default="claude-opus-5-5")
    p.add_argument("--skein", type=Path, default=ROOT / "build" / "skein")
    p.add_argument("--timeout", type=int, default=600, help="fh run, seconds")
    args = p.parse_args()
    args.out = args.out.resolve()
    args.skein = args.skein.resolve()
    card = ""
    if args.path == "skein":
        status, card = run([str(args.skein), "help", "--llm"], ROOT)
        if status != 0:
            sys.exit(f"workflow: skein help --llm failed:\n{card}")
    records = []
    for n in range(1, args.runs + 1):
        records.append(one_run(args.path, n, args, card))
        print(toon(records[-1:]).splitlines()[-1], flush=True)
        (args.out / f"results-{args.path}.toon").write_text(toon(records))
    print(toon(records), end="")


if __name__ == "__main__":
    main()
