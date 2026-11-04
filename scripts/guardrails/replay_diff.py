"""Compares the captured outputs of two skein builds (scripts/replay.sh).

    replay_diff.py [--accept] BASE_DIR HEAD_DIR
    replay_diff.py --changed BASE_DIR HEAD_DIR     list the models whose C differs

Each model has a .c (the generated C) or a .fail (why there is none), and,
when its C differs, a .toon (the results of skein run) or a .runfail. The C
must be the same, except for the line numbers of the source in the calls of
the runtime. A result number may differ by the tolerance that its row states.
"""

from __future__ import annotations

import csv
import difflib
import re
import sys
from pathlib import Path

SOURCE_LINE = re.compile(r"\b(gs_file\d+), (?:\d+|\(\w+\) \? \d+ : \d+)\)")
NUMBER = re.compile(r"^\s*((?:[A-Za-z_]\w*:\s*)?)(-?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*(.*)$")
TABLE = re.compile(r"^(\s*)([A-Za-z_]\w*)\[\d+\]\{([^}]*)\}:$")
WITHIN = re.compile(r"within\s+(\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*(%|\S*)")
LEFT_OUT = ("seconds",)
LEFT_OUT_LINES = ("peak_mb", "estimate_mb", "built ")


def without_source_lines(c: str) -> str:
    return SOURCE_LINE.sub(r"\1, LINE)", c)


def number_and_unit(field: str) -> tuple[float, str] | None:
    match = NUMBER.match(field)
    return (float(match.group(2)), match.group(1) + match.group(3)) if match else None


def tolerance(row: dict[str, str], value: float) -> float:
    within = WITHIN.search(row.get("limit", ""))
    if not within:
        return 0.0
    amount = float(within.group(1))
    return abs(value) * amount / 100 if within.group(2) == "%" else amount


def split_rows(text: str) -> list[tuple[dict[str, str], str]]:
    rows: list[tuple[dict[str, str], str]] = []
    columns: list[str] = []
    for line in text.splitlines():
        if line.startswith(LEFT_OUT_LINES):
            continue
        table = TABLE.match(line)
        if table:
            columns = table.group(3).split(",")
            rows.append(({}, line))
        elif columns and line.startswith("  "):
            values = next(csv.reader([line.strip()]))
            rows.append((dict(zip(columns, values)), line))
        else:
            columns = []
            rows.append(({}, line))
    return rows


def differences(base: str, head: str) -> list[str]:
    a, b = split_rows(base), split_rows(head)
    if len(a) != len(b):
        return [f"the results have {len(a)} rows in the base, {len(b)} in the branch"]
    found: list[str] = []
    for (ra, la), (rb, lb) in zip(a, b):
        if la == lb:
            continue
        if not ra or ra.keys() != rb.keys():
            found.append(f"- {la}\n+ {lb}")
            continue
        for column, value in ra.items():
            other = rb[column]
            if column in LEFT_OUT or value == other:
                continue
            x, y = number_and_unit(value), number_and_unit(other)
            if x and y and x[1] == y[1] and abs(x[0] - y[0]) <= tolerance(ra, max(abs(x[0]), abs(y[0]))):
                continue
            found.append(f"- {la}\n+ {lb}")
            break
    return found


def compare(base_dir: Path, head_dir: Path) -> tuple[int, list[str]]:
    reports: list[str] = []
    compared = 0
    for head_file in sorted(head_dir.rglob("*")):
        if not head_file.is_file():
            continue
        relative = head_file.relative_to(head_dir)
        base_file = base_dir / relative
        if not base_file.exists():
            if any(base_dir.glob(f"{relative.with_suffix('')}.*")):
                reports.append(f"{relative}: the base made other outputs than the branch")
            continue
        compared += 1
        base_text, head_text = base_file.read_text(), head_file.read_text()
        if relative.suffix == ".c":
            base_text, head_text = without_source_lines(base_text), without_source_lines(head_text)
        if relative.suffix == ".c" and base_text != head_text:
            diff = difflib.unified_diff(
                base_text.splitlines(), head_text.splitlines(), "base", "branch", lineterm="", n=1
            )
            reports.append(f"{relative}: the generated C differs\n" + "\n".join(list(diff)[:40]))
        elif relative.suffix == ".toon":
            found = differences(base_text, head_text)
            if found:
                reports.append(f"{relative}: results beyond the tolerance\n" + "\n".join(found[:10]))
        elif relative.suffix in (".fail", ".runfail") and base_text != head_text:
            reports.append(f"{relative}: the error differs\n- {base_text.strip()}\n+ {head_text.strip()}")
    return compared, reports


def changed(base_dir: Path, head_dir: Path) -> list[str]:
    models: list[str] = []
    for head_file in sorted(head_dir.rglob("*.c")):
        relative = head_file.relative_to(head_dir)
        base_file = base_dir / relative
        if base_file.exists() and (
            without_source_lines(base_file.read_text()) != without_source_lines(head_file.read_text())
        ):
            models.append(relative.with_suffix(".skn").as_posix())
    return models


def main() -> int:
    args = sys.argv[1:]
    if args[:1] == ["--changed"] and len(args) == 3:
        print("\n".join(changed(Path(args[1]), Path(args[2]))))
        return 0
    accept = "--accept" in args
    paths = [a for a in args if a != "--accept"]
    if len(paths) != 2:
        print("usage: replay_diff.py [--accept] BASE_DIR HEAD_DIR", file=sys.stderr)
        return 2
    compared, reports = compare(Path(paths[0]), Path(paths[1]))
    models = len({p.with_suffix("") for p in Path(paths[1]).rglob("*") if p.is_file() and p.suffix in (".c", ".fail")})
    differing = len({r.split(":", 1)[0].rsplit(".", 1)[0] for r in reports})
    for report in reports:
        print(report)
    if reports:
        print(f"replay: {differing} of {models} models differ")
    else:
        print(f"replay: ok, {models} models, {compared} outputs the same within tolerance")
    only_c = all("the generated C differs" in r.split("\n", 1)[0] for r in reports)
    return 0 if not reports or (accept and only_c) else 1


if __name__ == "__main__":
    sys.exit(main())
