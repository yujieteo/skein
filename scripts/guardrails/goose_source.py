"""Checks on the Goose source of Skein: comments, file and function size, dead code.

    python3 -I scripts/guardrails/goose_source.py check   fail on a difference from limits.txt
    python3 -I scripts/guardrails/goose_source.py write   write limits.txt from the source now

The lexer follows the lexical rules of the Goose spec (line comments, nested
block comments, strings, raw strings, character literals), so it never takes
text in a string for a comment.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Literal

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src"
LIMITS = Path(__file__).with_name("limits.txt")

IDENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
DECLARATION_KEYWORDS = ("fn", "let", "var", "struct", "type")

Kind = Literal["comment", "string", "char", "ident", "punct"]


@dataclass(frozen=True)
class Token:
    kind: Kind
    text: str
    line: int


def lex(text: str) -> Iterator[Token]:
    i = 0
    line = 1
    n = len(text)
    while i < n:
        c = text[i]
        if c == "\n":
            line += 1
            i += 1
        elif c in " \t\r":
            i += 1
        elif text.startswith("//", i):
            end = text.find("\n", i)
            end = n if end < 0 else end
            yield Token("comment", text[i:end], line)
            i = end
        elif text.startswith("/*", i):
            depth = 0
            start = i
            start_line = line
            while i < n:
                if text.startswith("/*", i):
                    depth += 1
                    i += 2
                elif text.startswith("*/", i):
                    depth -= 1
                    i += 2
                    if depth == 0:
                        break
                else:
                    line += text[i] == "\n"
                    i += 1
            yield Token("comment", text[start:i], start_line)
        elif text.startswith('"""', i):
            end = text.find('"""', i + 3)
            end = n if end < 0 else end + 3
            yield Token("string", text[i:end], line)
            line += text.count("\n", i, end)
            i = end
        elif c in "\"'":
            start = i
            i += 1
            while i < n and text[i] != c:
                i += 2 if text[i] == "\\" else 1
            i += 1
            yield Token("string" if c == '"' else "char", text[start:i], line)
        elif c.isalpha() or c == "_":
            match = IDENT.match(text, i)
            assert match is not None
            yield Token("ident", match.group(), line)
            i = match.end()
        else:
            yield Token("punct", c, line)
            i += 1


def source_files() -> list[Path]:
    return sorted(SOURCE.rglob("*.goose"))


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def comment_counts(files: list[Path]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for path in files:
        n = sum(1 for t in lex(path.read_text()) if t.kind == "comment")
        if n:
            counts[relative(path)] = n
    return counts


def function_sizes(path: Path) -> list[tuple[str, int, int]]:
    tokens = [t for t in lex(path.read_text()) if t.kind != "comment"]
    sizes: list[tuple[str, int, int]] = []
    for k, token in enumerate(tokens):
        if token.kind != "ident" or token.text != "fn" or k + 1 >= len(tokens):
            continue
        name = tokens[k + 1]
        if name.kind != "ident":
            continue
        parens = 0
        j = k + 2
        while j < len(tokens):
            t = tokens[j]
            if t.text == "(":
                parens += 1
            elif t.text == ")":
                parens -= 1
            elif parens == 0 and t.text in "{;":
                break
            j += 1
        if j >= len(tokens) or tokens[j].text != "{":
            continue
        depth = 0
        while j < len(tokens):
            if tokens[j].kind == "punct":
                depth += tokens[j].text == "{"
                depth -= tokens[j].text == "}"
                if depth == 0:
                    break
            j += 1
        sizes.append((name.text, token.line, tokens[j].line - token.line + 1))
    return sizes


def dead_declarations(files: list[Path]) -> list[str]:
    declared: dict[str, list[str]] = {}
    mentions: dict[str, int] = {}
    for path in files:
        tokens = [t for t in lex(path.read_text()) if t.kind != "comment"]
        depth = 0
        previous = ""
        for t in tokens:
            if t.kind == "punct":
                depth += t.text == "{"
                depth -= t.text == "}"
            elif t.kind == "ident":
                mentions[t.text] = mentions.get(t.text, 0) + 1
                if depth == 0 and previous in DECLARATION_KEYWORDS:
                    declared.setdefault(t.text, []).append(f"{relative(path)}:{t.line}")
            elif t.kind == "string":
                for word in IDENT.findall(t.text):
                    mentions[word] = mentions.get(word, 0) + 1
            previous = t.text if t.kind == "ident" else ""
    dead: list[str] = []
    for name, places in declared.items():
        if mentions.get(name, 0) <= len(places):
            dead += [f"{place}: {name} is declared and never named again" for place in places]
    return sorted(dead)


def measure() -> tuple[int, int]:
    files = source_files()
    biggest_file = max(len(p.read_text().splitlines()) for p in files)
    biggest_fn = max(size for p in files for _, _, size in function_sizes(p))
    return biggest_file, biggest_fn


def round_up(value: int, step: int) -> int:
    return -(-value // step) * step


def write() -> int:
    biggest_file, biggest_fn = measure()
    lines = [
        f"max_file_lines {round_up(biggest_file, 100)}",
        f"max_function_lines {round_up(biggest_fn, 50)}",
    ]
    LIMITS.write_text("\n".join(lines) + "\n")
    print(f"wrote {relative(LIMITS)}: {lines[0]}, {lines[1]}")
    return 0


def read_limits() -> tuple[int, int]:
    max_file = max_fn = 0
    for row in LIMITS.read_text().splitlines():
        word = row.split()
        if word[0] == "max_file_lines":
            max_file = int(word[1])
        elif word[0] == "max_function_lines":
            max_fn = int(word[1])
    return max_file, max_fn


def check() -> int:
    files = source_files()
    max_file, max_fn = read_limits()
    problems: list[str] = []
    for path in files:
        size = len(path.read_text().splitlines())
        if size > max_file:
            problems.append(f"{relative(path)}: {size} lines, the limit is {max_file}")
        for name, line, fn_size in function_sizes(path):
            if fn_size > max_fn:
                problems.append(
                    f"{relative(path)}:{line}: fn {name} has {fn_size} lines, the limit is {max_fn}"
                )
    for path_name, have in sorted(comment_counts(files).items()):
        problems.append(f"{path_name}: {have} comments. Give the code a name, a test or a check instead")
    problems += dead_declarations(files)
    for problem in problems:
        print(f"FAIL {problem}")
    if not problems:
        print(f"goose_source: ok ({len(files)} files)")
    return 1 if problems else 0


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) == 2 else ""
    if mode == "check":
        return check()
    if mode == "write":
        return write()
    print("usage: goose_source.py check|write", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
