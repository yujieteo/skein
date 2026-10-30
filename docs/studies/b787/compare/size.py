#!/usr/bin/env python3
"""Counts the lines and the tokens of each file given.

Lines are non-blank lines. "Code" lines also drop the full-line comments:
`//` in .skn, `**` in .inp and `#` in .py. Tokens use tiktoken o200k_base
on the whole file.

    size.py FILE...
"""

import sys
from pathlib import Path

import tiktoken

COMMENT = {".skn": "//", ".inp": "**", ".py": "#"}


def main():
    enc = tiktoken.get_encoding("o200k_base")
    total = [0, 0, 0]
    print("file,lines,code_lines,tokens")
    for name in sys.argv[1:]:
        path = Path(name)
        text = path.read_text()
        lines = [x.strip() for x in text.splitlines() if x.strip()]
        mark = COMMENT.get(path.suffix, None)
        code = [x for x in lines if not (mark and x.startswith(mark))]
        row = [len(lines), len(code), len(enc.encode(text))]
        total = [a + b for a, b in zip(total, row)]
        print(f"{path.name},{row[0]},{row[1]},{row[2]}")
    print(f"total,{total[0]},{total[1]},{total[2]}")


if __name__ == "__main__":
    main()
