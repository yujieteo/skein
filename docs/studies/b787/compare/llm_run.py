#!/usr/bin/env python3
"""Runs scripts/llm/workflow.py with the B787 request in place of its own.

The request is one line plus the problem sheet (../problem.md). All else
(the model, the card, the check, the repair rounds, the run) is the workflow.

    docs/studies/b787/compare/llm_run.py skein --out DIR
    docs/studies/b787/compare/llm_run.py inp --out DIR
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "scripts" / "llm"))
import workflow  # noqa: E402

LINE = ("Model the Boeing 787 forward fuselage pressurization of this "
        "problem sheet, and report the outputs and the checks it asks for.")
workflow.REQUEST = LINE + "\n\n" + (HERE.parent / "problem.md").read_text()

if __name__ == "__main__":
    workflow.main()
