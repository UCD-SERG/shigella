#!/usr/bin/env python3
"""
Quarto project post-render hook.

WHAT IT DOES.  Runs fix_table_breaks.py over the dissertation, so that a local
`quarto render` produces the same file CI produces.  Until now the fix existed
only as the "Keep tables whole" step in render-thesis-vignettes.yaml, so the
docx on your machine and the docx in the CI artifact differed.

WHY IT IS GUARDED.  _quarto.yml is shared with chapter2.qmd and chapter3.qmd,
and a post-render hook runs for every render in the project — including
`quarto render chapter2/chapter2.qmd`.  The guard means the hook does nothing
at all unless the dissertation was one of the outputs of this render, so the two
manuscripts build exactly as they did before.

The CI workflow still applies the script itself.  Running it twice is a no-op:
fix_table_breaks.py adds cantSplit only to rows that do not already carry it.
"""
import os
import subprocess
import sys
from pathlib import Path

TARGET = "thesis_ch1_ch2_ch3.docx"

outputs = os.environ.get("QUARTO_PROJECT_OUTPUT_FILES", "")
if not any(Path(line.strip()).name == TARGET for line in outputs.splitlines() if line.strip()):
    sys.exit(0)

here = Path(__file__).resolve().parent
subprocess.run(
    [sys.executable, str(here / "fix_table_breaks.py"), str(here / TARGET)],
    check=True,
)
