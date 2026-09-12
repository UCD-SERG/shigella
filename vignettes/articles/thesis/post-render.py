#!/usr/bin/env python3
"""
Quarto project post-render hook.

WHAT IT DOES.  Runs the two post-processing scripts over the dissertation, so
that a local `quarto render` produces the file that gets filed:

  fix_submission.py   the UC Davis submission mechanics — one w:pPr per
                      caption, table captions styled as tables, updateFields,
                      roman preliminary pages and an arabic body, no blank page
  fix_table_breaks.py cantSplit on every row, so a table or figure and its
                      caption stay on one page

The CI workflow runs fix_table_breaks.py itself, after this hook has already
run inside `quarto render`; both scripts are no-ops the second time.

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
for script in ("fix_submission.py", "fix_table_breaks.py"):
    subprocess.run(
        [sys.executable, str(here / script), str(here / TARGET)], check=True
    )
