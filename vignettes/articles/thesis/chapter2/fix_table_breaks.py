#!/usr/bin/env python3
"""
fix_table_breaks.py — stop tables splitting across pages in a Word file.

WHY THIS EXISTS.  Quarto wraps each captioned table in an outer single-row
table: the caption paragraph and the data table both sit inside one cell of it.
flextable's opts_word reaches the inner table only, so the outer wrapper can
still straddle a page boundary and take the caption with it.  A row with
cantSplit set cannot straddle one, and the wrapper has exactly one row — so
setting it there moves the caption and the table together to the next page.

WHAT IT SETS.
  w:cantSplit  on every <w:tr>, outer and inner.  Word's "Allow row to break
               across pages", off.
  w:keepNext   on every paragraph in a table except those in its last row, so
               Word keeps the rows together rather than filling the page.

WHAT IT CANNOT DO.  A table taller than a page still has to break; nothing
prevents that.  For those the repeated header from flextable's paginate() is
what keeps them readable.

USAGE
  python fix_table_breaks.py thesis_ch1_ch2_ch3.docx
  python fix_table_breaks.py in.docx out.docx      # keep the original
"""
import re
import shutil
import sys
import zipfile
from pathlib import Path

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def patch_document_xml(xml: str):
    """Return (patched_xml, n_rows_touched, n_paras_touched)."""
    rows = paras = 0

    # ---- 1. cantSplit on every row -------------------------------------
    # A row either has <w:trPr>…</w:trPr> already or it does not.
    def add_to_existing(m):
        nonlocal rows
        inner = m.group(1)
        if "cantSplit" in inner:
            return m.group(0)
        rows += 1
        return "<w:trPr><w:cantSplit/>" + inner + "</w:trPr>"

    xml = re.sub(r"<w:trPr>(.*?)</w:trPr>", add_to_existing, xml, flags=re.S)

    # rows with no trPr at all: trPr must be the first child of w:tr
    def add_missing(m):
        nonlocal rows
        rows += 1
        return "<w:tr><w:trPr><w:cantSplit/></w:trPr>"

    xml = re.sub(r"<w:tr>(?!\s*<w:trPr>)", add_missing, xml)

    # ---- 2. keepNext on paragraphs inside tables ------------------------
    # Walk the string tracking table depth and row index so the last row of a
    # table is left alone; keeping it with the next paragraph would drag the
    # following body text up against the table.
    out = []
    pos = 0
    depth = 0
    row_spans = []          # (start, end) of the final row of each open table

    # Find each table's last row so it can be skipped.
    last_rows = set()
    for tbl in re.finditer(r"<w:tbl>", xml):
        pass  # spans computed below with a proper scan

    tok = re.compile(r"<w:tbl>|</w:tbl>|<w:tr>|</w:tr>|<w:p>|<w:p/>|</w:p>")
    stack = []
    rows_of = {}
    order = []
    for m in tok.finditer(xml):
        t = m.group(0)
        if t == "<w:tbl>":
            stack.append(m.start())
            rows_of[m.start()] = []
            order.append(m.start())
        elif t == "</w:tbl>":
            if stack:
                stack.pop()
        elif t == "<w:tr>" and stack:
            rows_of[stack[-1]].append(m.start())
    for k in order:
        if rows_of[k]:
            last_rows.add(rows_of[k][-1])

    # Now insert keepNext into every <w:p> that is inside a table but not
    # inside a table's final row.
    result = []
    idx = 0
    depth = 0
    cur_row_start = None
    in_last_row = False
    for m in tok.finditer(xml):
        t = m.group(0)
        result.append(xml[idx:m.start()])
        idx = m.end()
        if t == "<w:tbl>":
            depth += 1
            result.append(t)
        elif t == "</w:tbl>":
            depth -= 1
            result.append(t)
        elif t == "<w:tr>":
            cur_row_start = m.start()
            in_last_row = m.start() in last_rows
            result.append(t)
        elif t == "</w:tr>":
            in_last_row = False
            result.append(t)
        elif t == "<w:p>" and depth > 0 and not in_last_row:
            paras += 1
            result.append("<w:p><w:pPr><w:keepNext/></w:pPr>")
        else:
            result.append(t)
    result.append(xml[idx:])
    xml = "".join(result)

    return xml, rows, paras


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    src = Path(sys.argv[1])
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else src
    if not src.exists():
        sys.exit(f"not found: {src}")

    tmp = src.with_suffix(".patched.docx")
    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(
        tmp, "w", zipfile.ZIP_DEFLATED
    ) as zout:
        rows = paras = 0
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/document.xml":
                xml = data.decode("utf-8")
                xml, rows, paras = patch_document_xml(xml)
                data = xml.encode("utf-8")
            zout.writestr(item, data)

    shutil.move(str(tmp), str(dst))
    print(f"  rows given cantSplit : {rows}")
    print(f"  paragraphs given keepNext : {paras}")
    print(f"  written: {dst}")


if __name__ == "__main__":
    main()
