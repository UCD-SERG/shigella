#!/usr/bin/env python3
"""
fix_table_breaks.py — stop tables splitting across pages in a Word file.

WHY THIS EXISTS.  Quarto wraps each captioned table and each captioned figure
in an outer single-row table: the caption paragraph and the data table (or the
image) both sit inside one cell of it.  flextable's opts_word reaches the inner
table only, so the outer wrapper can still straddle a page boundary and take
the caption with it.  A row with cantSplit set cannot straddle one, and the
wrapper has exactly one row — so setting it there moves the caption and the
table together to the next page.

WHAT IT SETS.
  w:cantSplit  on every <w:tr>, outer and inner.  Word's "Allow row to break
               across pages", off.
  w:keepNext   on every paragraph in a table except those in its last row, so
               Word keeps the rows together rather than filling the page.

BOTH ARE MERGED, NOT PREPENDED.  A <w:p> may carry at most one <w:pPr> and a
<w:tr> at most one <w:trPr>, so a property is added to the element that is
already there rather than in a second one.  Inside <w:pPr> the schema fixes the
order of the children: w:keepNext goes immediately after w:pStyle when there is
one, and first otherwise.  Adding a second <w:pPr> instead — which an earlier
version of this script did — produces markup Word may repair, and the repair
can drop the paragraph's real properties, its w:pStyle among them.

RUNNING IT TWICE IS A NO-OP.  A row that already has cantSplit and a paragraph
that already has keepNext are left alone, so it is safe for the render hook and
the CI step to both apply it.

WHAT IT CANNOT DO.  A table or a figure taller than the 9-inch text column still
has to break; nothing prevents that.  For tables the repeated header from
flextable's paginate() is what keeps them readable.

USAGE
  python fix_table_breaks.py thesis_ch1_ch2_ch3.docx
  python fix_table_breaks.py in.docx out.docx      # keep the original
"""
import re
import shutil
import sys
import zipfile
from pathlib import Path

# flextable writes <w:tbl> with a namespace declaration list; pandoc writes it
# bare.  Both have to be seen, or the table depth below is tracked wrongly.
TOKEN = re.compile(
    r"<w:tbl(?:\s[^>]*)?>|</w:tbl>"
    r"|<w:tr(?:\s[^>]*)?>|</w:tr>"
    r"|<w:p(?:\s[^>]*)?>"
)

TR_OPEN = re.compile(r"<w:tr(?:\s[^>]*)?>")
P_OPEN = re.compile(r"<w:p(?:\s[^>]*)?>")
PSTYLE = re.compile(r"<w:pStyle\b[^>]*/>")


def _row_spans(xml):
    """Return (start, end) of every <w:tr> …</w:tr>, and the set of starts that
    are the final row of the table they belong to."""
    stack = []
    rows_of = {}
    order = []
    spans = {}
    open_rows = []
    for m in TOKEN.finditer(xml):
        t = m.group(0)
        if t.startswith("<w:tbl"):
            stack.append(m.start())
            rows_of[m.start()] = []
            order.append(m.start())
        elif t == "</w:tbl>":
            if stack:
                stack.pop()
        elif TR_OPEN.fullmatch(t):
            if stack:
                rows_of[stack[-1]].append(m.start())
            open_rows.append(m.start())
        elif t == "</w:tr>":
            if open_rows:
                spans[open_rows.pop()] = m.end()
    last = {rows_of[k][-1] for k in order if rows_of[k]}
    return spans, last


def patch_document_xml(xml: str):
    """Return (patched_xml, n_rows_touched, n_paras_touched)."""
    spans, last_rows = _row_spans(xml)

    # Which paragraphs sit inside a table but not inside a table's final row?
    inside = []
    depth = 0
    in_last = 0
    for m in TOKEN.finditer(xml):
        t = m.group(0)
        if t.startswith("<w:tbl"):
            depth += 1
        elif t == "</w:tbl>":
            depth -= 1
        elif TR_OPEN.fullmatch(t):
            in_last += 1 if m.start() in last_rows else 0
        elif t == "</w:tr>":
            in_last = max(0, in_last - 1) if in_last else 0
        elif P_OPEN.fullmatch(t) and depth > 0 and not in_last:
            inside.append((m.start(), m.end()))

    # Build the patched string in one pass, back to front, so earlier offsets
    # stay valid.
    edits = []   # (position, old_len, replacement)
    rows = paras = 0

    for start, end in spans.items():
        seg = xml[start:end]
        m = re.match(r"(<w:tr(?:\s[^>]*)?>)\s*(<w:trPr>)", seg)
        if m:
            body_end = seg.index("</w:trPr>")
            if "cantSplit" in seg[:body_end]:
                continue                      # already set — leave it alone
            pos = start + m.end(2)
            edits.append((pos, 0, "<w:cantSplit/>"))
        else:
            m2 = TR_OPEN.match(seg)
            pos = start + m2.end()
            edits.append((pos, 0, "<w:trPr><w:cantSplit/></w:trPr>"))
        rows += 1

    for start, end in inside:
        rest = xml[end:end + 4000]
        m = re.match(r"\s*<w:pPr>", rest)
        if m:
            close = rest.index("</w:pPr>")
            head = rest[:close]
            if "<w:keepNext" in head:
                continue                      # already set — leave it alone
            # w:keepNext follows w:pStyle in the schema's element order.
            ps = PSTYLE.search(head)
            at = end + (ps.end() if ps else m.end())
            edits.append((at, 0, "<w:keepNext/>"))
        else:
            edits.append((end, 0, "<w:pPr><w:keepNext/></w:pPr>"))
        paras += 1

    for pos, old_len, text in sorted(edits, reverse=True):
        xml = xml[:pos] + text + xml[pos + old_len:]

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
