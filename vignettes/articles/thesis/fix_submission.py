#!/usr/bin/env python3
"""
fix_submission.py — the docx changes UC Davis Graduate Studies requires that
Quarto cannot express in the .qmd.

Everything here is layout or document plumbing.  Nothing touches prose.

WHAT IT DOES

1.  ONE w:pPr PER CAPTION.  Quarto emits every caption paragraph with two
    <w:pPr> elements — one holding w:jc, a second holding w:spacing, w:jc and
    w:pStyle.  A w:p permits exactly one.  LibreOffice merges them; Word may
    keep the first and discard the second, which is the one carrying
    w:pStyle.  The two are merged into one, with w:pStyle first and the rest in
    the order CT_PPrBase fixes, so the caption style survives whichever way a
    reader resolves it.

2.  TABLE CAPTIONS GET THE TABLE CAPTION STYLE.  Quarto gives all 43 captions —
    figures and tables alike — the ImageCaption style.  The List of Tables is a
    TOC field switched on `\\t "Table Caption"`, so as shipped it can never
    list anything, and the List of Figures lists the tables as well as the
    figures.  Captions whose text begins "Table" are moved to the TableCaption
    style, which reference.docx defines with the same formatting.

3.  THE THREE FRONT LISTS.  See front_lists() below — the two list field
    codes are repaired, and either list is inserted if it is absent.

4.  UPDATE FIELDS ON OPEN.  The three front lists are TOC fields with no cached
    result.  Without w:updateFields Word never offers to fill them and they
    open empty.  Pandoc writes its own settings.xml rather than copying the
    reference document's, so this cannot be set in reference.docx.

5.  ROMAN PRELIMINARY PAGES, ARABIC BODY.  UC Davis wants the preliminary pages
    in lowercase roman with the title page as i, and the body restarting at
    arabic 1.  That needs two sections, and a .qmd has no way to ask for one.
    The page-break paragraph before the Introduction heading is turned into a
    next-page section break carrying the roman numbering; the body's numbering
    comes from the sectPr in reference.docx.

    START_ROMAN_AT is 2 because the UC Davis title page is merged in
    separately and is page i.  If that title page turns out to be more than one
    page, raise this by one for each extra page.

6.  NO BLANK PAGES.  UC Davis prohibits them.  {{< pagebreak >}} compiles to an
    empty paragraph carrying <w:br w:type="page"/>.  That paragraph needs a
    line of room on the page it starts on; when the text before it ends at a
    page boundary there is none, so the paragraph moves to a fresh page, its
    break then moves to the one after, and the fresh page is left blank.

    Where a paragraph follows, the empty paragraph is dropped and that
    paragraph is given w:pageBreakBefore: the break happens in the same place
    and there is no paragraph left that could be stranded.  Where a table
    follows — Quarto wraps every captioned figure and table in one — the empty
    paragraph is kept but set to a one-point line, because w:pageBreakBefore
    inside a table cell is not honoured consistently.

    This is deliberately structural rather than a fix to the one blank page in
    the current build.  Filling the three front lists moves every page after
    them, and a spot fix would only move the problem.

7.  THE ABSTRACT HEADING IS CENTRED, and only that one.  See
    centre_abstract().

8.   THE QUARTO TITLE BLOCK IS DROPPED.  See strip_title_block().

9.   THE ABSTRACT COMES BEFORE THE THREE LISTS.  See front_matter_order().

10.  EACH FRONT LIST STARTS A NEW PAGE.  See break_before_front_lists().

RUNNING IT TWICE IS A NO-OP.
"""
import re
import shutil
import sys
import zipfile
from pathlib import Path

# The UC Davis title page is page i and is merged in separately, so this
# document's first page is ii.  One more per extra title page.
START_ROMAN_AT = 2

PAGEBREAK_PARA = re.compile(
    r'<w:p>\s*<w:r>\s*<w:br w:type="page"\s*/>\s*</w:r>\s*</w:p>'
)

# CT_PPrBase fixes the order of these children; anything not listed keeps the
# order it arrived in, after the ones that are.
PPR_ORDER = [
    "pStyle", "keepNext", "keepLines", "pageBreakBefore", "framePr",
    "widowControl", "numPr", "suppressLineNumbers", "pBdr", "shd", "tabs",
    "suppressAutoHyphens", "kinsoku", "wordWrap", "overflowPunct",
    "topLinePunct", "autoSpaceDE", "autoSpaceDN", "bidi", "adjustRightInd",
    "snapToGrid", "spacing", "ind", "contextualSpacing", "mirrorIndents",
    "suppressOverlap", "jc", "textDirection", "textAlignment",
    "textboxTightWrap", "outlineLvl", "divId", "cnfStyle", "rPr", "sectPr",
]


def _children(ppr_inner):
    """Split a pPr body into (tag, xml) pairs, later duplicates dropped."""
    seen = {}
    order = []
    for m in re.finditer(r"<w:(\w+)\b[^>]*?(?:/>|>.*?</w:\1>)", ppr_inner, re.S):
        tag = m.group(1)
        if tag not in seen:
            order.append(tag)
        seen[tag] = m.group(0)      # a later value wins: the second pPr is the
                                    # one Quarto means
    return [(t, seen[t]) for t in order]


def merge_caption_pprs(xml):
    """1 and 2: one pPr per paragraph, and table captions styled as tables.

    Paragraphs do not nest, so each one is taken as <w:p> up to the next
    </w:p> and rewritten in place.  Nothing outside a pPr is touched, and the
    caller asserts that.
    """
    merged = restyled = 0
    out = []
    pos = 0
    for m in re.finditer(r"<w:p>", xml):
        if m.start() < pos:
            continue
        close = xml.find("</w:p>", m.end())
        if close == -1:
            break
        para = xml[m.start():close + len("</w:p>")]
        pprs = re.findall(r"<w:pPr>.*?</w:pPr>", para, re.S)
        if len(pprs) < 2:
            continue
        # _children wants the pPr bodies, not the pPr elements themselves.
        inner = "".join(re.findall(r"<w:pPr>(.*?)</w:pPr>", para, re.S))
        kids = dict(_children(inner))
        text = re.search(r"<w:t[^>]*>([^<]*)</w:t>", para)
        if ("pStyle" in kids
                and 'w:val="ImageCaption"' in kids["pStyle"]
                and text and text.group(1).lstrip().startswith("Table")):
            kids["pStyle"] = '<w:pStyle w:val="TableCaption" />'
            restyled += 1
        ordered = [kids[t] for t in PPR_ORDER if t in kids]
        ordered += [v for t, v in kids.items() if t not in PPR_ORDER]
        rest = para
        for one in pprs:
            rest = rest.replace(one, "", 1)
        rest = rest[len("<w:p>"):]
        out.append(xml[pos:m.start()])
        out.append("<w:p><w:pPr>" + "".join(ordered) + "</w:pPr>" + rest)
        pos = close + len("</w:p>")
        merged += 1
    out.append(xml[pos:])
    return "".join(out), merged, restyled


SDT = ('<w:sdt><w:sdtPr><w:docPartObj>'
       '<w:docPartGallery w:val="Table of Contents" /><w:docPartUnique />'
       '</w:docPartObj></w:sdtPr><w:sdtContent>'
       '<w:p><w:pPr><w:pStyle w:val="TOCHeading" /></w:pPr>'
       '<w:r><w:t xml:space="preserve">{title}</w:t></w:r></w:p>'
       '<w:p><w:r><w:fldChar w:fldCharType="begin" w:dirty="true" />'
       '<w:instrText xml:space="preserve">TOC \\h \\z \\t &quot;{style},1&quot;</w:instrText>'
       '<w:fldChar w:fldCharType="separate" /><w:fldChar w:fldCharType="end" />'
       '</w:r></w:p></w:sdtContent></w:sdt>')

FRONT = [("List of Figures", "Image Caption"), ("List of Tables", "Table Caption")]


def front_lists(xml):
    """3: a table of contents, a list of figures and a list of tables, in that
    order, each with a field code Word will accept.

    Quarto emits all three, but the two list fields are switched
    `\t "Image Caption" \c` — a bare \c, which takes a SEQ identifier and is
    invalid without one — inside a content control whose docPartGallery is
    "List of Figures", which is not one of Word's galleries.  A Word that
    discards them on update leaves the table of contents alone, which is
    exactly what a docx that has been through Word looks like.  Both are
    rewritten to `\t "Style,1"` in a Table of Contents gallery, and inserted
    after the table of contents if they are not there at all.
    """
    repaired = inserted = 0
    for _, style in FRONT:
        bad = 'TOC \\h \\z \\t &quot;%s&quot; \\c' % style
        good = 'TOC \\h \\z \\t &quot;%s,1&quot;' % style
        if bad in xml:
            xml = xml.replace(bad, good)
            repaired += 1
    xml = xml.replace('<w:docPartGallery w:val="List of Figures" />',
                      '<w:docPartGallery w:val="Table of Contents" />')
    xml = xml.replace('<w:docPartGallery w:val="List of Tables" />',
                      '<w:docPartGallery w:val="Table of Contents" />')
    # Insert whichever is missing, after the last sdt already present.
    for title, style in FRONT:
        if ">%s</w:t>" % title in xml:
            continue
        end = xml.rfind("</w:sdt>")
        if end == -1:
            break
        end += len("</w:sdt>")
        xml = xml[:end] + SDT.format(title=title, style=style) + xml[end:]
        inserted += 1
    return xml, repaired, inserted


ABSTRACT_HEADING = re.compile(
    r'<w:p>\s*<w:pPr>((?:(?!</w:pPr>).)*?)</w:pPr>\s*<w:r>\s*'
    r'<w:t[^>]*>Abstract</w:t>\s*</w:r>\s*</w:p>', re.S
)


def centre_abstract(xml):
    """7: centre the word Abstract, and nothing else.

    Abstract shares the Heading1 style with Chapter 1, Chapter 2, Chapter 3 and
    the appendices, so centring the style would centre all of them.  Direct
    formatting on the one paragraph overrides the style for that paragraph
    alone, which is why this is a w:jc on the heading rather than a change to
    reference.docx.  w:jc comes late in the order CT_PPrBase fixes, so it goes
    at the end of the pPr.
    """
    m = ABSTRACT_HEADING.search(xml)
    if m is None or "Heading1" not in m.group(1):
        return xml, False
    if "<w:jc " in m.group(1):
        return xml, False                       # already centred
    at = m.start(1) + len(m.group(1))
    return xml[:at] + '<w:jc w:val="center"/>' + xml[at:], True



TITLE_STYLES = ("Title", "Subtitle", "Author", "Date")
SDT_BLOCK = re.compile(r"<w:sdt>.*?</w:sdt>", re.S)
TOC_HEADING = re.compile(
    r'(<w:p>\s*<w:pPr>\s*<w:pStyle w:val="TOCHeading" />)(\s*</w:pPr>)'
)


def strip_title_block(xml):
    """8: drop Quarto's title, subtitle, author and date paragraphs.

    The UC Davis title page is merged in separately and is the one that counts;
    this block would sit between it and the front lists.  Only the four
    paragraphs go — the YAML metadata stays, so pandoc still writes the title
    and author into docProps/core.xml, which is where ProQuest and Word's
    document properties read them from.
    """
    dropped = 0
    out, pos = [], 0
    for m in re.finditer(r"<w:p>", xml):
        if m.start() < pos:
            continue
        close = xml.find("</w:p>", m.end())
        if close == -1:
            break
        para = xml[m.start():close + len("</w:p>")]
        st = re.search(r'<w:pStyle w:val="(\w+)" />', para)
        if st and st.group(1) in TITLE_STYLES:
            out.append(xml[pos:m.start()])
            pos = close + len("</w:p>")
            dropped += 1
    out.append(xml[pos:])
    return "".join(out), dropped


def front_matter_order(xml):
    """9: Abstract first, then the three lists, then the Introduction.

    The lists are w:sdt blocks and the Abstract is ordinary content, so this
    lifts the three blocks out and puts them back immediately before the
    paragraph that carries the preliminary section's sectPr.  That paragraph is
    the last of the lowercase-roman section, so the Abstract and all three
    lists stay roman and the Introduction is still arabic 1.
    """
    blocks = SDT_BLOCK.findall(xml)
    if len(blocks) != 3:
        return xml, False
    anchor = xml.find("<w:sectPr>")           # the inserted preliminary sectPr
    if anchor == -1:
        return xml, False
    para = xml.rfind("<w:p>", 0, anchor)
    if para == -1 or para < xml.find(blocks[-1]):
        return xml, False                     # already after the lists
    for b in blocks:
        xml = xml.replace(b, "", 1)
    para = xml.rfind("<w:p>", 0, xml.find("<w:sectPr>"))
    return xml[:para] + "".join(blocks) + xml[para:], True


def break_before_front_lists(xml):
    """10: each front list starts a new page, like every other heading.

    w:pageBreakBefore goes after w:pStyle, which is where CT_PPrBase puts it.
    """
    xml, n = TOC_HEADING.subn(
        lambda m: m.group(0) if "pageBreakBefore" in m.group(0)
        else m.group(1) + "<w:pageBreakBefore/>" + m.group(2), xml)
    return xml, sum(1 for _ in re.finditer(
        r'<w:pStyle w:val="TOCHeading" /><w:pageBreakBefore/>', xml))


def section_break(xml):
    """4: turn the page break before the Introduction into a section break."""
    if "<w:pgNumType" in xml and 'w:fmt="lowerRoman"' in xml:
        return xml, False                       # already done
    body = re.search(r"<w:sectPr\b[^>]*>(.*?)</w:sectPr>", xml, re.S)
    if not body:
        return xml, False
    prelim = body.group(1)
    prelim = re.sub(r"<w:pgNumType[^>]*/>", "", prelim)
    prelim = prelim.replace(
        "<w:pgSz", '<w:type w:val="nextPage"/><w:pgSz', 1
    ).replace(
        "<w:cols",
        f'<w:pgNumType w:fmt="lowerRoman" w:start="{START_ROMAN_AT}"/><w:cols',
        1,
    )
    sect_para = ("<w:p><w:pPr>"
                 '<w:spacing w:before="0" w:after="0" w:line="20" w:lineRule="exact"/>'
                 "<w:sectPr>" + prelim + "</w:sectPr></w:pPr></w:p>")

    # The body begins at the thesis-only Introduction, not at Chapter 1.  Its
    # own {#sec-intro} label is the anchor: the heading text carries a section
    # number, so matching on the words would be brittle.
    heading = re.search(r'<w:bookmarkStart[^>]*w:name="sec-intro"\s*/>', xml)
    if heading is None:
        raise SystemExit("fix_submission: could not find the sec-intro bookmark")

    before = xml[:heading.start()]
    brk = None
    for m in PAGEBREAK_PARA.finditer(before):
        brk = m
    if brk is None:
        raise SystemExit("fix_submission: no page break before the Introduction")
    return xml[:brk.start()] + sect_para + xml[brk.end():], True


def no_blank_pages(xml):
    """5: no page-break paragraph can be left alone on a page."""
    # A one-point, zero-spaced empty paragraph: it always finds room.
    THIN = ('<w:pPr><w:spacing w:before="0" w:after="0" w:line="20" '
            'w:lineRule="exact"/><w:rPr><w:sz w:val="2"/></w:rPr></w:pPr>')
    SKIP = re.compile(r"(?:\s*<w:bookmark(?:Start|End)[^>]*/>)*\s*")
    out = []
    pos = 0
    hoisted = thinned = 0
    for m in PAGEBREAK_PARA.finditer(xml):
        if m.start() < pos:
            continue
        gap = SKIP.match(xml, m.end())
        nxt = xml[gap.end():]
        out.append(xml[pos:m.start()])
        if nxt.startswith("<w:p>"):
            lead = xml[m.end():gap.end()]
            if nxt.startswith("<w:p><w:pPr>"):
                # w:pStyle comes first in CT_PPrBase, so step over it.
                ps = re.match(r"<w:p><w:pPr>(<w:pStyle\b[^>]*/>)?", nxt)
                head = ps.group(0)
                out.append(lead + head + "<w:pageBreakBefore/>")
                pos = gap.end() + len(head)
            else:
                out.append(lead + "<w:p><w:pPr><w:pageBreakBefore/></w:pPr>")
                pos = gap.end() + len("<w:p>")
            hoisted += 1
        else:
            # a table follows: keep the paragraph, make it one point tall
            out.append("<w:p>" + THIN + m.group(0)[len("<w:p>"):])
            pos = m.end()
            thinned += 1
    out.append(xml[pos:])
    return "".join(out), hoisted, thinned


def patch_settings(xml):
    if "updateFields" in xml:
        return xml, False
    return xml.replace("<w:settings", "<w:settings", 1).replace(
        "</w:settings>", '<w:updateFields w:val="true"/></w:settings>', 1
    ), True


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    src = Path(sys.argv[1])
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else src
    if not src.exists():
        sys.exit(f"not found: {src}")

    tmp = src.with_suffix(".submission.docx")
    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(
        tmp, "w", zipfile.ZIP_DEFLATED
    ) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/document.xml":
                xml = data.decode("utf-8")
                strip = lambda t: re.sub(r"<w:pPr>.*?</w:pPr>", "", t, flags=re.S)
                before = strip(xml)
                xml, merged, restyled = merge_caption_pprs(xml)
                assert strip(xml) == before, "merge changed content outside w:pPr"
                xml, repaired, inserted = front_lists(xml)
                xml, dropped = strip_title_block(xml)
                xml, centred = centre_abstract(xml)
                xml, sectioned = section_break(xml)
                xml, reordered = front_matter_order(xml)
                xml, broken = break_before_front_lists(xml)
                xml, hoisted, thinned = no_blank_pages(xml)
                print(f"  caption pPr merged        : {merged}")
                print(f"  table captions restyled   : {restyled}")
                print(f"  list field codes repaired : {repaired}")
                print(f"  front lists inserted      : {inserted}")
                print(f"  title-block paras dropped : {dropped}")
                print(f"  abstract heading centred  : {centred}")
                print(f"  front matter reordered    : {reordered}")
                print(f"  front lists on a new page : {broken}")
                print(f"  section break inserted    : {sectioned}")
                print(f"  breaks -> pageBreakBefore  : {hoisted}")
                print(f"  breaks kept but thinned   : {thinned}")
                data = xml.encode("utf-8")
            elif item.filename == "word/settings.xml":
                xml = data.decode("utf-8")
                xml, added = patch_settings(xml)
                print(f"  updateFields added        : {added}")
                data = xml.encode("utf-8")
            zout.writestr(item, data)

    shutil.move(str(tmp), str(dst))
    print(f"  written: {dst}")


if __name__ == "__main__":
    main()
