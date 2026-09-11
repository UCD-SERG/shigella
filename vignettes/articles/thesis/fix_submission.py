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

7.  THE FRONT MATTER HEADINGS ARE CENTRED, and only those.  See
    centre_front_matter().

8.   THE QUARTO TITLE BLOCK IS DROPPED.  See strip_title_block().

9.   THE ABSTRACT COMES BEFORE THE THREE LISTS.  See front_matter_order().

10.  EACH FRONT LIST STARTS A NEW PAGE.  See break_before_front_lists().

11.  EVERY FIGURE CARRIES ALT TEXT.  See set_figure_alt_text(), and fig-alt.yml
     beside this script for the descriptions and for why they live there.
     Quarto leaves all seventeen alt-text slots empty in docx, so they are
     filled here, matched to their pictures by figure id rather than by
     position.  A mismatch writes nothing and exits non-zero.

12.  NO FIGURE IS WIDER THAN THE MARGINS ALLOW.  See cap_figure_width().
     UC Davis requires a one-inch margin on every side and says so for figures
     explicitly.  The text column is 6.5 in; a figure written {width=100%}
     resolves to exactly that, which leaves it nothing to spare.  Quarto wraps
     every captioned figure in a one-cell table, and Word lays the cell's
     content box inside the column rather than around it, so a 6.5 in image
     overflows to the right of the text boundary.  Eight of the seventeen are
     affected.  They are scaled to 6.35 in here rather than in the .qmd because
     those eight lines are byte-identical to chapter2.qmd and chapter3.qmd
     under the empty-diff invariant, and the journals have no margin rule to
     satisfy.  See deferred.md, which records the same reasoning for
     reference.docx.

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

# The widest an image may be drawn.  The text column is 6.5 in and the wrapper
# table's cell content box sits inside it, so 6.5 in does not fit; 6.35 in
# leaves room either side and centres like the nine narrower figures.  EMU,
# which is what wp:extent and a:ext are written in: 914400 to the inch.
MAX_FIG_EMU = round(6.35 * 914400)

# Both members of a drawing's size pair -- wp:extent on the inline, a:ext
# inside pic:spPr/a:xfrm -- carry the same cx and cy and must move together.
# Neither wp:effectExtent nor a:off matches this, since they carry l/t/r/b and
# x/y rather than cx/cy.
EXTENT = re.compile(
    r'<(?P<tag>wp:extent|a:ext)\s+cx="(?P<cx>\d+)"\s+cy="(?P<cy>\d+)"\s*/>'
)

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


FRONT_HEADING = re.compile(
    r'<w:p\b[^>]*>\s*<w:pPr>((?:(?!</w:pPr>).)*?)</w:pPr>.*?</w:p>', re.S
)
JC_AT = PPR_ORDER.index("jc")


def centre_front_matter(xml):
    """7: centre the front matter headings, and only those.

    The five are the Abstract, the Acknowledgments and the three list
    headings.  Abstract and Acknowledgments share the Heading1 style with
    Chapter 1, Chapter 2, Chapter 3 and the appendices, which are left-aligned
    and stay that way, so centring the style is not available; direct
    formatting on the paragraph overrides the style for that paragraph alone.

    What separates the five is position, not their titles: they are the
    headings of the lowercase-roman preliminary section, so they are the
    heading paragraphs before the sectPr that section_break() inserts.  That
    is why this runs after section_break() and front_matter_order() rather
    than on a list of titles -- "Abstract" was the only one that could be
    matched by its words, and a second document with a Dedication or a Vita
    would need no change here.

    w:jc comes late in the order CT_PPrBase fixes, so it goes after the last
    child that precedes it and before the first that follows.
    """
    end = xml.find("<w:sectPr")
    if end == -1:
        return xml, 0, 0
    out, pos, centred, found = [], 0, 0, 0
    for m in FRONT_HEADING.finditer(xml, 0, end):
        if m.start() < pos:
            continue
        inner = m.group(1)
        style = re.search(r'<w:pStyle w:val="([^"]+)"', inner)
        if not style or not (style.group(1).startswith("Heading")
                             or style.group(1) == "TOCHeading"):
            continue
        found += 1
        if "<w:jc " in inner:
            continue                            # already centred
        tags = [t for t, _ in _children(inner)]
        after = [t for t in tags
                 if t in PPR_ORDER and PPR_ORDER.index(t) > JC_AT]
        if after:
            at = m.start(1) + inner.find("<w:" + after[0])
        else:
            at = m.start(1) + len(inner)
        out.append(xml[pos:at])
        out.append('<w:jc w:val="center"/>')
        pos = at
        centred += 1
    out.append(xml[pos:])
    return "".join(out), centred, found



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


# ── Figure alt text (11) ────────────────────────────────────────────────────

ALT_PATH = Path(__file__).resolve().parent / "fig-alt.yml"

# `<id>: >` opening a block.  See fig-alt.yml for the grammar this accepts.
ALT_KEY = re.compile(r"^([A-Za-z0-9][A-Za-z0-9_.-]*):[ \t]*>[ \t]*$")

FIG_BOOKMARK = re.compile(r'<w:bookmarkStart\b[^>]*w:name="(fig-[^"]*)"[^>]*>')
DOCPR = re.compile(r"<wp:docPr\b[^>]*?/>")
DESCR_ATTR = re.compile(r'\sdescr="[^"]*"')


def _xml_attr(text):
    """Escape `text` for use inside a double-quoted XML attribute.

    Ampersand first, or the escapes introduced after it would be re-escaped.
    Apostrophes need no escape inside double quotes and are left alone so the
    stored text stays readable.
    """
    return (text.replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;").replace('"', "&quot;"))


def read_alt_text(path=ALT_PATH):
    """Read fig-alt.yml into {figure id: description}.

    Returns (mapping, problem).  A missing file is not a problem -- it means
    the feature is not in use and returns (None, None).  Anything the grammar
    does not cover is, and stops the run before the document is touched: this
    file exists to put words in front of a reader who cannot see the figure,
    and a half-read mapping is worse than none.

    Deliberately hand-parsed rather than handed to PyYAML.  The render workflow
    installs no Python packages and this script has no third-party imports;
    depending on one here would trade a filing-day render against a
    convenience.  The accepted subset is documented in fig-alt.yml.
    """
    if not path.exists():
        return None, None
    alt, key, buf = {}, None, []

    def close():
        if key is not None:
            alt[key] = " ".join(buf).strip()

    for n, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.rstrip()
        if not line.strip():
            close()
            key, buf = None, []
            continue
        if key is None:
            if line.strip() == "---" or line.lstrip().startswith("#"):
                continue
            m = ALT_KEY.match(line)
            if not m:
                return None, f"{path.name}:{n}: expected `<figure-id>: >`, got {line!r}"
            key = m.group(1)
            if key in alt:
                return None, f"{path.name}:{n}: duplicate entry for {key}"
            buf = []
        else:
            if not line.startswith("  "):
                return None, (f"{path.name}:{n}: continuation of {key} must be "
                              f"indented by two spaces, got {line!r}")
            buf.append(line.strip())
    close()
    if not alt:
        return None, f"{path.name}: no entries"
    return alt, None


def set_figure_alt_text(xml, alt):
    """11: give every picture the description its figure id maps to.

    Quarto emits seventeen `wp:docPr` elements with `descr=""`; Word reads that
    attribute as the alt text.  `fig-alt` does nothing here -- it reaches HTML
    only -- and the one source route that does reach `descr`, text in the
    `![...]` brackets, renders as a second visible caption on the pkgdown site.
    Hence after the render, and hence only for the dissertation.

    Matching is by name, never by position.  Quarto writes each figure's id
    into the document as a bookmark just before the picture, so a description
    can be tied to its own figure.  Counting would fail silently and by one,
    and a wrong description is invisible: the page is unchanged and a screen
    reader reads it out with confidence.  So the match has to be a bijection --
    every picture claims its own bookmark, every id in the document has an
    entry, every entry is used -- and anything short of that writes nothing.

    Returns (xml, matched, pictures, problem).
    """
    pics = list(DOCPR.finditer(xml))
    if alt is None:
        return xml, 0, len(pics), None

    marks = [(m.start(), m.group(1)) for m in FIG_BOOKMARK.finditer(xml)]
    names = [n for _, n in marks]
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        return xml, 0, len(pics), "duplicate figure bookmarks: " + ", ".join(dup)

    pairs, claimed = [], set()
    for m in pics:
        before = [n for o, n in marks if o < m.start()]
        if not before:
            return xml, 0, len(pics), "a picture has no figure bookmark before it"
        name = before[-1]
        if name in claimed:
            return xml, 0, len(pics), f"two pictures resolve to {name}"
        claimed.add(name)
        pairs.append((m, name))

    in_doc = {n for _, n in pairs}
    no_entry = sorted(in_doc - set(alt))
    no_figure = sorted(set(alt) - in_doc)
    if no_entry:
        return xml, 0, len(pics), ("figures with no entry in "
                                   f"{ALT_PATH.name}: " + ", ".join(no_entry))
    if no_figure:
        return xml, 0, len(pics), (f"entries in {ALT_PATH.name} matching no "
                                   "figure: " + ", ".join(no_figure))
    for name in sorted(in_doc):
        text = alt[name]
        if not text.strip():
            return xml, 0, len(pics), f"empty alt text for {name}"
        if any(ord(c) < 0x20 or ord(c) == 0x7F for c in text):
            return xml, 0, len(pics), f"control character in the alt text for {name}"

    # Back to front, so the offsets of the matches still ahead stay valid.
    for m, name in reversed(pairs):
        value = _xml_attr(alt[name])
        tag = m.group(0)
        if DESCR_ATTR.search(tag):
            # A lambda, not a replacement string: a backslash in the alt text
            # would otherwise be read as a regex escape.
            tag = DESCR_ATTR.sub(lambda _: f' descr="{value}"', tag, count=1)
        else:
            tag = tag.replace("<wp:docPr", f'<wp:docPr descr="{value}"', 1)
        xml = xml[:m.start()] + tag + xml[m.end():]
    return xml, len(pairs), len(pics), None


def cap_figure_width(xml):
    """Scale every image wider than MAX_FIG_EMU down to it, keeping its ratio.

    Returns the xml, how many size elements were changed, and the widths in
    inches before and after, so the caller can show what moved.  Idempotent:
    a second pass finds nothing above the cap.
    """
    changed = []

    def shrink(m):
        cx, cy = int(m.group("cx")), int(m.group("cy"))
        if cx <= MAX_FIG_EMU:
            return m.group(0)
        new_cy = round(cy * MAX_FIG_EMU / cx)
        changed.append((m.group("tag"), cx, cy, new_cy))
        return (f'<{m.group("tag")} cx="{MAX_FIG_EMU}" cy="{new_cy}" />')

    out = EXTENT.sub(shrink, xml)

    # Every drawing has exactly one wp:extent and one a:ext, so the two counts
    # must match; if they do not, the pair has been split and the picture would
    # disagree with its frame.
    wp = sum(1 for t, *_ in changed if t == "wp:extent")
    ax = sum(1 for t, *_ in changed if t == "a:ext")
    if wp != ax:
        return xml, 0, [], f"extent pair mismatch: {wp} wp:extent, {ax} a:ext"
    sizes = [(cx / 914400, cy / 914400, MAX_FIG_EMU / 914400, ncy / 914400)
             for t, cx, cy, ncy in changed if t == "wp:extent"]
    return out, len(changed), sizes, None


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    src = Path(sys.argv[1])
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else src
    if not src.exists():
        sys.exit(f"not found: {src}")

    alt_map, alt_problem = read_alt_text()
    if alt_problem:
        sys.exit(f"fix_submission: {alt_problem}\n"
                 f"  no alt text written; {dst} left unchanged.")

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
                xml, sectioned = section_break(xml)
                xml, reordered = front_matter_order(xml)
                xml, broken = break_before_front_lists(xml)
                xml, centred, found = centre_front_matter(xml)
                xml, hoisted, thinned = no_blank_pages(xml)
                xml, alt_ok, alt_n, alt_problem = set_figure_alt_text(
                    xml, alt_map)
                xml, capped, cap_sizes, cap_problem = cap_figure_width(xml)
                if cap_problem:
                    sys.exit(f"fix_submission: {cap_problem}")
                print(f"  caption pPr merged        : {merged}")
                print(f"  table captions restyled   : {restyled}")
                print(f"  list field codes repaired : {repaired}")
                print(f"  front lists inserted      : {inserted}")
                print(f"  title-block paras dropped : {dropped}")
                print(f"  front matter reordered    : {reordered}")
                print(f"  front lists on a new page : {broken}")
                print(f"  front matter centred      : {centred} of {found}")
                print(f"  section break inserted    : {sectioned}")
                print(f"  breaks -> pageBreakBefore  : {hoisted}")
                print(f"  breaks kept but thinned   : {thinned}")
                print(f"  figure alt text           : {alt_ok} of {alt_n} matched by id")
                print(f"  figures capped to 6.35 in : {len(cap_sizes)}"
                      f" ({capped} size elements)")
                for w, h, nw, nh in cap_sizes:
                    print(f"      {w:.3f} x {h:.3f} -> {nw:.3f} x {nh:.3f} in")
                data = xml.encode("utf-8")
            elif item.filename == "word/settings.xml":
                xml = data.decode("utf-8")
                xml, added = patch_settings(xml)
                print(f"  updateFields added        : {added}")
                data = xml.encode("utf-8")
            zout.writestr(item, data)

    if alt_problem:
        tmp.unlink()
        sys.exit(f"fix_submission: {alt_problem}\n"
                 f"  no alt text written; {dst} left unchanged.")

    shutil.move(str(tmp), str(dst))
    print(f"  written: {dst}")


if __name__ == "__main__":
    main()
