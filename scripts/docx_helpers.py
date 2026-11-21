"""python-docx helpers for the course paper.

Implements the common requirements of a Fuzhou University undergraduate thesis:
A4 paper, generous margins, SimSun body text at 小四 (12 pt) with 1.5 line
spacing, headings in 黑体, table captions above and figure captions below,
automatic table of contents field, page numbers in the footer.
"""
from __future__ import annotations

from pathlib import Path

import docx
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

BODY_FONT = "宋体"
HEAD_FONT = "黑体"
LATIN_FONT = "Times New Roman"


def set_run_font(run, name_cn: str = BODY_FONT, name_en: str = LATIN_FONT, size: float | None = None,
                 bold: bool | None = None, color: RGBColor | None = None) -> None:
    run.font.name = name_en
    run._element.rPr.rFonts.set(qn("w:eastAsia"), name_cn)
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.font.bold = bold
    if color is not None:
        run.font.color.rgb = color


def _set_style_font(style, name_cn: str, name_en: str, size: float | None = None, bold: bool | None = None):
    style.font.name = name_en
    style.element.rPr.rFonts.set(qn("w:eastAsia"), name_cn)
    if size is not None:
        style.font.size = Pt(size)
    if bold is not None:
        style.font.bold = bold


def _add_field(paragraph, instr: str) -> None:
    run = paragraph.add_run()
    fld_begin = OxmlElement("w:fldChar")
    fld_begin.set(qn("w:fldCharType"), "begin")
    instr_el = OxmlElement("w:instrText")
    instr_el.set(qn("xml:space"), "preserve")
    instr_el.text = instr
    fld_sep = OxmlElement("w:fldChar")
    fld_sep.set(qn("w:fldCharType"), "separate")
    text = OxmlElement("w:t")
    text.text = "右键更新域"
    fld_end = OxmlElement("w:fldChar")
    fld_end.set(qn("w:fldCharType"), "end")
    run._r.append(fld_begin)
    run._r.append(instr_el)
    run._r.append(fld_sep)
    run._r.append(text)
    run._r.append(fld_end)


def _page_number_footer(section) -> None:
    footer = section.footer
    p = footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run("— ")
    set_run_font(run, size=9)
    _add_field(p, "PAGE")
    run2 = p.add_run(" —")
    set_run_font(run2, size=9)


class Paper:
    """Thin wrapper around ``python-docx`` with thesis-ready defaults."""

    def __init__(self, title: str = ""):
        self.doc = Document()
        self.title = title
        self.fig_no = 0
        self.tab_no = 0
        self.heading_no = 0
        self._setup_page()
        self._setup_styles()

    # ------------------------------------------------------------------ setup
    def _setup_page(self) -> None:
        s = self.doc.sections[0]
        s.page_width, s.page_height = Cm(21.0), Cm(29.7)
        s.top_margin = Cm(2.8)
        s.bottom_margin = Cm(2.5)
        s.left_margin = Cm(3.0)
        s.right_margin = Cm(2.5)
        s.header_distance = Cm(1.5)
        s.footer_distance = Cm(1.6)
        _page_number_footer(s)
        # header text
        hp = s.header.paragraphs[0]
        hp.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = hp.add_run("福州大学本科生课程论文")
        set_run_font(r, name_cn=HEAD_FONT, size=9)

    def _setup_styles(self) -> None:
        st = self.doc.styles["Normal"]
        _set_style_font(st, BODY_FONT, LATIN_FONT, size=12)
        pf = st.paragraph_format
        pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
        pf.space_before = Pt(0)
        pf.space_after = Pt(0)
        pf.first_line_indent = Pt(24)
        pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        # body style used explicitly by paragraphs
        if "Body CN" not in [s.name for s in self.doc.styles]:
            body = self.doc.styles.add_style("Body CN", docx.enum.style.WD_STYLE_TYPE.PARAGRAPH)
            body.base_style = st
            _set_style_font(body, BODY_FONT, LATIN_FONT, size=12)
            body.paragraph_format.first_line_indent = Pt(24)
            body.paragraph_format.line_spacing = 1.5
            body.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        for lvl, size, before, after in ((1, 15, 18, 12), (2, 13, 12, 6), (3, 12, 6, 6)):
            name = f"Heading {lvl}"
            try:
                h = self.doc.styles[name]
            except KeyError:
                continue
            _set_style_font(h, HEAD_FONT, "Arial", size=size, bold=True)
            h.font.color.rgb = RGBColor(0, 0, 0)
            h.paragraph_format.space_before = Pt(before)
            h.paragraph_format.space_after = Pt(after)
            h.paragraph_format.first_line_indent = Pt(0)
            h.paragraph_format.line_spacing = 1.5

    # ---------------------------------------------------------------- content
    def para(self, text: str, indent: bool = True, bold_lead: str | None = None,
             size: float | None = 12, align=WD_ALIGN_PARAGRAPH.JUSTIFY,
             space_after: float = 0, first_indent: float | None = 24) -> None:
        p = self.doc.add_paragraph(style="Body CN")
        p.alignment = align
        p.paragraph_format.first_line_indent = Pt(first_indent if indent else 0)
        p.paragraph_format.space_after = Pt(space_after)
        if bold_lead:
            r = p.add_run(bold_lead)
            set_run_font(r, name_cn=HEAD_FONT, size=size)
        r = p.add_run(text)
        set_run_font(r, size=size)
        return p

    def heading(self, text: str, level: int = 1, numbered: bool = True) -> None:
        self.heading_no += 1
        p = self.doc.add_heading(level=level)
        p.paragraph_format.first_line_indent = Pt(0)
        r = p.add_run(text)
        set_run_font(r, name_cn=HEAD_FONT, name_en="Arial", size={1: 15, 2: 13, 3: 12}.get(level, 12), bold=True,
                     color=RGBColor(0, 0, 0))
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT

    def caption(self, text: str, kind: str = "figure", number: int | None = None) -> None:
        if kind == "figure":
            self.fig_no += 1
            label = f"图 {self.fig_no}  {text}"
            size = 10.5
        else:
            self.tab_no += 1
            n = number if number is not None else self.tab_no
            label = f"表 {n}  {text}"
            size = 10.5
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.first_line_indent = Pt(0)
        p.paragraph_format.space_before = Pt(3)
        p.paragraph_format.space_after = Pt(6)
        r = p.add_run(label)
        set_run_font(r, size=size, bold=(kind == "table"))

    def table(self, header: list[str], rows: list[list[str]], caption: str,
              widths: list[float] | None = None, font_size: float = 9,
              align_right: set[int] | None = None, number: int | None = None) -> None:
        self.caption(caption, kind="table", number=number)
        t = self.doc.add_table(rows=1, cols=len(header))
        t.style = "Table Grid"
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        t.autofit = False
        hdr = t.rows[0].cells
        for i, h in enumerate(header):
            hdr[i].text = ""
            p = hdr[i].paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.first_line_indent = Pt(0)
            r = p.add_run(h)
            set_run_font(r, name_cn=HEAD_FONT, size=font_size, bold=True)
        for row in rows:
            cells = t.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = ""
                p = cells[i].paragraphs[0]
                p.paragraph_format.first_line_indent = Pt(0)
                p.alignment = WD_ALIGN_PARAGRAPH.RIGHT if (align_right and i in align_right) else WD_ALIGN_PARAGRAPH.CENTER
                r = p.add_run(str(v))
                set_run_font(r, size=font_size)
        if widths:
            for i, w in enumerate(widths):
                for row in t.rows:
                    row.cells[i].width = Cm(w)
        # small gap after table
        sp = self.doc.add_paragraph()
        sp.paragraph_format.space_after = Pt(0)
        sp.paragraph_format.first_line_indent = Pt(0)

    def figure(self, path: str | Path, caption: str, width_cm: float = 14.0) -> None:
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.first_line_indent = Pt(0)
        p.paragraph_format.space_before = Pt(6)
        run = p.add_run()
        run.add_picture(str(path), width=Cm(width_cm))
        self.caption(caption, kind="figure")

    def page_break(self) -> None:
        self.doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)

    def toc(self) -> None:
        p = self.doc.add_paragraph()
        p.paragraph_format.first_line_indent = Pt(0)
        _add_field(p, r'TOC \o "1-3" \h \z \u')

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.doc.save(path)
        return path
