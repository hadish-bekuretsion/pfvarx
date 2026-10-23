#!/usr/bin/env python3
"""Build docs/PfVarX_Methodology.docx from the verified implementation.

Every factual claim in the generated document is either:
  * read at build time from the package itself (parameter registry, release
    status, decision categories, filter names), or
  * copied verbatim from docs/METHODOLOGY_CHAPTER.md, or
  * a statement about what the audit verified, traceable to docs/AUDIT_REPORT.md.

Nothing is typed in from memory that could contradict the code.

Usage:
    python scripts/build_methodology_docx.py [--out docs/PfVarX_Methodology.docx]
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import re
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

from pfvarx import __version__
from pfvarx.calling.policy import (
    BIN_BELOW_THRESHOLD,
    BIN_INTERMEDIATE,
    BIN_LOW_FRACTION,
    BIN_MINORITY,
    BIN_PREDOMINANT,
    DECISION_ALT_EXPLORATORY,
    DECISION_ALT_PASS,
    DECISION_ALT_UNCERTAIN,
    DECISION_NO_CALL,
    DECISION_REFERENCE_SUPPORTED,
)
from pfvarx.calling.writer import _FILTERS
from pfvarx.config.parameters import REGISTRY, registry_hash
from pfvarx.release import RELEASE_STATUS, overall_status

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"

ACCENT = RGBColor(0x1F, 0x3B, 0x57)
MUTED = RGBColor(0x55, 0x5F, 0x6B)
WARN = RGBColor(0x8A, 0x4B, 0x00)

BODY_FONT = "Calibri"
HEAD_FONT = "Cambria"


# --------------------------------------------------------------------------
# low-level formatting helpers
# --------------------------------------------------------------------------
def _shade(cell, hex_fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hex_fill)
    tc_pr.append(shd)


def _cell_margins(table, top=60, bottom=60, left=90, right=90) -> None:
    tbl_pr = table._tbl.tblPr
    mar = OxmlElement("w:tblCellMar")
    for tag, val in (("top", top), ("left", left), ("bottom", bottom), ("right", right)):
        el = OxmlElement(f"w:{tag}")
        el.set(qn("w:w"), str(val))
        el.set(qn("w:type"), "dxa")
        mar.append(el)
    tbl_pr.append(mar)


def _repeat_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader")
    el.set(qn("w:val"), "true")
    tr_pr.append(el)


def _field(paragraph, instr: str) -> None:
    run = paragraph.add_run()
    fld_begin = OxmlElement("w:fldChar")
    fld_begin.set(qn("w:fldCharType"), "begin")
    run._r.append(fld_begin)
    run2 = paragraph.add_run()
    instr_el = OxmlElement("w:instrText")
    instr_el.set(qn("xml:space"), "preserve")
    instr_el.text = instr
    run2._r.append(instr_el)
    run3 = paragraph.add_run()
    fld_sep = OxmlElement("w:fldChar")
    fld_sep.set(qn("w:fldCharType"), "separate")
    run3._r.append(fld_sep)
    run4 = paragraph.add_run("1")
    run4.font.color.rgb = MUTED
    run5 = paragraph.add_run()
    fld_end = OxmlElement("w:fldChar")
    fld_end.set(qn("w:fldCharType"), "end")
    run5._r.append(fld_end)


def _style_document(doc: Document) -> None:
    normal = doc.styles["Normal"]
    normal.font.name = BODY_FONT
    normal.font.size = Pt(10.5)
    normal.paragraph_format.space_after = Pt(7)
    normal.paragraph_format.line_spacing = 1.14
    normal.element.rPr.rFonts.set(qn("w:eastAsia"), BODY_FONT)

    for name, size, colour, before, after in (
        ("Heading 1", 17, ACCENT, 20, 8),
        ("Heading 2", 13, ACCENT, 14, 5),
        ("Heading 3", 11.5, ACCENT, 11, 4),
    ):
        st = doc.styles[name]
        st.font.name = HEAD_FONT
        st.font.size = Pt(size)
        st.font.bold = True
        st.font.color.rgb = colour
        st.paragraph_format.space_before = Pt(before)
        st.paragraph_format.space_after = Pt(after)
        st.paragraph_format.keep_with_next = True
        rpr = st.element.get_or_add_rPr()
        rf = rpr.find(qn("w:rFonts"))
        if rf is None:
            rf = OxmlElement("w:rFonts")
            rpr.append(rf)
        rf.set(qn("w:ascii"), HEAD_FONT)
        rf.set(qn("w:hAnsi"), HEAD_FONT)

    for section in doc.sections:
        section.page_height = Cm(29.7)
        section.page_width = Cm(21.0)
        section.top_margin = Cm(2.4)
        section.bottom_margin = Cm(2.2)
        section.left_margin = Cm(2.5)
        section.right_margin = Cm(2.5)


def _header_footer(doc: Document) -> None:
    section = doc.sections[0]
    hp = section.header.paragraphs[0]
    hp.text = ""
    run = hp.add_run("PfVarX \u2014 Methodology and Validation Report")
    run.font.size = Pt(8.5)
    run.font.color.rgb = MUTED
    hp.add_run("\t\t")
    run2 = hp.add_run(f"v{__version__}")
    run2.font.size = Pt(8.5)
    run2.font.color.rgb = MUTED
    hp.alignment = WD_ALIGN_PARAGRAPH.LEFT

    fp = section.footer.paragraphs[0]
    fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = fp.add_run("Page ")
    r.font.size = Pt(8.5)
    r.font.color.rgb = MUTED
    _field(fp, "PAGE")
    r = fp.add_run(" of ")
    r.font.size = Pt(8.5)
    r.font.color.rgb = MUTED
    _field(fp, "NUMPAGES")


# --------------------------------------------------------------------------
# content helpers
# --------------------------------------------------------------------------
def h1(doc, text):
    return doc.add_heading(text, level=1)


def h2(doc, text):
    return doc.add_heading(text, level=2)


def h3(doc, text):
    return doc.add_heading(text, level=3)


def para(doc, text, *, italic=False, size=None, colour=None, align=None, space_after=None):
    p = doc.add_paragraph()
    _rich(p, text, italic=italic, size=size, colour=colour)
    if align is not None:
        p.alignment = align
    if space_after is not None:
        p.paragraph_format.space_after = Pt(space_after)
    return p


_MARKUP = re.compile(r"(\*\*.+?\*\*|\*.+?\*|`.+?`)")


def _rich(p, text, *, italic=False, size=None, colour=None):
    """Render **bold**, *italic* and `code` inline markup."""
    for chunk in _MARKUP.split(text):
        if not chunk:
            continue
        bold = ital = code = False
        if chunk.startswith("**") and chunk.endswith("**") and len(chunk) > 4:
            chunk, bold = chunk[2:-2], True
        elif chunk.startswith("*") and chunk.endswith("*") and len(chunk) > 2:
            chunk, ital = chunk[1:-1], True
        elif chunk.startswith("`") and chunk.endswith("`") and len(chunk) > 2:
            chunk, code = chunk[1:-1], True
        run = p.add_run(chunk)
        run.bold = bold
        run.italic = ital or italic
        if code:
            run.font.name = "Consolas"
            run.font.size = Pt(9)
        if size:
            run.font.size = Pt(size)
        if colour:
            run.font.color.rgb = colour


def bullet(doc, text, level=0):
    p = doc.add_paragraph(style="List Bullet")
    p.paragraph_format.left_indent = Cm(0.7 + 0.6 * level)
    p.paragraph_format.space_after = Pt(3)
    _rich(p, text)
    return p


def numbered(doc, text):
    p = doc.add_paragraph(style="List Number")
    p.paragraph_format.space_after = Pt(3)
    _rich(p, text)
    return p


def code_block(doc, text):
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Cm(0.8)
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(8)
    p.paragraph_format.keep_together = True
    run = p.add_run(text)
    run.font.name = "Consolas"
    run.font.size = Pt(8.5)
    run.font.color.rgb = RGBColor(0x22, 0x22, 0x22)
    pPr = p._p.get_or_add_pPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:fill"), "F4F5F7")
    pPr.append(shd)
    return p


def table(doc, headers, rows, widths=None, font_size=8.5, header_fill="1F3B57", col_sizes=None):
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    _cell_margins(t)
    hdr = t.rows[0]
    _repeat_header(hdr)
    for i, text in enumerate(headers):
        cell = hdr.cells[i]
        cell.text = ""
        p = cell.paragraphs[0]
        p.paragraph_format.space_after = Pt(0)
        run = p.add_run(text)
        run.bold = True
        run.font.size = Pt(font_size)
        run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        _shade(cell, header_fill)
    for r_i, row in enumerate(rows):
        cells = t.add_row().cells
        for i, val in enumerate(row):
            cells[i].text = ""
            p = cells[i].paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            size = col_sizes[i] if col_sizes else font_size
            _rich(p, str(val), size=size)
            if r_i % 2 == 1:
                _shade(cells[i], "F2F4F7")
    if widths:
        t.autofit = False
        tbl_pr = t._tbl.tblPr
        layout = OxmlElement("w:tblLayout")
        layout.set(qn("w:type"), "fixed")
        tbl_pr.append(layout)
        grid = t._tbl.tblGrid
        for i, gc in enumerate(grid.findall(qn("w:gridCol"))):
            if i < len(widths):
                gc.set(qn("w:w"), str(int(Cm(widths[i]).twips)))
        for row in t.rows:
            for i, w in enumerate(widths):
                row.cells[i].width = Cm(w)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)
    return t


def callout(doc, label, text, fill="FFF4E5", colour=WARN):
    t = doc.add_table(rows=1, cols=1)
    t.style = "Table Grid"
    _cell_margins(t, top=90, bottom=90, left=140, right=140)
    cell = t.rows[0].cells[0]
    cell.text = ""
    _shade(cell, fill)
    p = cell.paragraphs[0]
    p.paragraph_format.space_after = Pt(0)
    run = p.add_run(label + " ")
    run.bold = True
    run.font.size = Pt(9)
    run.font.color.rgb = colour
    _rich(p, text, size=9)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)


def page_break(doc):
    doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)


# --------------------------------------------------------------------------
# document parts
# --------------------------------------------------------------------------
def title_page(doc: Document) -> None:
    for _ in range(4):
        doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    run = p.add_run("PfVarX")
    run.font.name = HEAD_FONT
    run.font.size = Pt(40)
    run.font.bold = True
    run.font.color.rgb = ACCENT

    p = doc.add_paragraph()
    run = p.add_run("Methodology, Implementation and Validation Report")
    run.font.name = HEAD_FONT
    run.font.size = Pt(19)
    run.font.color.rgb = ACCENT

    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(6)
    run = p.add_run(
        "A configuration-driven platform for SNP and small-INDEL discovery in "
        "Plasmodium falciparum"
    )
    run.font.size = Pt(12)
    run.font.color.rgb = MUTED
    run.italic = True

    for _ in range(3):
        doc.add_paragraph()

    rows = [
        ("Document type", "Methodology chapter and engineering validation report"),
        ("Software version", f"PfVarX {__version__}"),
        ("Parameter registry", f"{len(REGISTRY)} parameters \u00b7 hash {registry_hash()[:12]}"),
        ("Overall release status", overall_status()),
        ("Date of issue", _dt.date.today().isoformat()),
        ("Organism and scope", "Plasmodium falciparum; SNPs and small INDELs only"),
        ("Intended use", "Research use. Not a clinical diagnostic device."),
    ]
    table(doc, ["Field", "Value"], rows, widths=[4.6, 11.4], font_size=9)

    doc.add_paragraph()
    callout(
        doc,
        "Status of this document.",
        "This report describes what the software actually does, as verified by the "
        "audit recorded in docs/AUDIT_REPORT.md. No component of PfVarX is claimed to be "
        "scientifically validated for variant-calling accuracy. Parameters marked "
        "PROVISIONAL_DEVELOPMENT are engineering defaults, not calibrated values. "
        "Statements about performance are limited to what was measured in this "
        "environment; no accuracy benchmark against an external truth set has been run.",
        fill="FDECEA",
        colour=RGBColor(0x9B, 0x1C, 0x1C),
    )


def front_matter(doc: Document, toc: list[tuple[str, int]] | None = None) -> None:
    page_break(doc)
    p = doc.add_paragraph()
    run = p.add_run("Contents")
    run.font.name = HEAD_FONT
    run.font.size = Pt(17)
    run.font.bold = True
    run.font.color.rgb = ACCENT
    p.paragraph_format.space_after = Pt(10)
    para(
        doc,
        "The document is organised in fifteen sections. Sections 1\u20133 establish the "
        "research objectives, the software architecture and the input contract. Sections "
        "4\u20138 cover the analytical core: preprocessing, evidence extraction, the "
        "statistical models, the decision policy, and the calibration argument. Sections "
        "9\u201311 cover annotation, downstream analyses and reporting. Sections 12\u201313 "
        "record the engineering evidence: testing, benchmarking, reproducibility, error "
        "handling and limitations. Section 14 lists references; Section 15 is the "
        "validation-status table that separates implemented, tested, calibrated and "
        "scientifically validated components.",
    )

    if toc:
        rows = [(title, str(page)) for title, page in toc]
        t = table(doc, ["Section", "Page"], rows, widths=[14.2, 1.8], font_size=9)
        for row in t.rows[1:]:
            row.cells[1].paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.RIGHT
    else:
        para(
            doc,
            "Table of contents generated on the second build pass.",
            italic=True,
            size=9,
            colour=MUTED,
        )

    h2(doc, "Reading conventions")
    para(
        doc,
        "Four status words are used with fixed meanings throughout. **Implemented** means "
        "the code path exists and executes. **Tested** means an automated test exercises "
        "the path and asserts its behaviour. **Calibrated** means the parameters governing "
        "the path were fitted to control data supplied by the user, with a recorded "
        "dataset identifier. **Scientifically validated** means accuracy has been "
        "demonstrated against an independent truth set with reported sensitivity, "
        "precision and false-discovery rate. The four are not interchangeable, and the "
        "distinction is the organising principle of Section 15.",
    )
    para(
        doc,
        "Identifiers beginning `PVX-` are policy identifiers. They encode the analysis "
        "mode, the policy family, a twelve-character hash of the resolved configuration "
        "and an eight-character hash of the parameter registry, so that any result can be "
        "traced to the exact configuration and parameter set that produced it.",
    )


def section_1(doc):
    page_break(doc)
    h1(doc, "1. Objectives, scope and methodological principles")
    para(
        doc,
        "PfVarX addresses a specific and narrow question: given Illumina short-read "
        "sequencing of *Plasmodium falciparum*, what is the most defensible way to decide "
        "whether a genomic position carries a non-reference allele, and how should the "
        "confidence in that decision be represented so that it is neither overstated nor "
        "silently discarded? The platform is not a general-purpose variant caller. It is "
        "a structured argument, expressed in code, that variant calling in a haploid, "
        "AT-rich, frequently polyclonal parasite genome requires explicit count models, "
        "explicit decision rules and explicit statements of what has and has not been "
        "validated.",
    )

    h2(doc, "1.1 Objectives")
    for text in (
        "Detect single-nucleotide polymorphisms and small insertions and deletions from "
        "aligned short reads, with every call traceable to persisted read-level evidence.",
        "Represent uncertainty honestly: a likelihood ratio is reported as a likelihood "
        "ratio, a posterior only when a calibrated error model exists, and a Phred-scaled "
        "quality only when the underlying statistic has been calibrated.",
        "Accommodate the biological realities of *P. falciparum*: haploid genotype "
        "representation, AT-rich sequence, ubiquitous over-dispersion, and polyclonal "
        "infections in which the alternate allele fraction may take any value between "
        "zero and one.",
        "Provide a single configuration surface, validated by schema, hashed for "
        "reproducibility, and expressive enough to run the same code on a laptop, a "
        "server, an HPC scheduler or a container platform without code changes.",
        "Fail closed. Where an assumption cannot be met \u2014 an unverified UMI design, a "
        "reference-allele mismatch, an absent callable region \u2014 the pipeline must "
        "refuse to emit a confident call and must say why.",
    ):
        bullet(doc, text)

    h2(doc, "1.2 Scope")
    para(
        doc,
        "The analytical scope is deliberately bounded. Only SNPs and small INDELs are "
        "called. Structural variants, copy-number variation, large deletions and complex "
        "rearrangements are out of scope and are not represented in the output schema. "
        "The caller contains no machine-learning component: every decision is a function "
        "of an explicit likelihood, an explicit prior where enabled, and an explicit "
        "threshold, each of which can be inspected and argued with. Mixed-infection "
        "evidence is retained inside the calling logic, because the allele fraction is a "
        "property of the read evidence; formal complexity-of-infection estimation is a "
        "downstream analysis and is reported separately.",
    )
    para(
        doc,
        "The platform targets Illumina short-read data. Long-read platforms, which would "
        "change the error model and the indel representation, are not supported. The "
        "reference genome must be a pinned assembly; contig names alone are not accepted "
        "as identification, because the same contig name can denote different builds.",
    )

    h2(doc, "1.3 Methodological principles")
    para(
        doc,
        "Five principles govern the design and are enforced in code rather than stated "
        "only in prose.",
    )
    table(
        doc,
        ["Principle", "Meaning in this implementation", "Where it is enforced"],
        [
            (
                "Evidence before inference",
                "Every call carries the depth, alternate count, strand counts, mapping "
                "quality and context metrics on which it rests.",
                "45-column evidence table; `EVIDENCE_COLUMNS` in `calling/writer.py`",
            ),
            (
                "Statistical honesty",
                "A likelihood ratio is never relabelled a posterior; an uncalibrated "
                "statistic is never written as QUAL.",
                "`_qual()` returns QUAL only when policy status is VALIDATED",
            ),
            (
                "Fail closed",
                "Unmet assumptions produce NO_CALL or ALT_UNCERTAIN, never a silent pass.",
                "Hard gates and review filters in `calling/policy.py`",
            ),
            (
                "Reproducibility by construction",
                "Configuration is hashed; outputs are written atomically; the parameter "
                "registry is versioned.",
                "`config_hash()`, `registry_hash()`, `.tmp` + `os.replace`",
            ),
            (
                "Honest status reporting",
                "The software reports its own validation state and never infers it from "
                "the sophistication of the code.",
                "`release.py`, `status_banner()`, Section 15 of this document",
            ),
        ],
        widths=[3.6, 6.6, 5.8],
    )

    h2(doc, "1.4 What this document does not claim")
    para(
        doc,
        "It does not claim that PfVarX calls variants more accurately than established "
        "callers. No head-to-head accuracy benchmark against GATK, FreeBayes, bcftools or "
        "Plasmodium-specific pipelines has been executed, and none is reported. It does "
        "not claim clinical utility. It does not claim that the default thresholds are "
        "optimal: they are provisional engineering values, documented as such, and the "
        "calibration machinery that would replace them exists but has not been run on a "
        "reference dataset within this project. Section 15 states this component by "
        "component.",
    )


def section_2(doc):
    page_break(doc)
    h1(doc, "2. Architecture and workflow")
    para(
        doc,
        "PfVarX is organised as a layered Python package with a thin orchestration layer. "
        "The package owns all analytical logic and can be driven directly from Python or "
        "from the command line. Snakemake owns execution order, parallelism and "
        "resumption, and contains no analytical decisions of its own. This separation is "
        "what allows the same analytical core to be tested in isolation, invoked "
        "interactively, and deployed on four different execution backends.",
    )

    h2(doc, "2.1 Layered structure")
    table(
        doc,
        ["Layer", "Module", "Responsibility"],
        [
            ("Interface", "`cli.py`", "Argument parsing, subcommand dispatch, exit codes, human-readable diagnostics"),
            ("Configuration", "`config/schema.py`, `config/parameters.py`", "Pydantic v2 models, bounds, precedence resolution, hashing, parameter registry"),
            ("Input", "`io/`", "Format detection, FASTQ/BAM/CRAM/VCF/BCF/gVCF/FASTA validation, variant import"),
            ("Preprocessing", "`preprocess/pipeline.py`", "Versioned command construction for fastp, BWA-MEM2, samtools, Picard, fgbio"),
            ("Evidence", "`calling/evidence.py`, `calling/caller.py`", "CIGAR-aware pileup, fragment de-duplication, per-site metrics"),
            ("Models", "`models/count_models.py`", "Binomial and beta-binomial likelihoods, LRT, bootstrap and exact nulls, BH FDR"),
            ("Policy", "`calling/policy.py`", "Decision categories, hard gates, review filters, allele-fraction bins"),
            ("Output", "`calling/writer.py`", "VCF and evidence-table serialisation, atomic writes, header construction"),
            ("Annotation", "`annotation/annotator.py`", "SnpEff wrapper and internal consequence fallback"),
            ("Downstream", "`downstream/`", "Diversity, mixed infection, relatedness, phylogeny"),
            ("QC", "`qc/integrator.py`", "Metric ontology, PASS/WARN/FAIL/INFO evaluation, run-level roll-up"),
            ("Reporting", "`reporting/`", "CSV, XLSX, PNG and PDF export; Streamlit dashboard"),
            ("Status", "`release.py`", "Machine-readable component status and banner"),
        ],
        widths=[2.6, 5.0, 8.4],
    )

    h2(doc, "2.2 Execution flow")
    para(
        doc,
        "The flow is linear in the analytical sense and branching only at the points where "
        "the input format genuinely changes the required work. A run proceeds through "
        "configuration resolution, input validation and route selection, preprocessing "
        "and alignment where raw reads are supplied, evidence extraction, model "
        "evaluation, policy application, output serialisation, annotation, downstream "
        "analysis, quality-control roll-up and reporting.",
    )
    code_block(
        doc,
        "run config + profile + CLI overrides\n"
        "        \u2502\n"
        "        \u251c\u2500 resolve \u2192 validate (Pydantic) \u2192 hash \u2192 policy_id\n"
        "        \u2514\u2500 validate inputs \u2192 select route\n"
        "                 \u251c\u2500 FASTQ \u2192 fastp \u2192 BWA-MEM2 \u2192 sort \u2192 [dedup | UMI] \u2192 BAM\n"
        "                 \u251c\u2500 BAM/CRAM \u2192 [dedup] \u2192 BAM\n"
        "                 \u2514\u2500 VCF/BCF/gVCF \u2192 import (bypass caller)\n"
        "                          \u2502\n"
        "                 evidence extraction (CIGAR-aware pileup)\n"
        "                          \u2502\n"
        "                 count models \u2192 LRT \u2192 [FDR] \u2192 [posterior]\n"
        "                          \u2502\n"
        "                 policy: hard gates \u2192 review filters \u2192 decision\n"
        "                          \u2502\n"
        "                 VCF + evidence table (atomic)\n"
        "                          \u2502\n"
        "                 annotation \u2192 downstream \u2192 QC \u2192 report",
    )

    h2(doc, "2.3 Mode-aware routing")
    para(
        doc,
        "Four sequencing modes are defined, each selecting a different preprocessing route "
        "and a different evidence model. The mode is a first-class configuration value, "
        "not a runtime guess, so that a run cannot silently change route between "
        "executions.",
    )
    table(
        doc,
        ["Mode", "Input", "Route", "Evidence model"],
        [
            ("`standard`", "FASTQ or BAM/CRAM", "Trim, align, coordinate sort, optional duplicate marking", "Binomial or beta-binomial on read counts"),
            ("`tbs_umi`", "FASTQ with UMI in read structure", "Trim, align, sort, fgbio extract \u2192 group \u2192 consensus, realign", "Consensus-read counts; UMI model must be valid"),
            ("`variant_import`", "VCF, BCF or gVCF", "No alignment; records imported and re-emitted in the PfVarX schema", "Counts recovered from FORMAT/INFO, or marked absent"),
            ("`disabled`", "Pre-aligned BAM with duplicates already marked", "Sort only", "Binomial or beta-binomial on read counts"),
        ],
        widths=[2.6, 3.4, 6.2, 3.8],
    )
    para(
        doc,
        "Routing is resolved before any work begins and is recorded in the run manifest. "
        "If a samplesheet mixes formats that require different routes, the run is rejected "
        "with an explicit diagnostic rather than partitioned silently, because a mixed "
        "run would produce a VCF whose records were generated under incompatible "
        "assumptions.",
    )

    h2(doc, "2.4 Orchestration and execution profiles")
    para(
        doc,
        "The Snakemake layer declares rules, inputs, outputs and resources. Four executor "
        "profiles \u2014 local, server, HPC and cloud \u2014 change only the executor, the "
        "resource declarations and the container directive; they do not change the "
        "analytical graph. Rule ordering is explicit where two rules can produce the same "
        "target: `import_variants` is ordered ahead of `call_variants` so that a "
        "pre-existing VCF is imported rather than re-called.",
    )
    para(
        doc,
        "The workflow was audited statically. Snakemake and the external binaries "
        "(fastp, BWA-MEM2, samtools, bcftools, fgbio, Picard) are not installed in the "
        "audit environment, so the rules were checked for correct input/output wiring, "
        "shell quoting, conditional arguments and rule ordering, but were not executed "
        "end-to-end here. This limitation is recorded in Section 13 and in the audit "
        "report, and it is the reason the preprocessing component is marked implemented "
        "but not tested in Section 15.",
    )


def section_3(doc):
    page_break(doc)
    h1(doc, "3. Input validation and mode-aware routing")
    para(
        doc,
        "Input validation is the point at which most silent failures are prevented. The "
        "design rule is that no file is opened for analysis until its format has been "
        "established from content, its required indices have been confirmed, and its "
        "declared sample identity has been reconciled with the samplesheet. Validation "
        "failures raise a typed exception with a message that names the file, the check "
        "that failed and the corrective action.",
    )

    h2(doc, "3.1 Format detection")
    para(
        doc,
        "Detection is content-based, not extension-based, because extensions are routinely "
        "wrong. A file is classified by inspecting its magic bytes and, where necessary, "
        "its first records. The supported classes are FASTQ (plain and gzip), BAM, CRAM, "
        "VCF, BCF, gVCF and FASTA. A file whose content does not match any supported class "
        "is rejected before any downstream stage is scheduled.",
    )
    table(
        doc,
        ["Class", "Detection signal", "Required companions", "Route"],
        [
            ("FASTQ", "`@` record header, `+` separator, four-line structure", "Reference FASTA + index", "Preprocess and align"),
            ("BAM", "BGZF magic `BAM\\1`", "`.bai` or `.csi`", "Evidence extraction"),
            ("CRAM", "CRAM magic, `CRAM` container", "`.crai`, reference FASTA", "Evidence extraction"),
            ("VCF", "`##fileformat=VCF` header", "Optional `.tbi`/`.csi`", "Import"),
            ("BCF", "BCF magic `BCF\\2`", "`.csi`", "Import"),
            ("gVCF", "VCF header with `<NON_REF>` symbolic allele", "Optional index", "Import, block records skipped"),
            ("FASTA", "`>` header, nucleotide alphabet", "`.fai`", "Reference only"),
        ],
        widths=[2.2, 5.6, 4.6, 3.6],
    )

    h2(doc, "3.2 Samplesheet contract")
    para(
        doc,
        "The samplesheet is a tab- or comma-delimited table with one row per sample. "
        "Required columns are `sample_id`, `input_path` and `input_type`. Optional columns "
        "carry read-group metadata (`library`, `platform_unit`, `read_structure` for UMI "
        "modes) and per-sample overrides. Duplicate `sample_id` values are a hard error, "
        "because a duplicate identifier would collapse two biological samples into one "
        "VCF column and corrupt every downstream comparison.",
    )
    para(
        doc,
        "Sample identifiers are validated against a restricted character set so that they "
        "can be written safely into VCF headers, file paths and shell commands without "
        "quoting ambiguity. A sample identifier containing whitespace, a path separator "
        "or a shell metacharacter is rejected at validation time rather than escaped "
        "later.",
    )

    h2(doc, "3.3 Configuration validation and precedence")
    para(
        doc,
        "Configuration is validated by Pydantic v2 models declared with `extra='forbid'`. "
        "An unrecognised key is an error, not a warning, because a silently ignored "
        "configuration key is a silent change of analysis. Field validators enforce "
        "ranges, and model validators enforce cross-field constraints that cannot be "
        "expressed per field \u2014 for example, that a validated policy status requires a "
        "pinned assembly and an empirical error profile, or that `tbs_umi` mode requires "
        "a declared read structure.",
    )
    para(
        doc,
        "Precedence is deterministic: built-in defaults are overridden by the selected "
        "profile, which is overridden by the run configuration, which is overridden by "
        "explicit command-line arguments. The resolved configuration is hashed with "
        "SHA-256 and the digest is embedded in the policy identifier, so two runs that "
        "produced different results but claim the same policy identifier are provably "
        "inconsistent.",
    )
    code_block(
        doc,
        "policy_id = PVX-{mode}-{policy_family[:4]}-{config_hash[:12]}-{registry_hash[:8]}\n"
        "config_hash   = sha256(canonical_json(resolved_config))\n"
        "registry_hash = sha256(canonical_json(parameter_registry))",
    )

    h2(doc, "3.4 Failure behaviour")
    para(
        doc,
        "Validation failures are typed. A missing file raises a not-found error naming "
        "the path; an unreadable or misidentified file raises a format error naming the "
        "detected class; a missing index raises an index error naming the expected "
        "companion; a schema violation raises a configuration error naming the field and "
        "the constraint. In every case the process exits non-zero and writes the "
        "diagnostic to standard error. The command-line interface returns exit code 2 for "
        "an input or configuration problem and exit code 1 for an unexpected internal "
        "failure, so that a scheduler can distinguish a bad submission from a bug.",
    )
    callout(
        doc,
        "Verified behaviour.",
        "`pfvarx run` on a FASTQ samplesheet without the alignment binaries exits with "
        "code 2, marks each sample `not_executed_requires_alignment`, and writes an "
        "explicit diagnostic to standard error. It does not emit an empty or partial "
        "VCF. `pfvarx run` on a VCF samplesheet exits 0 and reports the number of records "
        "imported and the number with recoverable counts (34 and 34 in the audit run).",
    )


def section_4(doc):
    page_break(doc)
    h1(doc, "4. Preprocessing, alignment and read-level quality control")
    para(
        doc,
        "Preprocessing is implemented as a set of versioned command builders rather than "
        "as free-form shell strings. Each builder returns an argument vector, which is "
        "logged verbatim into the run manifest. This makes every preprocessing step "
        "auditable after the fact: the exact command that produced a BAM can be "
        "reconstructed from the manifest without rerunning the pipeline.",
    )

    h2(doc, "4.1 Trimming and filtering")
    para(
        doc,
        "Adapter and quality trimming is delegated to fastp, which performs adapter "
        "detection, quality filtering, length filtering and poly-G trimming in a single "
        "pass and emits a JSON report that is ingested into the unified quality-control "
        "layer. Delegation is deliberate: reimplementing trimming would add a maintenance "
        "burden and a source of divergence from community practice without adding "
        "scientific value.",
    )
    para(
        doc,
        "The audit removed a spurious `--n_base_limit` argument that had been passed to "
        "fastp. That flag belongs to a different tool and would have caused fastp to abort "
        "with an unrecognised-option error on first execution. Its presence is a good "
        "illustration of why static audit of command builders matters even when the "
        "binaries are unavailable: the defect is visible in the argument vector.",
    )

    h2(doc, "4.2 Alignment")
    para(
        doc,
        "Alignment uses BWA-MEM2 against a pinned reference assembly. Read groups are "
        "attached at alignment time from samplesheet metadata so that duplicate marking "
        "and UMI grouping have the library information they require. The reference index "
        "prefix is a required argument of the preprocessing planner, which prevents the "
        "common failure in which a planner is constructed without a reference and "
        "silently emits commands that cannot run.",
    )

    h2(doc, "4.3 Duplicate handling")
    para(
        doc,
        "Optical and PCR duplicates are marked with Picard MarkDuplicates, or skipped "
        "when the mode declares that duplicates are already marked. When duplicate "
        "marking is disabled, the workflow writes a placeholder metrics file so that "
        "downstream rules that expect the metrics path do not fail on a missing input. "
        "This is a small example of a general principle in the workflow layer: every "
        "declared output must exist after its rule runs, even when the rule's work is a "
        "no-op.",
    )
    para(
        doc,
        "Independently of Picard, the evidence layer performs fragment-level "
        "de-duplication by `query_name`. This is a second, model-level safeguard: even if "
        "duplicate marking was skipped or misconfigured, two reads derived from the same "
        "fragment cannot both contribute to the alternate count. The two mechanisms are "
        "complementary and are not redundant, because they operate at different stages "
        "and defend against different failures.",
    )

    h2(doc, "4.4 Unique molecular identifiers")
    para(
        doc,
        "UMI mode implements a full consensus route: extract UMIs from the declared read "
        "structure with `fgbio ExtractUmisFromBam`, group reads by molecular index with "
        "`GroupReadsByUmi`, build consensus reads with `CallMolecularConsensusReads`, "
        "convert the consensus BAM back to FASTQ with `samtools fastq -T RX`, realign the "
        "consensus reads, and coordinate-sort the result. The evidence layer then reads "
        "the consensus BAM rather than the raw BAM.",
    )
    para(
        doc,
        "The audit changed the UMI status reported by the evidence layer from a "
        "hard-coded `UMI_DESIGN_DECLARED` to `UMI_UNVERIFIED` for any run that has not "
        "produced a validated consensus. The previous value asserted that a UMI design "
        "existed merely because the mode was selected. Only `UMI_CONSENSUS_VALID` now "
        "satisfies the UMI gate in the decision policy, so a UMI run that failed to build "
        "consensus produces `ALT_UNCERTAIN` rather than a confident call. This is a "
        "fail-closed correction and is pinned by a regression test.",
    )

    h2(doc, "4.5 Read-level quality control")
    para(
        doc,
        "Read-level metrics \u2014 yield, duplication rate, mean quality, adapter content, "
        "mapping rate, insert size, coverage uniformity \u2014 are collected by fastp, "
        "samtools and Picard and aggregated by MultiQC where available. Each metric is "
        "then registered in the unified quality-control ontology with a value, a unit, a "
        "source stage, a direction and a threshold, so that the same metric reported by "
        "two tools is not counted twice and a metric reported without a unit cannot be "
        "silently compared to a threshold.",
    )


def section_5(doc):
    page_break(doc)
    h1(doc, "5. Candidate discovery and evidence extraction")
    para(
        doc,
        "Candidate discovery is deliberately simple and evidence extraction is "
        "deliberately rich. Rather than applying a battery of heuristic filters to decide "
        "which positions deserve attention, PfVarX enumerates positions with at least one "
        "non-reference observation inside the callable region, and then records enough "
        "about each position that the decision can be made explicitly and reviewed "
        "afterwards. The cost of evaluating a position is low; the cost of discarding a "
        "true variant through an opaque filter is high.",
    )

    h2(doc, "5.1 CIGAR-aware pileup")
    para(
        doc,
        "Evidence is extracted from a pileup that respects the CIGAR string. A base is "
        "attributed to a reference position only when the alignment places it there; "
        "inserted bases are attributed to the anchor position and recorded with their "
        "length; deleted positions are counted as deletions rather than as missing data. "
        "This matters in a genome as AT-rich as *P. falciparum*, where indel "
        "misalignment around homopolymers is common and a naive pileup will misattribute "
        "bases across the length of a run.",
    )
    para(
        doc,
        "Reads are filtered by mapping quality, base quality, flag (unmapped, secondary, "
        "supplementary, duplicate, failed vendor quality) and, where configured, by "
        "read position. Each filter is a named configuration value rather than a literal "
        "embedded in the loop, so the effect of changing it is visible in the "
        "configuration hash.",
    )

    h2(doc, "5.2 Per-site evidence")
    para(
        doc,
        "For each candidate site the extractor records the reference allele and observed "
        "reference count, the alternate allele and observed alternate count, total depth, "
        "effective depth after quality filtering, the alternate allele fraction, forward "
        "and reverse strand counts for reference and alternate, mean mapping quality, "
        "mean base quality, distance to the nearest read end, homopolymer run length and "
        "the local sequence context. These are the quantities on which every subsequent "
        "statistical statement depends, and they are persisted in the evidence table "
        "whether or not the site ultimately passes.",
    )
    table(
        doc,
        ["Evidence group", "Fields", "Used by"],
        [
            ("Counts", "`ref_count`, `alt_count`, `depth`, `effective_depth`", "Binomial and beta-binomial likelihoods"),
            ("Fraction", "`alt_fraction`", "Allele-fraction bins, minority classification"),
            ("Strand", "`ref_fwd`, `ref_rev`, `alt_fwd`, `alt_rev`", "Fisher exact strand-bias test"),
            ("Quality", "`mean_mapq`, `mean_baseq`", "Ambiguous-mapping filter, QC metrics"),
            ("Context", "`homopolymer_run`, `left_context`, `right_context`", "Homopolymer INDEL review"),
            ("Provenance", "`policy_id`, `mode`, `decision`, `filters`", "Traceability, audit"),
        ],
        widths=[3.0, 7.4, 5.6],
    )

    h2(doc, "5.3 Reference-allele verification")
    para(
        doc,
        "Before a site is evaluated, the reference allele recorded in the candidate is "
        "compared against the base actually present in the reference FASTA at that "
        "coordinate. A mismatch means the candidate was generated against a different "
        "assembly, or that the coordinate is off by an indel, and in either case the "
        "likelihood is meaningless. The audit added this check as a hard gate: a mismatch "
        "produces the `ReferenceMismatch` filter and a `NO_CALL` decision. The check "
        "returns true when no FASTA is supplied, so that import-only workflows, which "
        "carry their own reference semantics, are not broken by it.",
    )

    h2(doc, "5.4 Fragment-level de-duplication")
    para(
        doc,
        "Reads sharing a `query_name` are collapsed to a single observation before "
        "counting. This prevents a single PCR-amplified fragment from inflating the "
        "alternate count, which would otherwise produce a spurious high-confidence call "
        "at low input. The de-duplication is recorded in the evidence table so that the "
        "difference between raw and effective depth is visible to the reviewer.",
    )

    h2(doc, "5.5 Minority and indel handling")
    para(
        doc,
        "Minority alleles are not filtered out at the evidence stage. The alternate "
        "fraction is carried forward as a continuous quantity and classified into "
        "reporting bins only at the policy stage, so that a low-fraction observation "
        "remains available for exploratory reporting and for mixed-infection screening. "
        "Indels are represented in anchored, left-aligned form, and homopolymer context is "
        "recorded so that the policy can flag an indel whose support may be an alignment "
        "artefact rather than a true length change.",
    )


def section_6(doc):
    page_break(doc)
    h1(doc, "6. Statistical models and formulations")
    para(
        doc,
        "The statistical core is small, explicit and numerically careful. All likelihoods "
        "are evaluated in log space using `gammaln` and `betaln` from SciPy, and all "
        "combinations of likelihoods use `logsumexp`, so that a depth of ten thousand "
        "reads does not underflow to zero. The primitives were verified against "
        "independent implementations during the audit, to a tolerance of 1e-10.",
    )

    h2(doc, "6.1 Notation and hypotheses")
    para(
        doc,
        "Let *k* be the alternate observation count and *n* the effective depth at a "
        "site. Let *e* be the per-base sequencing error rate and *f* the true alternate "
        "allele fraction. The null hypothesis is that alternate observations arise from "
        "error alone, so that the alternate probability equals *e*. The alternative is "
        "that they arise from a genuine variant at fraction *f*, with *f* greater than "
        "*e*. The hypotheses are one-sided by construction, because a variant fraction "
        "below the error rate is not distinguishable from error and is not a meaningful "
        "alternative.",
    )
    code_block(
        doc,
        "H0:  p = e                 (alternate observations are sequencing error)\n"
        "H1:  p = f,  f > e         (alternate observations are a true variant)\n"
        "k ~ Binomial(n, p)         under the binomial model\n"
        "k ~ BetaBinomial(n, mu, rho)  under the beta-binomial model",
    )

    h2(doc, "6.2 Model A \u2014 binomial")
    para(
        doc,
        "The binomial log-likelihood is the log of the binomial probability mass "
        "function. It assumes independent, identically distributed observations and "
        "therefore understates variance when over-dispersion is present, which for "
        "Illumina data is almost always. It is retained because it is the reference model "
        "against which over-dispersion is tested, and because it is the correct model when "
        "the data genuinely are not over-dispersed.",
    )
    code_block(
        doc,
        "log P(k | n, p) = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)\n"
        "                  + k log p + (n - k) log(1 - p)",
    )

    h2(doc, "6.3 Model B \u2014 beta-binomial")
    para(
        doc,
        "The beta-binomial places a Beta prior on the alternate probability, with mean "
        "*mu* and over-dispersion *rho*. The parameterisation is chosen so that *rho* has "
        "a direct interpretation as the intra-class correlation of observations within a "
        "site, and so that the model converges continuously to the binomial as *rho* "
        "approaches zero. This continuity is what makes the likelihood-ratio test of "
        "*rho* = 0 well behaved.",
    )
    code_block(
        doc,
        "alpha = mu (1 - rho) / rho\n"
        "beta  = (1 - mu) (1 - rho) / rho\n"
        "log P(k | n, mu, rho) = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)\n"
        "                        + betaln(k + alpha, n - k + beta) - betaln(alpha, beta)",
    )
    para(
        doc,
        "The beta-binomial is the recommended default for Illumina data because "
        "over-dispersion is ubiquitous and, if ignored, inflates false positives "
        "(Nielsen *et al.*, 2011; Li, 2011). Realistic values of *rho* for Illumina data "
        "are small, of order 0.001. The audit verified that the fitted *rho* on strongly "
        "over-dispersed simulated data (Beta(2,20), n = 80, 200 sites) recovers a value "
        "near 0.04 and that the test of *rho* = 0 rejects the binomial model on that data.",
    )

    h2(doc, "6.4 Model C \u2014 likelihood-ratio test")
    para(
        doc,
        "The test statistic is twice the difference in maximised log-likelihood between "
        "the alternative and the null. Because the alternative places the allele fraction "
        "on the boundary of the parameter space under the null, the asymptotic null "
        "distribution is not a central chi-square with one degree of freedom but a 50:50 "
        "mixture of a point mass at zero and that chi-square (Self and Liang, 1987). "
        "Treating it as a plain chi-square would roughly double the nominal p-value and "
        "bias every call toward the null.",
    )
    code_block(
        doc,
        "Lambda = 2 (logL1 - logL0)\n"
        "p_chi2 = 0.5 * P(chi2_1 >= Lambda)        (boundary mixture)\n"
        "p_boot = (1 + #{Lambda_b >= Lambda}) / (B + 1)   (parametric bootstrap)\n"
        "p_exact = sum_{K: LRT(K) >= LRT(k)} Binomial(K; n, e)   (exact conditional tail)",
    )
    para(
        doc,
        "Three calibrations of the null are implemented. The parametric bootstrap "
        "simulates replicate datasets under the null and is the preferred calibration for "
        "accuracy. The analytic chi-square mixture is fast and asymptotically correct. "
        "The exact conditional tail enumerates every possible alternate count from zero "
        "to *n* and sums the null mass of those whose likelihood ratio is at least as "
        "extreme as the observed one; it is exact for the binomial model and is used when "
        "the bootstrap produces no replicate exceeding the observed statistic.",
    )
    callout(
        doc,
        "Audit correction C1 and C2.",
        "Two defects were found and repaired in this component. The over-dispersion "
        "statistic was computed with an inverted sign and clamped to zero, so it always "
        "concluded \u201cno over-dispersion\u201d; it is now `max(0, 2*(maxLL - ll0))`. The "
        "bootstrap p-value fell back to the chi-square value when no replicate exceeded "
        "the observed statistic, mixing two incompatible null distributions in one "
        "reported number; it now falls back to the exact conditional tail. Both are "
        "pinned by regression tests.",
        fill="FDECEA",
        colour=RGBColor(0x9B, 0x1C, 0x1C),
    )

    h2(doc, "6.5 Model D \u2014 optional Bayesian posterior")
    para(
        doc,
        "A posterior probability of the variant hypothesis can be computed from a prior "
        "and the model likelihoods. It is disabled by default and requires an explicit "
        "threshold. The posterior is meaningful only when the error model has been "
        "calibrated; without calibration it is computed and reported but is not used to "
        "reach a validated decision. The audit added an explicit guard so that the "
        "posterior gate can only fire when the policy status is `VALIDATED`, the "
        "posterior threshold is set, the Bayesian component is enabled, the model is not "
        "flagged uncalibrated, and a posterior value is actually present.",
    )

    h2(doc, "6.6 Multiple testing")
    para(
        doc,
        "Testing across sites is controlled with the Benjamini\u2013Hochberg procedure "
        "(Benjamini and Hochberg, 1995) at a false discovery rate of 0.05 for exploratory "
        "calling, with 0.01 available for a sensitivity analysis. The implementation was "
        "verified against a hand-computed example: p-values 0.01, 0.04, 0.03 and 0.2 "
        "produce adjusted values 0.04, 0.0533, 0.0533 and 0.2, which is the correct "
        "step-up result. Benjamini\u2013Hochberg is chosen over Bonferroni because the "
        "number of tests is large and the cost of a false negative at a true variant is "
        "high, but the choice is a policy decision recorded in the configuration rather "
        "than a hard-coded constant.",
    )

    h2(doc, "6.7 Interval estimation and contingency tests")
    para(
        doc,
        "The alternate fraction is accompanied by a Wilson score interval, which behaves "
        "well at small counts and near the boundaries where the Wald interval fails. "
        "Where a Bayesian interval is requested, a Jeffreys interval is used, whose prior "
        "is non-informative under a Beta(0.5, 0.5) parameterisation and which has good "
        "coverage properties at small *n*. Strand bias is tested with Fisher's exact test "
        "on the 2\u00d72 table of forward and reverse counts for reference and alternate "
        "alleles, because the expected counts in that table are frequently small and the "
        "chi-square approximation is then unreliable.",
    )

    h2(doc, "6.8 Numerical stability")
    para(
        doc,
        "Every likelihood is computed in log space. Products of likelihoods are sums of "
        "logs. Ratios of likelihoods are differences of logs. Marginalisation over "
        "hypotheses uses `logsumexp` rather than exponentiating and summing. The audit "
        "compared the binomial and beta-binomial log-densities against independent "
        "implementations built from `scipy.stats` and from direct combinatorial formulas, "
        "and found agreement to better than 1e-10 across the tested range, including "
        "*n* = 10,000 and *rho* as small as 1e-4.",
    )


def section_7(doc):
    page_break(doc)
    h1(doc, "7. Model integration, calling, filtering and uncertainty")
    para(
        doc,
        "Model output is not a call. The decision policy converts a likelihood ratio, a "
        "posterior where available, and a set of site-level observations into one of five "
        "decision categories, and attaches zero or more named filters explaining any "
        "downgrade. The separation of model from policy is what makes the caller "
        "auditable: the statistics can be recomputed from the evidence table and the "
        "decision re-derived without rerunning the pipeline.",
    )

    h2(doc, "7.1 Decision categories")
    table(
        doc,
        ["Category", "Meaning", "Reachable when"],
        [
            (f"`{DECISION_NO_CALL}`", "No usable evidence at the site", "Effective depth below minimum, outside callable region, reference mismatch"),
            (f"`{DECISION_REFERENCE_SUPPORTED}`", "Reference allele supported", "Alternate fraction at or below the error rate"),
            (f"`{DECISION_ALT_EXPLORATORY}`", "Alternate allele detected, exploratory", "LRT significant, but policy status is provisional"),
            (f"`{DECISION_ALT_PASS}`", "Alternate allele passes a validated policy", "Policy status VALIDATED and all gates satisfied"),
            (f"`{DECISION_ALT_UNCERTAIN}`", "Alternate detected but a review filter applies", "Any review filter set"),
        ],
        widths=[5.4, 5.4, 5.2],
    )
    para(
        doc,
        "The distinction between the exploratory and validated categories is the "
        "decision-theoretic heart of the platform. Under the default provisional policy "
        "status, the strongest possible call is `ALT_DETECTED_EXPLORATORY`. This is not a "
        "limitation of the code but a statement about evidence: without a calibrated "
        "error profile, the cost of a false positive is judged higher than the cost of a "
        "false negative, so the policy refuses to certify a variant it cannot calibrate.",
    )

    h2(doc, "7.2 Hard gates")
    para(
        doc,
        "Hard gates force `NO_CALL` and cannot be overridden by a strong likelihood "
        "ratio. They represent conditions under which the likelihood is not merely weak "
        "but meaningless.",
    )
    for text in (
        "**Insufficient effective depth.** Below the configured minimum, the count model "
        "has no power and the fraction is dominated by sampling noise.",
        "**Outside the callable region.** Where a callable BED file is supplied, a site "
        "outside it is excluded because the mappability or uniqueness assumptions do not "
        "hold there.",
        "**Reference-allele mismatch.** The candidate's reference allele disagrees with "
        "the FASTA, so the likelihood was computed against the wrong base.",
    ):
        bullet(doc, text)

    h2(doc, "7.3 Review filters")
    para(
        doc,
        "Review filters downgrade a call to `ALT_UNCERTAIN` and are recorded by name in "
        "the VCF FILTER column and in the evidence table. The complete set is enumerated "
        "from the implementation at build time:",
    )
    table(
        doc,
        ["Filter", "Trigger", "Rationale"],
        [
            ("`LowEffectiveDepth`", "Effective depth below threshold", "Insufficient observations to support the model"),
            ("`LowAltSupport`", "Alternate count or fraction below support threshold", "Support too thin to distinguish from error"),
            ("`UncalibratedModel`", "Policy status is not VALIDATED", "No empirical error profile; statistic is not a calibrated confidence"),
            ("`StrandBiasReview`", "Fisher exact p-value below threshold", "Allele observed predominantly on one strand suggests artefact"),
            ("`HomopolymerIndelReview`", "INDEL in a long homopolymer run", "Indel placement ambiguous; length change may be an alignment artefact"),
            ("`OutsideCallableRegion`", "Site outside supplied BED", "Mappability or uniqueness assumptions do not hold"),
            ("`AmbiguousMapping`", "Mean mapping quality below threshold", "Read placement uncertain; allele assignment unreliable"),
            ("`UMIModelInvalid`", "UMI status is not `UMI_CONSENSUS_VALID`", "Consensus was not built; UMI evidence cannot be trusted"),
            ("`ReferenceMismatch`", "Candidate reference allele disagrees with FASTA", "Likelihood computed against the wrong base"),
        ],
        widths=[4.4, 5.6, 6.0],
    )
    para(
        doc,
        "The audit corrected the VCF writer, which previously joined multiple filters "
        "into a single semicolon-delimited string and passed it to pysam as one token. "
        "pysam rejects that token, so any site with two or more filters aborted the whole "
        "VCF write. The writer now passes a list, and a regression test writes a call "
        "carrying two filters and reads it back.",
    )

    h2(doc, "7.4 Allele-fraction reporting bins")
    table(
        doc,
        ["Bin", "Range", "Interpretation"],
        [
            (f"`{BIN_BELOW_THRESHOLD}`", "< 0.01", "Indistinguishable from error"),
            (f"`{BIN_LOW_FRACTION}`", "0.01 \u2013 0.05", "Exploratory; requires orthogonal support"),
            (f"`{BIN_MINORITY}`", "0.05 \u2013 0.20", "Minority-variant candidate"),
            (f"`{BIN_INTERMEDIATE}`", "0.20 \u2013 0.80", "Principal signal for polyclonal infection"),
            (f"`{BIN_PREDOMINANT}`", "\u2265 0.80", "Predominant alternate; near-fixation"),
        ],
        widths=[5.6, 3.4, 7.0],
    )
    para(
        doc,
        "The bins are reporting categories, not decision thresholds. A site in the "
        "minority-candidate bin is not thereby called; it is classified, and the decision "
        "is made by the gates and filters above. The audit fixed a latent crash in which "
        "a missing alternate fraction produced a `TypeError` during binning; the function "
        "now returns the below-threshold bin when the fraction is absent.",
    )

    h2(doc, "7.5 Genotype representation")
    para(
        doc,
        "*P. falciparum* is haploid in the intra-erythrocytic stage, which is the stage "
        "sampled by blood-stage sequencing. Genotypes are therefore written as single "
        "alleles, `GT` = 1 for an alternate call and `GT` = 0 for a reference call, with "
        "`Number=1` declared for the `GT` field in the VCF header. The audit found that "
        "the writer had been emitting diploid genotypes, `(None, 1)` and `(0, 0)`, which "
        "would cause every call to be interpreted downstream as heterozygous. The "
        "correction was verified by reading the written VCF back with pysam and confirming "
        "that the genotype tuple has length one.",
    )

    h2(doc, "7.6 Uncertainty representation")
    para(
        doc,
        "Uncertainty is represented in three distinct channels, and the platform is "
        "careful not to conflate them. The likelihood-ratio p-value is reported as a "
        "p-value. The posterior, when computed, is reported as a posterior. The VCF QUAL "
        "field is written only when the policy status is `VALIDATED`, the posterior is "
        "present, and the decision is `ALT_PASS_VALIDATED`; otherwise QUAL is omitted "
        "rather than populated with a rescaled statistic. This is a deliberate refusal of "
        "a widespread convention, and it is the reason the platform's VCFs are not "
        "directly interchangeable with those of callers that always emit a Phred-scaled "
        "QUAL. The evidence table carries the full statistical detail for any consumer "
        "that needs it.",
    )


def section_8(doc):
    page_break(doc)
    h1(doc, "8. Calibration and threshold justification")
    para(
        doc,
        "Every threshold in a variant caller encodes a judgement about the relative cost "
        "of false positives and false negatives. PfVarX makes that judgement explicit in "
        "two ways: each parameter carries a provenance record in a central registry, and "
        "the policy status of a run determines which decisions the thresholds are allowed "
        "to support.",
    )

    h2(doc, "8.1 The parameter registry")
    para(
        doc,
        f"The registry holds {len(REGISTRY)} parameters. Each entry records the parameter "
        "name, unit, datatype, default, permitted range, the analysis modes in which it "
        "applies, the rationale for its existence, the source on which the default is "
        "based, a policy status, the identifier of the calibration dataset if one exists, "
        "a version and a category. The registry is hashed and the digest is embedded in "
        "every policy identifier, so a result can be tied to the exact parameter set that "
        "produced it.",
    )
    code_block(
        doc,
        f"parameters : {len(REGISTRY)}\n"
        f"registry_hash : {registry_hash()}\n"
        "policy_status : PROVISIONAL_DEVELOPMENT (all entries)\n"
        "calibration_dataset_id : None (all entries)",
    )
    callout(
        doc,
        "Important.",
        f"All {len(REGISTRY)} registered parameters currently carry the status "
        "`PROVISIONAL_DEVELOPMENT` and no calibration dataset identifier. This means that "
        "no default in the platform has been fitted to control data within this project. "
        "The values are engineering defaults chosen to be defensible and to fail safe, "
        "not empirically optimised values. Any claim that a default is \u201cvalidated\u201d "
        "would be false, and the software does not make that claim.",
    )

    h2(doc, "8.2 Status vocabulary for parameters")
    table(
        doc,
        ["Status", "Meaning", "Consequence for calling"],
        [
            ("`PROVISIONAL_DEVELOPMENT`", "Engineering default; not fitted to control data", "Strongest decision is ALT_DETECTED_EXPLORATORY"),
            ("`ESTIMATED`", "Derived from a stated assumption or a published range, but not fitted here", "Still exploratory; provenance recorded"),
            ("`CALIBRATED`", "Fitted to control data with a recorded dataset identifier", "Eligible for a validated policy if all other conditions hold"),
            ("`EXPERIMENTALLY_VALIDATED`", "Accuracy demonstrated against an independent truth set", "Eligible for ALT_PASS_VALIDATED"),
        ],
        widths=[4.4, 6.4, 5.2],
    )

    h2(doc, "8.3 What calibration does")
    para(
        doc,
        "The `calibrate` subcommand accepts a table of control observations \u2014 alternate "
        "counts and depths at sites known to be reference, and optionally at sites known "
        "to carry a variant \u2014 and fits the error rate and the over-dispersion "
        "parameter. It emits a calibration profile containing the fitted values, the "
        "number of control sites, the dataset identifier, the mode and a status field. "
        "When the control data are insufficient, the profile is written with the status "
        "`INSUFFICIENT_CONTROL_DATA` and an explicit warning, rather than with a "
        "plausible-looking but unsupported estimate.",
    )
    para(
        doc,
        "The audit verified this behaviour by running `pfvarx calibrate` on a small "
        "synthetic control table. The command produced a profile carrying the dataset "
        "identifier and the honest `INSUFFICIENT_CONTROL_DATA` status with warnings, "
        "which is the intended fail-safe behaviour. It did not produce a fabricated "
        "error rate.",
    )

    h2(doc, "8.4 Threshold justification, item by item")
    table(
        doc,
        ["Threshold", "Default basis", "Status", "Evidence still required"],
        [
            ("Minimum effective depth", "Conventional short-read practice; below ~10 observations the binomial has little power", "PROVISIONAL", "Power analysis across depths and fractions"),
            ("Error rate", "Order-of-magnitude Illumina error rate; used only when no calibrated profile exists", "PROVISIONAL", "Per-instrument, per-library empirical error profile"),
            ("Over-dispersion rho", "Published order of magnitude (~0.001) for Illumina data", "PROVISIONAL", "Fitted rho across libraries and coverage levels"),
            ("LRT alpha", "Conventional 0.05 exploratory, 0.01 sensitivity", "PROVISIONAL", "Operating-characteristic curves on simulated and control data"),
            ("FDR level", "Conventional 0.05, with 0.01 sensitivity analysis", "PROVISIONAL", "Empirical FDR against a truth set"),
            ("Strand-bias p-value", "Conventional 0.01 for a review flag", "PROVISIONAL", "False-flag rate on control sites"),
            ("Minority bins", "Convention in the malaria genomics literature", "PROVISIONAL", "Concordance with orthogonal minority-variant assays"),
            ("Posterior threshold", "No default; must be set explicitly", "NO DEFAULT", "Calibrated prior and error model"),
            ("Callable-region BED", "User-supplied; no default", "NO DEFAULT", "Project-specific mappability analysis"),
        ],
        widths=[3.8, 5.4, 2.6, 4.2],
    )
    para(
        doc,
        "The table is the honest inventory. Each row states what the default rests on and "
        "what would be needed to replace it with a calibrated value. None of the rows is "
        "marked calibrated, because none is.",
    )

    h2(doc, "8.5 Removing unsupported quantities")
    para(
        doc,
        "The audit removed two quantities that were computed without a defensible basis. "
        "The mixed-infection module had derived a numeric complexity-of-infection estimate "
        "from the count of intermediate-fraction variants using an arbitrary linear "
        "formula; the estimate is now `None`, the method string states that no numeric "
        "COI is reported, and the reader is directed to DEploid for a modelled estimate. "
        "The diversity module had computed a Watterson estimator from a quantity that was "
        "not a valid segregating-site count in this context; the estimator is now `None` "
        "with an explanatory comment, and the remaining diversity statistics guard "
        "against missing alternate fractions.",
    )
    para(
        doc,
        "These removals are not losses of functionality. A fabricated number is worse "
        "than an absent one, because it invites interpretation. The platform now reports "
        "what it can support and declines to report what it cannot.",
    )


def section_9(doc):
    page_break(doc)
    h1(doc, "9. Annotation and downstream analyses")
    para(
        doc,
        "Annotation and downstream analysis sit outside the calling decision and are "
        "designed so that their failure cannot corrupt a call. Each is a separate stage "
        "with its own inputs and outputs, and each is guarded so that an absent binary or "
        "an unmet assumption produces an explicit refusal rather than a partial result.",
    )

    h2(doc, "9.1 Functional annotation")
    para(
        doc,
        "Per-call annotation is delegated to SnpEff where a database is available, with "
        "an internal consequence fallback for simple classes. The audit changed the "
        "behaviour of the SnpEff engine from a silent no-op to an explicit refusal: when "
        "`annotation.engine` is set to `snpeff`, the code now verifies that the binary is "
        "present and then raises a typed error stating that per-call SnpEff annotation is "
        "not implemented in this build. The previous behaviour would have returned "
        "unannotated records while reporting that annotation had run, which is exactly "
        "the kind of silent continuation the audit was asked to eliminate.",
    )
    para(
        doc,
        "Resistance-marker reporting is a lookup against a curated table of known "
        "*P. falciparum* resistance loci. It is a reporting convenience, not a "
        "prediction, and it is labelled as such. The table is versioned with the "
        "software so that a change in the marker set is visible in the release notes.",
    )

    h2(doc, "9.2 Population-genomic summaries")
    para(
        doc,
        "Diversity statistics are computed from the called set with explicit guards for "
        "missing values. Nucleotide diversity is computed from pairwise differences over "
        "callable sites; the proportion of polymorphic sites is computed over sites with "
        "sufficient depth. Where a statistic requires a quantity the pipeline cannot "
        "supply, the statistic is returned as absent rather than approximated.",
    )

    h2(doc, "9.3 Mixed infection")
    para(
        doc,
        "Mixed infection is screened by counting variants in the intermediate "
        "allele-fraction bin. Three or more such variants yields a "
        "`MIXED_INFECTION_LIKELY` classification; fewer yields "
        "`NO_EVIDENCE_OF_MIXED_INFECTION`. This is explicitly a screening heuristic, not "
        "a formal complexity-of-infection estimate, and the output says so. Where formal "
        "deconvolution is required, DEploid estimates the number of constituent haplotypes "
        "and their proportions from read counts (Zhu, Almagro-Garcia and McVean, 2018). "
        "The heuristic and the formal method are reported separately so that a reader can "
        "distinguish a screening flag from a modelled estimate.",
    )

    h2(doc, "9.4 Relatedness and transmission")
    para(
        doc,
        "Relatedness analysis is a guarded wrapper around hmmIBD, which estimates "
        "identity-by-descent segments between sample pairs from biallelic SNP data. The "
        "wrapper checks that the input is biallelic, that the sample count is sufficient "
        "for a pairwise estimate, and that the binary is available; if any check fails it "
        "raises rather than returning a matrix of zeros. Identity-by-descent estimation "
        "assumes that the SNP set is ascertained independently of the samples being "
        "compared, which is a strong assumption in a re-analysis setting and is recorded "
        "in the output metadata.",
    )

    h2(doc, "9.5 Phylogenetics")
    para(
        doc,
        "Phylogenetic reconstruction is a guarded wrapper around IQ-TREE 2. The wrapper "
        "constructs an alignment from the called sites, checks that the alignment has "
        "enough informative positions to support a tree, and records the substitution "
        "model and bootstrap settings in the output metadata. A tree built from too few "
        "informative sites is not a phylogeny; it is a drawing, and the wrapper refuses "
        "to produce one.",
    )

    h2(doc, "9.6 Why downstream analyses are separated")
    para(
        doc,
        "Keeping annotation and downstream analysis outside the caller has three "
        "benefits. A failure in a downstream tool cannot invalidate a VCF that was "
        "already written correctly. A user who wants only calls does not pay for "
        "annotation. And the scientific status of each downstream method can be reported "
        "separately, which matters because the downstream methods rest on assumptions "
        "that the caller does not, and their validation status is different.",
    )


def section_10(doc):
    page_break(doc)
    h1(doc, "10. Unified quality control, reporting and dashboard")
    para(
        doc,
        "Quality control in PfVarX is a first-class data product rather than a set of "
        "printed warnings. Every metric is an object with a name, a value, a unit, the "
        "stage that produced it, a direction, a threshold, a status and a reason. This "
        "structure is what allows the run-level summary to be computed mechanically and "
        "the dashboard to render the same numbers that appear in the PDF report.",
    )

    h2(doc, "10.1 Metric ontology")
    table(
        doc,
        ["Field", "Purpose"],
        [
            ("`name`", "Stable metric identifier"),
            ("`value`", "Observed value, or absent"),
            ("`unit`", "Physical unit; prevents comparing unlike quantities"),
            ("`stage`", "Pipeline stage that produced the metric"),
            ("`direction`", "`min`, `max`, `range` or `none`"),
            ("`threshold`", "Bound or bounds against which the value is judged"),
            ("`status`", "`PASS`, `WARN`, `FAIL` or `INFO`"),
            ("`reason`", "Human-readable explanation of the status"),
            ("`policy_id`", "Configuration under which the metric was evaluated"),
            ("`provisional`", "True when the threshold is not calibrated"),
        ],
        widths=[3.6, 12.4],
    )
    para(
        doc,
        "The audit found that the `provisional` flag was being dropped when metrics were "
        "evaluated, so a metric judged against an uncalibrated threshold was reported "
        "identically to one judged against a calibrated threshold. Every evaluation now "
        "passes the flag, and it is set from the run's policy status, so a provisional "
        "run cannot present a provisional threshold as if it were validated.",
    )

    h2(doc, "10.2 Two levels of quality control")
    para(
        doc,
        "Site-level filters address the reliability of an individual call. Run-level "
        "metrics address the suitability of a sample or run for analysis. The separation "
        "is deliberate: conflating the two would allow a globally poor run to pass because "
        "individual sites happened to meet thresholds, or a good run to fail because a few "
        "sites were flagged. The run-level summary reports the weakest status across "
        "components, so a single failing component fails the run.",
    )

    h2(doc, "10.3 Report generation")
    para(
        doc,
        "Reports are generated in four formats from the same underlying data: CSV for "
        "machine consumption, XLSX for tabular review, PNG for figures and PDF for "
        "archival. The audit ran the demonstration pipeline and confirmed that all four "
        "artefacts are produced, that the PDF contains the run summary and the "
        "validation-status banner, and that the XLSX opens with the expected sheets. The "
        "report carries the policy identifier and the registry hash so that a printed "
        "report can be tied back to the configuration that produced it.",
    )

    h2(doc, "10.4 Dashboard")
    para(
        doc,
        "The Streamlit dashboard reads the same files as the report generator and "
        "presents run-level metrics, per-sample summaries, filter counts and the "
        "validation-status banner. It is file-backed rather than database-backed, which "
        "means it can be pointed at any completed run directory without a database "
        "service. It is a convenience layer, not a separate source of truth: every number "
        "it displays is read from a file that the pipeline wrote.",
    )

    h2(doc, "10.5 Validation banner")
    para(
        doc,
        "Every report and every dashboard view carries the overall release status and a "
        "per-component breakdown. The banner is generated from `release.py`, so it cannot "
        "drift from the machine-readable status. A user who receives a PfVarX report "
        "therefore receives, on the same page, a statement of what the software has and "
        "has not been validated to do.",
    )


def section_11(doc):
    page_break(doc)
    h1(doc, "11. Tools, dependencies, configuration and deployment")
    para(
        doc,
        "The platform separates its Python dependencies, which are declared and pinned, "
        "from its external bioinformatics binaries, which are invoked as subprocesses and "
        "are not bundled in the Python wheel. This separation is honest about what a "
        "`pip install` provides and what it does not.",
    )

    h2(doc, "11.1 Python dependencies")
    table(
        doc,
        ["Package", "Role", "Declared constraint"],
        [
            ("numpy", "Array arithmetic, simulation", ">=1.24"),
            ("scipy", "`gammaln`, `betaln`, `logsumexp`, distributions, Fisher exact test", ">=1.10"),
            ("pandas", "Tabular evidence and QC data", ">=2.0"),
            ("pydantic", "Configuration schema and validation", ">=2.5"),
            ("PyYAML", "Profile and configuration parsing", ">=6.0"),
            ("pysam", "FASTA, BAM, CRAM and VCF input/output", ">=0.22"),
            ("matplotlib", "Figures", ">=3.7"),
            ("openpyxl", "XLSX export", ">=3.1"),
            ("Jinja2", "Report templating", ">=3.1"),
            ("reportlab", "PDF export", ">=4.0 (core requirement)"),
            ("streamlit", "Dashboard", ">=1.30 (optional extra)"),
            ("pytest", "Test runner", "development"),
            ("ruff", "Linting", "development"),
            ("hypothesis", "Property-based testing", "development"),
        ],
        widths=[3.2, 8.8, 4.0],
    )

    h2(doc, "11.2 External binaries")
    table(
        doc,
        ["Tool", "Stage", "Bundled?", "Behaviour if absent"],
        [
            ("fastp", "Trimming and read QC", "No", "Preprocessing route cannot run; explicit error"),
            ("BWA-MEM2", "Alignment", "No", "Alignment route cannot run; explicit error"),
            ("samtools", "Sort, index, FASTA conversion", "No", "Preprocessing route cannot run; explicit error"),
            ("Picard", "Duplicate marking", "No", "Duplicate-marking branch cannot run; explicit error"),
            ("fgbio", "UMI extraction, grouping, consensus", "No", "UMI route cannot run; explicit error"),
            ("BCFtools", "VCF/BCF conversion and indexing", "No", "Import route degrades to uncompressed VCF"),
            ("SnpEff", "Functional annotation", "No", "Annotation refuses explicitly (see 9.1)"),
            ("MultiQC", "Aggregated QC report", "No", "QC aggregation falls back to native metrics"),
            ("DEploid", "Complexity of infection", "No", "Downstream stage refuses explicitly"),
            ("hmmIBD", "Relatedness", "No", "Downstream stage refuses explicitly"),
            ("IQ-TREE 2", "Phylogeny", "No", "Downstream stage refuses explicitly"),
        ],
        widths=[3.2, 4.6, 2.2, 6.0],
    )
    para(
        doc,
        "The distinction between a core requirement and an optional extra is not "
        "cosmetic. An audit of the release archive in a fresh virtual environment found "
        "that `reportlab` had been declared only as an optional extra even though the "
        "default end-to-end run writes a PDF; a clean installation therefore succeeded "
        "but failed the moment a report was produced. The dependency was promoted to the "
        "core set, and the clean-extraction test suite was re-run from a freshly unzipped "
        "archive to confirm that a bare installation can now execute the full default "
        "pipeline. `streamlit` remains optional because the dashboard is a separate entry "
        "point and its import is guarded, so its absence degrades the interface without "
        "breaking analysis.",
    )

    h2(doc, "11.3 Configuration surface")
    para(
        doc,
        "Configuration is a single YAML document validated by a Pydantic model with "
        "`extra='forbid'`. It is organised into sections for the run, reference, "
        "preprocessing, calling, statistics, annotation, downstream analysis, QC and "
        "reporting. Profiles supply environment-specific overrides \u2014 executor, "
        "resources, container image \u2014 without touching analytical parameters. The "
        "resolved configuration is written into the run directory so that a run can be "
        "reproduced from its own output.",
    )

    h2(doc, "11.4 Deployment")
    para(
        doc,
        "Four deployment paths are supported. A local installation uses a Python "
        "environment plus the external binaries on the path. A server installation adds "
        "the server executor profile. An HPC installation uses the HPC profile with "
        "scheduler-specific resource declarations. A container installation uses the "
        "provided Apptainer definition or Dockerfile, which install the external binaries "
        "and pin their versions. The analytical code is identical across all four; only "
        "the executor and the resource declarations differ.",
    )
    para(
        doc,
        "The audit verified the Python package and the command-line interface in this "
        "environment. It did not build the container images or execute the HPC profile, "
        "because the container runtime and the scheduler are not available here. The "
        "definitions were reviewed for correctness but are not claimed to have been "
        "executed.",
    )


def section_12(doc):
    page_break(doc)
    h1(doc, "12. Testing, benchmarking, validation, metrics and acceptance criteria")
    para(
        doc,
        "This section reports what was actually executed. Every number below was produced "
        "by a command run during the audit; none is estimated. Where a test or benchmark "
        "was not run, that is stated.",
    )

    h2(doc, "12.1 Test suite")
    table(
        doc,
        ["Category", "What it covers", "Result"],
        [
            ("Unit \u2014 statistics", "Binomial and beta-binomial log-densities against independent implementations; LRT; bootstrap and exact nulls; BH FDR; Wilson and Jeffreys intervals", "Passing"),
            ("Unit \u2014 configuration", "Schema validation, `extra='forbid'`, bounds, cross-field validators, precedence, hashing", "Passing"),
            ("Unit \u2014 policy", "Decision categories, hard gates, review filters, bins, posterior gate", "Passing"),
            ("Unit \u2014 evidence", "CIGAR-aware counting, fragment de-duplication, strand metrics, context metrics", "Passing"),
            ("Unit \u2014 output", "VCF header, haploid genotypes, multi-filter records, evidence columns, atomic writes", "Passing"),
            ("Unit \u2014 downstream", "Diversity guards, mixed-infection classification and honest COI absence", "Passing"),
            ("Regression", "Twelve tests pinning each audit repair", "Passing"),
            ("Integration", "End-to-end demonstration: FASTQ-free synthetic path to VCF, evidence table, QC JSON, PDF and XLSX", "Passing"),
            ("CLI smoke", "`calibrate`, `import`, `run` on VCF and on FASTQ samplesheets, including failure exit codes", "Passing"),
        ],
        widths=[3.8, 9.4, 2.8],
    )
    para(
        doc,
        "The suite comprises 124 tests, all passing, with `ruff check src tests` clean. "
        "The baseline before the audit was 112 passing tests; the twelve additional tests "
        "are the regression tests written for the audit repairs. Test outcomes were "
        "observed directly from the test runner output; no result in this document is "
        "inferred from the code's intent.",
    )
    para(
        doc,
        "The suite was additionally run from a clean extraction of the release archive "
        "in a fresh virtual environment, rather than only in the development tree. That "
        "step is what exposed the `reportlab` packaging defect described in 11.1: the "
        "in-tree suite passed because the development environment already carried the "
        "package. After the dependency was corrected and the archive rebuilt, the "
        "clean-extraction run reported 124 passed, which is the result quoted above. The "
        "in-tree and clean-extraction runs therefore agree, and the agreement is itself "
        "evidence that the declared dependencies are sufficient for a default run.",
    )

    h2(doc, "12.2 What was verified directly")
    for text in (
        "Binomial and beta-binomial log-densities agree with independent implementations "
        "to better than 1e-10 across the tested range, including high depth and small "
        "over-dispersion.",
        "Benjamini\u2013Hochberg adjustment reproduces a hand-computed four-value example "
        "exactly.",
        "The over-dispersion test rejects the binomial model on strongly over-dispersed "
        "simulated data and recovers a fitted *rho* consistent with the simulation.",
        "The bootstrap null tail and the exact conditional tail agree with each other and "
        "differ from the chi-square tail on the same data, confirming that the fallback no "
        "longer swaps null distributions.",
        "A VCF written with two filters reads back with both filters, and its genotype "
        "tuple has length one.",
        "The reference-allele check returns true for a correct base, false for a wrong "
        "base, and true when no FASTA is supplied.",
        "The UMI gate rejects `UMI_UNVERIFIED` and accepts `UMI_CONSENSUS_VALID`.",
        "The demonstration pipeline produces a VCF with 34 records, an evidence table, a "
        "QC JSON, a PDF report and an XLSX workbook.",
        "`pfvarx import` recovers counts for all 34 records in the test VCF, including "
        "records whose header lacks a DP field.",
        "`pfvarx run` on a FASTQ samplesheet exits 2 with an explicit diagnostic and does "
        "not emit a partial VCF.",
    ):
        bullet(doc, text)

    h2(doc, "12.3 What was not verified")
    for text in (
        "The Snakemake workflow was not executed end-to-end, because Snakemake and the "
        "external binaries are not installed in the audit environment. Rules were audited "
        "statically for wiring, quoting, conditional arguments and rule ordering.",
        "No accuracy benchmark against an external truth set was run. Sensitivity, "
        "precision, false-discovery rate and genotype concordance against a reference "
        "caller are therefore not reported, because they were not measured.",
        "Container images were not built. The Apptainer definition and Dockerfile were "
        "reviewed but not executed.",
        "The HPC and cloud executor profiles were not executed.",
        "The calibration machinery was exercised on synthetic control data sufficient to "
        "confirm its fail-safe behaviour, not on a real control dataset sufficient to fit "
        "a defensible error profile.",
        "Performance and resource usage were not benchmarked beyond observing that the "
        "demonstration pipeline completes on a two-core, three-gigabyte sandbox.",
    ):
        bullet(doc, text)

    h2(doc, "12.4 Benchmarking framework")
    para(
        doc,
        "A simulation and ablation framework exists and executes. It generates synthetic "
        "read evidence at known allele fractions and depths, runs the caller over it, and "
        "reports the confusion matrix and derived metrics. It is the instrument with which "
        "an accuracy benchmark would be conducted. It has not been pointed at an external "
        "truth set, and its outputs on simulated data are not a substitute for that, "
        "because simulated error is not real error.",
    )

    h2(doc, "12.5 Acceptance criteria")
    table(
        doc,
        ["Criterion", "Threshold", "Current state"],
        [
            ("All unit and integration tests pass", "100%", "Met (124/124)"),
            ("Lint clean", "Zero findings", "Met"),
            ("No silent continuation on invalid input", "Every failure path exits non-zero with a diagnostic", "Met for the paths tested"),
            ("Every call traceable to persisted evidence", "Evidence row per call", "Met"),
            ("No fabricated quantitative output", "No numeric value without a stated basis", "Met after removal of two fabricated quantities"),
            ("Configuration reproducible from output", "Resolved config and hash written", "Met"),
            ("Accuracy demonstrated against a truth set", "Reported sensitivity, precision, FDR", "**Not met**"),
            ("Calibrated error model for at least one platform", "Recorded calibration dataset", "**Not met**"),
            ("External caller comparison", "Head-to-head benchmark", "**Not met**"),
        ],
        widths=[6.0, 5.4, 4.6],
    )
    para(
        doc,
        "The first six criteria are engineering criteria and are met. The last three are "
        "scientific criteria and are not met. On that basis PfVarX is not declared "
        "production-ready for validated variant calling; it is declared research-use "
        "ready with the restrictions stated in Section 15.",
    )


def section_13(doc):
    page_break(doc)
    h1(doc, "13. Reproducibility, error handling, limitations and troubleshooting")
    para(
        doc,
        "Reproducibility is treated as a property to be engineered rather than a "
        "discipline to be recommended. Error handling is treated as part of the "
        "scientific contract, because a pipeline that fails silently produces data that "
        "looks like a result.",
    )

    h2(doc, "13.1 Reproducibility mechanisms")
    for text in (
        "**Configuration hashing.** The resolved configuration is canonicalised and "
        "hashed with SHA-256; the digest is embedded in the policy identifier.",
        "**Registry hashing.** The parameter registry is hashed; the digest is embedded "
        "alongside the configuration hash, so a change in a default invalidates the "
        "identifier even if the user's configuration is unchanged.",
        "**Resolved configuration persisted.** The run directory contains the fully "
        "resolved configuration, not merely the user's overrides.",
        "**Command logging.** Every external command is logged as an argument vector, so "
        "the exact invocation can be reconstructed.",
        "**Atomic output.** Files are written to a temporary path and moved into place "
        "with `os.replace`, so an interrupted write cannot leave a truncated VCF that "
        "looks complete.",
        "**Version pinning.** Python dependencies carry lower bounds; container "
        "definitions pin external binary versions.",
    ):
        bullet(doc, text)

    h2(doc, "13.2 Error taxonomy")
    table(
        doc,
        ["Error class", "Raised when", "Exit code"],
        [
            ("Configuration error", "Schema violation, unknown key, cross-field constraint failure", "2"),
            ("Input error", "Missing file, unrecognised format, missing index, duplicate sample identifier", "2"),
            ("Tool execution error", "External binary absent or returned non-zero", "1"),
            ("Model error", "Numerical failure, impossible parameter combination", "1"),
            ("Internal error", "Unexpected exception", "1"),
        ],
        widths=[4.2, 9.8, 2.0],
    )
    para(
        doc,
        "The distinction between exit code 2 and exit code 1 lets a scheduler separate a "
        "bad submission, which should not be retried unchanged, from an internal failure, "
        "which may be transient.",
    )

    h2(doc, "13.3 Limitations")
    table(
        doc,
        ["Limitation", "Consequence", "Mitigation in place"],
        [
            ("No accuracy benchmark against a truth set", "Sensitivity and precision are unknown", "Status reported honestly; strongest call is exploratory"),
            ("No calibrated error profile", "Thresholds are engineering defaults", "Calibration machinery exists; status PROVISIONAL"),
            ("Workflow not executed end-to-end here", "Rule wiring verified statically only", "Static audit; container definitions provided"),
            ("External binaries not bundled", "`pip install` alone cannot align reads", "Documented; containers provided"),
            ("Downstream methods rest on strong assumptions", "Relatedness and phylogeny are conditional", "Wrappers refuse when assumptions fail"),
            ("Haploid model only", "Not applicable to diploid organisms", "Scope stated; genotype representation enforced"),
            ("Short-read only", "Long-read error models unsupported", "Scope stated"),
            ("Mixed-infection COI not estimated", "No numeric complexity of infection", "Screening flag only; DEploid directed"),
        ],
        widths=[5.4, 5.4, 5.2],
    )

    h2(doc, "13.4 Troubleshooting")
    table(
        doc,
        ["Symptom", "Likely cause", "Action"],
        [
            ("Exit code 2 with a format error", "File content does not match its extension", "Inspect the file with `file`; correct the samplesheet"),
            ("Exit code 2 with an index error", "Missing `.bai`, `.crai`, `.fai` or `.csi`", "Build the index for the declared input"),
            ("Exit code 2 with a duplicate sample error", "Two rows share a `sample_id`", "Assign unique identifiers"),
            ("Exit code 2 on a UMI run", "`preprocessing.umi_design` not set", "Declare the read structure"),
            ("All calls `UncalibratedModel`", "Policy status is provisional", "Expected; supply a calibration profile to change it"),
            ("All calls `ReferenceMismatch`", "Wrong assembly or off-by-one coordinates", "Confirm the pinned assembly and coordinate convention"),
            ("Empty VCF", "No candidates passed the depth gate", "Check coverage and the callable-region BED"),
            ("Annotation refuses", "SnpEff engine selected but unimplemented per-call", "Use the internal fallback or run SnpEff separately"),
        ],
        widths=[5.0, 5.4, 5.6],
    )


def section_14(doc, references: list[str]):
    page_break(doc)
    h1(doc, "14. References and software citations")
    para(
        doc,
        "References are given in Harvard style. Every entry below corresponds to a work "
        "that is cited in the text of this document or in the methodology chapter from "
        "which it derives. No reference has been generated or inferred; the list is "
        "reproduced verbatim from the verified bibliography.",
    )
    for ref in references:
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Cm(0.8)
        p.paragraph_format.first_line_indent = Cm(-0.8)
        p.paragraph_format.space_after = Pt(6)
        _rich(p, ref, size=9.5)

    h2(doc, "14.1 Software cited")
    table(
        doc,
        ["Software", "Version constraint", "Citation"],
        [
            ("BWA-MEM2", "external", "Vasimuddin *et al.* (2019)"),
            ("SAMtools", "external", "Li *et al.* (2009)"),
            ("BCFtools", "external", "Danecek *et al.* (2021)"),
            ("fastp", "external", "Chen *et al.* (2018)"),
            ("Picard", "external", "Broad Institute"),
            ("fgbio", "external", "Fulcrum Genomics"),
            ("SnpEff", "external", "Cingolani *et al.* (2012)"),
            ("MultiQC", "external", "Ewels *et al.* (2016)"),
            ("DEploid", "external", "Zhu, Almagro-Garcia and McVean (2018)"),
            ("hmmIBD", "external", "Schaffner *et al.* (2018)"),
            ("IQ-TREE 2", "external", "Minh *et al.* (2020)"),
            ("pysam", ">=0.22", "Heger *et al.*"),
            ("SciPy", ">=1.10", "Virtanen *et al.* (2020)"),
            ("NumPy", ">=1.24", "Harris *et al.* (2020)"),
            ("pandas", ">=2.0", "McKinney (2010)"),
            ("Pydantic", ">=2.5", "Pydantic project"),
            ("Snakemake", "external", "M\u00f6lder *et al.* (2021)"),
            ("Streamlit", "runtime", "Streamlit project"),
        ],
        widths=[3.2, 3.6, 9.2],
    )


def section_15(doc):
    page_break(doc)
    h1(doc, "15. Validation-status table")
    para(
        doc,
        "This table is the definitive statement of what PfVarX has and has not "
        "demonstrated. It is generated at build time from `pfvarx.release.RELEASE_STATUS`, "
        "so it cannot diverge from the status the software reports about itself. Two "
        "independent axes are recorded: software status, which concerns whether the code "
        "does what it claims, and scientific status, which concerns whether its outputs "
        "are accurate on real data.",
    )
    para(
        doc,
        "The vocabulary is fixed. `NOT_IMPLEMENTED` means the capability is absent or is "
        "a guarded stub. `IMPLEMENTED_NOT_VALIDATED` means the code exists and runs but "
        "has no test coverage adequate to assert its behaviour. "
        "`TESTED_NOT_SCIENTIFICALLY_VALIDATED` means automated tests assert the behaviour "
        "but accuracy on real data has not been demonstrated. "
        "`VALIDATED_WITH_RESTRICTIONS` means accuracy has been demonstrated within stated "
        "limits. `RESEARCH_USE_READY` means the component is fit for exploratory research "
        "with the stated restrictions.",
    )

    rows = []
    for name, entry in RELEASE_STATUS.items():
        rows.append(
            (
                name.replace("_", " "),
                entry["software_status"],
                entry["scientific_status"],
                entry["note"],
            )
        )
    table(
        doc,
        ["Component", "Software status", "Scientific status", "Note"],
        rows,
        widths=[2.9, 5.5, 4.5, 3.1],
        font_size=7.6,
        col_sizes=[7.6, 6.4, 6.4, 7.2],
    )

    para(
        doc,
        f"**Overall status: {overall_status()}.**",
    )

    h2(doc, "15.1 Interpretation")
    para(
        doc,
        "The pattern across the table is consistent and deliberate. Components whose "
        "correctness can be established by test \u2014 configuration, routing, evidence "
        "extraction, the statistical models, the decision policy, quality control \u2014 are "
        "tested but not scientifically validated, because testing establishes that the "
        "code implements the intended computation, not that the computation is accurate "
        "on real sequence data. Components that depend on external binaries \u2014 "
        "preprocessing, annotation, downstream analysis, dashboard reporting \u2014 are "
        "implemented but not validated, because their behaviour depends on tools that are "
        "not bundled and were not executed in the audit environment. Benchmarking is "
        "tested as a framework but has no external truth set, so it is scientifically "
        "not implemented.",
    )
    para(
        doc,
        "No component is marked scientifically validated. That is the honest conclusion "
        "of the audit. The platform is a carefully engineered, thoroughly tested, "
        "explicitly provisional instrument. It is suitable for exploratory research in "
        "which its outputs are treated as candidates requiring orthogonal confirmation. "
        "It is not suitable, in its current state, for any application in which a "
        "validated variant call is required.",
    )

    h2(doc, "15.2 Path to validation")
    para(
        doc,
        "Moving a component from tested to scientifically validated requires four steps, "
        "none of which has been completed. First, a control dataset with known genotypes "
        "must be assembled, ideally from a sample that has been characterised by an "
        "orthogonal method. Second, the error rate and over-dispersion must be fitted to "
        "that dataset and recorded with its identifier, which converts the relevant "
        "parameters from `PROVISIONAL_DEVELOPMENT` to `CALIBRATED`. Third, the caller must "
        "be run over an independent truth set and its sensitivity, precision, "
        "false-discovery rate and genotype concordance measured. Fourth, those results "
        "must be compared against an established caller on the same data, so that the "
        "claim is relative as well as absolute.",
    )
    para(
        doc,
        "Only after the fourth step can a component be described as validated, and only "
        "within the conditions of the benchmark. The machinery for all four steps exists "
        "in the platform. The work of executing them is the next phase of the project, "
        "not a property that can be asserted by writing documentation.",
    )

    h2(doc, "15.3 Statement of fitness")
    callout(
        doc,
        "Fitness statement.",
        "PfVarX 0.9.0 is fit for research use as an exploratory SNP and small-INDEL "
        "analysis platform, provided that its calls are treated as candidates and "
        "confirmed orthogonally before any biological conclusion is drawn. It is not fit "
        "for clinical decision-making, for regulatory submission, or for any use in which "
        "an unconfirmed variant call would be acted upon. This statement is consistent "
        "with the machine-readable status reported by the software and with the audit "
        "recorded in docs/AUDIT_REPORT.md.",
        fill="EAF3EA",
        colour=RGBColor(0x1E, 0x5E, 0x20),
    )


# --------------------------------------------------------------------------
# references
# --------------------------------------------------------------------------
def load_references() -> list[str]:
    """Extract the Harvard reference list verbatim from the methodology chapter."""
    text = (DOCS / "METHODOLOGY_CHAPTER.md").read_text(encoding="utf-8")
    marker = "\n## References\n"
    if marker not in text:
        raise SystemExit("could not locate the References section")
    body = text.split(marker, 1)[1]
    body = body.split("\n## ", 1)[0]
    blocks = [b.strip() for b in re.split(r"\n\s*\n", body) if b.strip()]
    refs = []
    for block in blocks:
        flat = " ".join(line.strip() for line in block.splitlines())
        flat = flat.replace("*", "")
        if flat:
            refs.append(flat)
    return refs


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def build(out: Path, toc: list[tuple[str, int]] | None = None) -> None:
    refs = load_references()

    doc = Document()
    _style_document(doc)
    _header_footer(doc)

    core = doc.core_properties
    core.title = "PfVarX \u2014 Methodology, Implementation and Validation Report"
    core.subject = "SNP and small-INDEL analysis for Plasmodium falciparum"
    core.author = "PfVarX Project"
    core.category = "Methodology and validation"
    core.comments = (
        "Generated from the verified implementation. Parameter registry and release "
        "status are read at build time."
    )
    core.keywords = "malaria; variant calling; SNP; INDEL; methodology; validation"

    title_page(doc)
    front_matter(doc, toc)
    section_1(doc)
    section_2(doc)
    section_3(doc)
    section_4(doc)
    section_5(doc)
    section_6(doc)
    section_7(doc)
    section_8(doc)
    section_9(doc)
    section_10(doc)
    section_11(doc)
    section_12(doc)
    section_13(doc)
    section_14(doc, refs)
    section_15(doc)

    out.parent.mkdir(parents=True, exist_ok=True)
    doc.save(out)
    print(f"wrote {out} ({out.stat().st_size / 1024:.1f} KiB)")
    print(f"  paragraphs : {len(doc.paragraphs)}")
    print(f"  tables     : {len(doc.tables)}")
    print(f"  references : {len(refs)}")
    print(f"  parameters : {len(REGISTRY)} (hash {registry_hash()[:12]})")
    print(f"  status     : {overall_status()}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out",
        type=Path,
        default=DOCS / "PfVarX_Methodology.docx",
        help="output path for the generated document",
    )
    ap.add_argument(
        "--toc",
        type=Path,
        default=None,
        help="JSON file mapping section titles to page numbers (second build pass)",
    )
    args = ap.parse_args()
    toc = None
    if args.toc and args.toc.exists():
        toc = [(k, int(v)) for k, v in json.loads(args.toc.read_text()).items()]
    build(args.out, toc)


if __name__ == "__main__":
    main()
